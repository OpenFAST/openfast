   SUBROUTINE ElementMatrix_AM2(Nuu0,Nuuu,Nuuu0,Nrr0,Nrrr,Nrrr0,Nvvv,Nvvv0,&
                               &EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,DistrLoad_GL0,&
                               &ngp,dt,node_elem,dof_node,damp_flag,beta,&
                               &elf1,elf2,elm11,elm12,elm21,elm22)
                           
   !-------------------------------------------------------------------------------
   ! This subroutine total element forces and mass matrices
   !-------------------------------------------------------------------------------

   REAL(ReKi),INTENT(IN   )    :: Nuu0(:) ! Nodal initial position for each element
   REAL(ReKi),INTENT(IN   )    :: Nuuu(:) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi),INTENT(IN   )    :: Nuuu0(:) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi),INTENT(IN   )    :: Nrr0(:) ! Nodal rotation parameters for initial position
   REAL(ReKi),INTENT(IN   )    :: Nrrr(:) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi),INTENT(IN   )    :: Nrrr0(:) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi),INTENT(IN   )    :: Nvvv(:) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi),INTENT(IN   )    :: Nvvv0(:) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi),INTENT(IN   )    :: EStif0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),INTENT(IN   )    :: EMass0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),INTENT(IN   )    :: gravity(:) ! 
   REAL(ReKi),INTENT(IN   )    :: DistrLoad_GL(:,:) ! Nodal material properties for each element
   REAL(ReKi),INTENT(IN   )    :: DistrLoad_GL0(:,:) ! Nodal material properties for each element
   REAL(DbKi),INTENT(IN   )    :: dt ! Nodal material properties for each element
   REAL(ReKi),INTENT(  OUT)    :: elf1(:)  ! Total element force (Fd, Fc, Fb)
   REAL(ReKi),INTENT(  OUT)    :: elf2(:)  ! Total element force (Fd, Fc, Fb)
   REAL(ReKi),INTENT(  OUT)    :: elm11(:,:) ! Total element mass matrix
   REAL(ReKi),INTENT(  OUT)    :: elm12(:,:) ! Total element mass matrix
   REAL(ReKi),INTENT(  OUT)    :: elm21(:,:) ! Total element mass matrix
   REAL(ReKi),INTENT(  OUT)    :: elm22(:,:) ! Total element mass matrix
   INTEGER(IntKi),INTENT(IN   ):: ngp ! Number of Gauss points
   INTEGER(IntKi),INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi),INTENT(IN   ):: dof_node ! Degrees of freedom per node
   INTEGER(IntKi),INTENT(IN   ):: damp_flag ! Degrees of freedom per node
   REAL(ReKi),INTENT(IN   )    :: beta

!   REAL(ReKi),ALLOCATABLE      :: gp(:) ! Gauss points
!   REAL(ReKi),ALLOCATABLE      :: gw(:) ! Gauss point weights
!   REAL(ReKi),ALLOCATABLE      :: hhx(:) ! Shape function
!   REAL(ReKi),ALLOCATABLE      :: hpx(:) ! Derivative of shape function
!   REAL(ReKi),ALLOCATABLE      :: GLL_temp(:) ! Temp Gauss-Lobatto-Legendre points
!   REAL(ReKi),ALLOCATABLE      :: w_temp(:) ! Temp GLL weights
   REAL(ReKi)      :: gp(ngp) ! Gauss points
   REAL(ReKi)      :: gw(ngp) ! Gauss point weights
   REAL(ReKi)      :: hhx(node_elem) ! Shape function
   REAL(ReKi)      :: hpx(node_elem) ! Derivative of shape function
   REAL(ReKi)      :: GLL_temp(node_elem) ! Temp Gauss-Lobatto-Legendre points
   REAL(ReKi)      :: w_temp(node_elem) ! Temp GLL weights
   REAL(ReKi)                  :: uu0(6)
   REAL(ReKi)                  :: E10(3)
   REAL(ReKi)                  :: RR0(3,3)
   REAL(ReKi)                  :: RR00(3,3)
   REAL(ReKi)                  :: kapa(3)
   REAL(ReKi)                  :: kapa0(3)
   REAL(ReKi)                  :: E1(3)
   REAL(ReKi)                  :: E100(3)
   REAL(ReKi)                  :: Stif(6,6)
   REAL(ReKi)                  :: Stif0(6,6)
   REAL(ReKi)                  :: cet
   REAL(ReKi)                  :: cet0
   REAL(ReKi)                  :: uuu(6)
   REAL(ReKi)                  :: uuu0(6)
   REAL(ReKi)                  :: uup(3)
   REAL(ReKi)                  :: Jacobian
   REAL(ReKi)                  :: gpr
   REAL(ReKi)                  :: Fc(6)
   REAL(ReKi)                  :: Fc0(6)
   REAL(ReKi)                  :: Fd(6)
   REAL(ReKi)                  :: Fd0(6)
   REAL(ReKi)                  :: Fg(6)
   REAL(ReKi)                  :: Fg0(6)
   REAL(ReKi)                  :: vvv(6)
   REAL(ReKi)                  :: vvv0(6)
   REAL(ReKi)                  :: vvp(6)
   REAL(ReKi)                  :: vvp0(6)
   REAL(ReKi)                  :: mmm
   REAL(ReKi)                  :: mEta(3)
   REAL(ReKi)                  :: rho(3,3)
   REAL(ReKi)                  :: mEta0(3)
   REAL(ReKi)                  :: rho0(3,3)
   REAL(ReKi)                  :: Fb(6)
   REAL(ReKi)                  :: Fb0(6)
   REAL(ReKi)                  :: Mi(6,6)
   REAL(ReKi)                  :: Mi0(6,6)
!   REAL(ReKi)                  :: A1(6,6)
   REAL(ReKi)                  :: F1(6)
   REAL(ReKi)                  :: A2(6,6)
   REAL(ReKi)                  :: A3(6,6)
   REAL(ReKi)                  :: A4(6,6)
   REAL(ReKi)                  :: A5(6,6)
   REAL(ReKi)                  :: A6(6,6)
   REAL(ReKi)                  :: A7(6,6)
   REAL(ReKi)                  :: Oe(6,6)
   REAL(ReKi)                  :: Oe0(6,6)
   REAL(ReKi)                  :: Pe(6,6)
   REAL(ReKi)                  :: Pe0(6,6)
   REAL(ReKi)                  :: Qe(6,6)
   REAL(ReKi)                  :: Qe0(6,6)
   REAL(ReKi)                  :: temp_H(3,3)
   REAL(ReKi)                  :: temp_H0(3,3)
   REAL(ReKi)                  :: uup0(3)
   REAL(ReKi)                  :: Sd(6,6)
   REAL(ReKi)                  :: Sd0(6,6)
   REAL(ReKi)                  :: Od(6,6)
   REAL(ReKi)                  :: Od0(6,6)
   REAL(ReKi)                  :: Pd(6,6)
   REAL(ReKi)                  :: Pd0(6,6)
   REAL(ReKi)                  :: Qd(6,6)
   REAL(ReKi)                  :: Qd0(6,6)
   REAL(ReKi)                  :: betaC(6,6)
   REAL(ReKi)                  :: betaC0(6,6)
   REAL(ReKi)                  :: Gd(6,6)
   REAL(ReKi)                  :: Gd0(6,6)
   REAL(ReKi)                  :: Xd(6,6)
   REAL(ReKi)                  :: Xd0(6,6)
   REAL(ReKi)                  :: Yd(6,6)
   REAL(ReKi)                  :: Yd0(6,6)

   INTEGER(IntKi)              :: igp
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: m
   INTEGER(IntKi)              :: n
   INTEGER(IntKi)              :: temp_id1
   INTEGER(IntKi)              :: temp_id2
   INTEGER(IntKi)              :: allo_stat

! WRITE(*,*) "TEST2"
!i   ALLOCATE(gp(ngp), STAT = allo_stat)
!i   IF(allo_stat/=0) GOTO 9999
!i   gp = 0.0D0
!i
!i   ALLOCATE(gw(ngp), STAT = allo_stat)
!i   IF(allo_stat/=0) GOTO 9999
!i   gw = 0.0D0
!i
!i   ALLOCATE(hhx(node_elem), STAT = allo_stat)
!i   IF(allo_stat/=0) GOTO 9999
!i   hhx = 0.0D0
!i
!i   ALLOCATE(hpx(node_elem), STAT = allo_stat)
!i   IF(allo_stat/=0) GOTO 9999
!i   hpx = 0.0D0
!i   
!i   ALLOCATE(GLL_temp(node_elem), STAT = allo_stat)
!i   IF(allo_stat/=0) GOTO 9999
!i   GLL_temp = 0.0D0
!i   
!i   ALLOCATE(w_temp(node_elem), STAT = allo_stat)
!i   IF(allo_stat/=0) GOTO 9999
!i   w_temp = 0.0D0

   elf1(:) = 0.0D0
   elf2(:) = 0.0D0
   elm11(:,:) = 0.0D0
   elm12(:,:) = 0.0D0
   elm21(:,:) = 0.0D0
   elm22(:,:) = 0.0D0

   CALL BD_gen_gll_LSGL(ngp,GLL_temp,w_temp)
   CALL BldGaussPointWeight(ngp,gp,gw)

   DO igp=1,ngp
       gpr=gp(igp)
       CALL BldComputeJacobianLSGL(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian)

       CALL BldGaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10)
       Stif(:,:) = 0.0D0
       Stif(1:6,1:6) = EStif0_GL(1:6,1:6,igp)
       CALL BldGaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,uuu,uup,E1,RR0,kapa,Stif,cet)       
       Stif0(:,:) = 0.0D0
       Stif0(1:6,1:6) = EStif0_GL(1:6,1:6,igp)
       CALL BldGaussPointData(hhx,hpx,Nuuu0,Nrrr0,uu0,E10,node_elem,dof_node,uuu0,uup0,E100,RR00,kapa0,Stif0,cet0)       
       mmm  = 0.0D0
       mEta = 0.0D0
       rho  = 0.0D0
       mmm          = EMass0_GL(1,1,igp)
       mEta(2)      = -EMass0_GL(1,6,igp)
       mEta(3)      =  EMass0_GL(1,5,igp)
       rho(1:3,1:3) = EMass0_GL(4:6,4:6,igp)
       mEta0(:) = mEta(:)
       rho0(:,:) = rho(:,:)
       CALL BldGaussPointDataMass(hhx,hpx,Nvvv,RR0,node_elem,dof_node,vvv,vvp,mmm,mEta,rho)
       CALL BldGaussPointDataMass(hhx,hpx,Nvvv0,RR00,node_elem,dof_node,vvv0,vvp0,mmm,mEta0,rho0)
       CALL BeamDyn_StaticElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)
       CALL BeamDyn_StaticElasticForce(E100,RR00,kapa0,Stif0,cet0,Fc0,Fd0,Oe0,Pe0,Qe0)
       IF(damp_flag .NE. 0) THEN
           CALL DissipativeForce(beta,Stif,vvv,vvp,E1,Fc,Fd,Sd,Od,Pd,Qd,betaC,Gd,Xd,Yd)
           CALL DissipativeForce(beta,Stif0,vvv0,vvp0,E100,Fc0,Fd0,Sd0,Od0,Pd0,Qd0,betaC0,Gd0,Xd0,Yd0)
       ENDIF

       CALL MassMatrix(mmm,mEta,rho,Mi)
       CALL MassMatrix(mmm,mEta0,rho0,Mi0)
       CALL GyroForce(mEta,rho,uuu,vvv,Fb)
       CALL GyroForce(mEta0,rho0,uuu0,vvv0,Fb0)
       CALL GravityLoads(mmm,mEta,gravity,Fg)
       CALL GravityLoads(mmm,mEta0,gravity,Fg0)
       CALL AM2LinearizationMatrix(uuu,vvv,uuu0,vvv0,mmm,mEta,rho,A2,A3,A4,A5,A6,A7)

       Fd(:) = dt*(Fd + Fd0) + &
              &MATMUL(Mi,vvv-vvv0) + &
              &MATMUL(Mi0,vvv-vvv0) + &
              &dt*(Fb + Fb0) - &
              &dt*(DistrLoad_GL(:,igp)+Fg+DistrLoad_GL0(:,igp)+Fg0)

       Fc(:) = dt*(Fc + Fc0)

       CALL CrvMatrixH(uuu(4:6),temp_H)
       CALL CrvMatrixH(uuu0(4:6),temp_H0)
       F1(1:3) = -dt*(vvv(1:3)+vvv0(1:3)) - 2.0D0*(uuu0(1:3)-uuu(1:3))
       F1(4:6) = MATMUL(temp_H,uuu(4:6)-uuu0(4:6)) + &
                &MATMUL(temp_H0,uuu(4:6)-uuu0(4:6)) - &
                &dt*(vvv(4:6)+vvv0(4:6))

       DO i=1,node_elem
           DO j=1,node_elem
               DO m=1,dof_node
                   temp_id1 = (i-1)*dof_node+m
                   DO n=1,dof_node
                       temp_id2 = (j-1)*dof_node+n
                       elm11(temp_id1,temp_id2) = elm11(temp_id1,temp_id2) + hhx(i)*A6(m,n)*hhx(j)*Jacobian*gw(igp)
                       elm12(temp_id1,temp_id2) = elm12(temp_id1,temp_id2) - dt* hhx(i)*A7(m,n)*hhx(j)*Jacobian*gw(igp)
                       elm21(temp_id1,temp_id2) = elm21(temp_id1,temp_id2) + &
                                                 &hhx(i)*(A2(m,n)-A3(m,n)+dt*A5(m,n)+dt*Qe(m,n))*hhx(j)*Jacobian*gw(igp) + &
                                                 &hhx(i)*dt*Pe(m,n)*hpx(j)*Jacobian*gw(igp) + &
                                                 &hpx(i)*dt*Stif(m,n)*hpx(j)*Jacobian*gw(igp) + &
                                                 &hpx(i)*dt*Oe(m,n)*hhx(j)*Jacobian*gw(igp)
                       elm22(temp_id1,temp_id2) = elm22(temp_id1,temp_id2) + &
                                                 &hhx(i)*(Mi(m,n)+Mi0(m,n)+dt*A4(m,n))*hhx(j)*Jacobian*gw(igp)
                       IF(damp_flag .NE. 0) THEN
                           elm21(temp_id1,temp_id2) = elm21(temp_id1,temp_id2) + &
                                                     &hhx(i)*dt*Pd(m,n)*hpx(j)*Jacobian*gw(igp) + &
                                                     &hhx(i)*dt*Qd(m,n)*hhx(j)*Jacobian*gw(igp) + &
                                                     &hpx(i)*dt*Sd(m,n)*hpx(j)*Jacobian*gw(igp) + &
                                                     &hpx(i)*dt*Od(m,n)*hhx(j)*Jacobian*gw(igp)  
                           elm22(temp_id1,temp_id2) = elm22(temp_id1,temp_id2) + &
                                                     &hhx(i)*dt*Xd(m,n)*hhx(j)*Jacobian*gw(igp) + &
                                                     &hhx(i)*dt*Yd(m,n)*hpx(j)*Jacobian*gw(igp) + &
                                                     &hpx(i)*dt*Gd(m,n)*hhx(j)*Jacobian*gw(igp) + &
                                                     &hpx(i)*dt*betaC(m,n)*hpx(j)*Jacobian*gw(igp) 
                       ENDIF
                   ENDDO
               ENDDO
           ENDDO
       ENDDO
!DO i=1,18
!WRITE(*,*) "elm21",i,elm21(i,1)
!ENDDO
       DO i=1,node_elem
           DO j=1,dof_node
               temp_id1 = (i-1) * dof_node+j
               elf1(temp_id1) = elf1(temp_id1) - hhx(i)*F1(j)*Jacobian*gw(igp)
               elf2(temp_id1) = elf2(temp_id1) - hpx(i)*Fc(j)*Jacobian*gw(igp)
               elf2(temp_id1) = elf2(temp_id1) - hhx(i)*Fd(j)*Jacobian*gw(igp)
           ENDDO
       ENDDO
   ENDDO

!   DEALLOCATE(gp)
!   DEALLOCATE(gw)
!   DEALLOCATE(hhx)
!   DEALLOCATE(hpx)
!   DEALLOCATE(GLL_temp)
!   DEALLOCATE(w_temp)


!   9999 IF(allo_stat/=0) THEN
!            IF(ALLOCATED(gp))  DEALLOCATE(gp)
!            IF(ALLOCATED(gw))  DEALLOCATE(gw)
!            IF(ALLOCATED(hhx)) DEALLOCATE(hhx)
!            IF(ALLOCATED(hpx)) DEALLOCATE(hpx)
!            IF(ALLOCATED(GLL_temp)) DEALLOCATE(GLL_temp)
!            IF(ALLOCATED(w_temp)) DEALLOCATE(w_temp)
!        ENDIF

   END SUBROUTINE ElementMatrix_AM2
