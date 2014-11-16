   SUBROUTINE ElementMatrix_AM2(Nuu0,Nuuu,Nuuu0,Nrr0,Nrrr,Nrrr0,Nvvv,Nvvv0,Nuud0,Nvvd0,&
                               &EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
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
   REAL(ReKi),INTENT(IN   )    :: Nuud0(:) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi),INTENT(IN   )    :: Nvvd0(:) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi),INTENT(IN   )    :: EStif0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),INTENT(IN   )    :: EMass0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),INTENT(IN   )    :: gravity(:) ! 
   REAL(ReKi),INTENT(IN   )    :: DistrLoad_GL(:,:) ! Nodal material properties for each element
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
   REAL(ReKi),INTENT(IN   )    :: beta(:)

   REAL(ReKi)                  :: gp(ngp) ! Gauss points
   REAL(ReKi)                  :: gw(ngp) ! Gauss point weights
   REAL(ReKi)                  :: hhx(node_elem) ! Shape function
   REAL(ReKi)                  :: hpx(node_elem) ! Derivative of shape function
   REAL(ReKi)                  :: GLL_temp(node_elem) ! Temp Gauss-Lobatto-Legendre points
   REAL(ReKi)                  :: w_temp(node_elem) ! Temp GLL weights
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
   REAL(ReKi)                  :: uup0(3)
   REAL(ReKi)                  :: Jacobian
   REAL(ReKi)                  :: gpr
   REAL(ReKi)                  :: Fc(6)
   REAL(ReKi)                  :: Fd(6)
   REAL(ReKi)                  :: Fg(6)
   REAL(ReKi)                  :: vvv(6)
   REAL(ReKi)                  :: vvv0(6)
   REAL(ReKi)                  :: vvp(6)
   REAL(ReKi)                  :: vvp0(6)
   REAL(ReKi)                  :: uud0(6)
   REAL(ReKi)                  :: vvd0(6)
   REAL(ReKi)                  :: mmm
   REAL(ReKi)                  :: mEta(3)
   REAL(ReKi)                  :: rho(3,3)
   REAL(ReKi)                  :: Fb(6)
   REAL(ReKi)                  :: Mi(6,6)
   REAL(ReKi)                  :: A1(6,6)
   REAL(ReKi)                  :: F1(6)
   REAL(ReKi)                  :: A2(6,6)
   REAL(ReKi)                  :: A3(6,6)
   REAL(ReKi)                  :: A4(6,6)
   REAL(ReKi)                  :: A5(6,6)
   REAL(ReKi)                  :: A6(6,6)
   REAL(ReKi)                  :: A7(6,6)
   REAL(ReKi)                  :: Oe(6,6)
   REAL(ReKi)                  :: Pe(6,6)
   REAL(ReKi)                  :: Qe(6,6)
   REAL(ReKi)                  :: temp_H(3,3)
   REAL(ReKi)                  :: temp_H0(3,3)
   REAL(ReKi)                  :: Sd(6,6)
   REAL(ReKi)                  :: Od(6,6)
   REAL(ReKi)                  :: Pd(6,6)
   REAL(ReKi)                  :: Qd(6,6)
   REAL(ReKi)                  :: betaC(6,6)
   REAL(ReKi)                  :: Gd(6,6)
   REAL(ReKi)                  :: Xd(6,6)
   REAL(ReKi)                  :: Yd(6,6)
   INTEGER(IntKi)              :: igp
   INTEGER(IntKi)              :: i
   INTEGER(IntKi)              :: j
   INTEGER(IntKi)              :: m
   INTEGER(IntKi)              :: n
   INTEGER(IntKi)              :: temp_id1
   INTEGER(IntKi)              :: temp_id2


   elf1(:)    = 0.0D0
   elf2(:)    = 0.0D0
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
       CALL BldGaussPointData(hhx,hpx,Nuuu0,Nrrr0,uu0,E10,node_elem,dof_node,uuu0,uup0,E100,RR00,kapa0,Stif0,cet0)       
       CALL BldGaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,uuu,uup,E1,RR0,kapa,Stif,cet)       
!WRITE(*,*) uuu0(1:6)
       mmm  = 0.0D0
       mEta = 0.0D0
       rho  = 0.0D0
       mmm          =  EMass0_GL(1,1,igp)
       mEta(2)      = -EMass0_GL(1,6,igp)
       mEta(3)      =  EMass0_GL(1,5,igp)
       rho(1:3,1:3) =  EMass0_GL(4:6,4:6,igp)
       CALL BldGaussPointDataXdot(hhx,Nvvv0,node_elem,dof_node,vvv0)
       CALL BldGaussPointDataXdot(hhx,Nuud0,node_elem,dof_node,uud0)
       CALL BldGaussPointDataXdot(hhx,Nvvd0,node_elem,dof_node,vvd0)
       CALL BldGaussPointDataMass(hhx,hpx,Nvvv,RR0,node_elem,dof_node,vvv,vvp,mmm,mEta,rho)
       CALL BeamDyn_StaticElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)
       IF(damp_flag .NE. 0) THEN
           CALL DissipativeForce(beta,Stif,vvv,vvp,E1,Fc,Fd,Sd,Od,Pd,Qd,betaC,Gd,Xd,Yd)
       ENDIF

       CALL MassMatrix(mmm,mEta,rho,Mi)
       CALL GyroForce(mEta,rho,uuu,vvv,Fb)
       CALL GravityLoads(mmm,mEta,gravity,Fg)
       CALL AM2LinearizationMatrix(uuu,vvv,uuu0,vvv0,uud0,vvd0,mmm,mEta,rho,dt,A1,A2,A3,A4,A5,A6,A7)

       Fd(:) = MATMUL(Mi,vvv-vvv0-0.5D0*dt*vvd0) + &
              &0.5D0*dt*Fb + &
              &0.5D0*dt*Fd - &
              &0.5D0*dt*DistrLoad_GL(:,igp)

       Fc(:) = 0.5D0*dt*Fc

       CALL CrvMatrixH(uuu(4:6),temp_H)
       F1(1:3) = uuu(1:3) - uuu0(1:3) - 0.5D0*dt*(vvv(1:3)+vvv0(1:3))
       F1(4:6) = MATMUL(temp_H,uuu(4:6)-uuu0(4:6)) - &
                &0.5D0*dt*MATMUL(temp_H,uud0(4:6)) - &
                &0.5D0*dt*vvv(4:6)
!WRITE(*,*) F1(1:6)
       DO i=1,node_elem
           DO j=1,node_elem
               DO m=1,dof_node
                   temp_id1 = (i-1)*dof_node+m
                   DO n=1,dof_node
                       temp_id2 = (j-1)*dof_node+n
                       elm11(temp_id1,temp_id2) = elm11(temp_id1,temp_id2) + hhx(i)*A6(m,n)*hhx(j)*Jacobian*gw(igp)
                       elm12(temp_id1,temp_id2) = elm12(temp_id1,temp_id2) + hhx(i)*A7(m,n)*hhx(j)*Jacobian*gw(igp)
                       elm21(temp_id1,temp_id2) = elm21(temp_id1,temp_id2) + &
                                                 &hhx(i)*(A2(m,n)-A3(m,n)+0.5D0*dt*(-A1(m,n)+A5(m,n)+Qe(m,n)))*hhx(j)*Jacobian*gw(igp) + &
                                                 &hhx(i)*0.5D0*dt*Pe(m,n)*hpx(j)*Jacobian*gw(igp) + &
                                                 &hpx(i)*0.5D0*dt*Stif(m,n)*hpx(j)*Jacobian*gw(igp) + &
                                                 &hpx(i)*0.5D0*dt*Oe(m,n)*hhx(j)*Jacobian*gw(igp)
                       elm22(temp_id1,temp_id2) = elm22(temp_id1,temp_id2) + &
                                                 &hhx(i)*(Mi(m,n)+0.5D0*dt*A4(m,n))*hhx(j)*Jacobian*gw(igp)
                       IF(damp_flag .NE. 0) THEN
                           elm21(temp_id1,temp_id2) = elm21(temp_id1,temp_id2) + &
                                                     &hhx(i)*0.5D0*dt*Pd(m,n)*hpx(j)*Jacobian*gw(igp) + &
                                                     &hhx(i)*0.5D0*dt*Qd(m,n)*hhx(j)*Jacobian*gw(igp) + &
                                                     &hpx(i)*0.5D0*dt*Sd(m,n)*hpx(j)*Jacobian*gw(igp) + &
                                                     &hpx(i)*0.5D0*dt*Od(m,n)*hhx(j)*Jacobian*gw(igp)  
                           elm22(temp_id1,temp_id2) = elm22(temp_id1,temp_id2) + &
                                                     &hhx(i)*0.5D0*dt*Xd(m,n)*hhx(j)*Jacobian*gw(igp) + &
                                                     &hhx(i)*0.5D0*dt*Yd(m,n)*hpx(j)*Jacobian*gw(igp) + &
                                                     &hpx(i)*0.5D0*dt*Gd(m,n)*hhx(j)*Jacobian*gw(igp) + &
                                                     &hpx(i)*0.5D0*dt*betaC(m,n)*hpx(j)*Jacobian*gw(igp) 
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


   END SUBROUTINE ElementMatrix_AM2
