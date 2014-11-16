   SUBROUTINE ElementMatrix_CCSD(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,&
                                &EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                                &ngp,node_elem,dof_node,damp_flag,beta,&
                                &elf1,elf2,elm11,elm22)
                           
   !-------------------------------------------------------------------------------
   ! This subroutine total element forces and mass matrices
   !-------------------------------------------------------------------------------

   REAL(ReKi),INTENT(IN   )    :: Nuu0(:) ! Nodal initial position for each element
   REAL(ReKi),INTENT(IN   )    :: Nuuu(:) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi),INTENT(IN   )    :: Nrr0(:) ! Nodal rotation parameters for initial position
   REAL(ReKi),INTENT(IN   )    :: Nrrr(:) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi),INTENT(IN   )    :: Nvvv(:) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi),INTENT(IN   )    :: EStif0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),INTENT(IN   )    :: EMass0_GL(:,:,:) ! Nodal material properties for each element
   REAL(ReKi),INTENT(IN   )    :: gravity(:) ! 
   REAL(ReKi),INTENT(IN   )    :: DistrLoad_GL(:,:) ! Nodal material properties for each element
   REAL(ReKi),INTENT(  OUT)    :: elf1(:)  ! Total element force (Fd, Fc, Fb)
   REAL(ReKi),INTENT(  OUT)    :: elf2(:)  ! Total element force (Fd, Fc, Fb)
   REAL(ReKi),INTENT(  OUT)    :: elm11(:,:) ! Total element mass matrix
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
   REAL(ReKi)                  :: kapa(3)
   REAL(ReKi)                  :: E1(3)
   REAL(ReKi)                  :: Stif(6,6)
   REAL(ReKi)                  :: cet
   REAL(ReKi)                  :: uuu(6)
   REAL(ReKi)                  :: uup(3)
   REAL(ReKi)                  :: Jacobian
   REAL(ReKi)                  :: gpr
   REAL(ReKi)                  :: Fc(6)
   REAL(ReKi)                  :: Fd(6)
   REAL(ReKi)                  :: Fg(6)
   REAL(ReKi)                  :: vvv(6)
   REAL(ReKi)                  :: vvp(6)
   REAL(ReKi)                  :: mmm
   REAL(ReKi)                  :: mEta(3)
   REAL(ReKi)                  :: rho(3,3)
   REAL(ReKi)                  :: Fb(6)
   REAL(ReKi)                  :: Mi(6,6)
   REAL(ReKi)                  :: A1(6,6)
   REAL(ReKi)                  :: F1(6)
   REAL(ReKi)                  :: Oe(6,6)
   REAL(ReKi)                  :: Pe(6,6)
   REAL(ReKi)                  :: Qe(6,6)
   REAL(ReKi)                  :: temp_H(3,3)
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
   INTEGER(IntKi)              :: allo_stat


   elf1(:) = 0.0D0
   elf2(:) = 0.0D0
   elm11(:,:) = 0.0D0
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
       mmm  = 0.0D0
       mEta = 0.0D0
       rho  = 0.0D0
       mmm          = EMass0_GL(1,1,igp)
       mEta(2)      = -EMass0_GL(1,6,igp)
       mEta(3)      =  EMass0_GL(1,5,igp)
       rho(1:3,1:3) = EMass0_GL(4:6,4:6,igp)
       CALL BldGaussPointDataMass(hhx,hpx,Nvvv,RR0,node_elem,dof_node,vvv,vvp,mmm,mEta,rho)
       CALL MassMatrix(mmm,mEta,rho,Mi)
!DO i=1,6
!WRITE(*,*) i,Mi(i,i)
!ENDDO
       CALL GyroForce(mEta,rho,uuu,vvv,Fb)
       CALL GravityLoads(mmm,mEta,gravity,Fg)

       CALL BeamDyn_StaticElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)
       IF(damp_flag .NE. 0) THEN
           CALL DissipativeForce(beta,Stif,vvv,vvp,E1,Fc,Fd,Sd,Od,Pd,Qd,betaC,Gd,Xd,Yd)
       ENDIF

       A1(:,:) = 0.0D0
       CALL CrvMatrixH(uuu(4:6),temp_H)
       DO i = 1,3
           A1(i,i) = 1.0D0
       ENDDO
       A1(4:6,4:6) = temp_H(1:3,1:3)

       F1(1:6) = vvv(1:6)
       Fd(:) = Fd(:) + Fb(:) - DistrLoad_GL(:,igp) - Fg(:)

       DO i=1,node_elem
           DO j=1,node_elem
               DO m=1,dof_node
                   temp_id1 = (i-1)*dof_node+m
                   DO n=1,dof_node
                       temp_id2 = (j-1)*dof_node+n
                       elm11(temp_id1,temp_id2) = elm11(temp_id1,temp_id2) + hhx(i)*A1(m,n)*hhx(j)*Jacobian*gw(igp)
                       elm22(temp_id1,temp_id2) = elm22(temp_id1,temp_id2) + hhx(i)*Mi(m,n)*hhx(j)*Jacobian*gw(igp)
                   ENDDO
               ENDDO
           ENDDO
       ENDDO

       DO i=1,node_elem
           DO j=1,dof_node
               temp_id1 = (i-1) * dof_node+j
               elf1(temp_id1) = elf1(temp_id1) + hhx(i)*F1(j)*Jacobian*gw(igp)
               elf2(temp_id1) = elf2(temp_id1) - hpx(i)*Fc(j)*Jacobian*gw(igp)
               elf2(temp_id1) = elf2(temp_id1) - hhx(i)*Fd(j)*Jacobian*gw(igp)
           ENDDO
       ENDDO
DO i=1,18
!WRITE(*,*) i,elf2(i)
ENDDO
   ENDDO


   END SUBROUTINE ElementMatrix_CCSD
