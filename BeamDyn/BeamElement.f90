MODULE BeamElement

      USE GlobalDataFun
      IMPLICIT NONE

      PRIVATE

      PUBLIC

      CONTAINS
      
      
      SUBROUTINE ElementMatrx(Nuu0,Nuuu,Nrr0,Nrrr,hhp,Stif,cet,&
                              w,node_elem,nelem,norder,dof_node,elk,elf)

      REAL(ReKi),INTENT(IN)::Nuu0(:),Nuuu(:),Nrr0(:),Nrrr(:)
      REAL(ReKi),INTENT(IN)::hhp(:,:),Stif(:,:),cet
      REAL(ReKi),INTENT(IN)::w(:)
      INTEGER,INTENT(IN)::node_elem,nelem,norder,dof_node

      REAL(ReKi),INTENT(OUT),ALLOCATABLE::elk(:,:),elf(:)      

      REAL(ReKi),ALLOCATABLE::Fc_elem(:,:),Fd_elem(:,:),Oe_elem(:,:,:)
      REAL(ReKi),ALLOCATABLE::Pe_elem(:,:,:),Qe_elem(:,:,:),Se_elem(:,:,:)
      REAL(ReKi)::E10(3),RR0(3,3),kapa(3),E1(3)
      REAL(ReKi)::Fc(6),Fd(6),Oe(6,6),Pe(6,6),Qe(6,6)

      INTEGER::i,j,k,m,n,temp_id1,temp_id2
      
!      CALL NodalRelRot(Nuu0,node_elem,nelem,norder,dof_node,Nrr0)
!      CALL NodalRelRot(Nuuu,node_elem,nelem,norder,dof_node,Nrrr)

      DO i=1,node_elem
          CALL NodalDataAt0(node_elem,nelem,norder,dof_node,i,hhp,Nuu0,E10)
          CALL NodalData(Nuuu,Nrrr,Nuu0,Nrr0,E10,hhp,&
                         &node_elem,nelem,norder,dof_node,E1,RR0,kapa,Stif)
          CALL ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)
          Fc_elem(i,1:6) = Fc(1:6)
          Fd_elem(i,1:6) = Fd(1:6)
          Oe_elem(i,1:6,1:6) = Oe(1:6,1:6)
          Pe_elem(i,1:6,1:6) = Pe(1:6,1:6)
          Qe_elem(i,1:6,1:6) = Qe(1:6,1:6)
          Se_elem(i,1:6,1:6) = Stif(1:6,1:6)
      ENDDO

      DO i=1,node_elem
          DO j=1,6
              temp_id1 = (i-1)*dof_node+j
              DO k=1,6
                  temp_id2 = (i-1)*dof_node+k
                  elk(temp_id1,temp_id2) = w(i)*Qe_elem(i,j,k)*Jac
              ENDDO
          ENDDO
      ENDDO

      DO i=1,node_elem
          DO j=1,node_elem
              DO k=1,6
                  temp_id1=(i-1)*dof_node+k
                  DO m=1,6
                      temp_id2=(i-1)*dof_node+m
                      elk(temp_id1,temp_id2)=elk(temp_id1,temp_id2)+w(i)*Pe_elem(i,k,m)*hhp(j,i)
                      elk(temp_id1,temp_id2)=elk(temp_id1,temp_id2)+w(j)*Oe_elem(j,k,m)*hhp(i,j)
                  ENDDO
              ENDDO
          ENDDO
      ENDDO

      DO i=1,node_elem
          DO j=1,node_elem
              DO k=1,6
                  temp_id1=(i-1)*dof_node+k
                  DO m=1,6
                      temp_id2=(j-1)*dof_node+m
                      DO n=1,node_elem
                          elk(temp_id1,temp_id2)=elk(temp_id1,temp_id2)+w(n)*hhp(i,n)*Se_elem(n,k,m)*hhp(j,n)/Jac
                      ENDDO
                  ENDDO
              ENDDO
          ENDDO
      ENDDO

      DO i=1,node_elem
          DO j=1,6
              temp_id1 = (i-1)*dof_node+j
              elf(temp_id1) = -w(i)*Fd_elem(i,j)*Jac
              DO m=1,node_elem
                  elf(temp_id1) = elf(temp_id1)-w(m)*hhp(i,m)*Fc_elem(m,j)
              ENDDO
          ENDDO
      ENDDO

      END SUBROUTINE ElementMatrix


     

      SUBROUTINE InertialForce(m00,mEta,rho,vvv,aaa,Fi,Mi,Gi,Ki)

!------------------------------------------------------------------------------
!   INPUT VARIABLES
!------------------------------------------------------------------------------

      REAL(ReKi),INTENT(IN)::m00      ! this variable is 
      REAL(ReKi),INTENT(IN)::mEta(:)  ! this variable is 
      REAL(ReKi),INTENT(IN)::rho(:,:) ! this variable is 
      REAL(ReKi),INTENT(IN)::vvv(:)   ! this variable is
      REAL(ReKi),INTENT(IN)::aaa(:)   ! this variable is
      REAL(ReKi),INTENT(OUT)::Fi(6),Mi(6,6),Gi(6,6),Ki(6,6)

!------------------------------------------------------------------------------ 
!   LOCAL VARIABLES
!------------------------------------------------------------------------------

      REAL(ReKi)::beta(3),gama(3),nu(3),epsi(3,3),mu(3,3)
      REAL(ReKi)::ome(3),omd(3),tempV(3),tempA(3)

!------------------------------------------------------------------------------ 
!   LOGIC
!------------------------------------------------------------------------------

      ! Prepare
      ome(:) = vvv(4:6)
      omd(:) = aaa(4:6)
      tempV(:) = vvv(1:3)
      tempA(:) = aaa(1:3)

      beta = MATMUL(Tilde(ome),mEta)
      gama = MATMUL(rho,ome)
      nu = MATMUL(rho,omd)

      !Compute Fi
      Fi(1:3) = m00*tempA + MATMUL(Tilde(omd),mEta) + MATMUL(Tilde(ome),beta)
      Fi(4:6) = MATMUL(Tilde(mEta,tempA) + nu + MATMUL(Tilde(ome),gama) 

      !Mass Matrix
      Mi = ZERO
      Mi(1:3,1:3) = m00*I3(3,3)
      Mi(1:3,4:6) = TRANSPOSE(Tilde(mEta))
      Mi(4:6,1:3) = Tilde(mEta)
      Mi(4:6,4:6) = rho

      !Gyroscopic Matrix
      Gi = ZERO
      epsi = MATMUL(Tilde(ome),rho)
      mu = MATMUL(Tilde(ome),TRANSPOSE(Tilde(mEta))
      Gi(1:3,4:6) = TRANSPOSE(Tilde(beta)) + mu
      Gi(4:6,4:6) = epsi - Tilde(gama)

      !Stiffness Matrix
      Ki = ZERO
      Ki(1:3,4:6) = MATMUL(Tilde(omd),TRANSPOSE(Tilde(mEta))) +&
                   &MATMUL(Tilde(ome),mu)
      Ki(4:6,4:6) = MATMUL(Tilde(tempA),Tilde(mEta) + &
                   &MATMUL(rho,Tilde(ome)) - Tilde(nu) +&
                   &MATMUL(epsi,Tilde(ome)) - &
                   &MATMUL(Tilde(ome),Tilde(gama))

      END SUBROUTINE InertialForce


      SUBROUTINE NodalDataAt0(node_elem,nelem,norder,dof_node,nnode,hhp,Nuu0,E10)

      INTEGER,INTENT(IN)::node_elem,nelem,norder,dof_node,nnode
      REAL(ReKi),INTENT(IN)::hhp(:,:),Nuu0(:)
      REAL(ReKi),INTENT(OUT)::E10(:)
!      REAL(ReKi),INTENT(INOUT)::Nrr0(:)

      INTEGER:i,temp_id
!      REAL(ReKi)::Nrr0_temp(3),Nrr0_temp1(3)

      E10 = 0.0d0
      DO i=1,node_elem
          temp_id = ((nelem-1)*norder+i-1)*dof_node
          E10(1) = E10(1) + hhp(i,nnode)*Nuu0(temp_id+1)
          E10(2) = E10(2) + hhp(i,nnode)*Nuu0(temp_id+2)
          E10(3) = E10(3) + hhp(i,nnode)*Nuu0(temp_id+3)
      ENDDO


      END SUBROUTINE NodalDataAt0


      SUBROUTINE NodalRelRot(Nuuu,node_elem,nelem,norder,dof_node,Nrrr)

      REAL(ReKi),INTENT(IN)::Nuuu(:)
      INTEGER,INTENT(IN)::node_elem,nelem,norder,dof_node

      REAL(ReKi),INTENT(OUT)::Nrrr(:)

      INTEGER::j,k,temp_id,temp_id1
      REAL(ReKi)::Nuuu_temp1(3),Nuuu_temp(3),Nrrr_temp(3)
      
      DO j=1,node_elem
          temp_id = ((nelem-1)*norder+j-1)*dof_node
          DO k=1,3
              IF(j==1) Nuuu_temp1(k) = Nuuu(temp_id+k+3)
              Nuuu_temp(k) = Nuuu(temp_id+k+3)
          ENDDO
          CALL CrvCompose(Nrrr_temp,Nuuu_temp1,Nuuu_temp,1)
          DO k=1,3
              temp_id1 = (j-1)*3+k
              Nrrr(temp_id1) = Nrrr_temp(k)
          ENDDO
      ENDDO

      END SUBROUTINE NodalRelRot
  

      SUBROUTINE NodalData(Nuuu,Nrrr,Nuu0,Nrr0,E10,hhp,&
                           &node_elem,nelem,norder,dof_node,&
                           &E1,RR0,kapa,Stif)
 
      REAL(ReKi),INTENT(IN)::Nuuu(:),Nrrr(:),Nuu0(:),Nrr0(:)
      REAL(ReKi),INTENT(IN)::E10(:),hhp(:,:)
      INTEGER,INTENT(IN)::node_elem,nelem,norder,dof_node

      REAL(ReKi),INTENT(OUT)::E1(:),RR0(:,:),kapa(:)
      REAL(ReKi),INTENT(INOUT)::Stif(:,:)

      REAL(ReKi)::cc0(3),ccc(3),cc(3),rrr(3),tempR(3,3),tempR6(6,6)
      REAL(ReKi)::Nuuu_temp(3),tempH(3,3),rrp(3),Nuuu_temp1(3)
      INTEGER::i,temp_id
 
!      CALL NodalRelRot(Nuuu,node_elem,nelem,norder,dof_node,Nrrr) 
      DO i=1,node_elem
          temp_id = ((nelem-1)*norder+i-1)*dof_node
          E1(1) = E1(1) + hhp(i,nnode)*Nuuu(temp_id+1)
          E1(2) = E1(2) + hhp(i,nnode)*Nuuu(temp_id+2)
          E1(3) = E1(3) + hhp(i,nnode)*Nuuu(temp_id+3)
          
          temp_id = (i-1)*3
          rrp(1) = rrp(1) + hhp(i,nnode)*Nrrr(temp_id+1)
          rrp(2) = rrp(2) + hhp(i,nnode)*Nrrr(temp_id+2)
          rrp(3) = rrp(3) + hhp(i,nnode)*Nrrr(temp_id+3) 
      ENDDO
      E1 = E1 + E10

!      CALL NodalRelRot(Nuu0,node_elem,nelem,norder,dof_node,Nrr0)
      temp_id = ((nelem-1)*norder+nnode-1)*dof_node+3
      DO i=1,3
          cc0(i) = Nuu0(temp_id+i)
          ccc(i) = Nuuu(temp_id+i)
      ENDDO
      CALL CrvCompose(cc,ccc,cc0,0)
      CALL CrvMatrixR(cc,RR0)

      temp_R6 = 0.0d0
      DO i=1,3
          DO j=1,3
              tempR6(i,j) = RR0(i,j)
              tempR6(i+3,j+3) = RR0(i,j)
          ENDDO
      ENDDO

      Stif = MATMUL(tempR6,MATMUL(Stif,TRANSPOSE(tempR6))

      temp_id = (nnode-1)*3
      DO i=1,3
          rrr(i) = Nrrr(temp_id+i)
      ENDDO
     
      temp_id = ((nelem-1)*norder*dof_node+3
      DO i=1,3
          Nuuu_temp1(i) = Nuuu(temp_id+i)
      ENDDO
      CALL CrvmatrixH(rrr,tempH)
      cc = MATMUL(tempH,rrp)
      CALL CrvMatrixR(Nuuu_temp1,tempR)
      kapa = MATMUL(tempR,cc)
      
      END SUBROUTINE NodalData


      SUBROUTINE ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)

      REAL(ReKi),INTENT(IN)::E1(:),RR0(:,:),kapa(:)
      REAL(ReKi),INTENT(IN)::Stif(:,:),cet
      REAL(ReKi),INTENT(OUT)::Fc(:),Fd(:)    
      REAL(ReKi),INTENT(OUT)::Oe(:,:),Pe(:,:),Qe(:,:) 
 
      REAL(ReKi)::eee(6),fff(6)
      REAL(ReKi)::tempS(3),tempK(3)
      REAL(ReKi)::Wrk(3),e1s,k1s,Wrk33(3,3)
      REAL(ReKi)::C11(3,3),C12(3,3),C21(3,3),C22(3,3)

      INTEGER:: i,j

      eee = ZERO 
      DO i=1,3
          eee(i) = E1(i) - RR0(i,1)
          eee(i+3) = kapa(i)

          tempS(i) = eee(i)
          tempK(i) = eee(i+3)
      ENDDO
           
      fff = ZERO 
      fff = MATMUL(Stif,eee)

      Wrk = ZERO     
      Wrk = MATMUL(TRANSPOSE(RR0),tempS)
      e1s = Wrk(1)      !epsilon_{11} in material basis

      Wrk = ZERO
      Wrk = MATMUL(TRANSPOSE(RR0),tempK)
      k1s = Wrk(1)      !kapa_{1} in material basis
     
      DO i=1,3
          fff(i) = fff(i) + HALF*cet*k1s*k1s*RR0(i,1)
          fff(i+3) = fff(i+3) + cet*e1s*k1s*RR0(i,1)
      ENDDO 

      Fc = ZERO
      Fc = fff
      Wrk = ZERO 
      Wrk(1:3) = fff(1:3)
      Fd = ZERO 
      Fd(4:6) = MATMUL(Tilde(Wrk),E1)

      C11(:,:) = Stif(1:3,1:3)
      C12(:,:) = Stif(1:3,4:6)
      C21(:,:) = Stif(4:6,1:3)
      C22(:,:) = Stif(4:6,4:6)

      DO i=1,3
          Wrk(i) = RR0(i,1)
      ENDDO
      Wrk33 = OuterProduct(Wrk,Wrk) 
      C12 = C12 + cet*k1s*Wrk33
      C21 = C21 + cet*k1s*Wrk33
      C22 = C22 + cet*e1s*Wrk33

      epsi = ZERO 
      mu = ZERO
      epsi = MATMUL(C11,Tilde(E1))
      mu = MATMUL(C21,Tilde(E1))
   
      Wrk = ZERO

      Oe = ZERO
      Oe(1:3,4:6) = epsi(1:3,1:3)
      Oe(4:6,4:6) = mu(1:3,1:3)
      
      Wrk(:) = fff(1:3)
      Oe(1:3,4:6) = Oe(1:3,4:6) - Tilde(Wrk)
      Wrk(:) = fff(4:6)
      Oe(4:6,4:6) = Oe(4:6,4:6) - Tilde(Wrk)

      Pe = ZERO
      Wrk(:) = fff(1:3)
      Pe(4:6,1:3) = Tilde(Wrk) + TRANSPOSE(epsi)
      Pe(4:6,4:6) = TRANSPOSE(mu)

      Qe = ZERO
      Wrk33(1:3,1:3) = Oe(1:3,4:6)
      Qe(4:6,4:6) = MATMUL(TRANSPOSE(Tilde(E1)),Wrk33)

      END SUBROUTINE ElasticForce 
