      MODULE BeamElement

      USE GlobalDataFun
      IMPLICIT NONE

      PRIVATE

      PUBLIC

      CONTAINS
      
      
      SUBROUTINE ElementMatrx()

      CALL InertialForce(m00,mEta,rho,vvv,aaa,Fi,Mi,Gi,Ki)

      CALL ElasticForce()

      DO i=1,node_elem
          DO j=1,6
              DO m=1,6
                  temp_fun(j,m) = ShapeFun(i,j)
              ENDDO
              elf((i-1)*DOF_NODE+j) = -temp_fun(j,i)*Fi(j)*djw
              elf((i-1)*DOF_NODE+J) = elf - ShapeDer(
          ENDDO
      ENDDO



      SUBROUTINE InertialForce(m00,mEta,rho,vvv,aaa,Fi,Mi,Gi,Ki)

      REAL(DBL),INTENT(IN)::m00,mEta(:),rho(:,:)
      REAL(DBL),INTENT(IN)::vvv(:),aaa(:)
      REAL(DBL),INTENT(OUT)::Fi(6),Mi(6,6),Gi(6,6),Ki(6,6)

      REAL(DBL)::beta(3),gama(3),nu(3),epsi(3,3),mu(3,3)
      REAL(DBL)::ome(3),omd(3),tempV(3),tempA(3)

      ! Prepare
      ome(:) = vvv(4:6)
      omd(:) = aaa(4:6)
      tempV(:) = vvv(1:3)
      tempA(:) = aaa(1:3)

      beta = MATMUL(Tilde(ome),mEta)
      gama = MATMUL(rho,ome)
      nu = MATMUL(rho,omd)

      !Compute Fi
      Fi(1:3)= m00*tempA+MATMUL(Tilde(omd),mEta)+MATMUL(Tilde(ome),beta)
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





      SUBROUTINE ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)

      REAL(DBL),INTENT(IN)::E1(:),RR0(:,:),kapa(:)
      REAL(DBL),INTENT(IN)::Stif(:,:),cet
      REAL(DBL),INTENT(OUT)::Fc(:),Fd(:)    
      REAL(DBL),INTENT(OUT)::Oe(:,:),Pe(:,:),Qe(:,:) 
 
      REAL(DBL)::eee(6),fff(6)
      REAL(DBL)::tempS(3),tempK(3)
      REAL(DBL)::Wrk(3),e1s,k1s,Wrk33(3,3)
      REAL(DBL)::C11(3,3),C12(3,3),C21(3,3),C22(3,3)

      INTEGER:: i,j

      eee = ZERO 
      DO i=1,3
          eee(i) = E1(i) - RR0(i,1)
          eee(i+3) = kapa(i)

          tempS(i) = eee(i)
          tempK(i) = eee(i+3)
      ENDDO
           
 
      fff = MATMUL(Stif,eee)

      Wrk = ZERO     
      Wrk = MATMUL(TRANSPOSE(RR0),tempS)
      e1s = Wrk(1)

      Wrk = ZERO
      Wrk = MATMUL(TRANSPOSE(RR0),tempK)
      k1s = Wrk(1)
     
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

      Wrk33 = MATMUL(RR0,TRANSPOSE(RR0))
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

      Pe = 0.0d0
      Wrk(:) = fff(1:3)
      Pe(4:6,1:3) = Tilde(Wrk) + TRANSPOSE(epsi)
      Pe(4:6,4:6) = TRANSPOSE(mu)

      Qe = 0.0d0
      Wrk33(1:3,1:3) = Oe(1:3,4:6)
      Qe(4:6,4:6) = MATMUL(TRANSPOSE(Tilde(E1)),Wrk33)

      END SUBROUTINE ElasticForce 
