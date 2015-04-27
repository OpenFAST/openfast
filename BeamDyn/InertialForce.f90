   SUBROUTINE InertialForce(m00,mEta,rho,vvv,aaa,Fi,Mi,Gi,Ki)

   REAL(ReKi),INTENT(IN)::m00,mEta(:),rho(:,:)
   REAL(ReKi),INTENT(IN)::vvv(:),aaa(:)
   REAL(ReKi),INTENT(INOUT)::Fi(6),Mi(6,6),Gi(6,6),Ki(6,6)

   REAL(ReKi)::beta(3),gama(3),nu(3),epsi(3,3),mu(3,3)
   REAL(ReKi)::ome(3),omd(3),tempV(3),tempA(3)
   INTEGER(IntKi)::i

   ! Prepare
DO i=1,6
!WRITE(*,*) 'aaa:',i,aaa(i)
ENDDO
   ome(:) = vvv(4:6)
   omd(:) = aaa(4:6)
   tempV(:) = vvv(1:3)
   tempA(:) = aaa(1:3)

   beta = MATMUL(Tilde(ome),mEta)
   gama = MATMUL(rho,ome)
   nu = MATMUL(rho,omd)

   !Compute Fi
   Fi(1:3)= m00*tempA+MATMUL(Tilde(omd),mEta)+MATMUL(Tilde(ome),beta)
   Fi(4:6) = MATMUL(Tilde(mEta),tempA) + nu + MATMUL(Tilde(ome),gama) 

   !Mass Matrix
   Mi = 0.0D0
   DO i=1,3
       Mi(i,i) = m00
   ENDDO
   Mi(1:3,4:6) = TRANSPOSE(Tilde(mEta))
   Mi(4:6,1:3) = Tilde(mEta)
   Mi(4:6,4:6) = rho

   !Gyroscopic Matrix
   Gi = 0.0D0
   epsi = MATMUL(Tilde(ome),rho)
   mu = MATMUL(Tilde(ome),TRANSPOSE(Tilde(mEta)))
   Gi(1:3,4:6) = TRANSPOSE(Tilde(beta)) + mu
   Gi(4:6,4:6) = epsi - Tilde(gama)

   !Stiffness Matrix
   Ki = 0.0D0
   Ki(1:3,4:6) = MATMUL(Tilde(omd),TRANSPOSE(Tilde(mEta))) +&
                &MATMUL(Tilde(ome),mu)
   Ki(4:6,4:6) = MATMUL(Tilde(tempA),Tilde(mEta)) + &
                &MATMUL(rho,Tilde(omd)) - Tilde(nu) +&
                &MATMUL(epsi,Tilde(ome)) - &
                &MATMUL(Tilde(ome),Tilde(gama))
!   WRITE(*,*) "Ki at Node #"
!   DO i=1,6
!       WRITE(*,*) Ki(i,1),Ki(i,2),Ki(i,3),Ki(i,4),Ki(i,5),Ki(i,6)
!   ENDDO

   END SUBROUTINE InertialForce
