   SUBROUTINE ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)

   REAL(ReKi),INTENT(IN)::E1(:),RR0(:,:),kapa(:)
   REAL(ReKi),INTENT(IN)::Stif(:,:),cet
   REAL(ReKi),INTENT(OUT)::Fc(:),Fd(:)    
   REAL(ReKi),INTENT(OUT)::Oe(:,:),Pe(:,:),Qe(:,:) 
 
   REAL(ReKi)::eee(6),fff(6)
   REAL(ReKi)::tempS(3),tempK(3)
   REAL(ReKi)::Wrk(3),e1s,k1s,Wrk33(3,3)
   REAL(ReKi)::C11(3,3),C12(3,3),C21(3,3),C22(3,3)
   REAL(ReKi)::epsi(3,3),mu(3,3)

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