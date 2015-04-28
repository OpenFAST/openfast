   SUBROUTINE BD_ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)

   REAL(ReKi),INTENT(IN)::E1(:),RR0(:,:),kapa(:)
   REAL(ReKi),INTENT(IN)::Stif(:,:),cet
   REAL(ReKi),INTENT(OUT)::Fc(:),Fd(:)    
   REAL(ReKi),INTENT(OUT)::Oe(:,:),Pe(:,:),Qe(:,:) 
 
   REAL(ReKi)::eee(6),fff(6)
   REAL(ReKi)::tempS(3),tempK(3)
   REAL(ReKi)::Wrk(3),e1s,k1s,Wrk33(3,3)
   REAL(ReKi)::C11(3,3),C12(3,3),C21(3,3),C22(3,3)
   REAL(ReKi)::epsi(3,3),mu(3,3)

   INTEGER(IntKi):: i,j

   eee = 0.0D0 
   DO i=1,3
       eee(i) = E1(i) - RR0(i,1)
       eee(i+3) = kapa(i)

       tempS(i) = eee(i)
       tempK(i) = eee(i+3)
   ENDDO
   fff = 0.0D0 
   fff = MATMUL(Stif,eee)

   Wrk = 0.0D0     
   Wrk = MATMUL(TRANSPOSE(RR0),tempS)
   e1s = Wrk(1)      !epsilon_{11} in material basis

   Wrk = 0.0D0
   Wrk = MATMUL(TRANSPOSE(RR0),tempK)
   k1s = Wrk(1)      !kapa_{1} in material basis
     
   DO i=1,3
       fff(i) = fff(i) + 0.5D0*cet*k1s*k1s*RR0(i,1)
       fff(i+3) = fff(i+3) + cet*e1s*k1s*RR0(i,1)
   ENDDO 

   Fc = 0.0D0
   Fc = fff
   Wrk = 0.0D0 
   Wrk(1:3) = fff(1:3)
   Fd = 0.0D0 
   Fd(4:6) = MATMUL(TRANSPOSE(BD_Tilde(E1)),Wrk)

   C11 = 0.0D0
   C12 = 0.0D0
   C21 = 0.0D0
   C22 = 0.0D0
   C11(1:3,1:3) = Stif(1:3,1:3)
   C12(1:3,1:3) = Stif(1:3,4:6)
   C21(1:3,1:3) = Stif(4:6,1:3)
   C22(1:3,1:3) = Stif(4:6,4:6)

   Wrk = 0.0D0
   DO i=1,3
       Wrk(i) = RR0(i,1)
   ENDDO
   Wrk33 = 0.0D0
   Wrk33 = BD_OuterProduct(Wrk,Wrk) 
   C12 = C12 + cet*k1s*Wrk33
   C21 = C21 + cet*k1s*Wrk33
   C22 = C22 + cet*e1s*Wrk33

   epsi = 0.0D0 
   mu = 0.0D0
   epsi = MATMUL(C11,BD_Tilde(E1))
   mu = MATMUL(C21,BD_Tilde(E1))
   
   Wrk = 0.0D0

   Oe = 0.0D0
   Oe(1:3,4:6) = epsi(1:3,1:3)
   Oe(4:6,4:6) = mu(1:3,1:3)
   
   Wrk(1:3) = fff(1:3)
   Oe(1:3,4:6) = Oe(1:3,4:6) - BD_Tilde(Wrk)
   Wrk = 0.0D0
   Wrk(1:3) = fff(4:6)
   Oe(4:6,4:6) = Oe(4:6,4:6) - BD_Tilde(Wrk)

   Pe = 0.0D0
   Wrk = 0.0D0
   Wrk(1:3) = fff(1:3)
   Pe(4:6,1:3) = BD_Tilde(Wrk) + TRANSPOSE(epsi)
   Pe(4:6,4:6) = TRANSPOSE(mu)

   Qe = 0.0D0
   Wrk33 = 0.0D0
   Wrk33(1:3,1:3) = Oe(1:3,4:6)
   Qe(4:6,4:6) = MATMUL(TRANSPOSE(BD_Tilde(E1)),Wrk33)

   END SUBROUTINE BD_ElasticForce
