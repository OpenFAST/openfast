   SUBROUTINE ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd)
!----------------------------------------------------------------------------------------
! This subroutine computes Elastic Forces Fc and Fd.  
!----------------------------------------------------------------------------------------
   
   REAL(ReKi),INTENT(IN)::E1(:) ! E1 = x_0^\prime + u^\prime at Gauss point
   REAL(ReKi),INTENT(IN)::RR0(:,:) ! Rotation tensor at Gauss point
   REAL(ReKi),INTENT(IN)::kapa(:) ! Curvature starin vector at Gauss point
   REAL(ReKi),INTENT(IN)::Stif(:,:) ! C/S stiffness matrix resolved in inertial frame at Gauss point
   REAL(ReKi),INTENT(IN)::cet ! Extension-torsion coefficient at Gauss point
   REAL(ReKi),INTENT(OUT)::Fc(:) ! Elastic Force (N,M)
   REAL(ReKi),INTENT(OUT)::Fd(:) ! Elastic Force (E1,N)
 
   REAL(ReKi)::eee(6),fff(6)
   REAL(ReKi)::tempS(3),tempK(3)
   REAL(ReKi)::Wrk(3),e1s,k1s,Wrk33(3,3)

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
   Fd(4:6) = MATMUL(TRANSPOSE(Tilde(E1)),Wrk)

   END SUBROUTINE ElasticForce
