   SUBROUTINE CrvMatrixR(cc,Rr) 

   REAL(ReKi),INTENT(IN)::cc(:)
   REAL(ReKi),INTENT(OUT)::Rr(:,:)

   INTEGER(IntKi):: i, j
   REAL(ReKi):: c1,c2,c3,c0,tr0

   Rr = 0.0D0
      
   c1 = cc(1)/4.0D0
   c2 = cc(2)/4.0D0
   c3 = cc(3)/4.0D0

   c0 = 0.5D0*(1.0D0-c1*c1-c2*c2-c3*c3)
   tr0 = 1.0D0 - c0
   tr0 = 2.0D0/(tr0*tr0)

   Rr(1,1) = tr0*(c1*c1 + c0*c0) - 1.0D0 
   Rr(2,1) = tr0*(c1*c2 + c0*c3)
   Rr(3,1) = tr0*(c1*c3 - c0*c2)
   Rr(1,2) = tr0*(c1*c2 - c0*c3)
   Rr(2,2) = tr0*(c2*c2 + c0*c0) - 1.0D0 
   Rr(3,2) = tr0*(c2*c3 + c0*c1)
   Rr(1,3) = tr0*(c1*c3 + c0*c2)
   Rr(2,3) = tr0*(c2*c3 - c0*c1)
   Rr(3,3) = tr0*(c3*c3 + c0*c0) - 1.0D0

   END SUBROUTINE CrvMatrixR
