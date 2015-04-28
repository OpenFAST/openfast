   SUBROUTINE BD_CrvMatrixH(cc,Hh) 

   REAL(ReKi),INTENT(IN)::cc(:)
   REAL(ReKi),INTENT(OUT)::Hh(:,:)

   INTEGER(IntKi):: i, j
   REAL(ReKi):: cf1,cf2,cf3,cq,ocq,aa,cb0,cb1,cb2,cb3
   
   cf1 = cc(1)/4.0D0
   cf2 = cc(2)/4.0D0
   cf3 = cc(3)/4.0D0
   cq = cf1 * cf1 + cf2 * cf2 + cf3 * cf3
   ocq = 1.0D0 + cq
   aa = 2.0D0 * ocq * ocq
   cb0 = 2.0D0 * (1.0D0 - cq) / aa
   cb1 = cc(1)/aa
   cb2 = cc(2)/aa
   cb3 = cc(3)/aa
   
   Hh = 0.0D0
   
   Hh(1,1) = cb1 * cf1 + cb0
   Hh(2,1) = cb2 * cf1 + cb3
   Hh(3,1) = cb3 * cf1 - cb2
   Hh(1,2) = cb1 * cf2 - cb3
   Hh(2,2) = cb2 * cf2 + cb0
   Hh(3,2) = cb3 * cf2 + cb1
   Hh(1,3) = cb1 * cf3 + cb2
   Hh(2,3) = cb2 * cf3 - cb1
   Hh(3,3) = cb3 * cf3 + cb0

   END SUBROUTINE BD_CrvMatrixH
