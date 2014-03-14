   SUBROUTINE CrvMatrixHD(cc,cd,Hh,Hd)

   REAL(ReKi),INTENT(IN):: cc(:), cd(:), Hh(:,:)
   
   REAL(ReKi),INTENT(OUT):: Hd(:,:)

   REAL(ReKi):: cf1,cf2,cf3,cq,ocq,aa,cv0,cv1,cv2,cv3,ada

   cf1 = cc(1)/4.0D0
   cf2 = cc(2)/4.0D0
   cf3 = cc(3)/4.0D0

   cq = cf1*cf1 + cf2*cf2 + cf3*cf3
   ocq = 1.0D0 + cq
   aa = 2.0D0*ocq*ocq

   cv1 = cd(1)/aa
   cv2 = cd(2)/aa
   cv3 = cd(3)/aa
   cv0 = -(cf1*cv1 + cf2*cv2 + cf3*cv3)
   
   ada = -2.0D0*ocq*cv0

   Hd = 0.0D0

   Hd(1,1) = cv1*cf1 + cf1*cv1 + cv0 - ada*Hh(1,1)
   Hd(2,1) = cv2*cf1 + cf2*cv1 + cv3 - ada*Hh(2,1)
   Hd(3,1) = cv3*cf1 + cf3*cv1 - cv2 - ada*Hh(3,1)
   Hd(1,2) = cv1*cf2 + cf1*cv2 - cv3 - ada*Hh(1,2)
   Hd(2,2) = cv2*cf2 + cf2*cv2 + cv0 - ada*Hh(2,2)
   Hd(3,2) = cv3*cf2 + cf3*cv2 + cv1 - ada*Hh(3,2)
   Hd(1,3) = cv1*cf3 + cf1*cv3 + cv2 - ada*Hh(1,3)
   Hd(2,3) = cv2*cf3 + cf2*cv3 - cv1 - ada*Hh(2,3)
   Hd(3,3) = cv3*cf3 + cf3*cv3 + cv0 - ada*Hh(3,3)

   

   END SUBROUTINE CrvMatrixHD
