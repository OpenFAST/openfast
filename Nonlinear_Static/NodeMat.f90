subroutine NodeMat(dmat, dloc,h0, h1, b0, b1, node_total, blength,&
                     &Young,G1)
    
   implicit double precision (a-h,o-z)
    
   integer node_total
    
   double precision dmat(node_total,3), dloc(node_total)
   double precision h0, h1, blength
   double precision b0, b1,Young,G1
  
   integer i,j
   double precision htemp, btemp
   
   dmat = 0.0d0
    
   do i = 1, node_total

      ! h at dloc(i)
      htemp = h0 * (1. - dloc(i) / blength) + dloc(i) *h1 / blength
      btemp = b0 * (1. - dloc(i) / blength) + dloc(i) *b1 / blength
          
      !Area of cross-section at dolc(i): btemp*htemp
      !Second moment of inertia: btemp * htemp**3 / 12

      dmat(i,1) = Young * btemp * htemp 
      dmat(i,2) = G1 * btemp * htemp
      dmat(i,3) = Young * btemp * htemp**3 / 12.   !Second moment of inertia 
      
       write(*,*) "dmat"
       do j=1,3
           write(*,*) dmat(1,j)
       enddo
!      stop
!       dmat(i,1) = 70.0d+09 
!       dmat(i,2) = 22.4359d+09 
!       dmat(i,3) = 5.833d+09 
       
!       dmat(i,1) = 1.0d4
!       dmat(i,2) = 1.0d4
!       dmat(i,3) = 1.0d2                            
   enddo
    
   return
end subroutine
