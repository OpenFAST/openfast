subroutine NodeLoc(dloc, nquar, xmin, elem_length,& 
                       &xj, norder, elem_total, node_total, blength)

   implicit double precision (a-h,o-z)
    
   integer nquar, norder, node_total,elem_total
    
   double precision dloc(node_total), xj(norder+1)
   double precision xmin, elem_length, blength
   
   integer i, npos, nelem
   double precision temp
    
    
   dloc(1) = xmin
      
   npos = 1      
         
   do nelem = 1, elem_total 

      do j = 1, norder 

         npos = npos + 1

         dloc(npos) = (nelem - 1) * elem_length + elem_length*(xj(j+1)+1.d0)/2.d0               

      enddo
    
   enddo
      
   nquar = 1  

   temp = blength*2.0d0/4.0d0
   do i = 1, nmax 
      if ( abs(dloc(i)-temp) .lt. abs(dloc(nquar) -temp)) then
         nquar = i
      endif
   enddo
    
   return
end subroutine
