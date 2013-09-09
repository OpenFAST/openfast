   SUBROUTINE NodeLoc(dloc, xmin, elem_length,& 
                       &xj, norder, elem_total, node_total, blength)

   INTEGER(IntKi), INTENT(IN):: norder, node_total,elem_total
   REAL(ReKi), INTENT(OUT):: dloc(:), xj(:)   
 
   REAL(ReKi), INTENT(IN):: xmin, elem_length, blength
   
   INTEGER(IntKi):: i, j,npos, nelem
   INTEGER(IntKi):: temp
    
    
   dloc(1) = xmin
      
   npos = 1      
         
   do nelem = 1, elem_total 

      do j = 1, norder 

         npos = npos + 1

         dloc(npos) = (nelem - 1) * elem_length + elem_length*(xj(j+1)+1.d0)/2.d0               

      enddo
    
   enddo
      
   
   END SUBROUTINE NodeLoc
