SUBROUTINE lubksb(a,n,indx,b,ui)
   
   integer n,indx(n)
   double precision a(n,n),b(n),ui(n)
    
   integer i,ii,j,ll
   double precision sum
            
   ii = 0
   DO i = 1 , n
      ll = indx(i)
      sum = b(ll)
      b(ll) = b(i)
      IF (ii .ne. 0) THEN
         DO j = ii , i-1
            sum = sum - a(i,j) * b(j)
         ENDDO
      ELSE IF (sum .ne. 0.0) THEN
         ii = i
      ENDIF
      b(i) = sum
   ENDDO
      
   DO i = n,1,-1
      sum = b(i)
      DO j = i+1 , n
         sum = sum - a(i,j) * b(j)
      ENDDO
      b(i) = sum / a(i,i)
   ENDDO
    
   do i=1,n
      ui(i)=b(i)
   enddo
    
   RETURN
   
END subroutine
