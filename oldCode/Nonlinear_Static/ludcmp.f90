SUBROUTINE ludcmp(a,n,indx,d) 

   integer n,indx(n)
   double precision d,a(n,n),tiny
   PARAMETER (tiny=1.0e-20) 

   integer i,imax,j,k
   double precision aamax,dum,sum,vv(n)
   
   d=1.d0 
   do i=1,n 
      aamax=0. 
      do j=1,n 
        if (abs(a(i,j)) .gt. aamax) aamax=abs(a(i,j)) 
      enddo 
      if (aamax==0.) pause 'singular matrix.' 
      vv(i)=1./aamax 
   enddo 

   do j=1,n 
      do i=1,j-1 
         sum=a(i,j) 
         do k=1,i-1 
            sum=sum-a(i,k)*a(k,j) 
         enddo 
         a(i,j)=sum 
      enddo 
      aamax=0. 
      do i=j,n 
         sum=a(i,j) 
         do k=1,j-1 
            sum=sum-a(i,k)*a(k,j) 
         enddo 
         a(i,j)=sum 
         dum=vv(i)*abs(sum) 
         if (dum.ge.aamax) then 
            imax=i 
            aamax=dum 
         endif 
      enddo 
      if (j.ne.imax) then 
         do k=1,n 
            dum=a(imax,k) 
            a(imax,k)=a(j,k) 
            a(j,k)=dum 
         enddo 
         d=-d 
         vv(imax)=vv(j) 
      endif 
      indx(j)=imax 
      if (a(j,j).eq.0.) a(j,j)=tiny 
      if (j.ne.n) then 
        dum=1./a(j,j) 
        do i=j+1,n 
            a(i,j)=a(i,j)*dum 
        enddo 
      endif 
   enddo 

   return 
END SUBROUTINE ludcmp 
