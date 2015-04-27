SUBROUTINE ludcmp(a,n,indx,d) 
!***************************************************************************************
! This subroutine performs LU Decomposition
!***************************************************************************************

INTEGER(IntKi),INTENT(IN):: n ! DOF total - 6
INTEGER(IntKi),INTENT(OUT):: indx(:) 
REAL(ReKi),INTENT(INOUT):: a(:,:) ! Mass matrix
REAL(ReKi),INTENT(OUT):: d ! Determinate

REAL(ReKi),PARAMETER:: tolf = 1.0D-20

INTEGER(IntKi):: i,imax,j,k
REAL(ReKi):: aamax,dum,summ,vv(n)
	
	d=1.d0 
do i=1,n 
    aamax=0. 
    do j=1,n 
        if (abs(a(i,j)) .gt. aamax) aamax=abs(a(i,j)) 
    enddo 
    if (aamax==0.) then
        write(*,*) "singular matrix"
        stop
    endif
    vv(i)=1./aamax 
enddo 

do j=1,n 
	do i=1,j-1 
    	summ=a(i,j) 
        do k=1,i-1 
        	summ=summ-a(i,k)*a(k,j) 
        enddo 
        a(i,j)=summ 
    enddo 
    aamax=0. 
    do i=j,n 
        summ=a(i,j) 
        do k=1,j-1 
        	summ=summ-a(i,k)*a(k,j) 
    	enddo 
        a(i,j)=summ 
        dum=vv(i)*abs(summ) 
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
    if (a(j,j).eq.0.) a(j,j)=tolf 
    if (j.ne.n) then 
        dum=1./a(j,j) 
        do i=j+1,n 
            a(i,j)=a(i,j)*dum 
        enddo 
    endif 
enddo 
 
END SUBROUTINE ludcmp 
