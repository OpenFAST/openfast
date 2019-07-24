      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
            USE TIPrecision
      IMPLICIT                        NONE                                   
                                                                       
   ! Local variables. 
INTEGER (4)                  :: n,i
INTEGER (4), PARAMETER       :: NMAX=50
   
REAL(DbKi)                   :: eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
REAL(DbKi)                   :: errmax,h,xnew,yerr(NMAX),ytemp(NMAX)
REAL(DbKi),PARAMETER         :: SAFETY=0.9
REAL(DbKi),PARAMETER         :: PGROW=-.2
REAL(DbKi),PARAMETER         :: PSHRNK=-.25
REAL(DbKi),PARAMETER         :: ERRCON=1.89e-4

EXTERNAL          :: derivs
      
!U    USES derivs,rkck

      h=htry
1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.
      do i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
      enddo
      errmax=errmax/eps
      if(errmax.gt.1.)then
        h=SAFETY*h*(errmax**PSHRNK)
        if(h.lt.0.1*h)then
          h=.1*h
        endif
        xnew=x+h
        if(xnew.eq.x) then
           write(*,*)'ERROR:xfoil:rkqs: stepsize underflow in rkqs'
           STOP 1
         endif
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do i=1,n
          y(i)=ytemp(i)
        enddo
        return
      endif
      END
