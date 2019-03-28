      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
      USE TIPrecision
      IMPLICIT                        NONE                                   
                                                                       
   ! Local variables. 
INTEGER (4)                  :: n,i
INTEGER (4), PARAMETER       :: NMAX=50
   
REAL(DbKi)                   :: h,x,dydx(n),y(n),yerr(n),yout(n)
REAL(DbKi)                   :: ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX)
REAL(DbKi)                   :: ytemp(NMAX)
REAL(DbKi),PARAMETER         :: A2=.2
REAL(DbKi),PARAMETER         :: A3=.3
REAL(DbKi),PARAMETER         :: A4=.6 
REAL(DbKi),PARAMETER         :: A5=1.
REAL(DbKi),PARAMETER         :: A6=.875
REAL(DbKi),PARAMETER         :: B21=.2
REAL(DbKi),PARAMETER         :: B31=3./40.
REAL(DbKi),PARAMETER         :: B32=9./40.
REAL(DbKi),PARAMETER         :: B41=.3
REAL(DbKi),PARAMETER         :: B42=-.9
REAL(DbKi),PARAMETER         :: B43=1.2
REAL(DbKi),PARAMETER         :: B51=-11./54.
REAL(DbKi),PARAMETER         :: B52=2.5
REAL(DbKi),PARAMETER         :: B53=-70./27.
REAL(DbKi),PARAMETER         :: B54=35./27.
REAL(DbKi),PARAMETER         :: B61=1631./55296.
REAL(DbKi),PARAMETER         :: B62=175./512.
REAL(DbKi),PARAMETER         :: B63=575./13824.
REAL(DbKi),PARAMETER         :: B64=44275./110592.
REAL(DbKi),PARAMETER         :: B65=253./4096.
REAL(DbKi),PARAMETER         :: C1=37./378.
REAL(DbKi),PARAMETER         :: C3=250./621.
REAL(DbKi),PARAMETER         :: C4=125./594.
REAL(DbKi),PARAMETER         :: C6=512./1771.
REAL(DbKi),PARAMETER         :: DC1=C1-2825./27648.
REAL(DbKi),PARAMETER         :: DC3=C3-18575./48384.
REAL(DbKi),PARAMETER         :: DC4=C4-13525./55296.
REAL(DbKi),PARAMETER         :: DC5=-277./14336.
REAL(DbKi),PARAMETER         :: DC6=C6-.25

EXTERNAL          :: derivs

!U    USES derivs

      do i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
      enddo
      call derivs(x+A2*h,ytemp,ak2)
      do i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
      enddo
      call derivs(x+A3*h,ytemp,ak3)
      do i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
      enddo
      call derivs(x+A4*h,ytemp,ak4)
      do i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
      enddo
      call derivs(x+A5*h,ytemp,ak5)
      do i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+B65*ak5(i))
      enddo
      call derivs(x+A6*h,ytemp,ak6)
      do i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
      enddo
      do  i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*ak6(i))
      enddo
      return
      END
