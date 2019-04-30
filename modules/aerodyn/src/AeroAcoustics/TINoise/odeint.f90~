      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,rkqs)
      USE TINoisePATH
      USE TIPrecision 
!ODE solver "driver", using Runge-Kutta with adaptive stepsize control 
!From Numerical Recipes                                               
                                                                      
IMPLICIT                        NONE                                  
                                                                     
   ! Local variables.                                                 
                                                                      
INTEGER(4)                   :: nbad,nok,nvar                           
INTEGER(4)                   :: i,nstp                                                                     
REAL(DbKi)                   :: eps,h1,hmin,x1,x2,ystart(nvar)
REAL(DbKi)                   :: h,hdid,hnext,x,xsav,dydx(NMAX), y(NMAX),yscal(NMAX)

EXTERNAL         :: derivs,rkqs

      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) then
           pause 'stepsize smaller than minimum in odeint'
	     endif
        h=hnext
16    continue
      pause 'too many steps in odeint'
      return
      END
