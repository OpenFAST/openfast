!=======================================================================
!
      subroutine calll(x1,x2,y1,y2,dum,dum1,n1,n2,dipok,v1,v2,dv1d1,dv1d2,dv2d1,dv2d2)
      USE TIPrecision                                                
                                                                     
IMPLICIT                        NONE                                 
                                                                     
                                                                     
   ! Local variables.                                                
                                                                     
INTEGER(4)                   :: nvar, istr,nok,nbad                  
                                                                     
REAL(DbKi)                   :: n1,n2,small,pi,pi2,r1,r2,rsq,z,x,thet1
REAL(DbKi)                   :: dipok, v1_t,v2_t, dv1d1_t,dv1d2_t,dv2d1_t,dv2d2_t
REAL(DbKi)                   :: v1,v2,dv1d1,dv1d2,dv2d1,dv2d2,x1,x2,y1,y2,dum,dum1 
      

!                                        
      small = 1.0D-12                    
      pi = dacos(-1.0d0)                 
      pi2 = 2.0d0*dacos(-1.0d0)          
!      R1s = (x1-s1a)**2+(x2-s2a)**2     
!      R2s = (x1-s1b)**2+(x2-s2b)**2     
      r1 = (x1 - y1)                     
      r2 = (x2 - y2)                     
      rsq = r1**2+r2**2
      z = r1*n1+r2*n2
      x = r1*n2-r2*n1
      if(z.gt.small) then
        if(x.lt.small) then
          thet1 = datan(DBLE(z/x))+pi
        elseif(abs(x).le.small) then
          thet1 = 0.5d0*pi
        elseif(x.gt.small) then 
          thet1 = datan(DBLE(z/x))
        endif
      elseif(abs(z).le.small) then
        if(x.lt.small) then
          thet1 = 0.0d0
        elseif(abs(x).le.small) then
          write(*,*) 'SHITA!'
          stop
        elseif(x.gt.small) then 
!          write(*,*) 'SHITB!'
!          write(*,*) x1,x2,y1,y2
          thet1 = 0.0d0
!          stop
        endif
      else
        if(x.lt.small) then
          thet1 = datan(DBLE(z/x))-pi
        elseif(abs(x).le.small) then
          thet1 = -0.5d0*pi
        elseif(x.gt.small) then 
          thet1 = datan(DBLE(z/x))
        endif
      endif
!
      dipok = -thet1/pi2
      v1_t =  z/pi2*(1.0d0/rsq)
      v2_t =  -1.0d0/pi2*(x/rsq)
      dv1d1_t = -z*x/(pi*rsq**2)
      dv1d2_t = (rsq-2.0d0*z**2)/(2.0d0*pi*rsq**2)
      dv2d1_t = dv1d2_t
      dv2d2_t = -dv1d1_t
      v1 = +n1*v2_t+n2*v1_t
      v2 =  n2*v2_t-n1*v1_t
      dv1d1 = n2**2*dv1d1_t+n1*n2*(dv1d2_t+dv2d1_t)+n1**2*dv2d2_t
      dv1d2 = n2**2*dv1d2_t-n1*n2*(dv1d1_t-dv2d2_t)-n1**2*dv2d1_t
      dv2d1 = n2**2*dv2d1_t-n1*n2*(dv1d1_t-dv2d2_t)-n1**2*dv1d2_t
      dv2d2 = n2**2*dv2d2_t-n1*n2*(dv1d2_t+dv2d1_t)+n1**2*dv1d1_t



      return
      end
