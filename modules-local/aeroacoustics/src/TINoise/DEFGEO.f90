!-----subroutine---------------------------------------------g.guidati--
! 
subroutine DEFGEO  
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
! 
!     scope             read airfoil and mesh coordinates
! 
!.......................................................................
!               declarations
!.......................................................................
      USE TINoiseGeneric
      USE TINoiseInput
      USE TINoiseGeo
      USE TIPrecision   
      USE TICoords                                              
                                                                      
IMPLICIT                        NONE                                  
                                                                     
   ! Local variables.                                                 
                                                                      
INTEGER(4)                   :: klo,khi,i                                                                    
REAL(DbKi)                   :: n1,n2,nozzw,y1dum,y2dum,scp,d1y1,d1y2
REAL(DbKi)                   :: y1in(m_in),y2in(m_in),ssin(m_in),y1,y2
REAL(DbKi)                   :: d2y1in(m_in),d2y2in(m_in),y1wh,y2wh,ddss

!.......................................................................
!     Spline input coordinates
!.......................................................................

!-----read coordinates

      y1in = x_ti
      y2in = y_ti

!-----rotate coordinates
      do i=1,n_in
        y1dum =   (y1in(i)-0.25d0)*cos(alfa)+y2in(i)*sin(alfa) + 0.25d0
        y2dum = - (y1in(i)-0.25d0)*sin(alfa)+y2in(i)*cos(alfa)
	y1in(i) = y1dum
	y2in(i) = y2dum
        ssin(i) = dble(i-1)
      enddo

      call SPL_P (ssin, y1in, n_in, d2y1in)
      call SPL_P (ssin, y2in, n_in, d2y2in)       

      
!.......................................................................
!     Re-panel the airfoil contour for flow
!....................................................................... 
      do i=1, n+1
        swork(i) = dble(i-1)
      enddo

      do i=1, n
        ds(i) = swork(i+1) - swork(i)
      enddo
      

      do i=1,n+1

	scp = swork(i) * ssin(n_in) / swork(n+1)
	call SPL_EX  (ssin,n_in,scp,khi,klo)
	call SPL_EX1 (ssin,y1in,d2y1in,n_in,scp,yc1(i),d1y1,khi,klo)
	call SPL_EX1 (ssin,y2in,d2y2in,n_in,scp,yc2(i),d1y2,khi,klo)
        
        nc1(i) =  d1y2 / sqrt(d1y1**2+d1y2**2)
        nc2(i) = -d1y1 / sqrt(d1y1**2+d1y2**2)
        tc1(i) =  d1y1 / sqrt(d1y1**2+d1y2**2)
        tc2(i) =  d1y2 / sqrt(d1y1**2+d1y2**2)
!        write(20,*) real(yc1(i)),real(yc2(i)),real(nc1(i)),real(nc2(i))
      enddo

      call SPL_P (swork, yc1, n+1, d2yc1)
      call SPL_P (swork, yc2, n+1, d2yc2)       
      
!      close(20)
             

!.......................................................................
!     Determine wake point at 'infinity'
!.......................................................................
      
      y1wh = 0.5d0 * (yc1(2) + yc1(n))
      y2wh = 0.5d0 * (yc2(2) + yc2(n))
        
      ddss = sqrt((yc1(1)-y1wh)**2+(yc2(1)-y2wh)**2)
      ywinf1 = yc1(1) + (yc1(1) - y1wh) /ddss * 100.
      ywinf2 = yc2(1) + (yc2(1) - y2wh) /ddss * 100.
       
      ywn1 =  (yc2(1)-y2wh)/ddss
      ywn2 = -(yc1(1)-y1wh)/ddss
      
        
!.......................................................................
!     Re-panel the airfoil contour for acoustics
!....................................................................... 
      
      do i=1, na+1
	sworkat(i) = dble(i-1)
      enddo

      do i=1,na+1
        scp = sworkat(i) * dble(n) / dble(na)
        call SPL_EX  (swork,n+1,scp,khi,klo)
        call SPL_EX1 (swork,yc1,d2yc1,n+1,scp,y1,d1y1,khi,klo)
        call SPL_EX1 (swork,yc2,d2yc2,n+1,scp,y2,d1y2,khi,klo)
        yc1at(i) =  y1 
        yc2at(i) =  y2
      enddo
      
      call SPL_PP (sworkat, yc1at, na+1, 1.0d33,1.0d33,d2yc1at)
      call SPL_PP (sworkat, yc2at, na+1, 1.0d33,1.0d33,d2yc2at) 

      do i=1,na+1
        scp = sworkat(i) * dble(n) / dble(na)
        call SPL_EX  (sworkat,na+1,scp,khi,klo)
        call SPL_EX1 (sworkat,yc1at,d2yc1at,na+1,scp,y1,d1y1,khi,klo)
        call SPL_EX1 (sworkat,yc2at,d2yc2at,na+1,scp,y2,d1y2,khi,klo)
        nc1at(i) =  d1y2 / sqrt(d1y1**2+d1y2**2)
        nc2at(i) = -d1y1 / sqrt(d1y1**2+d1y2**2)
        tc1at(i) = d1y1 / sqrt(d1y1**2+d1y2**2)
        tc2at(i) = d1y2 / sqrt(d1y1**2+d1y2**2)
      enddo

!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 

