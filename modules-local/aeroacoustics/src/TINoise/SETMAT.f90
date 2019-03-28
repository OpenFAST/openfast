!-----subroutine---------------------------------------------g.guidati--
! 
     subroutine SETMAT  
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
! 
!     scope             set up aerodynamic influence 
!                       coefficients matrix  
! 
!.......................................................................
!               declarations
!.......................................................................
      USE TINoiseGeneric
      USE TINoiseDDD
      USE TINoiseGeo
      IMPLICIT                        NONE                                 
                                                                           
   ! Local variables.                                                      
INTEGER (4)                  :: i,j,ishift,jshift,ng,info                                 
                                                                           
REAL(DbKi)                   :: rhelp,solidangle,summ,herm1,herm2,herm3,herm4
REAL(DbKi)                   :: HERM(n+1,n+1,4),dipok,dv1,dv2,dv1d1,dv1d2
REAL(DbKi)                   :: dv2d1,dv2d2

      
!.......................................................................
!              
!.......................................................................
      ng = 4
      rhelp = acos(-tc1(1)*tc1(n+1) - tc2(1)*tc2(n+1))
      solidangle = (1.0 - 0.5*(2.0*rhelp/pi2)) 
      
      do i=1,n+2
        do j=1,n+2
	  Kern(i,j) = 0.0
	end do
      end do
      
      do i=1,N
!        write(*,*) i
        x1 = yc1(i)
        x2 = yc2(i)
!        write(55,*) real(x1),real(x2),real(nc1(i)),real(nc2(i)),
!     >              real(tc1(i)),real(tc2(i))   
        
        Kern(i,i) = 0.5d0
        if(i.eq.1) Kern(i,i) = solidangle

        summ = 0.0d0
        do j=1,N
          s1 = swork(j)   + 0.0000001d0
          s2 = swork(j+1) - 0.0000001d0
          call CDI0   (swork, yc1, yc2, d2yc1, d2yc2,n, ng, td, Ad, pi2i,herm1, herm2, herm3, herm4 )

          HERM(i,j,1) = herm1
          HERM(i,j,2) = herm2
          HERM(i,j,3) = herm3
          HERM(i,j,4) = herm4
          
          summ = summ + HERM(i,j,1) + HERM(i,j,2)
          
!          write(*,*) HERM(i,j,3),HERM(i,j,4)
          do jshift=0,1

            Kern(i,j+jshift) = Kern(i,j+jshift) - HERM(i,j,jshift+1)

            ishift=0
            if(j+jshift.eq.1     ) ishift= 2
            if(j+jshift.eq.2     ) ishift= 1
            if(j+jshift.eq.n     ) ishift=-1
            if(j+jshift.eq.n+1   ) ishift=-2
  	    Kern(i,j-2+ishift+jshift) = Kern(i,j-2+ishift+jshift) -dst(-2,ishift)*HERM(i,j,jshift+3)*.5

            Kern(i,j-1+ishift+jshift) = Kern(i,j-1+ishift+jshift) -dst(-1,ishift)*HERM(i,j,jshift+3)*.5

            Kern(i,j  +ishift+jshift) = Kern(i,j  +ishift+jshift) -dst( 0,ishift)*HERM(i,j,jshift+3)*.5

            Kern(i,j+1+ishift+jshift) = Kern(i,j+1+ishift+jshift) -dst( 1,ishift)*HERM(i,j,jshift+3)*.5

	    Kern(i,j+2+ishift+jshift) = Kern(i,j+2+ishift+jshift) -dst( 2,ishift)*HERM(i,j,jshift+3)*.5
      
          enddo
                
        enddo
!	if(i.eq.1) then
!  	  write(*,*) i,summ+1.0
!	  write(*,*) rhelp,solidangle
!	endif
        
        
                
        call calll(x1,x2,yc1(1),yc2(1),ywinf1,ywinf2,ywn1,ywn2,dipok,dv1,dv2,dv1d1,dv1d2,dv2d1,dv2d2)
                   
        Kern(i,n+2) = -dipok
        
      enddo

      if(.false.) then
      Kern(n+1,1)      = -25.0d0/12.0
      Kern(n+1,2)      =  48.0d0/12.0
      Kern(n+1,3)      = -36.0d0/12.0
      Kern(n+1,4)      =  16.0d0/12.0
      Kern(n+1,5)      =  -3.0d0/12.0
      Kern(n+1,n+1)    =  25.0d0/12.0
      Kern(n+1,n  )    = -48.0d0/12.0
      Kern(n+1,n-1)    =  36.0d0/12.0
      Kern(n+1,n-2)    = -16.0d0/12.0
      Kern(n+1,n-3)    =   3.0d0/12.0
      else
      Kern(n+1,1)      =  2.0
      Kern(n+1,2)      =  -3.0
      Kern(n+1,3)      =  1.0
      Kern(n+1,n+1)    =  -2.0
      Kern(n+1,n  )    =  3.0
      Kern(n+1,n-1)    =  -1.0
      endif
      
      
      Kern(n+2,1)   =  1.0
      Kern(n+2,n+1) = -1.0
      Kern(n+2,n+2) = -1.0
      Kern(1,n+2) = -solidangle/2.0

      call dgetrf (n+2,n+2,Kern,n+2,ipiv,info)

!.......................................................................
!               end of subroutine
!.......................................................................
     return
!***********************************************************************
      end 

