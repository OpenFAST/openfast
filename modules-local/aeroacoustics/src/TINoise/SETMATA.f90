!-----subroutine---------------------------------------------g.guidati--
! 
     subroutine SETMATA 
! 
!-----------------------------------------------------------------------
!               
!     franco guidati    IAG
!                       Uni Stuttgart
!
! 
!     scope             set up acoustic influence 
!                       coefficients matrix  
! 
!.......................................................................
!               declarations
!.......................................................................
      USE TINoiseGeneric
      USE TINoiseGeo
      USE TINoiseDDD
      USE TINoiseCancela
      IMPLICIT                        NONE                                 
                                                                           
   ! Local variables.                                                      
INTEGER (4)                  :: i,j,ishift,jshift,ng,info                                 
                                                                           
REAL(DbKi)                   :: rhelp,solidangle,summ
COMPLEX(DbKi)                :: HERM(na+1,na+1,4),ampi,value
COMPLEX(DbKi)                :: herm1, herm2, herm3, herm4, HANK1, HANK0
!.......................................................................
!              
!.......................................................................
      ng = 10

!      write(*,*) 'SETTING UP AIC-MATRIX'
      rhelp = acos(-tc1at(1)*tc1at(na+1) - tc2at(1)*tc2at(na+1))
      solidangle = 1.0 - 0.5*(2.0*rhelp/pi2)
      
      do i=1,na+2
        do j=1,na+2
	  Kerna(i,j) = 0.0d0
	end do
      end do
                 
      do i=1,na
      
!          write(*,*) i


        x1 = yc1at(i)
        x2 = yc2at(i)

        Kerna(i,i) = 0.5d0
        if(i.eq.1) Kerna(i,i) = solidangle

        do j=1,na
          s1 = sworkat(j)   + 0.0000001
          s2 = sworkat(j+1) - 0.0000001
          
          call CDA0   (sworkat, yc1at, yc2at, d2yc1at, d2yc2at,na, ng, td, Ad, pi2i, kwave,herm1, herm2, herm3, herm4 )
     
          HERM(i,j,1) = herm1 
          HERM(i,j,2) = herm2 
          HERM(i,j,3) = herm3 
          HERM(i,j,4) = herm4
          
          
          do jshift=0,1

            Kerna(i,j+jshift) = Kerna(i,j+jshift) - HERM(i,j,jshift+1)

            ishift=0
            if(j+jshift.eq.1     ) ishift= 2
            if(j+jshift.eq.2     ) ishift= 1
            if(j+jshift.eq.n     ) ishift=-1
            if(j+jshift.eq.n+1   ) ishift=-2
  	    Kerna(i,j-2+ishift+jshift) = Kerna(i,j-2+ishift+jshift) -dst(-2,ishift)*HERM(i,j,jshift+3)*.5

            Kerna(i,j-1+ishift+jshift) = Kerna(i,j-1+ishift+jshift) -dst(-1,ishift)*HERM(i,j,jshift+3)*.5

            Kerna(i,j  +ishift+jshift) = Kerna(i,j  +ishift+jshift) -dst( 0,ishift)*HERM(i,j,jshift+3)*.5

            Kerna(i,j+1+ishift+jshift) = Kerna(i,j+1+ishift+jshift) -dst( 1,ishift)*HERM(i,j,jshift+3)*.5

	    Kerna(i,j+2+ishift+jshift) = Kerna(i,j+2+ishift+jshift) -dst( 2,ishift)*HERM(i,j,jshift+3)*.5
      
          enddo
                
        enddo
	
      enddo
      


      Kerna(na+1,1)       =  1.0
      Kerna(na+1,na+1)    = -1.0
      Kerna(na+2,na+2)    =  1.0


      

      call zgetrf (na+2,na+2,Kerna,na+2,ipiva,info)
      if(info.ne.0) then
        write(*,*) 'oopsy - something is wrong'
        stop
      endif
!      write(*,*) 'AKUSTIK MATRIX INVERTIERT',info


!.......................................................................
!               end of subroutine
!.......................................................................
      return
!***********************************************************************
      end 

