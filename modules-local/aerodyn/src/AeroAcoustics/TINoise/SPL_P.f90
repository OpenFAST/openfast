!----------------------------------------------------j.ostertag---------
     subroutine SPL_P(x,y,n,y2)
!     
!     jasmin ostertag                       IAG
!                                           Universitaet Stuttgart
!     03.02.1998
!-----------------------------------------------------------------------
!
!     scope                     prepare spline interpolation
!                               (natural splines)
!-----------------------------------------------------------------------
!     declarations
!-----------------------------------------------------------------------
USE TIPrecision                                      
                                                     
IMPLICIT                        NONE                 
                                                     
                                                     
   ! Local variables.                                
                                                     
INTEGER(4)                   ::  n, i, k  
INTEGER(4),PARAMETER          ::  NMAX=5000        
                                                     
REAL(DbKi)                   :: yp1, ypn,p, qn, sig, un
REAL(DbKi)                   :: x(n), y(n), y2(n),u(NMAX)

               


      yp1 = -25.0d0 / 12.0d0 * y(1   ) &
            +48.0d0 / 12.0d0 * y(2   ) &
            -36.0d0 / 12.0d0 * y(3   ) &
            +16.0d0 / 12.0d0 * y(4   ) &
             -3.0d0 / 12.0d0 * y(5   )
      ypn =   3.0d0 / 12.0d0 * y(n-4) &
            -16.0d0 / 12.0d0 * y(n-3) &
            +36.0d0 / 12.0d0 * y(n-2) &
            -48.0d0 / 12.0d0 * y(n-1) &
           +25.0d0 / 12.0d0 * y(n  )

            
!-----------------------------------------------------------------------
!     prepare spline interpolation
!-----------------------------------------------------------------------
!     lower boundary condition is set either to be "natural"
!jo   erste Ableitung in Punkt 1      
      if (yp1.gt..99e30) then
        y2(1) = 0.
        u(1) = 0.
        
!     or else to have a first derivative        
      else
        y2(1) = -0.5
        u(1) = (3./(x(2) - x(1)))*((y(2) - y(1))/(x(2) - x(1)) - yp1)
      endif
      
!     this is the decomposition loop of the tridiagonal algorithm.
!     (y2 und u are used for temporary storage of the decomposed factors)
      do i=2,n-1
        sig = (x(i) - x(i-1))/(x(i+1) - x(i-1))
        p = sig*y2(i-1) + 2.0d0
        y2(i) = (sig - 1.0d0)/p
        u(i) = (6.0d0*((y(i+1) - y(i))/(x(i+1) - x(i)) - (y(i) - y(i-1))/ &
               (x(i) - x(i-1)))/(x(i+1) - x(i-1)) - sig*u(i-1))/p
      enddo
 
!     upper boundary condition is set either to be "natural"
!     erste Ableitung in Punkt n      
      if (ypn.gt..99e30) then
        qn = 0.
        un = 0.
        
!     or else to have a first derivative        
      else
        qn = 0.5
        un = (3.0d0/(x(n) - x(n-1)))*(ypn - (y(n) - y(n-1))/(x(n) - x(n-1)))  
      endif
     
      y2(n) = (un - qn*u(n-1))/(qn*y2(n-1) + 1.0d0)
          
!     backsubstitution loop of the tridiagonal algorithm
      do k=n-1,1,-1
        y2(k) = y2(k)*y2(k+1) + u(k)
      enddo
      
!-----------------------------------------------------------------------
!     end of subroutine
!-----------------------------------------------------------------------
      return
!-----------------------------------------------------------------------
      end
