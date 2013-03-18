subroutine CGSolver(RHS,KT,ui,bc,dof_total)
   
   integer dof_total
   double precision RHS(dof_total), KT(dof_total,dof_total)
   double precision ui(dof_total), bc(dof_total)
   
   integer lmax, i, j, l   
   double precision KTu(dof_total)
   double precision eps, alphatop, alphabot, alpha_cg, alphatop_new
   double precision beta_cg, p(dof_total), r(dof_total)
    
    
   eps = 1.0d-14
   lmax = 1000000
    
   ui = 0.0d0  ! use zero as initial condition
     
   alphatop = 0.d0  ! mas: this was missing

   KTu = 0.0d0

!  do i=1,dof_total
!     write(*,*) "Line=",i
!     do j=1,dof_total
         
!        write(*,*) KT(i,j)
!        KTu(i) = KTu(i) + KT(i,j)*ui(j)
!     enddo
!  enddo
    
!  do i=1,dof_total
!     write(*,*) "RHS Line=",i
!     write(*,*) RHS(i)
!  enddo
    
   do i = 1,dof_total
      r(i) = bc(i) * (RHS(i) - KTu(i))
      p(i) = r(i)
      alphatop = alphatop + r(i) * r(i)
   enddo
    
   do l = 2,lmax
    
      KTu=0.0d0
      do i=1,dof_total
         do j=1,dof_total
            KTu(i) = KTu(i) + KT(i,j)*p(j)
         enddo
      enddo
        
      alphabot = 0.0d0
      do i = 1,dof_total
         alphabot = alphabot + p(i) * KTu(i)
      enddo
        
      alpha_cg = alphatop / alphabot
        
        
      do i=1, dof_total
         ui(i) = bc(i) * (ui(i) + alpha_cg * p(i) )
         r(i) = bc(i) * (r(i) - alpha_cg * KTu(i))
      enddo
         
      alphatop_new = 0.d0
      do i=1,dof_total
         alphatop_new = alphatop_new + r(i) * r(i)
      enddo
        
      if(sqrt(alphatop_new) .LE. eps) then
         goto 20
      endif
          
!         if(sqrt(alphatop_new) .GT. sqrt(alphatop)) then
!             write(*,*) 'cg iteration diverging'
!             write(*,*) 'rsold = ', sqrt(alphatop)
!             write(*,*) 'rsnew = ', sqrt(alphatop_new)
!             write(*,*) 'its = ', l
!             stop
!         endif
        
      beta_cg = alphatop_new / alphatop
         
      do i=1, dof_total
         p(i) = bc(i) * (r(i) + beta_cg * p(i))
      enddo
          
      alphatop = alphatop_new
          
      if(l==lmax) then
         write(*,*) 'CG failed to converge in lmax = ', lmax
         write(*,*) 'eps=',sqrt(alphatop_new)
         write(*,*) 'stopping'
         stop
      endif

   enddo
         
20 write(*,*) 'CG iterations', l
 
   return
end subroutine
