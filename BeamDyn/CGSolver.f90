   SUBROUTINE CGSolver(RHS,KT,ui,bc,dof_total)
   
   INTEGER(IntKi), INTENT(IN):: dof_total
   REAL(ReKi), INTENT(IN):: RHS(:), KT(:,:), bc(:)
   REAL(ReKi), INTENT(INOUT):: ui(:)
   
   INTEGER(IntKi):: lmax, i, j, l   
   REAL(ReKi):: KTu(dof_total)
   REAL(ReKi):: eps, alphatop, alphabot, alpha_cg, alphatop_new
   REAL(ReKi):: beta_cg, p(dof_total), r(dof_total)
    
    
   eps = 1.0D-05
   lmax = 100000000
    
   ui = 0.0d0  ! use zero as initial condition
     
   alphatop = 0.d0  ! mas: this was missing

   KTu = 0.0d0

   DO i = 1,dof_total
      r(i) = bc(i) * (RHS(i) - KTu(i))
      p(i) = r(i)
      alphatop = alphatop + r(i) * r(i)
   ENDDO
    
   DO l = 2,lmax 
      KTu=0.0d0
      DO i=1,dof_total
         DO j=1,dof_total
            KTu(i) = KTu(i) + KT(i,j)*p(j)
         ENDDO
      ENDDO
        
      alphabot = 0.0d0
      DO i = 1,dof_total
         alphabot = alphabot + p(i) * KTu(i)
      ENDDO
        
      alpha_cg = alphatop / alphabot
        
        
      DO i=1, dof_total
         ui(i) = bc(i) * (ui(i) + alpha_cg * p(i) )
         r(i) = bc(i) * (r(i) - alpha_cg * KTu(i))
      ENDDO
         
      alphatop_new = 0.d0
      DO i=1,dof_total
         alphatop_new = alphatop_new + r(i) * r(i)
      ENDDO
        
      IF(SQRT(alphatop_new) .LE. eps) THEN
!         goto 20
!          WRITE(*,*) "eps =",SQRT(alphatop_new)
!          WRITE(*,*) "CG Iterations",l
          RETURN
      ENDIF
          
!         IF(sqrt(alphatop_new) .GT. sqrt(alphatop)) THEN
!             WRITE(*,*) 'cg iteration diverging'
!             WRITE(*,*) 'rsold = ', sqrt(alphatop)
!             WRITE(*,*) 'rsnew = ', sqrt(alphatop_new)
!             WRITE(*,*) 'its = ', l
!             STOP
!         ENDIF
        
      beta_cg = alphatop_new / alphatop
         
      DO i=1, dof_total
         p(i) = bc(i) * (r(i) + beta_cg * p(i))
      ENDDO
          
      alphatop = alphatop_new
          
      IF(l==lmax) THEN
         WRITE(*,*) 'CG failed to converge in lmax = ', lmax
         WRITE(*,*) 'eps=',SQRT(alphatop_new)
         WRITE(*,*) 'STOPPING'
         STOP
      ENDIF

   ENDDO
         
   20 WRITE(*,*) 'CG iterations', l
   RETURN 
   END SUBROUTINE CGSolver
