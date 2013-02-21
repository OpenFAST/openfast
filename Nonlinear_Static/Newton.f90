	subroutine NewtonRaphson(dof_node,dof_total,norder,node_total,elem_total,&
							&hhp,uf,dmat,wj,niter,Jacobian)
	
	integer dof_node,dof_total,norder,node_total,elem_total,niter
	double precision hhp(norder+1,norder+1),wj(norder+1),dmat(node_total,3)
	double precision uf(dof_total),Jacobian
	
	logical check
	double precision RHS(dof_total),KT(dof_total,dof_total),gradient(dof_total)
	double precision ui(dof_total),bc(dof_total)
	double precision fmin,fold,stpmax,sum,uold(dof_total)
	parameter(TOLF = 1.0d-12)
	parameter(STPMX=1000.0d0)
	integer i,j,m
	double precision temp1,temp2
	integer indx(dof_total)
	double precision d
	
	ui=0.0d0
	
	call AssembleRHS(RHS,dof_node,dof_total,uf,&
							&norder,hhp,wj,node_total,dmat,&
							&elem_total)
							
	if(MAXVAL(ABS(RHS))<TOLF) return
    
	fmin=0.5d0*DOT_PRODUCT(RHS,RHS)  ! calculate fnew=1/2*F.F
    
    temp1 = float(size(uf))
	
	temp2 = 0.0d0
	
	call Norm(dof_total,uf,temp2)
	
	stpmax=STPMX*MAX(temp2,temp1)   ! Calculate step length for line searches

	check=.FALSE.	


    
	do i=1,niter
   		write(*,*) "ITERATION=",i
   		! Assemble the nonlinear system
   		!---------------------------------------------------
   		call AssembleKT(KT,dof_node,dof_total,norder,node_total,elem_total,&
							&hhp,uf,dmat,wj,Jacobian)  

		
		do j=1,dof_total
			sum=0.0d0
			do m=1,dof_total
				sum=sum+KT(m,j)*(-RHS(m))
			enddo
			gradient(j)=sum
		enddo
		
		
!   gradient=MATMUL_sparse(-rhs,nsize,ne,irn,jcn,coef) ! compute gradient for the line search
   		uold=uf
   		fold=fmin

   ! Solve the linear system
   !====================================
   		bc=1.0d0
   		
   		bc(1) = 0.0d0
   		bc(2) = 0.0d0
   		bc(3) = 0.0d0
   
   		call CGSolver(RHS,KT,ui,bc,dof_total)
		
!		indx = 0
!		d=0.0d0
		
!		do j=1,dof_total
!			write(*,*) "KT Line=",j
!			do m=1,dof_total
!				write(*,*) KT(j,m)
!			enddo 
!		enddo
		
!		do j=1,dof_total
!			write(*,*) "RHS Line=",j
!			write(*,*) RHS(j)
!		enddo
!		stop
		
!		call ludcmp(KT,dof_total,indx,d)
		
!		do j=1,dof_total
!			write(*,*) "KT Line=",j
!			do m=1,dof_total
!				write(*,*) KT(j,m)
!			enddo 
!		enddo
!		stop
!		call lubksb(KT,dof_total,indx,RHS,ui)
		   		
!   		write(*,*) "ui"
!   		do j=1,dof_total
!   			write(*,*) ui(j)
!   		enddo
!   		stop
   		
   		   		
   		call LineSearch(dof_total,uold,fold,gradient,rhs,ui,fmin,stpmax,check,&
							&dof_node,uf,norder,hhp,wj,node_total,dmat,elem_total)
							
!		write(*,*) "test1",check
!		write(*,*) "TOLF",MAXVAL(ABS(RHS))<TOLF
							
   		IF(MAXVAL(ABS(RHS))<TOLF) THEN
       		check=.FALSE.  ! the function values converge to zero
	   		RETURN 
   		ENDIF 
   		
   		
   		IF(check) THEN  ! converges on deltaX
      	! Check for gradient of 1/2 F.F zero, i.e. spurious convergence
	  	!---------------------------------------------------------------------
	   		check=(MAXVAL(ABS(gradient)*MAX(ABS(uf),1.0D0)/MAX(fmin,0.5D0*SIZE(uf)))<TOLF) 
	   		
!	   		write(*,*) "test2",check
	   		EXIT
   		ENDIF   

   		IF(i==niter) write(*,*) "The solution does not converge after the maximum number of iterations" !The maximum number of iterations reached, the solution is still not converged. 
	enddo

	if(check) write(*,*) "The solution converges to a local minimum. Please restart &
	  &  with different number of load steps and number of iterations"
	  
	return
	end subroutine