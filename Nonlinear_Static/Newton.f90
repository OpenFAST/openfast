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

	ui=0.0d0
	call AssembleRHS(RHS,dof_node,dof_total,uf,&
							&norder,hhp,wj,node_total,dmat,&
							&elem_total)
							
	if(MAXVAL(ABS(RHS))<TOLF) return
    
	fmin=0.5d0*DOT_PRODUCT(RHS,RHS)  ! calculate fnew=1/2*F.F
    
    temp1 = float(size(UF))
	
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
   		call CGSolver(RHS,KT,ui,bc,dof_total)
   
   		call LineSearch(dof_total,xold,fold,gradient,rhs,ui,fmin,stpmax,check,&
							&dof_node,uf,norder,hhp,wj,node_total,dmat,elem_total)
							
   		IF(MAXVAL(ABS(RHS))<TOLF) THEN
       		check=.FALSE.  ! the function values converge to zero
	   		RETURN 
   		ENDIF 
   	
   		IF(check) THEN  ! converges on deltaX
      	! Check for gradient of 1/2 F.F zero, i.e. spurious convergence
	  	!---------------------------------------------------------------------
	   		check=(MAXVAL(ABS(gradient)*MAX(ABS(uf),1.0D0)/MAX(fmin,0.5D0*SIZE(uf)))<TOLF) 
	   		EXIT
   		ENDIF   

   		IF(i==niter) write(*,*) "The solution does not converge after the maximum number of iterations" !The maximum number of iterations reached, the solution is still not converged. 
	enddo

	if(check) write(*,*) "The solution converges to a local minimum. Please restart &
	  &  with different number of load steps and number of iterations"
	  
	return
	end subroutine