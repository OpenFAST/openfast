   SUBROUTINE StaticSolution(uuN0,uuNf,hhp,w,Jacobian,Stif0,F_ext,bc,&
                            &node_elem,dof_node,norder,elem_total,dof_total,node_total,dof_elem,niter)

   REAL(ReKi),INTENT(IN)::uuN0(:),hhp(:,:),w(:),Stif0(:,:,:),F_ext(:),bc(:),Jacobian

   REAL(Reki),INTENT(INOUT)::uuNf(:)
   INTEGER(IntKi),INTENT(IN)::node_elem,dof_node,norder,elem_total,dof_total,node_total,niter
   INTEGER(IntKi),INTENT(IN)::dof_elem
 
   REAL(ReKi)::errf,temp1,temp2,errx
   REAL(ReKi)::StifK(dof_total,dof_total), RHS(dof_total)
   REAL(ReKi)::ui(dof_total),rel_change(dof_total),ui_old(dof_total)
   REAL(ReKi),PARAMETER:: TOLF = 1.0d-5   

   INTEGER(IntKi)::i,j,k

   ui = 0.0D0

   DO i=1,niter
       WRITE(*,*) "N-R Iteration #:", i
       IF(i==6) STOP
       StifK = 0.0D0
       RHS = 0.0D0
       CALL BeamStatic(uuN0,uuNf,hhp,w,Jacobian,Stif0,F_ext,&
                      &node_elem,dof_node,norder,elem_total,dof_total,node_total,dof_elem,&
                      &StifK,RHS)
!       k=12
!       DO j=1,dof_total
!           WRITE(*,*) StifK(j,k+1),StifK(j,k+2),StifK(j,k+3),StifK(j,k+4),StifK(j,k+5),StifK(j,k+6)
!       ENDDO
!       STOP
!       RHS(dof_total-1) = RHS(dof_total-1) - 3.14159D+01
!       DO j=1,dof_total
!           WRITE(*,*) "j=",j
!           WRITE(*,*) RHS(j)
!       ENDDO
!       STOP
       RHS = RHS + F_ext
       errf = 0.0D0
       DO j=1,dof_node
           RHS(j) = 0.0D0   
       ENDDO
       CALL Norm(dof_total,RHS,errf)
       WRITE(*,*) "NORM(RHS) = ", errf
!       IF(errf .LE. TOLF) RETURN
       ui_old = ui
       
       CALL CGSolver(RHS,StifK,ui,bc,dof_total)
!       DO j=1,dof_total
!           WRITE(*,*) "j=",j
!           WRITE(*,*) "ui(j)=",ui(j)
!       ENDDO
!       STOP
       rel_change = ABS(ui - ui_old)
       CALL Norm(dof_total,uuNf,temp1)
       CALL Norm(dof_total,rel_change,temp2)
       errx = temp2
!       IF(errx .LT. TOLF*temp1) RETURN
           
       CALL UpdateConfiguration(ui,uuNf,node_total,dof_node)
!       IF(i==20) THEN    
       DO j=1,dof_total
!           WRITE(*,*) "j=",j
           WRITE(*,*) uuNf(j)
       ENDDO
!       ENDIF
       IF(i==niter) THEN
           WRITE(*,*) "Solution does not converge after the maximum number of iterations"
           STOP
       ENDIF
   ENDDO
   
   END SUBROUTINE StaticSolution
