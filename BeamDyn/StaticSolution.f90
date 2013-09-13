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

   INTEGER(IntKi)::i,j

   ui = 0.0D0

   DO i=1,niter
       StifK = 0.0D0
       RHS = 0.0D0
       CALL BeamStatic(uuN0,uuNf,hhp,w,Jacobian,Stif0,F_ext,&
                      &node_elem,dof_node,norder,elem_total,dof_total,node_total,dof_elem,&
                      &StifK,RHS)
       errf = 0.0D0
       CALL Norm(dof_total,RHS,errf)
       IF(errf .LE. TOLF) RETURN
       ui_old = ui
       WRITE(*,*) "# of N-R Iteration", i
       CALL CGSolver(RHS,StifK,ui,bc,dof_total)
       DO j=1,dof_total
           WRITE(*,*) "j=",j
           WRITE(*,*) ui(j)
       ENDDO
       rel_change = ABS(ui - ui_old)
       CALL Norm(dof_total,uuNf,temp1)
       CALL Norm(dof_total,rel_change,temp2)
       errx = temp2
       IF(errx .LT. TOLF*temp1) RETURN
           
       CALL UpdateConfiguration(ui,uuNf,node_total,dof_node)
           
       IF(i==niter) THEN
           WRITE(*,*) "Solution does not converge after the maximum number of iterations"
           STOP
       ENDIF
   ENDDO
   
   END SUBROUTINE StaticSolution
