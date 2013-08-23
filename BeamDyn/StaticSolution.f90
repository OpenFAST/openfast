   SUBROUTINE StaticSolution(uuN0,uuNf,hhp,w,Jacobian,Stif0,F_ext,&
                            &node_elem,dof_node,norder,elem_total,dof_total,node_total,niter)

   REAL(ReKi),INTENT(IN)::uuN0(:),hhp(:,:),w(:),Stif0(:,:,:),F_ext(:),Jacobian

   REAL(Reki),INTENT(INOUT)::uuNf(:)
   INTEGER,INTENT(IN)::node_elem,dof_node,norder,elem_total,dof_total,node_total,niter
 
   REAL(ReKi)::errf,temp1,temp2,errx
   REAL(ReKi)::ui(dof_total),bc(dof_total),rel_change(dof_total),ui_old(dof_total)
   PARAMETER(TOLF = 1.0d-5)   

   INTEGER::i,j

   DO i=1,niter
       CALL BeamGenerateStaticMatrix(uuN0,uuNf,hhp,w,Jacobian,Stif0,&
                                     &node_elem,dof_node,norder,elem_total,dof_total,node_total,&
                                     &StifK,RHS)
       CALL AppliedLoad(F_ext,hhp,w,Jacobian,elem_total,dof_total,RHS)
       CALL Norm(dof_total,RHS,errf)
       IF(errf .LE. TOLF) RETURN
       ui_old = ui
       bc = 1.0d0

       DO j=1,dof_node
           bc(j) = 0.0d0
       ENDDO
           
       CALL CGSolver(RHS,StifK,ui,bc,dof_total)
       rel_change = ABS(ui - ui_old)
       CALL Norm(dof_total,uf,temp1)
       CALL Norm(dof_total,rel_change,temp2)
       errx = temp2
       IF(errx .LT. TOLF*temp1) RETURN
           
       uf = uf + ui
           
       IF(i==niter) THEN
           WRITE(*,*) "Solution does not converge after the maximum number of iterations"
           STOP
       ENDIF
   ENDDO
   
   END SUBROUTINE StaticSolution
