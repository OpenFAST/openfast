   SUBROUTINE StaticSolution(uuN0,uuNf,hhp,w,Jacobian,Stif0,F_ext,bc,&
                            &node_elem,dof_node,norder,elem_total,dof_total,node_total,dof_elem,niter)

   REAL(ReKi),INTENT(IN)::uuN0(:),hhp(:,:),w(:),Stif0(:,:,:),F_ext(:),bc(:),Jacobian

   REAL(Reki),INTENT(INOUT)::uuNf(:)
   INTEGER(IntKi),INTENT(IN)::node_elem,dof_node,norder,elem_total,dof_total,node_total,niter
   INTEGER(IntKi),INTENT(IN)::dof_elem
 
   REAL(ReKi)::errf,temp1,temp2,errx
   REAL(ReKi)::StifK(dof_total,dof_total), RHS(dof_total)
   REAL(ReKi)::ui(dof_total)
   REAL(ReKi)::feqv(dof_total-6),Eref,ui_temp(dof_total-6),Enorm
   REAL(ReKi),PARAMETER:: TOLF = 1.0D-07   

   INTEGER(IntKi)::i,j,k

   ui = 0.0D0
   Eref = 0.0D0


   DO i=1,niter
       WRITE(*,*) "N-R Iteration #:", i
       IF(i==20) STOP
       StifK = 0.0D0
       RHS = 0.0D0
       CALL BeamStatic(uuN0,uuNf,hhp,w,Jacobian,Stif0,F_ext,&
                      &node_elem,dof_node,norder,elem_total,dof_total,node_total,dof_elem,&
                      &StifK,RHS)
       RHS = RHS + F_ext
       feqv = 0.0D0
       errf = 0.0D0
       DO j=1,dof_total-6
           feqv(j) = RHS(j+6) 
       ENDDO
       CALL Norm(dof_total-6,feqv,errf)
       WRITE(*,*) "NORM(feqv) = ", errf

       CALL CGSolver(RHS,StifK,ui,bc,dof_total)
       ui_temp = 0.0D0
       DO j=1,dof_total-6
           ui_temp(j) = ui(j+6)
       ENDDO
       IF(i==1) Eref = DOT_PRODUCT(ui_temp,feqv)*TOLF
       IF(i .GT. 1) THEN
           Enorm = 0.0D0 
           Enorm = DOT_PRODUCT(ui_temp,feqv)
           WRITE(*,*) "Enorm = ", Enorm
           WRITE(*,*) "Eref = ", Eref
           IF(Enorm .LE. Eref) RETURN
       ENDIF
           
       CALL UpdateConfiguration(ui,uuNf,node_total,dof_node)
       IF(i==niter) THEN
           WRITE(*,*) "Solution does not converge after the maximum number of iterations"
           STOP
       ENDIF
   ENDDO
   
   END SUBROUTINE StaticSolution
