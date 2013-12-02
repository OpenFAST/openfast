   SUBROUTINE StaticSolutionGL(uuN0,uuNf,Stif0,F_ext,bc,&
                            &node_elem,dof_node,elem_total,dof_total,node_total,ngp,niter)

   REAL(ReKi),INTENT(IN):: uuN0(:),F_ext(:),Stif0(:,:),bc(:)
   REAL(ReKi),INTENT(INOUT):: uuNf(:)
   INTEGER(IntKi),INTENT(IN):: niter,elem_total,node_elem,dof_node,ngp,dof_total,node_total

   REAL(ReKi):: StifK(dof_total,dof_total),RHS(dof_total)

   REAL(ReKi):: StifK_LU(dof_total-6,dof_total-6),RHS_LU(dof_total-6)

   REAL(ReKi):: ui(dof_total),ui_temp(dof_total-6)
   REAL(ReKi):: Eref,Enorm,feqv(dof_total-6),errf
   REAL(ReKi),PARAMETER:: TOLF = 1.0D-10
   REAL(ReKi):: d

   INTEGER(IntKi):: indx(dof_total-6)
   INTEGER(IntKi):: i,j
   INTEGER(IntKi):: temp_id,k !For Debug

   ui = 0.0D0
   Eref = 0.0D0

   DO i=1,niter
       WRITE(*,*) "N-R Iteration #:", i
!       IF(i==2) STOP
       StifK = 0.0D0
       RHS = 0.0D0
       CALL BldGenerateStaticElement(uuN0,uuNf,F_ext,Stif0,elem_total,node_elem,dof_node,ngp,StifK,RHS)
       
       RHS = RHS + F_ext
       
       feqv = 0.0D0
       errf = 0.0D0
       DO j=1,dof_total-6
           feqv(j) = RHS(j+6)
           RHS_LU(j) = RHS(j+6)
           DO k=1,dof_total-6
               StifK_LU(j,k) = StifK(j+6,k+6)
           ENDDO 
       ENDDO
       CALL Norm(dof_total-6,feqv,errf)
       WRITE(*,*) "NORM(feqv) = ", errf
   

!       CALL CGSolver(RHS,StifK,ui,bc,dof_total)
       CALL ludcmp(StifK_LU,dof_total-6,indx,d)
       CALL lubksb(StifK_LU,dof_total-6,indx,RHS_LU,ui_temp)       

       ui = 0.0D0
       DO j=1,dof_total-6
           ui(j+6) = ui_temp(j)
       ENDDO
!       ui_temp = 0.0D0
!       DO j=1,dof_total-6
!           ui_temp(j) = ui(j+6)
!       ENDDO
       IF(i==1) Eref = SQRT(DOT_PRODUCT(ui_temp,feqv)*TOLF)
       IF(i .GT. 1) THEN
           Enorm = 0.0D0 
           Enorm = SQRT(DOT_PRODUCT(ui_temp,feqv))
           WRITE(*,*) "Enorm = ", Enorm
           WRITE(*,*) "Eref = ", Eref
           IF(Enorm .LE. Eref) RETURN
       ENDIF
           
       CALL UpdateConfiguration(ui,uuNf,node_total,dof_node)
!       DO j=1,node_total
!           DO k=1,dof_node
!               WRITE(*,*) "uuNf",k," = ", uuNf((j-1)*dof_node+k)
!           ENDDO
!       ENDDO
       IF(i==niter) THEN
           WRITE(*,*) "Solution does not converge after the maximum number of iterations"
           STOP
       ENDIF
   ENDDO
   
   END SUBROUTINE StaticSolutionGL

