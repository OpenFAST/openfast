   SUBROUTINE BeamDyn_StaticSolution(uuN0,uuNf,Mass0,Stif0,gravity,u,&
                                     node_elem,dof_node,elem_total,&
                                     dof_total,node_total,ngp,niter,piter)

   REAL(ReKi),        INTENT(IN   ):: uuN0(:,:)
   REAL(ReKi),        INTENT(IN   ):: Mass0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: Stif0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: gravity(:)
   TYPE(BD_InputType),INTENT(IN   ):: u
   INTEGER(IntKi),    INTENT(IN   ):: niter
   INTEGER(IntKi),    INTENT(IN   ):: elem_total
   INTEGER(IntKi),    INTENT(IN   ):: node_elem
   INTEGER(IntKi),    INTENT(IN   ):: dof_node
   INTEGER(IntKi),    INTENT(IN   ):: ngp
   INTEGER(IntKi),    INTENT(IN   ):: dof_total
   INTEGER(IntKi),    INTENT(IN   ):: node_total
   REAL(ReKi),        INTENT(INOUT):: uuNf(:)
   INTEGER(IntKi),    INTENT(  OUT):: piter !! ADDED piter AS OUTPUT

   REAL(ReKi)                      :: StifK(dof_total,dof_total)
   REAL(ReKi)                      :: RHS(dof_total)
   REAL(ReKi)                      :: StifK_LU(dof_total-6,dof_total-6)
   REAL(ReKi)                      :: RHS_LU(dof_total-6)
   REAL(ReKi)                      :: ui(dof_total)
   REAL(ReKi)                      :: ui_temp(dof_total-6)
   REAL(ReKi)                      :: feqv(dof_total-6)
   REAL(ReKi)                      :: Eref
   REAL(ReKi)                      :: Enorm
   REAL(ReKi)                      :: temp
   REAL(ReKi),            PARAMETER:: TOLF = 1.0D-10
   REAL(ReKi)                      :: d
   INTEGER(IntKi)                  :: indx(dof_total-6)
   INTEGER(IntKi)                  :: i
   INTEGER(IntKi)                  :: j
   INTEGER(IntKi)                  :: k
   INTEGER(IntKi)                  :: temp_id
!   INTEGER(IntKi)                  :: temp_count

   Eref = 0.0D0
   DO i=1,niter
       StifK(:,:) = 0.0D0
       RHS(:)     = 0.0D0
       CALL BeamDyn_GenerateStaticElement(uuN0,uuNf,Mass0,Stif0,gravity,u,&
                                          elem_total,node_elem,dof_node,ngp,StifK,RHS)
WRITE(*,*) "niter = ",i

       DO j=1,node_total
           temp_id = (j-1)*dof_node
           RHS(temp_id+1:temp_id+3) = RHS(temp_id+1:temp_id+3) + u%Pointload%Force(1:3,j)
           RHS(temp_id+4:temp_id+6) = RHS(temp_id+4:temp_id+6) + u%Pointload%Moment(1:3,j)
       ENDDO

       feqv(:) = 0.0D0
       DO j=1,dof_total-6
           feqv(j) = RHS(j+6)
           RHS_LU(j) = RHS(j+6)
           DO k=1,dof_total-6
               StifK_LU(j,k) = StifK(j+6,k+6)
           ENDDO 
       ENDDO

       CALL ludcmp(StifK_LU,dof_total-6,indx,d)
       CALL lubksb(StifK_LU,dof_total-6,indx,RHS_LU,ui_temp)       

       ui(:) = 0.0D0
       DO j=1,dof_total-6
           ui(j+6) = ui_temp(j)
       ENDDO

       temp = Norm(feqv)
       WRITE(13,*) i,temp 
       piter=i !! ADDED BY NJ 3/18/14
!       IF(i==1) Eref = SQRT(DOT_PRODUCT(ui_temp,feqv)*TOLF)
!       IF(i .GT. 1) THEN
!           Enorm = 0.0D0 
!           Enorm = SQRT(DOT_PRODUCT(ui_temp,feqv))
!           IF(Enorm .LE. Eref) RETURN
!       ENDIF
       IF(i==1) Eref = TOLF * DOT_PRODUCT(ui_temp,feqv)
       IF(i==1) WRITE(*,*) "Eref = ", Eref
       IF(i .GT. 1) THEN
           Enorm = 0.0D0 
           Enorm = DOT_PRODUCT(ui_temp,feqv)
           WRITE(*,*) "Enorm = ", Enorm
           IF(Enorm .GT. Eref/TOLF) THEN
               WRITE(*,*) "Solution is diverging, exit N-R"
               piter=niter 
           ELSEIF(Enorm .LE. Eref) THEN
               RETURN
           ENDIF
       ENDIF
           
       CALL UpdateConfiguration(ui,uuNf,node_total,dof_node)
       IF(i==niter .OR. piter==niter) THEN
           WRITE(*,*) "Solution does not converge after the maximum number of iterations"
           EXIT !! USE EXIT INSTEAD OF STOP, NJ 3/18/2014
       ENDIF
   ENDDO
   
   END SUBROUTINE BeamDyn_StaticSolution

