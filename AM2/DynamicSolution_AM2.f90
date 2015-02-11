   SUBROUTINE DynamicSolution_AM2(uuN0,uuN,vvN,uuN00,vvN00,uud0,vvd0,&
                                  Stif0,Mass0,gravity,u,&
                                  damp_flag,beta,&
                                  node_elem,dof_node,elem_total,dof_total,&
                                  node_total,ngp,niter,nr_counter,dt,alpha)

   REAL(ReKi),        INTENT(IN   ):: uuN0(:,:)
   REAL(ReKi),        INTENT(IN   ):: uuN00(:)
   REAL(ReKi),        INTENT(IN   ):: vvN00(:)
   REAL(ReKi),        INTENT(IN   ):: uud0(:)
   REAL(ReKi),        INTENT(IN   ):: vvd0(:)
   REAL(ReKi),        INTENT(IN   ):: Stif0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: Mass0(:,:,:)
   REAL(ReKi),        INTENT(IN   ):: gravity(:)
   REAL(DbKi),        INTENT(IN   ):: dt
   REAL(DbKi),        INTENT(IN   ):: alpha
   TYPE(BD_InputType),INTENT(IN   ):: u
   INTEGER(IntKi),    INTENT(IN   ):: node_elem
   INTEGER(IntKi),    INTENT(IN   ):: dof_node
   INTEGER(IntKi),    INTENT(IN   ):: elem_total
   INTEGER(IntKi),    INTENT(IN   ):: dof_total
   INTEGER(IntKi),    INTENT(IN   ):: node_total
   INTEGER(IntKi),    INTENT(IN   ):: ngp
   INTEGER(IntKi),    INTENT(IN   ):: niter
   INTEGER(IntKi),    INTENT(IN   ):: damp_flag
   REAL(ReKi),        INTENT(IN   ):: beta(:)
   REAL(ReKi),        INTENT(INOUT):: uuN(:)
   REAL(ReKi),        INTENT(INOUT):: vvN(:)
   INTEGER(IntKi),    INTENT(INOUT):: nr_counter

   REAL(ReKi)                      :: MassM(dof_total*2,dof_total*2)
   REAL(ReKi)                      :: RHS(dof_total*2)
   REAL(ReKi)                      :: MassM_LU(dof_total*2-12,dof_total*2-12)
   REAL(ReKi)                      :: RHS_LU(dof_total*2-12)
   REAL(ReKi)                      :: F_PointLoad(dof_total)
   REAL(ReKi)                      :: F_PointLoad0(dof_total)
   REAL(ReKi)                      :: feqv(dof_total-6)
   REAL(ReKi)                      :: sol_temp(dof_total*2-12)
   REAL(ReKi)                      :: sol(dof_total*2)
   REAL(ReKi)                      :: d
   REAL(ReKi)                      :: temp
   REAL(ReKi)                      :: Enorm
   REAL(ReKi)                      :: Eref
   REAL(ReKi),            PARAMETER:: TOLF = 1.0D-03
   REAL(ReKi),            PARAMETER:: TOLF1 = 1.0D-04
   INTEGER(IntKi)                  :: indx(dof_total*2-12)
   INTEGER(IntKi)                  :: temp_id
   INTEGER(IntKi)                  :: i
   INTEGER(IntKi)                  :: j
   INTEGER(IntKi)                  :: k

   Eref = 0.0D0

   DO i=1,niter
       RHS(:) = 0.0D0
       MassM(:,:) = 0.0D0
!WRITE(*,*) "niter = ",i
       CALL GenerateDynamicElement_AM2(uuN0,uuN,vvN,uuN00,vvN00,uud0,vvd0,&
                                      &Stif0,Mass0,gravity,u,&
                                      &damp_flag,beta,&
                                      &elem_total,node_elem,dof_total,dof_node,ngp,dt,alpha,RHS,MassM)
!DO j=1,36
!   WRITE(*,*) j,RHS(j)
!ENDDO
       DO j=1,node_total
           temp_id = (j-1)*dof_node
           F_PointLoad(temp_id+1:temp_id+3) = u%PointLoad%Force(1:3,j)
           F_PointLoad(temp_id+4:temp_id+6) = u%PointLoad%Moment(1:3,j)
       ENDDO
       RHS(dof_total+1:dof_total*2) = RHS(dof_total+1:dof_total*2) + &
                                     &dt*alpha*F_PointLoad(1:dof_total)

       feqv(:) = 0.0D0
       DO j=1,dof_total-6
           feqv(j) = RHS(j+6)
           RHS_LU(j) = RHS(j+6)
           RHS_LU(j+dof_total-6) = RHS(j+6+dof_total)
           DO k=1,dof_total-6
               MassM_LU(j,k) = MassM(j+6,k+6)
               MassM_LU(j,k+dof_total-6) = MassM(j+6,k+dof_total+6)
               MassM_LU(j+dof_total-6,k) = MassM(j+dof_total+6,k+6)
               MassM_LU(j+dof_total-6,k+dof_total-6) = MassM(j+dof_total+6,k+dof_total+6)
           ENDDO
       ENDDO

       temp = Norm(RHS_LU)
WRITE(*,*) i, temp
IF(temp .LT. TOLF1) THEN
!IF(i .EQ. 3) THEN
   nr_counter = i
!   WRITE(*,*) i
   RETURN
ENDIF
       CALL ludcmp(MassM_LU,dof_total*2-12,indx,d)
       CALL lubksb(MassM_LU,dof_total*2-12,indx,RHS_LU,sol_temp)

!temp = MAXVAL(ABS(sol_temp))
!WRITE(*,*) i,temp
!DO j=1,24
!   WRITE(*,*) sol_temp(j)
!ENDDO

       sol(:) = 0.0D0
       DO j=1,dof_total-6
           sol(j+6) = sol_temp(j)
           sol(j+dof_total+6) = sol_temp(j+dof_total-6)
       ENDDO
!DO j=1,36
!WRITE(*,*) "sol(j) = ",j,sol(j)
!WRITE(*,*) "feqv(j) = ",j,feqv(j)
!ENDDO

!       IF(i==1) Eref = TOLF * DOT_PRODUCT(sol_temp(1:dof_total-6),feqv)
!       IF(i .GT. 1) THEN
!           Enorm = 0.0D0 
!           Enorm = DOT_PRODUCT(sol_temp(1:dof_total-6),feqv)
!           IF(Enorm .GT. Eref/TOLF) THEN
!           ELSEIF(Enorm .LE. Eref) THEN
!               RETURN
!           ENDIF
!       ENDIF
       CALL UpdateConfiguration(sol(1:dof_total),uuN,node_total,dof_node)
       CALL UpdateConfiguration_AM2(sol(dof_total+1:dof_total*2),vvN,node_total,dof_node)

       IF(i==niter) THEN
           WRITE(*,*) "Solution does not converge after the maximum number of iterations"
           STOP !! USE EXIT INSTEAD OF STOP, NJ 3/18/2014
       ENDIF
   ENDDO

   END SUBROUTINE DynamicSolution_AM2
