   SUBROUTINE DynamicSolution(uuN0,uuNi,vvNi,aaNi,xxNi,&
                            &uuNf,vvNf,aaNf,xxNf,hhp,w,Jacobian,Stif0,m00,mEta0,rho0,bc,&
                            &node_elem,dof_node,norder,elem_total,dof_total,node_total,dof_elem,coef,niter,&
                            &deltat,time)

   REAL(ReKi),INTENT(IN)::uuN0(:),hhp(:,:),w(:),Stif0(:,:,:),bc(:),Jacobian
   REAL(ReKi),INTENT(IN)::m00,mEta0(:),rho0(:,:)

   REAL(ReKi),INTENT(INOUT)::uuNi(:),vvNi(:),aaNi(:),xxNi(:)

   REAL(ReKi),INTENT(INOUT)::vvNf(:),aaNf(:),xxNf(:)
   REAL(Reki),INTENT(INOUT)::uuNf(:)

   REAL(DbKi),INTENT(IN):: coef(:), deltat,time

   INTEGER(IntKi),INTENT(IN)::node_elem,dof_node,norder,elem_total,dof_total,node_total,niter
   INTEGER(IntKi),INTENT(IN)::dof_elem
 
   REAL(ReKi)::errf,temp1,temp2,errx
   REAL(ReKi)::StifK(dof_total,dof_total), RHS(dof_total)
   REAL(ReKi)::MassM(dof_total,dof_total), DampG(dof_total,dof_total)
   REAL(ReKi)::F_ext(dof_total)

   REAL(ReKi)::ai(dof_total),rel_change(dof_total),ai_old(dof_total)
   REAL(ReKi),PARAMETER:: TOLF = 1.0d-5   

   INTEGER(IntKi)::i,j,k

   CALL TiSchmPredictorStep(uuNi,vvNi,aaNi,xxNi,coef,deltat,uuNf,vvNf,aaNf,xxNf,node_total,dof_node)
   CALL AppliedNodalLoad(F_ext,time,dof_total)
   ai = 0.0D0

   DO i=1,niter
       WRITE(*,*) "N-R Iteration #", i
!       IF(i==10) STOP
       StifK = 0.0D0
       RHS = 0.0D0
       MassM = 0.0D0
       DampG = 0.0D0

       CALL BeamDynamic(uuN0,uuNf,vvNf,aaNf,hhp,w,Jacobian,Stif0,m00,mEta0,rho0,&
                      &node_elem,dof_node,norder,elem_total,dof_total,node_total,dof_elem,&
                      &StifK,MassM,DampG,RHS)
!       k=0
!       DO j=1,dof_total
!           WRITE(*,*) "j=",j
!           WRITE(*,*) "RHS(j)=",RHS(j)
!           WRITE(*,*) StifK(j,1+k),StifK(j,2+k),StifK(j,3+k),StifK(j,4+k),StifK(j,5+k),StifK(j,6+k)
!           WRITE(*,*) MassM(j,1+k),MassM(j,2+k),MassM(j,3+k),MassM(j,4+k),MassM(j,5+k),MassM(j,6+k)
!           WRITE(*,*) DampG(j,1+k),DampG(j,2+k),DampG(j,3+k),DampG(j,4+k),DampG(j,5+k),DampG(j,6+k)
!       ENDDO
!       STOP
       StifK = MassM + coef(7) * DampG + coef(8) * StifK
       RHS = RHS + F_ext
!       k=0
!       DO j=1,dof_total
!           WRITE(*,*) "j=",j
!           WRITE(*,*) "RHS(j)=",RHS(j)
!           WRITE(*,*) StifK(j,1+k),StifK(j,2+k),StifK(j,3+k),StifK(j,4+k),StifK(j,5+k),StifK(j,6+k)
!       ENDDO
!       STOP
       errf = 0.0D0
       CALL Norm(dof_total,RHS,errf)
       WRITE(*,*) "Norm(RHS) = ", errf
       IF(errf .LE. TOLF) RETURN
       ai_old = ai
       
       CALL CGSolver(RHS,StifK,ai,bc,dof_total)
!       DO j=1,dof_total
!           WRITE(*,*) "j=",j
!           WRITE(*,*) "ai(j)=",ai(j)
!       ENDDO
!       STOP
       rel_change = ABS(ai - ai_old)
       CALL Norm(dof_total,aaNf,temp1)
       CALL Norm(dof_total,rel_change,temp2)
       errx = temp2
       IF(errx .LT. TOLF*temp1) RETURN
           
       CALL UpdateDynamic(ai,uuNf,vvNf,aaNf,xxNf,coef,node_total,dof_node)
           
!       DO j=1,dof_total
!           WRITE(*,*) "j=",j
!           WRITE(*,*) "xxNf(j)=",xxNf(j)
!       ENDDO
!       STOP
       IF(i==niter) THEN
           WRITE(*,*) "Solution does not converge after the maximum number of iterations"
           STOP
       ENDIF
   ENDDO
   
   END SUBROUTINE DynamicSolution
