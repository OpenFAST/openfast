   SUBROUTINE ElementMatrixDynLSGL(Nuu0,Nuuu,Nrr0,Nrrr,Next,Nvvv,Naaa,&
                                  &Stif0,m00,mEta0,rho0,ngp,node_elem,dof_node,elk,elf,elm,elg)

   REAL(ReKi),INTENT(IN):: Nuu0(:),Nuuu(:),Nrr0(:),Nrrr(:),Next(:)
   REAL(ReKi),INTENT(IN):: Nvvv(:),Naaa(:)
   REAL(ReKi),INTENT(IN):: Stif0(:,:)
   REAL(ReKi),INTENT(IN):: m00,mEta0(:),rho0(:,:)
   INTEGER(IntKi),INTENT(IN):: ngp,node_elem,dof_node

   REAL(ReKi),INTENT(OUT):: elk(:,:),elf(:)  
   REAL(ReKi),INTENT(OUT):: elg(:,:),elm(:,:)    

   REAL(ReKi),ALLOCATABLE:: gp(:),gw(:),hhx(:),hpx(:)
   REAL(ReKi),ALLOCATABLE:: GLL_temp(:),w_temp(:)
   
   REAL(ReKi):: uu0(6),E10(3),RR0(3,3),kapa(3),E1(3),Stif(6,6),cet
   REAL(ReKi):: uuu(6),uup(3),Jacobian,gpr
   REAL(ReKi):: Fc(6),Fd(6),Oe(6,6),Pe(6,6),Qe(6,6)
   REAL(ReKi):: vvv(6),aaa(6),mEta(3),rho(3,3)
   REAL(ReKi):: Fi(6),Mi(6,6),Ki(6,6),Gi(6,6)

   INTEGER(IntKi)::igp,i,j,m,n,temp_id1,temp_id2
   INTEGER(IntKi)::allo_stat


   ALLOCATE(gp(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gp = 0.0D0

   ALLOCATE(gw(ngp), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   gw = 0.0D0

   ALLOCATE(hhx(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hhx = 0.0D0

   ALLOCATE(hpx(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   hpx = 0.0D0
   
   ALLOCATE(GLL_temp(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   GLL_temp = 0.0D0
   
   ALLOCATE(w_temp(node_elem), STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   w_temp = 0.0D0

   elk = 0.0D0
   elf = 0.0D0
   elg = 0.0D0
   elm = 0.0D0
   
   CALL BDyn_gen_gll_LSGL(node_elem-1,GLL_temp,w_temp)
!   CALL BldSet1DGaussPointScheme(ngp,gp,gw)
   CALL BldGaussPointWeight(ngp,gp,gw)
   DO igp=1,ngp
       gpr = gp(igp)
       CALL BldComputeJacobianLSGL(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian)
!       WRITE(*,*) "Jacobian = ", Jacobian
       CALL BldGaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10)
       CALL BldGaussPointData(hhx,hpx,Nuuu,Nrrr,Stif0,uu0,E10,node_elem,dof_node,uuu,uup,E1,RR0,kapa,Stif,cet)
       CALL ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)
       
       CALL BldGaussPointDataMass(hhx,hpx,Nvvv,Naaa,RR0,mEta0,rho0,node_elem,dof_node,vvv,aaa,mEta,rho)
       CALL InertialForce(m00,mEta,rho,vvv,aaa,Fi,Mi,Gi,Ki)

       DO i=1,node_elem
           DO j=1,node_elem
               DO m=1,dof_node
                   temp_id1 = (i-1)*dof_node+m
                   DO n=1,dof_node
                       temp_id2 = (j-1)*dof_node+n
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hhx(i)*Qe(m,n)*hhx(j)*Jacobian*gw(igp)
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hhx(i)*Pe(m,n)*hpx(j)*Jacobian*gw(igp)
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hpx(i)*Oe(m,n)*hhx(j)*Jacobian*gw(igp)
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hpx(i)*Stif(m,n)*hpx(j)*Jacobian*gw(igp)
                       elk(temp_id1,temp_id2) = elk(temp_id1,temp_id2) + hhx(i)*Ki(m,n)*hhx(j)*Jacobian*gw(igp)
                       elm(temp_id1,temp_id2) = elm(temp_id1,temp_id2) + hhx(i)*Mi(m,n)*hhx(j)*Jacobian*gw(igp)
                       elg(temp_id1,temp_id2) = elg(temp_id1,temp_id2) + hhx(i)*Gi(m,n)*hhx(j)*Jacobian*gw(igp)
                   ENDDO
               ENDDO
           ENDDO
       ENDDO 

       DO i=1,node_elem
           DO j=1,dof_node
               temp_id1 = (i-1) * dof_node+j
               elf(temp_id1) = elf(temp_id1) - hhx(i)*Fd(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) - hpx(i)*Fc(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) - hhx(i)*Fi(j)*Jacobian*gw(igp)
           ENDDO
       ENDDO

   ENDDO

   DEALLOCATE(gp)
   DEALLOCATE(gw)
   DEALLOCATE(hhx)
   DEALLOCATE(hpx)
   DEALLOCATE(GLL_temp)
   DEALLOCATE(w_temp)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(gp))  DEALLOCATE(gp)
            IF(ALLOCATED(gw))  DEALLOCATE(gw)
            IF(ALLOCATED(hhx)) DEALLOCATE(hhx)
            IF(ALLOCATED(hpx)) DEALLOCATE(hpx)
            IF(ALLOCATED(GLL_temp)) DEALLOCATE(GLL_temp)
            IF(ALLOCATED(w_temp)) DEALLOCATE(w_temp)
        ENDIF

   END SUBROUTINE ElementMatrixDynLSGL
