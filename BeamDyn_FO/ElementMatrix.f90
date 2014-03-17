   SUBROUTINE ElementMatrix(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,NStif0,Nm00,NmEta0,Nrho0,&
                           &ngp,node_elem,dof_node,elf,elm)

   REAL(ReKi),INTENT(IN):: Nuu0(:),Nuuu(:),Nrr0(:),Nrrr(:)
   REAL(ReKi),INTENT(IN):: Nvvv(:)
   REAL(ReKi),INTENT(IN):: NStif0(:,:,:)
   REAL(ReKi),INTENT(IN):: Nm00(:),NmEta0(:,:),Nrho0(:,:,:)
   INTEGER(IntKi),INTENT(IN):: ngp,node_elem,dof_node

   REAL(ReKi),INTENT(OUT):: elf(:)  
   REAL(ReKi),INTENT(OUT):: elm(:,:)    

   REAL(ReKi),ALLOCATABLE:: gp(:),gw(:),hhx(:),hpx(:)
   REAL(ReKi),ALLOCATABLE:: GLL_temp(:),w_temp(:)
   
   REAL(ReKi):: uu0(6),E10(3),RR0(3,3),kapa(3),E1(3),Stif(6,6),cet
   REAL(ReKi):: uuu(6),uup(3),Jacobian,gpr
   REAL(ReKi):: Fc(6),Fd(6)
   REAL(ReKi):: vvv(6)
   REAL(ReKi):: mmm,mEta(3),rho(3,3)
   REAL(ReKi):: Fb(6),Mi(6,6)

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

   elf = 0.0D0
   elm = 0.0D0


   CALL BDyn_gen_gll_LSGL(node_elem-1,GLL_temp,w_temp)
   CALL BldGaussPointWeight(ngp,gp,gw)

   DO igp=1,ngp
       gpr=gp(igp)
       CALL BldComputeJacobianLSGL(gpr,Nuu0,node_elem,dof_node,gp,GLL_temp,ngp,igp,hhx,hpx,Jacobian)

       CALL BldGaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,NStif0,node_elem,dof_node,uu0,E10,Stif)
       CALL BldGaussPointData(hhx,hpx,Nuuu,Nrrr,uu0,E10,node_elem,dof_node,uuu,uup,E1,RR0,kapa,Stif,cet)       
       CALL BldGaussPointDataMass(hhx,hpx,Nvvv,RR0,Nm00,NmEta0,Nrho0,node_elem,dof_node,vvv,mmm,mEta,rho)

       CALL ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd)
       CALL MassMatrix(mmm,mEta,rho,uuu,Mi)
       CALL GyroForce(mEta,rho,uuu,vvv,Fb)

       DO i=1,node_elem
           DO j=1,node_elem
               DO m=1,dof_node
                   temp_id1 = (i-1)*dof_node+m
                   DO n=1,dof_node
                       temp_id2 = (j-1)*dof_node+n
                       elm(temp_id1,temp_id2) = elm(temp_id1,temp_id2) + hhx(i)*Mi(m,n)*hhx(j)*Jacobian*gw(igp)
                   ENDDO
               ENDDO
           ENDDO
       ENDDO

       DO i=1,node_elem
           DO j=1,dof_node
               temp_id1 = (i-1) * dof_node+j
               elf(temp_id1) = elf(temp_id1) - hhx(i)*Fd(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) - hpx(i)*Fc(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) - hhx(i)*Fb(j)*Jacobian*gw(igp)
           ENDDO
       ENDDO
   ENDDO

   DEALLOCATE(gp)
   DEALLOCATE(gw)
   DEALLOCATE(hhx)
   DEALLOCATE(hpx)
   DEALLOCATE(GLL_temp)
   DEALLOCATE(w_temp)

!   j=6
!   DO i=1,18
!       WRITE(*,*) elm(i,j+1), elm(i,j+2),elm(i,j+3),elm(i,j+4),elm(i,j+5),elm(i,j+6)
!       WRITE(*,*) elf(i)
!   ENDDO
!   STOP

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(gp))  DEALLOCATE(gp)
            IF(ALLOCATED(gw))  DEALLOCATE(gw)
            IF(ALLOCATED(hhx)) DEALLOCATE(hhx)
            IF(ALLOCATED(hpx)) DEALLOCATE(hpx)
            IF(ALLOCATED(GLL_temp)) DEALLOCATE(GLL_temp)
            IF(ALLOCATED(w_temp)) DEALLOCATE(w_temp)
        ENDIF

   END SUBROUTINE ElementMatrix
