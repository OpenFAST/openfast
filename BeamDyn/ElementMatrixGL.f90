   SUBROUTINE ElementMatrixGL(Nuu0,Nuuu,Nrr0,Nrrr,Next,Stif0,ngp,node_elem,dof_node,elk,elf)

   REAL(ReKi),INTENT(IN):: Nuu0(:),Nuuu(:),Nrr0(:),Nrrr(:),Next(:)
   REAL(ReKi),INTENT(IN):: Stif0(:,:)
   INTEGER(IntKi),INTENT(IN):: ngp,node_elem,dof_node

   REAL(ReKi),INTENT(INOUT):: elk(:,:),elf(:)      

   REAL(ReKi),ALLOCATABLE:: gp(:),gw(:),hhx(:),hpx(:)
   REAL(ReKi):: uu0(6),E10(3),RR0(3,3),kapa(3),E1(3),Stif(6,6),cet
   REAL(ReKi):: uuu(6),uup(3),Jacobian,gpr
   REAL(ReKi):: Fc(6),Fd(6),Oe(6,6),Pe(6,6),Qe(6,6)

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

   elk = 0.0D0
   elf = 0.0D0
   CALL BldSet1DGaussPointScheme(ngp,gp,gw)
   DO igp=1,ngp
       gpr = gp(igp)
       CALL BldComputeJacobian(gpr,Nuu0,node_elem,dof_node,hhx,hpx,Jacobian)
!       WRITE(*,*) "Jacobian = ", Jacobian
       CALL BldGaussPointDataAt0(hhx,hpx,Nuu0,Nrr0,node_elem,dof_node,uu0,E10)
       CALL BldGaussPointData(hhx,hpx,Nuuu,Nrrr,Stif0,uu0,E10,node_elem,dof_node,uuu,uup,E1,RR0,kapa,Stif,cet)
       CALL ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)

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
                   ENDDO
               ENDDO
           ENDDO
       ENDDO 

       DO i=1,node_elem
           DO j=1,dof_node
               temp_id1 = (i-1) * dof_node+j
               elf(temp_id1) = elf(temp_id1) - hhx(i)*Fd(j)*Jacobian*gw(igp)
               elf(temp_id1) = elf(temp_id1) - hpx(i)*Fc(j)*Jacobian*gw(igp)
           ENDDO
       ENDDO

   ENDDO

   DEALLOCATE(gp)
   DEALLOCATE(gw)
   DEALLOCATE(hhx)
   DEALLOCATE(hpx)

   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(gp))  DEALLOCATE(gp)
            IF(ALLOCATED(gw))  DEALLOCATE(gw)
            IF(ALLOCATED(hhx)) DEALLOCATE(hhx)
            IF(ALLOCATED(hpx)) DEALLOCATE(hpx)
        ENDIF

   END SUBROUTINE ElementMatrixGL
