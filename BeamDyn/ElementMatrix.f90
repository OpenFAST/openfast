   SUBROUTINE ElementMatrix(Nuu0,Nuuu,Nrr0,Nrrr,Next,hhp,Stif0,Jac,&
                            &w,node_elem,nelem,norder,dof_node,elk,elf)

   REAL(ReKi),INTENT(IN)::Nuu0(:),Nuuu(:),Nrr0(:),Nrrr(:),Next(:)
   REAL(ReKi),INTENT(IN)::hhp(:,:),Stif0(:,:,:),Jac
   REAL(ReKi),INTENT(IN)::w(:)
   INTEGER,INTENT(IN)::node_elem,nelem,norder,dof_node

   REAL(ReKi),INTENT(OUT)::elk(:,:),elf(:)      

   REAL(ReKi),ALLOCATABLE::Fc_elem(:,:),Fd_elem(:,:),Oe_elem(:,:,:)
   REAL(ReKi),ALLOCATABLE::Pe_elem(:,:,:),Qe_elem(:,:,:),Se_elem(:,:,:)
   REAL(ReKi)::E10(3),RR0(3,3),kapa(3),E1(3),Stif(6,6),cet
   REAL(ReKi)::Fc(6),Fd(6),Oe(6,6),Pe(6,6),Qe(6,6)

   INTEGER::i,j,k,m,n,temp_id1,temp_id2
   INTEGER::allo_stat
      

   ALLOCATE(Fc_elem(dof_node,node_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Fc_elem = 0.0D0
   
   ALLOCATE(Fd_elem(dof_node,node_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Fd_elem = 0.0D0
   
   ALLOCATE(Oe_elem(dof_node,dof_node,node_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Oe_elem = 0.0D0
   
   ALLOCATE(Pe_elem(dof_node,dof_node,node_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Pe_elem = 0.0D0
   
   ALLOCATE(Qe_elem(dof_node,dof_node,node_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Qe_elem = 0.0D0
   
   ALLOCATE(Se_elem(dof_node,dof_node,node_elem),STAT = allo_stat)
   IF(allo_stat/=0) GOTO 9999
   Se_elem = 0.0D0

   DO i=1,node_elem
       E10 = 0.0D0
       E1 = 0.0D0
       RR0 = 0.0D0
       kapa = 0.0D0
       Fc = 0.0D0
       Fd = 0.0D0
       Oe = 0.0D0
       Pe = 0.0D0
       Qe = 0.0D0
       Stif = 0.0D0
       cet = 0.0D0
       CALL NodalDataAt0(node_elem,nelem,norder,dof_node,i,hhp,Nuu0,E10)
       CALL NodalData(Nuuu,Nrrr,Nuu0,Nrr0,E10,hhp,Stif0,&
                      &node_elem,nelem,i,norder,dof_node,&
                      &E1,RR0,kapa,Stif,cet)
       CALL ElasticForce(E1,RR0,kapa,Stif,cet,Fc,Fd,Oe,Pe,Qe)
       Fc_elem(1:6,i) = Fc(1:6)
       Fd_elem(1:6,i) = Fd(1:6)
       Oe_elem(1:6,1:6,i) = Oe(1:6,1:6)
       Pe_elem(1:6,1:6,i) = Pe(1:6,1:6)
       Qe_elem(1:6,1:6,i) = Qe(1:6,1:6)
       Se_elem(1:6,1:6,i) = Stif(1:6,1:6)
   ENDDO

   DO i=1,node_elem
       DO j=1,6
           temp_id1 = (i-1)*dof_node+j
           DO k=1,6
               temp_id2 = (i-1)*dof_node+k
               elk(temp_id1,temp_id2) = w(i)*Qe_elem(j,k,i)*Jac
           ENDDO
       ENDDO
   ENDDO

   DO i=1,node_elem
       DO j=1,node_elem
           DO k=1,6
               temp_id1=(i-1)*dof_node+k
               DO m=1,6
                   temp_id2=(i-1)*dof_node+m
                   elk(temp_id1,temp_id2)=elk(temp_id1,temp_id2)+w(i)*Pe_elem(k,m,i)*hhp(j,i)
                   elk(temp_id1,temp_id2)=elk(temp_id1,temp_id2)+w(j)*Oe_elem(k,m,j)*hhp(i,j)
               ENDDO
           ENDDO
       ENDDO
   ENDDO

   DO i=1,node_elem
       DO j=1,node_elem
           DO k=1,6
               temp_id1=(i-1)*dof_node+k
               DO m=1,6
                   temp_id2=(j-1)*dof_node+m
                   DO n=1,node_elem
                       elk(temp_id1,temp_id2)=elk(temp_id1,temp_id2)+w(n)*hhp(i,n)*Se_elem(k,m,n)*hhp(j,n)/Jac
                   ENDDO
               ENDDO
           ENDDO
       ENDDO
   ENDDO

   DO i=1,node_elem
       DO j=1,6
           temp_id1 = (i-1)*dof_node+j
           elf(temp_id1) = -w(i)*Fd_elem(j,i)*Jac
           elf(temp_id1) = elf(temp_id1) + w(i) * Next((i-1)*dof_node+j) * Jac
           DO m=1,node_elem
              elf(temp_id1) = elf(temp_id1)-w(m)*hhp(i,m)*Fc_elem(m,j)
           ENDDO
       ENDDO
   ENDDO

   DEALLOCATE(Fc_elem)
   DEALLOCATE(Fd_elem)
   DEALLOCATE(Oe_elem)
   DEALLOCATE(Pe_elem)
   DEALLOCATE(Qe_elem)
   DEALLOCATE(Se_elem)
   
   9999 IF(allo_stat/=0) THEN
            IF(ALLOCATED(Fc_elem)) DEALLOCATE(Fc_elem)
            IF(ALLOCATED(Fd_elem)) DEALLOCATE(Fd_elem)
            IF(ALLOCATED(Oe_elem)) DEALLOCATE(Oe_elem)
            IF(ALLOCATED(Pe_elem)) DEALLOCATE(Pe_elem)
            IF(ALLOCATED(Qe_elem)) DEALLOCATE(Qe_elem)
            IF(ALLOCATED(Se_elem)) DEALLOCATE(Se_elem)
        ENDIF

   END SUBROUTINE ElementMatrix
