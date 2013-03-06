subroutine AssembleKT(KT,dof_node,dof_total,norder,node_total,elem_total,&
                     &hhp,uf,dmat,wj,Jacobian)
   
   integer dof_node,dof_total,norder,node_total,elem_total
   double precision hhp(norder+1,norder+1)
   double precision uf(dof_total),dmat(node_total),wj(norder+1)
   double precision KT(dof_total,dof_total)
   double precision Jacobian
   
   integer nelem,m,temp_id,temp_id2
   double precision KT_elem(dof_node*(norder+1),dof_node*(norder+1))
   
   KT = 0.0d0  
      
   do nelem=1,elem_total
      KT_elem = 0.0d0
      
      call KT_element(KT_elem,dof_node,dof_total,norder,hhp,uf,nelem,&
                     &node_total,dmat,wj,Jacobian)
                     
      do i=1,(norder+1)*dof_node
         temp_id = (nelem-1)*norder*dof_node+i
         do j=1,(norder+1)*dof_node
            temp_id2 = (nelem-1)*norder*dof_node+j
            KT(temp_id,temp_id2) = KT(temp_id,temp_id2) + KT_elem(i,j)
         enddo
      enddo
   enddo
   
!  write(*,*) "dof_total",dof_total
!  do i=1,dof_total
!     write(*,*) KT(9,i)
!  enddo
!  stop
   
   
   
!  do nelem = 1, elem_total
!     call KTu_element(KTu_elem,dof_node,dof_total,norder,hhp,uf,ui,nelem,&
!                    &nmax,dmat,wj)
!     do m=1,(norder+1)*dof_node
!        temp_id = (nelem-1)*norder*dof_node+m
!        KTu(temp_id) += KTu_elem(m)
!     enddo    
!  enddo

   return
   
end subroutine
