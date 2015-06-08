subroutine AssembleRHS(RHS, dof_node, dof_total, uf, &
                     & norder, hhp, wj, node_total, dmat, &
                     & elem_total, Jacobian, F_ext)
   
   integer dof_node,dof_total,norder,node_total,elem_total
   double precision dmat(node_total,3)
   double precision hhp(norder+1,norder+1),wj(norder+1)
   double precision uf(dof_total)
   double precision RHS(dof_total)
   
   integer nelem, m, temp_id, i
   double precision rhs_elem(dof_node*(norder+1))
   double precision F_ext
!   double precision FmL
   
 !  FmL=6.28d+00

 !  FmL=0.01*(6.28d+01)
  
   RHS = 0.0d0
   
   do nelem = 1, elem_total
   
      call RHS_element(rhs_elem, dof_node, dof_total,&
                     & uf, norder, hhp, wj, nelem, node_total,&
                     & dmat, Jacobian)
                     
      

      do m=1, (norder+1)*dof_node
         temp_id      = (nelem-1)*norder*dof_node + m
         RHS(temp_id) = RHS(temp_id) + rhs_elem(m)
      enddo    

   enddo
   
!  do i=1,dof_total
!     write(*,*) RHS(i)
!  enddo
!  stop
   
   RHS = - RHS
   
!  do nelem = 1, elem_max
!     call Extvector_elem
!     do m=1,(norder+1)*dof_node
!        temp_id = (nelem-1)*norder*dof_node+m
!        RHS(temp_id) += ext_elem(m)
!     enddo
!  enddo

   RHS(dof_total-1) = RHS(dof_total-1) + F_ext 
   
   return
   
end subroutine
