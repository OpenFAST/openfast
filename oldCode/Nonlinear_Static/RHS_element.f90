subroutine RHS_element(rhs_elem, dof_node, dof_total, &
                     & uf, norder, hhp, wj, nelem, node_total, &
                     & dmat,Jacobian)
   
   integer dof_node, norder, dof_total, nelem, node_total
   double precision rhs_elem(dof_node*(norder+1))
   double precision hhp(norder+1,norder+1), wj(norder+1), uf(dof_total)
   double precision dmat(node_total,3), Jacobian
   
   integer nnode, m, temp_id, i, j
   double precision tempATN(4,1), tempDChi(3,1)
   double precision Am(3,4), Dm(3,3), Chim(3,1)
   double precision rhs_temp(dof_node*(norder+1)), hhp_temp(norder+1,norder+1)
   
   
   rhs_elem = 0.0d0
   
   hhp_temp = TRANSPOSE(hhp)
   
   do nnode=1, norder+1
   
!   write(*,*) "dmat",dmat
   
      call Amatrix(Am, dof_total, dof_node, norder, hhp, uf, nelem, nnode, Jacobian)
      call Dmatrix(Dm, node_total, dmat, nelem, nnode)
      
      
      call Chimatrix(Chim, dof_total, dof_node, norder, hhp, uf, nelem, nnode, Jacobian)
         
      tempDChi = MATMUL(Dm, Chim)
      tempATN  = MATMUL(TRANSPOSE(Am), tempDChi)
   
      rhs_temp = 0.0d0

      do m=1, norder+1
         temp_id = (m-1)*dof_node
         rhs_temp(temp_id+1) = wj(nnode)*hhp_temp(nnode,m)*tempATN(1,1)
         rhs_temp(temp_id+2) = wj(nnode)*hhp_temp(nnode,m)*tempATN(2,1)
         rhs_temp(temp_id+3) = hhp_temp(nnode,m)*tempATN(3,1)
         if(nnode==m) rhs_temp(temp_id+3) = rhs_temp(temp_id+3) + tempATN(4,1)*Jacobian
         rhs_temp(temp_id+3) = wj(nnode)*rhs_temp(temp_id+3)   
      enddo
      
      rhs_elem = rhs_elem + rhs_temp
   
   enddo
   
   return
   
end subroutine
   
   
   
   
   
