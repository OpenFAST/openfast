subroutine Dmatrix(Dm,node_total,dmat,nelem,nnode)

   integer nelem, nnode,node_total
   double precision dmat(node_total,3) 
   double precision Dm(3,3)

   integer temp_id

   Dm=0.0d0

   temp_id = (nelem-1)*norder+nnode
   Dm(1,1) = dmat(temp_id,1)
   Dm(2,2) = dmat(temp_id,2)
   Dm(3,3) = dmat(temp_id,3)

   return

end subroutine

