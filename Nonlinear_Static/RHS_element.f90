	subroutine RHS_element(rhs_elem,dof_node,dof_total,&
								&uf,norder,hhp,wj,nelem,node_total,&
								&dmat)
	
	integer dof_node,norder,dof_total,nelem,node_total
	double precision rhs_elem(dof_node*(norder+1))
	double precision hhp(norder+1,norder+1),wj(norder+1),uf(dof_total)
	double precision dmat(node_total,3)
	
	integer nnode,m,temp_id
	double precision tempATN(4,1),tempDChi(3,1)
	double precision Am(3,4),Dm(3,3),Chim(3,1)
	
	
	rhs_elem = 0.0d0
	
	do nnode=1,norder+1
	
		call Amatrix(Am,dof_total,dof_node,norder,hhp,uf,nelem,nnode)
		call Dmatrix(Dm,node_total,dmat,nelem,nnode)
		call Chimatrix(Chim,dof_total,dof_node,norder,hhp,uf,nelem,nnode)
	
		tempDChi = MATMUL(Dm,Chim)
		tempATN = MATMUL(TRANSPOSE(Am),tempDChi)
	
		do m=1, norder+1
			temp_id = (m-1)*dof_node
			rhs_elem(temp_id+1) = hhp(nnode,m)*tempATN(1,1)
			rhs_elem(temp_id+2) = hhp(nnode,m)*tempATN(2,1)
			rhs_elem(temp_id+3) = hhp(nnode,m)*tempATN(3,1)
			if(nnode==m) rhs_elem(temp_id+3) = rhs_elem(temp_id+3) + tempATN(4,1)	
		enddo
	
		rhs_elem = wj(nnode)*rhs_elem
	
	enddo
	
	return
	
	end subroutine
	
	
	
	
	