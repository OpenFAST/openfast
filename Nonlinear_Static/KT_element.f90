	subroutine KT_element(KT_elem,dof_node,dof_total,norder,hhp,uf,nelem,&
							&node_total,dmat,wj,Jacobian)
	
	integer norder,dof_node,dof_total,nelem,node_total
	double precision hhp(norder+1,norder+1),uf(dof_total)
	double precision dmat(node_total,3),wj(norder+1)
	double precision Jacobian
	
	double precision KT_elem(dof_node*(norder+1),dof_node*(norder+1))
	double precision Am(3,4),Dm(3,3),Xim(4,4)
	double precision tempKT(4,4),KT_temp(dof_node*(norder+1),dof_node*(norder+1))
	integer nnode,m,j,i
	integer temp_id,temp_id2
	
	KT_elem = 0.0d0
	
	do nnode=1,norder+1
	
		call Amatrix(Am,dof_total,dof_node,norder,hhp,uf,nelem,nnode)
		call Dmatrix(Dm,node_total,dmat,nelem,nnode)
		call Ximatrix(Xim,dof_total,dof_node,norder,hhp,uf,nelem,nnode,&
						&node_total,dmat)
						
		
			
		tempKT = MATMUL(MATMUL(TRANSPOSE(Am),Dm),Am)
		
			
		tempKT = tempKT+Xim
		
		KT_temp = 0.0d0
		
		do m=1, norder+1
			temp_id = (m-1)*dof_node
			do j=1,norder+1
				temp_id2 = (j-1)*dof_node
				KT_temp(temp_id+1,temp_id2+1) = hhp(nnode,m)*hhp(nnode,j)*tempKT(1,1)
				KT_temp(temp_id+1,temp_id2+2) = hhp(nnode,m)*hhp(nnode,j)*tempKT(1,2)
				KT_temp(temp_id+1,temp_id2+3) = hhp(nnode,m)*hhp(nnode,j)*tempKT(1,3)
				KT_temp(temp_id+2,temp_id2+1) = hhp(nnode,m)*hhp(nnode,j)*tempKT(2,1)
				KT_temp(temp_id+2,temp_id2+2) = hhp(nnode,m)*hhp(nnode,j)*tempKT(2,2)
				KT_temp(temp_id+2,temp_id2+3) = hhp(nnode,m)*hhp(nnode,j)*tempKT(2,3)
				KT_temp(temp_id+3,temp_id2+1) = hhp(nnode,m)*hhp(nnode,j)*tempKT(3,1)
				KT_temp(temp_id+3,temp_id2+2) = hhp(nnode,m)*hhp(nnode,j)*tempKT(3,2)
				KT_temp(temp_id+3,temp_id2+3) = hhp(nnode,m)*hhp(nnode,j)*tempKT(3,3)
				if(nnode==j) then 
					KT_temp(temp_id+1,temp_id2+3) =KT_temp(temp_id+1,temp_id2+3) + hhp(nnode,m)*tempKT(1,4)
					KT_temp(temp_id+2,temp_id2+3) =KT_temp(temp_id+2,temp_id2+3) + hhp(nnode,m)*tempKT(2,4)
					KT_temp(temp_id+3,temp_id2+3) =KT_temp(temp_id+3,temp_id2+3) + hhp(nnode,m)*tempKT(3,4)
					if(nnode==m) KT_temp(temp_id+3,temp_id2+3) =KT_temp(temp_id+3,temp_id2+3) + tempKT(4,4)
				endif
				if(nnode==m) then
					KT_temp(temp_id+3,temp_id2+1) =KT_temp(temp_id+3,temp_id2+1) + hhp(nnode,j)*tempKT(4,1)
					KT_temp(temp_id+3,temp_id2+2) =KT_temp(temp_id+3,temp_id2+2) + hhp(nnode,j)*tempKT(4,2)
					KT_temp(temp_id+3,temp_id2+3) =KT_temp(temp_id+3,temp_id2+3) + hhp(nnode,j)*tempKT(4,3)
				endif								
			enddo	
		enddo
		
		KT_elem = KT_elem + wj(nnode)*KT_temp/Jacobian
				
	enddo
	
!	write(*,*) KT_elem
!	stop
		
!	temp_id2 = (nelem-1)*norder+1  !Starting node number in a element
!	do nnode=1,norder+1
!		do m=1,dof_elem
!			temp_id = (nnode-1)*dof_node 
!			KTu_elem(temp_id+1) = KT_elem(temp_id+1,m)*ui((temp_id2-1)*dof_node+m)
!		enddo
!	enddo
	return	
	end subroutine
	
	
	