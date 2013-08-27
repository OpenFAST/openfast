   PROGRAM DriverStatic

   IMPLICIT NONE
   INTEGER norder, elem_total, node_total
   INTEGER dof_node,dof_total,node_elem
   INTEGER niter, ReKi

   PARAMETER (norder = 3)
   PARAMETER (elem_total = 1)
   PARAMETER (node_elem = norder + 1)
   PARAMETER (node_total = elem_total * norder + 1)
   PARAMETER (dof_node = 6)
   PARAMETER (dof_total = dof_node * node_total)
   PARAMETER (ReKi = SELECTED_REAL_KIND(15,307))

   REAL(ReKi),PARAMETER:: ONE = 1.d0
   REAL(ReKi),PARAMETER:: TWO = 2.d0
   REAL(ReKi),PARAMETER:: THREE = 3.d0
   REAL(ReKi),PARAMETER:: FOUR = 4.d0
   REAL(ReKi),PARAMETER:: FIVE = 5.d0
   REAL(ReKi),PARAMETER:: SIX = 6.d0
   REAL(ReKi),PARAMETER:: SEVEN = 7.d0
   REAL(ReKi),PARAMETER:: EIGHT = 8.d0
   REAL(ReKi),PARAMETER:: NINE = 9.d0

   REAL(ReKi),PARAMETER:: PI = 3.1415926535897932d0
   REAL(ReKi),PARAMETER:: DEG_2_RAD = 1.7453292519943296d-02
   REAL(ReKi),PARAMETER:: RAD_2_DEG = 5.7295779513082321d01
   REAL(ReKi),PARAMETER:: ZERO = 0.0d0
   REAL(ReKi),PARAMETER:: HALF = 0.5d0

   REAL(ReKi),PARAMETER:: I3(3,3) = RESHAPE((/1.d0, 0.d0, 0.d0,&
                                              &0.d0, 1.d0, 0.d0,&
                                              &0.d0, 0.d0, 1.d0/),&
                                              &(/3,3/))

   REAL(ReKi)::xj(node_elem),w(node_elem),hhp(node_elem,node_elem)
   REAL(ReKi)::uuN0(dof_total),uuNf(dof_total), Jacobian
   REAL(ReKi)::F_ext(dof_total),StifK(dof_total,dof_total),RHS(dof_total)
   REAL(ReKi)::Stif0(dof_node,dof_node,node_total)
   REAL(ReKi)::dloc(node_total), xmin, xmax, blength elem_length
   INTEGER::nquar   


   CALL gen_gll(norder,xj,w)
   CALL gen_deriv(xj,hhp,node_elem)
      
   xmin = 0.0d0
   xmax = 10.0d0
   blength = xmax - xmin
   elem_length = blength / float(elem_total)
   Jacobian = elem_length / 2.0d0
   niter = 300

   CALL NodeLoc(dloc, nquar, xmin, elem_length,& 
                & xj, norder, elem_total, node_total, blength)

   uuN0 = ZERO
   Stif0 = ZERO
   F_ext = 0.0d0
   DO i=1,node_total
       uuN0((i-1)*dof_node+i) = dloc(i)
       Stif0(1,1,i) = 70.0d+09
       Stif0(2,2,i) = 70.0d+09
       Stif0(3,3,i) = 70.0d+09
       Stif0(4,4,i) = 
       Stif0(5,5,i) = 5.8333d+09
       Stif0(6,6,i) = 5.8333d+09  
   ENDDO
   F_ext(dof_total - 1) = (5.8333d+08) * 1.0 * PI * 0.1

   uuNf = ZERO

   CALL StaticSolution(uuN0,uuNf,hhp,w,Jacobian,Stif0,F_ext,&
                       &node_elem,dof_node,norder,elem_total,dof_total,node_total,niter)

   WRITE(*,*) "Finished Static Solution!"

   DO i=1, node_total
       DO j=1, dof_node
           temp_id = (i-1)*dof_node + j
           WRITE(30,*) i, uuNf(temp_id)
       ENDDO
   ENDDO





   END PROGRAM DriverStatic
