   SUBROUTINE NodeLocGL(dloc,xmin,elem_length,node_elem,elem_total)

   INTEGER(IntKi), INTENT(IN):: node_elem,elem_total
   REAL(ReKi), INTENT(OUT):: dloc(:)

   REAL(ReKi), INTENT(IN):: xmin, elem_length

   INTEGER(IntKi):: i, j,npos, nelem


   dloc(1) = xmin

   npos = 1

   do nelem = 1, elem_total

      do j = 1, node_elem - 1

         npos = npos + 1

         dloc(npos) = (nelem - 1) * elem_length + elem_length*j/(node_elem-1)

      enddo

   enddo


   END SUBROUTINE NodeLocGL   
