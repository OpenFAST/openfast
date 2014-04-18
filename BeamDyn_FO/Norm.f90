   SUBROUTINE Norm(dof_total,vector,norm_value)
   
   INTEGER(IntKi),INTENT(IN)::dof_total
   REAL(ReKi),INTENT(IN)::vector(:)
   REAL(ReKi),INTENT(OUT)::norm_value
   
   norm_value = SQRT(DOT_PRODUCT(vector,vector))
   
   END SUBROUTINE Norm
