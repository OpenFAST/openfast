subroutine Norm(dof_total,vector,norm_value)

   integer dof_total
   double precision vector(dof_total)
   double precision norm_value 

   norm_value=SQRT(DOT_PRODUCT(vector,vector))

   return
END subroutine
