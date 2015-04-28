   FUNCTION BD_Norm(vector)
   
   REAL(ReKi),INTENT(IN):: vector(:)
   REAL(ReKi)           :: BD_Norm
   
   BD_Norm = SQRT(DOT_PRODUCT(vector,vector))
   
   END FUNCTION BD_Norm
