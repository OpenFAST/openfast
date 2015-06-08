   FUNCTION Norm(vector)
   
   REAL(ReKi),INTENT(IN):: vector(:)
   REAL(ReKi)           :: Norm
   
   Norm = SQRT(DOT_PRODUCT(vector,vector))
   
   END FUNCTION Norm
