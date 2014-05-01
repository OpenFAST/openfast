   FUNCTION CrossProduct(a,b)

   REAL(ReKi),INTENT(IN):: a(3),b(3)
   REAL(ReKi)           :: CrossProduct(3)
      
   CrossProduct(1) = a(2) * b(3) - a(3) * b(2)
   CrossProduct(2) = a(3) * b(1) - a(1) * b(3)
   CrossProduct(3) = a(1) * b(2) - a(2) * b(1)

   END FUNCTION CrossProduct
