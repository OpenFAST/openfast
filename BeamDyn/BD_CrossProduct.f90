   FUNCTION BD_CrossProduct(a,b)

   REAL(ReKi),INTENT(IN):: a(3),b(3)
   REAL(ReKi)           :: BD_CrossProduct(3)
      
   BD_CrossProduct(1) = a(2) * b(3) - a(3) * b(2)
   BD_CrossProduct(2) = a(3) * b(1) - a(1) * b(3)
   BD_CrossProduct(3) = a(1) * b(2) - a(2) * b(1)

   END FUNCTION BD_CrossProduct
