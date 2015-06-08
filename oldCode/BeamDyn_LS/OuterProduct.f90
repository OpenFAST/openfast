!   SUBROUTINE OuterProduct(vec1,vec2,outer)

!   REAL(ReKi),INTENT(IN):: vec1(:),vec2(:)
!   REAL(ReKi),INTENT(OUT)::outer(SIZE(vec1),SIZE(vec2))

!   INTEGER(IntKi)::i,j,n1,n2

!   n1=SIZE(vec1)
!   n2=SIZE(vec2)

!   outer = 0.0d0
!   DO i=1,n1
!       DO j=1,n2
!           outer(i,j) = vec1(i) * vec2(j)
!       ENDDO
!   ENDDO
      
!   END SUBROUTINE OuterProduct

   FUNCTION OuterProduct(vec1,vec2)

   REAL(ReKi),INTENT(IN):: vec1(:),vec2(:)
   REAL(ReKi)::OuterProduct(SIZE(vec1),SIZE(vec2))

   INTEGER(IntKi)::i,j,n1,n2

   n1=SIZE(vec1)
   n2=SIZE(vec2)

   DO i=1,n1
       DO j=1,n2
           OuterProduct(i,j) = vec1(i) * vec2(j)
       ENDDO
   ENDDO
      
   END FUNCTION OuterProduct
