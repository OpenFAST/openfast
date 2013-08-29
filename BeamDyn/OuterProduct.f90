   SUBROUTINE OuterProduct(vec1,vec2,outer)

      DOUBLE PRECISION,INTENT(IN):: vec1(:),vec2(:)
      DOUBLE PRECISION,INTENT(INOUT)::outer(SIZE(vec1),SIZE(vec2))

      INTEGER::i,j,n1,n2

      n1=SIZE(vec1)
      n2=SIZE(vec2)

      outer = 0.0d0
      DO i=1,n1
          DO j=1,n2
              outer(i,j) = vec1(i) * vec2(j)
          ENDDO
      ENDDO
      
      END SUBROUTINE OuterProduct
