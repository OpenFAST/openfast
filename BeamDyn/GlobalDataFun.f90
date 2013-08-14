      MODULE GlobalDataFun
 
      IMPLICIT NONE

      PRIVATE


      PUBLIC ZERO,ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT,NINE,HALF
      PUBLIC ReKi,PI,DEG_2_RAD,RAD_2_DEG,I3
      PUBLIC Tilde, CrossProduct,OuterProduct
 

      INTEGER,PARAMETER:: ReKi = SELECTED_REAL_KIND(15,307)

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
      
      CONTAINS

      FUNCTION Tilde(vect)

      REAL(ReKi),INTENT(IN):: vect(3)
      REAL(ReKi):: Tilde(3,3)

      Tilde = 0.d0

      Tilde(1,2) = -vect(3)
      Tilde(1,3) = vect(2)
      Tilde(2,1) = vect(3)
      Tilde(2,3) = -vect(1)
      Tilde(3,1) = -vect(2)
      Tilde(3,2) = vect(1)

      END FUNCTION Tilde

      FUNCTION CrossProduct(a,b)

      REAL(ReKi),INTENT(IN):: a(3),b(3)
      REAL(ReKi):: CrossProduct(3)
      
      CrossProduct(1) = a(2) * b(3) - a(3) * b(2)
      CrossProduct(2) = a(3) * b(1) - a(1) * b(3)
      CrossProduct(3) = a(1) * b(2) - a(3) * b(1)

      END FUNCTION CrossProduct

      FUNCTION OuterProduct(vec1,vec2)

      REAL(ReKi),INTENT(IN):: vec1(:),vec2(:)
      REAL(ReKi)::OuterProduct(SIZE(vec1),SIZE(vec2))

      INTEGER::i,j,n1,n2

      n1=SIZE(vec1)
      n2=SIZE(vec2)

      DO i=1,n1
          DO j=1,n2
              OuterProduct(i,j) = vec1(i) * vec2(j)
          ENDDO
      ENDDO
      
      END FUNCTION OuterProduct
