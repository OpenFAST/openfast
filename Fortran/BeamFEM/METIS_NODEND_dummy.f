      SUBROUTINE METIS_NODEND(N,IPTR,IRN,METFTN,METOPT,INVPRM,PERM)
C Dummy routine that is called if MeTiS is not linked.
      INTEGER N
      INTEGER IPTR(N+1),IRN(*),METFTN,METOPT(8),INVPRM(N),PERM(N)
      PERM(1) = -1
      RETURN
      END

