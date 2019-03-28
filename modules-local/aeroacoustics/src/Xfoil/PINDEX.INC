C---------------------
C  Pointers for referencing polar force coefficients and parameters
C  First 4 pointers must be main polar plot variables.
C
      PARAMETER (
     &    IAL =  1,    ! alpha
     &    ICL =  2,    ! CL
     &    ICD =  3,    ! CD
     &    ICM =  4,    ! Cm
     &    ICW =  5,    ! CDwave
     &    ICV =  6,    ! CDvisc
     &    IMA =  7,    ! Mach
     &    IRE =  8,    ! Re
     &    INC =  9,    ! Ncrit
     &    ICH = 10,    ! Hinge moment
     &    IMC = 11 )   ! Minimum Cp on surface
      PARAMETER (IPTOT=11)
C
C---------------------
C  Pointers for referencing polar airfoil-side quantities
C
      PARAMETER (
     &    JTP = 1,     ! trip
     &    JTN = 2 )    ! transition
      PARAMETER (JPTOT=2)
