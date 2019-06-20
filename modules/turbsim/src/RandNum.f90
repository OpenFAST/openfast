!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2014, 2016  National Renewable Energy Laboratory
!
!    This file is part of TurbSim.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
MODULE TS_RandNum

   USE TurbSim_Types   
   USE Ran_Lux_Mod  

   IMPLICIT NONE


   INTEGER(IntKi), PARAMETER :: pRNG_RANLUX    = 1
   INTEGER(IntKi), PARAMETER :: pRNG_INTRINSIC = 2
   INTEGER(IntKi), PARAMETER :: pRNG_SNLW3     = 3   

   
   INTEGER, PARAMETER        :: LuxLevel       = 3       ! Luxury Level for RanLux RNG

   
CONTAINS

!=======================================================================
SUBROUTINE RandNum_Init(p, OtherSt, ErrStat, ErrMsg )

   ! Initialize the Random Number Generators


IMPLICIT                            NONE

TYPE(RandNum_ParameterType),  INTENT(IN   ) :: p           ! parameters for random number generation
TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt     ! other states for random number generation
INTEGER(IntKi)  ,             INTENT(OUT)   :: ErrStat     ! allocation status
CHARACTER(*) ,                INTENT(OUT)   :: ErrMsg      ! error message



REAL(ReKi)                       :: RN(1)
INTEGER                          :: I           ! loop counter
INTEGER                          :: NumSeeds    ! number of seeds in the intrinsic random number generator

ErrStat = ErrID_None
ErrMsg  = ""

IF (p%pRNG == pRNG_INTRINSIC) THEN ! RNG_type == 'NORMAL'


      ! determine the number of seeds necessary (gfortran needs 8 or 12 seeds, not just 2)

   CALL RANDOM_SEED ( SIZE = NumSeeds )

   IF ( NumSeeds /= 2 ) THEN
      CALL ProgWarn( ' The random number generator in use differs from the original code provided by NREL. This pRNG uses ' &
                        //TRIM(Int2LStr(NumSeeds))//' seeds instead of the 2 in the TurbSim input file.')
   END IF

   IF ( .NOT. ALLOCATED( OtherSt%nextSeed ) ) THEN
      CALL AllocAry( OtherSt%nextSeed, NumSeeds, 'nextSeed', ErrSTat, ErrMsg )
      IF (ErrStat >= AbortErrLev) RETURN
   END IF


         ! We'll just populate this with odd seeds = Seed(1) and even seeds = Seed(2)
   DO I = 1,NumSeeds,2
      OtherSt%nextSeed(I) = p%RandSeed(1)
   END DO
   DO I = 2,NumSeeds,2
      OtherSt%nextSeed(I) = p%RandSeed(2)
   END DO


   CALL RANDOM_SEED ( PUT=OtherSt%nextSeed )


ELSEIF (p%pRNG == pRNG_RANLUX) THEN ! RNG_type == 'RANLUX'

   CALL RLuxGo ( LuxLevel, ABS( p%RandSeed(1) ), 0, 0 )
   
   IF (.NOT. ALLOCATED( OtherSt%nextSeed ) ) THEN
      CALL AllocAry( OtherSt%nextSeed, 2, 'nextSeed', ErrStat, ErrMsg )
      IF (ErrStat >= AbortErrLev) RETURN 
   END IF
   
   

ELSE ! pRNG == pRNG_SNLW3

   
   IF (.NOT. ALLOCATED( OtherSt%nextSeed ) ) THEN
      CALL AllocAry( OtherSt%nextSeed, 3, 'nextSeed', ErrStat, ErrMsg )
      IF (ErrStat >= AbortErrLev) RETURN 
   END IF
   
   
      ! A quick and dirty way to get three random seeds for u, v, and w
      ! This implementation allows comparisons with Neil's SNLWIND-3D

   OtherSt%nextSeed = p%RandSeed  
      
   CALL ARand( OtherSt%nextSeed(1), RN,1,1)

   OtherSt%nextSeed(2) = OtherSt%nextSeed(1)+1
   CALL ARand( OtherSt%nextSeed(2), RN,1,1)

   OtherSt%nextSeed(3) = OtherSt%nextSeed(2)+1
   CALL ARand( OtherSt%nextSeed(3), RN,1,1)
      
   
ENDIF

END SUBROUTINE RandNum_Init
!=======================================================================
SUBROUTINE ARand(ix, RandNum_Ary,I, RNG_start)

IMPLICIT                     NONE

      ! Passed variables


REAL(ReKi),     INTENT(OUT)     :: RandNum_Ary(:)    ! Output: random numbers
INTEGER(IntKi), INTENT(IN)      :: I                 ! Input: Size of RandNum_Ary to be filled
INTEGER(IntKi), INTENT(INOUT)   :: ix                ! Input/Output: Seed  !BONNIE: should this be set to Integer(4), not default integer size?
INTEGER(IntKi), INTENT(IN)      :: RNG_start

      ! Local variables

INTEGER,    PARAMETER    :: B15  = 32768      ! = 2^15
INTEGER,    PARAMETER    :: B16  = 65536      ! = 2^16
INTEGER                  :: fHI
INTEGER                  :: K                 ! Loop counter
INTEGER                  :: leftLO
INTEGER,    PARAMETER    :: ranA = 16807      ! = 7^5
INTEGER                  :: rank
INTEGER,    PARAMETER    :: ranP = 2147483647 ! = 2^31 - 1 = huge(ranP)
INTEGER                  :: xHI
INTEGER                  :: xaLO

REAL(ReKi), PARAMETER    :: ranPR = 1.0 / REAL(ranP, ReKi)

!BONNIE: We should check that RandNum_Ary is dimensioned correctly....

DO K = RNG_start, RNG_start+I-1

   xHI    = ix / B16
   xaLO   = (ix - B16*xHI) * ranA     ! MOD( ix, B16 ) * ranA

   leftLO = xaLO / B16
   fHI    = xHI*ranA + leftLO
   rank   = fHI / B15

   ix     = (((xaLO - leftLO*B16) - ranP) + (fHI - rank*B15)*B16) + rank
            ! ranP is subtracted in order to avoid overflow
            ! MOD( xaLO, B16 ) - ranP + B16*MOD( fHI, B15 ) + (fHI / B15)

   IF (ix < 0) ix = ix + ranP

   RandNum_Ary(K) = REAL(ix) * ranPR

ENDDO

END SUBROUTINE ARand
!=======================================================================
SUBROUTINE Rnd4ParmLogNorm( p, OtherSt, RandNum, Parms, FnRange )
   ! This subroutine generates a random variate with a PDF of the form
   ! f(x) = A + B*exp(-0.5*(ln(x/C)/D)^2)
   ! a truely log-normal distribution has A = 0

IMPLICIT                            NONE
TYPE(RandNum_ParameterType),  INTENT(IN   ) :: p                   ! parameters for random number generation
TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt             ! other states for random number generation

REAL(ReKi), INTENT(IN)           :: Parms(4)            ! 1=a, 2=b, 3=c, 4=d
REAL(ReKi), INTENT(IN)           :: FnRange(2)

REAL(ReKi)                       :: fMAX                ! Max(f(x)) ! occurs at f(b)
REAL(ReKi)                       :: Gx                  ! The function g(x) = f(x)/fMAX
REAL(ReKi)                       :: MaxVALUE            ! Maximum value of returned variate
REAL(ReKi)                       :: MinVALUE            ! Minimum value of returned variate
REAL(ReKi), INTENT(OUT)          :: RandNum             ! numbers distributed with the pdf above
REAL(ReKi)                       :: RN                  ! A random number for the acceptance-rejection method

INTEGER                          :: Count
INTEGER, PARAMETER               :: MaxIter = 10000     ! Max number of iterations to converge (so we don't get an infinite loop)

fMAX     = Parms(1) + Parms(2) !See if this is correct.
MaxVALUE = MAX(FnRange(1),FnRange(2))
MinVALUE = MIN(FnRange(1),FnRange(2))

Count = 1

      ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )

DO WHILE (Count < MaxIter)

   CALL RndUnif( p, OtherSt, RN )        ! Generate RN from U(0,1)

   CALL RndUnif( p, OtherSt, RandNum )   ! Generate RandNum from h(y) = 1 / (MaxVALUE - MinVALUE)
   RandNum = RandNum*(MaxVALUE-MinVALUE)+MinVALUE;

   Gx = Parms(1) + Parms(2)*EXP(-0.5*(LOG(RandNum/Parms(3))/Parms(4))**2);
   Gx = Gx / fMAX

   IF ( RN <= Gx ) THEN
      Count = MaxIter      ! Let's keep this deviate
   ELSE
      Count = Count + 1    ! try again
      RandNum = -1
   ENDIF

ENDDO

END SUBROUTINE Rnd4ParmLogNorm
!=======================================================================
SUBROUTINE Rnd3ParmNorm( p, OtherSt, RandNum, A, B, C, xMin, xMax, ErrStat, ErrMsg )

! Calculates a deviate from a distribution with pdf:
! f(x) = A * EXP( -0.5 * ((x-B)/C)**2 )
! A 3-parameter normal distribution
! We assume the returned values are between -1 and 1, since this is for the Cross-Correlations, unless
! the optional values, xMin and xMax are used

   IMPLICIT                            NONE

   TYPE(RandNum_ParameterType),  INTENT(IN   ) :: p                   ! parameters for random number generation
   TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt             ! other states for random number generation

   REAL(ReKi),                   INTENT(IN)    :: A
   REAL(ReKi),                   INTENT(IN)    :: B
   REAL(ReKi),                   INTENT(IN)    :: C

   INTEGER(IntKi),                  intent(  out) :: ErrStat                         ! Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg                          ! Message describing error



REAL(ReKi)                                  :: fMAX                ! Max(f(x)) = f(B) = A
REAL(ReKi)                                  :: Gx                  ! The function g(x) = f(x)/fMAX
REAL(ReKi)                                  :: MaxVALUE            ! Maximum value of returned variate
REAL(ReKi)                                  :: MinVALUE            ! Minimum value of returned variate
REAL(ReKi),                   INTENT(OUT)   :: RandNum             ! numbers distributed with the pdf above
REAL(ReKi)                                  :: RN                  ! A random number for the acceptance-rejection method
REAL(ReKi),    OPTIONAL,      INTENT(IN)    :: xMax                ! The maximum returned iterate
REAL(ReKi),    OPTIONAL,      INTENT(IN)    :: xMin                ! The minimum returned iterate

INTEGER                                     :: Count
INTEGER, PARAMETER                          :: MaxIter = 10000     ! Max number of iterations to converge (so we don't get an infinite loop)


   ErrStat = ErrID_None
   ErrMsg  = ""

! If A < 0 then we have a minimum value in the center of the distribution, not a maximum -- this method won't work.
! If A < 1/(MaxVALUE-MinVALUE), then this acceptance-rejection method won't work.

IF ( PRESENT(xMax) ) THEN
   MaxVALUE = xMax
ELSE
   MaxVALUE = 1.0
ENDIF

IF ( PRESENT(xMin) ) THEN
   MinVALUE = xMin
ELSE
   MinVALUE = -MaxVALUE
ENDIF

RN = 1. / (MaxVALUE-MinVALUE)
IF (A < RN .OR. C==0. .OR. MaxVALUE <= MinVALUE) THEN
   ErrStat = ErrID_Fatal
   ErrMsg  = 'Rnd3ParmNorm: Parameter A must at least 1/(xMax-xMin) and parameter C cannot be zero in this 3-parameter normal distribution.'
   RETURN
ENDIF


fMAX  = A
Count = 1

      ! Generate a 3-parameter normal distribution on (-1,1) from a uniform distribution

DO WHILE (Count < MaxIter)

   CALL RndUnif( p, OtherSt, RN )        ! Generate RN from U(0,1)
   CALL RndUnif( p, OtherSt, RandNum )   ! Generate RandNum from h(y) = 1 / (MaxVALUE - MinVALUE)
   RandNum = RandNum*(MaxVALUE-MinVALUE)+MinVALUE;

   Gx = A * EXP( -0.5 * ((RandNum-B)/C)**2 )
   Gx = Gx / fMAX

   IF ( RN <= Gx ) THEN
      Count = MaxIter      ! Let's keep this deviate
   ELSE
      Count = Count + 1    ! try again
      RandNum = -100
   ENDIF

ENDDO

END SUBROUTINE Rnd3ParmNorm
!=======================================================================
SUBROUTINE RndExp( p, OtherSt, RandExpNum, mu )

   ! This subroutine computes an exponential distribution on (0,inf). If the
   ! number of random variates to return is large, a different algorithm will
   ! probably be faster (i.e. one that computes LOG(x) fewer times).
   ! mu must be positive, it defaults to 1.0 if the parameter is not included.
   ! RandNum has p.d.f.(x) = 1/mu * exp(-x/mu), x>=0
   ! The expected value of RandNum is mu.


   IMPLICIT                NONE

      ! Passed Variables

TYPE(RandNum_ParameterType),  INTENT(IN   ) :: p                   ! parameters for random number generation
TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt             ! other states for random number generation

REAL(ReKi),                   INTENT(OUT)   :: RandExpNum          ! The exponentially distributed numbers in (0,1)
REAL(ReKi), OPTIONAL,         INTENT(IN)    :: mu                  ! The exponential distribution parameter equal to the expected value of the distribution

      ! Local Variable

REAL(ReKi)                                  :: mu_use


IF ( PRESENT(mu) ) THEN
   IF (mu < 0.0) THEN
      CALL WrScr( 'RndExp: Parameter mu is negative. Using -mu for exponential distribution.')
   ENDIF
   mu_use = ABS( mu )
ELSE
   mu_use = 1.0
ENDIF


   ! Get a uniform distribution of random numbers

CALL RndUnif( p, OtherSt, RandExpNum )

IF ( RandExpNum == 0.0 )  THEN ! We shouldn't get two zeros in a row...
   CALL RndUnif( p, OtherSt, RandExpNum )
ENDIF

   ! Transform the uniform distribution to an exponential distribution

RandExpNum = - mu_use * LOG( RandExpNum )

END SUBROUTINE RndExp

!=======================================================================
SUBROUTINE RndJetHeight( p, OtherSt, RandNum )
! This function uses the Pearson IV equation

IMPLICIT                            NONE

TYPE(RandNum_ParameterType),  INTENT(IN   )  :: p                   ! parameters for random number generation
TYPE(RandNum_OtherStateType), INTENT(INOUT)  :: OtherSt             ! other states for random number generation

REAL(ReKi), PARAMETER                        :: a =   0.021548497
REAL(ReKi), PARAMETER                        :: b = -13.173289
REAL(ReKi), PARAMETER                        :: c =  13.43201034
REAL(ReKi), PARAMETER                        :: d =   0.896588964
REAL(ReKi), PARAMETER                        :: e =  -0.71128456
                                        
REAL(ReKi), PARAMETER                        :: MaxVALUE =  120     ! Maximum value of returned variate
REAL(ReKi), PARAMETER                        :: MinVALUE = -160     ! Minimum value of returned variate
REAL(ReKi), PARAMETER                        :: Parms(5) = (/ a, b, c, d, e /)
REAL(ReKi), INTENT(OUT)                      :: RandNum             ! numbers distributed with the pdf above
REAL(ReKi), PARAMETER                        :: RangeFn(2)  = (/ MinVALUE, MaxVALUE /)


   CALL RndPearsonIV( p, OtherSt, RandNum, Parms, RangeFn )


END SUBROUTINE RndJetHeight
!=======================================================================
SUBROUTINE RndModLogNorm( p, OtherSt, RandNum, Height )
   ! This subroutine generates a random variate with a PDF of the form
   ! f(x) = A + B*exp(-0.5*(ln(x/C)/D)^2)

!BJJ use Rnd4ParmLogNorm()
IMPLICIT                            NONE

TYPE(RandNum_ParameterType),  INTENT(IN   )  :: p              ! parameters for random number generation
TYPE(RandNum_OtherStateType), INTENT(INOUT)  :: OtherSt        ! other states for random number generation

REAL(ReKi),                   INTENT(OUT)    :: RandNum        ! Near-Log-Normally distributed numbers
REAL(ReKi), OPTIONAL,         INTENT(IN)     :: Height         ! height (in meters), determining what parameters to use

   ! Internal variables

REAL(ReKi), PARAMETER                        :: A(3) = (/-0.0041046,   -0.00566512,  -0.00216964  /)
REAL(ReKi), PARAMETER                        :: B(3) = (/ 0.162945643,  0.278246235,  0.113718973 /)
REAL(ReKi), PARAMETER                        :: C(3) = (/ 0.67493672,   0.203262077,  3.211606394 /)
REAL(ReKi), PARAMETER                        :: D(3) = (/ 2.391316782,  2.715789776,  1.700298642 /)
REAL(ReKi)                                   :: G              ! The function g(x) = f(x)/B
REAL(ReKi)                                   :: RN (2)         ! Two random numbers for the acceptance-rejection method


INTEGER                                      :: Count
INTEGER                                      :: Indx
INTEGER, PARAMETER                           :: MaxIter = 10000  ! Max number of iterations to converge (so we don't get an infinite loop)
INTEGER, PARAMETER                           :: MaxTime = 600    ! Maximum value of returned value (the data used to compute A,B,C,D is valid up to 600 s.)

   !Get the index closest to the station we want to use...
Indx = 2 ! Index 2 == NC-CC-SC stations (37 m)
IF ( PRESENT(Height) ) THEN
   IF (Height > 47) THEN
      Indx = 1    ! Index 1 == UC station (58 m)
   ELSEIF (Height < 26) THEN
      Indx = 3 ! Index 3 = LC station (15 m)
   ENDIF
ENDIF

Count = 1

      ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )

DO WHILE (Count < MaxIter)

   CALL RndUnif( p, OtherSt, RN(1) )
   CALL RndUnif( p, OtherSt, RN(2) )

   RandNum = RN(2)*MaxTime;

   g = A(Indx)/B(Indx) + EXP(-0.5*(LOG(RandNum/C(Indx))/D(Indx))**2);

   IF ( RN(1) <= g ) THEN
      Count = MaxIter      ! Let's keep this deviate
   ELSE
      Count = Count + 1    ! try again
      RandNum = -1
   ENDIF

ENDDO


END SUBROUTINE RndModLogNorm
!=======================================================================
SUBROUTINE RndNorm( p, OtherSt, RandNormNum, mu, sigma )

IMPLICIT                            NONE

TYPE(RandNum_ParameterType),  INTENT(IN   ) :: p                   ! parameters for random number generation
TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt             ! other states for random number generation

REAL(ReKi),                   INTENT(OUT)   :: RandNormNum         ! Normally distributed numbers
REAL(ReKi), OPTIONAL,         INTENT(IN)    :: mu                  ! mean of the distributed numbers - DEFAULT IS 0.0
REAL(ReKi), OPTIONAL,         INTENT(IN)    :: sigma               ! standard deviation of the distributed numbers - DEFAULT IS 1.0

   ! Internal variable

REAL(ReKi)                                  :: RN (2)              ! Two random numbers


      ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )

   CALL RndUnif( p, OtherSt, RN(1) )
   CALL RndUnif( p, OtherSt, RN(2) )

   RandNormNum = SQRT( PI / 8.0 ) * LOG( ( 1.0 + RN(1) ) / ( 1.0 - RN(1) ) )

   IF ( RN(2) < 0.5 ) THEN
      RandNormNum = -RandNormNum
   ENDIF


      ! Give the correct mean and standard deviation, if specified

   IF ( PRESENT( sigma ) ) THEN
      RandNormNum = RandNormNum * sigma
   ENDIF

   IF ( PRESENT( mu ) ) THEN
      RandNormNum = RandNormNum + mu
   ENDIF

END SUBROUTINE RndNorm
!=======================================================================
SUBROUTINE RndNWTCpkCTKE( p, OtherSt, RandNum )
   ! This subroutine generates a random variate with a PDF of the form
   ! f(x) = A + B * EXP( (-X + C + D - D*E*EXP(-( X + D*LOG(E) - C)/D)) / (D*E)
   ! Maximum, f(C) = A + B
   ! Uses the Acceptance-Rejection: f(x) = Cf*h(x)*g(x)
   !   where h(x) = 1 / (150-30) (for our domain)
   !         g(x) = f(x)/(A+B), and
   !         Cf   = (150-30)*(A+B)

IMPLICIT                            NONE

TYPE(RandNum_ParameterType),  INTENT(IN   ) :: p                   ! parameters for random number generation
TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt             ! other states for random number generation

REAL(ReKi),                   INTENT(OUT)   :: RandNum             ! numbers distributed with the pdf above

   ! Internal variables

REAL(ReKi), PARAMETER                       :: A =   0.000500609
REAL(ReKi), PARAMETER                       :: B =   0.286202317
REAL(ReKi), PARAMETER                       :: C = -38.4131676
REAL(ReKi), PARAMETER                       :: D = 244.6908697
REAL(ReKi), PARAMETER                       :: E =   0.02115063

REAL(ReKi), PARAMETER                       :: fMAX = A + B     ! Max(f(x))
REAL(ReKi)                                  :: Gx               ! The function g(x) = f(x)/B
REAL(ReKi), PARAMETER                       :: MaxVALUE = 150.0 ! Maximum value of returned variate
REAL(ReKi), PARAMETER                       :: MinVALUE =  30.0 ! Minimum value of returned variate
REAL(ReKi)                                  :: RN               ! A random number for the acceptance-rejection method
                                            
INTEGER                                     :: Count
INTEGER, PARAMETER                          :: MaxIter = 10000 ! Max number of iterations to converge (so we don't get an infinite loop)

Count = 1

      ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )

DO WHILE (Count < MaxIter)

   CALL RndUnif( p, OtherSt, RN )        ! Generate RN from U(0,1)

   CALL RndUnif( p, OtherSt, RandNum )   ! Generate RandNum from h(y) = 1 / (MaxVALUE - MinVALUE)
   RandNum = RandNum*(MaxVALUE-MinVALUE)+MinVALUE;

   Gx = A + B * EXP( (-RandNum + C + D - D*E*EXP(-( RandNum + D*LOG(E) - C)/D)) / (D*E) )
   Gx = Gx / fMAX

   IF ( RN <= Gx ) THEN
      Count = MaxIter      ! Let's keep this deviate
   ELSE
      Count = Count + 1    ! try again
      RandNum = -1
   ENDIF

ENDDO

END SUBROUTINE RndNWTCpkCTKE
!=======================================================================
SUBROUTINE RndNWTCuStar( p, OtherSt, RandNum ) 
!bjj 17-jul-2014: this isn't used anywhere....
   ! This subroutine generates a random variate with a PDF of the form
   ! f(x) = (A + Cx + Ex^2 + Gx^3) / (1 + Bx + Dx^2 + Fx^3 + Hx^4)
   ! using the acceptance/rejection method.

IMPLICIT                            NONE

TYPE(RandNum_ParameterType),  INTENT(IN   ) :: p                   ! parameters for random number generation
TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt             ! other states for random number generation

REAL(ReKi), INTENT(OUT)                     :: RandNum             ! numbers distributed with the pdf above

   ! Internal variables

REAL(ReKi), PARAMETER                       :: A =   4.50581       ! Scaling parameters for the pdf
REAL(ReKi), PARAMETER                       :: B =  -0.60722       ! Scaling parameters for the pdf
REAL(ReKi), PARAMETER                       :: C = -14.23826       ! Scaling parameters for the pdf
REAL(ReKi), PARAMETER                       :: D =  -0.96523       ! Scaling parameters for the pdf
REAL(ReKi), PARAMETER                       :: E =  15.92342       ! Scaling parameters for the pdf
REAL(ReKi), PARAMETER                       :: F =  14.41326       ! Scaling parameters for the pdf
REAL(ReKi), PARAMETER                       :: G =  -6.16188       ! Scaling parameters for the pdf
REAL(ReKi), PARAMETER                       :: H =  -4.82923       ! Scaling parameters for the pdf

REAL(DbKi)                                  :: Gx                  ! The function g(x) = f(x)/A
REAL(ReKi), PARAMETER                       :: MaxUstar = 1.0      ! Maximum value of returned value (the data used to compute A,B,C,D is valid up to 600 s.)
REAL(DbKi)                                  :: RandNum2            ! RandNum**2
REAL(DbKi)                                  :: RandNum3            ! RandNum**3
REAL(DbKi)                                  :: RandNum4            ! RandNum**4
REAL(ReKi)                                  :: RN                  ! A random number for the acceptance-rejection method
                                            
INTEGER                                     :: Count
INTEGER, PARAMETER                          :: MaxIter = 10000     ! Max number of iterations to converge (so we don't get an infinite loop)

Count = 1

      ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )

DO WHILE (Count < MaxIter)

   CALL RndUnif( p, OtherSt, RN )
   CALL RndUnif( p, OtherSt, RandNum )

   RandNum = RandNum*MaxUstar;

   RandNum2 = RandNum*RandNum
   RandNum3 = RandNum*RandNum2
   RandNum4 = RandNum*RandNum3

   Gx = (A + C*RandNum + E*RandNum2 + G*RandNum3) / &
        (1 + B*RandNum + D*RandNum2 + F*RandNum3  + H*RandNum4)
   Gx = Gx / A  ! This makes 0<Gx<=1

   IF ( RN <= Gx ) THEN
      Count = MaxIter      ! Let's keep this deviate
   ELSE
      Count = Count + 1    ! try again
      RandNum = -1
   ENDIF

ENDDO

END SUBROUTINE RndNWTCuStar
!=======================================================================
SUBROUTINE RndPearsonIV( p, OtherSt, RandNum, Parms, FnRange )
! This function uses the Pearson IV equation to generate a deviate from that distribution
! Equation 8186 in TableCurve

IMPLICIT                            NONE

TYPE(RandNum_ParameterType),  INTENT(IN   ) :: p                   ! parameters for random number generation
TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt             ! other states for random number generation

REAL(ReKi),                   INTENT(IN)    :: Parms(5)           ! 1=a, 2=b, 3=c, 4=d, 5=e
REAL(ReKi),                   INTENT(IN)    :: FnRange(2)
                                            
REAL(ReKi)                                  :: fMAX                ! Max(f(x)) ! occurs at f(b)
REAL(ReKi)                                  :: Gx                  ! The function g(x) = f(x)/fMAX
REAL(ReKi)                                  :: MaxVALUE            ! Maximum value of returned variate
REAL(ReKi)                                  :: MinVALUE            ! Minimum value of returned variate
REAL(ReKi)                                  :: n                   ! A temporary variable for calculating the function values
REAL(ReKi),                   INTENT(OUT)   :: RandNum             ! numbers distributed with the pdf above
REAL(ReKi)                                  :: RN                  ! A random number for the acceptance-rejection method
                                            
INTEGER                                     :: Count
INTEGER, PARAMETER                          :: MaxIter = 10000     ! Max number of iterations to converge (so we don't get an infinite loop)

fMAX     = Parms(1)
MaxVALUE = MAX(FnRange(1),FnRange(2))
MinVALUE = MIN(FnRange(1),FnRange(2))

Count = 1

      ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )

DO WHILE (Count < MaxIter)

   CALL RndUnif( p, OtherSt, RN )        ! Generate RN from U(0,1)

   CALL RndUnif( p, OtherSt, RandNum )   ! Generate RandNum from h(y) = 1 / (MaxVALUE - MinVALUE)
   RandNum = RandNum*(MaxVALUE-MinVALUE)+MinVALUE;

   n  = (RandNum - Parms(2))/Parms(3) - Parms(5)/(2.0*Parms(4))
   Gx = Parms(1) * (1+n**2)**(-Parms(4)) * EXP( -Parms(5)*(ATAN(n)+ATAN(Parms(5)/(2*Parms(4)))) ) / &
                                                  ( 1+Parms(5)**2/(4.*Parms(4)**2))**(-Parms(4) )
   Gx = Gx / fMAX

   IF ( RN <= Gx ) THEN
      Count = MaxIter      ! Let's keep this deviate
   ELSE
      Count = Count + 1    ! try again
      RandNum = -1
   ENDIF

ENDDO


END SUBROUTINE RndPearsonIV
!=======================================================================
SUBROUTINE RndpkCTKE_WFTA( p, OtherSt, RandNum )
! Calculates a deviate from a distribution with pdf:
! f(x) = (a + c*x^2 + e*x^4)/(1. + b*x^2 + d*x^4 + f*x^6), where

IMPLICIT                            NONE

TYPE(RandNum_ParameterType),  INTENT(IN   ) :: p                   ! parameters for random number generation
TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt             ! other states for random number generation

!67m and 85m
REAL(ReKi), PARAMETER                       :: A = 0.30985506
REAL(ReKi), PARAMETER                       :: B = 0.006902104
REAL(ReKi), PARAMETER                       :: C =-0.00206008
REAL(ReKi), PARAMETER                       :: D = 1.28884E-05
REAL(ReKi), PARAMETER                       :: E = 5.71475E-06
REAL(ReKi), PARAMETER                       :: F = 8.70606E-07

REAL(ReKi)                                  :: fMAX                ! Max(f(x)) = f(MinVALUE)
REAL(ReKi)                                  :: Gx                  ! The function g(x) = f(x)/fMAX
REAL(ReKi), PARAMETER                       :: MaxVALUE = 22.      ! Maximum value of returned variate
REAL(ReKi), PARAMETER                       :: MinVALUE =  4.      ! Minimum value of returned variate
REAL(ReKi),                   INTENT(OUT)   :: RandNum             ! numbers distributed with the pdf above
REAL(ReKi)                                  :: RN                  ! A random number for the acceptance-rejection method
REAL(ReKi)                                  :: x2                  ! A temp variable for x^2
                                            
INTEGER                                     :: Count
INTEGER,    PARAMETER                       :: MaxIter = 10000     ! Max number of iterations to converge (so we don't get an infinite loop)


x2   = MinVALUE**2
fMAX = (A + C*x2 + E*x2**2)/(1.0 + B*x2 + D*x2**2 + F*x2**3)

Count = 1

      ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )

DO WHILE (Count < MaxIter)

   CALL RndUnif( p, OtherSt, RN )        ! Generate RN from U(0,1)
   CALL RndUnif( p, OtherSt, RandNum )   ! Generate RandNum from h(y) = 1 / (MaxVALUE - MinVALUE)
   RandNum = RandNum*(MaxVALUE-MinVALUE)+MinVALUE;

   x2      = RandNum**2
   Gx      = (A + C*x2 + E*x2**2)/(1.0 + B*x2 + D*x2**2 + F*x2**3)
   Gx      = Gx / fMAX

   IF ( RN <= Gx ) THEN
      Count = MaxIter      ! Let's keep this deviate
   ELSE
      Count = Count + 1    ! try again
      RandNum = -1
   ENDIF

ENDDO

END SUBROUTINE RndpkCTKE_WFTA
!=======================================================================
SUBROUTINE RndPolyFit( p, OtherSt, RandNum, Coeffs, FnRange, fMAX)
!bjj: 15-jul-2014: no longer used
! Calculates a deviate from a distribution with pdf:
! f(x) = (Coeffs(1) + Coeffs(3)*x + Coeffs(5)*x^2 + Coeffs(7)*x^3 + Coeffs(9)*x^4 + Coeffs(11)*x^5) / &
!        (       1. + Coeffs(2)*x + Coeffs(4)*x^2 + Coeffs(6)*x^3 + Coeffs(8)*x^4 + Coeffs(10)*x^5)
! This equation covers the following (plus others) from Table Curve\SysStat:
! Eqn 7906, Eqn 7907, Eqn 7908, & Eqn 7909


IMPLICIT                            NONE

TYPE(RandNum_ParameterType),  INTENT(IN   ) :: p                   ! parameters for random number generation
TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt             ! other states for random number generation


REAL(ReKi), INTENT(IN)           :: Coeffs(11)
REAL(ReKi), INTENT(IN)           :: FnRange(2)

REAL(ReKi), INTENT(IN)           :: fMAX                ! Max(f(x))
REAL(ReKi)                       :: Gx                  ! The function g(x) = f(x)/fMAX
REAL(ReKi)                       :: MaxVALUE            ! Maximum value of returned variate
REAL(ReKi)                       :: MinVALUE            ! Minimum value of returned variate
REAL(ReKi), INTENT(OUT)          :: RandNum             ! numbers distributed with the pdf above
REAL(ReKi)                       :: RN                  ! A random number for the acceptance-rejection method

INTEGER                          :: Count
INTEGER, PARAMETER               :: MaxIter = 10000     ! Max number of iterations to converge (so we don't get an infinite loop)

MaxVALUE = MAX(FnRange(1),FnRange(2))
MinVALUE = MIN(FnRange(1),FnRange(2))


Count = 1

      ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )

DO WHILE (Count < MaxIter)

   CALL RndUnif( p, OtherSt, RN )        ! Generate RN from U(0,1)
   CALL RndUnif( p, OtherSt, RandNum )   ! Generate RandNum from h(y) = 1 / (MaxVALUE - MinVALUE)
   RandNum = RandNum*(MaxVALUE-MinVALUE)+MinVALUE;

   Gx = (Coeffs(1) + Coeffs(3)*RandNum + Coeffs(5)*RandNum**2 + Coeffs( 7)*RandNum**3 + &
                                         Coeffs(9)*RandNum**4 + Coeffs(11)*RandNum**5 ) / &
        (       1. + Coeffs(2)*RandNum + Coeffs(4)*RandNum**2 + Coeffs( 6)*RandNum**3 + &
                                         Coeffs(8)*RandNum**4 + Coeffs(10)*RandNum**5 )
   Gx = Gx / fMAX

   IF ( RN <= Gx ) THEN
      Count = MaxIter      ! Let's keep this deviate
   ELSE
      Count = Count + 1    ! try again
      RandNum = -1
   ENDIF

ENDDO


END SUBROUTINE RndPolyFit
!=======================================================================
SUBROUTINE RndTcohLLJ( p, OtherSt, RandNum, Height )
! Calculates a deviate from a distribution with pdf:
! f(x) = EXP(A + B*SQRT(x) + C*LOG(x) )  !Eqn 1376

IMPLICIT                            NONE

TYPE(RandNum_ParameterType),  INTENT(IN   ) :: p                   ! parameters for random number generation
TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt             ! other states for random number generation
 
REAL(ReKi),                   INTENT(  OUT) :: RandNum             ! numbers distributed with the pdf above


!67m and 85m
REAL(ReKi), PARAMETER                       :: A(2) = (/ -1.34064396 , -1.17577736  /)
REAL(ReKi), PARAMETER                       :: B(2) = (/ -0.26996911 , -0.23056567  /)
REAL(ReKi), PARAMETER                       :: C(2) = (/ -0.57793906 , -0.69871145  /)
                                            
REAL(ReKi)                                  :: fMAX                ! Max(f(x)) = f(MinValue)
REAL(ReKi)                                  :: Gx                  ! The function g(x) = f(x)/fMAX
REAL(ReKi)                                  :: Height              ! The height of the center of the billow, in meters
REAL(ReKi), PARAMETER                       :: MaxVALUE = 600.     ! Maximum value of returned variate
REAL(ReKi), PARAMETER                       :: MinVALUE = 2.5      ! Minimum value of returned variate
REAL(ReKi)                                  :: RN                  ! A random number for the acceptance-rejection method
                                            
INTEGER                                     :: Count
INTEGER                                     :: Indx
INTEGER,    PARAMETER                       :: MaxIter = 10000     ! Max number of iterations to converge (so we don't get an infinite loop)

IF (Height < 76) THEN
   Indx = 1
ELSE
   Indx = 2
ENDIF

fMAX = EXP(A(Indx) + B(Indx)*SQRT(MinVALUE) + C(Indx)*LOG(MinVALUE) )

Count = 1

      ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )

DO WHILE (Count < MaxIter)

   CALL RndUnif( p, OtherSt, RN )        ! Generate RN from U(0,1)
   CALL RndUnif( p, OtherSt, RandNum )   ! Generate RandNum from h(y) = 1 / (MaxVALUE - MinVALUE)
   RandNum = RandNum*(MaxVALUE-MinVALUE)+MinVALUE;

   Gx = EXP(A(Indx) + B(Indx)*SQRT(RandNum) + C(Indx)*LOG(RandNum) )
   Gx = Gx / fMAX

   IF ( RN <= Gx ) THEN
      Count = MaxIter      ! Let's keep this deviate
   ELSE
      Count = Count + 1    ! try again
      RandNum = -1
   ENDIF

ENDDO


END SUBROUTINE RndTcohLLJ
!=======================================================================
SUBROUTINE RndTcoh_WF( p, OtherSt, RandNum, SpecModel )


IMPLICIT                            NONE

TYPE(RandNum_ParameterType),  INTENT(IN   ) :: p                   ! parameters for random number generation
TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt             ! other states for random number generation
REAL(ReKi),                   INTENT(  OUT) :: RandNum             ! numbers distributed with the pdf above
INTEGER(IntKi),               INTENT(IN   ) :: SpecModel           ! unique number identifying the spectral model being used

REAL(ReKi)                                  :: FnRange( 2)         ! the min and max values of the returned deviate
REAL(ReKi)                                  :: ParmsLN( 4)         ! parameters for the 4-parameter Log-Norm function
REAL(ReKi)                                  :: ParmsPIV(5)         ! parameters for the Pearson IV function


SELECT CASE ( SpecModel )
   CASE ( SpecModel_WF_UPW )
      ParmsLN = (/ 0.0, 0.132537201, 0.348791907, 1.781668096 /)  !parm(1) = 0 b/c there's no offset... it's a 3-parameter function
      FnRange = (/ 0.4, 200.0 /)

      CALL Rnd4ParmLogNorm( p, OtherSt, RandNum, ParmsLN, FnRange )

   CASE ( SpecModel_WF_07D )
      ParmsPIV = (/ 0.108721975, 5.705449915, 2.769408844, 0.475906651, 1.067616671 /)
      FnRange  = (/ 0.6, 30.0 /)

      CALL RndPearsonIV( p, OtherSt, RandNum, ParmsPIV, FnRange )

   CASE ( SpecModel_WF_14D )
      ParmsPIV = (/ 0.080526074, 13.51637204, 6.391924365, 1.197332751, 0.390220799 /)
      FnRange  = (/ 5.0, 40.0 /)

      CALL RndPearsonIV( p, OtherSt, RandNum, ParmsPIV, FnRange )

END SELECT

END SUBROUTINE RndTcoh_WF
!=======================================================================
SUBROUTINE RndUnif( p, OtherSt, RandUnifNum )

!This subroutine produces uniformly distributed random numbers, based on
!the pRNG that is requested in TurbSim's input file.  This routine assumes
!that the random number generator has been initialized earlier in the main
!program.

IMPLICIT                            NONE

TYPE(RandNum_ParameterType),  INTENT(IN)    :: p                   ! parameters for random number generation
TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt             ! other states for random number generation

REAL(ReKi), INTENT(OUT)                     :: RandUnifNum         ! Uniformly distributed numbers
REAL(ReKi)                                  :: RN(1)

IF ( p%pRNG == pRNG_INTRINSIC ) THEN !RNG_type == 'NORMAL'

   CALL RANDOM_NUMBER( RN )

ELSEIF (p%pRNG == pRNG_RANLUX) THEN !RNG_type == 'RANLUX'

   CALL RanLux ( RN )

ELSE

   CALL ARand( OtherSt%NextSeed(1), RN, 1,  1)

ENDIF


RandUnifNum = RN(1)


END SUBROUTINE RndUnif

!======================================================================
SUBROUTINE RndPhases(p, OtherSt, PhaseAngles, NPoints, NumFreq, US, ErrStat, ErrMsg)


IMPLICIT NONE

TYPE(RandNum_ParameterType),  INTENT(IN   ) :: p                                !< parameters for random number generation
TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt                          !< other states for random number generation
INTEGER(IntKi)              , INTENT(IN)    ::  US                              !< unit number of file in which to print a summary of the scaling used. If < 1, will not print summary.
INTEGER(IntKi)  ,             INTENT(OUT)   :: ErrStat                          !< error level/status
CHARACTER(*) ,                INTENT(OUT)   :: ErrMsg                           !< error message
                                                                                
INTEGER(IntKi)              , INTENT(IN   ) :: NPoints                          !< number of points being simulated
INTEGER(IntKi)              , INTENT(IN   ) :: NumFreq                          !< number of frequencies being simulated

REAL(ReKi)                  , INTENT(  OUT) :: PhaseAngles(NPoints,NumFreq,3)   !< random phases

! local variables

REAL(ReKi), ALLOCATABLE                     :: RandNum(:)                       ! contains the uniformly-distributed random numbers for all the points and frequencies

INTEGER(IntKi)                              :: Indx                             ! holds the next index in the RandNum array for the SNLWnd3d generator
INTEGER(IntKi)                              :: IVec                             ! loop counter (=number of wind components = 3)
INTEGER(IntKi)                              :: IFreq                            ! loop counter (=number of frequencies)
INTEGER(IntKi)                              :: J                                ! loop counter (=number of points on grid)
INTEGER(IntKi)                              :: NumPointsFreq                    ! number of points * number of frequency, or 1/3 the size of RandNum

INTEGER                                     :: LuxLevelOut, InitSeed
CHARACTER(20)                               :: NextSeedText


! get the uniformly distributed random numbers:


   CALL AllocAry( RandNum, SIZE(PhaseAngles), 'RandNum', ErrStat, ErrMsg )


      ! Reinitialize the random number generator ( it was initialized when the
      ! seeds were read in ) so that the same seed always generates the same
      ! random phases, regardless of previous randomizations in this code.

   CALL RandNum_Init(p, OtherSt, ErrStat, ErrMsg)
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

   
      ! Let's go ahead and get all the random numbers we will need for the entire
      ! run.  This (hopefully) will be faster than getting them one at a time,
      ! but it will use more memory.
      ! These pRNGs have been initialized in the GetInput() subroutine

   IF ( p%pRNG == pRNG_INTRINSIC ) THEN  ! RNG_type == 'NORMAL'

         !The first two real numbers in the RandSeed array are used as seeds
         !The number of seeds needed are compiler specific, thus we can't assume only 2 seeds anymore

      CALL RANDOM_NUMBER ( RandNum )

      ! Let's harvest the random seeds so that they can be used for the next run if desired.
      ! Write them to the summary file.

      CALL RANDOM_SEED ( GET = OtherSt%nextSeed(:) )  ! bjj: 16-jul-2014: without the (:), I get an "internal compiler error" here using Intel(R) Visual Fortran Compiler XE 12.1.3.300 [Intel(R) 64]

      NextSeedText = ' Harvested seed #'
      
   ELSEIF ( p%pRNG == pRNG_RANLUX ) THEN ! RNG_type == 'RANLUX'

      CALL RanLux ( RandNum )

      CALL RLuxAT ( LuxLevelOut, InitSeed, OtherSt%nextSeed(1), OtherSt%nextSeed(2) ) !luxury level, seed, nextSeed
      
      !InitSeed = p%RandSeed(1)????
      NextSeedText = ' K'

   ELSE
 
      NumPointsFreq = NPoints*NumFreq

      Indx = 1
      DO IVec = 1,3  ! 3 wind components
         CALL ARand( OtherSt%nextSeed(IVec), RandNum, NumPointsFreq,  Indx) 
         Indx = Indx + NumPointsFreq
      ENDDO

      NextSeedText = ' Next seed #'
      
   ENDIF
         
! set them to the range 0-2pi and   
! sort them so we get the same random numbers as previous versions of TurbSim   
                                       
   DO IVec = 1,3
      DO IFreq = 1,NumFreq
         DO J=1,NPoints
            Indx = IFreq + (J-1)*NumFreq + (IVec-1)*NPoints*NumFreq  ! This sorts the random numbers as they were done previously

            PhaseAngles(J,IFreq,IVec)  = TwoPi*RandNum(Indx)
         ENDDO ! J
      ENDDO !IFreq
   ENDDO !IVec         
                           
   call cleanup()
   
   IF ( US > 0 ) THEN
   
      WRITE(US,"(//,'Harvested Random Seeds after Generation of the Random Numbers:',/)")

      DO J = 1,SIZE( OtherSt%nextSeed  )
         WRITE(US,"(I13,A,I2)")  OtherSt%nextSeed(J), TRIM(NextSeedText), J
      END DO
            
   END IF
   
contains
   subroutine cleanup()
      IF (ALLOCATED(RandNum)) DEALLOCATE(RandNum)
   end subroutine cleanup 
         
END SUBROUTINE RndPhases
!=======================================================================


END MODULE TS_RandNum
