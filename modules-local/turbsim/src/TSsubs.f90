MODULE TSSubs

   USE                     NWTC_Library
   IMPLICIT                NONE

      ! Create interface for a generic getWindSpeed that actually uses specific routines.

   INTERFACE getWindSpeed
      MODULE PROCEDURE getWindSpeedVal
      MODULE PROCEDURE getWindSpeedAry
   END INTERFACE


CONTAINS

!=======================================================================
SUBROUTINE ARand(ix, RandNum_Ary,I, RNG_start)

IMPLICIT                     NONE

      ! Passed variables

REAL(ReKi)               :: RandNum_Ary(:)    ! Output: random numbers
INTEGER, INTENT(IN)      :: I                 ! Input: Size of RandNum_Ary to be filled
INTEGER, INTENT(INOUT)   :: ix                ! Input/Output: Seed  !BONNIE: should this be set to 32-bits, not default integer size?
INTEGER, INTENT(IN)      :: RNG_start

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
FUNCTION Beta1( X, A0, A1, A2, A3 )

   ! This function is used in the calculation of the Wind Farm models' PSD

IMPLICIT                NONE

REAL(ReKi), INTENT(IN) :: X           ! Function input
REAL(ReKi), INTENT(IN) :: A0          ! Function input
REAL(ReKi), INTENT(IN) :: A1          ! Function input
REAL(ReKi), INTENT(IN) :: A2          ! Function input
REAL(ReKi), INTENT(IN) :: A3          ! Function input
REAL(ReKi)             :: Beta1       ! Function result

REAL(ReKi)             :: tmp1        ! temporary variable
REAL(ReKi)             :: tmp2        ! temporary variable


      tmp1  = X - A1
      tmp2  = A2 / 2.0

      Beta1 =          A0 / (1.0 + EXP( ( tmp1 + tmp2 ) / (-A3) )) * &
             (1.0 - ( 1.0 / (1.0 + EXP( ( tmp1 - tmp2 ) / (-A3) )) ))


RETURN
END FUNCTION Beta1
!=======================================================================
FUNCTION BETA2(X,A0,A1,A2)

   ! This function is used in the calculation of the Wind Farm models' PSD

IMPLICIT                NONE

REAL(ReKi),INTENT(IN) :: X           ! Function input
REAL(ReKi),INTENT(IN) :: A0          ! Function input
REAL(ReKi),INTENT(IN) :: A1          ! Function input
REAL(ReKi),INTENT(IN) :: A2          ! Function input
REAL(ReKi)            :: Beta2       ! Function output


      Beta2 = ( A0 / ( 2.50663 * A2 ) ) * EXP( -0.5 * ( (X-A1) / A2 )**2 )

RETURN
END FUNCTION Beta2
!=======================================================================
FUNCTION BETA3(X,A0,A1,A2)

   ! This function is used in the calculation of the Wind Farm models' PSD

IMPLICIT                NONE

REAL(ReKi)            :: Beta3       ! Function output
REAL(ReKi),INTENT(IN) :: X           ! Function input
REAL(ReKi),INTENT(IN) :: A0          ! Function input
REAL(ReKi),INTENT(IN) :: A1          ! Function input
REAL(ReKi),INTENT(IN) :: A2          ! Function input

      Beta3 = 0.5 * A0 * ( 1.0 + Beta10( (X-A1) / (1.414*A2) ) )

RETURN
END FUNCTION Beta3
!=======================================================================
FUNCTION Beta4(X,A0,A1,A2,A3)

   ! This function is used in the calculation of the Wind Farm models' PSD

IMPLICIT                NONE

REAL(ReKi)            :: Beta4       ! Function output
REAL(ReKi),INTENT(IN) :: X           ! Function input
REAL(ReKi),INTENT(IN) :: A0          ! Function input
REAL(ReKi),INTENT(IN) :: A1          ! Function input
REAL(ReKi),INTENT(IN) :: A2          ! Function input
REAL(ReKi),INTENT(IN) :: A3          ! Function input


      Beta4 = A0 * EXP( -X/A1 ) + A2 * EXP( -X/A3 )

RETURN
END FUNCTION Beta4
!=======================================================================
FUNCTION Beta5(X,A0,A1,A2)

   ! This function is used in the calculation of the Wind Farm models' PSD

IMPLICIT                NONE

REAL(ReKi)            :: Beta5       ! Function output
REAL(ReKi),INTENT(IN) :: X           ! Function input
REAL(ReKi),INTENT(IN) :: A0          ! Function input
REAL(ReKi),INTENT(IN) :: A1          ! Function input
REAL(ReKi),INTENT(IN) :: A2          ! Function input

      Beta5 = A0 / ( 1.0 + EXP( -(X-A1) / A2 ) )

RETURN
END FUNCTION Beta5
!=======================================================================
FUNCTION Beta6(A,X)

   ! This function is used in the calculation of the Wind Farm models' PSD

IMPLICIT                NONE

REAL(ReKi)            :: Beta6       ! Function output
REAL(ReKi),INTENT(IN) :: A           ! Function input
REAL(ReKi),INTENT(IN) :: X           ! Function input

   IF ( ( X < 0.0 ) .OR. ( A <= 0.0 ) ) THEN
     CALL TS_Abort( ' An error occurred in the Beta6 function. ' )
   ENDIF

   IF ( X < A + 1.0 ) THEN
      CALL Beta8( Beta6, A, X ) ! was CALL Beta8( GAMSER, A, X )
      ! Beta6 = GAMSER
   ELSE
      CALL Beta7( Beta6, A, X)  ! was CALL Beta7( GAMMCF, A, X )
      Beta6 = 1.0 - Beta6 ! was Beta6 = 1.0 - GAMMCF
   ENDIF

RETURN
END FUNCTION Beta6
!=======================================================================
SUBROUTINE Beta7(GAMMCF, A, X )

   ! This subroutine is used in the calculation of the Wind Farm models' PSD


IMPLICIT                   NONE

REAL(ReKi),INTENT(OUT)    :: GAMMCF        ! Subroutine Output
REAL(ReKi),INTENT(IN)     :: A             ! Subroutine Input
REAL(ReKi),INTENT(IN)     :: X             ! Subroutine Input

REAL(ReKi)                :: GLN
REAL(ReKi)                :: g
REAL(ReKi)                :: gOld
REAL(ReKi)                :: A0
REAL(ReKi)                :: A1
REAL(ReKi)                :: B0
REAL(ReKi)                :: B1
REAL(ReKi)                :: FAC
REAL(ReKi)                :: AN
REAL(ReKi)                :: ANA
REAL(ReKi)                :: ANF

REAL(ReKi), PARAMETER     :: eps  = 3.0E-7
REAL(ReKi), PARAMETER     :: ITmax = 100.0

LOGICAL                   :: continueIT


IF ( X <= 0.0 ) THEN
   CALL TS_Abort ( 'Input variable X must be positive in function Beta7()' )
ENDIF

   gOld = 0.0

   A0   = 1.0
   A1   = X
   B0   = 0.0
   B1   = 1.0
   FAC  = 1.0

   AN = 0.0
   continueIT = .TRUE.

   DO WHILE ( ( AN < ITmax ) .AND. continueIT )

      AN  = AN + 1.0

      ANA = AN - A
      A0  = ( A1 + A0*ANA )*FAC
      B0  = ( B1 + B0*ANA )*FAC

      ANF = AN*FAC
      A1  = X*A0 + ANF*A1
      B1  = X*B0 + ANF*B1

      IF ( A1 /= 0.0 ) THEN
         FAC = 1.0 / A1
         g   = B1*FAC

         IF( ABS( ( g - gOld ) / g ) < eps)  continueIT = .FALSE.

         gOld = g
      ENDIF

   ENDDO

   IF ( continueIT ) THEN
      CALL TS_Abort (' Value of A is too large or ITMAX is too small in BETA7. ' )
   ENDIF

   GLN  = Beta9( A )

   GAMMCF = EXP( -X + A*LOG( X ) - GLN ) * G

RETURN
END SUBROUTINE Beta7
!=======================================================================
SUBROUTINE Beta8(GAMSER,A,X)

   ! This subroutine is used in the calculation of the Wind Farm models' PSD

IMPLICIT                  NONE

REAL(ReKi),INTENT(OUT)  :: GAMSER        ! Subroutine Output
REAL(ReKi),INTENT(IN)   :: A             ! Subroutine Input
REAL(ReKi),INTENT(IN)   :: X             ! Subroutine Input

REAL(ReKi)              :: GLN
REAL(ReKi)              :: AP
REAL(ReKi)              :: Sum
REAL(ReKi)              :: del

REAL(ReKi),PARAMETER    :: eps = 3.0E-7  ! Tolerance
INTEGER,PARAMETER       :: ITmax = 100   ! Maximum loop iterations


INTEGER                 :: N             ! Loop counter

LOGICAL                 :: continueIT


   IF ( ( X > 0.0 ) .AND. ( A /= 0.0 ) ) THEN
      continueIT = .TRUE.

      AP  = A
      Sum = 1.0 / A
      del = Sum

      N = 1

      DO WHILE ( ( N <= ITmax ) .AND. ( continueIT ) )
         AP  = AP + 1.0
         del = del * X / AP
         Sum = Sum + del

         IF( ABS(del) < ABS(Sum) * eps ) continueIT = .FALSE.

         N = N + 1
      ENDDO

      IF ( continueIT ) THEN
         CALL TS_Abort (' Value of A is too large or ITMAX is too small in BETA8. ' )
      ENDIF

      GLN = Beta9( A )
      GAMSER = Sum * EXP( -X + A * LOG(X) - GLN)

   ELSEIF ( X == 0.0 ) THEN

         GAMSER = 0.0

   ELSE ! ( X < 0.0 )
        CALL TS_Abort( 'Error in Subroutine Beta8.  Invalid input.' )

   ENDIF

RETURN
END SUBROUTINE Beta8
!=======================================================================
FUNCTION Beta9(XX)

   ! This function is used in the calculation of the Wind Farm models' PSD

IMPLICIT                NONE

REAL(ReKi)             :: Beta9   ! Output value
REAL(ReKi),INTENT(IN)  :: XX      ! Input value

REAL(ReKi)             :: X
REAL(ReKi)             :: Tmp
REAL(ReKi)             :: SER

REAL(ReKi)             :: Cof(6) = (/ 76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.120858003E-2, -0.536382E-5 /)
REAL(ReKi), PARAMETER  :: STP  = 2.50662827465

INTEGER                :: J       ! Loop counter

IF ( XX <= -4.5 ) THEN
   CALL TS_Abort ( 'Input variable XX must be larger than -4.5 in function Beta9()' )
ENDIF


   X = XX - 1.0
   Tmp =   X + 5.5
   Tmp = ( X + 0.5 ) * LOG( Tmp ) - Tmp

   SER = 1.0

   DO J = 1,6
      X   = X + 1.0
      SER = SER + Cof(J) / X
   ENDDO


IF ( SER <= 0.0 ) THEN
   CALL TS_Abort ( 'Variable SER must be larger than 0.0 in function Beta9()' )
ENDIF


   Beta9 = Tmp + LOG( STP*SER )

RETURN
END FUNCTION Beta9
!=======================================================================
FUNCTION Beta10(X)

   ! This function is used in the calculation of the Wind Farm models' PSD

IMPLICIT                  NONE

REAL(ReKi)             :: Beta10
REAL(ReKi), INTENT(IN) :: X
REAL(ReKi), PARAMETER  :: Tmp = 0.5


   IF ( X < 0.0 ) THEN
      Beta10 = -Beta6(Tmp, X**2)
   ELSE
      Beta10 =  Beta6(Tmp, X**2)
   ENDIF

RETURN
END FUNCTION Beta10
!=======================================================================
SUBROUTINE CalcEvents( WindSpeed, MaxCTKE, Height )

      ! This subroutine calculates what events to use and when to use them.
      ! It computes the number of timesteps in the file, NumCTt.

   USE                         TSMods

   IMPLICIT                    NONE

      ! passed variables
REAL(ReKi), INTENT(IN)      :: WindSpeed           ! Hub height wind speed
REAL(ReKi), INTENT(OUT)     :: MaxCTKE             ! Maximum CTKE of events we've picked
REAL(ReKi), INTENT(IN)      :: Height              ! Height for expected length PDF equation

      ! local variables
REAL(ReKi)                  :: EventTimeSum = 0.0  ! Amount of time the coherent structure takes
REAL(ReKi)                  :: ExpectedTime        ! Amount of time the coherent structures should take
REAL(ReKi)                  :: iA                  ! Variable used to calculate IAT
REAL(ReKi)                  :: iB                  ! Variable used to calculate IAT
REAL(ReKi)                  :: iC                  ! Variable used to calculate IAT
REAL(ReKi)                  :: lambda              ! The expected value of interarrival times for the Poisson process
REAL(ReKi)                  :: rn                  ! random number
REAL(ReKi)                  :: TEnd                ! End time for the current event
REAL(ReKi)                  :: TStartNext   = 0.0  ! temporary start time for next event

INTEGER                     :: IStat               ! Status of memory allocation
INTEGER                     :: NewEvent            ! event number of the new event
INTEGER                     :: NumCompared         ! Number of events we've compared

LOGICAL(1)                  :: Inserted            ! Whether an event was inserted here

TYPE(Event), POINTER        :: PtrCurr  => NULL()  ! Pointer to the current event in the list
TYPE(Event), POINTER        :: PtrNew   => NULL()  ! A new event to be inserted into the list


      ! Compute the mean interarrival time and the expected length of events

   SELECT CASE ( TRIM(TurbModel) )

      CASE ( 'NWTCUP', 'NONE', 'USRVKM' )
         lambda = -0.000904*Rich_No + 0.000562*WindSpeed + 0.001389
         lambda = 1.0 / lambda

         IF ( TurbModel(1:4) == 'NONE' ) THEN
            ExpectedTime = 600.0
         ELSE
            CALL RndModLogNorm( ExpectedTime, Height )
         ENDIF

      CASE ( 'GP_LLJ', 'SMOOTH' , 'TIDAL', 'RIVER') ! HYDRO: added 'TIDAL' and 'RIVER' to the spectral models that get handled this way.
         iA     =        0.001797800 + (7.17399E-10)*Height**3.021144723
         iB     =  EXP(-10.590340100 - (4.92440E-05)*Height**2.5)
         iC     = SQRT(  3.655013599 + (8.91203E-06)*Height**3  )
         lambda = iA + iB*MIN( (WindSpeed**iC), HUGE(iC) )  ! lambda = iA + iB*(WindSpeed**iC)
         lambda = 1.0 / lambda

         CALL RndTcohLLJ( ExpectedTime, Height )

      CASE ( 'WF_UPW' )
        lambda = 0.000529*WindSpeed + 0.000365*Rich_No - 0.000596
        lambda = 1.0 / lambda

         CALL RndTcoh_WF( ExpectedTime )

      CASE ( 'WF_07D' )
         lambda = 0.000813*WindSpeed - 0.002642*Rich_No + 0.002676
         lambda = 1.0 / lambda

         CALL RndTcoh_WF( ExpectedTime )

      CASE ( 'WF_14D' )
         lambda = 0.001003*WindSpeed - 0.00254*Rich_No - 0.000984
         lambda = 1.0 / lambda

         CALL RndTcoh_WF( ExpectedTime )

      CASE DEFAULT
         !This should not happen

   END SELECT

   ExpectedTime = ExpectedTime * ( UsableTime - CTStartTime ) / 600.0  ! Scale for use with the amount of time we've been given


!BONNIE: PERHAPS WE SHOULD JUST PUT IN A CHECK THAT TURNS OFF THE COHERENT TIME STEP FILE IF THE
!        CTSTARTTIME IS LESS THAN THE USABLETIME... MAYBE WHEN WE'RE READING THE INPUT FILE...
ExpectedTime = MAX( ExpectedTime, REAL(0.0,ReKi) )  ! This occurs if CTStartTime = 0

      ! We start by adding events at random times

   NumCTEvents = 0                                    ! Number of events = length of our linked list
   NumCTt      = 0                                    ! Total number of time steps in the events we've picked
   MaxCTKE     = 0.0                                  ! Find the maximum CTKE for the events that we've selected


   CALL RndExp(rn, lambda)                            ! Assume the last event ended at time zero

   TStartNext = rn / 2.0

   IF ( KHtest ) THEN
      ExpectedTime = UsableTime   / 2                 ! When testing, add coherent events for half of the record
      TStartNext   = ExpectedTime / 2                 ! When testing, start about a quarter of the way into the record
   ENDIF

   IF ( TStartNext < CTStartTime ) THEN
      TStartNext = TStartNext + CTStartTime           ! Make sure the events start after time specified by CTStartTime
   ENDIF

   IF ( TStartNext > 0 ) NumCTt = NumCTt + 1          ! Add a point before the first event

   DO WHILE ( TStartNext < UsableTime .AND. EventTimeSum < ExpectedTime )

      CALL RndUnif( rn )

      NewEvent = INT( rn*( NumEvents - 1.0 ) ) + 1
      NewEvent = MAX( 1, MIN( NewEvent, NumEvents ) ) ! take care of possible rounding issues....


      IF ( .NOT. ASSOCIATED ( PtrHead ) ) THEN

         ALLOCATE ( PtrHead, STAT=IStat )             ! The pointer %Next is nullified in allocation

         IF ( IStat /= 0 ) THEN
            CALL TS_Abort ( 'Error allocating memory for new event.' )
         ENDIF

         PtrTail => PtrHead

      ELSE

         ALLOCATE ( PtrTail%Next, STAT=IStat )     ! The pointer PtrTail%Next%Next is nullified in allocation

         IF ( IStat /= 0 ) THEN
            CALL TS_Abort ( 'Error allocating memory for new event.' )
         ENDIF

         PtrTail => PtrTail%Next                   ! Move the pointer to point to the last record in the list

      ENDIF

      PtrTail%EventNum     = NewEvent
      PtrTail%TStart       = TStartNext
      PtrTail%delt         = EventLen( NewEvent ) / EventTS( NewEvent )          ! the average delta time in the event
      PtrTail%Connect2Prev = .FALSE.

      MaxCTKE              = MAX( MaxCTKE, pkCTKE( NewEvent ) )
      NumCTEvents          = NumCTEvents + 1

      TEnd = TStartNext + EventLen( NewEvent )


      IF ( KHtest ) THEN
         TStartNext   = UsableTime + TStartNext !TEnd + PtrTail%delt ! Add the events right after each other
      ELSE

         DO WHILE ( TStartNext <= TEnd )

            CALL RndExp(rn, lambda)                                                 ! compute the interarrival time
            TStartNext        = TStartNext + rn !+ EventLen( NewEvent )

         ENDDO

      ENDIF


      IF ( (TStartNext - TEnd) > PtrTail%delt ) THEN
         NumCTt = NumCTt + EventTS( NewEvent ) + 2                                  ! add a zero-line (essentially a break between events)
      ELSE
         NumCTt = NumCTt + EventTS( NewEvent ) + 1
      ENDIF

      EventTimeSum     = EventTimeSum + EventLen( NewEvent )

   ENDDO

      ! Write the number of separate events to the summary file

   IF (KHtest) THEN
      WRITE ( US,'(/)' )
   ELSE
      WRITE ( US,'(//A,F8.3," seconds")' ) 'Average expected time between events = ',lambda
   ENDIF

WRITE ( US, '(A,I8)'   )             'Number of coherent events            = ',NumCTEvents
WRITE ( US, '(A,F8.3," seconds")' )  'Predicted length of coherent events  = ',ExpectedTime

      ! Next, we start concatenating events until there is no space or we exceed the expected time

   IF ( TurbModel(1:4) /= 'NONE' ) THEN

      NumCompared = 0

      DO WHILE ( EventTimeSum < ExpectedTime .AND. NumCompared < NumCTEvents )

         CALL RndUnif( rn )

         NewEvent = INT( rn*( NumEvents - 1.0 ) ) + 1
         NewEvent = MAX( 1, MIN( NewEvent, NumEvents ) )    ! take care of possible rounding issues....

         NumCompared = 0
         Inserted    = .FALSE.

         DO WHILE ( NumCompared < NumCTEvents .AND. .NOT. Inserted )

            IF ( .NOT. ASSOCIATED ( PtrCurr ) ) THEN        ! Wrap around to the beginning of the list
               PtrCurr => PtrHead
            ENDIF


               ! See if the NewEvent fits between the end of event pointed to by PtrCurr and the
               ! beginning of the event pointed to by PtrCurr%Next

            IF ( ASSOCIATED( PtrCurr%Next ) ) THEN
               TStartNext = PtrCurr%Next%TStart
            ELSE !We're starting after the last event in the record
               TStartNext = UsableTime + 0.5 * EventLen( NewEvent )  ! We can go a little beyond the end...
            ENDIF

            IF ( TStartNext - (PtrCurr%TStart + EventLen( PtrCurr%EventNum ) + PtrCurr%delt) > EventLen( NewEvent ) ) THEN

               Inserted = .TRUE.

               ALLOCATE ( PtrNew, STAT=IStat )           ! The pointer %Next is nullified in allocation

               IF ( IStat /= 0 ) THEN
                  CALL TS_Abort ( 'Error allocating memory for new event.' )
               ENDIF

               PtrNew%EventNum      = NewEvent
               PtrNew%TStart        = PtrCurr%TStart + EventLen( PtrCurr%EventNum )
               PtrNew%delt          = EventLen( NewEvent ) / EventTS( NewEvent )          ! the average delta time in the event
               PtrNew%Connect2Prev  = .TRUE.

               PtrNew%Next  => PtrCurr%Next
               PtrCurr%Next => PtrNew
               PtrCurr      => PtrCurr%Next    ! Let's try to add the next event after the other events

               MaxCTKE              = MAX( MaxCTKE, pkCTKE( NewEvent ) )
               NumCTEvents          = NumCTEvents + 1
               NumCTt               = NumCTt + EventTS( NewEvent )  ! there is no break between events
                                   !(we may have one too many NumCTt here, so we'll deal with it when we write the file later)
               EventTimeSum         = EventTimeSum + EventLen( NewEvent )


            ELSE

               NumCompared = NumCompared + 1

            ENDIF

            PtrCurr => PtrCurr%Next

         ENDDO ! WHILE (NumCompared < NumCTEvents .AND. .NOT. Inserted)

      ENDDO ! WHILE (EventTimeSum < ExpectedTime .AND. NumCompared < NumCTEvents)

   ENDIF ! TurbModel /= 'NONE'

IF ( NumCTt > 0 ) THEN
   EventTimeStep = EventTimeSum / NumCTt                                          ! Average timestep of coherent event data
ELSE
   EventTimeStep = 0.0
ENDIF


WRITE ( US, '(A,F8.3," seconds")' )   'Length of coherent events            = ',EventTimeSum

END SUBROUTINE CalcEvents
!=======================================================================
SUBROUTINE ChebyshevVals(coeffs,x,y,MinX,MaxX)

   IMPLICIT NONE

      ! Passed variables

   REAL(ReKi), DIMENSION(:),   INTENT(IN)  :: coeffs  ! Coefficients defined on [-1,1]
   REAL(ReKi), DIMENSION(:),   INTENT(IN)  :: x       ! The x values where f(x)=y is desired
   REAL(ReKi), DIMENSION(:),   INTENT(OUT) :: y       ! The desired function values
   REAL(SiKi),                 INTENT(IN)  :: MinX    ! Min X of the interval the coeffs were originally calculated for
   REAL(SiKi),                 INTENT(IN)  :: MaxX    ! Max X of the interval the coeffs were originally calculated for

   INTEGER                                 :: i,j
   INTEGER                                 :: SC
   INTEGER                                 :: SX
   INTEGER                                 :: SY

   REAL(DbKi), DIMENSION(size(x))          :: x_hat
   REAL(DbKi), DIMENSION(size(coeffs))     :: BasisFn  ! The Chebyshev basis functions evaluated at x_hat

   SC = size(coeffs)
   SX = size(x)
   SY = size(y)

   IF (SX /= SY) THEN
      CALL TS_Abort( 'The x and y vectors in ChebyshevVals() must be the same size.' )
      SX = MIN(SX,SY)
      SY = SX
   ENDIF

   x_hat = (2.0*REAL(x(:),DbKi) - MaxX - MinX)/(MaxX - MinX)  ! Transform from [MinX,MaxX] to [-1,1]


   DO i=1,SX
      CALL ChebyshevFuncs( x_hat(i), BasisFn )

      y(i) = 0.
      DO j=1,SC
         y(i) = y(i) + coeffs(j)*REAL(BasisFn(j),ReKi)
      ENDDO

   ENDDO

   RETURN
   CONTAINS
     !-----------------------------------------------------------------------
      SUBROUTINE ChebyshevFuncs( x, Px )

         REAL(DbKi), INTENT(IN)                 ::  x
         REAL(DbKi), INTENT(OUT), DIMENSION(:)  ::  Px  ! The basis functions evaluated at x

         INTEGER                                ::  I
         INTEGER                                ::  S_Px  ! Size of Px, determines how many basis functions to use (i.e. order of the polynomial - 1)


            S_Px = SIZE(Px)

         !----------------------------
         ! Define the basis functions:
         !----------------------------
             Px(1) = 1

             IF (S_Px > 1) THEN

               Px(2) = x

                  ! Define Chebyshev polynomials recursively

               DO I=3,S_Px
                  Px(I) = 2.*x*Px(I-1) - Px(I-2)
               ENDDO

            ENDIF  !S_Px > 1

      END SUBROUTINE ChebyshevFuncs
END SUBROUTINE ChebyshevVals
!=======================================================================
SUBROUTINE CohSpecVMat( LC, NSize, Comp )

   ! This subroutine computes the coherence between two points on the grid.
   ! It stores the symmetric coherence matrix, packed into variable "Matrix"
   ! This replaces what formerly was the "ExCoDW" matrix.

USE                           TSMods

IMPLICIT                      NONE

   ! Passed variables

REAL(ReKi),   INTENT(IN)   :: LC             ! IEC coherency scale parameter
INTEGER,      INTENT(IN)   :: NSize          ! Size of dimension 2 of Matrix
CHARACTER(*), INTENT(IN)   :: Comp(3)        ! u, v, or w

   ! Internal variables

REAL(ReKi)                 :: CohDec         ! The coherence decrement
REAL(ReKi)                 :: CDistUM        ! Temporary variable for calculating coherence * dist
REAL(ReKi)                 :: CPh            ! Cosine of the random phase
REAL(ReKi), ALLOCATABLE    :: Dist(:)        ! The distance between points
REAL(ReKi)                 :: DistLC
REAL(ReKi), ALLOCATABLE    :: DistU(:)
REAL(ReKi), ALLOCATABLE    :: DistZMExp(:)
REAL(ReKi)                 :: dY             ! the lateral distance between two points
REAL(ReKi)                 :: Ph             ! Phase angle.
REAL(ReKi)                 :: SPh            ! Sine of the random phase
REAL(ReKi)                 :: UM             ! The mean wind speed of the two points
REAL(ReKi)                 :: ZM             ! The mean height of the two points

INTEGER                    :: J
INTEGER                    :: JJ             ! Index of point J
INTEGER                    :: JJ1
INTEGER                    :: JY             ! Index of y-value of point J
INTEGER                    :: JZ             ! Index of z-value of point J
INTEGER                    :: I
INTEGER                    :: IF1            ! Index to real part of vector
INTEGER                    :: IF2            ! Index to complex part of vector
INTEGER                    :: IFreq
INTEGER                    :: II             ! The index of point I
INTEGER                    :: Indx
INTEGER                    :: IRand          ! Index of Random Number of Array
INTEGER                    :: IVec           ! wind component, 1=u, 2=v, 3=w
INTEGER                    :: IY             ! The index of the y-value of point I
INTEGER                    :: IZ             ! Index of the z-value of point I
INTEGER                    :: Stat

LOGICAL                    :: IdentityCoh

! ------------ arrays allocated -------------------------
Stat = 0.

IF ( .NOT. ALLOCATED( Dist ) ) ALLOCATE( Dist(NSize),      STAT=Stat )

IF ( Stat /= 0 ) THEN
   CALL TS_Abort('Error allocating '//TRIM( Int2LStr( ReKi*NSize/1024**2 ) )//' MB for the Dist coherence array.')
ENDIF

IF ( .NOT. ALLOCATED( DistU     ) ) ALLOCATE( DistU(NSize),     STAT=Stat )

IF ( Stat /= 0 )  THEN
   CALL TS_Abort('Error allocating '//TRIM( Int2LStr( ReKi*NSize/1024**2 ) )//' MB for the turbulence DistU coherence array.')
ENDIF

IF ( .NOT. ALLOCATED( DistZMExp ) ) ALLOCATE( DistZMExp(NSize), STAT=Stat )

IF ( Stat /= 0 )  THEN
   CALL TS_Abort('Error allocating '//TRIM( Int2LStr( ReKi*NSize/1024**2 ) )//' MB for the turbulence DistZMExp coherence array.')
ENDIF

!! Initialize the velocity matrix  !V(:,:,:) = 0.0 - done outside the subroutine

!--------------------------------------------------------------------------------
! Calculate the distances and other parameters that don't change with frequency
!---------------------------------------------------------------------------------

IF ( (TurbModel(1:3) == 'IEC') .OR. (TurbModel == 'MODVKM') ) THEN
   InCohB(:)    = 0.12/LC
   DistZMExp(:) = 1.0
ENDIF
   
         II = 0
POINT_I: DO IZ=1,ZLim   !NumGrid_Z
            DO IY=1,IYmax(IZ) !NumGrid_Y
      
               II = II + 1                            ! Index of point I: S(I)  !equivalent to II = ( IZ - 1 )*NumGrid_Y + IY
               IF (II > NTot) EXIT POINT_I            ! Don't go past the end of the array; this exits the IY loop

               JJ = 0
POINT_J:       DO JZ=1,IZ
                  DO JY=1,IYmax(JZ) !NumGrid_Y

                     JJ = JJ + 1                      ! Index of point J: S(J)  !equivalent to JJ = ( JZ - 1 )*NumGrid_Y + JY

                     IF ( JJ > II )  EXIT POINT_J     ! The coherence matrix is symmetric

                     IF ( IZ > NumGrid_Z ) THEN       ! Get the correct location if we're using an extra point for the hub
                        I = YLim
                        IF ( JZ > NumGrid_Z ) THEN
                           J = YLim
                        ELSE
                           J = JY
                        ENDIF
                     ELSE
                        I = IY
                        J = JY
                     ENDIF

                     JJ1       = JJ - 1
                     Indx      = NTot*JJ1 - JJ*JJ1/2 + II   !Index of matrix ExCoDW (now Matrix), coherence between points I & J

                     IF ( .NOT. PeriodicY ) THEN
                        Dist(Indx)= SQRT( ( Y(I) - Y(J) )**2  + ( Z(IZ) - Z(JZ) )**2 )
                     ELSE
                        dY = Y(I) - Y(J)
                        IF (dY > 0.5*GridWidth ) THEN
                           dY = dY - GridWidth - GridRes_Y
                        ELSE IF (dY < -0.5*GridWidth ) THEN 
                           dY = dY + GridWidth + GridRes_Y
                        END IF

                        Dist(Indx)= SQRT( ( dY )**2  + ( Z(IZ) - Z(JZ) )**2 )
                     END IF

                     IF ( (TurbModel(1:3) == 'IEC')  .OR. (TurbModel == 'MODVKM')) THEN
                        DistU(Indx) = Dist(Indx)/UHub
!                           TRH(Indx) = EXP( -InCDec(IVec)*SQRT( ( Freq(IFreq) * Dist / UHub )**2 + (0.12*Dist/LC)**2 ) )
                     ELSE
                        UM       = 0.5*( U(IZ) + U(JZ) )
                        ZM       = 0.5*( Z(IZ) + Z(JZ) )
                        
                        DistU(Indx)     = Dist(Indx)/UM                           
                        DistZMExp(Indx) = ( Dist(Indx)/ZM )**COHEXP     ! Note: 0**0 = 1

!                       TRH(Indx) = EXP( -InCDec(IVec) * DistZMExp*SQRT( ( Freq(IFreq)* DistU )**2 + (InCohB(IVec)*Dist)**2 ) )
                     ENDIF !TurbModel

                  ENDDO    ! JY
            ENDDO POINT_J  ! JZ

      ENDDO             ! IY
   ENDDO POINT_I        ! IZ

IF ( COH_OUT ) THEN

      ! Write the coherence for three frequencies, for debugging purposes

      CALL OpenFOutFile( UC, TRIM(RootName)//'.coh' )
!      WRITE( UC, '(A4,X,3(A16),2(A11))' ) 'Comp','Distance', 'Avg Height', 'Avg Wind Speed', 'c(F1)','c(F2)'
      
      WRITE( UC, '(A4,X,A16,1X,'//INT2LSTR(NSize)//'(G10.4,1X))' ) 'Comp','Freq',(I,I=1,NSize)
      WRITE( UC,   '(5X,A16,1X,'//INT2LSTR(NSize)//'(G10.4,1X))' ) 'Distance',     Dist(:)
      WRITE( UC,   '(5X,A16,1X,'//INT2LSTR(NSize)//'(G10.4,1X))' ) 'Distance/U',   DistU(:)
      WRITE( UC,   '(5X,A16,1X,'//INT2LSTR(NSize)//'(G10.4,1X))' ) '(r/Z)^CohExp', DistZMExp(:)
            
ENDIF

DO IVec = 1,3

   IF ( ( (TurbModel(1:3) == 'IEC') .OR. (TurbModel == 'MODVKM')) .AND. ( IVec /= 1 ) ) THEN
         ! There is no coherence defined for the v or w component of the IEC spectral models
      IdentityCoh = .TRUE.
   ELSE
      IdentityCoh = .FALSE.
   ENDIF

   CALL WrScr ( '    '//Comp(IVec)//'-component matrices' )
   
   !--------------------------------------------------------------------------------
   ! Calculate the coherence, Veers' H matrix (CSDs), and the fourier coefficients
   !---------------------------------------------------------------------------------

   IF2 = 0
   DO IFREQ = 1,NumFreq

      !---------------------------------------------------
      ! Calculate the coherence and Veers' H matrix (CSDs)
      !---------------------------------------------------
      IF (IdentityCoh) THEN

            ! -----------------------------------------------------------------------------------      
            !  The coherence is the Identity (as is Cholesky); the Veers' H matrix is as follows:   
            ! -----------------------------------------------------------------------------------      

         Indx = 1
         DO J = 1,NTot ! The column number

               ! The diagonal entries of the matrix:

            TRH(Indx) = SQRT( ABS( S(IFreq,J,IVec) ) )

               ! The off-diagonal values:
            Indx = Indx + 1
            DO I = J+1,NTot ! The row number
               TRH(Indx) = 0.0
               Indx = Indx + 1
            ENDDO ! I
         ENDDO ! J
            
      ELSE 
      
            ! -----------------------------------------------      
            ! Create the coherence matrix for this frequency      
            ! -----------------------------------------------      
            
         DO II=1,NTot
            DO JJ=1,II

                  JJ1       = JJ - 1
                  Indx      = NTot*JJ1 - JJ*JJ1/2 + II   !Index of matrix ExCoDW (now Matrix), coherence between points I & J

                  TRH(Indx) = EXP( -1.0 * InCDec(IVec) * DistZMExp(Indx)* & 
                                      SQRT( (Freq(IFreq)*DistU(Indx) )**2 + (InCohB(IVec)*Dist(Indx))**2 ) )

            ENDDO ! JJ
         ENDDO ! II                  
         
         IF (COH_OUT) THEN
!            IF (IFreq == 1 .OR. IFreq == NumFreq) THEN
               WRITE( UC, '(I3,2X,F15.5,1X,'//INT2LSTR(NSize)//'(G10.4,1X))' ) IVec, Freq(IFreq), TRH(1:NSize)
!            ENDIF
         ENDIF

            ! -------------------------------------------------------------      
            ! Calculate the Cholesky factorization for the coherence matrix
            ! -------------------------------------------------------------      
         
         CALL SPPTRF( 'L', NTot, TRH, Stat )  ! 'L'ower triangular 'TRH' matrix (packed form), of order 'NTot'; returns Stat

         IF ( Stat /= 0 ) THEN
            CALL TS_Abort('Error '//TRIM(Int2LStr(Stat))//' in the Cholesky factorization occurred at frequency '//&
                           TRIM(Int2LStr(IFreq))//' ('//TRIM(Flt2LStr(Freq(IFreq)))//' Hz)'//&
                        '. The '//Comp(IVec)//'-component coherence matrix cannot be factored.  '//&
                        'Check the input file for invalid physical properties or modify the coherence exponent '//&
                        'or grid spacing.')
         ENDIF
         
            ! -------------------------------------------------------------      
            ! Create the lower triangular matrix, H, from Veer's method
            ! -------------------------------------------------------------      
         
         Indx = 1
         DO J = 1,NTot  ! Column
            DO I = J,NTot ! Row

                  ! S(IFreq,I,IVec) should never be less than zero, but the ABS makes sure...

               TRH(Indx) = TRH(Indx) * SQRT( ABS( S(IFreq,I,IVec) ) )

               Indx = Indx + 1

            ENDDO !I
         ENDDO !J
         
      ENDIF !IdentityCoh
      
      ! -------------------------------------------------------------      
      ! Calculate the correlated fourier coefficients.   
      ! -------------------------------------------------------------

      IF2      = IF2 + 2
      IF1      = IF2 - 1
      
      DO J=1,NTot

            ! Apply a random phase to each of the columns of H to
            ! produce random phases in the wind component.
            ! Then sum each of the rows into the vector V.
            
         IRand = IFreq + (J-1)*NumFreq + (IVec-1)*NTot*NumFreq  ! This sorts the random numbers as they were done previously

         Ph    = TwoPi*RandNum(IRand)
         CPh   = COS( Ph )
         SPh   = SIN( Ph )

         Indx  = NTot*(J-1) - J*(J-1)/2 + J !Index of H(I,J)
         DO I=J,NTot

            V(IF1,I,IVec) = V(IF1,I,IVec) + TRH(Indx)*CPh  !Real part
            V(IF2,I,IVec) = V(IF2,I,IVec) + TRH(Indx)*SPh  !Imaginary part

            Indx = Indx + 1      !H(I,J)

         ENDDO ! I  
      ENDDO ! J

   ENDDO !IFreq

ENDDO !IVec

IF (COH_OUT)  CLOSE( UC )

IF ( ALLOCATED( Dist      ) ) DEALLOCATE( Dist      )
IF ( ALLOCATED( DistU     ) ) DEALLOCATE( DistU     )
IF ( ALLOCATED( DistZMExp ) ) DEALLOCATE( DistZMExp )


RETURN

END SUBROUTINE CohSpecVMat
!=======================================================================
SUBROUTINE GetChebCoefs(URef, RefHt)

   ! This subroutine determines what Chebyshev Coefficients will be used
   ! for the jet wind speed and wind direction profiles

USE                        TSMods

IMPLICIT                   NONE

REAL(ReKi)              :: UH_coef(4,11)     ! The coefficients that Neil developed for calculating the Chebyshev coefficients
REAL(ReKi)              :: WD_coef(4,11)     ! The coefficients that Neil developed for calculating the Chebyshev coefficients
REAL(ReKi)              :: ChebyCoef_tmp(11)
REAL(ReKi),INTENT(IN)   :: URef              ! The reference wind speed (i.e. target value at hub height)
REAL(ReKi)              :: UTmp1             !
REAL(ReKi)              :: UTmp2             !
REAL(ReKi),INTENT(IN)   :: RefHt             ! The height of the reference wind speed

INTEGER                 :: I                 ! A loop counter


      ! Let's calculate the wind speed at the jet height

   CALL get_coefs(ZJetMax, UH_coef, WD_coef)


   IF ( RefHt == ZJetMax ) THEN

      UJetMax = URef

      DO I=1,11
         ChebyCoef_WS(I) = UJetMax*UH_coef(1,I) + Rich_No*UH_coef(2,I) &
                           + Ustar*UH_coef(3,I) +         UH_coef(4,I)
      ENDDO

   ELSE

         ! Calculate the coefficients without UJetMax

      DO I=1,11
         ChebyCoef_WS(I) = Rich_No*UH_coef(2,I) + Ustar*UH_coef(3,I) + UH_coef(4,I) ! +UJetMax*UH_coef(1,I)
      ENDDO

      Utmp1              = getWindSpeed(URef, RefHt, RefHt, RotorDiameter, PROFILE='JET')

         ! Now calculate the coefficients with just UJetMax term

      ChebyCoef_tmp(:)   = ChebyCoef_WS(:)
      ChebyCoef_WS(:)    = UH_coef(1,:)

      Utmp2              = getWindSpeed(URef, RefHt, RefHt, RotorDiameter, PROFILE='JET')       ! This uses the ChebyCoef_WS values, & ignores the first 2 inputs
      UJetMax            = (Uref - Utmp1)/Utmp2

         ! Get the final coefficients, using the computed UJetMax
      ChebyCoef_WS(:)    = UJetMax*ChebyCoef_WS(:) + ChebyCoef_tmp(:)

   ENDIF

   DO I=1,11
      ChebyCoef_WD(I)    = UJetMax*WD_coef(1,I) + Rich_No*WD_coef(2,I) &
                           + Ustar*WD_coef(3,I) +         WD_coef(4,I)
   ENDDO

!print *, 'UJetMax wind speed at ', ZJetMax, ' m: ', UJetMax, 'm/s'
!Utmp1 = getWindSpeed(URef, RefHt, ZJetMax, 999.9, PROFILE='JET')
!print *, "Calc'd  wind speed at ", ZJetMax, ' m: ', Utmp1, 'm/s'

RETURN
END SUBROUTINE GetChebCoefs
!=======================================================================
SUBROUTINE GetDefaultCoh(WS,Ht)
   ! These numbers come from Neil's analysis

   USE                     TSMods, ONLY : InCDec
   USE                     TSMods, ONLY : InCohB
   USE                     TSMods, ONLY : RICH_NO
   USE                     TSMods, ONLY : TurbModel
!   USE                     TSMods, ONLY : UHub

!   REAL(ReKi), PARAMETER                :: a =  0.007697495  !coeffs for WF_xxD best-fit equations
!   REAL(ReKi), PARAMETER                :: b =  0.451759656  !coeffs for WF_xxD best-fit equations
!   REAL(ReKi), PARAMETER                :: c =  6.559106387  !coeffs for WF_xxD best-fit equations
!   REAL(ReKi), PARAMETER                :: d = -0.10471942   !coeffs for WF_xxD best-fit equations
!   REAL(ReKi), PARAMETER                :: e = -1.19488521   !coeffs for WF_xxD best-fit equations
!   REAL(ReKi), PARAMETER                :: f =  0.005529328  !coeffs for WF_xxD best-fit equations
!   REAL(ReKi), PARAMETER                :: g =  0.059157163  !coeffs for WF_xxD best-fit equations
   REAL(ReKi)                           :: Coeffs(10,3)      ! coeffs for WS category coherence decrements
   REAL(ReKi), INTENT(IN)               :: Ht         !Height, usually hub height
   REAL(ReKi)                           :: Ht1        !Height, set to bounds of the individual models
   REAL(ReKi)                           :: Ht2        !Height squared
   REAL(ReKi)                           :: Ht3        !Height cubed
   REAL(ReKi), INTENT(IN)               :: WS         !Wind speed, usually = UHub
   REAL(ReKi)                           :: WS1        !Wind speed, set to bounds of individual models
   REAL(ReKi)                           :: RI1        !RICH_NO, set to bounds of individual models
   REAL(ReKi)                           :: RI2        !RICH_NO squared
   REAL(ReKi)                           :: RI3        !RICH_NO  cubed

   INTEGER                              :: I
   INTEGER                              :: Ri_Cat


      IF (RICH_NO <= 0.00 ) THEN
         IF ( RICH_NO <= - 1.0 ) THEN
            Ri_Cat = 1
         ELSE
            Ri_Cat = 2
         ENDIF
      ELSEIF (RICH_NO <= 0.25 ) THEN
         IF (RICH_NO <= 0.10 ) THEN
            Ri_Cat = 3
         ELSE
            Ri_Cat = 4
         ENDIF
      ELSE
            Ri_Cat = 5
      ENDIF

      SELECT CASE ( TurbModel )

         CASE ( 'GP_LLJ' )
            HT1 = MAX( 60.0, MIN( Ht, 100.0 ) )
            IF ( WS <= 14.0 ) THEN
               IF ( WS <= 8.0 ) THEN
                  IF     ( WS <= 6.0 ) THEN
                     coeffs(:,3) = (/  3.1322E+00,  2.2819E-03,  2.9214E+00, -5.2203E-04,  1.1877E+00, &
                                      -5.7605E-02,  3.7233E-06, -3.5021E-01, -1.7555E-03,  3.9712E-04 /)    !W  5
                     IF  ( WS <= 4.0 ) THEN !      WS <=  4
                        RI1 = MAX( 0.0, MIN( RICH_NO, 1.0 ) )
                        coeffs(:,1) = (/  4.8350E+00, -4.0113E-02,  7.8134E+00, -2.0069E-05, -1.9518E-01, &
                                         -1.4009E-01,  2.3195E-06,  8.2029E-02, -7.4979E-04,  6.1186E-04 /) !U  3
                        coeffs(:,2) = (/  3.2587E+00, -5.9086E-02,  9.7426E+00,  5.7360E-04,  2.1274E-01, &
                                         -1.6398E-01, -8.3786E-07,  6.6896E-02, -3.5254E-03,  6.4833E-04 /) !V  3
                     ELSE                   !  4 < WS <=  6
                        RI1 = MAX( -0.5, MIN( RICH_NO, 1.0 ) )
                        coeffs(:,1) = (/  9.2474E+00, -4.9849E-02,  6.0887E+00, -5.9124E-04,  4.4312E-02, &
                                         -1.1966E-01,  5.2652E-06, -1.0373E-01,  4.0480E-03,  5.5761E-04 /) !U  5
                        coeffs(:,2) = (/  3.6355E+00,  1.7701E-02,  4.2165E+00, -5.8828E-04,  9.5592E-02, &
                                         -6.5313E-02,  3.3875E-06, -1.7981E-02, -1.6375E-03,  3.0423E-04 /) !V  5
                     ENDIF
                  ELSE                      ! 6  < WS <=  8
                     RI1 = MAX( -0.5, MIN( RICH_NO, 1.0 ) )
                     coeffs(:,1) = (/  1.1795E+01, -7.5393E-02,  9.5279E+00, -3.4922E-04, -5.8973E-01, &
                                      -1.6753E-01,  4.4267E-06,  2.1797E-01,  7.7887E-04,  7.4912E-04 /)    !U  7
                     coeffs(:,2) = (/  1.7730E+00,  9.6577E-02,  8.1310E+00, -1.2028E-03,  3.0145E-02, &
                                      -1.2282E-01,  4.6866E-06,  3.5748E-02, -2.9013E-03,  4.8368E-04 /)    !V  7
                     coeffs(:,3) = (/  9.1695E-01,  9.1488E-02,  6.7163E+00, -1.2938E-03,  1.0315E+00, &
                                      -1.1976E-01,  5.6039E-06, -2.0416E-01, -3.4698E-03,  6.0175E-04 /)    !W  7
                  ENDIF
               ELSE ! 8.0 < WS <= 14.0
                  IF     (WS <= 10.0) THEN  !  8 < WS <= 10
                     RI1 = MAX( -0.5, MIN( RICH_NO, 1.0 ) )
                     coeffs(:,1) = (/  8.4674E+00,  1.2922E-01,  8.6170E+00, -3.3048E-03, -3.1928E-02, &
                                      -1.2515E-01,  1.8209E-05,  2.9087E-01, -9.3031E-03,  5.0706E-04 /)    !U  9
                     coeffs(:,2) = (/  2.8145E+00,  1.0257E-01,  4.2987E+00, -1.4901E-03,  4.9698E-02, &
                                      -3.9964E-02,  6.7640E-06,  2.2980E-01, -1.0046E-02,  1.3037E-04 /)    !V  9
                     coeffs(:,3) = (/  2.4952E+00,  5.8000E-02,  1.9851E+00, -9.4027E-04, -4.0135E-02, &
                                      -1.8377E-02,  4.3320E-06, -1.0441E-01,  3.6831E-03,  8.6637E-05 /)    !W  9
                  ELSEIF (WS <= 12.0) THEN  ! 10 < WS <= 12
                     RI1 = MAX( -0.5, MIN( RICH_NO, 1.0 ) )
                     coeffs(:,1) = (/  1.2473E+01,  3.2270E-02,  1.4508E+01, -2.2856E-03, -1.4652E+00, &
                                      -2.4114E-01,  1.4919E-05,  5.5578E-01, -8.5528E-04,  1.0273E-03 /)    !U  11
                     coeffs(:,2) = (/  1.0882E+00,  1.9425E-01,  8.1533E+00, -2.5574E-03,  4.3113E-01, &
                                      -8.0465E-02,  1.0478E-05,  1.1640E-01, -1.1717E-02,  1.6476E-04 /)    !V  11
                     coeffs(:,3) = (/  5.0280E-01,  1.1637E-01,  4.0130E+00, -1.2034E-03, -2.7592E-01, &
                                      -3.8744E-02,  3.4213E-06, -1.5144E-02,  2.4042E-03,  4.7818E-05 /)    !W  11
                  ELSE                      ! 12 < WS <= 14.0
                     RI1 = MAX( -1.0, MIN( RICH_NO, 1.0 ) )
                     coeffs(:,1) = (/  8.6311E+00,  2.5614E-01,  1.1165E+01, -5.1685E-03,  3.0895E+00, &
                                      -1.9190E-01,  2.7162E-05, -2.6513E-01, -3.6479E-02,  8.8431E-04 /)    !U  13
                     coeffs(:,2) = (/  1.2842E+00,  2.4007E-01,  5.3653E+00, -3.2589E-03,  3.4715E+00, &
                                      -6.8865E-02,  1.3756E-05, -4.8465E-01, -4.0608E-02,  3.8578E-04 /)    !V  13
                     coeffs(:,3) = (/  4.3681E+00,  1.2251E-02,  1.3826E+00, -1.1592E-04,  3.3654E+00, &
                                      -5.2367E-02, -4.4086E-08, -3.5254E-01, -1.6780E-02,  3.9048E-04 /)    !W  13
                  ENDIF
               ENDIF
            ELSE ! WS > 14
               IF (WS <= 20.0 ) THEN
                  IF     (WS <= 16.0) THEN  ! 14 < WS <= 16
                     RI1 = MAX( -1.0, MIN( RICH_NO, 1.0 ) )
                     coeffs(:,1) = (/  1.3972E-01,  6.3486E-01,  1.7576E+01, -1.0017E-02,  2.8458E+00, &
                                      -2.5233E-01,  4.6539E-05, -1.8899E-01, -2.6717E-02,  9.5173E-04 /)    !U  15
                     coeffs(:,2) = (/ -7.1243E+00,  5.6768E-01,  1.2886E+01, -7.3277E-03,  3.7880E+00, &
                                      -1.4733E-01,  3.0898E-05, -1.5056E-01, -2.9500E-02,  3.6703E-04 /)    !V  15
                     coeffs(:,3) = (/ -1.1004E+01,  5.3470E-01,  5.3118E+00, -5.8999E-03,  1.9009E+00, &
                                      -2.4063E-02,  2.1755E-05, -4.5798E-01,  1.6885E-02, -3.9974E-04 /)    !W  15
                  ELSEIF (WS <= 18.0) THEN  ! 16 < WS <= 18
                     RI1 = MAX( -0.5, MIN( RICH_NO, 1.0 ) )
                     coeffs(:,1) = (/ -6.9650E+00,  8.8636E-01,  2.3467E+01, -1.1973E-02, -4.3750E+00, &
                                      -3.5519E-01,  5.0414E-05,  9.1789E-01,  9.8340E-03,  1.5885E-03 /)    !U  17
                     coeffs(:,2) = (/  5.5495E-03,  3.2906E-01,  1.4609E+01, -4.1635E-03, -2.1246E+00, &
                                      -1.8887E-01,  1.6964E-05,  3.7805E-01,  1.1880E-03,  8.8265E-04 /)    !V  17
                     coeffs(:,3) = (/ -1.3195E+00,  2.0022E-01,  2.3490E+00, -2.1308E-03,  3.5582E+00, &
                                       1.4379E-02,  7.6830E-06, -7.6155E-01, -2.4660E-02, -2.0199E-04 /)    !W  17
                  ELSE                      ! 18 < WS <= 20
                     RI1 = MAX( -0.5, MIN( RICH_NO, 1.0 ) )
                     coeffs(:,1) = (/ -1.3985E+01,  1.3161E+00,  3.4773E+01, -1.9237E-02, -1.9845E+00, &
                                      -5.5817E-01,  8.8310E-05,  1.7142E+00, -4.2907E-02,  2.3932E-03 /)    !U  19
                     coeffs(:,2) = (/ -1.2400E+01,  8.6854E-01,  1.9923E+01, -1.1557E-02, -1.0441E+00, &
                                      -2.4593E-01,  4.9813E-05,  2.7861E-01, -8.6189E-03,  9.4314E-04 /)    !V  19
                     coeffs(:,3) = (/ -9.3436E+00,  6.4950E-01,  1.5316E+01, -8.7208E-03,  1.7329E+00, &
                                      -2.2411E-01,  3.6288E-05, -8.0006E-01, -2.6439E-03,  7.9293E-04 /)    !W  19
                  ENDIF
               ELSE ! WS > 20
                  IF     (WS <= 22.0) THEN  ! 20 < WS <= 22
                     RI1 = MAX( -0.5, MIN( RICH_NO, 1.0 ) )
                     coeffs(:,1) = (/ -2.4317E+01,  1.8176E+00,  5.3359E+01, -2.5973E-02,  6.0349E+00, &
                                      -7.9927E-01,  1.1558E-04,  1.5926E+00, -1.5005E-01,  3.1688E-03 /)    !U  21
                     coeffs(:,2) = (/  8.0459E+00,  1.8058E-01,  1.9426E+01, -3.6730E-03, -9.9717E-01, &
                                      -1.8249E-01,  1.9237E-05,  4.9173E-01, -1.8255E-02,  6.9371E-04 /)    !V  21
                     coeffs(:,3) = (/ -2.3544E+01,  1.1403E+00,  8.3526E+00, -1.4511E-02,  7.2014E+00, &
                                       5.0216E-02,  5.9947E-05, -1.0659E+00, -7.4769E-02, -9.8390E-04 /)    !W  21
                  ELSEIF (WS <= 24.0) THEN  ! 22 < WS <= 24
                     RI1 = MAX( 0.0, MIN( RICH_NO, 1.0 ) )
                     coeffs(:,1) = (/ -3.5790E+01,  1.5374E+00,  1.1322E+02, -1.6884E-02, -1.7767E+01, &
                                      -1.8122E+00,  6.8247E-05,  7.2101E+00,  3.5536E-02,  7.9269E-03 /)    !U  23
                     coeffs(:,2) = (/ -7.2883E+01,  2.8210E+00,  8.6392E+01, -3.1084E-02, -2.4938E+01, &
                                      -1.5898E+00,  1.0997E-04,  7.1972E+00,  1.2624E-01,  9.3084E-03 /)    !V  23
                     coeffs(:,3) = (/ -3.2844E+01,  1.2683E+00,  3.2032E+01, -1.3197E-02, -1.1129E+01, &
                                      -3.6741E-01,  4.2852E-05,  4.1336E+00,  2.4775E-02,  1.8431E-03 /)    !W  23
                  ELSE                      ! 24 < WS
                     RI1 = MAX( -0.5, MIN( RICH_NO, 1.0 ) )
                     coeffs(:,1) = (/  2.2906E+01,  9.3209E-02,  1.5448E+01, -5.7421E-03, -8.9114E+00, &
                                      -3.1547E-02,  4.0144E-05,  5.4544E-01,  5.3557E-02, -3.1299E-04 /)    !U  25
                     coeffs(:,2) = (/ -1.1903E+01,  1.1104E+00,  1.7962E+01, -1.6045E-02, -9.2458E+00, &
                                      -4.4526E-02,  6.9880E-05,  2.8017E+00, -2.7211E-02, -8.4099E-04 /)    !V  25
                     coeffs(:,3) = (/  6.1054E-01,  7.1841E-03,  4.2996E+00,  2.9071E-04, -2.0002E+00, &
                                      -7.0403E-02, -2.8931E-06,  2.3943E-02,  1.8395E-02,  5.0406E-04 /)    !W  25
                  ENDIF
               ENDIF
            ENDIF


            HT2 = HT1*HT1
            HT3 = HT1*HT2
            RI2 = RI1*RI1
            RI3 = RI1*RI2

            DO I = 1,3
               InCDec(I) = coeffs( 1,I) + coeffs(2,I)*Ht1 + coeffs(3,I)*RI1 &
                                        + coeffs(4,I)*Ht2 + coeffs(5,I)*RI2 + coeffs( 6,I)*Ht1*RI1 &
                                        + coeffs(7,I)*Ht3 + coeffs(8,I)*RI3 + coeffs( 9,I)*Ht1*RI2 &
                                                                            + coeffs(10,I)*Ht2*RI1
            ENDDO

            WS1 = MAX(  2.0, WS )
            SELECT CASE ( Ri_Cat )
               CASE ( 1, 2)
!                 InCDec   = (/            1.744591004*WS1**0.593219225, &
!                              -0.58750092+1.937230512*WS1**0.400548383, &
!                              -0.57833219+1.450654739*WS1**0.443191083 /)
                  InCohB   = (/-0.00014115+0.006826264/WS1, &
                                           0.014025749/WS1, &
                               0.000480386+0.020982336/WS1 /)

               CASE ( 3 )
!                 InCDec   = (/            1.962126171*WS1**0.575523536, &
!                              -2.79495117+3.698342796*WS1**0.305415750, &
!                                          0.887573173*WS1**0.498317195 /)
                  InCohB   = (/-0.00016838+0.009764148/WS1, &
                                           0.018582932/WS1, &
                               0.001865953+0.061952454/WS1 /)

               CASE ( 4 )
!                 InCDec   = (/            0.817085986*WS1**1.045777184, &
!                                          0.599696362*WS1**1.038373995, &
!                                          1.327586050*WS1**0.590370871 /)
                  InCohB   = (/0.000175033+0.004195814/WS1, &
                                           0.008479460/WS1, &
                               0.002318082+0.027820652/WS1 /)

               CASE ( 5 )
!                 InCDec   = (/            0.959999473*WS1**0.972466847, &
!                              0.082701643+0.867230846*WS1**0.925895412, &
!                                          1.524380209*WS1**0.548060899 /)
                  InCohB   = (/0.000241808+0.004267702/WS1, &
                                           0.005408592/WS1, &
                               0.001150319+0.010744459/WS1 /)
               END SELECT


         CASE ( 'NWTCUP', 'USRVKM' )
            HT1 = MAX( 25.0, MIN( Ht, 50.0 ) )

            IF ( WS <= 14.0 ) THEN
               RI1 = MAX( -1.0, MIN( RICH_NO, 1.0 ) )
               IF ( WS <= 8.0 ) THEN
                  IF     (WS <= 4.0 ) THEN  !      WS <=  4
                     coeffs(:,1) = (/  8.1767E+00, -3.1018E-01,  3.3055E-01,  4.4232E-03,  4.6550E-01, &
                                      -2.4582E-02, -5.8568E-06, -8.7873E-02,  1.3070E-02,  3.1871E-04 /)   !U  3
                     coeffs(:,2) = (/  5.8003E+00, -2.0838E-01,  2.8727E-01,  2.8669E-03,  6.9669E-01, &
                                      -8.2249E-03, -2.4732E-06, -1.0826E-01,  9.9973E-03,  1.8546E-05 /)   !V  3
                     coeffs(:,3) = (/  5.9625E+00, -2.9247E-01, -9.3269E-01,  4.4089E-03,  1.3779E-01, &
                                       2.6993E-02, -6.1784E-06, -7.2920E-02,  1.7028E-02, -3.3753E-04 /)   !W  3
                  ELSEIF (WS <= 6.0 ) THEN  !  4 < WS <=  6
                     coeffs(:,1) = (/  1.2891E+01, -4.8265E-01,  3.5549E+00,  6.6099E-03,  8.2275E-01, &
                                      -1.5913E-01, -7.9740E-06, -1.2357E-02,  3.2084E-03,  1.7145E-03 /)   !U  5
                     coeffs(:,2) = (/  8.0267E+00, -2.5275E-01,  1.3801E+00,  3.2447E-03,  1.6004E+00, &
                                      -3.2592E-02, -5.1265E-06, -9.8552E-02, -1.3513E-02,  2.8075E-04 /)   !V  5
                     coeffs(:,3) = (/  7.9593E+00, -3.6336E-01,  1.4974E+00,  5.4012E-03,  9.5041E-01, &
                                      -1.0152E-01, -1.0865E-05,  4.3121E-02, -3.2447E-03,  1.3797E-03 /)   !W  5
                  ELSE                      ! 6  < WS <=  8
                     coeffs(:,1) = (/  1.3702E+01, -4.4674E-01,  3.7943E+00,  5.9350E-03,  9.6026E-01, &
                                      -1.7425E-01, -7.2917E-06, -8.8426E-02,  5.1530E-03,  2.0554E-03 /)   !U  7
                     coeffs(:,2) = (/  9.2471E+00, -2.6247E-01,  1.4504E+00,  3.2436E-03,  1.8823E+00, &
                                      -3.2180E-02, -5.9491E-06, -2.0100E-01, -1.7619E-02,  3.8519E-04 /)   !V  7
                     coeffs(:,3) = (/  8.9439E+00, -3.8885E-01,  2.2175E+00,  5.6207E-03,  7.6040E-01, &
                                      -1.3502E-01, -9.2514E-06,  1.9269E-02,  3.8862E-03,  1.7674E-03 /)   !W  7
                  ENDIF
               ELSE ! 8.0 < WS <= 14.0
                  IF     (WS <= 10.0) THEN  !  8 < WS <= 10
                     coeffs(:,1) = (/  1.9061E+01, -4.5354E-01,  7.5961E+00,  5.2422E-03,  1.5158E+00, &
                                      -2.4908E-01, -2.5277E-06, -1.6660E-01,  1.1369E-02,  3.0156E-03 /)   !U  9
                     coeffs(:,2) = (/  1.3362E+01, -3.3806E-01,  7.0401E+00,  4.5349E-03,  2.6798E+00, &
                                      -2.3637E-01, -9.9075E-06, -2.2373E-01, -1.6644E-03,  2.3879E-03 /)   !V  9
                     coeffs(:,3) = (/  8.8401E+00, -2.9945E-01,  3.7883E+00,  4.4581E-03,  2.0417E+00, &
                                      -2.7852E-01, -7.0750E-06, -6.2618E-02,  1.4646E-02,  3.8512E-03 /)   !W  9
                  ELSEIF (WS <= 12.0) THEN  ! 10 < WS <= 12
                     coeffs(:,1) = (/  3.4011E+01, -1.2590E+00,  1.6320E+01,  1.9225E-02,  6.8346E+00, &
                                      -8.8950E-01, -6.2453E-05, -2.4945E-01, -4.3892E-02,  1.2078E-02 /)   !U  11
                     coeffs(:,2) = (/  1.7135E+01, -4.0754E-01,  1.0282E+01,  5.7832E-03,  6.3056E+00, &
                                      -2.8536E-01, -3.0216E-05, -5.3170E-01, -5.7090E-02,  2.8463E-03 /)   !V  11
                     coeffs(:,3) = (/  1.3002E+01, -4.8326E-01,  3.2819E+00,  7.8800E-03,  2.7094E+00, &
                                      -2.5714E-01, -3.0117E-05, -2.1404E-01, -4.2711E-03,  4.1067E-03 /)   !W  11
                  ELSE                      ! 12 < WS <= 14
                     coeffs(:,1) = (/  2.6682E+01, -9.7229E-01,  1.3191E+01,  1.7604E-02, -1.3537E+00, &
                                      -6.4082E-01, -7.8242E-05,  1.7548E-01,  9.7417E-02,  1.0259E-02 /)   !U  13
                     coeffs(:,2) = (/  1.7083E+01, -4.7346E-01,  1.3515E+01,  7.7832E-03,  5.8633E-01, &
                                      -6.1815E-01, -3.3752E-05, -1.7300E-01,  4.3584E-02,  8.9289E-03 /)   !V  13
                     coeffs(:,3) = (/  1.6015E+01, -6.3912E-01,  1.3137E+01,  9.4757E-03,  2.5549E+00, &
                                      -8.1438E-01, -1.5565E-05,  2.9244E-02,  2.2779E-02,  1.1982E-02 /)   !W  13
                  ENDIF
               ENDIF
            ELSE ! WS > 14
               IF (WS <= 20.0 ) THEN
                  IF     (WS <= 16.0) THEN  ! 14 < WS <= 16
                     RI1 = MAX( -1.0, MIN( RICH_NO, 1.0 ) )
                     coeffs(:,1) = (/  2.9459E+01, -7.3181E-01,  9.4613E+00,  9.2172E-03,  6.1086E+00, &
                                      -4.9990E-01, -2.9994E-05, -6.9606E-01, -8.5076E-03,  8.1330E-03 /)   !U  15
                     coeffs(:,2) = (/  1.7540E+01, -2.6071E-01,  9.3639E+00,  1.3341E-03,  9.4294E+00, &
                                      -4.2565E-01, -2.7836E-06, -6.7708E-01, -6.9127E-02,  6.2290E-03 /)   !V  15
                     coeffs(:,3) = (/  1.2792E+01, -4.6469E-01,  4.6350E+00,  1.0633E-02,  1.8523E+00, &
                                      -3.2417E-01, -8.5038E-05, -2.2253E-01, -7.3351E-04,  5.4781E-03 /)   !W  15
                  ELSEIF (WS <= 18.0) THEN  ! 16 < WS <= 18
                     RI1 = MAX( -1.0, MIN( RICH_NO, 1.0 ) )
                     coeffs(:,1) = (/  1.7775E+01,  4.5287E-01,  1.6417E+01, -2.3724E-02,  5.8998E+00, &
                                      -5.3502E-01,  2.6202E-04, -9.9466E-02,  4.1386E-02,  4.5663E-03 /)   !U  17
                     coeffs(:,2) = (/  1.2022E+01,  2.4246E-01,  1.3875E+01, -1.1725E-02,  5.1917E+00, &
                                      -5.4329E-01,  1.1893E-04, -2.0308E-01,  6.5256E-02,  5.6597E-03 /)   !V  17
                     coeffs(:,3) = (/  1.2680E+01, -1.4768E-01,  7.1498E+00, -3.0341E-03,  1.9747E+00, &
                                      -3.8374E-01,  7.0412E-05,  2.2297E-01,  5.9943E-02,  5.3514E-03 /)   !W  17
                  ELSE                      ! 18 < WS <= 20
                     RI1 = MAX( -0.5, MIN( RICH_NO, 1.0 ) )
                     coeffs(:,1) = (/  3.1187E+01, -6.8540E-01,  7.1288E+00,  1.1923E-02,  8.8547E+00, &
                                       6.3133E-02, -9.4673E-05, -2.5710E+00, -5.4077E-02, -1.2797E-04 /)   !U  19
                     coeffs(:,2) = (/  1.2664E+01,  9.1858E-02,  1.9050E+01, -2.8868E-03,  7.2969E+00, &
                                      -4.4573E-01, -6.1033E-06, -2.0960E+00, -1.9913E-02,  4.9023E-03 /)   !V  19
                     coeffs(:,3) = (/  2.2146E+01, -7.6940E-01,  1.1948E+01,  1.0400E-02,  5.0034E+00, &
                                      -4.3958E-01, -2.5936E-05, -3.0848E-01, -6.3381E-02,  5.1204E-03 /)   !W  19
                  ENDIF
               ELSE ! WS > 20
                  RI1 = MAX( -0.5, MIN( RICH_NO, 1.0 ) )
                  IF     (WS <= 22.0) THEN  ! 20 < WS <= 22
                     coeffs(:,1) = (/  2.5165E+01, -7.7660E-02,  1.9692E+01, -1.1794E-02,  9.8635E+00, &
                                      -2.5520E-01,  2.0573E-04, -4.9850E+00,  1.1272E-01,  1.3267E-03 /)   !U  21
                     coeffs(:,2) = (/  2.1691E+01, -3.1787E-01,  3.2327E+01, -4.5546E-03,  1.1194E+01, &
                                      -8.0823E-01,  1.4306E-04, -4.3418E+00,  7.3163E-02,  6.3637E-03 /)   !V  21
                     coeffs(:,3) = (/  1.4634E+01, -3.9394E-01,  1.1617E+01,  5.6387E-03,  5.4799E+00, &
                                      -3.9011E-01, -1.0420E-05, -2.4279E+00,  6.6452E-02,  4.9504E-03 /)   !W  21
                  ELSEIF (WS <= 24.0) THEN  ! 22 < WS <= 24
                     coeffs(:,1) = (/  7.3816E+00,  1.0538E+00,  2.1578E+01, -3.3487E-02, -6.4986E+00, &
                                      -8.6782E-01,  3.2397E-04,  1.1412E+00,  2.2982E-01,  1.4660E-02 /)   !U  23
                     coeffs(:,2) = (/  6.5302E+00,  1.0524E+00,  2.4596E+01, -4.1648E-02,  4.0584E+00, &
                                      -6.1130E-01,  4.5468E-04, -3.6547E+00,  2.3176E-01,  8.4385E-03 /)   !V  23
                     coeffs(:,3) = (/  1.3424E+01,  2.6104E-02,  7.6014E+00, -1.2744E-02,  1.0735E+01, &
                                       2.2086E-01,  1.9309E-04, -5.9548E+00,  8.6483E-02, -3.9550E-03 /)   !W  23
                  ELSE                      ! 24 < WS
                     coeffs(:,1) = (/ -1.6629E+01,  1.3094E+00, -4.4183E+00, -8.4860E-03, -1.3800E+01, &
                                      -5.5221E-01, -5.6659E-05,  8.1834E+00, -8.2497E-03,  1.8383E-02 /)   !U  25
                     coeffs(:,2) = (/  3.4796E+00,  7.1144E-01,  1.2153E+01, -2.7309E-02,  1.0003E+00, &
                                      -6.3570E-01,  3.4424E-04, -8.5038E-01,  1.2822E-01,  1.3181E-02 /)   !V  25
                     coeffs(:,3) = (/  2.7014E+00,  1.1794E-01,  2.1378E+00,  4.5539E-03,  1.6899E+00, &
                                       1.2254E-01, -9.6940E-05, -2.3430E-01, -2.3826E-02,  5.5964E-05 /)   !W  25
                  ENDIF
               ENDIF
            ENDIF

            HT2 = HT1*HT1
            HT3 = HT1*HT2
            RI2 = RI1*RI1
            RI3 = RI1*RI2

            DO I = 1,3
               InCDec(I) = coeffs( 1,I) + coeffs(2,I)*Ht1 + coeffs(3,I)*RI1 &
                                        + coeffs(4,I)*Ht2 + coeffs(5,I)*RI2 + coeffs( 6,I)*Ht1*RI1 &
                                        + coeffs(7,I)*Ht3 + coeffs(8,I)*RI3 + coeffs( 9,I)*Ht1*RI2 &
                                                                            + coeffs(10,I)*Ht2*RI1
            ENDDO

            WS1 = MAX(  2.0, WS )
            SELECT CASE ( Ri_Cat )
               CASE ( 1 )
!                 InCDec   = (/            1.623224368*WS1**1.015099356, &
!                                          0.884720872*WS1**1.192553093, &
!                                          1.338245093*WS1**0.841757461 /)
                  InCohB   = (/ -2.524e-05+0.002122544/WS1, &
                                           0.004367773*WS1**(-1.14945936), &
                                           0.031284497*WS1**(-0.72509517) /)

               CASE ( 2 )
!                 InCDec   = (/            1.478475074*WS1**0.752442176, &
!                                          1.310684825*WS1**0.624122449, &
!                                          0.849106068*WS1**0.627688235 /)
                  InCohB   = (/            0.003320615*WS1**(-1.18592214), &
                                           0.005402681*WS1**(-0.98637053), &
                                           0.091649927*WS1**(-1.48835650) /)

               CASE ( 3 )
!                 InCDec   = (/            1.596175944*WS1**0.674743966, &
!                                          1.114069218*WS1**0.638049141, &
!                                          0.473225245*WS1**0.784331891 /)
                  InCohB   = (/            0.002387997*WS1**(-0.85956868), &
                                           0.009481901*WS1**(-1.02518835), &
                                           0.052147706*WS1**(-0.88949864) /)

               CASE ( 4 )
!                 InCDec   = (/            1.293345620*WS1**0.955639280, &
!                                          1.296399839*WS1**0.838281755, &
!                                          0.333750239*WS1**1.103784094 /)
                  InCohB   = (/            0.002870978*WS1**(-1.07398490), &
                                           0.002435238*WS1**(-0.68685045), &
                                           0.125356016*WS1**(-1.34791890) /)

               CASE ( 5 )
!                 InCDec   = (/            1.325256941*WS1**1.039629269, &
!                                          1.014004299*WS1**1.082810576, &
!                                          0.206383058*WS1**1.435200799 /)
                  InCohB   = (/            0.003545043*WS1**(-1.03669585), &
                                           0.003996215*WS1**(-0.95313438), &
                                           0.125103070*WS1**(-1.02886635) /)
               END SELECT

         CASE ( 'WF_UPW' )
            HT1 = MAX( 5.0, MIN( Ht, 35.0 ) )
            IF ( WS <= 14.0 ) THEN
               IF ( WS <= 10 ) THEN
                  RI1 = MAX( -0.5, MIN( RICH_NO, 0.15 ) )
                  IF  ( WS <=  8.0 ) THEN   !      WS <= 8
                     coeffs(:,1) = (/  1.6715E+01, -3.8639E-01,  7.1817E+00,  1.5550E-03, -1.4293E+00, &
                                      -2.0350E-01,  8.5532E-06, -3.4710E+00, -1.9743E-02, -3.9949E-04 /) !Upw_U 7
                     coeffs(:,2) = (/  8.4145E+00, -4.7610E-02,  3.9097E+00, -7.1412E-04,  1.8295E+01, &
                                       2.2583E-01, -1.6965E-05,  2.0769E+01, -9.1670E-02, -8.0300E-03 /) !Upw_V 7
                  ELSE                      !  8 < WS <= 10
                     coeffs(:,1) = (/  1.5432E+01, -2.1254E-01,  5.3075E+00, -2.9928E-03,  2.1647E+00, &
                                       1.1787E-02,  6.7458E-05, -9.0445E-01, -7.5941E-02, -4.7053E-03 /) !Upw_U 9
                     coeffs(:,2) = (/  7.5921E+00,  3.3520E-02,  1.2231E+01, -7.0018E-03,  6.0889E+01, &
                                       2.1810E-01,  1.1718E-04,  7.7287E+01, -1.3828E-01, -9.6568E-03 /) !Upw_V 9
                  ENDIF
               ELSE
                  RI1 = MAX( -0.5, MIN( RICH_NO, 0.05 ) )
                  IF  ( WS <= 12.0 ) THEN   ! 10 < WS <= 12
                     coeffs(:,1) = (/  1.3539E+01, -8.4892E-02, -1.9237E+00, -1.1485E-03, -4.0840E-01, &
                                       3.0956E-01,  2.4048E-05, -1.1523E+00,  9.6877E-03, -4.0606E-03 /) !Upw_U 11
                     coeffs(:,2) = (/  7.7451E+00, -1.3818E-01, -9.5197E-01,  3.9610E-03,  8.3255E-01, &
                                       7.2166E-02, -4.5012E-05, -2.0948E-01, -2.1400E-02, -2.9788E-04 /) !Upw_V 11
                  ELSE                      ! 12 < WS <= 14
                     coeffs(:,1) = (/  1.2857E+01, -7.9408E-03, -1.5310E+00, -4.1077E-03,  1.0496E+00, &
                                       1.9473E-01,  7.2808E-05,  1.8380E-01, -1.6559E-02, -2.0872E-03 /) !Upw_U 13
                     coeffs(:,2) = (/  7.2452E+00, -6.2662E-02, -2.4865E+00,  3.2123E-03, -1.0281E-01, &
                                       1.9698E-01, -7.5745E-05, -1.1637E+00, -4.6458E-02, -2.7037E-03 /) !Upw_V 13
                  ENDIF
               ENDIF
            ELSE
               RI1 = MAX( -0.5, MIN( RICH_NO, 0.05 ) )
               IF  ( WS  <= 18.0 ) THEN
                  IF ( WS <= 16.0 ) THEN   ! 14 < WS <= 16
                     coeffs(:,1) = (/  1.4646E+01, -1.5023E-01, -9.7543E-01, -3.5607E-03,  4.8663E+00, &
                                      -9.4360E-03,  1.4932E-04,  5.9503E+00,  7.4028E-02,  5.2698E-03 /) !Upw_U 15
                     coeffs(:,2) = (/  1.0133E+01, -3.1417E-01,  2.5400E+00,  6.6777E-03,  3.0790E+00, &
                                      -2.5801E-01, -4.9501E-05,  2.8879E+00, -1.6722E-02,  4.8297E-03 /) !Upw_V 15
                  ELSE                     ! 16 < WS <= 18
                     coeffs(:,1) = (/  1.5282E+01, -2.7642E-01,  2.5903E+00,  9.8716E-03,  5.9314E-01, &
                                      -4.2790E-01, -1.6474E-04, -7.0065E-01, -3.2694E-02,  2.4583E-03 /) !Upw_U 17
                     coeffs(:,2) = (/  1.2464E+01, -3.4306E-01,  3.6261E+00,  5.8254E-03,  2.2592E+00, &
                                      -1.1498E-01, -6.6196E-05,  1.3610E+00, -1.3345E-02,  1.0932E-03 /) !Upw_V 17
                  ENDIF
               ELSE
                  IF ( WS <= 20.0 ) THEN   ! 18 < WS <= 20
                     coeffs(:,1) = (/  1.5059E+01, -8.0478E-02,  8.7088E+00, -1.7854E-03,  3.9922E+00, &
                                      -6.0268E-01,  4.3906E-05,  3.3463E+00, -6.6490E-02,  1.2290E-02 /) !Upw_U 19
                     coeffs(:,2) = (/  1.0672E+01, -2.8104E-01,  7.8021E+00,  6.6360E-03,  2.4345E+00, &
                                      -4.9103E-01, -8.3745E-05,  4.4084E-01, -9.2432E-02,  8.3096E-03 /) !Upw_V 19
                  ELSE                     ! 20 < WS
                     coeffs(:,1) = (/  1.8592E+01,  1.3888E-01,  1.6732E+01, -1.1880E-02,  2.3622E+01, &
                                       6.8199E-01,  7.3664E-05,  4.1289E+00, -3.8604E-01, -3.0381E-02 /) !Upw_U 21
                     coeffs(:,2) = (/  7.7137E+00,  1.2732E-01,  1.3477E+01,  1.9164E-03,  3.7133E+01, &
                                       3.8975E-01, -2.2818E-04,  1.8816E+01, -7.5304E-01, -2.1856E-02 /) !Upw_V 21
                  ENDIF
               ENDIF
            ENDIF

            HT2 = HT1*HT1
            HT3 = HT1*HT2
            RI2 = RI1*RI1
            RI3 = RI1*RI2

            DO I = 1,2
               InCDec(I) = coeffs( 1,I) + coeffs(2,I)*Ht1 + coeffs(3,I)*RI1 &
                                        + coeffs(4,I)*Ht2 + coeffs(5,I)*RI2 + coeffs( 6,I)*Ht1*RI1 &
                                        + coeffs(7,I)*Ht3 + coeffs(8,I)*RI3 + coeffs( 9,I)*Ht1*RI2 &
                                                                            + coeffs(10,I)*Ht2*RI1
            ENDDO

            WS1 = MAX(  3.0, WS )
!           InCDec(1:2)   = (/             5.640176786*WS1**0.269850341, &
!                              6.059554513+18.44124731/WS1**1.5 /)
            InCohB(1:2)   = (/ 0.000448295+0.002502915/WS1, &
                               0.001539069+0.005954785/WS1 /)


            InCDec(3)     =  0.4*InCDec(1)  !cohA(w) = cohA(u)/2.5, number derived from histograms of u/w for NWTC and LLLJP data
            InCohB(3)     = 10.0*InCohB(1)  !cohB(w) = cohB(u)*10, number derived from histograms of w/u for NWTC and LLLJP data

         CASE ( 'WF_07D', 'WF_14D' )
            HT1 = MAX( 5.0, MIN( Ht, 35.0 ) )
            IF ( WS <= 12.0 ) THEN
               IF     ( WS <=  8.0 ) THEN  !      WS <= 8
                  RI1 = MAX( -0.5, MIN( RICH_NO, 0.15 ) )
                  coeffs(:,1) = (/  1.0310E+01, -6.4824E-03, -1.3258E+00, -2.7238E-03, -6.8515E+00, &
                                    3.1602E-02,  5.5982E-05, -8.4777E+00,  2.1506E-02,  4.9745E-04 /) !Dwn_U 7
                  coeffs(:,2) = (/  6.9491E+00, -1.3378E-01,  1.7961E-01, -4.9439E-04, -1.8140E+00, &
                                   -4.2321E-02,  4.4962E-05, -3.6939E+00, -8.9465E-03,  4.7867E-04 /) !Dwn_V 7
               ELSEIF ( WS <= 10.0 ) THEN  !  8 < WS <= 10
                  RI1 = MAX( -0.5, MIN( RICH_NO, 0.05 ) )
                  coeffs(:,1) = (/  9.7420E+00,  6.1610E-02,  5.6636E-02, -5.5949E-03, -1.3014E+00, &
                                    2.0655E-01,  8.9989E-05, -1.9837E+00,  5.4957E-03, -3.5496E-03 /) !Dwn_U 9
                  coeffs(:,2) = (/  7.1063E+00, -1.7021E-01,  1.2560E+00, -4.2616E-04,  9.0937E-01, &
                                   -1.3022E-01,  4.7976E-05,  2.1302E-01, -4.3159E-04,  1.5443E-03 /) !Dwn_V 9
               ELSE                        ! 10 < WS <= 12
                  RI1 = MAX( -0.5, MIN( RICH_NO, 0.05 ) )
                  coeffs(:,1) = (/  1.0869E+01, -9.1393E-03, -1.1695E+00, -3.3725E-03,  3.2199E-01, &
                                    7.2692E-02,  7.0565E-05,  6.9573E-01,  2.5360E-02,  1.0187E-03 /) !Dwn_U 11
                  coeffs(:,2) = (/  6.9882E+00, -1.3517E-01, -3.0492E-01, -4.6775E-04,  4.6897E-01, &
                                   -2.0102E-03,  3.3908E-05,  1.4604E-02,  1.1729E-02, -6.2775E-05 /) !Dwn_V 11
               ENDIF
            ELSE
               RI1 = MAX( -0.5, MIN( RICH_NO, 0.05 ) )
               IF     ( WS <= 14.0 ) THEN  ! 12 < WS <= 14
                  coeffs(:,1) = (/  1.1105E+01,  5.3789E-02, -9.4253E-02, -5.4203E-03, -1.0114E+00, &
                                    1.1421E-01,  7.6110E-05, -1.2654E+00,  1.5121E-02, -2.9055E-03 /) !Dwn_U 13
                  coeffs(:,2) = (/  7.5741E+00, -8.3945E-02,  3.7020E+00, -6.0317E-03,  3.1339E-01, &
                                   -2.1921E-01,  1.5598E-04,  6.2478E-01,  5.9490E-02,  3.4785E-03 /) !Dwn_V 13
               ELSE                        ! 14 < WS
                  coeffs(:,1) = (/  1.2256E+01,  2.0131E-02,  1.9465E+00, -7.6608E-03,  1.5031E+00, &
                                   -1.0916E-01,  1.3634E-04,  1.3451E+00, -1.6458E-02,  3.8312E-03 /) !Dwn_U 15
                  coeffs(:,2) = (/  7.7749E+00, -2.2712E-01,  1.3675E+00,  6.7944E-03,  4.2033E-02, &
                                   -6.8887E-02, -9.6117E-05, -1.5526E+00, -2.2357E-02, -1.5311E-03 /) !Dwn_V 15
               ENDIF
            ENDIF

            HT2 = HT1*HT1
            HT3 = HT1*HT2
            RI2 = RI1*RI1
            RI3 = RI1*RI2

            DO I = 1,2
               InCDec(I) = coeffs( 1,I) + coeffs(2,I)*Ht1 + coeffs(3,I)*RI1 &
                                        + coeffs(4,I)*Ht2 + coeffs(5,I)*RI2 + coeffs( 6,I)*Ht1*RI1 &
                                        + coeffs(7,I)*Ht3 + coeffs(8,I)*RI3 + coeffs( 9,I)*Ht1*RI2 &
                                                                            + coeffs(10,I)*Ht2*RI1
            ENDDO

            WS1 = MAX(  3.0, WS )
!           WS2 = WS1*WS1
!           WS3 = WS2*WS1
!           InCDec(1:2)   = (/ (a+c*WS1+e*WS2+g*WS3)/(1+b*WS1+d*WS2+f*WS3), &
!                                               3.357892649*WS1**0.1198781 /)
            InCohB(1:2)   = (/ 4.49289e-05+0.004933460/WS1, &
                                0.00158053+0.014268899/WS1 /)
            InCDec(3)     =  0.4*InCDec(1)  !cohA(w) = cohA(u)/2.5, number derived from histograms of u/w for NWTC and LLLJP data
            InCohB(3)     = 10.0*InCohB(1)  !cohB(w) = cohB(u)*10, number derived from histograms of w/u for NWTC and LLLJP data

         CASE ( 'USRINP' )
            InCDec = (/   WS, HUGE(InCohB(1)), HUGE(InCohB(1)) /)
            InCohB = (/0.0  , 0.0            , 0.0             /)         

         CASE DEFAULT   ! includes CASE ( 'SMOOTH' )

            InCDec = (/   WS, 0.75*WS, 0.75*WS /)  ! The davenport exponential parameter indicates that coh(v) ~ coh(w) in NWTC and LLLJP data
            InCohB = (/0.0  , 0.0    , 0.0     /)

      END SELECT

END SUBROUTINE GetDefaultCoh
!=======================================================================
SUBROUTINE GetDefaultRS( UW, UV, VW )
   ! This subroutine is used to get the default values of the Reynolds
   !  stresses.

   USE                     TSMods, ONLY : GridHeight
   USE                     TSMods, ONLY : HubHt
   USE                     TSMods, ONLY : H_ref                           ! This is needed for the HYDRO spectral models.
   USE                     TSMods, ONLY : RICH_NO
   USE                     TSMods, ONLY : RotorDiameter
   USE                     TSMods, ONLY : TurbModel
   USE                     TSMods, ONLY : UHub
   USE                     TSMods, ONLY : Ustar
   USE                     TSMods, ONLY : UVskip                                   ! Flag to determine if UV cross-feed term should be skipped or used
   USE                     TSMods, ONLY : UWskip                                   ! Flag to determine if UW cross-feed term should be skipped or used
   USE                     TSMods, ONLY : VWskip                                   ! Flag to determine if VW cross-feed term should be skipped or used
   USE                     TSMods, ONLY : WindProfileType
   USE                     TSMods, ONLY : zL

   REAL(ReKi)                          :: Coef(11)    !coefficients for pdf in NWTCUP random residual calculation
   REAL(ReKi)                          :: FnRange(2)  !min/max values for range in NWTCUP random residual calculation
   REAL(ReKi)                          :: fMax        !max value for pdf in NWTCUP random residual calculation
   REAL(ReKi), INTENT(INOUT)           :: UW  ! PC_UW
   REAL(ReKi), INTENT(OUT)             :: UV  ! PC_UV
   REAL(ReKi), INTENT(OUT)             :: VW  ! PC_VW
   REAL(ReKi)                          :: rnd
   REAL(ReKi)                          :: rndSgn
   REAL(ReKi)                          :: SignProb
   REAL(ReKi)                          :: Shr
   REAL(ReKi)                          :: Ustar2
   REAL(ReKi)                          :: V(2)
   REAL(ReKi)                          :: Z(2)
   REAL(ReKi)                          :: ZLtmp

   Z(2) = HubHt + 0.5*RotorDiameter    ! top of the grid
   Z(1) = Z(2) - GridHeight            ! bottom of the grid
   V(:) = getWindSpeed(UHub, HubHt, Z, RotorDiameter, PROFILE=WindProfileType)

   Shr = ( V(2)-V(1) ) / GridHeight    ! dv/dz

!BJJ: check the ranges of our best-fit parameters, using domains of measured values

   SELECT CASE ( TurbModel )
      CASE ( 'GP_LLJ' )
         ZLtmp  = MIN( MAX( ZL,    REAL(-1.00,ReKi) ), REAL(1.0,ReKi) )  !Limit the observed values of z/L
         UStar2 = MIN( MAX( UStar, REAL( 0.15,ReKi) ), REAL(1.0,ReKi) )  !Limit the observed values of u*
         Ustar2 = Ustar2*Ustar2
      CASE ( 'NWTCUP' )
         ZLtmp  = MIN( MAX( ZL,    REAL(-0.5,ReKi) ), REAL(3.5,ReKi) )  !Limit the observed values of z/L
         UStar2 = MIN( MAX( UStar, REAL( 0.2,ReKi) ), REAL(1.4,ReKi) )  !Limit the observed values of u*
         Ustar2 = Ustar2*Ustar2
!      CASE ( 'WF_UPW' )
!      CASE ( 'WF_07D' )
!      CASE ( 'WF_14D' )

      CASE DEFAULT
         ZLtmp  = ZL
         Ustar2 = UStar*Ustar
   END SELECT
   
   !-------------------------------------------------------------------------------------------------
   ! default UW Reynolds stress
   !-------------------------------------------------------------------------------------------------

   CALL  RndUnif( rndSgn )
   SELECT CASE ( TRIM(TurbModel) )

      CASE ( 'GP_LLJ' )

         IF (UW <= 0) THEN  !We don't have a local u* value to tie it to; otherwise, assume UW contains magnitude of value we want
            IF ( HubHt >= 100.5 ) THEN     ! 116m
               UW =  0.0399 - 0.00371*Uhub - 0.00182*RICH_NO + 0.00251*ZLtmp - 0.402*Shr + 1.033*Ustar2
            ELSEIF ( HubHt >= 76.0 ) THEN  ! 85 m
               UW = 0.00668 - 0.00184*Uhub + 0.000709*RICH_NO  + 0.264*Shr + 1.065*Ustar2  !magnitude
            ELSEIF ( HubHt >= 60.5 ) THEN  ! 67 m
               UW = -0.0216 + 0.00319*Uhub  - 0.00205*ZLtmp + 0.206*Shr + 0.963*Ustar2    !magnitude
            ELSE                           ! 54 m
               UW = -0.0373 + 0.00675*Uhub  - 0.00277*ZLtmp + 0.851*Ustar2                !magnitude
            ENDIF
            UW = MAX(UW,0.0)
            
         ENDIF
                     
         IF (UW > 0) THEN
            SignProb = 0.765 + 0.57/PI * ATAN( 0.78511*LOG(UW)+3.42584)
            IF (rndSgn <= SignProb) UW = -UW
         ENDIF

      CASE ( 'NWTCUP' )

         IF ( HubHt > 47.0 ) THEN      ! 58m data
            UW = 0.165 - 0.0232*UHub - 0.0129*RICH_NO + 1.337*Ustar2 - 0.758*SHR
         ELSEIF ( HubHt >= 26.0 ) THEN ! 37m data
            UW = 0.00279 - 0.00139*UHub + 1.074*Ustar2 + 0.179*SHR
         ELSE                          ! 15m data
            UW = -0.1310 + 0.0239*UHub + 0.556*Ustar2
         ENDIF
         UW = MAX(UW,0.0)
         
         IF (UW > 0) THEN !i.e. not equal to zero
            SignProb = 0.765 + 0.57/PI * ATAN( 0.88356*LOG(UW)+2.47668)
            IF (rndSgn <= SignProb) UW = -UW
         ENDIF

      CASE ( 'WF_14D' )

         UW = -Ustar2
         IF ( rndSgn > 0.9937 )  UW = -UW

      CASE ( 'USRINP' )
         UW = 0.0
         UWskip = .true.

      CASE ( 'TIDAL', 'RIVER' ) ! HYDROTURBSIM specific.
         UW = -Ustar2*(1-HubHt/H_ref)
      CASE DEFAULT

         UW = -Ustar2

   END SELECT

   !-------------------------------------------------------------------------------------------------
   ! default UV Reynolds stress
   !-------------------------------------------------------------------------------------------------

   CALL  RndUnif( rndSgn )
   SELECT CASE ( TurbModel )

      CASE ( 'GP_LLJ' )

         IF ( HubHt >= 100.5 ) THEN     ! 116m
            UV = 0.199 - 0.0167*Uhub + 0.0115*ZLtmp + 1.143*Ustar2
            UV = MAX(UV,0.0)
            IF ( rndSgn < 0.6527 ) UV = -UV
         ELSEIF ( HubHt >= 76.0 ) THEN  ! 85 m
            UV = 0.190 - 0.0156*Uhub + 0.00931*ZLtmp + 1.101*Ustar2
            UV = MAX(UV,0.0)
            IF ( rndSgn < 0.6394 ) UV = -UV
         ELSEIF ( HubHt >= 60.5 ) THEN  ! 67 m
            UV = 0.178 - 0.0141*Uhub + 0.00709*ZLtmp + 1.072*Ustar2
            UV = MAX(UV,0.0)
            IF ( rndSgn < 0.6326 ) UV = -UV
         ELSE                           ! 54 m
            UV = 0.162 - 0.0123*Uhub + 0.00784*RICH_NO + 1.024*Ustar2
            UV = MAX(UV,0.0)
            IF ( rndSgn < 0.6191 ) UV = -UV
         ENDIF

      CASE ( 'NWTCUP' )

            ! Get the magnitude and add the sign
         IF ( HubHt > 47.0 ) THEN      ! 58m data
            UV = 0.669 - 0.0300*UHub - 0.0911*RICH_NO + 1.421*Ustar2 - 1.393*SHR
         ELSEIF ( HubHt >= 26.0 ) THEN ! 37m data
            UV = 1.521 - 0.00635*UHub - 0.2200*RICH_NO + 3.214*Ustar2 - 3.858*SHR
         ELSE                          ! 15m data
            UV = 0.462 - 0.01400*UHub + 1.277*Ustar2
         ENDIF
         UV = MAX(UV,0.0)
         IF (UV > 0) THEN !i.e. not equal to zero
            SignProb = 0.33 + 0.64/PI * ATAN( -0.374775*LOG(UV)-0.205681)
            IF (rndSgn <= SignProb) UV = -UV
         ENDIF

      CASE ( 'WF_UPW' )

         UV = 0.0202 + 0.890*Ustar2 - 2.461*Shr
         UV = MAX(UV,0.0)
         IF ( rndSgn < 0.7315 ) UV = -UV

      CASE ( 'WF_07D' )

         UV = 0.5040 + 0.177*Ustar2
         UV = MAX(UV,0.0)
         IF ( rndSgn < 0.7355 ) UV = -UV

      CASE ( 'WF_14D' )

         UV = 0.0430 + 0.258*Ustar2
         UV = MAX(UV,0.0)
         IF ( rndSgn < 0.4423 ) UV = -UV

      CASE DEFAULT

         UV  = 0.0
         UVskip = .TRUE.  !use whatever comes our way from the random phases

   END SELECT


   !-------------------------------------------------------------------------------------------------
   ! default VW Reynolds stress
   !-------------------------------------------------------------------------------------------------

   CALL  RndUnif( rndSgn )
   SELECT CASE ( TurbModel )

      CASE ( 'GP_LLJ' )

         IF ( HubHt >= 100.5 ) THEN     ! 116m
            VW =  0.0528  - 0.00210*Uhub - 0.00531*RICH_NO - 0.519*Shr + 0.283*Ustar2
            VW = MAX(VW,0.0)
            IF ( rndSgn < 0.2999 ) VW = -VW
         ELSEIF ( HubHt >= 76.0 ) THEN  ! 85 m
            VW =  0.0482  - 0.00264*Uhub - 0.00391*RICH_NO - 0.240*Shr + 0.265*Ustar2
            VW = MAX(VW,0.0)
            IF ( rndSgn < 0.3061 ) VW = -VW
         ELSEIF ( HubHt >= 60.5 ) THEN  ! 67 m
            VW =  0.0444  - 0.00249*Uhub - 0.00403*RICH_NO - 0.141*Shr + 0.250*Ustar2
            VW = MAX(VW,0.0)
            IF ( rndSgn < 0.3041 ) VW = -VW
         ELSE                           ! 54 m
            VW =  0.0443  - 0.00261*Uhub - 0.00371*RICH_NO - 0.107*Shr + 0.226*Ustar2
            VW = MAX(VW,0.0)
            IF ( rndSgn < 0.3111 ) VW = -VW
         ENDIF

      CASE ( 'NWTCUP' )

         IF ( HubHt > 47.0 ) THEN      ! 58m data
            VW = 0.174 + 0.00154*UHub - 0.0270*RICH_NO + 0.380*Ustar2 - 1.131*Shr - 0.00741*ZLtmp
         ELSEIF ( HubHt >= 26.0 ) THEN ! 37m data
            VW = 0.120 + 0.00283*UHub - 0.0227*RICH_NO + 0.306*Ustar2 - 0.825*Shr
         ELSE                          ! 15m data
            VW = 0.0165 + 0.00833*UHub                 + 0.224*Ustar2
         ENDIF
         VW = MAX(VW,0.0)
         IF (VW > 0) THEN !i.e. not equal to zero
            SignProb = 0.725 + 0.65/PI * ATAN( 0.654886*LOG(VW)+1.777198)
            IF (rndSgn <= SignProb) VW = -VW
         ENDIF

      CASE ( 'WF_UPW' )

         VW = 0.0263 + 0.273*Ustar2 - 0.684*Shr
         VW = MAX(VW,0.0)
         IF ( rndSgn < 0.3139 ) VW = -VW

      CASE ( 'WF_07D' )

         VW = 0.241 + 0.118*Ustar2
         VW = MAX(VW,0.0)
         IF ( rndSgn < 0.0982 ) VW = -VW

      CASE ( 'WF_14D' )

         VW =-0.0224 + 0.159*Ustar2
         VW = MAX(VW,0.0)
         IF ( rndSgn < 0.8436 ) VW = -VW

      CASE DEFAULT

         VW  = 0.0
         VWskip = .TRUE.  !use whatever comes our way from the random phases

   END SELECT


RETURN
END SUBROUTINE GetDefaultRS
!=======================================================================
SUBROUTINE GetFiles

  ! This subroutine is used to open the summary output file.

USE              TSMods

IMPLICIT         NONE


CALL GetRoot( InFile, RootName )


   ! Open summary file.

CALL OpenFOutFile( US, TRIM( RootName )//'.sum' ) ! Formatted output file


   ! Write the program name and version, date and time into the summary file.

   ! Let's make sure the binary file and the full-field file have the same date and time.
DescStr = 'generated by '//TRIM( ProgName )//TRIM( ProgVer )//' on '//CurDate()//' at '//CurTime()//'.'

FormStr = "( / 'This summary file was ', A / )"
WRITE (US,FormStr)  TRIM(DescStr)

   ! Capitalize the first letter of the string.

DescStr = 'This full-field file was '//TRIM(DescStr)


RETURN
END SUBROUTINE GetFiles
!=======================================================================
SUBROUTINE GetInput


  ! This subroutine is used to read parameters from the input file.

USE                  TSMods

IMPLICIT             NONE

   ! Local variables

REAL(ReKi)        :: InCVar     (2)                           ! Contains the coherence parameters (used for input)
REAL(ReKi)        :: RefHt           ! Height for reference wind speed.
REAL(ReKi)        :: tmp             ! variable for estimating Ustar
REAL(ReKi)        :: TmpUary (3)     !Temporary vector to store windSpeed(z) values
REAL(ReKi)        :: TmpUstar(3)     !Temporary vector to store ustar(z) values
REAL(ReKi)        :: TmpUstarD       !Temporary ustarD value
REAL(ReKi)        :: TmpZary (3)     !Temporary vector to store height(z) values
REAL(ReKi)        :: TmpZLary(3)     !Temporary vector to store zL(z) values
REAL(ReKi)        :: URef            ! Wind speed at the reference height.

INTEGER           :: IOS             ! I/O status
INTEGER           :: TmpIndex        ! Contains the index number when searching for substrings

LOGICAL           :: getPLExp        ! Whether the power law exponent needs to be calculated
LOGICAL           :: Randomize       ! Whether to randomize the coherent turbulence scaling
LOGICAL           :: UseDefault      ! Whether or not to use a default value


CHARACTER(99)     :: Line            ! An input line
CHARACTER(1)      :: Line1           ! The first character of an input line


UnEc = US
Echo = .FALSE.       ! Do not echo the input into a separate file

   !==========================================================
   ! Initialize temporary variables   
   !==========================================================
TmpUary    = (/ 0.0, 0.0, 0.0 /)
TmpUstarD  = 0.0
TmpUstar   = (/ 0.0, 0.0, 0.0 /)
UstarOffset= 0.0
UstarSlope = 1.0
zlOffset   = 0.0

   ! Open input file.
CALL OpenFInpFile( UI, InFile )

CALL WrScr1(' Reading the input file "'//TRIM(InFile)//'".' )

   !==========================================================
   ! Read the runtime options.
   !==========================================================

FormStr = "( / 'Runtime Options:' / )"
WRITE (US,FormStr)
CALL ReadCom( UI, InFile, "File Heading Line 1" )
CALL ReadCom( UI, InFile, "File Heading Line 2" )
CALL ReadCom( UI, InFile, "Runtime Options Heading" )
!READ (UI,'(//)')


   ! ---------- Read Random seed 1 -----------------------

CALL ReadIVar( UI, InFile, RandSeed(1), "RandSeed(1)", "Random seed #1")

FormStr = "( I10 , 2X , 'Random seed #1' )"
WRITE (US,FormStr)  RandSeed(1)

RandSeedTmp = RandSeed(1)

   ! ---------- Read Random seed 2 -----------------------

CALL ReadCVar( UI, InFile, Line, "RandSeed(2)", "Random seed #2")

   ! Check if alternate random number generator is to be used

   READ (Line,*,IOSTAT=IOS) Line1

   CALL Conv2UC( Line1 )

   IF ( (Line1 == 'T') .OR.  (Line1 == 'F') ) THEN
      CALL TS_Abort (' Invalid RNG type. ')
   ENDIF

   READ (Line,*,IOSTAT=IOS)  RandSeed(2)

   IF (IOS == 0) THEN

      FormStr = "( I10 , 2X , 'Random seed #2' )"
      WRITE (US,FormStr)  RandSeed(2)
      RNG_type = "NORMAL"

   ELSE

      RNG_type = ADJUSTL( Line )
      CALL Conv2UC( RNG_type )

      IF ( RNG_type /= "RANLUX" .AND. RNG_type /= "RNSNLW" ) THEN
         CALL TS_Abort( 'Invalid alternative random number generator.' )
      ENDIF

      FormStr  = "( 4X, A6, 2X, 'Type of random number generator' )"
      WRITE (US,FormStr)  RNG_type

   ENDIF

   ! Initialize the RNG

CALL RndInit()


   ! --------- Read the flag for writing the binary HH (GenPro) turbulence parameters. -------------

CALL ReadLVar( UI, InFile, WrBHHTP, "the flag for writing the binary HH turbulence parameters", &
                                    "Output binary HH turbulence parameters?")

FormStr = "( L10 , 2X , 'Output binary HH turbulence parameters?' )"
WRITE (US,FormStr)  WrBHHTP


   ! --------- Read the flag for writing the formatted turbulence parameters. ----------------------

CALL ReadLVar( UI, InFile, WrFHHTP, "the flag for writing the formatted turbulence parameters", &
                                    "Output formatted turbulence parameters?")

FormStr = "( L10 , 2X , 'Output formatted turbulence parameters?' )"
WRITE (US,FormStr)  WrFHHTP


   ! ---------- Read the flag for writing the AeroDyn HH files. -------------------------------------

CALL ReadLVar( UI, InFile, WrADHH, "the flag for writing AeroDyn HH files", "Output AeroDyn HH files?")

FormStr = "( L10 , 2X , 'Output AeroDyn HH files?' )"
WRITE (US,FormStr)  WrADHH

   ! ---------- Read the flag for writing the AeroDyn FF files. ---------------------------------------

CALL ReadLVar( UI, InFile, WrADFF, "the flag for writing AeroDyn FF files", "Output AeroDyn FF files?")

FormStr = "( L10 , 2X , 'Output AeroDyn FF files?' )"
WRITE (US,FormStr)  WrADFF

   ! ---------- Read the flag for writing the BLADED FF files. -----------------------------------------

CALL ReadLVar( UI, InFile, WrBLFF, "the flag for writing BLADED FF files", "Output BLADED FF files?")

FormStr = "( L10 , 2X , 'Output BLADED FF files?' )"
WRITE (US,FormStr)  WrBLFF


   ! ---------- Read the flag for writing the AeroDyn tower files. --------------------------------------

CALL ReadLVar( UI, InFile, WrADTWR, "the flag for writing tower data", "Output tower data?")

FormStr = "( L10 , 2X , 'Output tower data?' )"
WRITE (US,FormStr)  WrADTWR


   ! ---------- Read the flag for writing the formatted FF files. ---------------------------------------

CALL ReadLVar( UI, InFile, WrFmtFF, "the flag for writing formatted FF files", "Output formatted FF files?")

FormStr = "( L10 , 2X , 'Output formatted FF files?' )"
WRITE (US,FormStr)  WrFmtFF


   ! ---------- Read the flag for writing coherent time series files. --------------------------------------

CALL ReadLVar( UI, InFile, WrACT, "the flag for writing coherent time series files", "Output coherent time series files?")

FormStr = "( L10 , 2X , 'Output coherent turbulence time step file?' )"
WRITE (US,FormStr)  WrACT


   ! ---------- Read the flag for turbine rotation. -----------------------------------------------------------

CALL ReadLVar( UI, InFile, Clockwise, "the flag for direction of turbine rotation", "Clockwise rotation when looking downwind?")

FormStr = "( L10 , 2X , 'Clockwise rotation when looking downwind?' )"
WRITE (US,FormStr)  Clockwise

   ! ---------- Read the flag for determining IEC scaling -----------------------------------------------------
CALL ReadLIVar( UI, InFile, ScaleIEC, "ScaleIEC, the switch for scaling IEC turbulence", &
               "Scale IEC turbulence models to specified standard deviation?")

FormStr = "( I2, ' - ', A5, 2X , 'IEC turbulence models scaled to exact specified standard deviation' )"
   SELECT CASE ( ScaleIEC )
      CASE (0)
         WRITE (US,FormStr)  ScaleIEC, "NONE"
      CASE (1, -1)   ! included the -1 for reading t/f on other systems
         WRITE (US,FormStr)  ScaleIEC, "HUB"
         ScaleIEC = 1;
      CASE (2)
         WRITE (US,FormStr)  ScaleIEC, "ALL"
      CASE DEFAULT
         CALL TS_Abort ( 'The value for parameter ScaleIEC must be 0, 1, or 2.' )   
   ENDSELECT
         

   !==================================================================================
   ! Read the turbine/model specifications.
   !===================================================================================

FormStr = "( // 'Turbine/Model Specifications:' / )"
WRITE (US,FormStr)
CALL ReadCom( UI, InFile, "Turbine/Model Specifications Heading Line 1" )
CALL ReadCom( UI, InFile, "Turbine/Model Specifications Heading Line 2" )
!READ (UI,'(/)')


   ! ------------ Read in the vertical matrix dimension. ---------------------------------------------

CALL ReadIVar( UI, InFile, NumGrid_Z, "the vertical matrix dimension", "Vertical grid-point matrix dimension")

   IF ( NumGrid_Z < 2 )  THEN
      CALL TS_Abort ( 'The matrix must be >= 2x2.' )
   ENDIF

FormStr = "( I10 , 2X , 'Vertical grid-point matrix dimension' )"
WRITE (US,FormStr)  NumGrid_Z


   ! ------------ Read in the lateral matrix dimension. ---------------------------------------------

CALL ReadIVar( UI, InFile, NumGrid_Y, "the horizontal matrix dimension", "Horizontal grid-point matrix dimension")

   IF ( NumGrid_Y < 2 )  THEN
      CALL TS_Abort ( 'The matrix must be >= 2x2.' )
   ENDIF

FormStr = "( I10 , 2X , 'Horizontal grid-point matrix dimension' )"
WRITE (US,FormStr)  NumGrid_Y


   ! ------------ Read in the time step. ---------------------------------------------

CALL ReadRVar( UI, InFile, TimeStep, "the time step", "Time step [seconds]")

   IF ( TimeStep <=  0.0 )  THEN
      CALL TS_Abort ( 'The time step must be greater than zero.' )
   ENDIF

FormStr = "( F10.3 , 2X , 'Time step [seconds]' )"
WRITE (US,FormStr)  TimeStep


   ! ------------ Read in the analysis time. ---------------------------------------------

CALL ReadRVar( UI, InFile, AnalysisTime, "the analysis time", "Analysis time [seconds]")

   IF ( AnalysisTime <=  0.0 )  THEN
      CALL TS_Abort ( 'The analysis time must be greater than zero.' )
   ENDIF

FormStr = "( F10.3 , 2X , 'Analysis time [seconds]' )"
WRITE (US,FormStr)  AnalysisTime


   ! ------------ Read in the usable time. ---------------------------------------------

CALL ReadRVar( UI, InFile, UsableTime, "the usable output time", "Usable output time [seconds]")

   IF ( UsableTime <=  0.0 )  THEN
      CALL TS_Abort ( 'The usable output time must be greater than zero.' )
   ENDIF

FormStr = "( F10.3 , 2X , 'Usable output time [seconds]' )"
WRITE (US,FormStr)  UsableTime


   ! ------------ Read in the hub height. ---------------------------------------------

CALL ReadRVar( UI, InFile, HubHt, "the hub height", "Hub height [m]")

   IF ( HubHt <=  0.0 )  THEN
      CALL TS_Abort ( 'The hub height must be greater than zero.' )
   ENDIF

FormStr = "( F10.3 , 2X , 'Hub height [m]' )"
WRITE (US,FormStr)  HubHt


   ! ------------ Read in the grid height. ---------------------------------------------

CALL ReadRVar( UI, InFile, GridHeight, "the grid height", "Grid height [m]")

   IF ( 0.5*GridHeight > HubHt  )THEN
      CALL TS_Abort( 'The hub must be higher than half of the grid height.')
   ENDIF

FormStr = "( F10.3 , 2X , 'Grid height [m]' )"
WRITE (US,FormStr)  GridHeight


   ! ------------ Read in the grid width. ---------------------------------------------

CALL ReadRVar( UI, InFile, GridWidth, "the grid width", "Grid width [m]")

   IF ( GridWidth <=  0.0 )  THEN
      CALL TS_Abort ( 'The grid width must be greater than zero.' )
   ENDIF

FormStr = "( F10.3 , 2X , 'Grid width [m]' )"
WRITE (US,FormStr)  GridWidth


   ! ***** Calculate the diameter of the rotor disk *****

RotorDiameter = MIN( GridWidth, GridHeight )


   ! ------------ Read in the vertical flow angle. ---------------------------------------------

CALL ReadRVar( UI, InFile, VFlowAng, "the vertical flow angle", "Vertical flow angle [degrees]")

   IF ( ABS( VFlowAng ) > 45.0 )  THEN
      CALL TS_Abort ( 'The vertical flow angle must not exceed +/- 45 degrees.' )
   ENDIF

FormStr = "( F10.3 , 2X , 'Vertical flow angle [degrees]' )"
WRITE (US,FormStr)  VFlowAng


   ! ------------ Read in the horizontal flow angle. ---------------------------------------------

CALL ReadRVar( UI, InFile, HFlowAng, "the horizontal flow angle", "Horizontal flow angle [degrees]")

FormStr = "( F10.3 , 2X , 'Horizontal flow angle [degrees]' )"
WRITE (US,FormStr)  HFlowAng


   !==========================================================
   ! Read the meteorological boundary conditions.
   !==========================================================

FormStr = "( // 'Meteorological Boundary Conditions:' / )"
WRITE (US,FormStr)
!READ (UI,'(/)')
CALL ReadCom( UI, InFile, "Meteorological Boundary Conditions Heading Line 1" )
CALL ReadCom( UI, InFile, "Meteorological Boundary Conditions Heading Line 2" )


   ! ------------ Read in the turbulence model. ---------------------------------------------

CALL ReadCVar( UI, InFile, TurbModel, "the spectral model", "spectral model")

   TurbModel = ADJUSTL( TurbModel )
   CALL Conv2UC( TurbModel )

   SELECT CASE ( TRIM(TurbModel) )
      CASE ( 'IECKAI' )
         TMName = 'IEC Kaimal'
      CASE ( 'IECVKM' )
         TMName = 'IEC von Karman'
      CASE ( 'TIDAL' )
         TMName = 'Tidal Channel Turbulence'
      CASE ( 'RIVER' )
         TMName = 'River Turbulence'
      CASE ( 'IECMAN' )
         TMName = 'IEC Mann'
      CASE ( 'SMOOTH' )
         TMName = 'RISO Smooth Terrain'
      CASE ( 'WF_UPW' )
         TMName = 'NREL Wind Farm Upwind'
      CASE ( 'WF_07D' )
         TMName = 'NREL 7D Spacing Wind Farm'
      CASE ( 'WF_14D' )
         TMName = 'NREL 14D Spacing Wind Farm'
      CASE ( 'NONE'   )
         TMName = 'No fluctuating wind components'
      CASE ( 'MODVKM' )
         TMName = 'Modified von Karman'
      CASE ( 'NWTCUP' )
         TMName = 'NREL National Wind Technology Center'
      CASE ( 'GP_LLJ' )
         TMName = 'Great Plains Low-Level Jet'
      CASE ( 'USRVKM' )
         TMName = 'von Karman model with user-defined specifications'
      CASE ( 'USRINP' )
         TMName = 'User-input '
         CALL GetUSRspec("UsrSpec.inp")      ! bjj: before documenting, please fix this hard-coded name!
      CASE DEFAULT
!BONNIE: add the UsrVKM model to this list when the model is complete as well as USRINP
         CALL TS_Abort ( 'The turbulence model must be "IECKAI", "IECVKM", "SMOOTH",' &
                    //' "WF_UPW", "WF_07D", "WF_14D", "NWTCUP", "GP_LLJ", "TIDAL", or "NONE".' )

   END SELECT  ! TurbModel

FormStr = "( 4X , A6 , 2X , '"//TRIM( TMName )//" spectral model' )"
WRITE (US,FormStr)  TurbModel

   ! ------------ Read in the IEC standard and edition numbers. ---------------------------------------------

CALL ReadCVar( UI, InFile, Line, "the number of the IEC standard", "Number of the IEC standard")

   IF ( TurbModel(1:3) == 'IEC' ) THEN

      CALL Conv2UC( LINE )

      IF ( (Line(1:1) == 'T') .OR.  (Line(1:1) == 'F') ) THEN
         CALL TS_Abort ( 'The number of the IEC standard must be either "1", "2", or "3"' &
                          // ' with an optional IEC 61400-1 edition number ("1-ED2").' )
      ENDIF

      TmpIndex = INDEX(Line, "-ED")

      IF ( TmpIndex > 0 ) THEN
         READ ( Line(TmpIndex+3:),*,IOSTAT=IOS ) IECedition

         CALL CheckIOS( IOS, InFile, 'the IEC edition number', NumType )

         IF ( IECedition < 1 .OR. IECedition > SIZE(IECeditionSTR) ) THEN
            CALL TS_Abort( 'Invalid IEC edition number.' )
         ENDIF

         Line = Line(1:TmpIndex-1)
      ELSE
         IECedition = 0
      ENDIF

      READ ( Line,*,IOSTAT=IOS ) IECstandard

      SELECT CASE ( IECstandard )
         CASE ( 1 )
            IF (IECedition < 1 ) THEN  ! Set up the default
               IF ( TurbModel(4:6) == 'VKM' ) THEN
                  IECedition = 2   ! The von Karman model is not specified in edition 3 of the -1 standard
               ELSE
                  IECedition = 3
               ENDIF
            ELSE
               IF ( IECedition < 2 ) THEN
                  CALL TS_Abort( 'The IEC edition number must be 2 or 3' )
               ENDIF
            ENDIF

         CASE ( 2 )
               ! The scaling is the same as 64100-1, Ed. 2 with "A" or user-specified turbulence
            IF (IECedition < 1 ) THEN  ! Set up the default
               IECedition = 2
            ELSE
               CALL TS_Abort( ' The edition number cannot be specified for the 61400-2 standard.')
            ENDIF
            IECeditionSTR(IECedition) = 'IEC 61400-2 Ed. 2: 2005'

         CASE ( 3 )
               ! The scaling for 61400-3 (Offshore) is the same as 61400-1 except it has a different power law exponent
            IF (IECedition < 1 ) THEN  ! Set up the default

               IF ( TurbModel /= 'IECKAI' ) THEN
                  CALL TS_Abort( ' The von Karman model is not valid for the 61400-3 standard.')
               ENDIF
               IECedition = 3   ! This is the edition of the -1 standard

            ELSE
               CALL TS_Abort( ' The edition number cannot be specified for the 61400-3 standard.')
            ENDIF
            IECeditionSTR(IECedition) = 'IEC 61400-3 Ed. 1: 2006'

         CASE DEFAULT
            CALL TS_Abort( 'The number of the IEC 61400-x standard must be 1, 2, or 3.')

      END SELECT

      FormStr = "( 7X, I3, 2X, 'IEC standard: ', A )"
      WRITE (US,FormStr)  IECstandard, TRIM(IECeditionSTR(IECedition))

   ELSE ! NOT IEC
      IECstandard = 0
      IECedition  = 0

      FormStr = "( 7X, A3, 2X, 'IEC standard' )"
      WRITE (US,FormStr)  'N/A'

   ENDIF ! IEC


   ! ------------ Read in the IEC turbulence characteristic. ---------------------------------------------

CALL ReadCVar( UI, InFile, Line, "the IEC turbulence characteristic", "IEC turbulence characteristic")
!!$! -- begin block --
!!$! This block reads turbulence intensity, and stores it in the variable TurbIntH20.  This variable is not currently used, but will be soon for user-specified turbulence intensity for the HYDRO spectral models.
!!$! This code is copied from TurbSim.f90
!!$READ (Line,*,IOSTAT=IOS)  PerTurbInt ! This is to read the 
!!$IECTurbC = ADJUSTL( Line )
!!$IF ( IOS /= 0 ) THEN
!!$   CALL Conv2UC(IECTurbC)
!!$   IF ( IECTurbC == 'A' ) THEN
!!$      TurbIntH20  = 0.16
!!$   ELSEIF ( IECTurbC == 'B' ) THEN
!!$      TurbIntH20  = 0.14
!!$   ELSEIF ( IECTurbC == 'C' ) THEN
!!$      TurbIntH20  = 0.12
!!$   ELSEIF ( IECTurbC == 'K' ) THEN
!!$      TurbIntH20  = 0.16
!!$   ELSE   ! We should never get here, but just to be complete...
!!$      !print *, IECTurbC
!!$      CALL TS_Abort( ' Invalid IEC turbulence characteristic.' )
!!$   ENDIF
!!$ELSE
!!$   TurbIntH20 = 0.01*PerTurbInt
!!$ENDIF
!!$! -- end block --
   


   IF ( ( TurbModel(1:3) == 'IEC' ) .OR. ( TurbModel == 'MODVKM' ) ) THEN

      READ (Line,*,IOSTAT=IOS) Line1

      CALL Conv2UC( Line1 )

      IF ( (Line1 == 'T') .OR.  (Line1 == 'F') ) THEN
         CALL TS_Abort ( 'The IEC turbulence characteristic must be either "A", "B", "C", or a real number.')
      ENDIF

      ! Check to see if the entry was a number.

      READ (Line,*,IOSTAT=IOS)  PerTurbInt

      IF ( IOS == 0 )  THEN

         ! Let's use turbulence value.

         NumTurbInp = .TRUE.
         FormStr    = "( F10.3 , 2X , 'Percent turbulence intensity, ', A )"
         WRITE (US,FormStr)  PerTurbInt, TRIM(IECeditionSTR(IECedition))

      ELSE

         ! Let's use one of the standard turbulence values (A or B or C).

         NumTurbInp = .FALSE.
         IECTurbC   = ADJUSTL( Line )
         CALL Conv2UC(  IECTurbC )

         SELECT CASE ( IECTurbC )
            CASE ( 'A' )
            CASE ( 'B' )
               IF ( IECstandard == 2 ) THEN
                  CALL TS_Abort (' The IEC 61400-2 turbulence characteristic must be either "A" or a real number.' )
               ENDIF
            CASE ( 'C' )
               IF ( IECstandard == 2 ) THEN
                  CALL TS_Abort (' The IEC 61400-2 turbulence characteristic must be either "A" or a real number.' )
               ELSEIF ( IECedition < 3 ) THEN
                  CALL TS_Abort (' The turbulence characteristic for '//TRIM(IECeditionSTR(IECedition) )// &
                                  ' must be either "A", "B", or a real number.' )
               ENDIF
            CASE DEFAULT
               CALL TS_Abort ( 'The IEC turbulence characteristic must be either "A", "B", "C", or a real number.' )
         END SELECT  ! IECTurbC

         FormStr  = "( 9X , A1 , 2X , 'IEC turbulence characteristic' )"
         WRITE (US,FormStr)  IECTurbC

      ENDIF

   ELSE  !not IECKAI and not IECVKM and not MODVKM

      Line = ADJUSTL( Line )
      CALL Conv2UC( Line )

      IF ( Line(1:6) == 'KHTEST' ) THEN
         KHtest = .TRUE.
         FormStr = "( 4X, A6, 2X, 'Kelvin-Helmholtz billow test case' )"
         WRITE (US,FormStr)  'KHTEST'


         IF ( .NOT. TurbModel == 'NWTCUP' ) THEN
            CALL TS_Abort( 'The KH test can be used with the "NWTCUP" spectral model only.' )
         ENDIF

         IF ( .NOT. WrACT ) THEN
            CALL WRScr( ' Coherent turbulence time step files must be generated when using the "KHTEST" option.' )
            WRACT  = .TRUE.
         ENDIF

      ELSE
         FormStr = "( 7X, A3, 2X, 'IEC turbulence characteristic' )"
         WRITE (US,FormStr)  'N/A'
      ENDIF


         ! These variables are not used for non-IEC turbulence

      IECedition = 0
      NumTurbInp = .TRUE.
      PerTurbInt = 0.0

   ENDIF

   ! ------------ Read in the IEC wind turbulence type ---------------------------------------------

CALL ReadCVar( UI, InFile, Line, "the IEC turbulence type", "IEC turbulence type")

   IF ( ( TurbModel(1:3) == 'IEC' ) ) THEN

      CALL Conv2UC( Line )

      IECTurbE   = Line(1:1)

      ! Let's see if the first character is a number (for the ETM case)
      SELECT CASE ( IECTurbE )
         CASE ('1')
            Vref = 50.0
            Line = Line(2:)
         CASE ('2')
            Vref = 42.5
            Line = Line(2:)
         CASE ('3')
            Vref = 37.5
            Line = Line(2:)
         CASE DEFAULT
               ! There's no number at the start of the string so let's move on (it's NTM).
            Vref = -999.9
            IECTurbE = ' '
      END SELECT

      SELECT CASE ( TRIM( Line ) )
         CASE ( 'NTM'   )
            IEC_WindType = IEC_NTM
            IEC_WindDesc = 'Normal Turbulence Model'
         CASE ( 'ETM'   )
            IEC_WindType = IEC_ETM
            IEC_WindDesc = 'Extreme Turbulence Model'
         CASE ( 'EWM1'  )
            IEC_WindType = IEC_EWM1
            IEC_WindDesc = 'Extreme 1-Year Wind Speed Model'
         CASE ( 'EWM50' )
            IEC_WindType = IEC_EWM50
            IEC_WindDesc = 'Extreme 50-Year Wind Speed Model'
         CASE DEFAULT
            CALL TS_Abort ( ' Valid entries for the IEC wind turbulence are "NTM", "xETM", "xEWM1", or "xEWM50", '// &
                             'where x is the wind turbine class (1, 2, or 3).' )
      END SELECT

      IF ( IEC_WindType /= IEC_NTM ) THEN

         IF (IECedition /= 3 .OR. IECstandard == 2) THEN
            CALL TS_Abort ( ' The extreme turbulence and extreme wind speed models are available with '// &
                         'the IEC 61400-1 Ed. 3 or 61400-3 scaling only.')
         ENDIF

         IF (Vref < 0. ) THEN
            CALL TS_Abort ( ' A wind turbine class (1, 2, or 3) must be specified with the '// &
                         'extreme turbulence and extreme wind types. (i.e. "1ETM")')
         ENDIF

         IF ( NumTurbInp ) THEN
            CALL TS_Abort ( ' When the turbulence intensity is entered as a percent, the IEC wind type must be "NTM".' )
         ENDIF
         
      ELSE

         IECTurbE = ' '

      ENDIF

      FormStr = "( 4X, A6 , 2X , 'IEC ', A )"
      WRITE (US,FormStr)  TRIM(IECTurbE)//TRIM(Line), TRIM(IEC_WindDesc)

   ELSE
      IEC_WindType = IEC_NTM

      FormStr = "( A10 , 2X , 'IEC turbulence type' )"
      WRITE (US,FormStr)  'N/A'
   ENDIF

   ! ------------ Read in the ETM c parameter (IEC 61400-1, Ed 3: Section 6.3.2.3, Eq. 19) ----------------------
UseDefault = .TRUE.
ETMc = 2;
CALL ReadRVarDefault( UI, InFile, ETMc, "the ETM c parameter", 'IEC Extreme Turbulence Model (ETM) "c" parameter [m/s]', &
                      UseDefault, IGNORE=(IEC_WindType /= IEC_ETM ))

   IF ( ETMc <= 0. ) THEN
      CALL TS_Abort('The ETM "c" parameter must be a positive number');
   ENDIF

   ! ------------ Read in the wind profile type -----------------------------------------------------------------

UseDefault = .TRUE.         ! Calculate the default value
SELECT CASE ( TRIM(TurbModel) )
   CASE ( 'GP_LLJ' )
      WindProfileType = 'JET'
   CASE ( 'IECKAI','IECVKM','IECMAN','MODVKM' )
      WindProfileType = 'IEC'
   CASE ( 'TIDAL' )
      WindProfileType = 'H2L'
   CASE ( 'USRVKM' )
      WindProfileType = 'USR'
   CASE DEFAULT
      WindProfileType = 'IEC'
END SELECT

CALL ReadCVarDefault( UI, InFile, WindProfileType, "the wind profile type", "Wind profile type", UseDefault )

      ! Make sure the variable is valid for this turbulence model

   SELECT CASE ( TRIM(WindProfileType) )
      CASE ( 'JET','J' )
         IF ( TurbModel /= 'GP_LLJ' ) THEN
            CALL TS_Abort( 'The jet wind profile is available with the GP_LLJ spectral model only.')
         ENDIF
      CASE ( 'LOG','L' )
         IF (IEC_WindType /= IEC_NTM ) THEN
            CALL TS_Abort( ' The IEC turbulence type must be NTM for the logarithmic wind profile.' )
!bjj check that IEC_WindType == IEC_NTM for non-IEC
         ENDIF
      CASE ( 'PL',  'P' )
      CASE ( 'H2L', 'H' )
         IF ( TRIM(TurbModel)/='TIDAL' ) THEN
            CALL TS_Abort(  'The "H2L" mean profile type should be used only with the "TIDAL" spectral model.' )
         ENDIF
      CASE ( 'IEC', 'N/A' )
      CASE ( 'USR', 'U' )
      CASE DEFAULT
         CALL TS_Abort( 'The wind profile type must be "JET", "LOG", "PL", "IEC", "USR", "H2L", or default.' )
   END SELECT

   IF ( TRIM(TurbModel)=='TIDAL' .AND. TRIM(WindProfileType) /= "H2L" ) THEN
      WindProfileType = 'H2L'
      CALL TS_Warn  ( 'Overwriting wind profile type to "H2L" for the "TIDAL" spectral model.', .TRUE.)
   ENDIF

   IF ( KHTest ) THEN
      IF ( TRIM(WindProfileType) /= 'IEC' .AND. TRIM(WindProfileType) /= 'PL' ) THEN
         WindProfileType = 'IEC'
         CALL TS_Warn  ( 'Overwriting wind profile type for the KH test.', .FALSE.)
      ENDIF
   ENDIF


   ! ------------ Read in the height for the reference wind speed. ---------------------------------------------

CALL ReadRVar( UI, InFile, RefHt, "the reference height", "Reference height [m]")

   IF ( RefHt <=  0.0 .AND. WindProfileType(1:1) /= 'U' )  THEN
      CALL TS_Abort ( 'The reference height must be greater than zero.' )
   ENDIF

FormStr = "( F10.3 , 2X , 'Reference height [m]' )"
WRITE (US,FormStr)  RefHt
H_ref = RefHt ! Define the variable H_ref, for later use in HYDRO spectral models (RefHt gets modified later in the code)

   ! ------------ Read in the reference wind speed. -----------------------------------------------------

UseDefault = .FALSE.
URef       = -999.9

! If we specify a Ustar (i.e. if Ustar /= "default") then we can enter "default" here,
! otherwise, we get circular logic...

CALL ReadRVarDefault( UI, InFile, URef, "the reference wind speed", "Reference wind speed [m/s]", UseDefault, &
                     IGNORE=(IEC_WindType == IEC_EWM1 .OR. IEC_WindType == IEC_EWM50 .OR. WindProfileType(1:1) == 'U') )

   NumUSRz = 0  ! initialize the number of points in a user-defined wind profile

   IF ( ( WindProfileType(1:1) /= 'J' .OR. .NOT. UseDefault) .AND. &
        (IEC_WindType /= IEC_EWM1 .AND. IEC_WindType /= IEC_EWM50 .AND. WindProfileType(1:1) /= 'U') ) THEN
      IF ( URef <=  0.0 )  THEN
         CALL TS_Abort ( 'The reference wind speed must be greater than zero.' )
      ENDIF

   ELSEIF ( WindProfileType(1:1) == 'U' ) THEN
      RefHt = HubHt
      CALL GetUSR( UI, InFile, 37 ) !Read the last several lines of the file, then return to line 37
      URef = getWindSpeed(URef, RefHt, RefHt, RotorDiameter, PROFILE=WindProfileType) !This is UHub

   ENDIF   ! Otherwise, we're using a Jet profile with default wind speed (for now it's -999.9)


   ! ------------ Read in the jet height -------------------------------------------------------------

UseDefault = .FALSE.
ZJetMax    = -999.9

CALL ReadRVarDefault( UI, InFile, ZJetMax, "the jet height", "Jet height [m]", UseDefault, IGNORE=WindProfileType(1:1) /= 'J')

   IF ( WindProfileType(1:1) == 'J' .AND. .NOT. UseDefault ) THEN
      IF ( ZJetMax <  70.0 .OR. ZJetMax > 490.0 )  THEN
         CALL TS_Abort ( 'The height of the maximum jet wind speed must be between 70 and 490 m.' )
      ENDIF
   ENDIF


   ! ------------ Read in the power law exponent, PLExp ---------------------------------------------

SELECT CASE ( TurbModel )
   CASE ('WF_UPW','WF_07D','WF_14D','NWTCUP')
      IF ( KHtest ) THEN
         UseDefault = .TRUE.
         PLExp      = 0.3
      ELSE
         UseDefault = .FALSE.             ! This case needs the Richardson number to get a default
         PLExp      = 0.
      ENDIF

   CASE DEFAULT
      UseDefault = .TRUE.
      PLExp      = PowerLawExp( RICH_NO )  ! These cases do not use the Richardson number to get a default

END SELECT
getPLExp = .NOT. UseDefault

CALL ReadRVarDefault( UI, InFile, PLExp, "the power law exponent", "Power law exponent", UseDefault, &
               IGNORE=( ((INDEX( 'JLUH', WindProfileType(1:1) ) > 0.)) .OR. &
                         IEC_WindType == IEC_EWM1 .OR. IEC_WindType == IEC_EWM50 ) )

   IF ( .NOT. UseDefault ) THEN  ! We didn't use a default (we entered a number on the line)
      getPLExp = .FALSE.

      IF ( KHtest ) THEN
         IF ( PLExp /= 0.3 ) THEN
            PLExp = 0.3
            CALL TS_Warn  ( 'Overwriting the power law exponent for KH test.', .FALSE. )
         ENDIF
      ENDIF
   ENDIF


   ! ------------ Read in the surface roughness length, Z0 (that's z-zero) ---------------------------------------------

UseDefault = .TRUE.
SELECT CASE ( TurbModel )
   CASE ('SMOOTH')
      Z0 = 0.010
   CASE ('GP_LLJ')
      Z0 = 0.005
   CASE ('WF_UPW')
      Z0 = 0.018
   CASE ('NWTCUP')
      Z0 = 0.021
   CASE ('WF_07D')
      Z0 = 0.233
   CASE ('WF_14D')
      Z0 = 0.064
   CASE DEFAULT !IEC values
      Z0 = 0.030 ! Represents smooth, homogenous terrain
END SELECT

CALL ReadRVarDefault( UI, InFile, Z0, "the roughness length", "Surface roughness length [m]", UseDefault, &
                       IGNORE=TRIM(TurbModel)=='TIDAL')

   IF ( Z0 <= 0.0 ) THEN
      CALL TS_Abort ( 'The surface roughness length must be a positive number or "default".')
   ENDIF


   !=================================================================================
   ! Read the meteorological boundary conditions for non-IEC models. !
   !=================================================================================

IF ( TurbModel /= 'IECKAI' .AND. TurbModel /= 'IECVKM' ) THEN

   FormStr = "( // 'Non-IEC Meteorological Boundary Conditions:' / )"
   WRITE (US,FormStr)
   !READ (UI,'(/)')
   CALL ReadCom( UI, InFile, "Non-IEC Meteorological Boundary Conditions Heading Line 1" )
   CALL ReadCom( UI, InFile, "Non-IEC Meteorological Boundary Conditions Heading Line 2" )
   


      ! ------------ Read in the site latitude, LATITUDE. ---------------------------------------------

   UseDefault = .TRUE.
   Latitude   = 45.0

   CALL ReadRVarDefault( UI, InFile, Latitude, "the site latitude", "Site latitude [degrees]", UseDefault)

      IF ( ABS(Latitude) < 5.0 .OR. ABS(Latitude) > 90.0 ) THEN
         CALL TS_Abort( 'The latitude must be between -90 and 90 degrees but not between -5 and 5 degrees.' )
      ENDIF

   Fc = 2.0 * Omega * SIN( ABS(Latitude*D2R) )  ! Calculate Coriolis parameter from latitude

ELSE

   Latitude = 0.0                               !Not used in IEC specs
   Fc = 0.0

ENDIF    ! Not IECKAI and Not IECVKM

IF ( TurbModel /= 'IECKAI' .AND. TurbModel /= 'IECVKM' .AND. TurbModel /= 'MODVKM' .AND. TurbModel /= 'IECMAN') THEN


      ! ------------ Read in the gradient Richardson number, RICH_NO. ---------------------------------------------

   CALL ReadRVar( UI, InFile, RICH_NO, "the gradient Richardson number", "Gradient Richardson number")

   IF ( KHtest ) THEN
      IF ( RICH_NO /= 0.02 ) THEN
         RICH_NO = 0.02
         CALL TS_Warn ( 'Overwriting the Richardson Number for KH test.', .FALSE. )
      ENDIF
   ENDIF

   IF ( TurbModel(1:1) == 'U' ) THEN
      IF ( RICH_NO /= 0.0 ) THEN
         RICH_NO = 0.0
         CALL TS_Warn ( 'Overwriting the Richardson Number for the '//TRIM(TurbModel)//' model.', .FALSE. )
      ENDIF
   ENDIF

   IF ( TRIM(TurbModel)=='TIDAL' ) THEN
      WRITE (US,"( A10, 2X, A )")  "N/A", 'Gradient Richardson number'
      RICH_NO = 0
   ELSE      
      FormStr = "( F10.3 , 2X , 'Gradient Richardson number' )"
      WRITE (US,FormStr)  RICH_NO
   END IF

      ! ***** Calculate M-O z/L parameter   :  z/L is a number in (-inf, 1] *****

   IF ( TurbModel == 'NWTCUP' ) THEN
         ! Calculate disk averaged Z/L from turbine layer Ri for NWTC/LIST experiment

      RICH_NO = MIN( MAX( RICH_NO, REAL(-1.0,ReKi) ), REAL(1.0,ReKi) )  ! Ensure that: -1 <= RICH_NO <= 1

      IF ( RICH_NO <= -0.1 ) THEN
         ZL = -0.254 + 1.047*RICH_NO
      ELSEIF ( RICH_NO < 0 ) THEN
         ZL = 10.369*RICH_NO/(1.0 - 19.393*RICH_NO)
      ELSE  !( RICH_NO < 0.155 ) THEN
         ZL = 2.535*MIN( RICH_NO, REAL(0.155, ReKi) ) / (1.0 - 6.252*MIN( RICH_NO, REAL(0.155,ReKi) ))
      ENDIF


   ELSEIF (TurbModel == 'GP_LLJ') THEN

      RICH_NO = MIN( MAX( RICH_NO, REAL(-1.0,ReKi) ), REAL(1.0,ReKi) )  ! Ensure that: -1 <= RICH_NO <= 1

      IF ( RICH_NO <= -0.1 ) THEN
         ZL = -0.047 + 1.054*RICH_NO
      ELSEIF ( RICH_NO < 0 ) THEN
         ZL = 2.213*RICH_NO/(1.0 - 4.698*RICH_NO)
      ELSE  !( RICH_NO < 0.1367 ) THEN
         ZL = 3.132*MIN( RICH_NO, REAL(0.1367,ReKi) ) / (1.0 - 6.762*MIN( RICH_NO, REAL(0.1367,ReKi) ))
      ENDIF

   ELSE ! see Businger, J.A.; Wyngaard, J.C.; Izumi, Y.; Bradley, E.F. (1971). "Flux-Profile Relationships in the Atmospheric Surface Layer." Journal of the Atmospheric Sciences (28); pp.181-189.

      IF ( RICH_NO <= 0.0 ) THEN
         ZL = RICH_NO
         !PhiM = (1.0 - 16.0*ZL)**-0.25
      ELSEIF ( RICH_NO < 0.16667 ) THEN
         ZL = MIN( RICH_NO / ( 1.0 - 5.0*RICH_NO ), REAL(1.0,ReKi) )  ! The MIN() will take care of rounding issues.
         !PhiM = (1.0 + 5.0*ZL)
      ELSE
         ZL = 1.0
      ENDIF

   ENDIF !TurbModels

   ZL = MIN( ZL, REAL(1.0,ReKi) )

      ! ***** Calculate M-O length scale, L [meters] *****
      ! L should be constant in the surface layer

   IF ( ZL /= 0.0 ) THEN  
        L = HubHt / ZL ! Since ZL is the average ZL over the rotor disk, we should use HubHt to estimate L instead
   ELSE
      L = HUGE( L )
   ENDIF

      ! ***** Calculate power law exponent, if needed *****

   IF ( getPLExp ) THEN
      PLExp = PowerLawExp( RICH_NO )
   ENDIF

      ! ------------ Read in the shear/friction velocity, Ustar, first calculating UstarDiab ------------------------

         ! Set up the heights for the zl- and ustar-profile averages across the rotor disk
      TmpZary  = (/ HubHt-RotorDiameter/2., HubHt, HubHt+RotorDiameter/2. /)
      IF (TmpZary(3) .GE. profileZmin .AND. TmpZary(1) .LE. profileZmax ) THEN  !set height limits so we don't extrapolate too far
         DO TmpIndex = 1,3
            TmpZary(TmpIndex) = MAX( MIN(TmpZary(TmpIndex), profileZmax), profileZmin)
         ENDDO
      ENDIF
   
   
   UseDefault = .TRUE.

   UstarDiab  = getUstarDiab(URef, RefHt)
   Ustar      = UstarDiab

   SELECT CASE ( TurbModel )

      CASE ('WF_UPW')

         IF ( ZL < 0.0 ) THEN
            Ustar = 1.162 * UstarDiab**( 2.0 / 3.0 )
         ELSE ! Include the neutral case to avoid strange discontinuities
            Ustar = 0.911 * UstarDiab**( 2.0 / 3.0 )
         ENDIF

      CASE ( 'WF_07D', 'WF_14D'  )

         IF ( ZL < 0.0 ) THEN
            Ustar = 1.484 * UstarDiab**( 2.0 / 3.0 )
         ELSE ! Include the neutral case with the stable one to avoid strange discontinuities
            Ustar = 1.370 * UstarDiab**( 2.0 / 3.0 )
         ENDIF

      CASE ('GP_LLJ')
         TmpUstarD  = -1.0
         IF ( URef < 0 ) THEN ! (1) We can't get a wind speed because Uref was default
            UseDefault = .FALSE. ! We'll calculate the default value later              
         ELSE
            Ustar = 0.17454 + 0.72045*UstarDiab**1.36242
         ENDIF
                                           
      CASE ( 'NWTCUP' )
         UStar = 0.2716 + 0.7573*UstarDiab**1.2599

      CASE ( 'TIDAL', 'RIVER' )
         ! Use a constant drag coefficient for the HYDRO spectral models.
         UStar = Uref*0.05 ! This corresponds to a drag coefficient of 0.0025.
         !UStar = Uref*0.04 ! This corresponds to a drag coefficient of 0.0016.

   END SELECT

   CALL ReadRVarDefault( UI, InFile, UStar, "the friction or shear velocity", "Friction or shear velocity [m/s]", UseDefault )

   IF ( Uref < 0.0 .AND. UseDefault ) THEN  ! This occurs if "default" was entered for both GP_LLJ wind speed and UStar
      CALL TS_Abort( 'The reference wind speed and friction velocity cannot both be "default."')
   ELSEIF (UStar <= 0) THEN
      CALL TS_Abort( 'The friction velocity must be a positive number.')
   ENDIF            


         ! ***** Calculate wind speed at hub height *****

   IF ( WindProfileType(1:1) == 'J' ) THEN
      IF ( ZJetMax < 0 ) THEN ! Calculate a default value
         ZJetMax = -14.820561*Rich_No + 56.488123*zl + 166.499069*uStar + 188.253377
         ZJetMax = 1.9326*ZJetMax - 252.7267  ! Correct with the residual

         CALL RndJetHeight( tmp ) ! Add a random amount

         ZJetMax = MIN( MAX(ZJetMax + tmp, REAL(70.,ReKi) ), REAL(490.,ReKi) )
      ENDIF

      IF ( URef < 0 ) THEN ! Calculate a default value

         UJetMax = MAX( -21.5515 + 6.6827*LOG(ZJetMax), REAL(5.0,ReKi) ) !Jet max must be at least 5 m/s (occurs ~50 m); shouldn't happen, but just in case....

         CALL Rnd3ParmNorm( tmp, REAL(0.1076,ReKi), REAL(-0.1404,ReKi), REAL(3.6111,ReKi),  REAL(-15.,ReKi), REAL(20.,ReKi) )

         IF (UJetMax + tmp > 0 ) UJetMax = UJetMax + tmp

         CALL GetChebCoefs( UJetMax, ZJetMax ) ! These coefficients are a function of UJetMax, ZJetMax, RICH_NO, and uStar        

         URef = getWindSpeed(UJetMax, ZJetMax, RefHt, RotorDiameter, PROFILE=WindProfileType)

      ELSE
         CALL GetChebCoefs(URef, RefHt)
      ENDIF

   ENDIF !Jet wind profile

   UHub = getWindSpeed(URef, RefHt, HubHt, RotorDiameter, PROFILE=WindProfileType)
   
         ! ***** Get uStar- and zl-profile values, if required, and determine offsets *****
      IF ( TmpUstarD == 0.0 ) THEN
         TmpUstarD = Ustar
      ELSE
         IF ( TmpUstarD < 0.0 ) THEN  ! If > 0, we've already calculated these things...
            UstarDiab = getUstarDiab(URef, RefHt) !bjj: is this problematic for anything else?
            TmpUary   = getWindSpeed(URef, RefHt, TmpZary, RotorDiameter, PROFILE=WindProfileType)
            TmpUstar  = getUstarARY( TmpUary,     TmpZary )         
            TmpUstarD = SUM(TmpUstar) / SIZE(TmpUstar)    ! The average of those 3 points
         ENDIF
         UstarOffset = Ustar - TmpUstarD
         TmpUstar(:) = TmpUstar(:) + UstarOffset
      ENDIF  
      
      TmpZLary = getZLARY(TmpUary, TmpZary)
      zlOffset = zL - SUM(TmpZLary) / SIZE(TmpUstar) 
              

      ! ------------- Read in the mixing layer depth, ZI ---------------------------------------------

   UseDefault = .TRUE.
   IF ( ZL >= 0.0 .AND. TurbModel /= 'GP_LLJ' ) THEN  !We must calculate ZI for stable GP_LLJ. z/L profile can change signs so ZI must be defined for spectra.
      ZI = 0.0
   ELSE
      IF ( Ustar < UstarDiab ) THEN
         ZI = ( 0.04 * Uref ) / ( 1.0E-4 * LOG10( RefHt / Z0 ) )  !for "very" windy days
      ELSE
         !Should Wind Farm models use the other definition since that was what was used in creating those models?
         ZI = Ustar / (6.0 * Fc)
      ENDIF
   ENDIF

   CALL ReadRVarDefault( UI, InFile, ZI, "the mixing layer depth", "Mixing layer depth [m]", UseDefault, &
                                                                  IGNORE=(ZL>=0. .AND. TurbModel /= 'GP_LLJ') )

      IF ( ( ZL < 0.0 ) .AND. ( ZI <= 0.0 ) ) THEN
         CALL TS_Abort ( 'The mixing layer depth must be a positive number for unstable flows.')
      ENDIF


      ! Get the default mean Reynolds stresses
   UWskip     = .FALSE.
   UVskip     = .FALSE.
   VWskip     = .FALSE.   
   PC_UW      = TmpUstar(2)**2   ! Used only for GP-LLJ in GetDefaultRS()
   
   CALL GetDefaultRS(  PC_UW, PC_UV, PC_VW )


       ! ----------- Read in the mean hub u'w' Reynolds stress, PC_UW ---------------------------------------------          
   UseDefault = .TRUE.
   CALL ReadRVarDefault( UI, InFile, PC_UW, "the mean hub u'w' Reynolds stress", &
                                            "Mean hub u'w' Reynolds stress", UseDefault, IGNORESTR = UWskip )
   IF (.NOT. UWskip) THEN  
      TmpUstarD = ( TmpUstar(1)- 2.0*TmpUstar(2) + TmpUstar(3) )
   
      IF ( TmpUstarD /= 0.0 ) THEN
         UstarSlope  = 3.0*(Ustar -  SQRT( ABS(PC_UW) ) ) / TmpUstarD
         UstarOffset = SQRT( ABS(PC_UW) ) - UstarSlope*(TmpUstar(2) - UstarOffset)         
      ELSE
         UstarSlope  = 0.0
         UstarOffset = SQRT( ABS(PC_UW) )
      ENDIF
   ENDIF   

      ! ------------ Read in the mean hub u'v' Reynolds stress, PC_UV ---------------------------------------------
   UseDefault = .TRUE.
   CALL ReadRVarDefault( UI, InFile, PC_UV, "the mean hub u'v' Reynolds stress", &
                                            "Mean hub u'v' Reynolds stress", UseDefault, IGNORESTR = UVskip )

      ! ------------ Read in the mean hub v'w' Reynolds stress, PC_VW ---------------------------------------------
   UseDefault = .TRUE.
   CALL ReadRVarDefault( UI, InFile, PC_VW, "the mean hub v'w' Reynolds stress", &
                                            "Mean hub v'w' Reynolds stress", UseDefault, IGNORESTR = VWskip )

      ! ------------ Read in the u component coherence decrement (coh-squared def), InCDec(1) = InCDecU ------------
   CALL GetDefaultCoh( UHub, HubHt )

   UseDefault = .TRUE.
   InCVar(1) = InCDec(1)
   InCVar(2) = InCohB(1)

   CALL ReadRAryDefault( UI, InFile, InCVar, "the u-component coherence parameters", &
                                             "u-component coherence parameters", UseDefault)
   InCDec(1) = InCVar(1)
   InCohB(1) = InCVar(2)

      IF ( InCDec(1) <= 0.0 ) THEN
         CALL TS_Abort ( 'The u-component coherence decrement must be a positive number.')
      ENDIF

      ! ------------ Read in the v component coherence decrement (coh-squared def), InCDec(2) = InCDecV ----------

   UseDefault = .TRUE.
   InCVar(1) = InCDec(2)
   InCVar(2) = InCohB(2)

   CALL ReadRAryDefault( UI, InFile, InCVar, "the v-component coherence parameters", &
                                             "v-component coherence parameters", UseDefault)

   InCDec(2) = InCVar(1)
   InCohB(2) = InCVar(2)

      IF ( InCDec(2) <= 0.0 ) THEN
         CALL TS_Abort ( 'The v-component coherence decrement must be a positive number.')
      ENDIF

      ! ------------ Read in the w component coherence decrement (coh-squared def), InCDec(3) = InCDecW -------

   UseDefault = .TRUE.
   InCVar(1) = InCDec(3)
   InCVar(2) = InCohB(3)

   CALL ReadRAryDefault( UI, InFile, InCVar, "the w-component coherence parameters", &
                                             "w-component coherence parameters", UseDefault)

   InCDec(3) = InCVar(1)
   InCohB(3) = InCVar(2)

      IF ( InCDec(3) <= 0.0 ) THEN
         CALL TS_Abort ( 'The w-component coherence decrement must be a positive number.')
      ENDIF

         ! ------------ Read in the coherence exponent, COHEXP -----------------------------------

   UseDefault = .TRUE.
   CohExp     = 0.0    ! was 0.25
   CALL ReadRVarDefault( UI, InFile, CohExp, "the coherence exponent", "Coherence exponent", UseDefault)

      IF ( COHEXP < 0.0 ) THEN
         CALL TS_Abort ( 'The coherence exponent must be non-negative.')
      ENDIF


         !=================================================================================
         ! Read the Coherent Turbulence Scaling Parameters, if necessary.  !
         !=================================================================================

   IF ( WrACT ) THEN

      FormStr = "( // 'Coherent Turbulence Scaling Parameters:' / )"
      WRITE (US,FormStr)
      !READ (UI,'(/)')                        ! Read header line
      CALL ReadCom( UI, InFile, "Coherent Turbulence Scaling Parameters Heading Line 1" )
      CALL ReadCom( UI, InFile, "Coherent Turbulence Scaling Parameters Heading Line 2" )


         ! ------------ Read the name of the path containg event file definitions, CTEventPath --------------------------

      CALL ReadCVar( UI, InFile, CTEventPath, "the path of the coherent turbulence events", "Coherence events path")

      IF ( LEN( TRIM(CTEventPath) ) <= 10 )  THEN
         FormStr = "( A10 , 2X , 'Name of the path containing the coherent turbulence data files' )"
      ELSE
         FormStr = "( A, /, 12X , 'Name of the path containing the coherent turbulence data files' )"
      ENDIF

      WRITE (US,FormStr)  TRIM(CTEventPath)

      CALL ReadCVar( UI, InFile, Line, "the event file type", "Event file type")


      IF ( KHtest ) THEN

         CText = 'les'
         CTEventFile = TRIM(CTEventPath)//'\Events.xtm'

         CALL WrScr( ' LES events will be used for the KH test.' )

      ELSE

         CText = Line  !This will preserve the case formatting, in case it matters.

         CALL Conv2UC( Line )

         IF (Line(1:6) == "RANDOM") THEN
             CALL RndUnif( tmp )

             IF ( tmp <= 0.5 ) THEN
                 CText = 'les'
             ELSE
                 CText = 'dns'
             ENDIF
         ENDIF

         CTEventFile = TRIM(CTEventPath)//'\Events.'//TRIM(CText)

      ENDIF

      FormStr = "( 7X, A3, 2X, 'Type of coherent turbulence data files' )"

      WRITE (US,FormStr) TRIM(CText)


         ! ------------ Read the Randomization Flag, Randomize -----------------------------------

      CALL ReadLVar( UI, InFile, Randomize, "the randomization flag", "Randomize CT Scaling")


      IF ( KHtest ) THEN
         Randomize = .FALSE.
         CALL WrScr( ' Billow will cover rotor disk for KH test. ' )
      ENDIF

      FormStr = "( L10 , 2X , 'Randomize the disturbance scale and location?' )"
      WRITE (US,FormStr)  Randomize


         ! ------------ Read the Disturbance Scale, DistScl ---------------------------------------------

      CALL ReadRVar( UI, InFile, DistScl, "the disturbance scale", "Disturbance scale")


         ! ------------ Read the Lateral Fractional Location of tower centerline in wave, CTLy ----------

      CALL ReadRVar( UI, InFile, CTLy, "the fractional location of tower centerline from right", &
                           "Location of tower centerline")

         ! ------------ Read the Vertical Fraction Location of hub in wave, CTLz ------------------------

      CALL ReadRVar( UI, InFile, CTLz, "the fractional location of hub height from the bottom", &
                     "Location of hub height")

      IF ( KHtest ) THEN
            DistScl = 1.0
            CTLy    = 0.5
            CTLz    = 0.5
      ELSEIF ( Randomize ) THEN

         CALL RndUnif( tmp )

            ! Assume a 75% chance of coherent turbulence being the size of the rotor
            ! If the rotor is too small, assume a 100% chance.
            ! If the turbulence is not the size of the rotor, assume it's half the size
            ! of the disk, with equal probablilty of being in the lower or upper half.

         IF ( tmp > 0.25 .OR. RotorDiameter <= 30.0 ) THEN

            DistScl = 1.0
            CTLy    = 0.5
            CTLz    = 0.5

         ELSE

            DistScl = 0.5
            CTLy    = 0.5

            IF ( tmp < 0.125 ) THEN
               CTLz = 0.0 ! The hub is on the bottom of the dataset (i.e. the billow is on the top of the disk)
            ELSE
               CTLz = 1.0 ! The hub is on the top of the dataset
            ENDIF

         ENDIF

      ELSE  !Don't randomize:

         IF ( DistScl < 0.0 ) THEN
            CALL TS_Abort ('The disturbance scale must be a positive.')
         ELSEIF ( RotorDiameter <= 30.0 .AND. DistScl < 1.0 ) THEN
            CALL TS_Abort ('The disturbance scale must be at least 1.0 for rotor diameters less than 30.')
         ELSEIF ( RotorDiameter*DistScl <= 15.0  ) THEN
            CALL TS_Abort ('The coherent turbulence must be greater than 15 meters in height.  '//&
                        'Increase the rotor diameter or the disturbance scale. ')
         ENDIF

      ENDIF


      FormStr = "( F10.3 , 2X , 'Disturbance scale (ratio of wave height to rotor diameter)' )"
      WRITE (US,FormStr)  DistScl

      FormStr = "( F10.3 , 2X , 'Fractional location of tower centerline from right' )"
      WRITE (US,FormStr)  CTLy

      FormStr = "( F10.3 , 2X , 'Fractional location of hub height from the bottom of the dataset' )"
      WRITE (US,FormStr)  CTLz

         ! ---------- Read the Minimum event start time, CTStartTime --------------------------------------------

      CALL ReadRVar( UI, InFile, CTStartTime, "the minimum start time for coherent structures", "CTS Start Time")

      CTStartTime = MAX( CTStartTime, REAL(0.0,ReKi) ) ! A Negative start time doesn't really make sense...

      FormStr = "( F10.3 , 2X , 'Minimum start time for coherent structures [seconds]' )"
      WRITE (US,FormStr)  CTStartTime

   ENDIF    ! WrACT


ELSE  ! IECVKM or IECKAI models

   RICH_NO = 0.0                       ! Richardson Number in neutral conditions
   COHEXP  = 0.0                       ! Coherence exponent

   IF ( IECedition == 3 ) THEN
      IncDec = (/ 12.00, 0.0, 0.0 /)   ! u-, v-, and w-component coherence decrement for IEC Ed. 3
   ELSE
      IncDec = (/  8.80, 0.0, 0.0 /)   ! u-, v-, and w-component coherence decrement
   ENDIF


      ! The following variables are not used in the IEC calculations

   ZL      = 0.0                       ! M-O z/L parameter
   L       = 0.0                       ! M-O length scale
   Ustar   = 0.0                       ! Shear or friction velocity
   ZI      = 0.0                       ! Mixing layer depth
   PC_UW   = 0.0                       ! u'w' x-axis correlation coefficient
   PC_UV   = 0.0                       ! u'v' x-axis correlation coefficient
   PC_VW   = 0.0                       ! v'w' x-axis correlation coefficient
   DistScl = 0.0                       ! coherent turbulence disturbance scale
   CTLy    = 0.0                       ! coherent turbulence scaling values
   CTLz    = 0.0                       ! coherent turbulence scaling values
   Ym_max  = 0.0                       ! coherent turbulence scaling values
   Zm_max  = 0.0                       ! coherent turbulence scaling values

   IF ( NumTurbInp .AND. PerTurbInt == 0.0 ) THEN    ! This will produce constant winds, instead of an error when the transfer matrix is singular
      TurbModel = 'NONE'
   ENDIF

      ! Calculate wind speed at hub height

   UHub    = getWindSpeed(URef, RefHt, HubHt, RotorDiameter, PROFILE=WindProfileType)


ENDIF


   ! Done reading the input file.

CLOSE (UI)


RETURN
END SUBROUTINE GetInput
!=======================================================================
FUNCTION getUstarARY(WS, Ht)

   USE                                  TSMods, ONLY: UstarDiab      ! Diabatic u*0
   USE                                  TSMods, ONLY: RICH_NO        ! Richardson number
   USE                                  TSMods, ONLY: UstarOffset
   USE                                  TSMods, ONLY: UstarSlope
   USE                                  TSMods, ONLY: profileZmax
   USE                                  TSMods, ONLY: profileZmin
   

   IMPLICIT                              NONE

   REAL(ReKi),   INTENT(IN)           :: Ht(:)                       ! Height at which ustar is defined
   REAL(ReKi),   INTENT(IN)           :: WS(:)                       ! Wind speed(s) at heights, Ht
   
   REAL(ReKi)                         :: tmpZ                        ! a temporary value
   REAL(ReKi)                         :: getUstarARY(SIZE(Ht))       ! the array of ustar values
   
   INTEGER                            :: IZ
   INTEGER                            :: Zindx
   INTEGER                            :: Zindx_mn (1)
   INTEGER                            :: Zindx_mx (1)
   
   LOGICAL                            :: mask(SIZE(Ht))
   
   mask = Ht.GE.profileZmin
   IF ( ANY(mask) ) THEN
      Zindx_mn = MINLOC( Ht, MASK=mask ) 
            
      mask = Ht.LE.profileZmax
      IF ( ANY(mask) ) THEN
         Zindx_mx = MAXLOC( Ht, MASK=mask )    
            
         DO IZ = 1,SIZE(Ht)
            IF ( Ht(IZ) < profileZmin ) THEN
               Zindx = Zindx_mn(1)
            ELSEIF ( Ht(IZ) > profileZmax ) THEN
               Zindx = Zindx_mx(1)
            ELSE
               Zindx = IZ
            ENDIF
            
            tmpZ = Ht(Zindx)      !ustar is constant below 50 meters, and we don't want to extrapolate too high (last measurement is at 116 m)
               
            getUstarARY(  IZ) = ( 0.045355367 +  4.47275E-8*tmpZ**3)                                                &
                              + ( 0.511491978 -  0.09691157*LOG(tmpZ) - 199.226951/tmpZ**2           ) * WS(Zindx)  &
                              + (-0.00396447  - 55.7818832/tmpZ**2                                   ) * RICH_NO    &
                              + (-5.35764429  +  0.102002162*tmpZ/LOG(tmpZ) + 25.30585136/SQRT(tmpZ) ) * UstarDiab
         ENDDO         
            
      ELSE ! All are above the max height so we'll use the old relationship at all heights
         getUstarARY(:) = 0.17454 + 0.72045*UstarDiab**1.36242
      ENDIF
            
   ELSE ! All are below the min height so we'll use the diabatic Ustar value 
      getUstarARY(:) = UstarDiab
   ENDIF
   
   getUstarARY = UstarSlope * getUstarARY(:) + UstarOffset  ! These terms are used to make the ustar profile match the rotor-disk averaged value and input hub u'w'
      
END FUNCTION
!=======================================================================
FUNCTION getUstarDiab(URef, RefHt)

   USE                                  TSMods, ONLY: z0             ! Surface roughness length -- It must be > 0 (which we've already checked for)
   USE                                  TSMods, ONLY: ZJetMax        ! Height of the low-level jet
   USE                                  TSMods, ONLY: ZL             ! M-O stability parameter

   IMPLICIT                              NONE

   REAL(ReKi),   INTENT(IN)           :: URef                        ! Wind speed at reference height
   REAL(ReKi),   INTENT(IN)           :: RefHt                       ! Reference height
   
   REAL(ReKi)                         :: tmp                         ! a temporary value
   REAL(ReKi)                         :: psiM                        
   REAL(ReKi)                         :: getUstarDiab                ! the diabatic u* value (u*0)
   
   IF ( ZL >= 0 ) THEN !& ZL < 1
      psiM = -5.0*MIN(ZL, REAL(1.0,ReKi) )
   ELSE
      tmp = (1.0 - 15.0*ZL)**0.25  
            
      !psiM = -2.0*LOG( (1.0 + tmp)/2.0 ) - LOG( (1.0 + tmp*tmp)/2.0 ) + 2.0*ATAN( tmp ) - 0.5 * PI
      psiM = -LOG( 0.125 * ( (1.0 + tmp)**2 * (1.0 + tmp*tmp) ) ) + 2.0*ATAN( tmp ) - 0.5 * PI

   ENDIF

   getUstarDiab = ( 0.4 * Uref ) / ( LOG( RefHt / z0 ) - psiM )   

END FUNCTION
!=======================================================================
SUBROUTINE GetUSR(U_in, FileName, NLines)

   USE                                   TSMods, ONLY: HubHt             ! Hub Height
   USE                                   TSMods, ONLY: L_USR             ! von Karman length scale
   USE                                   TSMods, ONLY: NumUSRz           ! Number of user-defined heights for the profiles
   USE                                   TSMods, ONLY: Sigma_USR         ! standard deviation profile
   USE                                   TSMods, ONLY: StdScale          ! scaling for standard deviation profile
   USE                                   TSMods, ONLY: TurbModel         ! Type of turbulence model
   USE                                   TSMods, ONLY: U_USR             ! wind speed profile
   USE                                   TSMods, ONLY: WindDir_USR       ! wind direction profile
   USE                                   TSMods, ONLY: Z_USR             ! heights corresponding to the profiles
   IMPLICIT                              NONE

   CHARACTER(*), INTENT(IN)           :: FileName                       ! Name of the input file
   CHARACTER(200)                     :: LINE

   REAL(ReKi)                         :: L_Usr_Tmp
   REAL(ReKi)                         :: Sigma_USR_Tmp
   REAL(ReKi)                         :: U_USR_Tmp
   REAL(ReKi)                         :: WindDir_USR_Tmp
   REAL(ReKi)                         :: Z_USR_Tmp

   INTEGER                            :: I
   INTEGER                            :: Indx
   INTEGER                            :: J
   INTEGER                            :: IOAstat                        ! Input/Output/Allocate status
   INTEGER, INTENT(IN), OPTIONAL      :: NLines                         ! Number of lines to be skipped, if the file must be rewound
   INTEGER, INTENT(IN)                :: U_in                           ! Input unit.  This file is assumed to be already open

   LOGICAL                            :: ReadSigL                       ! Whether or not to read the last 2 columns

      ! Find the end of the input file, where the "User-Defined Variables" are located
   READ ( U_in, '(A)', IOSTAT=IOAstat ) LINE

   IF ( IOAstat /= 0 ) THEN
      CALL TS_Abort( 'Could not read entire input file for user-defined variables.' )
   ENDIF

   CALL Conv2UC ( LINE )

   DO WHILE ( INDEX( LINE, 'USER-DEFINED' ) == 0 )

      READ ( U_in, '(A)', IOSTAT=IOAstat ) LINE

      IF ( IOAstat /= 0 ) THEN
         CALL TS_Abort( 'Could not read entire input file for user-defined variables.' )
      ENDIF

      CALL Conv2UC( LINE )

   ENDDO


      ! ---------- Read the size of the arrays --------------------------------------------
   CALL ReadIVar( U_in, FileName, NumUSRz, "NumUSRz", "Number of heights in the user-defined profiles" )

   IF ( NumUSRz < 1 ) THEN
      CALL TS_Abort( 'The number of heights specified in the user-defined profiles must be at least 1.')
   ENDIF

   DO I=1,3
         ! ---------- Read the size of the arrays --------------------------------------------
      CALL ReadRVar( U_in, FileName, StdScale(I), "StdScale", "Scaling value for user-defined standard deviation profile" )

      IF ( StdScale(I) <= 0. ) THEN
         CALL TS_Abort( 'The scaling value for the user-defined standard deviation profile must be positive.')
      ENDIF
   ENDDO

      ! Allocate the data arrays
   ALLOCATE ( Z_USR(NumUSRz) , STAT=IOAstat )

   IF ( IOAstat /= 0 )  THEN
      CALL TS_Abort ( 'Error allocating '//TRIM( Int2LStr( ReKi*NumUSRz/1024**2 ) )//' MB for the user-defined height array.' )
   ENDIF

   ALLOCATE ( U_USR(NumUSRz) , STAT=IOAstat )

   IF ( IOAstat /= 0 )  THEN
      CALL TS_Abort ( 'Error allocating '//TRIM( Int2LStr( ReKi*NumUSRz/1024**2 ) )//' MB for the user-defined wind speed array.' )
   ENDIF

   ALLOCATE ( WindDir_USR(NumUSRz) , STAT=IOAstat )

   IF ( IOAstat /= 0 )  THEN
      CALL TS_Abort ( 'Error allocating '//TRIM( Int2LStr( ReKi*NumUSRz/1024**2 ) )// &
                       ' MB for the user-defined wind direction array.' )
   ENDIF

   IF ( TurbModel == 'USRVKM' ) THEN
      ReadSigL = .TRUE.

      ALLOCATE ( Sigma_USR(NumUSRz) , STAT=IOAstat )

      IF ( IOAstat /= 0 )  THEN
         CALL TS_Abort ( 'Error allocating '//TRIM( Int2LStr( ReKi*NumUSRz/1024**2 ) )//' MB for the user-defined sigma array.' )
      ENDIF

      ALLOCATE ( L_USR(NumUSRz) , STAT=IOAstat )

      IF ( IOAstat /= 0 )  THEN
         CALL TS_Abort ( 'Error allocating '//TRIM( Int2LStr( ReKi*NumUSRz/1024**2 ) )// &
                      ' MB for the user-defined length scale array.' )
      ENDIF
   ELSE
      ReadSigL = .FALSE.
   ENDIF

      ! ---------- Skip 4 lines --------------------------------------------
   DO I=1,4
      CALL ReadCom( U_in, FileName, "Headers for user-defined variables" )
   ENDDO

   DO I=1,NumUSRz

      IF ( ReadSigL ) THEN
         READ( U_in, *, IOSTAT=IOAstat ) Z_USR(I), U_USR(I), WindDir_USR(I), Sigma_USR(I), L_USR(I)
      ELSE
         READ( U_in, *, IOSTAT=IOAstat ) Z_USR(I), U_USR(I), WindDir_USR(I)
      ENDIF

      IF ( IOAstat /= 0 ) THEN
         CALL TS_Abort( 'Could not read entire user-defined variable list on line '//Int2LStr(I)//'.' )
      ENDIF

      IF ( ReadSigL ) THEN
         IF ( Sigma_USR(I) <= REAL( 0., ReKi ) ) THEN
            CALL TS_Abort( 'The standard deviation must be a positive number.' );
         ELSEIF ( L_USR(I) <= REAL( 0., ReKi ) ) THEN
            CALL TS_Abort( 'The length scale must be a positive number.' );
         ENDIF
      ENDIF

      IF ( WindDir_USR(I) > 360. ) THEN
         J = INT ( WindDir_USR(I) / 360. )
         WindDir_USR(I) = WindDir_USR(I) - J * 360.
      ELSEIF ( WindDir_USR(I) < 0. ) THEN
         J = INT ( -WindDir_USR(I) / 360. ) +1
         WindDir_USR(I) = WindDir_USR(I) + J * 360.
      ENDIF
   ENDDO

      ! Sort the arrays
   DO I=2,NumUSRz
      IF ( Z_USR(I) < Z_USR(I-1) ) THEN

         Indx = 1
         DO J=I-2,1,-1
            IF ( Z_USR(I) > Z_USR(J) ) THEN
               Indx = J+1
               EXIT
            ELSEIF ( Z_USR(I) == Z_USR(J) ) THEN
               CALL TS_Abort( 'Error: user-defined values must contain unique heights.' )
            ENDIF
         ENDDO

         Z_USR_Tmp       = Z_USR(I)
         U_USR_Tmp       = U_USR(I)
         WindDir_USR_Tmp = WindDir_USR(I)

         DO J=I,Indx+1,-1
            Z_USR(J)       = Z_USR(J-1)
            U_USR(J)       = U_USR(J-1)
            WindDir_USR(J) = WindDir_USR(J-1)
         ENDDO

         Z_USR(Indx)       = Z_USR_Tmp
         U_USR(Indx)       = U_USR_Tmp
         WindDir_USR(Indx) = WindDir_USR_Tmp

         IF ( ReadSigL ) THEN
            Sigma_USR_Tmp   = Sigma_USR(I)
            L_USR_Tmp       = L_USR(I)

            DO J=I,Indx+1,-1
               Sigma_USR(J)   = Sigma_USR(J-1)
               L_USR(J)       = L_USR(J-1)
            ENDDO

            Sigma_USR(Indx)   = Sigma_USR_Tmp
            L_USR(Indx)       = L_USR_Tmp
         ENDIF ! ReadSigL

      ENDIF
   ENDDO

      ! Rewind the file, if necessary.
   IF ( PRESENT(NLines) ) THEN
      REWIND( U_in , IOSTAT=IOAstat )

      IF ( IOAstat /= 0 ) THEN
         CALL TS_Abort( 'Error rewinding the file '//TRIM(FileName)//'.' )
      ENDIF

      DO I = 1,NLines
         CALL ReadCom( U_in, FileName, "Line "//Int2LStr(I) )
      ENDDO
   ENDIF

END SUBROUTINE GetUSR
!=======================================================================
SUBROUTINE GetUSRSpec(FileName)

   USE                                   TSMods, ONLY: TurbModel         ! Type of turbulence model
   USE                                   TSMods, ONLY: NumUSRf           ! Length of user-defined spectra (# frequencies)
   USE                                   TSMods, ONLY: USpec             ! Unit number for the user-defined spectra file
   USE                                   TSMods, ONLY: Freq_USR          ! Frequencies for the user-defined spectra
   USE                                   TSMods, ONLY: Uspec_USR         ! User-defined u-component spectra
   USE                                   TSMods, ONLY: Vspec_USR         ! User-defined v-component spectra
   USE                                   TSMods, ONLY: Wspec_USR         ! User-defined w-component spectra
   
   IMPLICIT                              NONE

   CHARACTER(*), INTENT(IN)           :: FileName                       ! Name of the input file
   CHARACTER(200)                     :: LINE
   
   REAL(ReKi)                         :: Freq_USR_Tmp
   REAL(ReKi)                         :: U_USR_Tmp
   REAL(ReKi)                         :: V_USR_Tmp
   REAL(ReKi)                         :: W_USR_Tmp
   REAL(ReKi)                         :: SpecScale (3)

   INTEGER                            :: I
   INTEGER                            :: Indx
   INTEGER                            :: J
   INTEGER                            :: IOAstat                        ! Input/Output/Allocate status


      ! --------- Open the file ---------------

   CALL OpenFInpFile( USpec, FileName )

   CALL WrScr1(' Reading the user-defined spectra input file "'//TRIM(FileName)//'".' )


      ! --------- Read the comment lines at the beginning of the file ---------------
   DO I=1,3
      READ ( USpec, '(A)', IOSTAT=IOAstat ) LINE

      IF ( IOAstat /= 0 ) THEN
         CALL TS_Abort( 'Could not read entire input file for user-defined spectra.' )
      ENDIF
   ENDDO


      ! ---------- Read the size of the arrays --------------------------------------------
   CALL ReadIVar( USpec, FileName, NumUSRf, "NumUSRf", "Number of frequencies in the user-defined spectra" )

   IF ( NumUSRf < 3 ) THEN
      CALL TS_Abort( 'The number of frequencies specified in the user-defined spectra must be at least 3.')
   ENDIF

   DO I=1,3
         ! ---------- Read the scaling for the arrays --------------------------------------------
      CALL ReadRVar( USpec, FileName, SpecScale(I), "SpecScale", "Scaling value for user-defined standard deviation profile" )

      IF ( SpecScale(I) <= 0. ) THEN
         CALL TS_Abort( 'The scaling value for the user-defined spectra must be positive.')
      ENDIF
   ENDDO

      ! Allocate the data arrays
   ALLOCATE ( Freq_USR(NumUSRf) , STAT=IOAstat )

   IF ( IOAstat /= 0 )  THEN
      CALL TS_Abort ( 'Error allocating '//TRIM( Int2LStr( ReKi*NumUSRf/1024**2 ) )//' MB for the user-defined frequency array.')
   ENDIF

   ALLOCATE ( Uspec_USR(NumUSRf) , STAT=IOAstat )

   IF ( IOAstat /= 0 )  THEN
      CALL TS_Abort ( 'Error allocating '//TRIM( Int2LStr( ReKi*NumUSRf/1024**2 ) )//' MB for the user-defined u spectra array.')
   ENDIF

   ALLOCATE ( Vspec_USR(NumUSRf) , STAT=IOAstat )

   IF ( IOAstat /= 0 )  THEN
      CALL TS_Abort ( 'Error allocating '//TRIM( Int2LStr( ReKi*NumUSRf/1024**2 ) )//' MB for the user-defined v spectra array.')
   ENDIF
   
   ALLOCATE ( Wspec_USR(NumUSRf) , STAT=IOAstat )

   IF ( IOAstat /= 0 )  THEN
      CALL TS_Abort ( 'Error allocating '//TRIM( Int2LStr( ReKi*NumUSRf/1024**2 ) )//' MB for the user-defined w spectra array.')
   ENDIF
   


      ! ---------- Skip 4 lines --------------------------------------------
   DO I=1,4
      CALL ReadCom( USpec, FileName, "Headers for user-defined variables" )
   ENDDO

      ! ---------- Read the data lines --------------------------------------
   DO I=1,NumUSRf

      READ( USpec, *, IOSTAT=IOAstat ) Freq_USR(I), Uspec_USR(I), Vspec_USR(I), Wspec_USR(I)

      IF ( IOAstat /= 0 ) THEN
         CALL TS_Abort( 'Could not read entire user-defined spectra on line '//Int2LStr(I)//'.' )
      ENDIF

      IF ( ( Uspec_USR(I) <= REAL( 0., ReKi ) ) .OR. &
           ( Vspec_USR(I) <= REAL( 0., ReKi ) ) .OR. & 
           ( Wspec_USR(I) <= REAL( 0., ReKi ) ) ) THEN
         CALL TS_Abort( 'The spectra must contain positive numbers.' );
!      ELSEIF ( Freq_USR(I) <= REAL( 0., ReKi ) ) THEN
!         CALL TS_Abort( 'The frequencies must be a positive number.' );
      ENDIF
      
         ! Scale by the factors earlier in the input file
      
      Uspec_USR(I) = Uspec_USR(I)*SpecScale(1)
      Vspec_USR(I) = Vspec_USR(I)*SpecScale(2)
      Wspec_USR(I) = Wspec_USR(I)*SpecScale(3)
      
   ENDDO

      ! ------- Sort the arrays by frequency -----------------------------------
   DO I=2,NumUSRf
      IF ( Freq_USR(I) < Freq_USR(I-1) ) THEN

         Indx = 1
         DO J=I-2,1,-1
            IF ( Freq_USR(I) > Freq_USR(J) ) THEN
               Indx = J+1
               EXIT
            ELSEIF ( Freq_USR(I) == Freq_USR(J) ) THEN
               CALL TS_Abort( 'Error: user-defined spectra must contain unique frequencies.' )
            ENDIF
         ENDDO

         Freq_USR_Tmp    = Freq_USR(I)
         U_USR_Tmp       = Uspec_USR(I)
         V_USR_Tmp       = Vspec_USR(I)
         W_USR_Tmp       = Wspec_USR(I)

         DO J=I,Indx+1,-1
            Freq_USR(J)    = Freq_USR(J-1)
            Uspec_USR(J)   = Uspec_USR(J-1)
            Vspec_USR(J)   = Vspec_USR(J-1)
            Wspec_USR(J)   = Wspec_USR(J-1)
         ENDDO

         Freq_USR(Indx)    = Freq_USR_Tmp
         Uspec_USR(I)      = U_USR_Tmp
         Vspec_USR(I)      = V_USR_Tmp
         Wspec_USR(I)      = W_USR_Tmp

      ENDIF
   ENDDO

      ! --------- Close the file ---------------------------------------

   CLOSE( USpec )      

END SUBROUTINE GetUSRSpec
!=======================================================================

FUNCTION getWindSpeedAry(URef, RefHt, Ht, RotorDiam, profile, UHangle)

   ! Determine the wind speed at a given height, with reference wind speed.
   ! Use power law if given height and reference height are within rotor disk.
   ! Use log profile if reference height is below rotor disk.

   USE                                  TSMods, ONLY: ChebyCoef_WD   ! Chebyshev coefficients (must be defined before calling this function!)
   USE                                  TSMods, ONLY: ChebyCoef_WS   ! Chebyshev coefficients (must be defined before calling this function!)
   USE                                  TSMods, ONLY: HFlowAng       ! The horizontal flow angle of the mean wind speed at hub height
   USE                                  TSMods, ONLY: HubHt          ! The hub height (used with HFlowAng)
   USE                                  TSMods, ONLY: IEC_EWM1       ! IEC Extreme 1-yr wind speed model
   USE                                  TSMods, ONLY: IEC_EWM50      ! IEC Extreme 50-yr wind speed model
   USE                                  TSMods, ONLY: IEC_WindType   ! Type of IEC wind
   USE                                  TSMods, ONLY: NumUSRz        ! Number of user-defined heights
   USE                                  TSMods, ONLY: PLExp          ! Power law exponent
   USE                                  TSMods, ONLY: U_Usr          ! User-defined wind speeds
   USE                                  TSMods, ONLY: Ustar          ! Friction or shear velocity
   USE                                  TSMods, ONLY: Vref           ! IEC Extreme wind value
   USE                                  TSMods, ONLY: WindDir_USR    ! User-defined wind directions
   USE                                  TSMods, ONLY: Z_Usr          ! User-defined heights
   USE                                  TSMods, ONLY: z0             ! Surface roughness length -- It must be > 0 (which we've already checked for)
   USE                                  TSMods, ONLY: ZJetMax        ! Height of the low-level jet
   USE                                  TSMods, ONLY: ZL             ! M-O stability parameter

   IMPLICIT                              NONE

   REAL(ReKi),   INTENT(IN)           :: URef                        ! Wind speed at reference height
   REAL(ReKi),   INTENT(IN)           :: RefHt                       ! Reference height
   REAL(ReKi),   INTENT(IN)           :: Ht(:)                       ! Height where wind speed should be calculated
   REAL(ReKi),   INTENT(IN)           :: RotorDiam                   ! Diameter of rotor disk (meters)
   REAL(ReKi),   INTENT(OUT),OPTIONAL :: UHangle(SIZE(Ht))           ! Horizontal wind angle
   REAL(ReKi)                         :: getWindSpeedAry(SIZE(Ht))   ! This function, approximate wind speed at Ht
   
   REAL(SiKi),   PARAMETER            :: MinZ = 3.                   ! lower bound (height) for Cheby polynomial
   REAL(SiKi),   PARAMETER            :: MaxZ = 500.                 ! upper bound (height) for Cheby polynomial

   CHARACTER(*), INTENT(IN), OPTIONAL :: profile                     ! String that determines what profile is to be used
   CHARACTER(3)                       :: profile_type

   REAL(ReKi)                         :: psiM                        ! The diabatic term for the log wind profile
   REAL(ReKi)                         :: tmp                         ! A temporary variable for calculating psiM
   REAL(ReKi)                         :: tmpHt(2)
   REAL(ReKi)                         :: tmpWS(2)

   INTEGER                            :: I
   INTEGER                            :: Indx
   INTEGER                            :: J

   IF ( IEC_WindType == IEC_EWM50 ) THEN
      getWindSpeedAry(:) = VRef*( Ht(:)/HubHt )**0.11                ! [IEC 61400-1 6.3.2.1 (14)]  !bjj: this PLExp should be set to 0.11 already, why is this hard-coded?
      RETURN
   ELSEIF ( IEC_WindType == IEC_EWM1 ) THEN
      getWindSpeedAry(:) = 0.8*VRef*( Ht(:)/HubHt )**0.11            ! [IEC 61400-1 6.3.2.1 (14), (15)]
      RETURN
   ENDIF


   IF ( PRESENT( profile ) ) THEN
      profile_type = profile
      CALL Conv2UC( profile_type )
   ELSE
      profile_type = 'IEC'
   ENDIF

   SELECT CASE ( TRIM(profile_type) )

      CASE ( 'JET', 'J' )

         CALL ChebyshevVals( ChebyCoef_WS, Ht, getWindSpeedAry, MinZ, MaxZ ) ! We originally calculated the coeffs for 3-500 m in height

         IF ( PRESENT( UHangle ) ) THEN

               ! Calculate the wind direction at this height
            CALL ChebyshevVals( ChebyCoef_WD, Ht(:), UHangle(:), MinZ, MaxZ )

               ! Compute the wind direction at hub height & the jet height
            tmpHt(1) = HubHt
            tmpHt(2) = ZJetMax
            CALL ChebyshevVals( ChebyCoef_WD, tmpHt, tmpWS(1:2), MinZ, MaxZ )

               ! Make sure none of the directions are more than 45 degrees from the direction at the jet height
            IF ( ABS(tmpWS(1) - tmpWS(2) ) > 45. ) THEN  ! The direction at the hub height
               tmpWS(1) = tmpWS(2) + SIGN(REAL(45.,ReKi), tmpWS(1) - tmpWS(2))
            ENDIF

            DO I = 1,SIZE(UHangle) ! The directions at all the heights
               IF ( ABS(UHangle(I) - tmpWS(2) ) > 45. ) THEN
                  UHangle(I) = tmpWS(2) + SIGN(REAL(45.,ReKi), UHangle(I) - tmpWS(2))
               ENDIF

               ! Remove the hub height direction so that we have a relative direction, then
               ! add the mean flow angle. (Note that the Chebyshev profile is cw looking upwind,
               ! but the horizontal angle is ccw looking upwind)

               UHangle(I) = HFlowAng - (UHangle(I) - tmpWS(1)) ! This is the counter-clockwise angle of the wind
            ENDDO

         ENDIF


      CASE ( 'LOG', 'L' )

!         IF ( Ht > 0.0 .AND. RefHt > 0.0 .AND. RefHt /= Z0 ) THEN

            IF ( ZL >= 0 ) THEN !& ZL < 1
               psiM = -5.0*ZL
            ELSE
               tmp = (1.0 - 15.0*ZL)**0.25

               !psiM = -2.0*LOG( (1.0 + tmp)/2.0 ) - LOG( (1.0 + tmp*tmp)/2.0 ) + 2.0*ATAN( tmp ) - 0.5 * PI
               psiM = -LOG( 0.125 * ( (1.0 + tmp)**2 * (1.0 + tmp*tmp) ) ) + 2.0*ATAN( tmp ) - 0.5 * PI
            ENDIF

!            IF ( Ustar > 0. ) THEN
!               getWindSpeedAry(:) = ( UstarDiab / 0.4 ) * ( LOG( Ht(:) / Z0 ) - psiM )
!            ELSE
               !In neutral conditions, psiM is 0 and we get the IEC log wind profile:
               getWindSpeedAry(:) = URef*( LOG( Ht(:) / Z0 ) - psiM )/( LOG( RefHt / Z0 ) - psiM )
!            ENDIF
               
!         ENDIF
      CASE ( 'H2L', 'H' )

               ! Calculate the windspeed.
               ! RefHt and URef both get modified consistently, therefore RefHt is used instead of H_ref.
               !print *, RefHt,H_ref,URef,Ustar ! Need to include H_ref from TSmods for this print statement to work.
               getWindSpeedAry(:) = LOG(Ht(:)/RefHt)*Ustar/0.41+URef
            
      CASE ( 'PL', 'P' )

!         IF ( RefHt > 0.0 .AND. Ht > 0.0 ) THEN
            getWindSpeedAry(:) = URef*( Ht(:)/RefHt )**PLExp
!         ENDIF

      CASE ( 'USR', 'U' )

         DO J = 1,SIZE(Ht)
            IF ( Ht(J) <= Z_USR(1) ) THEN
               getWindSpeedAry(J) = U_USR(1)
            ELSEIF ( Ht(J) >= Z_USR(NumUSRz) ) THEN
               getWindSpeedAry(J) = U_USR(NumUSRz)
            ELSE
               ! Find the two points between which the height lies

               DO I=2,NumUSRz
                  IF ( Ht(J) <= Z_USR(I) ) THEN
                     Indx = I-1

                     ! Let's just do a linear interpolation for now
                     getWindSpeedAry(J) = (Ht(J) - Z_USR(Indx)) * ( U_USR(Indx) - U_USR(I) ) / ( Z_USR(Indx) - Z_USR(I) ) &
                                        + U_USR(Indx)
                     EXIT
                  ENDIF
               ENDDO

            ENDIF

            IF ( PRESENT( UHangle ) ) THEN
                  ! Calculate the wind direction at this height

               IF ( Ht(J) <= Z_USR(1) ) THEN
                  UHangle(J) = WindDir_USR(1)
               ELSEIF ( Ht(J) >= Z_USR(NumUSRz) ) THEN
                  UHangle(J) = WindDir_USR(NumUSRz)
               ELSE
                  I = Indx + 1

                     ! Let's just do a linear interpolation for now
                  !we need to check if the angle goes through 360, before we do the interpolation
                  tmpWS(1) = WindDir_USR(Indx) - WindDir_USR(I)
                  IF ( tmpWS(1) > 180. ) THEN
                     tmpWS(1) = WindDir_USR(Indx)
                     tmpWS(2) = WindDir_USR(I   ) + 360.
                  ELSEIF ( tmpWS(1) < -180. ) THEN
                     tmpWS(1) = WindDir_USR(Indx) + 360.
                     tmpWS(2) = WindDir_USR(I   )
                  ELSE
                     tmpWS(1) = WindDir_USR(Indx)
                     tmpWS(2) = WindDir_USR(I   )
                  ENDIF

                 UHangle(J) = (Ht(J) - Z_USR(Indx)) * ( tmpWS(1) - tmpWS(2) ) / ( Z_USR(Indx) - Z_USR(I) ) + tmpWS(1)

               ENDIF

               UHangle(J) = HFlowAng + UHangle(J)  ! This is the counter-clockwise angle of the wind

            ENDIF
         ENDDO

      CASE DEFAULT   ! This is how it worked before

         DO I=1,SIZE(getWindSpeedAry)
            IF ( Ht(I) == RefHt ) THEN
               getWindSpeedAry(I) = URef
            ELSEIF ( ABS( Ht(I)-RefHt ) <= 0.5*RotorDiam ) THEN
               getWindSpeedAry(I) = URef*( Ht(I)/RefHt )**PLExp
            ELSEIF ( Ht(I) > 0.0 .AND. RefHt > 0.0 .AND. RefHt /= Z0 ) THEN !Check that we don't have an invalid domain
               getWindSpeedAry(I) = URef*LOG( Ht(I)/Z0 )/LOG( RefHt/Z0 )
            ELSE
               getWindSpeedAry(I) = 0.0
            ENDIF
         ENDDO

   END SELECT

RETURN
END FUNCTION getWindSpeedAry
!=======================================================================
FUNCTION getWindSpeedVal(URef, RefHt, Ht, RotorDiam, profile, UHangle)

   ! Determine the wind speed at a given height, with reference wind speed.
   ! Use power law if given height and reference height are within rotor disk.
   ! Use log profile if reference height is below rotor disk.

   USE                                  TSMods, ONLY: ChebyCoef_WD   ! Chebyshev coefficients (must be defined before calling this function!)
   USE                                  TSMods, ONLY: ChebyCoef_WS   ! Chebyshev coefficients (must be defined before calling this function!)
   USE                                  TSMods, ONLY: HFlowAng       ! The horizontal flow angle of the mean wind speed at hub height
   USE                                  TSMods, ONLY: HubHt          ! The hub height (used with HFlowAng)
   USE                                  TSMods, ONLY: IEC_EWM1       ! IEC Extreme 1-yr wind speed model
   USE                                  TSMods, ONLY: IEC_EWM50      ! IEC Extreme 50-yr wind speed model
   USE                                  TSMods, ONLY: IEC_WindType   ! Type of IEC wind
   USE                                  TSMods, ONLY: NumUSRz        ! Number of user-defined heights
   USE                                  TSMods, ONLY: PLExp          ! Power law exponent
   USE                                  TSMods, ONLY: U_Usr          ! User-defined wind speeds
   USE                                  TSMods, ONLY: Ustar          ! Friction or shear velocity
   USE                                  TSMods, ONLY: Vref           ! IEC Extreme wind value
   USE                                  TSMods, ONLY: WindDir_USR    ! User-defined wind directions
   USE                                  TSMods, ONLY: Z_Usr          ! User-defined heights
   USE                                  TSMods, ONLY: z0             ! Surface roughness length -- It must be > 0 (which we've already checked for)
   USE                                  TSMods, ONLY: ZJetMax        ! Height of the low-level jet
   USE                                  TSMods, ONLY: ZL             ! M-O stability parameter

   IMPLICIT                              NONE

   REAL(ReKi),   INTENT(IN)           :: URef                        ! Wind speed at reference height
   REAL(ReKi),   INTENT(IN)           :: RefHt                       ! Reference height
   REAL(ReKi),   INTENT(IN)           :: Ht                          ! Height where wind speed should be calculated
   REAL(ReKi),   INTENT(IN)           :: RotorDiam                   ! Diameter of rotor disk (meters)
   REAL(ReKi),   INTENT(OUT),OPTIONAL :: UHangle                     ! Horizontal wind angle
   REAL(ReKi)                         :: getWindSpeedVal             ! This function, approximate wind speed at Ht

   REAL(SiKi),   PARAMETER            :: MinZ = 3.                   ! lower bound (height) for Cheby polynomial
   REAL(SiKi),   PARAMETER            :: MaxZ = 500.                 ! upper bound (height) for Cheby polynomial

   CHARACTER(*), INTENT(IN), OPTIONAL :: profile                     ! String that determines what profile is to be used
   CHARACTER(3)                       :: profile_type

   REAL(ReKi)                         :: psiM                        ! The diabatic term for the log wind profile
   REAL(ReKi)                         :: tmp                         ! A temporary variable for calculating psiM
   REAL(ReKi)                         :: tmpHt(2)
   REAL(ReKi)                         :: tmpWS(2)

   INTEGER                            :: I
   INTEGER                            :: Indx

   ! IF Z0 <= 0.0    CALL ProgAbort('The surface roughness must be a positive number')

   IF ( IEC_WindType == IEC_EWM50 ) THEN
      getWindSpeedVal = VRef*( Ht/HubHt )**0.11                      ! [IEC 61400-1 6.3.2.1 (14)]
      RETURN
   ELSEIF ( IEC_WindType == IEC_EWM1 ) THEN
      getWindSpeedVal = 0.8*VRef*( Ht/HubHt )**0.11                  ! [IEC 61400-1 6.3.2.1 (14), (15)]
      RETURN
   ENDIF

   IF ( PRESENT( profile ) ) THEN
      profile_type = profile
      CALL Conv2UC( profile_type )
   ELSE
      profile_type = 'IEC'
   ENDIF

   SELECT CASE ( TRIM(profile_type) )

      CASE ( 'JET', 'J' )

         tmpHt(1) = Ht
         CALL ChebyshevVals( ChebyCoef_WS, tmpHt(1:1), tmpWS(1:1), MinZ, MaxZ ) ! We originally calculated the coeffs for 3-500 m in height
         getWindSpeedVal = tmpWS(1)

         IF ( PRESENT( UHangle ) ) THEN
               ! Calculate the wind direciton at this height
            CALL ChebyshevVals( ChebyCoef_WD, tmpHt(1:1), tmpWS(1:1), MinZ, MaxZ )
            UHangle = tmpWS(1)

               ! Compute the wind direction at hub height & the jet height
            tmpHt(1) = HubHt
            tmpHt(2) = ZJetMax
            CALL ChebyshevVals( ChebyCoef_WD, tmpHt(1:2), tmpWS(1:2), MinZ, MaxZ )

               ! Make sure none of the directions are more than 45 degrees from the direction at the jet height
            IF ( ABS(tmpWS(1) - tmpWS(2) ) > 45. ) THEN  ! The direction at the hub height
               tmpWS(1) = tmpWS(2) + SIGN(REAL(45.,ReKi), tmpWS(1) - tmpWS(2))
            ENDIF

            IF ( ABS(UHangle - tmpWS(2) ) > 45. ) THEN
               UHangle  = tmpWS(2) + SIGN(REAL(45.,ReKi), UHangle  - tmpWS(2))
            ENDIF

               ! Remove the hub height direction so that we have a relative direction, then
               ! add the mean flow angle. (Note that the Chebyshev profile is clockwise looking
               ! from above, but the horizontal angle is counter-clockwise looking from above.)

            UHangle = HFlowAng - (UHangle - tmpWS(1)) ! This is the counter-clockwise angle of the wind

         ENDIF


      CASE ( 'LOG', 'L' ) !Panofsky, H.A.; Dutton, J.A. (1984). Atmospheric Turbulence: Models and Methods for Engineering Applications. New York: Wiley-Interscience; 397 pp.

         IF ( Ht > 0.0 .AND. RefHt > 0.0 .AND. RefHt /= Z0 ) THEN

            IF ( ZL >= 0 ) THEN !& ZL < 1
               psiM = -5.0*ZL
            ELSE
               tmp = (1.0 - 15.0*ZL)**0.25

               !psiM = -2.0*LOG( (1.0 + tmp)/2.0 ) - LOG( (1.0 + tmp*tmp)/2.0 ) + 2.0*ATAN( tmp ) - 0.5 * PI
               psiM = -LOG( 0.125 * ( (1.0 + tmp)**2 * (1.0 + tmp*tmp) ) ) + 2.0*ATAN( tmp ) - 0.5 * PI
            ENDIF

!            IF ( Ustar > 0. ) THEN
!               getWindSpeedVal = ( UstarDiab / 0.4 ) * ( LOG( Ht / Z0 ) - psiM )
!            ELSE
               !In neutral conditions, psiM is 0 and we get the IEC log wind profile:
               getWindSpeedVal = URef*( LOG( Ht / Z0 ) - psiM )/( LOG( RefHt / Z0 ) - psiM )
!            ENDIF

         ELSE
            getWindSpeedVal = 0.0
         ENDIF

      CASE ( 'H2L', 'H' )
               ! Calculate the windspeed.
               ! RefHt and URef both get modified consistently, therefore RefHt is used instead of H_ref.
               !print *, RefHt,H_ref,URef,Ustar ! need to include H_ref from TSmods for this print statement to work.
               getWindSpeedVal = LOG(Ht/RefHt)*Ustar/0.41+URef


      CASE ( 'PL', 'P' )

         IF ( RefHt > 0.0 .AND. Ht > 0.0 ) THEN
            getWindSpeedVal = URef*( Ht/RefHt )**PLExp      ! [IEC 61400-1 6.3.1.2 (10)]
         ELSE
            getWindSpeedVal = 0.0
         ENDIF

      CASE ( 'USR', 'U' )

         IF ( Ht <= Z_USR(1) ) THEN
            getWindSpeedVal = U_USR(1)
         ELSEIF ( Ht >= Z_USR(NumUSRz) ) THEN
            getWindSpeedVal = U_USR(NumUSRz)
         ELSE
            ! Find the two points between which the height lies

            DO I=2,NumUSRz
               IF ( Ht <= Z_USR(I) ) THEN
                  Indx = I-1

                  ! Let's just do a linear interpolation for now
                  getWindSpeedVal = (Ht - Z_USR(Indx)) * ( U_USR(Indx) - U_USR(I) ) / ( Z_USR(Indx) - Z_USR(I) ) + U_USR(Indx)
                  EXIT
               ENDIF
            ENDDO

         ENDIF

         IF ( PRESENT( UHangle ) ) THEN
               ! Calculate the wind direction at this height

            IF ( Ht <= Z_USR(1) ) THEN
               UHangle = WindDir_USR(1)
            ELSEIF ( Ht >= Z_USR(NumUSRz) ) THEN
               UHangle = WindDir_USR(NumUSRz)
            ELSE
               I = Indx + 1

                  ! Let's just do a linear interpolation for now
               !we need to check if the angle goes through 360, before we do the interpolation
               tmpWS(1) = WindDir_USR(Indx) - WindDir_USR(I)
               IF ( tmpWS(1) > 180. ) THEN
                  tmpWS(1) = WindDir_USR(Indx)
                  tmpWS(2) = WindDir_USR(I   ) + 360.
               ELSEIF ( tmpWS(1) < -180. ) THEN
                  tmpWS(1) = WindDir_USR(Indx) + 360.
                  tmpWS(2) = WindDir_USR(I   )
               ELSE
                  tmpWS(1) = WindDir_USR(Indx)
                  tmpWS(2) = WindDir_USR(I   )
               ENDIF

               UHangle = (Ht - Z_USR(Indx)) * ( tmpWS(1) - tmpWS(2) ) / ( Z_USR(Indx) - Z_USR(I) ) + tmpWS(1)

            ENDIF

            UHangle = HFlowAng + UHangle  ! This is the counter-clockwise angle of the wind

         ENDIF

      CASE DEFAULT   ! This is how it worked before

         IF ( Ht == RefHt ) THEN
            getWindSpeedVal = URef
         ELSEIF ( ABS( Ht-RefHt ) <= 0.5*RotorDiam ) THEN
            getWindSpeedVal = URef*( Ht/RefHt )**PLExp                ! [IEC 61400-1 6.3.1.2 (10)]
         ELSEIF ( Ht > 0.0 .AND. RefHt > 0.0 .AND. RefHt /= Z0 ) THEN !Check that we don't have an invalid domain
            getWindSpeedVal = URef*LOG( Ht/Z0 )/LOG( RefHt/Z0 )
         ELSE
            getWindSpeedVal = 0.0
         ENDIF

   END SELECT


RETURN
END FUNCTION getWindSpeedVal
!=======================================================================
FUNCTION getZLARY(WS, Ht)

   USE                                  TSMods, ONLY: RICH_NO        ! Richardson number
   USE                                  TSMods, ONLY: L              ! Rotor-disk averaged L
   USE                                  TSMods, ONLY: profileZmax
   USE                                  TSMods, ONLY: profileZmin
   USE                                  TSMods, ONLY: UstarDiab      ! Diabatic u*0
   USE                                  TSMods, ONLY: WindProfileType
   USE                                  TSMods, ONLY: ZL             ! Rotor-disk averaged z/L
   USE                                  TSMods, ONLY: ZLOffset       ! Offset to align profile with rotor-disk averaged z/L

   IMPLICIT                              NONE

   REAL(ReKi),   INTENT(IN)           :: Ht(:)                       ! Height at which local z/L is defined
   REAL(ReKi),   INTENT(IN)           :: WS(:)                       ! Wind speed(s) at heights, Ht
   
   REAL(ReKi)                         :: tmpZ                        ! a temporary value
   REAL(ReKi)                         :: getZLary(SIZE(Ht))          ! the array of z/L values
   
   INTEGER                            :: IZ
   INTEGER                            :: Zindx
   INTEGER                            :: Zindx_mn (1)
   INTEGER                            :: Zindx_mx (1)
   
   LOGICAL                            :: mask(SIZE(Ht))
   
   mask = Ht.GE.profileZmin
   IF ( ANY(mask) ) THEN
      Zindx_mn = MINLOC( Ht, MASK=mask ) 
            
      mask = Ht.LE.profileZmax
      IF ( ANY(mask) ) THEN
         Zindx_mx = MAXLOC( Ht, MASK=mask )    
            
         DO IZ = 1,SIZE(Ht)
            IF ( Ht(IZ) < profileZmin ) THEN
               Zindx = Zindx_mn(1)
               tmpZ  = Ht(IZ) / Ht(Zindx)    ! This keeps L constant below 50 m
            ELSEIF ( Ht(IZ) > profileZmax ) THEN
               Zindx = Zindx_mx(1)
               tmpZ  = 1.0                   ! L changes above measurement height, but since we don't know how much, we're going to keep z/L constant
            ELSE
               Zindx = IZ
               tmpZ  = 1.0
            ENDIF  !L is constant below 50 meters, and we don't want to extrapolate too high (last measurement is at 116 m)

            IF ( INDEX( 'JU', WindProfileType(1:1) ) > 0 ) THEN
               IF ( Rich_No >= 0 ) THEN
                  getZLary( IZ) =                     - 0.352464*Rich_No + 0.005272*WS(Zindx) + 0.465838
               ELSE
                  getZLary( IZ) =  0.004034*Ht(Zindx) + 0.809494*Rich_No - 0.008298*WS(Zindx) - 0.386632 
               ENDIF !Rich_NO
            ELSE   
               IF ( Rich_No >= 0 ) THEN
                  getZLary( IZ) =  0.003068*Ht(Zindx) + 1.140264*Rich_No + 0.036726*WS(Zindx) - 0.407269 
               ELSE
                  getZLary( IZ) =  0.003010*Ht(Zindx) + 0.942617*Rich_No                     - 0.221886 
               ENDIF   
            ENDIF
            getZLary( IZ) = MIN( getZLary( IZ), 1.0 )
            getZLary( IZ) = getZLary(IZ) * tmpZ
            
         ENDDO         
            
      ELSE ! All are above the max height so instead of extrapolating, we'll use ZL at all heights
         getZLary(:) = ZL
      ENDIF
            
   ELSE ! All are below the min height so we'll keep L constant (as is the case in the surface layer)
      getZLary(:) = Ht(:) / L
   ENDIF
   
   getZLARY = getZLARY(:) + ZLOffset  ! This offset term is used to make the zl profile match the rotor-disk averaged value


END FUNCTION getZLARY
!=======================================================================
SUBROUTINE GP_Turb ( Ht, Ucmp, ZL_tmp, UStar_tmp, Spec )
   ! This subroutine defines the 3-D turbulence spectrum that can be expected over terrain
   ! and heights similiar to the LLLJP project as developed by Neil Kelley & Bonnie Jonkman at NREL.
   ! The use of this subroutine requires that variables have the units of meters and seconds.

USE                     TSMods

IMPLICIT                NONE

   ! Passed variables

REAL(ReKi),INTENT(IN) :: Ht                      ! Height (local)
REAL(ReKi),INTENT(IN) :: Ucmp                    ! Longitudinal Velocity (local)
REAL(ReKi),INTENT(IN) :: Ustar_tmp               ! Local ustar
REAL(ReKi),INTENT(IN) :: ZL_tmp                  ! Local z/l
REAL(ReKi)            :: Spec     (:,:)          ! Working array for PSD

   ! Internal variables

REAL(ReKi), PARAMETER :: Exp53  = 5.0 / 3.0
REAL(ReKi), PARAMETER :: Exp23  = 2.0 / 3.0
REAL(ReKi), PARAMETER :: Exp32  = 3.0 / 2.0
REAL(ReKi)            :: fi       ! Temporary variable for calculation of Spec
REAL(ReKi)            :: fr       ! Temporary variable for calculation of Spec
REAL(ReKi)            :: Freq2    ! Temporary variable for the reduced frequency squared
REAL(ReKi)            :: fr_ih(3) ! Scaling for high-frequency peak location
REAL(ReKi)            :: fr_il(3) ! Scaling for low-frequency peak location
REAL(ReKi)            :: HtZI     ! Temporary variable for calculation of Spec
REAL(ReKi)            :: HtZI2    ! Temporary variable for calculation of Spec
REAL(ReKi)            :: phiE
REAL(ReKi)            :: phiM     ! Non-Dimensional Wind Shear
REAL(ReKi)            :: Pr_ih(3) ! Scaling for magnitude of high-frequency peak
REAL(ReKi)            :: Pr_il(3) ! Scaling for magnitude of low-frequency peak
REAL(ReKi)            :: ps_h
REAL(ReKi)            :: ps_l
REAL(ReKi), PARAMETER :: Scales(2,3) = RESHAPE( (/ 79.0, 13.0, 3.5,    &
                                                  263.0, 32.0, 8.6 /), &
                                         SHAPE=(/2,3/), ORDER=(/2,1/) )
REAL(ReKi)            :: tmpF     ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpFw    ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpX     ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpPhi   ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpZIL   ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpZIU   ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpZU    ! Temporary variable for calculation of Spec
REAL(ReKi)            :: UDen
REAL(ReKi)            :: Ustar_loc ! Local ustar
REAL(ReKi)            :: uStar2   ! Temporary variable holding Ustar-squared
REAL(ReKi)            :: Ustar2F
REAL(ReKi)            :: VDen
REAL(ReKi)            :: X_h
REAL(ReKi)            :: X_l
REAL(ReKi)            :: ZL_loc   ! Local z/l

INTEGER               :: I        ! DO LOOP counter
INTEGER               :: IC       ! DO LOOP counter

uStar2    = Ustar_tmp * Ustar_tmp  ! We don't use ustar_loc here b/c this ustar_loc was used to calculate non-dimensional spectral; this is to scale to dimensional values

ustar_loc = MAX( MIN(ustar_tmp, REAL(1.0,ReKi) ), REAL( 0.15,ReKi) )  ! make sure ustar does not go beyond the observed range that the values were calcualted over
zl_loc    = MAX( MIN(zl_tmp,    REAL(1.0,ReKi) ), REAL(-1.00,ReKi) )  ! make sure z/l does not go beyond the calculated range


IF (zL_loc >= 0) THEN

   phiM   =  1.0 + 4.7*(zL_loc)                   ! = q
   phiE   = (1.0 + 2.5*(zL_loc)**0.6)**Exp32
   
   zl_loc = MAX( zl_loc, REAL(0.025,ReKi) ) !This will prevent 0**-x from becoming infinite.  ustar_loc has this built in already.  This value is the observed min here anyway.

      ! Calculate NEUTRAL/STABLE spectral estimates

   fr_il(1) =  0.014746*(  zl_loc**(-0.37495232))*(ustar_loc**(-0.6167086) )*exp(-0.994591040*zl_loc+1.676298830*ustar_loc)
   fr_ih(1) =  0.043108*(  zl_loc**(-0.39311528))*(ustar_loc**(-2.1719048) )*exp( 0.152732100*zl_loc+2.939119120*ustar_loc)
   Pr_il(1) =  0.003043*(  zl_loc**(-0.60526081))*(ustar_loc**(-2.4348077) )*exp( 1.386013230*zl_loc+2.185372290*ustar_loc)
   Pr_ih(1) = 15.468066*(  zl_loc**( 0.27375765))*(ustar_loc**( 1.8091998) )*exp(-0.266223760*zl_loc-3.091731900*ustar_loc)

   fr_il(2) =  0.0008437*( zl_loc**(-0.79592929))*(ustar_loc**(-1.78297586))*exp( 1.316511335*zl_loc+0.175154746*ustar_loc)
   fr_ih(2) =  1.5278523*( zl_loc**(-0.14197939))*(ustar_loc**( 0.02684469))*exp(-0.261902952*zl_loc-0.672772974*ustar_loc)
   Pr_il(2) =  0.0222952*( zl_loc**( 0.18448738))*(ustar_loc**(-2.23473414))*exp(-1.216594402*zl_loc+1.491864128*ustar_loc)
   Pr_ih(2) =  1.6568440*( zl_loc**(-0.03919916))*(ustar_loc**( 0.57537263))*exp(+0.282805584*zl_loc-1.199845489*ustar_loc)

   fr_il(3) =  1.
   fr_ih(3) =  0.97627403*(zl_loc**(-0.05470045))*(ustar_loc**(0.09666427) )*exp(-0.301255210*zl_loc-0.063122900*ustar_loc)
   Pr_il(3) =  0.
   Pr_ih(3) =  0.69547455*(zl_loc**(-0.00800265))*(ustar_loc**(-0.1352012) )*exp( 0.041784840*zl_loc-0.003785870*ustar_loc)

   fr_il(1) = MAX( MIN( fr_il(1), REAL(0.30,ReKi) ), REAL(0.015,ReKi) )
   fr_ih(1) = MAX( MIN( fr_ih(1), REAL(2.5 ,ReKi) ), REAL(1.25 ,ReKi) )
   Pr_il(1) = MAX( MIN( Pr_il(1), REAL(0.75,ReKi) ), REAL(0.1  ,ReKi) )
   Pr_ih(1) = MAX( MIN( Pr_ih(1), REAL(0.75,ReKi) ), REAL(0.25 ,ReKi) )

   fr_il(2) = MAX( MIN( fr_il(2), REAL(0.3 ,ReKi) ), REAL(0.005,ReKi) )
   fr_ih(2) = MAX( MIN( fr_ih(2), REAL(2.5 ,ReKi) ), REAL(0.75 ,ReKi) )
   Pr_il(2) = MAX( MIN( Pr_il(2), REAL(1.4 ,ReKi) ), REAL(0.05 ,ReKi) )
   Pr_ih(2) = MAX( MIN( Pr_ih(2), REAL(1.0 ,ReKi) ), REAL(0.5  ,ReKi) )

   fr_ih(3) = MAX( MIN( fr_ih(3), REAL(1.4 ,ReKi) ), REAL(0.5  ,ReKi) )
   Pr_ih(3) = MAX( MIN( Pr_ih(3), REAL(1.1 ,ReKi) ), REAL(0.6  ,ReKi) )

   tmpPhi = ( (phiE / phiM)**Exp23 )
   tmpF   = Ht / (Ucmp * phiM)


   DO IC = 1,3  ! Wind components
      DO I = 1,NumFreq
         tmpX  = Freq(I)*tmpF             ! reduced frequency divided by q (q = phiM here)
         X_l   = tmpX/fr_il(ic)
         X_h   = tmpX/fr_ih(ic)

         ps_l  = (Pr_il(ic)*scales(1,ic)*X_l*tmpPhi) / (1.0 + scales(2,ic)*X_l**Exp53);
         ps_h  = (Pr_ih(ic)*scales(1,ic)*X_h*tmpPhi) / (1.0 + scales(2,ic)*X_h**Exp53);

         Spec(I,IC) = (ps_l + ps_h)*uStar2/Freq(I)
      ENDDO
   ENDDO


ELSE
      ! Calculate UNSTABLE spectral estimates
      fr_il(:) = 1.
      fr_ih(:) = 1.
      Pr_il(:) = 1.
      Pr_ih(:) = 1.

! THESE VALUES ARE BASED ON A SMALL AMOUNT OF DATA AND DON'T SEEM TO BEHAVE VERY WELL FOR THE GENERAL CASE.
! Using 1 for each of these values creates the spectral estimates for the SMOOTH model.
!
!      nzl      = -zl_loc
!
!      fr_il(1) = MIN(10.0,  0.0443117*(   (nzl)**(-0.42429))*(ustar_loc**(- 2.03969))*exp( 7.18271*(nzl)+ 1.11017*ustar_loc))
!      fr_ih(1) = MIN( 5.0,  1.10957*(     (nzl)**( 0.18200))*(ustar_loc**(- 0.13968))*exp( 2.48651*(nzl)+ 0.88788*ustar_loc))
!      Pr_il(1) = MIN(20.0,  1.08387e-004*((nzl)**( 0.32784))*(ustar_loc**(- 6.69897))*exp(-8.25590*(nzl)+14.46554*ustar_loc))
!      Pr_ih(1) = MIN( 5.0,  0.0870653*(   (nzl)**(-0.55618))*(ustar_loc**(- 0.85499))*exp( 3.66686*(nzl)- 0.34810*ustar_loc))
!
!      fr_il(2) = MIN( 5.0,  2.8412e-013*( (nzl)**(-0.43587))*(ustar_loc**(-14.62097))*exp( 2.41002*(nzl)+31.59745*ustar_loc))
!      fr_ih(2) = MIN( 5.0,  0.12219003 *( (nzl)**(-0.20010))*(ustar_loc**(- 1.11780))*exp( 1.66314*(nzl)+ 1.74815*ustar_loc))
!      Pr_il(2) = MIN(10.0,  6.6853e-018*( (nzl)**(-1.48280))*(ustar_loc**(-18.80570))*exp( 9.92010*(nzl)+41.12724*ustar_loc))
!      Pr_ih(2) = MIN( 5.0,  2.47627547 *( (nzl)**( 0.04305))*(ustar_loc**(- 0.01287))*exp(-2.74234*(nzl)- 0.95780*ustar_loc))
!
!      fr_il(3) = MIN(30.0,  2.66408e-004*((nzl)**(-0.65260))*(ustar_loc**(- 4.82119))*exp( 7.08116*(nzl)+ 5.85913*ustar_loc))
!      fr_ih(3) = MIN( 5.0,  0.0118916*  ( (nzl)**( 0.09544))*(ustar_loc**(- 2.82943))*exp(-3.21429*(nzl)+ 5.95403*ustar_loc))
!      Pr_il(3) = MIN(10.0,  3.6709e-011*( (nzl)**(-0.96751))*(ustar_loc**(-11.48936))*exp( 5.06644*(nzl)+26.26320*ustar_loc))
!      Pr_ih(3) = MIN( 5.0, 13.53430*    ( (nzl)**(-0.14450))*(ustar_loc**(  1.32560))*exp( 1.66323*(nzl)- 4.28085*ustar_loc))

   tmpZIL = ( ABS(ZI / L) )**Exp23
   HtZI   = Ht / ZI

   tmpZIU = ZI / Ucmp
   tmpZU  = Ht / Ucmp
   HtZI2  = (1.0 - HtZI)**2
   UDen   = 1.0 + 15.0*HtZI
   VDen   = 1.0 +  2.8*HtZI

   DO I=1,NumFreq
      fi      = Freq(I)*tmpZIU
      tmpF    = Freq(I)*tmpZU                ! reduced frequency
      Ustar2F = uStar2/Freq(I)               ! Normalizing term

         ! u component

      fr    = tmpF/UDen
      X_l   = fi/fr_il(1)
      X_h   = fr/fr_ih(1)

      ps_l = (Pr_il(1)* tmpZIL              *  0.50*X_l)/( 1.0 +  2.2*X_l**Exp53 )
      ps_h = (Pr_ih(1)*(HtZI2/(UDen**Exp23))*105.00*X_h)/((1.0 + 33.0*X_h)**Exp53)

      Spec(I,1) = (ps_l + ps_h) * Ustar2F


         ! v component

      fr    = tmpF/VDen
      X_l   = fi/fr_il(2)
      X_h   = fr/fr_ih(2)

      ps_l = (Pr_il(2)* tmpZIL              *  0.95*X_l)/((1.0 +  2.0*X_l)**Exp53)
      ps_h = (Pr_ih(2)*(HtZI2/(VDen**Exp23))* 17.00*X_h)/((1.0 +  9.5*X_h)**Exp53)

      Spec(I,2) = (ps_l + ps_h) * Ustar2F


         ! w component

      Freq2 = tmpF**2
      tmpFw = SQRT( (Freq2 + (0.3*HtZI)**2 ) / (Freq2 + 0.15**2) )
      X_l   = fi  /fr_il(3)
      X_h   = tmpF/fr_ih(3)

      ps_l = tmpFw*(Pr_il(3)*tmpZIL*0.95*X_l)/((1.0 +  2.0*X_l)**Exp53)
      ps_h =       (Pr_ih(3)*HtZI2 *2.00*X_h)/( 1.0 +  5.3*X_h**Exp53 )

      Spec(I,3) = (ps_l + ps_h) * Ustar2F

   ENDDO
ENDIF


RETURN
END SUBROUTINE GP_Turb
!=======================================================================
SUBROUTINE IEC_Kaim ( Ht, Ucmp, Spec )


   ! This subroutine defines the Kaimal PSD model as specified by IEC 61400-1, 2nd Ed. & 3rd Ed.
   ! the use of this subroutine requires that all variables have the units of meters and seconds.

USE                     TSMods

IMPLICIT                NONE

      ! Passed variables

REAL(ReKi),INTENT(IN) :: Ht                      ! Input: Height (Should be HubHt), value ignored
REAL(ReKi),INTENT(IN) :: Ucmp                    ! Input: Velocity (Should be UHub), value ignored
REAL(ReKi)            :: Spec   (:,:)            ! Output: target spectrum

      ! Internal variables

REAL(ReKi),PARAMETER  :: Exp1    = 5.0/3.0
REAL(ReKi)            :: Lambda
REAL(ReKi)            :: LambdaU

REAL(ReKi)            :: LU      (3)
REAL(ReKi)            :: SigmaLU (3)
REAL(ReKi)            :: SigmaU
REAL(ReKi)            :: SigmaV
REAL(ReKi)            :: SigmaW

INTEGER               :: I
INTEGER               :: IVec


   ! Setup scaling values.


SigmaU = SigmaIEC
SigmaV = 0.8*SigmaU
SigmaW = 0.5*SigmaU


   ! Define integral scales.

IF (IECedition == 2) THEN
   IF ( HubHt < 30.0 )  THEN
      Lambda = 0.7*HubHt
   ELSE
      Lambda = 21.0
   ENDIF
ELSE !ED 3
   IF ( HubHt < 60.0 )  THEN
      Lambda = 0.7*HubHt
   ELSE
      Lambda = 42.0
   ENDIF
ENDIF

LambdaU = Lambda/Uhub
LU(1) = 8.10*LambdaU
LU(2) = 2.70*LambdaU
LU(3) = 0.66*LambdaU

   ! Create the spectrum.

SigmaLU(1) = 4.0*SigmaU**2*LU(1)
SigmaLU(2) = 4.0*SigmaV**2*LU(2)
SigmaLU(3) = 4.0*SigmaW**2*LU(3)


DO IVec = 1,3

   LU(IVec) = 6.0*LU(IVec)

   DO I = 1,NumFreq
      Spec(I,IVec) = SigmaLU(IVec) / ( 1.0 + LU(IVec)*Freq(I) )**Exp1
   ENDDO !I

ENDDO !IVec


RETURN
END SUBROUTINE IEC_Kaim
!=======================================================================
SUBROUTINE IEC_vKrm ( Ht, Ucmp, Spec )


   ! This subroutine defines the von Karman PSD model as specified by IEC 61400-1 (2nd Ed).
   ! The use of this subroutine requires that all variables have the units of meters and seconds.


USE                     TSMods

IMPLICIT                NONE

      ! Passed variables

REAL(ReKi),INTENT(IN) :: Ht   ! Must be HubHt, value ignored
REAL(ReKi),INTENT(IN) :: Ucmp ! Must be UHub,  value ignored
REAL(ReKi)            :: Spec     (:,:)

      ! Internal variables

REAL(ReKi),PARAMETER  :: Exp1 =  5.0/6.0
REAL(ReKi),PARAMETER  :: Exp2 = 11.0/6.0
REAL(ReKi)            :: FLU2
REAL(ReKi)            :: L1_U
REAL(ReKi)            :: Lambda
REAL(ReKi)            :: SigmaL1_U
REAL(ReKi)            :: Tmp

INTEGER               :: I


   ! Set up scaling values.

IF ( HubHt  <  30.0 )  THEN
  Lambda = 0.7*HubHt
ELSE
  Lambda = 21.0
ENDIF

   ! Define u-component integral scale.

L1_U   = 3.5*Lambda/UHub
SigmaL1_U = 2.0*SigmaIEC*SigmaIEC*L1_U

DO I=1,NumFreq

   FLU2      = ( Freq(I)*L1_U )**2
   Tmp       = 1.0 + 71.0*FLU2

   Spec(I,1) = 2.0*SigmaL1_U/Tmp**Exp1
   Spec(I,2) = SigmaL1_U*( 1.0 + 189.0*FLU2 )/Tmp**Exp2
   Spec(I,3) = Spec(I,2)

ENDDO ! I

RETURN
END SUBROUTINE IEC_vKrm
!=======================================================================
SUBROUTINE InF_Turb ( Ht, Ucmp, Spec )

   ! This subroutine defines the 3-D turbulence spectrum that can be expected to exist upstream of a large, multi-row
   ! wind park.  It is based on the smooth or homogeneous terrain models of Hojstrup, Olesen, and Larsen of RISO
   ! National Laboratory in Denmark.  The RISO model has been adjusted to reflect the different spectral scaling present
   ! in the flow upwind of a large wind park.  The scaling is based on measurements made by the National Renewable Energy
   ! Laboratory (NREL) in San Gorgonio Pass, California.

USE                     TSMods

IMPLICIT                NONE

   ! Passed variables

REAL(ReKi),INTENT(IN) :: Ht                      ! Height   ( input )
REAL(ReKi),INTENT(IN) :: Ucmp                    ! Velocity ( input )
REAL(ReKi)            :: Spec     (:,:)          ! Working array for PSD ( output )

   ! Internal variables

REAL(ReKi)            :: den                     ! Denominator (replaces Pum_ih, Pum_il, fum_ih, fum_il)
REAL(ReKi), PARAMETER :: Exp1  = 5.0 / 3.0
REAL(ReKi), PARAMETER :: Exp2  = 2.0 / 3.0
REAL(ReKi), PARAMETER :: Exp3  = 3.0 / 2.0
REAL(ReKi)            :: F                       ! Reduced frequency
REAL(ReKi)            :: Fi
REAL(ReKi)            :: Fq                      ! reduced frequency / q
REAL(ReKi)            :: fur_ih
REAL(ReKi)            :: fur_il
REAL(ReKi)            :: fvr_ih
REAL(ReKi)            :: fvr_il
REAL(ReKi)            :: fwr_ih
REAL(ReKi)            :: fwr_il
REAL(ReKi)            :: Fw
REAL(ReKi)            :: HtU                     ! Height / Ucmp
REAL(ReKi)            :: HtZI                    ! Height / ZI     -- used to avoid recalculation
REAL(ReKi)            :: HtZI2                   ! (1.0 - Height / ZI)^2
REAL(ReKi)            :: num                     ! Numerator   (replaces Puo_ih, Puo_il, fuo_ih, fuo_il)
REAL(ReKi)            :: phiE
REAL(ReKi)            :: phiEQ                   ! temp variable
REAL(ReKi)            :: phiM
REAL(ReKi)            :: Ps_h
REAL(ReKi)            :: Ps_l
REAL(ReKi)            :: Pur_ih
REAL(ReKi)            :: Pur_il
REAL(ReKi)            :: Pvr_ih
REAL(ReKi)            :: Pvr_il
REAL(ReKi)            :: Pwr_ih
REAL(ReKi)            :: Pwr_il
REAL(ReKi)            :: q
REAL(ReKi)            :: UDen                    !
REAL(ReKi)            :: UDen2                   !
REAL(ReKi)            :: Ustar2                  ! Ustar**2
REAL(ReKi)            :: Ustar2F                 ! Ustar**2 / Frequency
REAL(ReKi)            :: VDen                    !
REAL(ReKi)            :: VDen2                   !
REAL(ReKi)            :: X                       ! Temporary variable
REAL(ReKi)            :: ZInL                    ! ZI / -L         -- used to avoid recalculation
REAL(ReKi)            :: ZIU                     ! ZI / Ucmp
REAL(ReKi), PARAMETER :: ZI_UVlimit = 1350.0
REAL(ReKi), PARAMETER :: ZI_Wlimit  = 1600.0
REAL(ReKi), PARAMETER :: ZL_MaxObs  =  0.15
REAL(ReKi), PARAMETER :: ZL_MinObs  = -1.00


INTEGER               :: I                       ! Loop counter


Ustar2 = Ustar*Ustar

IF ( ZL < 0 ) THEN
      ! BEGIN UNSTABLE FLOW LOOP

   ! Unstable high-frequency range scaling...

   X = - MAX( ZL, ZL_MinObs)

   Num    = 0.691114  + 0.0791666*X    ! was "Original" Puo_ih =
   Den    = 0.77991   + 0.1761624  / ( 1.0 + EXP( -(X - 0.0405364) / (-0.0184402) ) ) ! was "Measured" Pum_ih =
   Pur_ih = 0.10 * ( Num / Den )
   IF (ZI > ZI_UVlimit) Pur_ih = (ZI / ZI_UVlimit) * Pur_ih

   Num    = 0.421958 * EXP( 0.20739895*X )
   Den    = 0.5247865 + 0.0419204 / ( 1.0 + EXP( -(X - 0.0434172) / (-0.0179269) ) )
   Pvr_ih = Num / Den

   Num    = 0.222875  + 0.1347188*X
   Den    = 0.3542331 + 0.0168806 / ( 1.0 + EXP( -(X - 0.0388899) / (-0.0220998) ) )
   Pwr_ih = 0.80 * ( Num / Den )

   Num    = 0.047465  + 0.0132692*X
   Den    = 0.0599494 - 0.0139033*EXP(-X / 0.02603846)
   fur_ih = 1.75 * ( Num / Den )
   IF (ZI > ZI_UVlimit) fur_ih = (ZI / ZI_UVlimit)*fur_ih

   Num    = 0.18377384 * EXP( 0.2995136*X )
   Den    = 0.1581509  + 0.09501906*X
   fvr_ih = 1.50 * ( Num / Den )
   IF (ZI > ZI_UVlimit) fvr_ih = (ZI / ZI_UVlimit)*fvr_ih

   Num    = 0.3419874 + 0.24985029 * EXP(-X / 0.02619489)
   Den    = 0.451295  + 0.2355227*X
   fwr_ih = 2.0 * ( Num / Den )
   IF (ZI > ZI_Wlimit) fwr_ih = 0.35*(ZI / ZI_Wlimit)*fwr_ih


   ! Unstable low-frequency range scaling...

   Num    = -0.436922  + 2.784789 / ( 1.0 + EXP( -(X - 0.104094) / 0.136708 ) )
   Den    =  0.1392684 + 1.7396251*X
   Pur_il = 2.00 * ( Num / Den )

   Num    = 0.467006  + (5.3032075*X)**1.1713260
   Den    = 0.1425146 + 2.2011562*X
   Pvr_il = 0.25 * ( Num / Den )

   Num    = 0.086908   + (2.3719755 *X)**1.3106297
   Den    = 0.00251981 + (0.50642167*X)**0.6607754
   Pwr_il = Num / Den

   Num    = 0.467962 + 0.9270681*EXP( -X / 0.02039003 )
   Den    = 0.759259 - 0.1448362*X        ! X < 5.24
   fur_il = Num / Den

   Num    = 0.369625 + 1.0772852*EXP( -X / 0.0210098 )
   !Den    = 0.759259 - 0.1448362*X calculated previously
   fvr_il = 2.25 * ( Num / Den )
   IF (ZI > ZI_UVlimit) fvr_il = (ZI / ZI_UVlimit)*fvr_il

   Num    = 3.39482 * EXP( 0.279914*X )
   Den    = 4.59769 + 12.58881*EXP( -X / 0.03351852 )
   fwr_il = 2.25 * ( Num / Den )
   IF (ZI > ZI_Wlimit) fwr_il=4.0*(ZI / ZI_Wlimit)*fwr_il

   HtZI  = Ht / ZI
   HtZI2 = (1.0 - HtZI)**2
   ZInL  = ( ZI / ( -L ) )**Exp2
   HtU   = Ht / Ucmp
   ZIU   = ZI / Ucmp
   UDen  = 1.0 + 15.0*HtZI
   VDen  = 1.0 +  2.8*HtZI
   UDen2 = HtZI2 / UDen**Exp2
   VDen2 = HtZI2 / VDen**Exp2


   DO I = 1,NumFreq

         F   = Freq(I) * HtU
         Fi  = Freq(I) * ZIU

         ! Bonnie: These () around 0.3 HtZI are incorrect as compared to the original SMOOTH model. (For now, leave as is since parameters were-supposedly-calculated with this formulation)
         Fw  = SQRT( (F**2 + (0.3*HtZI**2) ) / (F**2 + 0.0225) )

         Ustar2F = Ustar2 / Freq(I)

      ! CALCULATE UNSTABLE LONGITUDINAL SPECTRAL COMPONENT, nSu(n)/(u*)^2, then multiply by (u*)^2/n

         X    = Fi / fur_il
         Ps_l = Pur_il * ( (0.5*X) / (1.0 + 2.2 * X**Exp1) ) * ZInL

         X     = F / (UDen * fur_ih)     !  Fru = F / UDen
         Ps_h = ( (105.0 * X) / (1.0 + 33.0 * X )**Exp1 ) * UDen2
         Ps_h = Ps_h*Pur_ih

         Spec(I,1) = ( Ps_l + Ps_h ) * Ustar2F

      ! CALCULATE UNSTABLE CROSSWIND SPECTRAL COMPONENT, nSv(n)/(u*)^2, then multiply by (u*)^2/n

         X    = Fi / fvr_il
         Ps_l = ( (0.95*X) / (1.0 + 2.0 * X)**Exp1 ) * ZInL
         Ps_l = Ps_l*Pvr_il

         X    = F / (VDen * fvr_ih)      ! Frv = F / VDen
         Ps_h = ( (17.0 * X) / (1.0 + 9.5*X)**Exp1 ) * VDen2
         Ps_h = Ps_h*Pvr_ih

         Spec(I,2) = ( Ps_l + Ps_h ) * Ustar2F

      ! CALCULATE UNSTABLE VERTICAL SPECTRAL COMPONENT, nSw(n)/(u*)^2, then multiply by (u*)^2/n

         X    = Fi / fwr_il
         Ps_l = Fw * ( (0.95*X) / (1.0 + 2.0*X)**Exp1 ) * ZInL
         Ps_l = Ps_l*Pwr_il

         X    = F / fwr_ih
         Ps_h = ( (2.0*X) / (1.0 + 5.3 * X**Exp1) ) * HtZI2
         Ps_h = Ps_h * Pwr_ih

         Spec(I,3) = ( Ps_l + Ps_h ) * Ustar2F

      ENDDO

ELSE ! ZL >= 0    ! BEGIN STABLE FLOW LOOP

      X = MIN(ZL, ZL_MaxObs)

   ! Get stable spectral peaks

   ! Calculate smooth terrain scaling functions

      phiE = (1.0 + 2.5 * X**0.6) **Exp3
      phiM = 1.0 + 4.7*X

      q = phiM

   ! Stable high-frequency (shear) range scaling ...

      Num    = 0.8029768 + ( 1.708247*X )**3.669245
      Den    = 1.5431 * EXP( 1.6379*X )
      Pur_ih = 0.01*( Num / Den )

      Num    = 0.419234  + ( 2.759119*X )**1.4483715
      Den    = 0.89717 * EXP( 1.67034*X )
      Pvr_ih = 1.30 * ( Num / Den )

      Num    = 0.239692  + ( 2.3531204*X )**1.062937
      Den    = 0.5324 * EXP( 1.6314*X )
      Pwr_ih = Num / Den
      Pwr_ih = Pwr_ih - 2.0*X
      IF ( Pwr_ih <= 0.0) Pwr_ih = 1.0
      Pwr_ih = 1.5*Pwr_ih

      Num    = 0.042393 + ( 1.28175*X )**1.409066
      Den    = 0.045 + 0.21137*X
      fur_ih = 3.5 * ( Num / Den )

      Num    = 0.220831 + (0.630632*X)**0.8120686
      Den    = 0.160 + 0.74876*X
      fvr_ih = 1.25 * ( Num / Den )

      Num    = 0.382558 + (1.3640485*X)**1.524565
      Den    = 0.350 + 1.638806*X
      fwr_ih = 1.5 * ( Num / Den )

   ! Low-frequency range scaling...

      Num    = 0.88418 + (11.665367*X)**0.794753
      Den    = 1.55288 * EXP( 1.56925*X )
      Pur_il = 1.50 * ( Num / Den )

      Num    = 0.4671733 + 4.3093084 * X**(0.859202)
      Den    = 0.90382 * EXP( 1.59076*X )
      Pvr_il = 0.75 * ( Num / Den )

      Num    = 0.076136 + 2.644456 * X**(1.207014)
      Den    = 0.533202 * EXP( 1.51415*X )
      Pwr_il = Num / Den
      Pwr_il = Pwr_il - 1.75*X

      Num    = 0.009709 + ( 0.4266236*X )**1.644925
      Den    = 0.045 + 0.212038*X
      fur_il = 2.00 * ( Num / Den )
      fur_il = ABS(fur_il - 3.0*X)

      Num    = 0.0220509 + ( 0.93256713*X )**1.719292
      Den    = 0.160 + 0.74985*X
      fvr_il = 1.15 * ( Num / Den )

      Num    = 0.0351474 + ( 1.4410838*X )**1.833043
      Den    = 0.350 + 1.645667*X
      fwr_il = Num / Den


      phiEQ = (phiE / q)**Exp2

      DO I = 1,NumFreq

         ! CALCULATE Reduced Frequency, f

         f  = Freq(I) * Ht / Ucmp
         fq = f / q     ! was XU = f/qu, XV = f/qv, XW = f/qw

         Ustar2F = Ustar2 / Freq(I)

         ! CALCULATE NEUTRAL/STABLE LONGITUDINAL SPECTRAL COMPONENT, nSu(n)/(u*)^2, then multiply by (u*)^2/n

         X    = fq / fur_ih
         Ps_h = ( (79.0 * X) / (1.0 + 263.0 * X**Exp1) ) * phiEQ
         Ps_h = Ps_h*Pur_ih

         X    = fq / fur_il
         Ps_l = ( (79.0 * X) / (1.0 + 263.0 * X**Exp1) ) * phiEQ
         Ps_l = Ps_l*Pur_il

         Spec(I,1) = ( Ps_l + Ps_h ) * Ustar2F

         ! CALCULATE NEUTRAL/STABLE CROSSWIND SPECTRAL COMPONENT, nSv(n)/(u*)^2, then multiply by (u*)^2/n

         X    = fq / fvr_ih
         Ps_h = ( (13.0 * X) / (1.0 + 32.0 * X**Exp1) ) * phiEQ
         Ps_h = Ps_h*Pvr_ih

         X    = fq / fvr_il
         Ps_l = ( (13.0 * X) / (1.0 + 32.0 * X**Exp1) ) * phiEQ
         Ps_l = Ps_l*Pvr_il

         Spec(I,2) = ( Ps_h + Ps_l ) * Ustar2F

         ! CALCULATE NEUTRAL/STABLE VERTICAL SPECTRAL COMPONENT, nSw(n)/(u*)^2, then multiply by (u*)^2/n

         X    = fq / fwr_ih
         Ps_h = ( ( 3.5 * X) / (1.0 +  8.6 * X**Exp1) ) * phiEQ
         Ps_h = Ps_h*Pwr_ih

         X    = fq / fwr_il
         Ps_l = ( ( 3.5 * X) / (1.0 +  8.6 * X**Exp1) ) * phiEQ
         Ps_l = Ps_l*Pwr_il

         Spec(I,3) = ( Ps_l + Ps_h ) * Ustar2F

   ENDDO ! I

ENDIF  ! ZL < 0


RETURN
END SUBROUTINE InF_Turb
!=======================================================================
SUBROUTINE NWTC_Turb ( Ht, Ucmp, Spec )

   ! This subroutine defines the 3-D turbulence spectrum that can be expected
   ! over terrain and heights similiar to the NWTC LIST project as developed
   ! by Neil Kelley & Bonnie Jonkman at NREL. The use of this subroutine
   ! requires that variables have the units of meters and seconds.

USE                     TSMods

IMPLICIT                NONE

   ! Passed variables

REAL(ReKi),INTENT(IN) :: Ht                      ! Height (local)
REAL(ReKi),INTENT(IN) :: Ucmp                    ! Longitudinal Velocity (local)
REAL(ReKi)            :: Spec     (:,:)          ! Working array for PSD

   ! Internal variables

REAL(ReKi), PARAMETER :: Exp53  = 5.0 / 3.0
REAL(ReKi), PARAMETER :: Exp23  = 2.0 / 3.0
REAL(ReKi), PARAMETER :: Exp32  = 3.0 / 2.0
REAL(ReKi)            :: fi       ! Temporary variable for calculation of Spec
REAL(ReKi)            :: fr       ! Temporary variable for calculation of Spec
REAL(ReKi)            :: Freq2    ! Temporary variable for the reduced frequency squared
REAL(ReKi)            :: fr_ih(3) ! Scaling for high-frequency peak location
REAL(ReKi)            :: fr_il(3) ! Scaling for low-frequency peak location
REAL(ReKi)            :: HtZI     ! Temporary variable for calculation of Spec
REAL(ReKi)            :: HtZI2    ! Temporary variable for calculation of Spec
REAL(ReKi)            :: phiE
REAL(ReKi)            :: phiM     ! Non-Dimensional Wind Shear
REAL(ReKi)            :: Pr_ih(3) ! Scaling for magnitude of high-frequency peak
REAL(ReKi)            :: Pr_il(3) ! Scaling for magnitude of low-frequency peak
REAL(ReKi)            :: ps_h
REAL(ReKi)            :: ps_l
REAL(ReKi), PARAMETER :: Scales(2,3) = RESHAPE( (/ 79.0, 13.0, 3.5,    &
                                                  263.0, 32.0, 8.6 /), &
                                         SHAPE=(/2,3/), ORDER=(/2,1/) )
REAL(ReKi)            :: tmpF     ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpFw    ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpX     ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpPhi   ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpZIL   ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpZIU   ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpZU    ! Temporary variable for calculation of Spec
REAL(ReKi)            :: UDen
REAL(ReKi)            :: Ustar_tmp ! Disk-averaged ustar, limited by the observed range of values for fitting these emperical functions
REAL(ReKi)            :: uStar2   ! Temporary variable holding Ustar-squared
REAL(ReKi)            :: Ustar2F
REAL(ReKi)            :: VDen
REAL(ReKi)            :: X_h
REAL(ReKi)            :: X_l
REAL(ReKi)            :: ZL_tmp   ! Disk-averaged z/l, limited by the observed range of z/l for fitting these emperical functions

INTEGER               :: I        ! DO LOOP counter
INTEGER               :: IC       ! DO LOOP counter

uStar2 = Ustar * Ustar

IF (zL >= 0) THEN

   zl_tmp = max( min(zl, REAL(3.5,ReKi) ), REAL(0.005,ReKi) )

      ! Calculate NEUTRAL/STABLE spectral estimates

   fr_il(1)  = 0.096376774*(zl_tmp**(-0.315715361)) * exp(-0.385026736*zl_tmp)
   fr_ih(1)  = 1.690996304*(zl_tmp**(-0.340366943)) * exp(-0.132661086*zl_tmp)
   Pr_il(1)  = 1.209487882*(zl_tmp**( 0.052273494)) * exp( 0.189014328*zl_tmp)
   Pr_ih(1)  = 0.224103219*(zl_tmp**( 0.169561956)) * exp( 0.222723480*zl_tmp)

   fr_il(2)  = 0.032285308*(zl_tmp**(-0.387804427)) * exp(-0.388660410*zl_tmp)
   fr_ih(2)  = 0.473438689*(zl_tmp**(-0.441450751)) * exp( 0.290697895*zl_tmp)
   Pr_il(2)  = 1.285421087*(zl_tmp**( 0.006644801)) * exp( 0.354496483*zl_tmp)
   Pr_ih(2)  = 0.991251080*(zl_tmp**( 0.343831230)) * exp(-0.605373943*zl_tmp)

   fr_il(3)  = 0.097156827*(zl_tmp**(-0.096412942)) * exp(-0.616256651*zl_tmp)
   fr_ih(3)  = 0.469904415*(zl_tmp**(-0.218253779)) * exp(-0.157526974*zl_tmp)
   Pr_il(3)  = 0.368138932*(zl_tmp**( 0.093776256)) * exp( 0.109020969*zl_tmp)
   Pr_ih(3)  = 0.638868926*(zl_tmp**( 0.035396647)) * exp(-0.031884105*zl_tmp)

   fr_il(1)  = MAX( MIN( fr_il(1),REAL( 0.40,ReKi) ), REAL(0.015,ReKi) )
   fr_ih(1)  = MAX( MIN( fr_ih(1),REAL(10.0 ,ReKi) ), REAL(0.35 ,ReKi) )
   Pr_il(1)  = MAX( MIN( Pr_il(1),REAL( 2.25,ReKi) ), REAL(0.8  ,ReKi) )
   Pr_ih(1)  = MAX( MIN( Pr_ih(1),REAL( 0.8 ,ReKi) ), REAL(0.05 ,ReKi) )

   fr_il(2)  = MAX( MIN( fr_il(2),REAL( 0.23,ReKi) ), REAL(0.003,ReKi) )
   fr_ih(2)  = MAX( MIN( fr_ih(2),REAL( 3.0 ,ReKi) ), REAL(0.25 ,ReKi) )
   Pr_il(2)  = MAX( MIN( Pr_il(2),REAL( 2.25,ReKi) ), REAL(0.95 ,ReKi) )
   Pr_ih(2)  = MAX( MIN( Pr_ih(2),REAL( 1.0 ,ReKi) ), REAL(0.2  ,ReKi) )

   fr_il(3)  = MAX( MIN( fr_il(3),REAL( 0.175,ReKi)), REAL(0.006,ReKi) )
   fr_ih(3)  = MAX( MIN( fr_ih(3),REAL( 1.25 ,ReKi)), REAL(0.2  ,ReKi) )
   Pr_il(3)  = MAX( MIN( Pr_il(3),REAL( 0.75 ,ReKi)), REAL(0.2  ,ReKi) )
   Pr_ih(3)  = MAX( MIN( Pr_ih(3),REAL( 1.0  ,ReKi)), REAL(0.25 ,ReKi) )

   phiM   =  1.0 + 4.7*(zl_tmp)                   ! = q
   phiE   = (1.0 + 2.5*(zl_tmp)**0.6)**Exp32

   tmpPhi = ( (phiE / phiM)**Exp23 )
   tmpF   = Ht / (Ucmp * phiM)


   DO IC = 1,3  ! Wind components
      DO I = 1,NumFreq
         tmpX  = Freq(I)*tmpF             ! reduced frequency divided by q (q = phiM here)
         X_l   = tmpX/fr_il(ic)
         X_h   = tmpX/fr_ih(ic)

         ps_l  = (Pr_il(ic)*scales(1,ic)*X_l*tmpPhi) / (1.0 + scales(2,ic)*X_l**Exp53);
         ps_h  = (Pr_ih(ic)*scales(1,ic)*X_h*tmpPhi) / (1.0 + scales(2,ic)*X_h**Exp53);

         Spec(I,IC) = (ps_l + ps_h)*uStar2/Freq(I)
      ENDDO
   ENDDO

!print *, fr_il
!print *, fr_ih
!print *, Pr_il
!print *, Pr_ih
!print *, zl_tmp, Ustar

ELSE
      ! Calculate UNSTABLE spectral estimates

   zl_tmp    = abs( min( max( zl  ,REAL(-0.5,ReKi) ),REAL( -0.025,ReKi) ) )
   ustar_tmp =      max( min(ustar,REAL( 1.4,ReKi) ),REAL(  0.2  ,ReKi) )

   fr_il(1)  =   0.08825035*(zl_tmp**(-0.08806865))*(ustar_tmp**(-0.26295052))*exp( 1.74135233*zl_tmp + 1.86785832*ustar_tmp)
   fr_ih(1)  =   1.34307411*(zl_tmp**(-0.55126969))*(ustar_tmp**(-0.07034031))*exp( 0.40185202*zl_tmp - 0.55083463*ustar_tmp)
   Pr_il(1)  =  57.51578485*(zl_tmp**(-1.89080610))*(ustar_tmp**( 4.03260796))*exp( 6.09158000*zl_tmp - 7.47414385*ustar_tmp)
   Pr_ih(1)  =   4.52702491*(zl_tmp**( 0.72447070))*(ustar_tmp**(-0.10602646))*exp(-3.73265876*zl_tmp - 0.49429015*ustar_tmp)

   fr_il(2)  =   0.58374913*(zl_tmp**(-0.53220033))*(ustar_tmp**( 1.49509302))*exp( 3.61867635*zl_tmp - 0.98540722*ustar_tmp)
   fr_ih(2)  =   4.30596626*(zl_tmp**( 0.31302745))*(ustar_tmp**(-0.26457011))*exp(-1.41513284*zl_tmp + 0.91503248*ustar_tmp)
   Pr_il(2)  =  32.06436225*(zl_tmp**(-1.43676866))*(ustar_tmp**( 3.57797045))*exp( 5.31617813*zl_tmp - 5.76800891*ustar_tmp)
   Pr_ih(2)  =   3.93109762*(zl_tmp**( 0.57974534))*(ustar_tmp**(-0.20510478))*exp(-4.85367443*zl_tmp - 0.61610914*ustar_tmp)

   fr_il(3)  =   0.81092087*(zl_tmp**(-0.03483105))*(ustar_tmp**( 0.58332966))*exp(-0.10731274*zl_tmp - 0.16463702*ustar_tmp)
   fr_ih(3)  =   1.05515450*(zl_tmp**(-0.25002535))*(ustar_tmp**( 0.14528047))*exp( 1.00641958*zl_tmp - 0.67569359*ustar_tmp)
   Pr_il(3)  =   6.60003543*(zl_tmp**(-0.45005503))*(ustar_tmp**( 1.35937877))*exp( 2.45632937*zl_tmp - 1.98267575*ustar_tmp)
   Pr_ih(3)  =  16.56290180*(zl_tmp**( 0.40464339))*(ustar_tmp**( 0.82276250))*exp(-3.92300971*zl_tmp - 1.82957067*ustar_tmp)


   fr_il(1)  = MAX( MIN( fr_il(1), REAL(1.50,ReKi) ), REAL(0.2 ,ReKi)  )
   fr_ih(1)  = MAX( MIN( fr_ih(1), REAL(8.0 ,ReKi) ), REAL(0.1 ,ReKi)  )
   Pr_il(1)  = MAX( MIN( Pr_il(1), REAL(8.0 ,ReKi) ), REAL(1.0 ,ReKi)  )
   Pr_ih(1)  = MAX( MIN( Pr_ih(1), REAL(1.2 ,ReKi) ), REAL(0.1 ,ReKi)  )

   fr_il(2)  = MAX( MIN( fr_il(2), REAL(2.3 ,ReKi) ), REAL(0.12,ReKi)  )
   fr_ih(2)  = MAX( MIN( fr_ih(2), REAL(7.5 ,ReKi) ), REAL(1.8 ,ReKi)  )
   Pr_il(2)  = MAX( MIN( Pr_il(2), REAL(8.0 ,ReKi) ), REAL(0.2 ,ReKi)  )
   Pr_ih(2)  = MAX( MIN( Pr_ih(2), REAL(0.9 ,ReKi) ), REAL(0.2 ,ReKi)  )

   fr_il(3)  = MAX( MIN( fr_il(3), REAL(1.4 ,ReKi) ), REAL(0.2 ,ReKi)  )
   fr_ih(3)  = MAX( MIN( fr_ih(3), REAL(1.75,ReKi) ), REAL(0.95,ReKi)  )
   Pr_il(3)  = MAX( MIN( Pr_il(3), REAL(7.0 ,ReKi) ), REAL(1.0 ,ReKi)  )
   Pr_ih(3)  = MAX( MIN( Pr_ih(3), REAL(1.0 ,ReKi) ), REAL(0.3 ,ReKi)  )

   tmpZIL = (-ZI / L)**Exp23
   HtZI   = Ht / ZI

   tmpZIU = ZI / Ucmp
   tmpZU  = Ht / Ucmp
   HtZI2  = (1.0 - HtZI)**2
   UDen   = 1.0 + 15.0*HtZI
   VDen   = 1.0 +  2.8*HtZI

   DO I=1,NumFreq
      fi      = Freq(I)*tmpZIU

      tmpF    = Freq(I)*tmpZU                ! reduced frequency
      Ustar2F = uStar2/Freq(I)               ! Normalizing term

         ! u component

      fr    = tmpF/UDen
      X_l   = fi/fr_il(1)
      X_h   = fr/fr_ih(1)

      ps_l = (Pr_il(1)* tmpZIL              *  0.50*X_l)/( 1.0 +  2.2*X_l**Exp53 )
      ps_h = (Pr_ih(1)*(HtZI2/(UDen**Exp23))*105.00*X_h)/((1.0 + 33.0*X_h)**Exp53)

      Spec(I,1) = (ps_l + ps_h) * Ustar2F


         ! v component

      fr    = tmpF/VDen
      X_l   = fi/fr_il(2)
      X_h   = fr/fr_ih(2)

      ps_l = (Pr_il(2)* tmpZIL              *  0.95*X_l)/((1.0 +  2.0*X_l)**Exp53)
      ps_h = (Pr_ih(2)*(HtZI2/(VDen**Exp23))* 17.00*X_h)/((1.0 +  9.5*X_h)**Exp53)

      Spec(I,2) = (ps_l + ps_h) * Ustar2F


         ! w component

      Freq2 = tmpF**2
      tmpFw = SQRT( (Freq2 + (0.3*HtZI)**2 ) / (Freq2 + 0.15**2) )
      X_l   = fi  /fr_il(3)
      X_h   = tmpF/fr_ih(3)

      ps_l = tmpFw*(Pr_il(3)*tmpZIL*0.95*X_l)/((1.0 +  2.0*X_l)**Exp53)
      ps_h =       (Pr_ih(3)*HtZI2 *2.00*X_h)/( 1.0 +  5.3*X_h**Exp53 )

      Spec(I,3) = (ps_l + ps_h) * Ustar2F

   ENDDO
ENDIF


RETURN

END SUBROUTINE NWTC_Turb
!=======================================================================
SUBROUTINE OutF_Turb ( Ht, Ucmp, Spec )

USE                     TSMods

IMPLICIT                NONE

   !Passed variables

REAL(ReKi),INTENT(IN) :: Ht                      ! Height
REAL(ReKi),INTENT(IN) :: Ucmp                    ! Velocity
REAL(ReKi)            :: Spec     (:,:)          ! Working array for PSD

   ! Internal variables

REAL(ReKi)            :: A0
REAL(ReKi)            :: A1
REAL(ReKi)            :: A2
REAL(ReKi)            :: A3
REAL(ReKi), PARAMETER :: Exp1  = 5.0 / 3.0
REAL(ReKi), PARAMETER :: Exp2  = 2.0 / 3.0
REAL(ReKi), PARAMETER :: Exp3  = 3.0 / 2.0
REAL(ReKi)            :: den                     ! Denominator (replaces Pum_oh, Pum_ol, fum_oh, fum_ol, Pvm_oh)
REAL(ReKi)            :: F                       ! Reduced frequency
REAL(ReKi)            :: Fi
REAL(ReKi)            :: fur_oh
REAL(ReKi)            :: fur_ol
REAL(ReKi)            :: fvr_oh
REAL(ReKi)            :: fvr_ol
REAL(ReKi)            :: fvr_wk
REAL(ReKi)            :: Fw
REAL(ReKi)            :: fwr_oh
REAL(ReKi)            :: fwr_ol
REAL(ReKi)            :: Fq                      ! reduced frequency / q
REAL(ReKi)            :: HtZI                    ! Height / ZI     -- used to avoid recalculation
REAL(ReKi)            :: HtZI2                   ! (1.0 - Height / ZI)^2
REAL(ReKi)            :: num                     ! Numerator   (replaces Puo_oh, Puo_ol, fuo_oh, fuo_ol, Pvo_wk)
REAL(ReKi)            :: phiE
REAL(ReKi)            :: phiM
REAL(ReKi)            :: Ps_h
REAL(ReKi)            :: Ps_l
REAL(ReKi)            :: Ps_wk
REAL(ReKi)            :: Pur_oh                  ! High Frequency Range
REAL(ReKi)            :: Pur_ol                  ! Low Frequency Range
REAL(ReKi)            :: Pvr_oh
REAL(ReKi)            :: Pvr_ol
REAL(ReKi)            :: Pvr_wk
REAL(ReKi)            :: Pwr_oh
REAL(ReKi)            :: Pwr_ol
REAL(ReKi)            :: q
REAL(ReKi)            :: tmp                     ! holds calculation common to several formulae
REAL(ReKi)            :: UDen                    !
REAL(ReKi)            :: UDen2                   !
REAL(ReKi)            :: Ustar2                  ! Ustar**2
REAL(ReKi)            :: Ustar2F                 ! Ustar**2 / Frequency
REAL(ReKi)            :: VDen                    !
REAL(ReKi)            :: VDen2                   !
REAL(ReKi)            :: X                       ! Temporary variable
REAL(ReKi)            :: ZInL                    ! ZI / -L         -- used to avoid recalculation
REAL(ReKi), PARAMETER :: ZL_MaxObs =  0.4        ! The largest z/L value where the spectral peak scaling should work.
REAL(ReKi), PARAMETER :: ZL_MinObs = -1.0        ! The smallest z/L value where the spectral peak scaling should work.

INTEGER               :: I                       ! Loop counter



Ustar2 = Ustar*Ustar

IF (ZL < 0) THEN


      ! Get Unstable spectral peaks

   ! Unstable high-frequency range scaling...

   X = - MAX( ZL, ZL_MinObs )

   Num    = 0.598894  + 0.282106  * EXP(-X / 0.0594047)
   Den    = 0.600977  + 9.137681  / (1.0 + EXP( -(X + 0.830756) / (-0.252026) ))
   Pur_oh = 0.1 * (Num / Den)

   Num    = 0.4830249 + 0.3703596 * EXP(-X / 0.0553952)
   Den    = 0.464604  + 1.900294  / (1.0 + EXP( -(X + 0.928719) / (-0.317242) ))
   Pvr_oh = 5.0 * (Num / Den)

   Num    = 0.320112  + 0.229540  * EXP(-X / 0.0126555)
   Den    = 0.331887  + 1.933535  / (1.0 + EXP( -(X + 1.19018 ) / (-0.3011064) ))
   Pwr_oh = 1.25 * (Num / Den)

   Num    =  0.049279 + EXP(0.245214 * X * 2.478923) ! was Num    = 0.049279  + EXP(0.245214 * X)**2.478923
   Den    = -2.333556 + 2.4111804 / (1.0 + EXP( -(X + 0.623439) / 0.1438076))
   fur_oh =  0.3 * (Num / Den)

   Num    = -2.94362   + 3.155970 / (1.0 + EXP( -(X + 0.872698) / 0.245246))
   Den    =  0.0171463 + 0.188081 / (1.0 + EXP( -(X + 0.711851) / 0.688910))
   fvr_oh =  2.0 * (Num / Den)

   Num    = 0.7697576 * EXP( -X / 3.8408779 ) - 0.561527 * EXP( -X / 0.1684403 ) ! was Num = Beta4(X,A0,A1,A2,A3)
   Den    = 0.512356  - 0.044946  / (1.0 + EXP( -(X - 0.066061) / (-0.0121168) ))
   fwr_oh = 1.75 * (Num / Den)
   IF (ZI < 1350.0 ) fwr_oh = (ZI / 1350.0) * fwr_oh

      ! Unstable low-frequency range scaling ...

   Num    = 0.796264 + 0.316895 / (1.0 + EXP( -(X - 0.082483) / 0.027480 ))
   Den    = 0.07616  + EXP(0.303919 * X * 0.390906)   ! was Den = 0.07616 + EXP(0.303919*X)**0.390906
   Pur_ol = 4.0 * (Num / Den)
   IF (ZI < 1600.0) Pur_ol = (ZI / 1600.0) * Pur_ol

   Num    = 0.812483 + 0.1332134 * X
   Den    = 0.104132 + EXP(0.714674 * X * 0.495370)   ! was Den = 0.104132 + EXP(0.714674*X)**0.495370
   Pvr_ol = Num / Den
   Pvr_ol = (ZI / 1600.0)*Pvr_ol

   Num    = 0.371298  + 0.0425447 * X
   Den    = 0.0004375 + EXP(0.4145751 * X * 0.6091557)   ! was Den = 0.0004375 + EXP(0.4145751*X)**0.6091557
   Pwr_ol = 0.75 * (Num / Den)

   Num    = 0.859809 * EXP(0.157999 * X)
   Den    = 0.81459 + 0.021942 * X
   fur_ol = 1.5 * (Num / Den)
   IF (ZI > 1850.0) fur_ol = 2.6 * (ZI / 1850.0) * fur_ol

   !A0 =  0.8121775
   !A1 =  4.122E+15
   !A2 = -0.594909
   !A3 =  0.0559581
   Num    = 0.8121775 * EXP( -X / 4.122E+15 ) - 0.594909 * EXP( -X / 0.0559581 ) ! was Num = BETA4(X,A0,A1,A2,A3)
   Den    = 0.72535  - 0.0256291 * X
   fvr_ol = 3.0 * (Num / Den)
   fvr_ol = (ZI / 1600.0) * fvr_ol

   Num    = 6.05669  * EXP(-0.97418 * X)
   Den    = 3.418386 + 9.58012 / (1.0 + EXP( -(X - 0.0480283) / (-0.022657) ))
   fwr_ol = 0.9 * (Num / Den)

      ! Unstable Wake Range Scaling for v-component only

   Num    = 0.247754 + 0.16703142 * EXP(-X / 0.1172513)
   Den    = 0.464604 + 1.900294 / (1.0 + EXP( -(X + 0.928719) / (-0.317242) ))
   Pvr_wk = 0.05 * (Num / Den)

   !A0 = 0.72435
   !A1 = 0.0436448
   !A2 = 0.08527
   Num    = 0.72435 / (1.0 + EXP( -(X - 0.0436448) / 0.08527 ))    ! was Num = BETA5(X,A0,A1,A2)
   Den    = 0.0171463 + 0.188081 / (1.0 + EXP( -(X + 0.711851) / 0.688910))
   fvr_wk = 3.0 * (Num / Den)

   HtZI  = Ht / ZI
   HtZI2 = (1.0 - HtZI)**2
   ZInL  = ( ZI / (-L) )**Exp2
   UDen  = 1.0 + 15.0 * HtZI
   VDen  = 1.0 +  2.8 * HtZI
   UDen2 = HtZI2 / UDen**Exp2
   VDen2 = HtZI2 / VDen**Exp2

   DO I = 1,NumFreq

      ! Calculate f,fi,fru,frv

      F   = Freq(I)*Ht / Ucmp
      Fi  = Freq(I)*ZI / Ucmp
      Fw  = SQRT( (F**2 + (0.3*HtZI**2)) / (F**2 + 0.0225) )

      Ustar2F = Ustar2 / Freq(I)

         ! CALCULATE UNSTABLE LONGITUDINAL SPECTRAL COMPONENT, nSu(n)/(u*)^2, then multiply by (u*)^2/n

         ! No identifiable wake contribution was found in u-component

      X    = Fi / fur_ol
      Ps_l = ( (0.5*X) / (1.0 + 2.2 * X**Exp1) ) * ZInL
      Ps_l = ABS(Ps_l*Pur_ol)

      X    = F / (UDen * fur_oh)
      Ps_h = ( (105.0 * X) / (1.0 + 33.0 * X)**Exp1 ) * UDen2
      Ps_h = Ps_h * Pur_oh
      Spec(I,1) = ( Ps_l + Ps_h ) * Ustar2F

         ! CALCULATE UNSTABLE CROSSWIND SPECTRAL COMPONENT, nSv(n)/(u*)^2, then multiply by (u*)^2/n

      X    = Fi / fvr_ol
      Ps_l = ( (0.95*X) / (1.0 + 2.0*X)**Exp1 ) * ZInL
!     Ps_l = ABS(Psv_l)*Pvr_ol

      X    = F / (VDen * fvr_oh)
      Ps_h = ( (17.0 * X) / (1.0 + 9.5*X)**Exp1 ) * VDen2
!     Ps_h = Ps_h*Pvr_oh

         ! Wake contribution for v-component only
      X     = F / (VDen * fvr_wk)
      Ps_wk = ( (17.0 * X) / (1.0 + 9.5 * X)**Exp1 ) * VDen2
      Ps_wk = Ps_wk*Pvr_wk
      Spec(I,2) = ( Ps_l + Ps_h + Ps_wk ) * Ustar2F

         ! CALCULATE UNSTABLE VERTICAL SPECTRAL COMPONENT, nSw(n)/(u*)^2, then multiply by (u*)^2/n

         ! No identifiable wake contribution was found in w-component

      X    = Fi / fwr_ol
      Ps_l = Fw*( (0.95 * X) / (1.0 + 2.0 * X)**Exp1 ) * ZInL
      Ps_l = ABS(Ps_l)*Pwr_ol

      X    = F / fwr_oh
      Ps_h = ( (2.0 * X) / (1.0 + 5.3 * X**Exp1) ) * HtZI2
      Ps_h = Ps_h*Pwr_oh

      Spec(I,3) = ( Ps_l + Ps_h ) * Ustar2F

   ENDDO

ELSE  ! ZL >= 0

         ! BEGIN STABLE FLOW LOOP...

   ! Get Stable spectral peaks

   ! Stable high-frequency (wake) range scaling...

   X      = MIN( ZL, ZL_MaxObs )

   Num    = 0.149471 + 0.028528 * &
            EXP( -EXP( -( (X - 0.003580) / 0.0018863 ) ) - ( (X - 0.0035802) / 0.0018863) + 1.0)
   Den    = 1.563166  * EXP(1.137965 * X)
   Pur_oh = 0.35 * (Num / Den)

   A0     = 2.66666062
   A1     = 0.0034082
   A2     = 0.0229827
   Num    = Beta3(X, A0, A1, A2)
   A0     = 33.942268
   A1     =  0.0160732
   A2     = -0.008654
   A3     =  0.0053586
   Num    = Num + Beta1(X, A0, A1, A2, A3)
   Num    = 1.0 / Num
   Den    = 0.9170783 * EXP(1.152393 * X)
   Pvr_oh = 2.25 * (Num / Den)

   A0     = 0.1990569
   A1     = 0.0286048
   A2     = 0.006751
   Num    = Beta3(X, A0, A1, A2)
   A0     = 0.0435354
   A1     = 0.0599214
   A2     = 0.0520877
   Num    = Num + Beta2(X, A0, A1, A2)
   Den    = 0.539112  * EXP(1.124104 * X)
   Pwr_oh = 0.9  * (Num / Den)

   tmp    = -(X - 0.037003738) / 0.01612278
   Num    = 0.764910145 + 0.654370025 * EXP( -EXP(tmp) + tmp + 1.0 )
   Den    = 0.045 + 0.209305  * X
   fur_oh = Num / Den

   A0     = 0.5491507
   A1     = 0.0099211
   A2     = 0.0044011
   Num    = Beta3(X, A0, A1, A2)
   A0     = 0.0244484
   A1     = 0.0139515
   A2     = 0.0109543
   Num    = Num + Beta2(X, A0, A1, A2)
   Den    = 0.160 + 0.7496606 * X
   fvr_oh = 0.5 * (Num / Den)

   Num    = 0.391962642 + 0.546722344*EXP( -0.5* ( (X - 0.023188588) / 0.018447575)**2 )
   Den    = 0.350 + 1.6431833 * X
   fwr_oh = 2.0 * (Num / Den)

      ! Stable low-frequency range scaling ...

   Num    = 0.894383 + (1.55915 * X)**3.111778
   Den    = 1.563317 * EXP(1.137965 * X)
   Pur_ol = 0.9 * (Num / Den)

   Num    = 0.747514 + (1.57011 * X)**1.681581
   Den    = 0.910783 * EXP(1.1523931 * X)
   Pvr_ol = 0.60 * (Num / Den - 1.75 * X)

   Num    = 0.376008 * EXP(1.4807733* X)
   Den    = 0.539112 * EXP(1.124104 * X)
   Pwr_ol = 0.6 * (Num / Den - 2.0 * X)

   Num    = 0.023450 + (0.3088194 * X)**1.24710
   Den    = 0.045    +  0.209305  * X
   fur_ol = 1.5 * (Num / Den - X)
   fur_ol = MAX( fur_ol, REAL( 0.1,ReKi ) )   ! We divide by this number so it should not get too small.

   Num    = 0.051616 + (0.8950263 * X)**1.37514
   Den    = 0.160    +  0.749661  * X
   fvr_ol = 0.5 * (Num / Den)

   Num    = 0.250375 - 0.690491  * X + 2.4329342 * X**2
   Den    = 0.350    + 1.6431833 * X
   fwr_ol = Num / Den

      ! Calculate smooth terrain scaling functions

   phiE   = (1.0 + 2.5 * X**0.6)**Exp3
   phiM   =  1.0 + 4.7 * X
   q      = phiM

   tmp    = (phiE / q)**Exp2

   DO I = 1,NumFreq

      F       = Freq(I) * Ht / Ucmp    ! Reduced frequency

      Fq      = F / q
      Ustar2F = Ustar2 / Freq(I)

         ! CALCULATE NEUTRAL/STABLE LONGITUDINAL SPECTRAL COMPONENT, nSu(n)/(u*)^2, then multiply by (u*)^2/n

      X    = Fq / fur_ol
      Ps_l = ( (79.0 * X) / (1.0 + 263.0 * X**Exp1) ) * tmp
      Ps_l = ABS(Ps_l * Pur_ol)

      X    = Fq / fur_oh
      Ps_h = ( (79.0 * X) / (1.0 + 263.0 * X**Exp1) ) * tmp
      Ps_h = Ps_h * Pur_oh

      Spec(I,1) = (Ps_l + Ps_h) * Ustar2F

         ! CALCULATE NEUTRAL/STABLE CROSSWIND SPECTRAL COMPONENT, nSv(n)/(u*)^2, then multiply by (u*)^2/n

      X    = Fq / fvr_ol
      Ps_l = ( (13.0 * X) / (1.0 + 32.0 * X**Exp1) ) * tmp
      Ps_l = ABS(Ps_l * Pvr_ol)

      X    = Fq / fvr_oh
      Ps_h = ( (13.0 * X) / (1.0 + 32.0 * X**Exp1) ) * tmp
      Ps_h = Ps_h * Pvr_oh

      Spec(I,2) = (Ps_h + Ps_l) * Ustar2F

         ! CALCULATE NEUTRAL/STABLE VERTICAL SPECTRAL COMPONENT, nSw(n)/(u*)^2, then multiply by (u*)^2/n

      X    = Fq / fwr_ol
      Ps_l = ( (3.5 * X) / (1.0 + 8.6 * X**Exp1) ) * tmp
      Ps_l = ABS(Ps_l * Pwr_ol)

      X    = Fq / fwr_oh
      Ps_h = ( (3.5 * X) / (1.0 + 8.6 * X**Exp1) ) * tmp
      Ps_h = Ps_h * Pwr_oh

      Spec(I,3) = (Ps_l + Ps_h) * Ustar2F

   ENDDO

ENDIF    ! ZL < 0

RETURN
END SUBROUTINE OutF_Turb
!=======================================================================
FUNCTION pkCTKE_LLJ(Ht)

USE              TSMods, ONLY: ZL
USE              TSMods, ONLY: UStar

IMPLICIT                 NONE

REAL(ReKi)               :: A                                        ! A constant/offset term in the pkCTKE calculation   
REAL(ReKi)               :: A_uSt                                    ! The scaling term for Ustar 
REAL(ReKi)               :: A_zL                                     ! The scaling term for z/L
REAL(ReKi), INTENT(IN)   :: Ht                                       ! The height at the billow center
REAL(ReKi)               :: pkCTKE_LLJ                               ! The max CTKE expected for LLJ coh structures
REAL(ReKi)               :: rndCTKE                                  ! The random residual
REAL(ReKi), PARAMETER    :: RndParms(5) = (/0.252510525, -0.67391279, 2.374794977, 1.920555797, -0.93417558/) ! parameters for the Pearson IV random residual
REAL(ReKi)               :: z                                        ! The height
REAL(ReKi), PARAMETER    :: z_Ary(4)    = (/54., 67., 85., 116./)    ! Aneomoeter heights

INTEGER                  :: Zindx_mn (1)


   Zindx_mn = MINLOC( ABS(z_Ary-Ht) )
   
   SELECT CASE ( Zindx_mn(1) )
      CASE ( 1 )  ! 54 m
         A     = -0.051   
         A_zL  = -0.0384  
         A_uSt =  9.9710
         
      CASE ( 2 )  ! 67 m
         A     = -0.054
         A_zL  = -0.1330  
         A_uSt = 10.2460
         
      CASE ( 3 )  ! 85 m
         A     = -0.062   
         A_zL  = -0.1320  
         A_uSt = 10.1660
         
      CASE ( 4 )  !116 m
         A     = -0.092   
         A_zL  = -0.3330  
         A_uSt = 10.7640

      CASE DEFAULT !This should not occur
         CALL TS_Abort( 'Error in pkCTKE_LLJ():: Height index is invalid.' )
   END SELECT   
    
   CALL RndPearsonIV( rndCTKE, RndParms, (/REAL(-10.,ReKi), REAL(17.5,ReKi)/) )

   pkCTKE_LLJ = MAX(0.0, A + A_uSt*UStar + A_zL*ZL + rndCTKE)
   
END FUNCTION pkCTKE_LLJ
!=======================================================================
FUNCTION PowerLawExp( Ri_No )

   ! This function calculates the power law exponent for the wind turbulence models
   ! WF_UPW, WF_07D, and WF_14D

USE                      TSMods

IMPLICIT                 NONE

REAL(ReKi), INTENT(IN) :: Ri_No                ! Richardson Number
REAL(ReKi)             :: PowerLawExp          ! Power Law exponent for particular model


IF ( KHtest ) THEN
   PowerLawExp = 0.3
   RETURN
ENDIF

SELECT CASE ( TRIM(TurbModel) )

   CASE ('WF_UPW', 'NWTCUP')
      IF ( Ri_No > 0.0 ) THEN
         PowerLawExp = 0.14733
      ELSE
         PowerLawExp = 0.087687698 + 0.059641545*EXP(Ri_No/0.04717783)
      ENDIF

   CASE ( 'WF_07D', 'WF_14D' )
      IF ( Ri_No > 0.04 ) THEN
         PowerLawExp = 0.17903
      ELSE
         PowerLawExp = 0.127704032 + 0.031228952*EXP(Ri_No/0.0805173)
      ENDIF

   CASE ('SMOOTH', 'GP_LLJ', 'TIDAL', 'RIVER')
      ! A 1/7 power law seems to work ok for HYDRO spectral models also...
      PowerLawExp = 0.143

   CASE DEFAULT
      IF ( IEC_WindType == IEC_EWM1 .OR. IEC_WindType == IEC_EWM50 ) THEN
         PowerLawExp = 0.11         ! [IEC 61400-1 6.3.2.1 (14)]
      ELSEIF ( IECstandard == 3 ) THEN
         PowerLawExp = 0.14         ! [IEC 61400-3 Page 22 (3)]
      ELSE
         PowerLawExp = 0.2          ! [IEC 61400-1 6.3.1.2 (10)]
      ENDIF

END SELECT

RETURN
END FUNCTION PowerLawExp
!=======================================================================
SUBROUTINE PSDcal ( Ht, Ucmp, ZL_loc, Ustar_loc )

   ! This routine calls the appropriate spectral model.

USE                                TSMods

IMPLICIT                            NONE


REAL(ReKi), INTENT(IN)           :: Ht             ! Height
REAL(ReKi), INTENT(IN)           :: Ucmp           ! Velocity
REAL(ReKi), INTENT(IN), OPTIONAL :: Ustar_loc      ! Local ustar
REAL(ReKi), INTENT(IN), OPTIONAL :: ZL_loc         ! Local Z/L


SELECT CASE ( TRIM(TurbModel) )
   CASE ( 'GP_LLJ' )
      IF ( PRESENT(Ustar_loc) .AND. PRESENT(ZL_loc) ) THEN
         CALL GP_Turb   ( Ht, Ucmp, ZL_loc, Ustar_loc, Work )
      ELSE
         CALL GP_Turb   ( Ht, Ucmp, ZL,     Ustar,     Work )
      ENDIF
   CASE ( 'IECKAI' )
!bjj TEST: CALL TestSpec ( Ht, Ucmp, Work )
      CALL IEC_Kaim  ( Ht, Ucmp, Work )
   CASE ( 'IECVKM' )
      CALL IEC_vKrm  ( Ht, Ucmp, Work )
   CASE ('NWTCUP')
      CALL NWTC_Turb  ( Ht, Ucmp, Work )
   CASE ( 'SMOOTH' )
      CALL St_Turb   ( Ht, Ucmp, Work )
   CASE ( 'TIDAL', 'RIVER' )
      CALL Tidal_Turb  ( Ucmp, Work )
   CASE ( 'USRINP' )
      CALL UsrSpec   ( Ht, Ucmp, Work )
   CASE ( 'USRVKM' )
      CALL vonKrmn   ( Ht, Ucmp, Work )
   CASE ('WF_UPW')
      CALL InF_Turb  ( Ht, Ucmp, Work )
   CASE ( 'WF_07D', 'WF_14D' )
      CALL OutF_Turb ( Ht, Ucmp, Work )
   CASE ( 'NONE  ' )
      Work(:,:) = 0.0

   CASE ( 'MODVKM' )
      IF (MVK) THEN
         CALL Mod_vKrm( Ht, Ucmp, Work )
      ELSE
         CALL TS_Abort ( 'Specified turbulence PSD, "'//TRIM( TurbModel )//'", not availible.' )
      ENDIF

   CASE DEFAULT
      CALL TS_Abort ( 'Specified turbulence PSD, "'//TRIM( TurbModel )//'", not availible.' )
END SELECT

IF ( PSD_OUT ) THEN
   !IF ( ABS(Ht - HubHt) < Tolerance ) THEN
      WRITE( UP, FormStr ) 1, Ht, Work(:,1)
      WRITE( UP, FormStr ) 2, Ht, Work(:,2)
      WRITE( UP, FormStr ) 3, Ht, Work(:,3)
   !ENDIF
ENDIF


RETURN
END SUBROUTINE PSDcal
!=======================================================================
SUBROUTINE ReadEventFile( Un, ScaleWid, ScaleVel, CTKE )

      ! This subroutine reads the events definitions from the event data file

   USE                     TSMods


   IMPLICIT                NONE


      ! Passed Variables

INTEGER,    INTENT(IN)  :: Un             ! I/O Unit
REAL(ReKi),INTENT(IN)   :: CTKE           ! Predicted maximum CTKE
REAL(ReKi),INTENT(INOUT):: ScaleVel       ! The shear we're scaling for
REAL(ReKi),INTENT(IN)   :: ScaleWid       ! The height of the wave we're scaling with

      ! Local variables
REAL(ReKi)              :: MaxEvtCTKE        ! The maximum CTKE in the dataset of events

INTEGER                 :: AllocStat      ! Array allocation status
INTEGER                 :: I              ! DO loop counter
INTEGER                 :: IOS            ! I/O Status


   MaxEvtCTKE = 0.0  ! initialize the MAX variable


         ! Read the nondimensional lateral width of the dataset, Ym_max

   CALL ReadRVar( Un, CTEventFile, Ym_max, "the nondimensional lateral width of the dataset", &
                                                   "Nondimensional lateral dataset width")

         ! Read the nondimensional vertical height of the dataset, Zm_max

   CALL ReadRVar( Un, CTEventFile, Zm_max, "the nondimensional vertical height of the dataset", &
                                          "Nondimensional vertical dataset height")


         ! Read the rest of the header

   CALL ReadIVar( Un, CTEventFile, NumEvents, "NumEvents", "the number of coherent structures.")


   IF ( NumEvents > 0 ) THEN

            ! Allocate memory for coherent event start times and lengths

      ALLOCATE ( EventName(NumEvents) , STAT=AllocStat )

      IF ( AllocStat /= 0 )  THEN
         CALL TS_Abort ( 'Error allocating memory for the coherent event name.' )
      ENDIF

      ALLOCATE ( EventTS(NumEvents) , STAT=AllocStat )

      IF ( AllocStat /= 0 )  THEN
         CALL TS_Abort ( 'Error allocating memory for the coherent event timestep lengths.' )
      ENDIF

      ALLOCATE ( EventLen(NumEvents) , STAT=AllocStat )

      IF ( AllocStat /= 0 )  THEN
         CALL TS_Abort ( 'Error allocating memory for the coherent event lengths.' )
      ENDIF

      ALLOCATE ( pkCTKE(NumEvents) , STAT=AllocStat )

      IF ( AllocStat /= 0 )  THEN
         CALL TS_Abort ( 'Error allocating memory for the coherent event peak CTKE.' )
      ENDIF


            ! Read the last header lines

      CALL ReadCom( Un, CTEventFile, 'the fourth header line')  ! A blank line
      CALL ReadCom( Un, CTEventFile, 'the fifth header line')   ! The column heading lines


            ! Read the event definitions and scale times by TScale

      DO I=1,NumEvents

         READ ( Un, *, IOSTAT=IOS )  EventName(I),  EventTS(I), EventLen(I), pkCTKE(I)

         IF ( IOS /= 0 )  THEN
            CALL TS_Abort ( 'Error reading event '//TRIM( Int2LStr( I ) )//' from the coherent event data file.' )
         ENDIF
         MaxEvtCTKE = MAX( MaxEvtCTKE, pkCTKE(I) )

      ENDDO
      
      IF ( MaxEvtCTKE > 0.0 ) THEN          
            ScaleVel = MAX( ScaleVel, SQRT( CTKE / MaxEvtCTKE ) )
            ! Calculate the Velocity Scale Factor, based on the requested maximum CTKE
      ENDIF
         
         ! Calculate the TimeScaleFactor, based on the Zm_max in the Events file.
         
      TSclFact = ScaleWid / (ScaleVel * Zm_max)
      
         ! Scale the time based on TSclFact
         
      DO I=1,NumEvents
         EventLen(I) = EventLen(I)*TSclFact
      ENDDO
      
   ELSE
   
      TSclFact = ScaleWid / (ScaleVel * Zm_max) 

   ENDIF  ! FileNum > 0


END SUBROUTINE ReadEventFile
!=======================================================================
SUBROUTINE ReadCVarDefault ( UnIn, Fil, CharVar, VarName, VarDescr, Def, IGNORE )


      ! This routine reads a single character variable from the next line of the input file.
      ! The input is allowed to be "default"

      ! Argument declarations:


   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.

   LOGICAL, INTENT(INOUT)       :: Def                                             ! - on input whether or not to use the default - on output, whether a default was used
   LOGICAL, INTENT(IN), OPTIONAL:: IGNORE                                          ! whether to ignore this input

   CHARACTER(250)               :: CharLine                                        ! Character string being read.
   CHARACTER(*), INTENT(INOUT)  :: CharVar                                         ! Character variable being read.
   CHARACTER( *), INTENT(IN)    :: Fil                                             ! Name of the input file.
   CHARACTER( *), INTENT(IN)    :: VarDescr                                        ! Text string describing the variable.
   CHARACTER( *), INTENT(IN)    :: VarName                                         ! Text string containing the variable name.

      ! Local declarations:

!  CHARACTER(38)                :: Frmt = "( 2X, ES11.4e2, 2X, A, T30, ' - ', A )" ! Output format for real parameters
   CHARACTER(38)                :: Frmt = "( A10, 2X, A )"                         ! Output format for real parameters


   CALL ReadCVar( UnIn, Fil, CharLine, VarName, VarDescr )

   IF ( PRESENT(IGNORE) ) THEN
      IF ( IGNORE ) THEN
         WRITE (UnEc,Frmt)  "N/A", VarDescr
         Def = .TRUE.
         RETURN
      ENDIF

   ENDIF

   CALL Conv2UC( CharLine )

   IF ( TRIM(CharLine) == 'DEFAULT' ) THEN

      CALL WrScr ( '    A default value will be used for '//TRIM(VarName)//'.' )

      IF ( Def ) THEN  ! use the value as a default
         WRITE (UnEc,Frmt)  CharVar, VarDescr
      ELSE
         WRITE (UnEc,Frmt)  "CALCULATED", VarDescr
      ENDIF

      Def = .TRUE.

   ELSE

      CharVar = CharLine

      WRITE (UnEc,Frmt)  CharVar, VarDescr

      Def = .FALSE.

   ENDIF

   RETURN
END SUBROUTINE ReadCVarDefault ! ( UnIn, Fil, RealVar, VarName, VarDescr )
!=======================================================================
SUBROUTINE ReadLIVar( UnIn, Fil, IntVar, VarName, VarDescr )

      ! This routine reads a single integer variable from the next line of the input file.
      ! BJJ modified from NWTC_subs ReadIVar(), with just the call to ReadNum() removed.  This 
      !     will allow users to keep their T/F values in their input files if they want to.

      ! Argument declarations:

   INTEGER, INTENT(OUT)         :: IntVar                                          ! Integer variable being read.
   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.

   CHARACTER(*), INTENT(IN)     :: Fil                                             ! Name of the input file.
   CHARACTER(*), INTENT(IN)     :: VarDescr                                        ! Text string describing the variable.
   CHARACTER(*), INTENT(IN)     :: VarName                                         ! Text string containing the variable name.


      ! Local declarations:

   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

   CHARACTER(33)                :: Frmt = "( 2X, I11, 2X, A, T30, ' - ', A )"      ! Output format for integer parameters.


   READ (UnIn,*,IOSTAT=IOS)  IntVar

   CALL CheckIOS ( IOS, Fil, VarName, NumType )

   IF ( Echo )  THEN
      WRITE (UnEc,Frmt)  IntVar, VarName, VarDescr
   END IF


   RETURN


END SUBROUTINE ReadLIVar
!=======================================================================
SUBROUTINE ReadRAryDefault ( UnIn, Fil, RealAry, VarName, VarDescr, Def, IGNORE )

      ! This routine reads a real array from the next line of the input file.
      ! The input is allowed to be "default"

      ! Argument declarations:

   REAL(ReKi), INTENT(INOUT)    :: RealAry (:)                                     ! Real variable being read.

   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.

   LOGICAL, INTENT(INOUT)       :: Def                                             ! - on input whether or not to use the default - on output, whether a default was used
   LOGICAL, INTENT(IN), OPTIONAL:: IGNORE                                          ! whether or not to ignore this input

   CHARACTER(250)               :: CharLine                                        ! Character string being read.
   CHARACTER( *), INTENT(IN)    :: Fil                                             ! Name of the input file.
   CHARACTER( *), INTENT(IN)    :: VarDescr                                        ! Text string describing the variable.
   CHARACTER( *), INTENT(IN)    :: VarName                                         ! Text string containing the variable name.

      ! Local declarations:

   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

!  CHARACTER(38)                :: Frmt = "( 2X, ES11.4e2, 2X, A, T30, ' - ', A )" ! Output format for real parameters
   CHARACTER(38)                :: Frmt = "( F10.3, 2X, , A )"                     ! Output format for real parameters

   WRITE(Frmt,"(I2)") SIZE(RealAry)-1

   Frmt = "( '(',F9.3,"//TRIM(Frmt)//"(',',G10.3),')',2X , A )"

   CALL ReadCVar( UnIn, Fil, CharLine, VarName, VarDescr ) !Maybe I should read this in explicitly...

   IF ( PRESENT(IGNORE) ) THEN
      IF ( IGNORE ) THEN
         WRITE (UnEc,"( A10, 2X, A )")  "N/A", VarDescr
         Def = .TRUE.
         RETURN
      ENDIF

   ENDIF

   CALL Conv2UC( CharLine )

   IF ( TRIM(CharLine) == 'DEFAULT' ) THEN

      CALL WrScr ( '    A default value will be used for '//TRIM(VarName)//'.' )

      IF ( Def ) THEN  ! use the value as a default
         WRITE (UnEc,Frmt)  RealAry, VarDescr
      ELSE
         WRITE (UnEc,"( A10, 2X, A )")  "CALCULATED", VarDescr
      ENDIF

      Def = .TRUE.

   ELSE

      IF ( INDEX( CharLine(1:1), 'TF') > 0 ) THEN    ! We don't want 'T' or 'F' read as -1 or 0.
         CALL WrScr1 ( ' Invalid numerical input for "'//TRIM( VarName )//'".' )
      ENDIF

      READ (CharLine,*,IOSTAT=IOS)  RealAry

      IF (IOS /=0) THEN
         READ (CharLine,*,IOSTAT=IOS)  RealAry(1)  ! Try reading only the first element
      ENDIF

      CALL CheckIOS ( IOS, Fil, VarName, NumType )

      WRITE (UnEc,Frmt)  RealAry, VarDescr

      Def = .FALSE.

   ENDIF


   RETURN

END SUBROUTINE ReadRAryDefault
!=======================================================================
SUBROUTINE ReadRVarDefault ( UnIn, Fil, RealVar, VarName, VarDescr, Def, IGNORE, IGNORESTR )

      ! This routine reads a single real variable from the next line of the input file.
      ! The input is allowed to be "default"

      ! Argument declarations:

   REAL(ReKi), INTENT(INOUT)    :: RealVar                                         ! Real variable being read.

   INTEGER, INTENT(IN)          :: UnIn                                            ! I/O unit for input file.

   LOGICAL, INTENT(INOUT)       :: Def                                             ! - on input whether or not to use the default - on output, whether a default was used
   LOGICAL, INTENT(IN), OPTIONAL:: IGNORE                                          ! whether or not to ignore this input
   LOGICAL, INTENT(OUT),OPTIONAL:: IGNORESTR                                       ! whether or not user requested to ignore this input

   CHARACTER(250)               :: CharLine                                        ! Character string being read.
   CHARACTER( *), INTENT(IN)    :: Fil                                             ! Name of the input file.
   CHARACTER( *), INTENT(IN)    :: VarDescr                                        ! Text string describing the variable.
   CHARACTER( *), INTENT(IN)    :: VarName                                         ! Text string containing the variable name.

      ! Local declarations:

   INTEGER                      :: IOS                                             ! I/O status returned from the read statement.

!  CHARACTER(38)                :: Frmt = "( 2X, ES11.4e2, 2X, A, T30, ' - ', A )" ! Output format for real parameters
   CHARACTER(38)                :: Frmt = "( F10.3, 2X, A )"                       ! Output format for real parameters


   CALL ReadCVar( UnIn, Fil, CharLine, VarName, VarDescr )

   IF ( PRESENT(IGNORE) ) THEN

      IF ( IGNORE ) THEN
         WRITE (UnEc,"( A10, 2X, A )")  "N/A", VarDescr
         Def = .TRUE.
         RETURN
      ENDIF

   ENDIF


   CALL Conv2UC( CharLine )

   IF ( PRESENT(IGNORESTR) ) THEN
      IF ( TRIM( CharLine ) == 'NONE' ) THEN
         IGNORESTR = .TRUE.
         WRITE (UnEc,"( A10, 2X, A )")  "N/A", VarDescr
         Def = .TRUE.
         RealVar = 0.0  ! This is set for the Reynolds stress inputs, but if IGNORESTR is used for other inputs, it may need to be changed
         RETURN
      ENDIF
   ENDIF

   IF ( TRIM(CharLine) == 'DEFAULT' ) THEN

      CALL WrScr ( '    A default value will be used for '//TRIM(VarName)//'.' )

      IF ( PRESENT(IGNORESTR) ) THEN
         IF ( IGNORESTR ) THEN  !We've told it to ignore this input, by default
            WRITE (UnEc,"( A10, 2X, A )")  "N/A", VarDescr
            RealVar = 0.0  ! This is set for the Reynolds stress inputs, but if IGNORESTR is used for other inputs, it may need to be changed
            RETURN
         ENDIF
      ENDIF

      IF ( Def ) THEN  ! use the value as a default
         WRITE (UnEc,Frmt)  RealVar, VarDescr
      ELSE
         WRITE (UnEc,"( A10, 2X, A )")  "CALCULATED", VarDescr
      ENDIF

      Def = .TRUE.

   ELSE

      IF ( INDEX( CharLine(1:1), 'TF') > 0 ) THEN    ! We don't want 'T' or 'F' read as -1 or 0.
         CALL WrScr1 ( ' Invalid numerical input for "'//TRIM( VarName )//'".' )
      ENDIF

      READ (CharLine,*,IOSTAT=IOS)  RealVar

      CALL CheckIOS ( IOS, Fil, VarName, NumType )

      WRITE (UnEc,Frmt)  RealVar, VarDescr

      Def = .FALSE.
      
      IF ( PRESENT(IGNORESTR) ) THEN
         IGNORESTR = .FALSE.
      ENDIF        

   ENDIF


   RETURN
END SUBROUTINE ReadRVarDefault
!=======================================================================
SUBROUTINE Rnd4ParmLogNorm( RandNum, Parms, FnRange )
   ! This subroutine generates a random variate with a PDF of the form
   ! f(x) = A + B*exp(-0.5*(ln(x/C)/D)^2)
   ! a truely log-normal distribution has A = 0

IMPLICIT                            NONE

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

   CALL RndUnif( RN )        ! Generate RN from U(0,1)

   CALL RndUnif( RandNum )   ! Generate RandNum from h(y) = 1 / (MaxVALUE - MinVALUE)
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
SUBROUTINE Rnd3ParmNorm( RandNum, A, B, C, xMin, xMax )

! Calculates a deviate from a distribution with pdf:
! f(x) = A * EXP( -0.5 * ((x-B)/C)**2 )
! A 3-parameter normal distribution
! We assume the returned values are between -1 and 1, since this is for the Cross-Correlations, unless
! the optional values, xMin and xMax are used

IMPLICIT                            NONE

REAL(ReKi), INTENT(IN)           :: A
REAL(ReKi), INTENT(IN)           :: B
REAL(ReKi), INTENT(IN)           :: C

REAL(ReKi)                       :: fMAX                ! Max(f(x)) = f(B) = A
REAL(ReKi)                       :: Gx                  ! The function g(x) = f(x)/fMAX
REAL(ReKi)                       :: MaxVALUE            ! Maximum value of returned variate
REAL(ReKi)                       :: MinVALUE            ! Minimum value of returned variate
REAL(ReKi), INTENT(OUT)          :: RandNum             ! numbers distributed with the pdf above
REAL(ReKi)                       :: RN                  ! A random number for the acceptance-rejection method
REAL(ReKi), INTENT(IN), OPTIONAL :: xMax                ! The maximum returned iterate
REAL(ReKi), INTENT(IN), OPTIONAL :: xMin                ! The minimum returned iterate

INTEGER                          :: Count
INTEGER, PARAMETER               :: MaxIter = 10000     ! Max number of iterations to converge (so we don't get an infinite loop)


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
   CALL TS_Abort('Parameter A must at least 1/(xMax-xMin) and parameter C cannot be zero in this 3-parameter normal distribution.')
ENDIF


fMAX  = A
Count = 1

      ! Generate a 3-parameter normal distribution on (-1,1) from a uniform distribution

DO WHILE (Count < MaxIter)

   CALL RndUnif( RN )        ! Generate RN from U(0,1)
   CALL RndUnif( RandNum )   ! Generate RandNum from h(y) = 1 / (MaxVALUE - MinVALUE)
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
SUBROUTINE RndExp( RandExpNum, mu )

   ! This subroutine computes an exponential distribution on (0,inf). If the
   ! number of random variates to return is large, a different algorithm will
   ! probably be faster (i.e. one that computes LOG(x) fewer times).
   ! mu must be positive, it defaults to 1.0 if the parameter is not included.
   ! RandNum has p.d.f.(x) = 1/mu * exp(-x/mu), x>=0
   ! The expected value of RandNum is mu.


   IMPLICIT                NONE

      ! Passed Variables

REAL(ReKi), INTENT(OUT)          :: RandExpNum             ! The exponentially distributed numbers in (0,1)
REAL(ReKi), INTENT(IN), OPTIONAL :: mu                     ! The exponential distribution parameter equal to the expected value of the distribution

      ! Local Variable

REAL(ReKi)                       :: mu_use


IF ( PRESENT(mu) ) THEN
   IF (mu < 0.0) THEN
      CALL TS_Abort ( 'Invalid mu parameter in exponential distribution.' )
   ENDIF
   mu_use = mu
ELSE
   mu_use = 1.0
ENDIF


   ! Get a uniform distribution of random numbers

CALL RndUnif( RandExpNum )

IF ( RandExpNum == 0.0 )  THEN ! We shouldn't get two zeros in a row...
   CALL RndUnif( RandExpNum )
ENDIF

   ! Transform the uniform distribution to an exponential distribution

RandExpNum = - mu_use * LOG( RandExpNum )

END SUBROUTINE RndExp
!=======================================================================
SUBROUTINE RndInit

   ! Initialize the Random Number Generators

USE                                 Ran_Lux_Mod
USE                                 TSMods

IMPLICIT                            NONE

REAL(ReKi)                       :: RN(1)
INTEGER                          :: I           ! loop counter
INTEGER                          :: NumSeeds    ! number of seeds in the intrinsic random number generator 
INTEGER                          :: Sttus       ! allocation status


IF (RNG_type == 'NORMAL') THEN


      ! determine the number of seeds necessary (gfortran needs 8 or 12 seeds, not just 2)
      
   CALL RANDOM_SEED ( SIZE = NumSeeds )
      
   IF ( NumSeeds /= 2 ) THEN
      CALL ProgWarn( ' The random number generator in use differs from the original code provided by NREL. This pRNG uses ' &
                        //TRIM(Int2LStr(NumSeeds))//' seeds instead of the 2 in the TurbSim input file.')
   END IF

   IF ( .NOT. ALLOCATED( RandSeedAry ) ) THEN
      ALLOCATE ( RandSeedAry ( NumSeeds ), STAT=Sttus)
      IF (Sttus/= 0 ) THEN
         CALL ProgAbort( ' Error allocating space for RandSeedAry array.' )
         RETURN
      END IF   
   END IF


         ! We'll just populate this with odd seeds = Seed(1) and even seeds = Seed(2)
   DO I = 1,NumSeeds,2
      RandSeedAry(I) = RandSeed(1)
   END DO
   DO I = 2,NumSeeds,2
      RandSeedAry(I) = RandSeed(2)
   END DO
                     
                  
   CALL RANDOM_SEED ( PUT=RandSeedAry )
   

ELSEIF (RNG_type == 'RANLUX') THEN

   CALL RLuxGo ( LuxLevel, ABS( RandSeed(1) ), 0, 0 )

ELSE

      ! A quick and dirty way to get three random seeds for u, v, and w
      ! This implementation allows comparisons with Neil's SNLWIND-3D

   CALL ARand( RandSeed(1), RN,1,1)

   RandSeed(2) = RandSeed(1)+1
   CALL ARand( RandSeed(2), RN,1,1)

   RandSeed(3) = RandSeed(2)+1
   CALL ARand( RandSeed(3), RN,1,1)

ENDIF

END SUBROUTINE RndInit
!=======================================================================
SUBROUTINE RndJetHeight( RandNum )
! This function uses the Pearson IV equation

IMPLICIT                            NONE

REAL(ReKi), PARAMETER            :: a =   0.021548497
REAL(ReKi), PARAMETER            :: b = -13.173289
REAL(ReKi), PARAMETER            :: c =  13.43201034
REAL(ReKi), PARAMETER            :: d =   0.896588964
REAL(ReKi), PARAMETER            :: e =  -0.71128456

REAL(ReKi), PARAMETER            :: MaxVALUE =  120     ! Maximum value of returned variate
REAL(ReKi), PARAMETER            :: MinVALUE = -160     ! Minimum value of returned variate
REAL(ReKi), PARAMETER            :: Parms(5) = (/ a, b, c, d, e /)
REAL(ReKi), INTENT(OUT)          :: RandNum             ! numbers distributed with the pdf above
REAL(ReKi), PARAMETER            :: RangeFn(2)  = (/ MinVALUE, MaxVALUE /)


   CALL RndPearsonIV( RandNum, Parms, RangeFn )


END SUBROUTINE RndJetHeight
!=======================================================================
SUBROUTINE RndModLogNorm( RandNum, Height )
   ! This subroutine generates a random variate with a PDF of the form
   ! f(x) = A + B*exp(-0.5*(ln(x/C)/D)^2)

!BJJ use Rnd4ParmLogNorm()
IMPLICIT                            NONE

REAL(ReKi), INTENT(OUT)          :: RandNum        ! Near-Log-Normally distributed numbers
REAL(ReKi), INTENT(IN), OPTIONAL :: Height         ! height (in meters), determining what parameters to use

   ! Internal variables

REAL(ReKi), PARAMETER            :: A(3) = (/-0.0041046,   -0.00566512,  -0.00216964 /)
REAL(ReKi), PARAMETER            :: B(3) = (/0.162945643,  0.278246235,  0.113718973 /)
REAL(ReKi), PARAMETER            :: C(3) = (/0.67493672,   0.203262077,  3.211606394 /)
REAL(ReKi), PARAMETER            :: D(3) = (/2.391316782,  2.715789776,  1.700298642 /)
REAL(ReKi)                       :: G              ! The function g(x) = f(x)/B
REAL(ReKi)                       :: RN (2)         ! Two random numbers for the acceptance-rejection method


INTEGER                          :: Count
INTEGER                          :: Indx
INTEGER, PARAMETER               :: MaxIter = 10000  ! Max number of iterations to converge (so we don't get an infinite loop)
INTEGER, PARAMETER               :: MaxTime = 600    ! Maximum value of returned value (the data used to compute A,B,C,D is valid up to 600 s.)

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

   CALL RndUnif( RN(1) )
   CALL RndUnif( RN(2) )

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
SUBROUTINE RndNorm( RandNormNum, mu, sigma )

IMPLICIT                            NONE

REAL(ReKi), INTENT(OUT)          :: RandNormNum    ! Normally distributed numbers
REAL(ReKi), INTENT(IN), OPTIONAL :: mu             ! mean of the distributed numbers - DEFAULT IS 0.0
REAL(ReKi), INTENT(IN), OPTIONAL :: sigma          ! standard deviation of the distributed numbers - DEFAULT IS 1.0

   ! Internal variable

REAL(ReKi)                       :: RN (2)         ! Two random numbers


      ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )

   CALL RndUnif( RN(1) )
   CALL RndUnif( RN(2) )

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
SUBROUTINE RndNWTCpkCTKE( RandNum )
   ! This subroutine generates a random variate with a PDF of the form
   ! f(x) = A + B * EXP( (-X + C + D - D*E*EXP(-( X + D*LOG(E) - C)/D)) / (D*E)
   ! Maximum, f(C) = A + B
   ! Uses the Acceptance-Rejection: f(x) = Cf*h(x)*g(x)
   !   where h(x) = 1 / (150-30) (for our domain)
   !         g(x) = f(x)/(A+B), and
   !         Cf   = (150-30)*(A+B)

IMPLICIT                            NONE

REAL(ReKi), INTENT(OUT)          :: RandNum        ! numbers distributed with the pdf above

   ! Internal variables

REAL(ReKi), PARAMETER            :: A =   0.000500609
REAL(ReKi), PARAMETER            :: B =   0.286202317
REAL(ReKi), PARAMETER            :: C = -38.4131676
REAL(ReKi), PARAMETER            :: D = 244.6908697
REAL(ReKi), PARAMETER            :: E =   0.02115063

REAL(ReKi), PARAMETER            :: fMAX = A + B     ! Max(f(x))
REAL(ReKi)                       :: Gx               ! The function g(x) = f(x)/B
REAL(ReKi), PARAMETER            :: MaxVALUE = 150.0 ! Maximum value of returned variate
REAL(ReKi), PARAMETER            :: MinVALUE =  30.0 ! Minimum value of returned variate
REAL(ReKi)                       :: RN               ! A random number for the acceptance-rejection method

INTEGER                          :: Count
INTEGER, PARAMETER               :: MaxIter = 10000 ! Max number of iterations to converge (so we don't get an infinite loop)

Count = 1

      ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )

DO WHILE (Count < MaxIter)

   CALL RndUnif( RN )        ! Generate RN from U(0,1)

   CALL RndUnif( RandNum )   ! Generate RandNum from h(y) = 1 / (MaxVALUE - MinVALUE)
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
SUBROUTINE RndNWTCuStar( RandNum )
   ! This subroutine generates a random variate with a PDF of the form
   ! f(x) = (A + Cx + Ex^2 + Gx^3) / (1 + Bx + Dx^2 + Fx^3 + Hx^4)
   ! using the acceptance/rejection method.

IMPLICIT                            NONE

REAL(ReKi), INTENT(OUT)          :: RandNum        ! numbers distributed with the pdf above

   ! Internal variables

REAL(ReKi), PARAMETER            :: A =   4.50581      ! Scaling parameters for the pdf
REAL(ReKi), PARAMETER            :: B =  -0.60722      ! Scaling parameters for the pdf
REAL(ReKi), PARAMETER            :: C = -14.23826      ! Scaling parameters for the pdf
REAL(ReKi), PARAMETER            :: D =  -0.96523      ! Scaling parameters for the pdf
REAL(ReKi), PARAMETER            :: E =  15.92342      ! Scaling parameters for the pdf
REAL(ReKi), PARAMETER            :: F =  14.41326      ! Scaling parameters for the pdf
REAL(ReKi), PARAMETER            :: G =  -6.16188      ! Scaling parameters for the pdf
REAL(ReKi), PARAMETER            :: H =  -4.82923      ! Scaling parameters for the pdf

REAL(DbKi)                       :: Gx              ! The function g(x) = f(x)/A
REAL(ReKi), PARAMETER            :: MaxUstar = 1.0  ! Maximum value of returned value (the data used to compute A,B,C,D is valid up to 600 s.)
REAL(DbKi)                       :: RandNum2        ! RandNum**2
REAL(DbKi)                       :: RandNum3        ! RandNum**3
REAL(DbKi)                       :: RandNum4        ! RandNum**4
REAL(ReKi)                       :: RN              ! A random number for the acceptance-rejection method

INTEGER                          :: Count
INTEGER, PARAMETER               :: MaxIter = 10000 ! Max number of iterations to converge (so we don't get an infinite loop)

Count = 1

      ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )

DO WHILE (Count < MaxIter)

   CALL RndUnif( RN )
   CALL RndUnif( RandNum )

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
SUBROUTINE RndPearsonIV( RandNum, Parms, FnRange )
! This function uses the Pearson IV equation to generate a deviate from that distribution
! Equation 8186 in TableCurve

IMPLICIT                            NONE

REAL(ReKi), INTENT(IN)           :: Parms(5)           ! 1=a, 2=b, 3=c, 4=d, 5=e
REAL(ReKi), INTENT(IN)           :: FnRange(2)

REAL(ReKi)                       :: fMAX                ! Max(f(x)) ! occurs at f(b)
REAL(ReKi)                       :: Gx                  ! The function g(x) = f(x)/fMAX
REAL(ReKi)                       :: MaxVALUE            ! Maximum value of returned variate
REAL(ReKi)                       :: MinVALUE            ! Minimum value of returned variate
REAL(ReKi)                       :: n                   ! A temporary variable for calculating the function values
REAL(ReKi), INTENT(OUT)          :: RandNum             ! numbers distributed with the pdf above
REAL(ReKi)                       :: RN                  ! A random number for the acceptance-rejection method

INTEGER                          :: Count
INTEGER, PARAMETER               :: MaxIter = 10000     ! Max number of iterations to converge (so we don't get an infinite loop)

fMAX     = Parms(1)
MaxVALUE = MAX(FnRange(1),FnRange(2))
MinVALUE = MIN(FnRange(1),FnRange(2))

Count = 1

      ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )

DO WHILE (Count < MaxIter)

   CALL RndUnif( RN )        ! Generate RN from U(0,1)

   CALL RndUnif( RandNum )   ! Generate RandNum from h(y) = 1 / (MaxVALUE - MinVALUE)
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
SUBROUTINE RndpkCTKE_WFTA( RandNum )
! Calculates a deviate from a distribution with pdf:
! f(x) = (a + c*x^2 + e*x^4)/(1. + b*x^2 + d*x^4 + f*x^6), where

IMPLICIT                            NONE

!67m and 85m
REAL(ReKi), PARAMETER            :: A = 0.30985506
REAL(ReKi), PARAMETER            :: B = 0.006902104
REAL(ReKi), PARAMETER            :: C =-0.00206008
REAL(ReKi), PARAMETER            :: D = 1.28884E-05
REAL(ReKi), PARAMETER            :: E = 5.71475E-06
REAL(ReKi), PARAMETER            :: F = 8.70606E-07

REAL(ReKi)                       :: fMAX                ! Max(f(x)) = f(MinVALUE)
REAL(ReKi)                       :: Gx                  ! The function g(x) = f(x)/fMAX
REAL(ReKi), PARAMETER            :: MaxVALUE = 22.      ! Maximum value of returned variate
REAL(ReKi), PARAMETER            :: MinVALUE =  4.      ! Minimum value of returned variate
REAL(ReKi), INTENT(OUT)          :: RandNum             ! numbers distributed with the pdf above
REAL(ReKi)                       :: RN                  ! A random number for the acceptance-rejection method
REAL(ReKi)                       :: x2                  ! A temp variable for x^2

INTEGER                          :: Count
INTEGER, PARAMETER               :: MaxIter = 10000     ! Max number of iterations to converge (so we don't get an infinite loop)


x2   = MinVALUE**2
fMAX = (A + C*x2 + E*x2**2)/(1.0 + B*x2 + D*x2**2 + F*x2**3)

Count = 1

      ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )

DO WHILE (Count < MaxIter)

   CALL RndUnif( RN )        ! Generate RN from U(0,1)
   CALL RndUnif( RandNum )   ! Generate RandNum from h(y) = 1 / (MaxVALUE - MinVALUE)
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
SUBROUTINE RndPolyFit( RandNum, Coeffs, FnRange, fMAX)
! Calculates a deviate from a distribution with pdf:
! f(x) = (Coeffs(1) + Coeffs(3)*x + Coeffs(5)*x^2 + Coeffs(7)*x^3 + Coeffs(9)*x^4 + Coeffs(11)*x^5) / &
!        (       1. + Coeffs(2)*x + Coeffs(4)*x^2 + Coeffs(6)*x^3 + Coeffs(8)*x^4 + Coeffs(10)*x^5)
! This equation covers the following (plus others) from Table Curve\SysStat:
! Eqn 7906, Eqn 7907, Eqn 7908, & Eqn 7909

IMPLICIT                            NONE


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

   CALL RndUnif( RN )        ! Generate RN from U(0,1)
   CALL RndUnif( RandNum )   ! Generate RandNum from h(y) = 1 / (MaxVALUE - MinVALUE)
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
SUBROUTINE RndTcohLLJ( RandNum, Height )
! Calculates a deviate from a distribution with pdf:
! f(x) = EXP(A + B*SQRT(x) + C*LOG(x) )  !Eqn 1376

IMPLICIT                            NONE

!67m and 85m
REAL(ReKi), PARAMETER            :: A(2) = (/ -1.34064396 , -1.17577736  /)
REAL(ReKi), PARAMETER            :: B(2) = (/ -0.26996911 , -0.23056567  /)
REAL(ReKi), PARAMETER            :: C(2) = (/ -0.57793906 , -0.69871145  /)

REAL(ReKi)                       :: fMAX                ! Max(f(x)) = f(MinValue)
REAL(ReKi)                       :: Gx                  ! The function g(x) = f(x)/fMAX
REAL(ReKi)                       :: Height              ! The height of the center of the billow, in meters
REAL(ReKi), PARAMETER            :: MaxVALUE = 600.     ! Maximum value of returned variate
REAL(ReKi), PARAMETER            :: MinVALUE = 2.5      ! Minimum value of returned variate
REAL(ReKi), INTENT(OUT)          :: RandNum             ! numbers distributed with the pdf above
REAL(ReKi)                       :: RN                  ! A random number for the acceptance-rejection method

INTEGER                          :: Count
INTEGER                          :: Indx
INTEGER, PARAMETER               :: MaxIter = 10000     ! Max number of iterations to converge (so we don't get an infinite loop)

IF (Height < 76) THEN
   Indx = 1
ELSE
   Indx = 2
ENDIF

fMAX = EXP(A(Indx) + B(Indx)*SQRT(MinVALUE) + C(Indx)*LOG(MinVALUE) )

Count = 1

      ! Generate a normal distribution on (0,1) from a uniform distribution ( ACTUALLY [0,1) )

DO WHILE (Count < MaxIter)

   CALL RndUnif( RN )        ! Generate RN from U(0,1)
   CALL RndUnif( RandNum )   ! Generate RandNum from h(y) = 1 / (MaxVALUE - MinVALUE)
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
SUBROUTINE RndTcoh_WF( RandNum )

USE                                 TSMods, ONLY: TurbModel ! The name of the turbulence model (valid for WF_UPW or WF_07D)

IMPLICIT                            NONE
REAL(ReKi)                       :: FnRange( 2)         ! the min and max values of the returned deviate
REAL(ReKi)                       :: ParmsLN( 4)         ! parameters for the 4-parameter Log-Norm function
REAL(ReKi)                       :: ParmsPIV(5)         ! parameters for the Pearson IV function
REAL(ReKi), INTENT(OUT)          :: RandNum             ! numbers distributed with the pdf above


SELECT CASE ( TurbModel )
   CASE ( 'WF_UPW' )
      ParmsLN = (/ 0.0, 0.132537201, 0.348791907, 1.781668096 /)  !parm(1) = 0 b/c there's no offset... it's a 3-parameter function
      FnRange = (/ 0.4, 200.0 /)

      CALL Rnd4ParmLogNorm( RandNum, ParmsLN, FnRange )

   CASE ( 'WF_07D' )
      ParmsPIV = (/ 0.108721975, 5.705449915, 2.769408844, 0.475906651, 1.067616671 /)
      FnRange  = (/ 0.6, 30.0 /)

      CALL RndPearsonIV( RandNum, ParmsPIV, FnRange )

   CASE ( 'WF_14D' )
      ParmsPIV = (/ 0.080526074, 13.51637204, 6.391924365, 1.197332751, 0.390220799 /)
      FnRange  = (/ 5.0, 40.0 /)

      CALL RndPearsonIV( RandNum, ParmsPIV, FnRange )

END SELECT

END SUBROUTINE RndTcoh_WF
!=======================================================================
SUBROUTINE RndUnif( RandUnifNum )

!This subroutine produces uniformly distributed random numbers, based on
!the pRNG that is requested in TurbSim's input file.  This routine assumes
!that the random number generator has been initialized earlier in the main
!program.

USE                                 Ran_Lux_Mod
USE                                 TSMods

IMPLICIT                            NONE

REAL(ReKi), INTENT(OUT)          :: RandUnifNum        ! Uniformly distributed numbers
REAL(ReKi)                       :: RN(1)

IF (RNG_type == 'NORMAL') THEN

   CALL RANDOM_NUMBER( RN )

ELSEIF (RNG_type == 'RANLUX') THEN

   CALL RanLux ( RN )

ELSE

   CALL ARand( RandSeed(1), RN, 1,  1)

ENDIF


RandUnifNum = RN(1)


END SUBROUTINE RndUnif
!=======================================================================
SUBROUTINE ST_Turb ( Ht, Ucmp, Spec )

   ! This subroutine defines the 3-D turbulence spectrum that can be expected over flat,
   ! homogeneous terrain as developed by RISO authors Hojstrup, Olesen, and Larsen.
   ! The use of this subroutine requires that variables have the units of meters and seconds.

USE                     TSMods

IMPLICIT                NONE

   ! Passed variables

REAL(ReKi),INTENT(IN) :: Ht                      ! Height
REAL(ReKi),INTENT(IN) :: Ucmp                    ! Longitudinal Velocity
REAL(ReKi)            :: Spec     (:,:)          ! Working array for PSD

   ! Internal variables

REAL(ReKi), PARAMETER :: Exp1  = 5.0 / 3.0
REAL(ReKi), PARAMETER :: Exp2  = 2.0 / 3.0
REAL(ReKi), PARAMETER :: Exp3  = 3.0 / 2.0
REAL(ReKi)            :: fi       ! Temporary variable for calculation of Spec
REAL(ReKi)            :: fr       ! Temporary variable for calculation of Spec
REAL(ReKi)            :: HtZI     ! Temporary variable for calculation of Spec
REAL(ReKi)            :: HtZI2    ! Temporary variable for calculation of Spec
REAL(ReKi)            :: phiE
REAL(ReKi)            :: phiM     ! Non-Dimensional Wind Shear
REAL(ReKi)            :: ps_h
REAL(ReKi)            :: ps_l
REAL(ReKi)            :: tmpF     ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpN     ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpX     ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpXX    ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpPhi   ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpZIL   ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpZIU   ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpZU    ! Temporary variable for calculation of Spec
REAL(ReKi)            :: UDen
REAL(ReKi)            :: uStar2   ! Temporary variable holding Ustar-squared
REAL(ReKi)            :: Ustar2F
REAL(ReKi)            :: VDen

INTEGER               :: I        ! DO LOOP counter

uStar2 = Ustar * Ustar

IF (zL >= 0) THEN

   ! Calculate NEUTRAL/STABLE spectral estimates

   phiM   =  1.0 + 4.7*(zL)            ! = q
   phiE   = (1.0 + 2.5*(zL)**0.6)**Exp3

   tmpPhi = uStar2 * ( (phiE / phiM)**Exp2 )
   tmpF   = Ht / (Ucmp * phiM)

   DO I = 1,NumFreq
      tmpX  = Freq(I)*tmpF             ! reduced frequency divided by q (q = phiM here)
      tmpXX = tmpX**Exp1
      tmpN  = tmpPhi / Freq(I) * tmpX  ! normalization factor used to obtain power spectrum components

      Spec(I,1) = tmpN * (79.0) / (1.0 + 263.0*tmpXX)
      Spec(I,2) = tmpN * (13.0) / (1.0 +  32.0*tmpXX)
      Spec(I,3) = tmpN * ( 3.5) / (1.0 +   8.6*tmpXX)
   ENDDO

ELSE
   ! Calculate UNSTABLE spectral estimates
   tmpZIL = (- ZI / L)**Exp2
   HtZI   = Ht / ZI

   HtZI2  = (1.0 - HtZI)**2
   tmpZU  = Ht / Ucmp
   tmpZIU = ZI / Ucmp
   UDen   = 1.0 + 15.0*HtZI
   VDen   = 1.0 +  2.8*HtZI

   DO I = 1,NumFreq

      Fi   = Freq(I)*tmpZIU
      tmpF = Freq(I)*tmpZU                ! reduced frequency

      Ustar2F = uStar2/Freq(I)

      ! u component
      Fr   = tmpF / UDen
      ps_l = ( (  0.5*Fi) / (1.0 +  2.2*  Fi**Exp1)) * tmpZIL
      ps_h = ( (105.0*Fr) / (1.0 + 33.0*Fr )**Exp1 ) * HtZI2 / UDen**Exp2

      Spec(I,1) = (ps_l + ps_h) * Ustar2F

      ! v component
      Fr   = tmpF / VDen
      ps_l = ( ( 0.95*Fi) / (1.0 +  2.0*Fi)**Exp1 ) * tmpZIL
      ps_h = ( (17.00*Fr) / (1.0 +  9.5*Fr)**Exp1 ) * HtZI2 / VDen**Exp2

      Spec(I,2) = (ps_l + ps_h) * Ustar2F

      ! w component
      tmpN = SQRT( (tmpF**2 + (0.3*HtZI)**2) / (tmpF**2 + 0.0225) )

      ps_l = tmpN * ( (0.95*Fi  ) / (1.0 +  2.0*Fi )**Exp1 ) * tmpZIL
      ps_h =        ( (2.00*tmpF) / (1.0 +  5.3*tmpF**Exp1)) * HtZI2

      Spec(I,3) = (ps_l + ps_h) * Ustar2F

   ENDDO

ENDIF


RETURN
END SUBROUTINE ST_Turb
!=======================================================================
SUBROUTINE Tidal_Turb ( Shr_DuDz, Spec )
   ! HYDROTURBSIM specific.

   ! This subroutine defines the 3-D turbulence expected in a tidal channel.
   ! It is similar to the 'smooth' spectral model (RISO; Hojstrup, Olesen and Larsen) for wind,
   ! but is scaled by the TKE (SigmaU**2), and du/dz rather than Ustar and u/z.  
   ! The fit is based on data from Puget Sound, estimated by L. Kilcher.
   ! The use of this subroutine requires that variables have the units of meters and seconds.
   ! Note that this model does not require height.
   

USE                     TSMods

IMPLICIT                NONE

   ! Passed variables

REAL(ReKi),INTENT(IN) :: Shr_DuDz                ! Shear (du/dz)
REAL(ReKi)            :: Spec     (:,:)          ! Working array for PSD

   ! Internal variables

REAL(ReKi), PARAMETER :: Exp1  = 5.0 / 3.0
REAL(ReKi)            :: tmpX                   ! Temporary variable for calculation of Spec
REAL(ReKi)            :: tmpvec(3)              ! Temporary vector for calculation of Spec
REAL(ReKi)            :: tmpa (3)               ! Spectra coefficients
REAL(ReKi)            :: tmpb (3)               ! Spectra coefficients
INTEGER               :: I                      ! DO LOOP counter


SELECT CASE ( TRIM(TurbModel) )
   CASE ( 'TIDAL' )
      tmpa = (/ 0.193, 0.053 , 0.0362 /)*TwoPi ! These coefficients were calculated using Shr_DuDz in units of 'radians', so we multiply these coefficients by 2*pi.
      tmpb = (/ 0.201, 0.0234, 0.0124 /)*(TwoPi**Exp1)
   CASE ( 'RIVER' )
      ! THESE ARE NOT VERIFIED YET!!!, therefore they are undocumented.
      tmpa = (/ 0.081, 0.056 , 0.026 /)*TwoPi
      tmpb = (/ 0.16, 0.025, 0.020 /)*(TwoPi**Exp1)
END SELECT

tmpvec = tmpa*(/Sigma_U2, Sigma_V2, Sigma_W2/)/Shr_DuDz

DO I = 1,NumFreq  
   tmpX  = (Freq(I)/Shr_DuDz)**Exp1
   Spec(I,1) = tmpvec(1) / (1.0 + tmpb(1)*tmpX)
   Spec(I,2) = tmpvec(2) / (1.0 + tmpb(2)*tmpX)
   Spec(I,3) = tmpvec(3) / (1.0 + tmpb(3)*tmpX)
ENDDO

RETURN
END SUBROUTINE Tidal_Turb
!=======================================================================
!!!!bjj Start of proposed change
!!!SUBROUTINE TestSpec ( Ht, Ucmp, Spec )
!!!
!!!
!!!USE                     TSMods
!!!
!!!IMPLICIT                NONE
!!!
!!!      ! Passed variables
!!!
!!!REAL(ReKi),INTENT(IN) :: Ht                      ! Input: Height (Should be HubHt)
!!!REAL(ReKi),INTENT(IN) :: Ucmp                    ! Input: Velocity (Should be UHub)
!!!REAL(ReKi)            :: Spec   (:,:)            ! Output: target spectrum
!!!
!!!INTEGER               :: I
!!!INTEGER               :: IVec
!!!
!!!
!!!   ! Create the spectrum.
!!!
!!!DO IVec = 1,3
!!!
!!!   DO I = 1,NumFreq
!!!      Spec(I,IVec) = 0.0
!!!   ENDDO !I
!!!   !I = INT( NumFreq/2 )
!!!   I = INT( 100 )
!!!   Spec( I, IVec ) = 1/Freq(1)
!!!
!!!print *, 'Test Spectra: sine wave with frequency ', Freq(I), ' Hz.'
!!!
!!!ENDDO !IVec
!!!
!!!
!!!RETURN
!!!END SUBROUTINE TestSpec
!!!!=======================================================================
!!!!bjj End of proposed change
SUBROUTINE TS_Abort ( Message )
   ! This routine outputs fatal warning messages and ends the program.
      
   USE                              TSMods, ONLY:   US                           ! unit number for summary file

      ! Argument declarations.

   CHARACTER(*), INTENT(IN)      :: Message                                      ! Warning message.


      ! Write the message to the summary file
      
   WRITE (US, "(/'ERROR:  ', A / )") Message
   WRITE (US, "('ABORTING PROGRAM.')" )

      ! Write the message to the screen
   CALL ProgAbort ( Message )

RETURN

END SUBROUTINE TS_Abort
!=======================================================================
SUBROUTINE TS_Warn ( Message, WrSum )

   ! This routine outputs non-fatal warning messages and returns to the calling routine.
      
   USE                              TSMods, ONLY:   US                           ! unit number for summary file

      ! Argument declarations.

   CHARACTER(*), INTENT(IN)      :: Message                                      ! Warning message.
   LOGICAL,      INTENT(IN)      :: WrSum                                        ! Whether to print a message in the sum file.

      ! Write the message to the screen

   CALL WrScr( '' )
   CALL ProgWarn ( Message )
   CALL WrScr( '' )


      ! Write the message to the summary file if requested

   IF ( WrSum ) THEN
      WRITE (US, "(/'WARNING:  ', A / )") Message
   ENDIF


RETURN
END SUBROUTINE TS_Warn ! ( Message )
!=======================================================================
SUBROUTINE UsrSpec ( Ht, Ucmp, Spec )


   USE                     TSMods

   IMPLICIT                NONE

         ! Passed variables

   REAL(ReKi),INTENT(IN) :: Ht    ! local height
   REAL(ReKi),INTENT(IN) :: Ucmp  ! local wind speed
   REAL(ReKi)            :: Spec     (:,:)

      ! Internal variables

   REAL(ReKi)            :: Tmp
   
   
   INTEGER               :: I
   INTEGER               :: Indx
   INTEGER               :: J

      ! --------- Interpolate to the desired frequencies ---------------
   
   Indx = 1;

   DO I=1,NumFreq

      IF ( Freq(I) <= Freq_USR(1) ) THEN
         Spec(I,1) = Uspec_USR(1)
         Spec(I,2) = Vspec_USR(1)
         Spec(I,3) = Wspec_USR(1)
      ELSEIF ( Freq(I) >= Freq_USR(NumUSRf) ) THEN
         Spec(I,1) = Uspec_USR(NumUSRf)
         Spec(I,2) = Vspec_USR(NumUSRf)
         Spec(I,3) = Wspec_USR(NumUSRf)
      ELSE
      
            ! Find the two points between which the frequency lies

         DO J=(Indx+1),NumUSRf
            IF ( Freq(I) <= Freq_USR(J) ) THEN
               Indx = J-1
               
                  ! Let's just do a linear interpolation for now

               Tmp  = (Freq(I) - Freq_USR(Indx)) / ( Freq_USR(Indx) - Freq_USR(J) )

               Spec(I,1) = Tmp * ( Uspec_USR(Indx) - Uspec_USR(J) ) + Uspec_USR(Indx)
               Spec(I,2) = Tmp * ( Vspec_USR(Indx) - Vspec_USR(J) ) + Vspec_USR(Indx)
               Spec(I,3) = Tmp * ( Wspec_USR(Indx) - Wspec_USR(J) ) + Wspec_USR(Indx)  
               
               EXIT
            ENDIF
         ENDDO ! J

      ENDIF

   ENDDO ! I

   RETURN


END SUBROUTINE UsrSpec
!=======================================================================
SUBROUTINE vonKrmn ( Ht, Ucmp, Spec )


   ! This subroutine defines the von Karman PSD model.
   ! The use of this subroutine requires that all variables have the units of meters and seconds.


USE                     TSMods

IMPLICIT                NONE

      ! Passed variables

REAL(ReKi),INTENT(IN) :: Ht    ! local height
REAL(ReKi),INTENT(IN) :: Ucmp  ! local wind speed
REAL(ReKi)            :: Spec     (:,:)

      ! Internal variables

REAL(ReKi),PARAMETER  :: Exp1 =  5.0/6.0
REAL(ReKi),PARAMETER  :: Exp2 = 11.0/6.0
REAL(ReKi)            :: FLU2
REAL(ReKi)            :: L1_U
REAL(ReKi)            :: Lambda
REAL(ReKi)            :: Lvk        ! von Karman length scale
REAL(ReKi)            :: Sigma      ! Standard deviation
REAL(ReKi)            :: SigmaL1_U
REAL(ReKi)            :: Tmp

INTEGER               :: I

   ! Define isotropic integral scale.
IF ( ALLOCATED( L_USR ) ) THEN
   IF ( Ht <= Z_USR(1) ) THEN
      Lvk = L_USR(1)   ! Extrapolation: nearest neighbor for heights below minimum height specified
   ELSEIF ( Ht >= Z_USR(NumUSRz) ) THEN
      Lvk = L_USR(NumUSRz)  ! Extrapolation: nearest neighbor for heights above maximum height specified
   ELSE !Interpolation: linear between user-defined height/integral scale curves
      DO I=2,NumUSRz
         IF ( Ht <= Z_USR(I) ) THEN
            Lvk = (Ht - Z_USR(I-1)) * ( L_USR(I-1) - L_USR(I) ) / ( Z_USR(I-1) - Z_USR(I) ) + L_USR(I-1)
            EXIT
         ENDIF
      ENDDO
   ENDIF
ELSE
   IF ( Ht  <  150.0 )  THEN
      Lambda = 0.7*Ht
   ELSE
      Lambda = 105.0
   ENDIF
   Lvk = 3.5*Lambda
ENDIF

   ! Define isotropic integral scale.
IF ( ALLOCATED( Sigma_USR ) ) THEN
   IF ( Ht <= Z_USR(1) ) THEN
      Sigma = Sigma_USR(1)
   ELSEIF ( Ht >= Z_USR(NumUSRz) ) THEN
      Sigma = Sigma_USR(NumUSRz)
   ELSE
      DO I=2,NumUSRz
         IF ( Ht <= Z_USR(I) ) THEN
            Sigma = (Ht - Z_USR(I-1)) * ( Sigma_USR(I-1) - Sigma_USR(I) ) / ( Z_USR(I-1) - Z_USR(I) ) + Sigma_USR(I-1)
            EXIT
         ENDIF
      ENDDO
   ENDIF
ELSE
    Sigma = Ustar*2.15 !bjj: BONNIE, make sure this is defined, or else define ustar for this model...
ENDIF


L1_U   = Lvk/Ucmp
SigmaL1_U = 2.0*Sigma*Sigma*L1_U

DO I=1,NumFreq

   FLU2      = ( Freq(I)*L1_U )**2
   Tmp       = 1.0 + 71.0*FLU2

   Spec(I,1) = (StdScale(1)**2)*2.0*SigmaL1_U/Tmp**Exp1
   Spec(I,2) = SigmaL1_U*( 1.0 + 189.0*FLU2 )/Tmp**Exp2
   Spec(I,3) = Spec(I,2)

   Spec(I,2) = (StdScale(2)**2)*Spec(I,2)
   Spec(I,3) = (StdScale(3)**2)*Spec(I,3)

ENDDO ! I

RETURN
END SUBROUTINE vonKrmn
!=======================================================================
SUBROUTINE WrBinBLADED(USig, VSig, WSig)

USE                         TSMods


IMPLICIT                    NONE

REAL(ReKi)                  :: U_C1                          ! Scale for converting BLADED U data
REAL(ReKi)                  :: U_C2                          ! Offset for converting BLADED U data
REAL(ReKi),INTENT(INOUT)    :: USig                          ! Standard deviation of U
REAL(ReKi)                  :: V_C                           ! Scale for converting BLADED V data
REAL(ReKi),INTENT(INOUT)    :: VSig                          ! Standard deviation of V
REAL(ReKi)                  :: W_C                           ! Scale for converting BLADED W data
REAL(ReKi),INTENT(INOUT)    :: WSig                          ! Standard deviation of W
REAL(ReKi)                  :: TI(3)                         ! Turbulence intensity for scaling data
REAL(ReKi)                  :: ZTmp                          ! This is the vertical center of the grid

INTEGER(B4Ki)               :: CFirst
INTEGER(B4Ki)               :: CLast
INTEGER(B4Ki)               :: CStep
INTEGER(B4Ki)               :: I
INTEGER(B4Ki)               :: II
INTEGER(B4Ki)               :: IT
INTEGER(B4Ki)               :: IY
INTEGER(B4Ki)               :: IZ

INTEGER(B4Ki)               :: IP
INTEGER(B2Ki)               :: TmpVarray(3*NumGrid_Y*NumGrid_Z) ! This array holds the normalized velocities before being written to the binary file
INTEGER(B2Ki),ALLOCATABLE   :: TmpTWRarray(:)                   ! This array holds the normalized tower velocities

INTEGER                     :: AllocStat



      ! Put normalizing factors into the summary file.  The user can use them to
      ! tell a simulation program how to rescale the data.

   USig  = MAX(100.0*Tolerance, USig)
   VSig  = MAX(100.0*Tolerance, VSig)
   WSig  = MAX(100.0*Tolerance, WSig)

   TI(1) = USig / UHub
   TI(2) = VSig / UHub
   TI(3) = WSig / UHub

   FormStr = "(//,'Normalizing Parameters for Binary Data (approximate statistics):',/)"
   WRITE (US,FormStr)

   FormStr = "(3X,A,' =',F9.4,A)"
   WRITE (US,FormStr)  'UBar ', UHub, ' m/s'
   WRITE (US,FormStr)  'TI(u)', 100.0*TI(1), ' %'
   WRITE (US,FormStr)  'TI(v)', 100.0*TI(2), ' %'
   WRITE (US,FormStr)  'TI(w)', 100.0*TI(3), ' %'

   Ztmp  = ( HubHt - GridHeight / 2.0 - Z(1) )  ! This is the grid offset

   WRITE (US,'()')
   WRITE (US,FormStr)  'Height Offset', Ztmp, ' m'
   WRITE (US,FormStr)  'Grid Base    ', Z(1), ' m'

      ! Calculate some numbers for normalizing the data.

   U_C1 = 1000.0/( UHub*TI(1) )
   U_C2 = 1000.0/TI(1)
   V_C  = 1000.0/( UHub*TI(2) )
   W_C  = 1000.0/( UHub*TI(3) )


   ZTmp     = Z(1) + GridHeight/2.0  !This is the vertical center of the grid

   IF ( WrBLFF )  THEN

      CALL WrScr ( ' Generating BLADED binary time-series file "'//TRIM( RootName )//'.wnd"' )

               ! Put header information into the binary data file.

      WRITE (UBFFW)   INT(  -99          , B2Ki )               ! -99 = New Bladed format
      WRITE (UBFFW)   INT(    4          , B2Ki )               ! 4 = improved von karman (but needed for next 7 inputs)
      WRITE (UBFFW)   INT(    3          , B4Ki )               ! 3 = number of wind components
      WRITE (UBFFW)  REAL( Latitude      , SiKi )               ! Latitude (degrees)
      WRITE (UBFFW)  REAL(   z0          , SiKi )               ! Roughness length (m)
      WRITE (UBFFW)  REAL( Ztmp          , SiKi )               ! Reference Height (m) ( Z(1) + GridHeight / 2.0 )
      WRITE (UBFFW)  REAL( 100.0*TI(1)   , SiKi )               ! Longitudinal turbulence intensity (%)
      WRITE (UBFFW)  REAL( 100.0*TI(2)   , SiKi )               ! Lateral turbulence intensity (%)
      WRITE (UBFFW)  REAL( 100.0*TI(3)   , SiKi )               ! Vertical turbulence intensity (%)

      WRITE (UBFFW)  REAL( GridRes_Z     , SiKi )               ! grid spacing in vertical direction, in m
      WRITE (UBFFW)  REAL( GridRes_Y     , SiKi )               ! grid spacing in lateral direction, in m
      WRITE (UBFFW)  REAL( TimeStep*UHub , SiKi )               ! grid spacing in longitudinal direciton, in m
      WRITE (UBFFW)   INT( NumOutSteps/2 , B4Ki )               ! half the number of points in alongwind direction
      WRITE (UBFFW)  REAL( UHub          , SiKi )               ! the mean wind speed in m/s
      WRITE (UBFFW)  REAL( 0             , SiKi )               ! the vertical length scale of the longitudinal component in m
      WRITE (UBFFW)  REAL( 0             , SiKi )               ! the lateral length scale of the longitudinal component in m
      WRITE (UBFFW)  REAL( 0             , SiKi )               ! the longitudinal length scale of the longitudinal component in m
      WRITE (UBFFW)   INT( 0             , B4Ki )               ! an unused integer
      WRITE (UBFFW)   INT( RandSeed(1)   , B4Ki )               ! the random number seed
      WRITE (UBFFW)   INT( NumGrid_Z     , B4Ki )               ! the number of grid points vertically
      WRITE (UBFFW)   INT( NumGrid_Y     , B4Ki )               ! the number of grid points laterally
      WRITE (UBFFW)   INT( 0             , B4Ki )               ! the vertical length scale of the lateral component, not used
      WRITE (UBFFW)   INT( 0             , B4Ki )               ! the lateral length scale of the lateral component, not used
      WRITE (UBFFW)   INT( 0             , B4Ki )               ! the longitudinal length scale of the lateral component, not used
      WRITE (UBFFW)   INT( 0             , B4Ki )               ! the vertical length scale of the vertical component, not used
      WRITE (UBFFW)   INT( 0             , B4Ki )               ! the lateral length scale of the vertical component, not used
      WRITE (UBFFW)   INT( 0             , B4Ki )               ! the longitudinal length scale of the vertical component, not used


         ! Compute parameters for ordering output for FF AeroDyn files. (This is for BLADED compatibility.)

      IF ( Clockwise )  THEN
         CFirst = NumGrid_Y
         CLast  = 1
         CStep  = -1
      ELSE
         CFirst = 1
         CLast  = NumGrid_Y
         CStep  = 1
      ENDIF


         ! Loop through time.

      DO IT=1,NumOutSteps  !Use only the number of timesteps requested originally

            ! Write out grid data in binary form.
         IP = 1
         DO IZ=1,NumGrid_Z
            DO IY=CFirst,CLast,CStep

               II = ( IZ - 1 )*NumGrid_Y + IY

               TmpVarray(IP)   = NINT( U_C1*V(IT,II,1) - U_C2, B2Ki )  ! Put the data into a temp array so that the WRITE() command works faster
               TmpVarray(IP+1) = NINT( V_C *V(IT,II,2)       , B2Ki )
               TmpVarray(IP+2) = NINT( W_C *V(IT,II,3)       , B2Ki )

               IP = IP + 3;
            ENDDO ! IY
         ENDDO ! IZ

         WRITE ( UBFFW )  TmpVarray(:) ! bjj: We cannot write the array including time because of stack overflow errors.. otherwise use compile option to put this on the heap instead of the stack

      ENDDO ! IT

      CLOSE ( UBFFW )


   ENDIF ! WrBLFF


   IF ( WrADTWR .AND. ( WrBLFF .OR. .NOT. WrADFF ) ) THEN
      IF ( ExtraHubPT ) THEN
         IZ = ZLim - NumGrid_Z - 1
         I  = NumGrid_Z*NumGrid_Y + 2
      ELSE
         IZ = ZLim - NumGrid_Z
         I  = NumGrid_Z*NumGrid_Y + 1
      ENDIF

      IF ( ExtraTwrPt ) THEN
         IY = I
         I  = I + 1
      ELSE
         IZ = IZ + 1
         IY = (NumGrid_Y / 2) + 1      !The grid location of the top tower point
      ENDIF


      CALL WrScr ( ' Generating tower binary time-series file "'//TRIM( RootName )//'.twr"' )


      WRITE (UATWR)  REAL( GridRes_Z ,       SiKi )         ! grid spacing in vertical direction, in m
      WRITE (UATWR)  REAL( TimeStep*UHub ,   SiKi )         ! grid spacing in longitudinal direciton, in m
      WRITE (UATWR)  REAL( Z(1) ,            SiKi )         ! The vertical location of the highest tower grid point in m
      WRITE (UATWR)   INT( NumOutSteps ,     B4Ki )         ! The number of points in alongwind direction
      WRITE (UATWR)   INT( IZ ,              B4Ki )         ! the number of grid points vertically
      WRITE (UATWR)  REAL( UHub ,            SiKi )         ! the mean wind speed in m/s
      WRITE (UATWR)  REAL( 100.0*TI(1),      SiKi )         ! Longitudinal turbulence intensity
      WRITE (UATWR)  REAL( 100.0*TI(2),      SiKi )         ! Lateral turbulence intensity
      WRITE (UATWR)  REAL( 100.0*TI(3),      SiKi )         ! Vertical turbulence intensity


      ALLOCATE ( TmpTWRarray( 3*(NTot-I+2) ) , STAT=AllocStat )

      IF ( AllocStat /= 0 )  THEN
         CALL TS_Abort ( 'Error allocating memory for temporary tower wind speed array.' )
      ENDIF


      DO IT=1,NumOutSteps

         TmpTWRarray(1) = NINT( U_C1*V(IT,IY,1) - U_C2 , B2Ki )
         TmpTWRarray(2) = NINT( V_C *V(IT,IY,2)        , B2Ki )
         TmpTWRarray(3) = NINT( W_C *V(IT,IY,3)        , B2Ki )

         IP = 4
         DO II = I,NTot    ! This assumes we have points in a single line along the center of the hub
            TmpTWRarray(IP  ) = NINT( U_C1*V(IT,II,1) - U_C2 , B2Ki )
            TmpTWRarray(IP+1) = NINT( V_C *V(IT,II,2)        , B2Ki )
            TmpTWRarray(IP+2) = NINT( W_C *V(IT,II,3)        , B2Ki )

            IP = IP + 3
         ENDDO    ! II

         WRITE (UATWR) TmpTWRarray(:)

      ENDDO ! IT

      CLOSE ( UATWR )

      IF ( ALLOCATED( TmpTWRarray) ) DEALLOCATE( TmpTWRarray )

   ENDIF !WrADWTR


END SUBROUTINE WrBinBLADED
!=======================================================================
SUBROUTINE WrBinTURBSIM

USE                       TSMods


IMPLICIT                  NONE

REAL(SiKi), PARAMETER     :: IntMax   =  32767.0
REAL(SiKi), PARAMETER     :: IntMin   = -32768.0
REAL(SiKi), PARAMETER     :: IntRng   = IntMax - IntMin ! Max Range of 2-byte integer

REAL(SiKi)                :: UOff                       ! Offset for the U component
REAL(SiKi)                :: UScl                       ! Slope  for the U component
REAL(ReKi)                :: VMax(3)                    ! Maximum value of the 3 wind components
REAL(ReKi)                :: VMin(3)                    ! Minimum value of the 3 wind components
REAL(SiKi)                :: VOff                       ! Offset for the V component
REAL(SiKi)                :: VScl                       ! Slope  for the V component
REAL(SiKi)                :: WOff                       ! Offset for the W component
REAL(SiKi)                :: WScl                       ! Slope  for the W component

INTEGER, PARAMETER        :: DecRound  = 3              ! Number of decimal places to round to
INTEGER                   :: IC                         ! counter for the velocity component of V
INTEGER                   :: II                         ! counter for the point on the grid/tower
INTEGER                   :: IT                         ! counter for the timestep
INTEGER(B4Ki)             :: LenDesc                    ! Length of the description string
INTEGER(B4Ki)             :: NumGrid                    ! Number of points on the grid
INTEGER(B4Ki)             :: NumTower                   ! Number of points on the tower
INTEGER(B4Ki)             :: TwrStart                   ! First index of a tower point
INTEGER(B4Ki)             :: TwrTop                     ! The index of top of the tower (it could be on the grid instead of at the end)

INTEGER(B4Ki)             :: IP
INTEGER(B2Ki),ALLOCATABLE :: TmpVarray(:)                ! This array holds the normalized velocities before being written to the binary file

INTEGER                   :: AllocStat

      ! Find the range of our velocity

   DO IC=1,3

         ! Initialize the Min/Max values

      VMin(IC) = V(1,1,IC)
      VMax(IC) = V(1,1,IC)

      DO II=1,NTot   ! Let's check all of the points
         DO IT=1,NumOutSteps  ! Use only the number of timesteps requested originally

            IF ( V(IT,II,IC) > VMax(IC) ) THEN

               VMax(IC) = V(IT,II,IC)

            ELSEIF ( V(IT,II,IC) < VMin(IC) ) THEN

               VMin(IC) = V(IT,II,IC)

            ENDIF

         ENDDO !IT
      ENDDO !II

   ENDDO !IC


      ! Calculate the scaling parameters for each component


   IF ( VMax(1) == VMin(1) ) THEN
      UScl = 1
   ELSE
      UScl = IntRng/REAL( VMax(1) - VMin(1) , SiKi )
   ENDIF

   IF ( VMax(2) == VMin(2) ) THEN
      VScl = 1
   ELSE
      VScl = IntRng/REAL( VMax(2) - VMin(2) , SiKi )
   ENDIF

   IF ( VMax(3) == VMin(3) ) THEN
      WScl = 1
   ELSE
      WScl = IntRng/REAL( VMax(3) - VMin(3) , SiKi )
   ENDIF


   UOff = IntMin - UScl*REAL( VMin(1)    , SiKi )
   VOff = IntMin - VScl*REAL( VMin(2)    , SiKi )
   WOff = IntMin - WScl*REAL( VMin(3)    , SiKi )


      ! Find the first tower point

   NumGrid  = NumGrid_Y*NumGrid_Z

   IF ( WrADTWR ) THEN

      TwrStart = NumGrid + 1

      IF ( ExtraHubPT ) THEN
         TwrStart = TwrStart + 1
      ENDIF

      IF ( ExtraTwrPt ) THEN
         TwrTop   = TwrStart
         TwrStart = TwrStart + 1
      ELSE
         TwrTop = INT(NumGrid_Y / 2) + 1      ! The top tower point is on the grid where Z = 1
      ENDIF

      NumTower = Ntot - TwrStart + 2

   ELSE

      NumTower = 0

   ENDIF


   LenDesc = LEN_TRIM( DescStr )             ! Length of the string that contains program name, version, date, and time

   CALL WrScr ( ' Generating AeroDyn binary time-series file "'//TRIM( RootName )//'.bts"' )


      ! Write the header


   WRITE (UAFFW)   INT(   7                , B2Ki )          ! TurbSim format

   WRITE (UAFFW)   INT( NumGrid_Z          , B4Ki )          ! the number of grid points vertically
   WRITE (UAFFW)   INT( NumGrid_Y          , B4Ki )          ! the number of grid points laterally
   WRITE (UAFFW)   INT( NumTower           , B4Ki )          ! the number of tower points
   WRITE (UAFFW)   INT( NumOutSteps        , B4Ki )          ! the number of time steps

   WRITE (UAFFW)  REAL( GridRes_Z          , SiKi )          ! grid spacing in vertical direction, in m
   WRITE (UAFFW)  REAL( GridRes_Y          , SiKi )          ! grid spacing in lateral direction, in m
   WRITE (UAFFW)  REAL( TimeStep           , SiKi )          ! grid spacing in delta time, in m/s
   WRITE (UAFFW)  REAL( UHub               , SiKi )          ! the mean wind speed in m/s at hub height
   WRITE (UAFFW)  REAL( HubHt              , SiKi )          ! the hub height, in m
   WRITE (UAFFW)  REAL( Z(1)               , SiKi )          ! the height of the grid bottom, in m

   WRITE (UAFFW)  REAL( UScl               , SiKi )          ! the U-component slope for scaling
   WRITE (UAFFW)  REAL( UOff               , SiKi )          ! the U-component offset for scaling
   WRITE (UAFFW)  REAL( VScl               , SiKi )          ! the V-component slope for scaling
   WRITE (UAFFW)  REAL( VOff               , SiKi )          ! the V-component offset for scaling
   WRITE (UAFFW)  REAL( WScl               , SiKi )          ! the W-component slope for scaling
   WRITE (UAFFW)  REAL( WOff               , SiKi )          ! the W-component offset for scaling

   WRITE (UAFFW)   INT( LenDesc            , B4Ki )          ! the number of characters in the string, max 200

   DO II=1,LenDesc

      WRITE (UAFFW)  INT( IACHAR( DescStr(II:II) ), B1Ki )   ! converted ASCII characters

   ENDDO

      ALLOCATE ( TmpVarray( 3*(NumGrid_Z*NumGrid_Y + NumTower) ) , STAT=AllocStat )

      IF ( AllocStat /= 0 )  THEN
         CALL TS_Abort ( 'Error allocating memory for temporary wind speed array.' )
      ENDIF

      ! Loop through time.

   DO IT=1,NumOutSteps  !Use only the number of timesteps requested originally

         ! Write out grid data in binary form. II = (IZ - 1)*NumGrid_Y + IY, IY varies most rapidly

      IP = 1

      DO II=1,NumGrid

         TmpVarray(IP)   =  NINT( Max( Min( REAL(UScl*V(IT,II,1) + UOff, SiKi), IntMax ),IntMin) , B2Ki )
         TmpVarray(IP+1) =  NINT( Max( Min( REAL(VScl*V(IT,II,2) + VOff, SiKi), IntMax ),IntMin) , B2Ki )
         TmpVarray(IP+2) =  NINT( Max( Min( REAL(WScl*V(IT,II,3) + Woff, SiKi), IntMax ),IntMin) , B2Ki )

         IP = IP + 3
      ENDDO ! II


      IF ( WrADTWR ) THEN

            ! Write out the tower data in binary form

            ! Value at the top of the tower (bottom of grid)
         TmpVarray(IP)   =  NINT( Max( Min( REAL(UScl*V(IT,TwrTop,1) + UOff, SiKi), IntMax ),IntMin) , B2Ki )
         TmpVarray(IP+1) =  NINT( Max( Min( REAL(VScl*V(IT,TwrTop,2) + VOff, SiKi), IntMax ),IntMin) , B2Ki )
         TmpVarray(IP+2) =  NINT( Max( Min( REAL(WScl*V(IT,TwrTop,3) + Woff, SiKi), IntMax ),IntMin) , B2Ki )

         IP = IP + 3
         DO II=TwrStart,NTot
                ! Values of tower data
            TmpVarray(IP)   =  NINT( Max( Min( REAL(UScl*V(IT,II,1) + UOff, SiKi), IntMax ),IntMin) , B2Ki )
            TmpVarray(IP+1) =  NINT( Max( Min( REAL(VScl*V(IT,II,2) + VOff, SiKi), IntMax ),IntMin) , B2Ki )
            TmpVarray(IP+2) =  NINT( Max( Min( REAL(WScl*V(IT,II,3) + Woff, SiKi), IntMax ),IntMin) , B2Ki )

            IP = IP + 3
         ENDDO ! II

      ENDIF

      WRITE ( UAFFW ) TmpVarray(:)
   ENDDO ! IT

   CLOSE ( UAFFW )

   IF ( ALLOCATED( TmpVarray ) ) DEALLOCATE( TmpVarray )


END SUBROUTINE WrBinTURBSIM
!=======================================================================
SUBROUTINE WrFormattedFF(HubIndx)

USE                     TSMods

IMPLICIT                NONE

REAL(ReKi), ALLOCATABLE      :: ZRow       (:)                           ! The horizontal locations of the grid points (NumGrid_Y) at each height.

INTEGER                      :: AllocStat
INTEGER, INTENT(IN)          :: HubIndx
INTEGER                      :: II
INTEGER                      :: IT
INTEGER                      :: IVec
INTEGER                      :: IY
INTEGER                      :: IZ

CHARACTER(1)                 :: Comp (3) = (/ 'u', 'v', 'w' /)           ! The names of the wind components
CHARACTER(200)               :: FormStr3                                 ! String used to store format specifiers.
CHARACTER(200)               :: FormStr4                                 ! String used to store format specifiers.
CHARACTER(200)               :: FormStr5                                 ! String used to store format specifiers.
CHARACTER(200)               :: FormStr6                                 ! String used to store format specifiers.



   FormStr  = "( / 'This full-field turbulence file was generated by ' , A , A , ' on ' , A , ' at ' , A , '.' / )"
   FormStr1 = "( ' | ', A,'-comp |  Y  x  Z  | Grid Resolution (Y x Z) | Time-step | Hub Elev | Mean U |')"
   FormStr2 = "(I14,I6,F11.3,F11.3,F15.3,F11.2,F10.2)"
   FormStr3 = "(/,' Z Coordinates (m):')"
   FormStr4 = "(/,' Y Coordinates (m):')"
   FormStr5 = "(1X,98(F8.3),:)"
   FormStr6 = "(/,1X,2(F8.3))"

      ! Allocate the array of wind speeds.

   ALLOCATE ( ZRow(NumGrid_Y) , STAT=AllocStat )

   IF ( AllocStat /= 0 )  THEN
      CALL TS_Abort ( 'Error allocating memory for array of wind speeds.' )
   ENDIF


   DO IVec=1,3

      CALL WrScr ( ' Generating full-field formatted file "'//TRIM(RootName)//'.'//Comp(IVec)//'".' )

      CALL OpenFOutFile ( UFFF, TRIM( RootName )//'.'//Comp(IVec) )


         ! Create file header.

      WRITE (UFFF,FormStr )  TRIM(ProgName), TRIM( ProgVer ), CurDate(), CurTime()

      WRITE (UFFF,FormStr1)  Comp(IVec)

      WRITE (UFFF,FormStr2)  NumGrid_Y, NumGrid_Z, GridRes_Y, GridRes_Z, TimeStep, HubHt, UHub
      WRITE (UFFF,FormStr3)
      WRITE (UFFF,FormStr5)  ( Z(IZ)-HubHt, IZ=1,NumGrid_Z )
      WRITE (UFFF,FormStr4)
      WRITE (UFFF,FormStr5)  ( Y(IY), IY=1,NumGrid_Y )

         ! Write out elapsed time & hub-level value before component grid.

      DO IT=1,NumOutSteps

         WRITE(UFFF,FormStr6)  TimeStep*( IT - 1 ), V(IT,HubIndx,IVec)

         DO IZ=1,NumGrid_Z  ! From the top to the bottom

            II = ( NumGrid_Z - IZ )*NumGrid_Y

            DO IY=1,NumGrid_Y  ! From the left to the right
               ZRow(IY) = V(IT,II+IY,IVec)
            ENDDO ! IY

            WRITE (UFFF,FormStr5)  ( ZRow(IY), IY=1,NumGrid_Y )

         ENDDO ! IZ

      ENDDO ! IT

      CLOSE ( UFFF )

   ENDDO ! IVec

      ! Deallocate the array of wind speeds.

   IF ( ALLOCATED( ZRow ) )  DEALLOCATE( ZRow )

END SUBROUTINE WrFormattedFF
!=======================================================================
SUBROUTINE WriteEvents( UnOut, UnIn, TScale )

    ! This subroutine writes the events as calculated in CalcEvents.

USE                     TSMods

IMPLICIT                NONE

      ! Passed Variables

REAL(ReKi), INTENT(IN)  :: TScale                 ! Time scaling factor
INTEGER,    INTENT(IN)  :: UnIn                   ! I/O Unit for input file
INTEGER,    INTENT(IN)  :: UnOut                  ! I/O Unit for output file


      ! Local Variables

REAL(ReKi)              :: CurrentTime = 0.0      ! the current time (in seconds)
REAL(ReKi)              :: CTTime                 ! Time from beginning of event file
REAL(ReKi)              :: deltaTime = 0.0        ! difference between two time steps in the event files

INTEGER                 :: FileNum                ! File Number in the event file
INTEGER                 :: IE                     ! Loop counter for event number
INTEGER                 :: IStat                  ! Status of file read
INTEGER                 :: IT                     ! Loop counter for time step


CHARACTER(200)          :: InpFile                ! Name of the input file
TYPE (Event), POINTER   :: PtrCurr  => NULL()     ! Pointer to the current event


IF (DEBUG_OUT) THEN
   WRITE (UD,'(/,A)' ) 'Computed Coherent Events'
   WRITE (UD,*) 'Event#     Start Time        End Time'
ENDIF


   PtrCurr => PtrHead


   DO IE = 1,NumCTEvents

      IF ( .NOT. ASSOCIATED ( PtrCurr ) ) EXIT     ! This shouldn't be necessary, given the way we created the list


      IF ( .NOT. PtrCurr%Connect2Prev ) THEN

         IF ( CurrentTime < PtrCurr%TStart ) THEN

            WRITE ( UnOut, '(G14.7,1x,I5.5)') CurrentTime, 0     ! Print end of previous event

            NumCTt = NumCTt - 1  ! Let's make sure the right number of points have been written to the file.

            IF ( CurrentTime < PtrCurr%TStart - PtrCurr%delt ) THEN  !This assumes a ramp of 1 delta t for each structure....

               WRITE ( UnOut, '(G14.7,1x,I5.5)') MAX(PtrCurr%TStart - PtrCurr%delt, REAL(0.0, ReKi) ), 0
               NumCTt = NumCTt - 1

            ENDIF

         ENDIF

      ENDIF  ! NOT Connect2Prev


      WRITE ( InpFile, '(I5.5)' ) EventName( PtrCurr%EventNum )
      InpFile = TRIM( CTEventPath )//'\Event'//TRIM( InpFile)//'.dat'

      CALL OpenFInpFile( UnIn, InpFile )


      DO IT = 1,EventTS( PtrCurr%EventNum )

         READ  ( UnIn, *, IOSTAT=IStat ) FileNum, CTTime, deltaTime

         IF (IStat /= 0) THEN
            CALL TS_Abort( 'Error reading event file'//TRIM( InpFile ) )
         ENDIF

         CurrentTime = PtrCurr%TStart + CTTime*TScale

         WRITE ( UnOut, '(G14.7,1x,I5.5)') CurrentTime, FileNum
         NumCTt = NumCTt - 1

      ENDDO    ! IT: Event timestep


      CLOSE ( UnIn )


         ! Add one (delta time) space between events

      CurrentTime = CurrentTime + deltaTime*TScale

IF (DEBUG_OUT) THEN
   WRITE ( UD,'(I6, 2(2x,F14.5))' ) PtrCurr%EventNum, PtrCurr%TStart, CurrentTime
ENDIF

      PtrCurr => PtrCurr%Next
      
!bjj deallocate the linked list!!!!

   ENDDO !IE: number of events

   WRITE ( UnOut, '(G14.7,1x,I5.5)') CurrentTime, 0     !Add the last line
   NumCTt = NumCTt - 1

      ! Let's append zero lines at the end of the file if we haven't output NumCTt lines, yet.
      ! We've subtracted from NumCTt every time we wrote a line so now the number in NumCTt is
      ! how many lines short we are.

   IF ( deltaTime > 0 ) THEN
      deltaTime = deltaTime*TScale
   ELSE
      deltaTime = 0.5
   ENDIF

   DO IE = 1, NumCTt  ! Write zeros at the end if we happened to insert an event that overwrote one of our zero lines
      CurrentTime = CurrentTime + deltaTime
      WRITE ( UnOut, '(G14.7,1x,I5.5)') CurrentTime, 0
   ENDDO

END SUBROUTINE WriteEvents
!=======================================================================
END MODULE TSSubs
