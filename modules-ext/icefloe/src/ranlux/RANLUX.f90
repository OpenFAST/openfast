! Borrowed from the v1.50 distribution of Turbsim
 Module Ran_Lux_Mod
! Subtract-and-borrow random number generator proposed by Marsaglia and Zaman, implemented by F. James with the name RCARRY in 1991,
! and later improved by Martin Luescher in 1993 to produce "Luxury Pseudorandom Numbers". Fortran 77 coded by F. James, 1993.
! Converted to Fortran 90 [and Lahey Elf90 subset] by Loren Meissner, 1995.
!
! References: M. Luscher, Computer Physics Communications 79 (1994) 100; F. James, Computer Physics Communications 79 (1994) 111
!
! LUXURY LEVELS. -- The available luxury levels are:
! level 0 (p = 24) : equivalent to the original RCARRY of Marsaglia and Zaman, very long period, but fails many tests.
! level 1 (p = 48) : considerable improvement in quality over level 0, now passes the gap test, but still fails spectral test.
! level 2 (p = 97) : passes all known tests, but theoretically still defective.
! level 3 (p = 223) : DEFAULT VALUE. Any theoretically possible correlations have very small chance of being observed.
! level 4 (p = 389) : highest possible luxury, all 24 bits chaotic.
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
! Calling sequences for RanLux:
! call RanLux (RVec)
!   Returns a vector RVec of Len(RVec) 32-bit random floating point numbers X, such that 0.0 < X < 1.0 .
! call RLuxGo (Lux, Int, K1, K2)
!   Initializes the generator from one 32-bit integer INT and sets Luxury Level LUX which is integer between zero and MaxLev, or
!   if Lux > 24, it sets p = Lux directly. K1 and K2 should be set to zero unless restarting at a break point given by output of
!   RLuxAt (see RLuxAt).
! call RLuxAt (Lux, Int, K1, K2)
!   Gets the values of four integers which can be used to restart the RanLux generator at the current point by calling RLuxGo.
!   K1 and K2 specify how many numbers were generated since the initialization with Lux and Int. The restarting skips over
!   K1 + K2 * 1E9 numbers, so it can be long. A more efficient but less convenient way of restarting is by:
!     call RLuxIn (ISVec) ! Restart the generator from vector ISVec of 25 32-bit integers (see RLXUt)
!     call RLuxUt (ISVec) ! Output the current values of the 25 32-bit integer Seeds, to be used for restarting.
! The array argument to RLuxIn or RLuxUt must be dimensioned 25 in the calling program
! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
!                                     default
!       Luxury Level    0     1     2   * 3 *   4
!                       0    24    73   199   365
!   corresponds to p = 24    48    97   223   389
!           time factor 1     2     3     6    10 on slow workstation
!                       1   1.5     2     3     5 on fast mainframe
!
! NotYet is .TRUE. if no initialization has been performed yet.
!Start bjj:  We want to write to the screen instead of "print *"
!   use     NWTC_IO
!End bjj:
  use precision
  implicit none

  integer, parameter :: NSeeds = 25, MaxLev = 4, LxDflt = 3
  real(ReKi), parameter :: TwoP12 = 4096.0
  integer, parameter :: IGiga = 1000000000, JSDFlt = 314159265, ITwo24 = 2 ** 24, ICons = 2147483563
  integer :: I_ranLux
  integer, parameter :: Next(NSeeds - 1) = (/ NSeeds - 1, (I_ranLux, I_ranLux = 1, NSeeds - 2) /)  ! Table look-up (faster than Mod function).

  integer :: I24 = 24, J24 = 10, In24 = 0, Kount = 0, LuxLev = LxDflt, MKount = 0 ! Initialized variables are automatically saved.
  integer, dimension(0: MaxLev) :: NDSkip = (/ 0, 24, 73, 199, 365 /) ! Initialized variables are automatically saved.
  integer, save :: NSkip, InSeed
  real(ReKi) :: Carry = 0.0 ! Initialized variables are automatically saved.
  real(ReKi), save :: Seeds(NSeeds - 1), TwoM24, TwoM12
  logical, save :: NotYet = .TRUE.

  real(ReKi) :: Uni

!bjj
  character(300) :: RanLux_str
!bjj

  private :: RCarry

contains
!============================================================================
  subroutine RanLux (RVec)
    ! Default Initialization by Multiplicative Congruential
    real(ReKi), intent(out) :: RVec(:)
    integer :: ISeeds(NSeeds - 1), I, IVec, JSeed, K, LEnv, LP

    real(ReKi) :: tmpTwoM24, tmpTwoM24Seed

! start subroutine RanLux
    LEnv = Size (RVec)
    if (NotYet) then
      NotYet = .FALSE.
      JSeed = JSDFlt
      InSeed = JSeed
!begin bjj
!      print *, " RanLux default initialization: ", JSeed
!      write( RanLux_str, '(I12)' ) JSeed
!      CALL WrScr( " RanLux default initialization: "//TRIM( ADJUSTL( RanLux_str ) ) )
!end bjj
      LuxLev = LxDflt
      NSkip = NDSkip(LuxLev)
      LP = NSkip + NSeeds - 1
      In24 = 0
      Kount = 0
      MKount = 0
!begin bjj
!      print *, " RanLux default luxury level = ", LuxLev, " p = ", LP
!      write( RanLux_str, '(A,I5,A,I12)' ) " RanLux default luxury level = ", LuxLev, " p = ", LP
!      CALL WrScr( TRIM( RanLux_str ) )
!end bjj

      TwoM24 = 1.0
      do I = 1, NSeeds - 1
        TwoM24 = TwoM24 * 0.5
        K = JSeed / 53668
        JSeed = 40014 * (JSeed - K * 53668) - K * 12211
        if (JSeed < 0) JSeed = JSeed + ICons
        ISeeds(I) = Mod (JSeed, ITwo24)
      end do
      TwoM12 = TwoM24 * 4096.0
      Seeds = Real (ISeeds) * TwoM24
      I24 = NSeeds - 1
      J24 = 10
      Carry = Merge (TwoM24, 0.0_ReKi, Seeds(NSeeds - 1) == 0.0)
    end if

      !bjj added to speed up later calculations (b/c I had to fix the "where" statement)
    tmpTwoM24Seed = TwoM24 * Seeds(J24)
    tmpTwoM24 = TwoM24 * TwoM24
      !bjj end of modifications

    do IVec = 1, LEnv
      RVec(IVec) = RCarry (1)
    ! Skipping to Luxury. As proposed by Martin Luscher.
      In24 = In24 + 1
      if (In24 == NSeeds - 1) then
        In24 = 0
        Kount = Kount + NSkip
        Uni = RCarry (NSkip)
      end if

         !bjj modified code to eliminate code crashing on "where" statement for large arrays
      ! "Pad" small numbers (with less than 12 "significant" bits) and eliminate zero values (in case someone takes a logarithm)
      if ( RVec(IVec) < TwoM12 ) RVec(IVec) = RVec(IVec) + tmpTwoM24Seed
      if ( Rvec(IVec) == 0.0 )  RVec(IVec) = tmpTwoM24
         !bjj end of modifications

   end do
      !bjj removed to eliminate crashing in SNwind
   ! "Pad" small numbers (with less than 12 "significant" bits) and eliminate zero values (in case someone takes a logarithm)
    !where (RVec < TwoM12) RVec = RVec + TwoM24 * Seeds(J24)
    !where (Rvec == 0.0) RVec = TwoM24 * TwoM24
      !bjj end of modifications

    Kount = Kount + LEnv
    if (Kount >= IGiga) then
      MKount = MKount + 1
      Kount = Kount - IGiga
    end if
    return
  end subroutine RanLux
!============================================================================
! Input and float integer Seeds from previous run
  subroutine RLuxIn (ISDext)
    integer, intent(in) :: ISDext(:)
    integer :: I, ISD
! start subroutine RLuxIn
    if (Size(ISDext) /= NSeeds) then
!begin bjj
!      print *, " Array size for RLuxIn must be ", NSeeds
!      write( RanLux_str, '(I5)' ) NSeeds
!      CALL WrScr( " Array size for RLuxIn must be "//TRIM( ADJUSTL(RanLux_str) ) )
!end bjj

      return
    end if
    ! The following IF block added by Phillip Helbig, based on conversation with Fred James;
    ! an equivalent correction has been published by James.
    if (NotYet) then
!begin bjj
!      print *, " Proper results only with initialisation from 25 integers obtained with RLuxUt"
!      CALL WrScr( " Proper results only with initialisation from 25 integers obtained with RLuxUt" )
!end bjj
      NotYet = .FALSE.
    end if
    TwoM24 = 1.0
    do I = 1, NSeeds - 1
      TwoM24 = TwoM24 * 0.5
    end do
    TwoM12 = TwoM24 * 4096.0
!Start bjj
!    print *, " Full initialization of RanLux with 25 integers:"
!    print *, ISDext
!    CALL WrScr ( " Full initialization of RanLux with 25 integers:" )
!    write( RanLux_str, '(25(I11,1x))' ) ISDext
!    CALL WrScr ( TRIM( RanLux_str ) )
!End bjj
    Seeds = Real (ISDext(: NSeeds - 1)) * TwoM24
    Carry = 0.0
    if (ISDext(NSeeds) < 0) Carry = TwoM24
    ISD = Abs (ISDext(NSeeds))
    I24 = Mod (ISD, 100)
    ISD = ISD / 100
    J24 = Mod (ISD, 100)
    ISD = ISD / 100
    In24 = Mod (ISD, 100)
    ISD = ISD / 100
    LuxLev = ISD

!start bjj
   write( RanLux_str, "(I5)" ) LuxLev
!end bjj

    if (LuxLev <= MaxLev) then
      NSkip = NDSkip(LuxLev)
!start bjj
!      print *, " RanLux luxury level set by RLuxIn to: ", LuxLev\
!       CALL WrScr( " RanLux luxury level set by RLuxIn to: "//TRIM(ADJUSTL(RanLux_str) ))
!end bjj
    else if (LuxLev >= NSeeds - 1) then
      NSkip = LuxLev - NSeeds + 1
!start bjj
!      print *, " RanLux p-value set by RLuxIn to:", LuxLev
!      CALL WrScr( " RanLux p-value set by RLuxIn to: "//TRIM(ADJUSTL(RanLux_str) ))
!end bjj
    else
      NSkip = NDSkip(MaxLev)
!start bjj
!      print *, " RanLux illegal luxury RLuxIn: ", LuxLev
!      CALL WrScr( " RanLux illegal luxury RLuxIn: "//TRIM(ADJUSTL(RanLux_str) ))
!end bjj
      LuxLev = MaxLev
    end if
    InSeed = - 1
    return
  end subroutine RLuxIn
!============================================================================
! Ouput Seeds as integers
  subroutine RLuxUt (ISDext)
    integer, intent(out) :: ISDext(:)
! start subroutine RLuxUt
    if (Size(ISDext) /= NSeeds) then
      ISDext = 0
!start bjj
!      print *, " Array size for RLuxUt must be ", NSeeds
!      write( RanLux_str, '(I20)' ) NSeeds
!      CALL WrScr( " Array size for RLuxUt must be "//TRIM( ADJUSTL(RanLux_str )))
!end bjj
      return
    end if
    ISDext(: NSeeds - 1) = Int (Seeds * TwoP12 * TwoP12)
    ISDext(NSeeds) = Merge (-ISDext(NSeeds), I24 + 100 * J24 + 10000 * In24 + 1000000 * LuxLev, Carry > 0.0)
    return
  end subroutine RLuxUt
!============================================================================
! Output the "convenient" restart point
  subroutine RLuxAt (LOut, InOut, K1, K2)
    integer, intent(out) :: LOut, InOut, K1, K2
! start subroutine RLuxAt
    LOut = LuxLev
    InOut = InSeed
    K1 = Kount
    K2 = MKount
    return
  end subroutine RLuxAt
!============================================================================
! Initialize from one or three integers
  subroutine RLuxGo (Lux, Int, K1, K2)
    integer, intent(in) :: Lux, Int, K1, K2
    integer :: ISeeds(NSeeds - 1), ILx, I, IOuter, IZip, IZip2, JSeed, K
! start subroutine RLuxGo
    if (Lux < 0) then
      LuxLev = LxDflt
    else if (Lux <= MaxLev) then
      LuxLev = Lux
    else if (Lux < NSeeds - 1 .or. Lux > 2000) then
      LuxLev = MaxLev
!start bjj
!      print *, " RanLux illegal luxury level in RLuxGo: ", Lux
!       write( RanLux_str, '(I20)' ) Lux
!       Call WrScr( " RanLux illegal luxury level in RLuxGo: "//TRIM( ADJUSTL(RanLux_str ) ))
!end bjj
    else
      LuxLev = Lux
      do ILx = 0, MaxLev
        if (Lux == NDSkip(ILx) + NSeeds - 1) then
          LuxLev = ILx
        end if
      end do
    end if
    if (LuxLev <= MaxLev) then
      NSkip = NDSkip(LuxLev)
!start bjj
!      print *, " RanLux luxury level set by RLuxGo :", LuxLev, " p = ", NSkip + NSeeds - 1
!      write (RanLux_str, '(A,I5)') " RanLux luxury level set by RLuxGo :", LuxLev
!      write (RanLux_str, '(A,I12)') TRIM(RanLux_str)//" p = ", NSkip + NSeeds - 1
!      CALL WrScr( TRIM(RanLux_str) )
!end bjj
    else
      NSkip = LuxLev - 24
!start bjj
!      print *, " RanLux p-value set by RLuxGo to:", LuxLev
!      write( RanLux_str, '(I20)' ) LuxLev
!      CALL WrScr( " RanLux p-value set by RLuxGo to: "//TRIM( ADJUSTL(RanLux_str ) ))
!end bjj
    end if
    In24 = 0
    if (Int < 0) then
!start bjj
!      print *, " Illegal initialization by RLuxGo, negative input seed"
!       CALL WrScr( " Illegal initialization by RLuxGo, negative input seed" )
!end bjj
   else if (Int > 0) then
      JSeed = Int
!start bjj
!      print *, " RanLux initialized by RLuxGo from Seeds", JSeed, K1, K2
!      write( RanLux_str, '(3(I12))' ) JSeed, K1, K2
!      CALL WrScr( " RanLux initialized by RLuxGo from Seeds"//TRIM( RanLux_str ) )
!end bjj
    else
      JSeed = JSDFlt
!start bjj
!      print *, " RanLux initialized by RLuxGo from default seed"
!      CALL WrScr( " RanLux initialized by RLuxGo from default seed" )
!end bjj
    end if
    InSeed = JSeed
    NotYet = .FALSE.
    TwoM24 = 1.0
    do I = 1, NSeeds - 1
      TwoM24 = TwoM24 * 0.5
      K = JSeed / 53668
      JSeed = 40014 * (JSeed - K * 53668) - K * 12211
      if (JSeed < 0) JSeed = JSeed + ICons
      ISeeds(I) = Mod (JSeed, ITwo24)
    end do
    TwoM12 = TwoM24 * 4096.0
    Seeds = Real (ISeeds) * TwoM24
    I24 = NSeeds - 1
    J24 = 10
    Carry = Merge (TwoM24, 0.0_ReKi, Seeds(NSeeds - 1) == 0.0)

    ! If restarting at a break point, skip K1 + IGIGA * K2
    ! Note that this is the number of numbers delivered to the user PLUS the number skipped (if Luxury > 0) .
    Kount = Abs (K1)
    MKount = Abs (K2)
    if (Kount + MKount /= 0) then
      do IOuter = 1, MKount + 1
        Uni = RCarry (Merge (Kount, IGiga, IOuter == MKount + 1))
      end do
      ! Get the right value of IN24 by direct calculation
      In24 = Mod (Kount, NSkip + NSeeds - 1)
      if (MKount > 0) then
        IZip = Mod (IGiga, NSkip + NSeeds - 1)
        IZip2 = MKount * IZip + In24
        In24 = Mod (IZip2, NSkip + NSeeds - 1)
      end if
      ! Now IN24 had better be between zero and 23 inclusive
      if ((In24 < 1) .or. (In24 >= NSeeds - 1)) then
!start bjj
!        print *, " Error in restarting with RLuxGo: the values", Int, K1, K2, " cannot occur at luxury level", LuxLev
!        write( RanLux_str, '(A,3(I12),A,I5)' ) " Error in restarting with RLuxGo: the values ", Int, K1, K2, &
!                                            " cannot occur at luxury level ", LuxLev
!        CALL WrScr( TRIM(RanLux_str ) )
!end bjj
        In24 = 0
      end if
    end if
    return
  end subroutine RLuxGo
!============================================================================
  function RCarry (N) result (Uni)  ! Private (in module); generates a sequence of N uniform random numbers; returns the last one.
    real(ReKi) :: Uni
    integer, intent(in) :: N
    integer :: Many
! start function RCarry
    do Many = 1, N
    ! The Generator proper: "Subtract-with-borrow", as proposed by Marsaglia and Zaman, Florida State University, March, 1989
      Uni = Seeds(J24) - Seeds(I24) - Carry
      if (Uni < 0.0) then
        Uni = Uni + 1.0
        Carry = TwoM24
      else
        Carry = 0.0
      end if
      Seeds(I24) = Uni
      I24 = Next(I24)
      J24 = Next(J24)
    end do
    return
  end function RCarry
!============================================================================
 end Module Ran_Lux_Mod
