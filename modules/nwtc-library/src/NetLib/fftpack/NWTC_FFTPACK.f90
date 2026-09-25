   ! NOTE: This MODULE isused in HydroDyn and TurbSim.
   ! BJJ: 02/22/2008: Updated to work with NWTC_Library v1.01.09
   !      all Abort() functions changed to ProgAbort()
   ! BJJ: 12/03/2010: Updated to add optional ErrStat return values
   !                  instead of aborting on errors (note not all changes documented)
   ! BJJ: 12/20/2010: Updated to add defined type and remove global data variables
   !                  Also updated to check that transform has been initialized for the
   !                    correct type (to avoid having wSave too small)
   ! ADP: 07/28/2014: Added in the complex FFT routines from fftpack v. 4.1
   ! ADP: 08/15/2026: upgraded from fftpack v4.1 to v5.1 and added interfaces for 2D fft
!=======================================================================
MODULE NWTC_FFTPACK
!-----------------------------------------------------------------------
! DESCRIPTION OF THE INVERSE FOURIER TRANSFORM ROUTINE:
!
! Given an array, H, of N complex numbers, calculate an array, h, of N real
! numbers:
!     h(J) = the sum from K=1,...,N of [ H(K)*exp(i*(J-1)*(K-1)*2*pi/N) ]
!     for J = 1,...,N
!
! where:
!     i=sqrt(-1)
!
! In order for h to be real, the real components of H must be an even function
! of frequency and the imaginary components of H must be an odd function of
! frequency.  Thus, only the first N/2 + 1 values of H are unique.  (The first
! N/2 + 1 are the positive frequencies including zero; the last N/2 - 1 values
! are the negative frequencies.)
!
! We accomplish this by defining a real array, R, where:
!     R(1) = REAL( H(1) ),
!     R(2) = REAL( H(2) ), R(3) = IMAG( H(2) ),
!     R(4) = REAL( H(3) ), R(5) = IMAG( H(3) ),...
!     R(N) = REAL( H(N/2 + 1) ).
! Note that the values of IMAG( H(1) ) and IMAG( H(N/2 + 1) ) must be zero for
! the result to be real, else the routine will abort.
!
! We return the array, r = h, of real numbers as calculated by:
!     h(J) = r(J) = R(1) + (-1)**(J-1)*R(N)
!                 + the sum from K=2,...,N/2 of
!                   [  2*R(2*K-2)*COS((K-1)*(J-1)*2*PI/N)
!                     -2*R(2*K-1)*SIN((K-1)*(J-1)*2*PI/N) ]
!     for J = 1,...,N, where N is an even number
!
! The routine is most effecient when N is a product of small primes.
!
! If the Normalization flag is set to "TRUE" in the initialization, we
! normalize the result by 1/N.
!------------------------------------------------------------------------
! DESCRIPTION OF THE COSINE TRANSFORM ROUTINE:
!
! Given an array, X, of N real numbers, calculate an array, x, of N real
! numbers:
!     x(J) = X(1) + (-1)**(J-1)*X(N)
!          + the sum from K=2,...,N-1 of [ 2*X(K)*COS((K-1)*(J-1)*PI/(N-1)) ]
!     for J = 1,...,N, where N is an odd number
!
! The routine is most effecient when N-1 is a product of small primes.
!
! If the Normalization flag is set to "TRUE" in the initialization, we
! normalize the result by 1/(N-1).
!------------------------------------------------------------------------
! DESCRIPTION OF THE SINE TRANSFORM ROUTINE:
!
! Given an array, X, of N real numbers, calculate an array, x, of N real
! numbers:
!     x(1) = X(1) = 0
!     x(J) = the sum from K=2,...,N-1 of [ 2*X(K)*SIN((K-1)*(J-1)*PI/(N-1)) ]
!     for J = 2,...,N-1, where N is an odd number
!     x(N) = X(N) = 0
!
! Note that the values of X(1) and X(N) must be zero, else the routine will
! abort.
!
! The routine is most effecient when N-1 is a product of small primes.
!
! If the Normalization flag is set to "TRUE" in the initialization, we
! normalize the result by 1/(N-1).
!------------------------------------------------------------------------

! We need the Precision module and the Abort() and
! Int2LStr() functions from the NWTC_IO module.

   USE                                    NWTC_Library


   IMPLICIT                               NONE

   INTEGER, PARAMETER, PRIVATE         :: Undef_trans   = -1            ! transformation type is undefined
   INTEGER, PARAMETER, PRIVATE         :: COS_trans     = 1             ! COSINE transformation
   INTEGER, PARAMETER, PRIVATE         :: Fourier_trans = 2             ! FAST FOURIER transformation
   INTEGER, PARAMETER, PRIVATE         :: SIN_trans     = 3             ! SINE transformation
   INTEGER, PARAMETER, PRIVATE         :: Fourier2D_trans = 4           ! 2D FAST FOURIER transformation
   INTEGER, PARAMETER, PRIVATE         :: CFourier2D_trans = 5          ! 2D complex FAST FOURIER transformation

   TYPE, PUBLIC :: FFT_DataType
      PRIVATE
      REAL(SiKi)                       :: InvN          = 0.0_SiKi      ! Normalization constant
      REAL(SiKi), ALLOCATABLE          :: wSave(:)                      ! Trig/factor table for FFTPACK 5.1
      INTEGER                          :: N             = -1            ! Number of steps
      INTEGER                          :: LenWork       = 0             ! Required scratch workspace size
      LOGICAL                          :: Normalize     = .FALSE.       ! Whether or not to normalize
      INTEGER                          :: TransformType = Undef_trans   ! the type of transfer function this is for
   END TYPE FFT_DataType      

   TYPE, PUBLIC :: FFT2D_DataType
      PRIVATE
      REAL(SiKi)                       :: InvN          = 0.0_SiKi      ! Normalization constant = 1/(L*M)
      REAL(SiKi), ALLOCATABLE          :: wSave(:)                      ! Trig/factor table for FFTPACK 5.1
      INTEGER                          :: L             = -1            ! Number of rows
      INTEGER                          :: M             = -1            ! Number of columns
      INTEGER                          :: LenWork       = 0             ! Required scratch workspace size
      LOGICAL                          :: Normalize     = .FALSE.       ! Whether or not to normalize
      INTEGER                          :: TransformType = Undef_trans   ! the type of transfer function this is for
   END TYPE FFT2D_DataType


!------------------------------------------------------------------------
CONTAINS

   SUBROUTINE ApplyCOST( TRH, FFT_Data, ErrStat )
   
         ! Perform cosine transform.

      IMPLICIT                         NONE

      REAL(SiKi), INTENT(INOUT)     :: TRH(:)
      TYPE(FFT_DataType), INTENT(IN):: FFT_Data             ! the handle to this instance of the FFT Module
      
      INTEGER, INTENT(OUT), OPTIONAL:: ErrStat
      
      REAL(SiKi)                    :: wWork(FFT_Data%LenWork)
      LOGICAL                       :: TrapErrors
      INTEGER                       :: IER
      

         
      IF ( PRESENT(ErrStat) ) THEN
         TrapErrors = .TRUE.
         ErrStat = ErrID_None
      ELSE
         TrapErrors = .FALSE.         
      END IF


        ! Make sure the array isn't too small

      IF ( SIZE(TRH) < FFT_Data%N )  THEN
          CALL ProgAbort( 'Error in call to cosine transform.  Array size is not large enough.', TrapErrors )
          ErrStat = ErrID_Fatal         ! The code can't get here unless PRESENT(ErrStat)
          RETURN
      END IF
      
      IF ( FFT_Data%TransformType /= COS_trans ) THEN
          CALL ProgAbort( 'Error in call to cosine transform. FFT_Data not initialized for cosine transform.', TrapErrors )
          ErrStat = ErrID_Fatal
          RETURN
      END IF      
      

        ! Perform the cosine transform with a FFTPACK 5.1 routine.
        ! COST1F computes a normalized DCT-I: endpoints get 1/(2*(N-1)), interior 1/(N-1).
        ! The old COST formula was un-normalized, so post-scale to recover old behavior.

      CALL COST1F(FFT_Data%N, 1, TRH, SIZE(TRH), FFT_Data%wSave, SIZE(FFT_Data%wSave), &
                  wWork, FFT_Data%LenWork, IER)
      IF (IER /= 0) THEN
          CALL ProgAbort( 'Error in cosine transform (COST1F).', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      TRH(1) = 2.0_SiKi * REAL(FFT_Data%N - 1, SiKi) * TRH(1)
      TRH(2:FFT_Data%N-1) = REAL(FFT_Data%N - 1, SiKi) * TRH(2:FFT_Data%N-1)
      TRH(FFT_Data%N) = 2.0_SiKi * REAL(FFT_Data%N - 1, SiKi) * TRH(FFT_Data%N)

      IF (FFT_Data%Normalize) THEN
          TRH(1:FFT_Data%N) = FFT_Data%InvN * TRH(1:FFT_Data%N)
      ENDIF

   END SUBROUTINE ApplyCOST
!------------------------------------------------------------------------
   SUBROUTINE ApplyCFFT( TRH_complex_return, TRH_complex, FFT_Data, ErrStat )
         ! Perform Backward complex FFT: given TRH_complex, an array of complex numbers,
         ! return an array TRH_complex_return, of complex numbers

         ! TRH_complex is of size FFT_Data%N/2 and represents the complex amplitude in frequency
         !  space of only the positive frequencies.  This is padded with zeros for the upper half
         !  the frequency domain.

         ! TRH_complex_return is of size FFT_Data%N and represents the complex amplitude in the
         !  time domain.

      IMPLICIT                         NONE

      COMPLEX(SiKi), INTENT(OUT)    :: TRH_complex_return(:)
      COMPLEX(SiKi), INTENT(IN)     :: TRH_complex(:)
      TYPE(FFT_DataType), INTENT(IN):: FFT_Data             ! the handle to this instance of the FFT Module
      INTEGER, INTENT(OUT), OPTIONAL:: ErrStat

      REAL(SiKi)                    :: wWork(FFT_Data%LenWork)
      INTEGER                       :: I
      INTEGER                       :: IER
      LOGICAL                       :: TrapErrors

         
      IF ( PRESENT(ErrStat) ) THEN
         TrapErrors = .TRUE.
         ErrStat = ErrID_None
      ELSE
         TrapErrors = .FALSE.         
      END IF


        ! Make sure the arrays aren't too small
      IF ( ( SIZE(TRH_complex_return) < FFT_Data%N ) .OR. ( SIZE(TRH_complex) < ( FFT_Data%N/2 + 1 ) ) )  THEN
          CALL ProgAbort( 'Error in call to FFT.  Array size is not large enough.', TrapErrors )
          ErrStat = ErrID_Fatal         ! The code can't get here unless PRESENT(ErrStat)
          RETURN
      END IF
      
      IF ( FFT_Data%TransformType /= Fourier_trans ) THEN
          CALL ProgAbort( 'Error in call to FFT. FFT_Data not initialized for Fourier transform.', TrapErrors )
          ErrStat = ErrID_Fatal
          RETURN
      END IF      


        ! Make sure that the imaginary components at the zeroeth and largest
        ! positive frequency are zero, else abort.

      IF ( .NOT. EqualRealNos( 0.0_SiKi, AIMAG( TRH_complex(1    ) ) ) ) THEN
          CALL ProgAbort( 'Error in call to FFT.  The imaginary component at the zeroeth frequency must be zero.', TrapErrors )
          ErrStat = ErrID_Fatal         ! The code can't get here unless PRESENT(ErrStat)
          RETURN      
      ELSE IF ( .NOT. EqualRealNos( 0.0_SiKi, AIMAG( TRH_complex(FFT_Data%N/2+1) ) ) )  THEN
          CALL ProgAbort( 'Error in call to FFT. '// &
                          'The imaginary component at the largest positive frequency must be zero.', TrapErrors )
          ErrStat = ErrID_Fatal         ! The code can't get here unless PRESENT(ErrStat)
          RETURN
      END IF


         ! Populate the array for the frequency information.  Only the first half is populated (note that
         ! this algorithm does not make any assumptions about double sided conjugate pairing)
      TRH_complex_return = CMPLX(0.0_SiKi,0.0_SiKi)
      DO I=1,FFT_Data%N/2
         TRH_complex_return(I) = TRH_complex(I)
      ENDDO

        ! FFTPACK 5.1 CFFT1B takes COMPLEX arrays directly
      CALL CFFT1B(FFT_Data%N, 1, TRH_complex_return, SIZE(TRH_complex_return), &
                  FFT_Data%wSave, SIZE(FFT_Data%wSave), wWork, FFT_Data%LenWork, IER)
      IF (IER /= 0) THEN
          CALL ProgAbort( 'Error in complex FFT (CFFT1B).', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

         ! Apply normalization, if any

      IF (FFT_Data%Normalize) THEN
          TRH_complex_return(1:FFT_Data%N) = FFT_Data%InvN * TRH_complex_return(1:FFT_Data%N)
      ENDIF

   END SUBROUTINE ApplyCFFT
  !------------------------------------------------------------------------
   SUBROUTINE ApplyCFFT_f( TRH_complex, FFT_Data, ErrStat )
!FIXME: THIS ROUTINE HAS NOT BEEN TESTED!!!!!
         ! Perform Forward complex FFT: 
         ! give an array TRH, of complex amplitudes in the time domain,
         ! return TRH, a complex array in the frequency domain

      IMPLICIT                         NONE

      COMPLEX(SiKi), INTENT(INOUT)  :: TRH_complex(:)
      TYPE(FFT_DataType), INTENT(IN):: FFT_Data             ! the handle to this instance of the FFT Module
      INTEGER, INTENT(OUT), OPTIONAL:: ErrStat
      
      REAL(SiKi)                    :: wWork(FFT_Data%LenWork)
      INTEGER                       :: IER
      LOGICAL                       :: TrapErrors
      
      
      IF ( PRESENT(ErrStat) ) THEN
         TrapErrors = .TRUE.
         ErrStat = ErrID_None
      ELSE
         TrapErrors = .FALSE.         
      END IF


        ! Make sure the array isn't too small

      IF ( SIZE(TRH_complex) < FFT_Data%N )  THEN
          CALL ProgAbort( 'Error in call to FFT.  Array size is not large enough.', TrapErrors )
          ErrStat = ErrID_Fatal         ! The code can't get here unless PRESENT(ErrStat)
          RETURN
      END IF
      
      IF ( FFT_Data%TransformType /= Fourier_trans ) THEN
          CALL ProgAbort( 'Error in call to FFT. FFT_Data not initialized for Fourier transform.', TrapErrors )
          ErrStat = ErrID_Fatal
          RETURN
      END IF            

        ! FFTPACK 5.1 CFFT1F takes COMPLEX arrays directly

      CALL CFFT1F(FFT_Data%N, 1, TRH_complex, SIZE(TRH_complex), &
                  FFT_Data%wSave, SIZE(FFT_Data%wSave), wWork, FFT_Data%LenWork, IER)
      IF (IER /= 0) THEN
          CALL ProgAbort( 'Error in complex FFT (CFFT1F).', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      IF (FFT_Data%Normalize) THEN
          TRH_complex(1:FFT_Data%N) = FFT_Data%InvN * TRH_complex(1:FFT_Data%N)
      ENDIF

   END SUBROUTINE ApplyCFFT_f
  !------------------------------------------------------------------------
   SUBROUTINE ApplyFFT( TRH, FFT_Data, ErrStat )
         ! Perform Backward FFT: given TRH, a REAL array representing complex numbers,
         ! return an array TRH, of real numbers.
         !     CALL FOURTH ( TRH, NumSteps, 1, WorkT, NumSteps+2 ) ! Sandia

      IMPLICIT                         NONE

      REAL(SiKi), INTENT(INOUT)     :: TRH(:)
      TYPE(FFT_DataType), INTENT(IN):: FFT_Data             ! the handle to this instance of the FFT Module
      INTEGER, INTENT(OUT), OPTIONAL:: ErrStat
      
      REAL(SiKi)                    :: wWork(FFT_Data%LenWork)
      INTEGER                       :: IER
      LOGICAL                       :: TrapErrors
      
         
      IF ( PRESENT(ErrStat) ) THEN
         TrapErrors = .TRUE.
         ErrStat = ErrID_None
      ELSE
         TrapErrors = .FALSE.         
      END IF



        ! Make sure the array isn't too small

      IF ( SIZE(TRH) < FFT_Data%N )  THEN
          CALL ProgAbort( 'Error in call to FFT.  Array size is not large enough.', TrapErrors )
          ErrStat = ErrID_Fatal         ! The code can't get here unless PRESENT(ErrStat)
          RETURN
      END IF
      
      IF ( FFT_Data%TransformType /= Fourier_trans ) THEN
          CALL ProgAbort( 'Error in call to FFT. FFT_Data not initialized for Fourier transform.', TrapErrors )
          ErrStat = ErrID_Fatal
          RETURN
      END IF            

        ! Perform the FFT with a FFTPACK 5.1 routine
        ! FFTPACK 5.1 RFFT1B includes internal normalization (HALF/HALFM) that
        ! 4.1's RFFTB did not. Pre-scale to get 4.1-equivalent un-normalized result.

      TRH(2:FFT_Data%N-1:2) = 2.0_SiKi * TRH(2:FFT_Data%N-1:2)
      TRH(3:FFT_Data%N-1:2) = -2.0_SiKi * TRH(3:FFT_Data%N-1:2)

      CALL RFFT1B(FFT_Data%N, 1, TRH, SIZE(TRH), FFT_Data%wSave, SIZE(FFT_Data%wSave), &
                  wWork, FFT_Data%LenWork, IER)
      IF (IER /= 0) THEN
          CALL ProgAbort( 'Error in FFT (RFFT1B).', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      IF (FFT_Data%Normalize) THEN
          TRH(1:FFT_Data%N) = FFT_Data%InvN * TRH(1:FFT_Data%N)
      ENDIF

   END SUBROUTINE ApplyFFT
  !------------------------------------------------------------------------
   SUBROUTINE ApplyFFT_f( TRH, FFT_Data, ErrStat )
         ! Perform Forward FFT: 
         ! give an array TRH, of real numbers,
         ! return TRH, a REAL array representing complex numbers (magnitude and phase).

      IMPLICIT                         NONE

      REAL(SiKi), INTENT(INOUT)     :: TRH(:)
      TYPE(FFT_DataType), INTENT(IN):: FFT_Data             ! the handle to this instance of the FFT Module
      INTEGER, INTENT(OUT), OPTIONAL:: ErrStat
      
      REAL(SiKi)                    :: wWork(FFT_Data%LenWork)
      INTEGER                       :: IER
      LOGICAL                       :: TrapErrors
      
         
      IF ( PRESENT(ErrStat) ) THEN
         TrapErrors = .TRUE.
         ErrStat = ErrID_None
      ELSE
         TrapErrors = .FALSE.         
      END IF


        ! Make sure the array isn't too small

      IF ( SIZE(TRH) < FFT_Data%N )  THEN
          CALL ProgAbort( 'Error in call to FFT.  Array size is not large enough.', TrapErrors )
          ErrStat = ErrID_Fatal         ! The code can't get here unless PRESENT(ErrStat)
          RETURN
      END IF
      
      IF ( FFT_Data%TransformType /= Fourier_trans ) THEN
          CALL ProgAbort( 'Error in call to FFT. FFT_Data not initialized for Fourier transform.', TrapErrors )
          ErrStat = ErrID_Fatal
          RETURN
      END IF            

        ! Perform the FFT with a FFTPACK 5.1 routine
        ! FFTPACK 5.1 RFFT1F includes internal normalization (SN=1/N, TSN=2/N, TSNM=-2/N)
        ! that 4.1's RFFTF did not. Post-scale to recover 4.1-equivalent un-normalized output.

      CALL RFFT1F(FFT_Data%N, 1, TRH, SIZE(TRH), FFT_Data%wSave, SIZE(FFT_Data%wSave), &
                  wWork, FFT_Data%LenWork, IER)
      IF (IER /= 0) THEN
          CALL ProgAbort( 'Error in FFT (RFFT1F).', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      TRH(1) = REAL(FFT_Data%N, SiKi) * TRH(1)
      TRH(2:FFT_Data%N-1:2) = REAL(FFT_Data%N, SiKi) / 2.0_SiKi * TRH(2:FFT_Data%N-1:2)
      TRH(3:FFT_Data%N-1:2) = -REAL(FFT_Data%N, SiKi) / 2.0_SiKi * TRH(3:FFT_Data%N-1:2)
      IF (MOD(FFT_Data%N, 2) == 0) TRH(FFT_Data%N) = REAL(FFT_Data%N, SiKi) * TRH(FFT_Data%N)

      IF (FFT_Data%Normalize) THEN
          TRH(1:FFT_Data%N) = FFT_Data%InvN * TRH(1:FFT_Data%N)
      ENDIF

   END SUBROUTINE ApplyFFT_f
!------------------------------------------------------------------------
   SUBROUTINE ApplyFFT_cx( TRH, TRH_complex, FFT_Data, ErrStat )
         ! Perform Backward FFT: given TRH, a REAL array representing complex numbers,
         ! return an array TRH, of real numbers.

      IMPLICIT                         NONE

      REAL(SiKi),    INTENT(OUT)    :: TRH(:)
      COMPLEX(SiKi), INTENT(IN)     :: TRH_complex(:)
      TYPE(FFT_DataType), INTENT(IN):: FFT_Data             ! the handle to this instance of the FFT Module
      INTEGER, INTENT(OUT), OPTIONAL:: ErrStat

      REAL(SiKi)                    :: wWork(FFT_Data%LenWork)
      INTEGER                       :: I
      INTEGER                       :: Indx
      INTEGER                       :: IER
      
      LOGICAL                       :: TrapErrors

         
      IF ( PRESENT(ErrStat) ) THEN
         TrapErrors = .TRUE.
         ErrStat = ErrID_None
      ELSE
         TrapErrors = .FALSE.         
      END IF



        ! Make sure the arrays aren't too small

      IF ( ( SIZE(TRH) < FFT_Data%N ) .OR. ( SIZE(TRH_complex) < ( FFT_Data%N/2 + 1 ) ) )  THEN
          CALL ProgAbort( 'Error in call to FFT.  Array size is not large enough.', TrapErrors )
          ErrStat = ErrID_Fatal         ! The code can't get here unless PRESENT(ErrStat)
          RETURN
      END IF
      
      IF ( FFT_Data%TransformType /= Fourier_trans ) THEN
          CALL ProgAbort( 'Error in call to FFT. FFT_Data not initialized for Fourier transform.', TrapErrors )
          ErrStat = ErrID_Fatal
          RETURN
      END IF      

        ! Make sure that the imaginary components at the zeroeth and largest
        ! positive frequency are zero, else abort.

      IF ( .NOT. EqualRealNos( 0.0_SiKi, AIMAG( TRH_complex(1    ) ) ) ) THEN
          CALL ProgAbort( 'Error in call to FFT.  The imaginary component at the zeroeth frequency must be zero.', TrapErrors )
          ErrStat = ErrID_Fatal         ! The code can't get here unless PRESENT(ErrStat)
          RETURN      
      ELSE IF ( .NOT. EqualRealNos( 0.0_SiKi, AIMAG( TRH_complex(FFT_Data%N/2+1) ) ) )  THEN
          CALL ProgAbort( 'Error in call to FFT. '// &
                          'The imaginary component at the largest positive frequency must be zero.', TrapErrors )
          ErrStat = ErrID_Fatal         ! The code can't get here unless PRESENT(ErrStat)
          RETURN
      END IF

        ! Initialize the TRH array with Complex numbers

      TRH(1) = REAL( TRH_complex(1    ) )

      Indx = 1
      DO I=2,FFT_Data%N-2, 2
        Indx     = Indx + 1  ! I/2 + 1

        TRH(I)   =  REAL( TRH_complex(Indx) )
        TRH(I+1) = AIMAG( TRH_complex(Indx) )
      ENDDO

      TRH(FFT_Data%N) = REAL( TRH_complex(FFT_Data%N/2+1) )


        ! Perform the FFT with a FFTPACK 5.1 routine
        ! FFTPACK 5.1 RFFT1B includes internal normalization (HALF/HALFM) that
        ! 4.1's RFFTB did not. Pre-scale to get 4.1-equivalent un-normalized result.

      TRH(2:FFT_Data%N-1:2) = 2.0_SiKi * TRH(2:FFT_Data%N-1:2)
      TRH(3:FFT_Data%N-1:2) = -2.0_SiKi * TRH(3:FFT_Data%N-1:2)

      CALL RFFT1B(FFT_Data%N, 1, TRH, SIZE(TRH), FFT_Data%wSave, SIZE(FFT_Data%wSave), &
                  wWork, FFT_Data%LenWork, IER)
      IF (IER /= 0) THEN
          CALL ProgAbort( 'Error in FFT (RFFT1B).', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      IF (FFT_Data%Normalize) THEN
          TRH(1:FFT_Data%N) = FFT_Data%InvN * TRH(1:FFT_Data%N)
      ENDIF


   END SUBROUTINE ApplyFFT_cx
  !------------------------------------------------------------------------
   SUBROUTINE ApplySINT( TRH, FFT_Data, ErrStat )
         ! Perform sine transform.

      IMPLICIT                         NONE

      REAL(SiKi), INTENT(INOUT)     :: TRH(:)
      TYPE(FFT_DataType), INTENT(IN):: FFT_Data             ! the handle to this instance of the FFT Module
      INTEGER, INTENT(OUT), OPTIONAL:: ErrStat
      
      REAL(SiKi)                    :: wWork(FFT_Data%LenWork)
      INTEGER                       :: IER
      LOGICAL                       :: TrapErrors
      
         
      IF ( PRESENT(ErrStat) ) THEN
         TrapErrors = .TRUE.
         ErrStat = ErrID_None
      ELSE
         TrapErrors = .FALSE.         
      END IF


        ! Make sure the array isn't too small

      IF ( SIZE(TRH) < FFT_Data%N )  THEN
          CALL ProgAbort( 'Error in call to sine transform.  Array size is not large enough.', TrapErrors )
          ErrStat = ErrID_Fatal
          RETURN
      END IF
      
      IF ( FFT_Data%TransformType /= SIN_trans ) THEN
          CALL ProgAbort( 'Error in call to sine transform. FFT_Data not initialized for sine transform.', TrapErrors )
          ErrStat = ErrID_Fatal
          RETURN
      END IF

        ! Make sure that the value at the zeroeth and largest positive
        ! frequency are zero, else abort.

      IF ( TRH(1) /= 0.0 )  THEN
          CALL ProgAbort( 'Error in call to FFT.  The value at the zeroeth frequency must be zero.', TrapErrors )
          ErrStat = ErrID_Fatal
          RETURN
      ELSE IF ( TRH(FFT_Data%N) /= 0.0 ) THEN
          CALL ProgAbort( 'Error in call to FFT.  The value at the largest positive frequency must be zero.', TrapErrors )
          ErrStat = ErrID_Fatal
          RETURN
      END IF

        ! Perform the sine transform with a FFTPACK 5.1 routine.
        ! SINT1B produces half the old SINT result, so post-scale by 2.

      CALL SINT1B(FFT_Data%N-2, 1, TRH(2:), FFT_Data%N-1, FFT_Data%wSave, SIZE(FFT_Data%wSave), &
                  wWork, FFT_Data%LenWork, IER)
      IF (IER /= 0) THEN
          CALL ProgAbort( 'Error in sine transform (SINT1B).', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      TRH(2:FFT_Data%N-1) = 2.0_SiKi * TRH(2:FFT_Data%N-1)

      IF (FFT_Data%Normalize) THEN
          TRH(1:FFT_Data%N) = FFT_Data%InvN * TRH(1:FFT_Data%N)
      ENDIF

   END SUBROUTINE ApplySINT
  !------------------------------------------------------------------------
   SUBROUTINE ExitCFFT(FFT_Data, ErrStat)
   
      TYPE(FFT_DataType), INTENT(INOUT) :: FFT_Data             ! the handle to this instance of the FFT Module
      INTEGER, INTENT(OUT), OPTIONAL    :: ErrStat
      INTEGER                           :: Alloc_Stat 

        ! This subroutine cleans up the backward FFT working space

      FFT_Data%N = -1
      FFT_Data%TransformType = Undef_trans
      
      Alloc_Stat = 0
      IF ( ALLOCATED (FFT_Data%wSave)    ) DEALLOCATE( FFT_Data%wSave, STAT=Alloc_Stat )

      IF ( PRESENT( ErrStat ) ) ErrStat = Alloc_Stat

   END SUBROUTINE ExitCFFT
  !------------------------------------------------------------------------
   SUBROUTINE ExitCOST(FFT_Data, ErrStat)
   
      TYPE(FFT_DataType), INTENT(INOUT) :: FFT_Data             ! the handle to this instance of the FFT Module
      INTEGER, INTENT(OUT), OPTIONAL    :: ErrStat
      INTEGER                           :: Alloc_Stat


        ! This subroutine cleans up the cosine transform working space

      FFT_Data%N = -1
      FFT_Data%TransformType = Undef_trans

      Alloc_Stat = 0
      IF ( ALLOCATED (FFT_Data%wSave)    ) DEALLOCATE( FFT_Data%wSave, STAT=Alloc_Stat )

      IF ( PRESENT( ErrStat ) ) ErrStat = Alloc_Stat

   END SUBROUTINE ExitCOST
  !------------------------------------------------------------------------
   SUBROUTINE ExitFFT(FFT_Data, ErrStat)
   
      TYPE(FFT_DataType), INTENT(INOUT) :: FFT_Data             ! the handle to this instance of the FFT Module
      INTEGER, INTENT(OUT), OPTIONAL    :: ErrStat
      INTEGER                           :: Alloc_Stat 

        ! This subroutine cleans up the backward FFT working space

      FFT_Data%N = -1
      FFT_Data%TransformType = Undef_trans
      
      Alloc_Stat = 0
      IF ( ALLOCATED (FFT_Data%wSave)    ) DEALLOCATE( FFT_Data%wSave, STAT=Alloc_Stat )

      IF ( PRESENT( ErrStat ) ) ErrStat = Alloc_Stat

   END SUBROUTINE ExitFFT
  !------------------------------------------------------------------------
   SUBROUTINE ExitSINT(FFT_Data, ErrStat)
   
      TYPE(FFT_DataType), INTENT(INOUT) :: FFT_Data             ! the handle to this instance of the FFT Module
      INTEGER, INTENT(OUT), OPTIONAL    :: ErrStat
      INTEGER                           :: Alloc_Stat 

        ! This subroutine cleans up the sine transform working space

      FFT_Data%N = -1
      FFT_Data%TransformType = Undef_trans
      
      Alloc_Stat = 0
      IF ( ALLOCATED (FFT_Data%wSave)    ) DEALLOCATE( FFT_Data%wSave, STAT=Alloc_Stat )

      IF ( PRESENT( ErrStat ) ) ErrStat = Alloc_Stat            

   END SUBROUTINE ExitSINT
  !------------------------------------------------------------------------
   SUBROUTINE CheckFFTPACKRealKind( ErrStat )

        ! This subroutine verifies that fftpack5.1.f was compiled with a default REAL
        ! kind of SiKi.  FFTPACK 5.1 declares its arrays as bare REAL/COMPLEX, but this
        ! wrapper hands it explicitly kinded REAL(SiKi)/COMPLEX(SiKi) buffers and passes
        ! their element counts as LENSAV/LENWRK.  If the build promotes the default REAL
        ! to 8 bytes in fftpack5.1.f (DOUBLE_PRECISION does this by default, via
        ! -fdefault-real-8 for GNU or -real-size 64 for Intel), FFTPACK writes twice as
        ! many bytes as those buffers hold and silently corrupts the heap and stack.
        ! See the FFTPACK_SOURCES block in modules/nwtc-library/CMakeLists.txt.

      IMPLICIT                         NONE

      INTEGER, INTENT(OUT),OPTIONAL :: ErrStat        ! returns non-zero if an error occurred

      INTEGER, EXTERNAL             :: FFTPACK_REALKIND  ! from src/NetLib/fftpack/fftpack_kind.f


      IF ( PRESENT(ErrStat) ) ErrStat = ErrID_None

      IF ( FFTPACK_REALKIND() /= SiKi ) THEN
         CALL ProgAbort ( 'FFTPACK 5.1 was compiled with a default REAL kind of '// &
                          TRIM(Num2LStr(FFTPACK_REALKIND()))//', but NWTC_FFTPACK requires '// &
                          TRIM(Num2LStr(SiKi))//'.  The build must suppress default-real promotion '// &
                          'for fftpack5.1.f and fftpack_kind.f.', PRESENT(ErrStat) )
         IF ( PRESENT(ErrStat) ) ErrStat = ErrID_Fatal
      ENDIF


   END SUBROUTINE CheckFFTPACKRealKind
  !------------------------------------------------------------------------

   SUBROUTINE InitCOST( NumSteps, FFT_Data, NormalizeIn, ErrStat )

        ! This subroutine initializes the cosine transform working space

      IMPLICIT                         NONE

      INTEGER, INTENT(IN)           :: NumSteps       ! Number of steps in the array
      INTEGER                       :: Sttus          ! Array allocation status
      INTEGER                       :: IER            ! FFTPACK error return
      INTEGER                       :: LenSav         ! Size of wSave array
      INTEGER                       :: LenWrk         ! Size of wWork array

      TYPE(FFT_DataType),INTENT(OUT):: FFT_Data       ! the handle to this instance of the FFT Module
      LOGICAL, INTENT(IN), OPTIONAL :: NormalizeIn    ! Whether or not to normalize
      INTEGER, INTENT(OUT),OPTIONAL :: ErrStat        ! returns non-zero if an error occurred



      IF ( PRESENT(ErrStat) ) ErrStat = ErrID_None

        ! Verify FFTPACK's default REAL kind matches this wrapper's (SiKi)

      CALL CheckFFTPACKRealKind( ErrStat )
      IF ( PRESENT(ErrStat) ) THEN
         IF ( ErrStat >= AbortErrLev ) RETURN
      ENDIF


        ! Number of timesteps in the time series returned from the cosine transform
        ! N should be odd:

      FFT_Data%N  = NumSteps

      IF ( MOD(FFT_Data%N,2) /= 1 ) THEN
         CALL ProgAbort ( 'The number of steps in the cosine transform must be odd', PRESENT(ErrStat) )
         ErrStat = ErrID_Fatal
         RETURN
      ENDIF

        ! Determine if we should normalize the cosine transform:

      IF ( PRESENT( NormalizeIn ) ) THEN
          FFT_Data%Normalize = NormalizeIn
          FFT_Data%InvN      = 1. / ( FFT_Data%N - 1 )
      ELSE
          FFT_Data%Normalize = .FALSE.
      ENDIF

        ! FFTPACK 5.1 wsave: 2*N + log2(N) + 4;  work: N-1

      LenSav = 2*FFT_Data%N + INT(LOG(REAL(FFT_Data%N, SiKi))/LOG(2.0_SiKi)) + 4
      LenWrk = FFT_Data%N - 1

      ALLOCATE ( FFT_Data%wSave(LenSav) , STAT=Sttus )

      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( 'Error allocating memory for the cosine transform working array.', PRESENT(ErrStat) )
         ErrStat = ErrID_Fatal
         RETURN
      ENDIF

      FFT_Data%LenWork = LenWrk


        ! Initialize the FFTPACK 5.1 working space

      CALL COST1I(FFT_Data%N, FFT_Data%wSave, LenSav, IER)
      IF (IER /= 0) THEN
         CALL ProgAbort ( 'Error initializing cosine transform (COST1I).', PRESENT(ErrStat) )
         IF ( PRESENT(ErrStat) ) ErrStat = ErrID_Fatal
         RETURN
      ENDIF

      FFT_Data%TransformType = COS_trans


   END SUBROUTINE InitCOST
  !------------------------------------------------------------------------
   SUBROUTINE InitCFFT( NumSteps, FFT_Data, NormalizeIn, ErrStat )

        ! This subroutine initializes the backward FFT working space

      IMPLICIT                         NONE

      INTEGER, INTENT(IN)           :: NumSteps       ! Number of steps in the array
      INTEGER                       :: Sttus          ! Array allocation status
      INTEGER                       :: IER            ! FFTPACK error return
      INTEGER                       :: LenSav         ! Size of wSave array
      INTEGER                       :: LenWrk         ! Size of wWork array

      TYPE(FFT_DataType),INTENT(OUT):: FFT_Data       ! the handle to this instance of the FFT Module
      LOGICAL, INTENT(IN), OPTIONAL :: NormalizeIn    ! Whether or not to normalize the FFT
      INTEGER, INTENT(OUT),OPTIONAL :: ErrStat        ! returns non-zero if an error occurred


      IF ( PRESENT(ErrStat) ) ErrStat = ErrID_None

        ! Verify FFTPACK's default REAL kind matches this wrapper's (SiKi)

      CALL CheckFFTPACKRealKind( ErrStat )
      IF ( PRESENT(ErrStat) ) THEN
         IF ( ErrStat >= AbortErrLev ) RETURN
      ENDIF


        ! Number of timesteps in the time series returned from the backward FFT
        ! N should be even:

      FFT_Data%N  = NumSteps

      IF ( MOD(FFT_Data%N,2) /= 0 ) THEN
         CALL ProgAbort ( 'The number of steps in the FFT must be even', PRESENT(ErrStat) ) ! For this Real FFT
         ErrStat = ErrID_Fatal
         RETURN
      ENDIF

        ! Determine if we should normalize the FFT

      IF ( PRESENT( NormalizeIn ) ) THEN
          FFT_Data%Normalize = NormalizeIn
          FFT_Data%InvN      = 1. / FFT_Data%N
      ELSE
          FFT_Data%Normalize = .FALSE.
      ENDIF

        ! FFTPACK 5.1 wsave: 2*N + log2(N) + 4;  work: 2*N

      LenSav = 2*FFT_Data%N + INT(LOG(REAL(FFT_Data%N, SiKi))/LOG(2.0_SiKi)) + 4
      LenWrk = 2*FFT_Data%N

      ALLOCATE ( FFT_Data%wSave(LenSav) , STAT=Sttus )

      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( 'Error allocating memory for the complex FFT working array.', PRESENT(ErrStat) )
         ErrStat = ErrID_Fatal
         RETURN
      ENDIF

      FFT_Data%LenWork = LenWrk


        ! Initialize the FFTPACK 5.1 working space

      CALL CFFT1I(FFT_Data%N, FFT_Data%wSave, LenSav, IER)
      IF (IER /= 0) THEN
         CALL ProgAbort ( 'Error initializing complex FFT (CFFT1I).', PRESENT(ErrStat) )
         IF ( PRESENT(ErrStat) ) ErrStat = ErrID_Fatal
         RETURN
      ENDIF


      FFT_Data%TransformType = Fourier_trans
 
   END SUBROUTINE InitCFFT
  !------------------------------------------------------------------------
   SUBROUTINE InitFFT( NumSteps, FFT_Data, NormalizeIn, ErrStat )

        ! This subroutine initializes the backward FFT working space

      IMPLICIT                         NONE

      INTEGER, INTENT(IN)           :: NumSteps       ! Number of steps in the array
      INTEGER                       :: Sttus          ! Array allocation status
      INTEGER                       :: IER            ! FFTPACK error return
      INTEGER                       :: LenSav         ! Size of wSave array
      INTEGER                       :: LenWrk         ! Size of wWork array

      TYPE(FFT_DataType),INTENT(OUT):: FFT_Data       ! the handle to this instance of the FFT Module
      LOGICAL, INTENT(IN), OPTIONAL :: NormalizeIn    ! Whether or not to normalize the FFT
      INTEGER, INTENT(OUT),OPTIONAL :: ErrStat        ! returns non-zero if an error occurred


      IF ( PRESENT(ErrStat) ) ErrStat = ErrID_None

        ! Verify FFTPACK's default REAL kind matches this wrapper's (SiKi)

      CALL CheckFFTPACKRealKind( ErrStat )
      IF ( PRESENT(ErrStat) ) THEN
         IF ( ErrStat >= AbortErrLev ) RETURN
      ENDIF


        ! Number of timesteps in the time series returned from the backward FFT
        ! N should be even:

      FFT_Data%N  = NumSteps

      IF ( MOD(FFT_Data%N,2) /= 0 ) THEN
         CALL ProgAbort ( 'The number of steps in the FFT must be even', PRESENT(ErrStat) ) ! For this Real FFT
         ErrStat = ErrID_Fatal
         RETURN
      ENDIF

        ! Determine if we should normalize the FFT

      IF ( PRESENT( NormalizeIn ) ) THEN
          FFT_Data%Normalize = NormalizeIn
          FFT_Data%InvN      = 1. / FFT_Data%N
      ELSE
          FFT_Data%Normalize = .FALSE.
          FFT_Data%InvN      = 1.
      ENDIF

        ! FFTPACK 5.1 wsave: N + log2(N) + 4;  work: N

      LenSav = FFT_Data%N + INT(LOG(REAL(FFT_Data%N, SiKi))/LOG(2.0_SiKi)) + 4
      LenWrk = FFT_Data%N

      ALLOCATE ( FFT_Data%wSave(LenSav) , STAT=Sttus )

      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( 'Error allocating memory for the FFT working array.', PRESENT(ErrStat) )
         ErrStat = ErrID_Fatal
         RETURN
      ENDIF

      FFT_Data%LenWork = LenWrk


        ! Initialize the FFTPACK 5.1 working space

      CALL RFFT1I(FFT_Data%N, FFT_Data%wSave, LenSav, IER)
      IF (IER /= 0) THEN
         CALL ProgAbort ( 'Error initializing FFT (RFFT1I).', PRESENT(ErrStat) )
         IF ( PRESENT(ErrStat) ) ErrStat = ErrID_Fatal
         RETURN
      ENDIF

      FFT_Data%TransformType = Fourier_trans
 
   END SUBROUTINE InitFFT
  !------------------------------------------------------------------------
   SUBROUTINE InitSINT( NumSteps, FFT_Data, NormalizeIn, ErrStat )

        ! This subroutine initializes the sine transform working space

      IMPLICIT                         NONE

      INTEGER, INTENT(IN)           :: NumSteps       ! Number of steps in the array
      INTEGER                       :: Sttus          ! Array allocation status
      INTEGER                       :: IER            ! FFTPACK error return
      INTEGER                       :: N_sint         ! Transform length (N-2)
      INTEGER                       :: LenSav         ! Size of wSave array
      INTEGER                       :: LenWrk         ! Size of wWork array

      TYPE(FFT_DataType),INTENT(OUT):: FFT_Data       ! the handle to this instance of the FFT Module
      LOGICAL, INTENT(IN), OPTIONAL :: NormalizeIn    ! Whether or not to normalize
      INTEGER, INTENT(OUT),OPTIONAL :: ErrStat        ! returns non-zero if an error occurred


      IF ( PRESENT(ErrStat) ) ErrStat = ErrID_None

        ! Verify FFTPACK's default REAL kind matches this wrapper's (SiKi)

      CALL CheckFFTPACKRealKind( ErrStat )
      IF ( PRESENT(ErrStat) ) THEN
         IF ( ErrStat >= AbortErrLev ) RETURN
      ENDIF


        ! Number of timesteps in the time series returned from the sine transform
        ! N should be odd:

      FFT_Data%N  = NumSteps

      IF ( MOD(FFT_Data%N,2) /= 1 ) THEN
         CALL ProgAbort ( 'The number of steps in the sine transform must be odd.', PRESENT(ErrStat) )
         ErrStat = ErrID_Fatal
         RETURN
      ENDIF

        ! Determine if we should normalize the sine transform:

      IF ( PRESENT( NormalizeIn ) ) THEN
          FFT_Data%Normalize = NormalizeIn
          FFT_Data%InvN      = 1. / ( FFT_Data%N - 1 )
      ELSE
          FFT_Data%Normalize = .FALSE.
      ENDIF

        ! FFTPACK 5.1: sine transform length is N-2 (interior points only)

      N_sint = FFT_Data%N - 2
      LenSav = N_sint/2 + N_sint + INT(LOG(REAL(N_sint, SiKi))/LOG(2.0_SiKi)) + 4
      LenWrk = 2*N_sint + 2

      ALLOCATE ( FFT_Data%wSave(LenSav) , STAT=Sttus )

      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( 'Error allocating memory for the sine transform working array.', PRESENT(ErrStat) )
         ErrStat = ErrID_Fatal
         RETURN
      ENDIF

      FFT_Data%LenWork = LenWrk


        ! Initialize the FFTPACK 5.1 working space

      CALL SINT1I(N_sint, FFT_Data%wSave, LenSav, IER)
      IF (IER /= 0) THEN
         CALL ProgAbort ( 'Error initializing sine transform (SINT1I).', PRESENT(ErrStat) )
         IF ( PRESENT(ErrStat) ) ErrStat = ErrID_Fatal
         RETURN
      ENDIF


      FFT_Data%TransformType = SIN_trans
      

   END SUBROUTINE InitSINT
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  ! 2D FFT ROUTINES
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  ! DESCRIPTION OF THE 2D REAL FOURIER TRANSFORM:
  !
  ! Given a real L x M array R, RFFT2F computes the normalized 2D DFT.
  ! The output is stored in a packed format in the same L x M array.
  ! RFFT2B computes the inverse, recovering the original signal exactly.
  !
  ! Forward * Backward = Identity (FFTPACK 5.1 normalizes internally)
  ! With NormalizeIn=.TRUE., each call additionally divides by L*M.
  !------------------------------------------------------------------------
  ! DESCRIPTION OF THE 2D COMPLEX FOURIER TRANSFORM:
  !
  ! Given a complex L x M array C, CFFT2F computes the normalized 2D DFT.
  ! CFFT2B computes the inverse, recovering the original signal exactly.
  !
  ! Forward * Backward = Identity (FFTPACK 5.1 normalizes internally)
  ! With NormalizeIn=.TRUE., each call additionally divides by L*M.
  !------------------------------------------------------------------------

   SUBROUTINE InitFFT2D( L, M, FFT_Data, NormalizeIn, ErrStat )

      IMPLICIT                         NONE

      INTEGER, INTENT(IN)           :: L              ! Number of rows
      INTEGER, INTENT(IN)           :: M              ! Number of columns
      INTEGER                       :: Sttus
      INTEGER                       :: IER
      INTEGER                       :: LenSav
      INTEGER                       :: LenWrk

      TYPE(FFT2D_DataType),INTENT(OUT) :: FFT_Data
      LOGICAL, INTENT(IN), OPTIONAL :: NormalizeIn
      INTEGER, INTENT(OUT),OPTIONAL :: ErrStat


      IF ( PRESENT(ErrStat) ) ErrStat = ErrID_None

        ! Verify FFTPACK's default REAL kind matches this wrapper's (SiKi)

      CALL CheckFFTPACKRealKind( ErrStat )
      IF ( PRESENT(ErrStat) ) THEN
         IF ( ErrStat >= AbortErrLev ) RETURN
      ENDIF


      FFT_Data%L = L
      FFT_Data%M = M

      IF ( PRESENT( NormalizeIn ) ) THEN
          FFT_Data%Normalize = NormalizeIn
          FFT_Data%InvN      = 1.0_SiKi / REAL(L*M, SiKi)
      ELSE
          FFT_Data%Normalize = .FALSE.
      ENDIF

      ! RFFT2I requires: L+log2(L)+4 + 2*M+log2(M)+4 + M+log2(M)+4
      LenSav = L + INT(LOG(REAL(L, SiKi))/LOG(2.0_SiKi)) + 4 &
             + 2*M + INT(LOG(REAL(M, SiKi))/LOG(2.0_SiKi)) + 4 &
             + M + INT(LOG(REAL(M, SiKi))/LOG(2.0_SiKi)) + 4
      LenWrk = (L+1)*M

      ALLOCATE ( FFT_Data%wSave(LenSav), STAT=Sttus )

      IF ( Sttus /= 0 ) THEN
         CALL ProgAbort ( 'Error allocating memory for the 2D FFT working array.', PRESENT(ErrStat) )
         IF ( PRESENT(ErrStat) ) ErrStat = ErrID_Fatal
         RETURN
      ENDIF

      FFT_Data%LenWork = LenWrk

      CALL RFFT2I(L, M, FFT_Data%wSave, LenSav, IER)
      IF (IER /= 0) THEN
         CALL ProgAbort ( 'Error initializing 2D FFT (RFFT2I).', PRESENT(ErrStat) )
         IF ( PRESENT(ErrStat) ) ErrStat = ErrID_Fatal
         RETURN
      ENDIF

      FFT_Data%TransformType = Fourier2D_trans

   END SUBROUTINE InitFFT2D
  !------------------------------------------------------------------------
   SUBROUTINE ApplyFFT2D( R, FFT_Data, ErrStat )
         ! Perform backward 2D real FFT: spectral -> spatial.
         ! Recovers original signal exactly (internally normalized).

      IMPLICIT                         NONE

      REAL(SiKi), INTENT(INOUT)     :: R(:,:)
      TYPE(FFT2D_DataType), INTENT(IN) :: FFT_Data
      INTEGER, INTENT(OUT), OPTIONAL:: ErrStat

      REAL(SiKi)                    :: wWork(FFT_Data%LenWork)
      INTEGER                       :: IER
      LOGICAL                       :: TrapErrors


      IF ( PRESENT(ErrStat) ) THEN
         TrapErrors = .TRUE.
         ErrStat = ErrID_None
      ELSE
         TrapErrors = .FALSE.
      END IF

      IF ( SIZE(R,1) < FFT_Data%L .OR. SIZE(R,2) < FFT_Data%M ) THEN
          CALL ProgAbort( 'Error in call to 2D FFT. Array size is not large enough.', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      IF ( FFT_Data%TransformType /= Fourier2D_trans ) THEN
          CALL ProgAbort( 'Error in call to 2D FFT. FFT_Data not initialized for 2D Fourier transform.', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      CALL RFFT2B(SIZE(R,1), FFT_Data%L, FFT_Data%M, R, FFT_Data%wSave, SIZE(FFT_Data%wSave), &
                  wWork, FFT_Data%LenWork, IER)
      IF (IER /= 0) THEN
          CALL ProgAbort( 'Error in 2D FFT (RFFT2B).', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      IF (FFT_Data%Normalize) THEN
          R(1:FFT_Data%L, 1:FFT_Data%M) = FFT_Data%InvN * R(1:FFT_Data%L, 1:FFT_Data%M)
      ENDIF

   END SUBROUTINE ApplyFFT2D
  !------------------------------------------------------------------------
   SUBROUTINE ApplyFFT2D_f( R, FFT_Data, ErrStat )
         ! Perform forward 2D real FFT: spatial -> spectral.

      IMPLICIT                         NONE

      REAL(SiKi), INTENT(INOUT)     :: R(:,:)
      TYPE(FFT2D_DataType), INTENT(IN) :: FFT_Data
      INTEGER, INTENT(OUT), OPTIONAL:: ErrStat

      REAL(SiKi)                    :: wWork(FFT_Data%LenWork)
      INTEGER                       :: IER
      LOGICAL                       :: TrapErrors


      IF ( PRESENT(ErrStat) ) THEN
         TrapErrors = .TRUE.
         ErrStat = ErrID_None
      ELSE
         TrapErrors = .FALSE.
      END IF

      IF ( SIZE(R,1) < FFT_Data%L .OR. SIZE(R,2) < FFT_Data%M ) THEN
          CALL ProgAbort( 'Error in call to 2D FFT. Array size is not large enough.', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      IF ( FFT_Data%TransformType /= Fourier2D_trans ) THEN
          CALL ProgAbort( 'Error in call to 2D FFT. FFT_Data not initialized for 2D Fourier transform.', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      CALL RFFT2F(SIZE(R,1), FFT_Data%L, FFT_Data%M, R, FFT_Data%wSave, SIZE(FFT_Data%wSave), &
                  wWork, FFT_Data%LenWork, IER)
      IF (IER /= 0) THEN
          CALL ProgAbort( 'Error in 2D FFT (RFFT2F).', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      IF (FFT_Data%Normalize) THEN
          R(1:FFT_Data%L, 1:FFT_Data%M) = FFT_Data%InvN * R(1:FFT_Data%L, 1:FFT_Data%M)
      ENDIF

   END SUBROUTINE ApplyFFT2D_f
  !------------------------------------------------------------------------
   SUBROUTINE ExitFFT2D(FFT_Data, ErrStat)

      TYPE(FFT2D_DataType), INTENT(INOUT) :: FFT_Data
      INTEGER, INTENT(OUT), OPTIONAL :: ErrStat

      IF ( PRESENT(ErrStat) ) ErrStat = ErrID_None

      IF ( ALLOCATED(FFT_Data%wSave) ) DEALLOCATE( FFT_Data%wSave )

   END SUBROUTINE ExitFFT2D
  !------------------------------------------------------------------------
   SUBROUTINE InitCFFT2D( L, M, FFT_Data, NormalizeIn, ErrStat )

      IMPLICIT                         NONE

      INTEGER, INTENT(IN)           :: L              ! Number of rows
      INTEGER, INTENT(IN)           :: M              ! Number of columns
      INTEGER                       :: Sttus
      INTEGER                       :: IER
      INTEGER                       :: LenSav
      INTEGER                       :: LenWrk

      TYPE(FFT2D_DataType),INTENT(OUT) :: FFT_Data
      LOGICAL, INTENT(IN), OPTIONAL :: NormalizeIn
      INTEGER, INTENT(OUT),OPTIONAL :: ErrStat


      IF ( PRESENT(ErrStat) ) ErrStat = ErrID_None

        ! Verify FFTPACK's default REAL kind matches this wrapper's (SiKi)

      CALL CheckFFTPACKRealKind( ErrStat )
      IF ( PRESENT(ErrStat) ) THEN
         IF ( ErrStat >= AbortErrLev ) RETURN
      ENDIF


      FFT_Data%L = L
      FFT_Data%M = M

      IF ( PRESENT( NormalizeIn ) ) THEN
          FFT_Data%Normalize = NormalizeIn
          FFT_Data%InvN      = 1.0_SiKi / REAL(L*M, SiKi)
      ELSE
          FFT_Data%Normalize = .FALSE.
      ENDIF

      ! CFFT2I requires: 2*L+log2(L)+2*M+log2(M)+8
      LenSav = 2*L + INT(LOG(REAL(L, SiKi))/LOG(2.0_SiKi)) &
             + 2*M + INT(LOG(REAL(M, SiKi))/LOG(2.0_SiKi)) + 8
      LenWrk = 2*L*M

      ALLOCATE ( FFT_Data%wSave(LenSav), STAT=Sttus )

      IF ( Sttus /= 0 ) THEN
         CALL ProgAbort ( 'Error allocating memory for the 2D complex FFT working array.', PRESENT(ErrStat) )
         IF ( PRESENT(ErrStat) ) ErrStat = ErrID_Fatal
         RETURN
      ENDIF

      FFT_Data%LenWork = LenWrk

      CALL CFFT2I(L, M, FFT_Data%wSave, LenSav, IER)
      IF (IER /= 0) THEN
         CALL ProgAbort ( 'Error initializing 2D complex FFT (CFFT2I).', PRESENT(ErrStat) )
         IF ( PRESENT(ErrStat) ) ErrStat = ErrID_Fatal
         RETURN
      ENDIF

      FFT_Data%TransformType = CFourier2D_trans

   END SUBROUTINE InitCFFT2D
  !------------------------------------------------------------------------
   SUBROUTINE ApplyCFFT2D( C, FFT_Data, ErrStat )
         ! Perform backward 2D complex FFT: spectral -> spatial.
         ! Recovers original signal exactly (internally normalized).

      IMPLICIT                         NONE

      COMPLEX(SiKi), INTENT(INOUT)  :: C(:,:)
      TYPE(FFT2D_DataType), INTENT(IN) :: FFT_Data
      INTEGER, INTENT(OUT), OPTIONAL:: ErrStat

      REAL(SiKi)                    :: wWork(FFT_Data%LenWork)
      INTEGER                       :: IER
      LOGICAL                       :: TrapErrors


      IF ( PRESENT(ErrStat) ) THEN
         TrapErrors = .TRUE.
         ErrStat = ErrID_None
      ELSE
         TrapErrors = .FALSE.
      END IF

      IF ( SIZE(C,1) < FFT_Data%L .OR. SIZE(C,2) < FFT_Data%M ) THEN
          CALL ProgAbort( 'Error in call to 2D complex FFT. Array size is not large enough.', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      IF ( FFT_Data%TransformType /= CFourier2D_trans ) THEN
          CALL ProgAbort( 'Error in call to 2D complex FFT. FFT_Data not initialized for 2D complex transform.', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      CALL CFFT2B(SIZE(C,1), FFT_Data%L, FFT_Data%M, C, FFT_Data%wSave, SIZE(FFT_Data%wSave), &
                  wWork, FFT_Data%LenWork, IER)
      IF (IER /= 0) THEN
          CALL ProgAbort( 'Error in 2D complex FFT (CFFT2B).', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      IF (FFT_Data%Normalize) THEN
          C(1:FFT_Data%L, 1:FFT_Data%M) = FFT_Data%InvN * C(1:FFT_Data%L, 1:FFT_Data%M)
      ENDIF

   END SUBROUTINE ApplyCFFT2D
  !------------------------------------------------------------------------
   SUBROUTINE ApplyCFFT2D_f( C, FFT_Data, ErrStat )
         ! Perform forward 2D complex FFT: spatial -> spectral.

      IMPLICIT                         NONE

      COMPLEX(SiKi), INTENT(INOUT)  :: C(:,:)
      TYPE(FFT2D_DataType), INTENT(IN) :: FFT_Data
      INTEGER, INTENT(OUT), OPTIONAL:: ErrStat

      REAL(SiKi)                    :: wWork(FFT_Data%LenWork)
      INTEGER                       :: IER
      LOGICAL                       :: TrapErrors


      IF ( PRESENT(ErrStat) ) THEN
         TrapErrors = .TRUE.
         ErrStat = ErrID_None
      ELSE
         TrapErrors = .FALSE.
      END IF

      IF ( SIZE(C,1) < FFT_Data%L .OR. SIZE(C,2) < FFT_Data%M ) THEN
          CALL ProgAbort( 'Error in call to 2D complex FFT. Array size is not large enough.', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      IF ( FFT_Data%TransformType /= CFourier2D_trans ) THEN
          CALL ProgAbort( 'Error in call to 2D complex FFT. FFT_Data not initialized for 2D complex transform.', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      CALL CFFT2F(SIZE(C,1), FFT_Data%L, FFT_Data%M, C, FFT_Data%wSave, SIZE(FFT_Data%wSave), &
                  wWork, FFT_Data%LenWork, IER)
      IF (IER /= 0) THEN
          CALL ProgAbort( 'Error in 2D complex FFT (CFFT2F).', TrapErrors )
          IF (PRESENT(ErrStat)) ErrStat = ErrID_Fatal
          RETURN
      END IF

      IF (FFT_Data%Normalize) THEN
          C(1:FFT_Data%L, 1:FFT_Data%M) = FFT_Data%InvN * C(1:FFT_Data%L, 1:FFT_Data%M)
      ENDIF

   END SUBROUTINE ApplyCFFT2D_f
  !------------------------------------------------------------------------
   SUBROUTINE ExitCFFT2D(FFT_Data, ErrStat)

      TYPE(FFT2D_DataType), INTENT(INOUT) :: FFT_Data
      INTEGER, INTENT(OUT), OPTIONAL :: ErrStat

      IF ( PRESENT(ErrStat) ) ErrStat = ErrID_None

      IF ( ALLOCATED(FFT_Data%wSave) ) DEALLOCATE( FFT_Data%wSave )

   END SUBROUTINE ExitCFFT2D
  !------------------------------------------------------------------------


END MODULE NWTC_FFTPACK
!=======================================================================
