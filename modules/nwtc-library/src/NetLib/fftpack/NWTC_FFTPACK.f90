   ! NOTE: This MODULE isused in HydroDyn and TurbSim.
   ! BJJ: 02/22/2008: Updated to work with NWTC_Library v1.01.09
   !      all Abort() functions changed to ProgAbort()
   ! BJJ: 12/03/2010: Updated to add optional ErrStat return values
   !                  instead of aborting on errors (note not all changes documented)
   ! BJJ: 12/20/2010: Updated to add defined type and remove global data variables
   !                  Also updated to check that transform has been initialized for the
   !                    correct type (to avoid having wSave too small)
   ! ADP: 07/28/2014: Added in the complex FFT routines from fftpack v. 4.1
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

   TYPE, PUBLIC :: FFT_DataType
      PRIVATE
      REAL(SiKi)                       :: InvN          = 0.0_SiKi      ! Normalization constant
      REAL(SiKi), ALLOCATABLE          :: wSave(:)                      ! Working array for performing transforms
      INTEGER                          :: N             = -1            ! Number of steps
      LOGICAL                          :: Normalize     = .FALSE.       ! Whether or not to normalize
      INTEGER                          :: TransformType = Undef_trans   ! the type of transfer function this is for
   END TYPE FFT_DataType      


!------------------------------------------------------------------------
CONTAINS

   SUBROUTINE ApplyCOST( TRH, FFT_Data, ErrStat )
   
         ! Perform cosine transform.

      IMPLICIT                         NONE

      REAL(SiKi), INTENT(INOUT)     :: TRH(:)
      TYPE(FFT_DataType), INTENT(IN):: FFT_Data             ! the handle to this instance of the FFT Module
      
      INTEGER, INTENT(OUT), OPTIONAL:: ErrStat
      
      LOGICAL                       :: TrapErrors
      

         
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
      

        ! Perform the cosine transform with a FFTpack routine

      CALL COST(FFT_Data%N, TRH, FFT_Data%wSave) ! FFTpack routine

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

      REAL(SiKi), ALLOCATABLE       :: TRH(:)   ! real array to help process the complex-array that fftpack defines as IMPLICIT (real)

      INTEGER                       :: I
     
      INTEGER(IntKi)                :: ErrStatTmp 
      LOGICAL                       :: TrapErrors
      character(ErrMsgLen)          :: ErrMsg

      ErrStatTmp  = ErrID_None
         
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


      CALL AllocAry( TRH, 2*size(TRH_complex_return,1), 'ApplyCFFT:TRH', ErrStat, ErrMsg  ) !allocate two real for each complex variable
         IF (ErrStat >= AbortErrLev) THEN
            CALL WrScr( TRIM(ErrMsg) )
            RETURN
         END IF

         !TRH = TRANSFER(TRH_complex_return, TRH)  ! this function apparently uses stack space and is causing stack overflow on large models
      do i=1,size(TRH_complex_return,1)
         TRH(2*i-1)   = REAL(TRH_complex_return(i))
         TRH(2*i  ) = AIMAG(TRH_complex_return(i))
      end do
      

      CALL CFFTB(FFT_Data%N, TRH, FFT_Data%wSave)

      ! put real values back into complex array
         !TRH = TRH_complex_return = TRANSFER(TRH, TRH_complex)  ! this function apparently uses stack space and is causing stack overflow on large models
      do i=1,size(TRH_complex_return,1)
         TRH_complex_return(i) = CMPLX(TRH(2*i-1),TRH(2*i))
      end do
      
      DEALLOCATE(TRH)


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
      
      
      REAL(SiKi), ALLOCATABLE       :: TRH(:)
      
      LOGICAL                       :: TrapErrors
      character(ErrMsgLen)          :: ErrMsg   
      
      
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

        ! Perform the FFT with a FFTpack routine

      CALL AllocAry( TRH, 2*size(TRH,1), 'ApplyCFFT_f:TRH', ErrStat, ErrMsg )
         IF (ErrStat >= AbortErrLev) THEN
            CALL WrScr( TRIM(ErrMsg) )
            RETURN
         END IF
         
      TRH = TRANSFER(TRH_complex, TRH)
        
      CALL CFFTF(FFT_Data%N, TRH, FFT_Data%wSave) ! FFTpack routine
   
      ! put real values back into complex array
      TRH_complex = TRANSFER(TRH, TRH_complex)
      DEALLOCATE(TRH)
      
      
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

        ! Perform the FFT with a FFTpack routine

      CALL RFFTB(FFT_Data%N, TRH, FFT_Data%wSave) ! FFTpack routine

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

        ! Perform the FFT with a FFTpack routine

      CALL RFFTF(FFT_Data%N, TRH, FFT_Data%wSave) ! FFTpack routine
   
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

      INTEGER                       :: I
      INTEGER                       :: Indx
      
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


        ! Perform the FFT with a FFTpack routine

      CALL RFFTB(FFT_Data%N, TRH, FFT_Data%wSave)

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

        ! Perform the sine transform with a FFTpack routine

      CALL SINT(FFT_Data%N-2, TRH(2:FFT_Data%N-1), FFT_Data%wSave) ! FFTpack routine

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
   SUBROUTINE InitCOST( NumSteps, FFT_Data, NormalizeIn, ErrStat )

        ! This subroutine initializes the cosine transform working space

      IMPLICIT                         NONE

      INTEGER, INTENT(IN)           :: NumSteps       ! Number of steps in the array
      INTEGER                       :: Sttus          ! Array allocation status

      TYPE(FFT_DataType),INTENT(OUT):: FFT_Data       ! the handle to this instance of the FFT Module
      LOGICAL, INTENT(IN), OPTIONAL :: NormalizeIn    ! Whether or not to normalize
      INTEGER, INTENT(OUT),OPTIONAL :: ErrStat        ! returns non-zero if an error occurred



      IF ( PRESENT(ErrStat) ) ErrStat = ErrID_None

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

        ! According to FFTPACK documentation, the working array must be at
        ! least size 3N+15

      ALLOCATE ( FFT_Data%wSave(3*FFT_Data%N + 15) , STAT=Sttus )

      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( 'Error allocating memory for the cosine transform working array.', PRESENT(ErrStat) )
         ErrStat = Sttus
         RETURN
      ENDIF


        ! Initialize the FFTPACK working space

      CALL COSTI(FFT_Data%N, FFT_Data%wSave)

      FFT_Data%TransformType = COS_trans


   END SUBROUTINE InitCOST
  !------------------------------------------------------------------------
   SUBROUTINE InitCFFT( NumSteps, FFT_Data, NormalizeIn, ErrStat )

        ! This subroutine initializes the backward FFT working space

      IMPLICIT                         NONE

      INTEGER, INTENT(IN)           :: NumSteps       ! Number of steps in the array
      INTEGER                       :: Sttus          ! Array allocation status

      TYPE(FFT_DataType),INTENT(OUT):: FFT_Data       ! the handle to this instance of the FFT Module
      LOGICAL, INTENT(IN), OPTIONAL :: NormalizeIn    ! Whether or not to normalize the FFT
      INTEGER, INTENT(OUT),OPTIONAL :: ErrStat        ! returns non-zero if an error occurred


      IF ( PRESENT(ErrStat) ) ErrStat = ErrID_None

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

        ! According to FFTPACK documentation, the working array must be at
        ! least size 4N+15

      ALLOCATE ( FFT_Data%wSave(4*FFT_Data%N + 15) , STAT=Sttus )

      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( 'Error allocating memory for the complex FFT working array.', PRESENT(ErrStat) )
         ErrStat = ErrID_Fatal
         RETURN
      ENDIF


        ! Initialize the FFTPACK working space

      CALL CFFTI(FFT_Data%N, FFT_Data%wSave)


      FFT_Data%TransformType = Fourier_trans
 
   END SUBROUTINE InitCFFT
  !------------------------------------------------------------------------
   SUBROUTINE InitFFT( NumSteps, FFT_Data, NormalizeIn, ErrStat )

        ! This subroutine initializes the backward FFT working space

      IMPLICIT                         NONE

      INTEGER, INTENT(IN)           :: NumSteps       ! Number of steps in the array
      INTEGER                       :: Sttus          ! Array allocation status

      TYPE(FFT_DataType),INTENT(OUT):: FFT_Data       ! the handle to this instance of the FFT Module
      LOGICAL, INTENT(IN), OPTIONAL :: NormalizeIn    ! Whether or not to normalize the FFT
      INTEGER, INTENT(OUT),OPTIONAL :: ErrStat        ! returns non-zero if an error occurred


      IF ( PRESENT(ErrStat) ) ErrStat = ErrID_None

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

        ! According to FFTPACK documentation, the working array must be at
        ! least size 2N+15

      ALLOCATE ( FFT_Data%wSave(2*FFT_Data%N + 15) , STAT=Sttus )

      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( 'Error allocating memory for the FFT working array.', PRESENT(ErrStat) )
         ErrStat = ErrID_Fatal
         RETURN
      ENDIF


        ! Initialize the FFTPACK working space

      CALL RFFTI(FFT_Data%N, FFT_Data%wSave)

      FFT_Data%TransformType = Fourier_trans
 
   END SUBROUTINE InitFFT
  !------------------------------------------------------------------------
   SUBROUTINE InitSINT( NumSteps, FFT_Data, NormalizeIn, ErrStat )

        ! This subroutine initializes the sine transform working space

      IMPLICIT                         NONE

      INTEGER, INTENT(IN)           :: NumSteps       ! Number of steps in the array
      INTEGER                       :: Sttus          ! Array allocation status

      TYPE(FFT_DataType),INTENT(OUT):: FFT_Data       ! the handle to this instance of the FFT Module
      LOGICAL, INTENT(IN), OPTIONAL :: NormalizeIn    ! Whether or not to normalize
      INTEGER, INTENT(OUT),OPTIONAL :: ErrStat        ! returns non-zero if an error occurred


      IF ( PRESENT(ErrStat) ) ErrStat = ErrID_None

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

        ! According to FFTPACK documentation, the working array must be at
        ! least size 2.5N+15; however, our N is +2 greater than their N

      ALLOCATE ( FFT_Data%wSave( CEILING( 2.5*(FFT_Data%N-2) ) + 15 ) , STAT=Sttus )

      IF ( Sttus /= 0 )  THEN
         CALL ProgAbort ( 'Error allocating memory for the sine transform working array.', PRESENT(ErrStat) )
         ErrStat = ErrID_Fatal
         RETURN
      ENDIF


        ! Initialize the FFTPACK working space

      CALL SINTI(FFT_Data%N-2, FFT_Data%wSave)


      FFT_Data%TransformType = SIN_trans
      

   END SUBROUTINE InitSINT
  !------------------------------------------------------------------------


END MODULE NWTC_FFTPACK
!=======================================================================
