   
   
! Module TSSubs is part of the TurbSim program.
!
!**************************************************************************************************************
MODULE TSSubs

   USE                     NWTC_Library   
   USE                     ModifiedvKrm_mod

   
use TS_Profiles  
use TS_RandNum
use ts_errors
use TS_VelocitySpectra

   IMPLICIT                NONE



CONTAINS


!=======================================================================
SUBROUTINE CohSpecVMat( LC, NSize )

   ! This subroutine computes the coherence between two points on the grid.
   ! It stores the symmetric coherence matrix, packed into variable "Matrix"
   ! This replaces what formerly was the "ExCoDW" matrix.

USE                           TSMods
USE NWTC_LAPACK

IMPLICIT                      NONE

   ! Passed variables

REAL(ReKi),   INTENT(IN)   :: LC             ! IEC coherency scale parameter
INTEGER,      INTENT(IN)   :: NSize          ! Size of dimension 2 of Matrix

   ! Internal variables

REAL(ReKi), ALLOCATABLE    :: Dist(:)        ! The distance between points
REAL(ReKi), ALLOCATABLE    :: DistU(:)
REAL(ReKi), ALLOCATABLE    :: DistZMExp(:)
REAL(ReKi)                 :: dY             ! the lateral distance between two points
REAL(ReKi)                 :: UM             ! The mean wind speed of the two points
REAL(ReKi)                 :: ZM             ! The mean height of the two points

INTEGER                    :: J
INTEGER                    :: JJ             ! Index of point J
INTEGER                    :: JJ1
INTEGER                    :: JY             ! Index of y-value of point J
INTEGER                    :: JZ             ! Index of z-value of point J
INTEGER                    :: I
INTEGER                    :: IFreq
INTEGER                    :: II             ! The index of point I
INTEGER                    :: Indx
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

V(:,:,:) = 0.0_ReKi    ! initialize the velocity matrix

!--------------------------------------------------------------------------------
! Calculate the distances and other parameters that don't change with frequency
!---------------------------------------------------------------------------------

IF ( SpecModel == SpecModel_IECKAI .or. SpecModel == SpecModel_IECVKM .OR. SpecModel == SpecModel_MODVKM .OR. SpecModel == SpecModel_API ) THEN
   InCohB(:)    = 0.12/LC
   DistZMExp(:) = 1.0
ENDIF

         II = 0
POINT_I: DO IZ=1,p_grid%ZLim   !NumGrid_Z
            DO IY=1,p_grid%IYmax(IZ) !NumGrid_Y

               II = II + 1                            ! Index of point I: S(I)  !equivalent to II = ( IZ - 1 )*NumGrid_Y + IY
               IF (II > p_grid%NPoints) EXIT POINT_I            ! Don't go past the end of the array; this exits the IZ loop

               JJ = 0
POINT_J:       DO JZ=1,IZ
                  DO JY=1,p_grid%IYmax(JZ) !NumGrid_Y

                     JJ = JJ + 1                      ! Index of point J: S(J)  !equivalent to JJ = ( JZ - 1 )*NumGrid_Y + JY

                     IF ( JJ > II )  EXIT POINT_J     ! The coherence matrix is symmetric

                     IF ( IZ > p_grid%NumGrid_Z ) THEN       ! Get the correct location if we're using an extra point for the hub
                        I = p_grid%YLim
                        IF ( JZ > p_grid%NumGrid_Z ) THEN
                           J = p_grid%YLim
                        ELSE
                           J = JY
                        ENDIF
                     ELSE
                        I = IY
                        J = JY
                     ENDIF

                     JJ1       = JJ - 1
                     Indx      = p_grid%NPoints*JJ1 - JJ*JJ1/2 + II   !Index of matrix ExCoDW (now Matrix), coherence between points I & J

                     IF ( .NOT. PeriodicY ) THEN
                        Dist(Indx)= SQRT( ( p_grid%Y(I) - p_grid%Y(J) )**2  + ( p_grid%Z(IZ) - p_grid%Z(JZ) )**2 )
                     ELSE
                        dY = p_grid%Y(I) - p_grid%Y(J)
                        IF (dY > 0.5*p_grid%GridWidth ) THEN
                           dY = dY - p_grid%GridWidth - p_grid%GridRes_Y
                        ELSE IF (dY < -0.5*p_grid%GridWidth ) THEN
                           dY = dY + p_grid%GridWidth + p_grid%GridRes_Y
                        END IF

                        Dist(Indx)= SQRT( ( dY )**2  + ( p_grid%Z(IZ) - p_grid%Z(JZ) )**2 )
                     END IF

                     IF ( SpecModel == SpecModel_IECKAI .or. SpecModel == SpecModel_IECVKM  .OR. SpecModel == SpecModel_MODVKM .OR. SpecModel == SpecModel_API ) THEN
                        DistU(Indx) = Dist(Indx)/UHub
!                           TRH(Indx) = EXP( -InCDec(IVec)*SQRT( ( p_grid%Freq(IFreq) * Dist / UHub )**2 + (0.12*Dist/LC)**2 ) )
                     ELSE
                        UM       = 0.5*( U(IZ) + U(JZ) )
                        ZM       = 0.5*( p_grid%Z(IZ) + p_grid%Z(JZ) )

                        DistU(Indx)     = Dist(Indx)/UM
                        DistZMExp(Indx) = ( Dist(Indx)/ZM )**COHEXP     ! Note: 0**0 = 1

!                       TRH(Indx) = EXP( -InCDec(IVec) * DistZMExp*SQRT( ( p_grid%Freq(IFreq)* DistU )**2 + (InCohB(IVec)*Dist)**2 ) )
                     ENDIF !SpecModel

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

   IF ( ( SpecModel == SpecModel_IECKAI .or. SpecModel == SpecModel_IECVKM .OR. SpecModel == SpecModel_MODVKM .OR. SpecModel == SpecModel_API ) .AND. ( IVec /= 1 ) ) THEN
         ! There is no coherence defined for the v or w component of the IEC spectral models
      IdentityCoh = .TRUE.
   ELSE
      IdentityCoh = .FALSE.
   ENDIF

   CALL WrScr ( '    '//Comp(IVec)//'-component matrices' )

   !--------------------------------------------------------------------------------
   ! Calculate the coherence, Veers' H matrix (CSDs), and the fourier coefficients
   !---------------------------------------------------------------------------------

   DO IFREQ = 1,p_grid%NumFreq

      !---------------------------------------------------
      ! Calculate the coherence and Veers' H matrix (CSDs)
      !---------------------------------------------------
      IF (.NOT. IdentityCoh) THEN

            ! -----------------------------------------------
            ! Create the coherence matrix for this frequency
            ! -----------------------------------------------

         DO II=1,p_grid%NPoints
            DO JJ=1,II

                  JJ1       = JJ - 1
                  Indx      = p_grid%NPoints*JJ1 - JJ*JJ1/2 + II   !Index of matrix ExCoDW (now Matrix), coherence between points I & J

                  TRH(Indx) = EXP( -1.0 * InCDec(IVec) * DistZMExp(Indx)* &
                              SQRT( (p_grid%Freq(IFreq)*DistU(Indx) )**2 + (InCohB(IVec)*Dist(Indx))**2 ) )

            ENDDO ! JJ
         ENDDO ! II

      END IF      
      
      CALL Coh2Coeffs( IdentityCoh, IVec, IFreq, TRH, S, PhaseAngles, V, NSize )

   ENDDO !IFreq

ENDDO !IVec

IF (COH_OUT)  CLOSE( UC )

IF ( ALLOCATED( Dist      ) ) DEALLOCATE( Dist      )
IF ( ALLOCATED( DistU     ) ) DEALLOCATE( DistU     )
IF ( ALLOCATED( DistZMExp ) ) DEALLOCATE( DistZMExp )


RETURN

END SUBROUTINE CohSpecVMat
!=======================================================================
SUBROUTINE CohSpecVMat_API( LC, NSize)

   ! This subroutine computes the coherence between two points on the grid.
   ! It stores the symmetric coherence matrix, packed into variable "Matrix"
   ! This replaces what formerly was the "ExCoDW" matrix.

USE                           TSMods
USE NWTC_LAPACK

IMPLICIT                      NONE

   ! Passed variables

REAL(ReKi),   INTENT(IN)   :: LC             ! IEC coherency scale parameter
INTEGER,      INTENT(IN)   :: NSize          ! Size of dimension 2 of Matrix

   ! Internal variables

!REAL(ReKi), ALLOCATABLE    :: Dist(:)        ! The distance between points
REAL(ReKi), ALLOCATABLE    :: Dist_Y(:)        ! The Y distance between points
REAL(ReKi), ALLOCATABLE    :: Dist_Z(:)        ! The Z distance between points
!mlb REAL(ReKi), ALLOCATABLE    :: Dist_Z12(:)        ! The distance between points (not really a distance!)
REAL(ReKi), ALLOCATABLE    :: Z1Z2(:)        ! Z(IZ)*Z(JZ)

INTEGER                    :: J
INTEGER                    :: JJ             ! Index of point J
INTEGER                    :: JJ1
INTEGER                    :: JY             ! Index of y-value of point J
INTEGER                    :: JZ             ! Index of z-value of point J
INTEGER                    :: I
INTEGER                    :: IFreq
INTEGER                    :: II             ! The index of point I
INTEGER                    :: Indx
INTEGER                    :: IVec           ! wind component, 1=u, 2=v, 3=w
INTEGER                    :: IY             ! The index of the y-value of point I
INTEGER                    :: IZ             ! Index of the z-value of point I
INTEGER                    :: Stat

LOGICAL                    :: IdentityCoh


!REAL :: Coef_QX=1.00
REAL :: Coef_QY=1.00
REAL :: Coef_QZ=1.25
!REAL :: Coef_PX=0.40
REAL :: Coef_PY=0.40
REAL :: Coef_PZ=0.50
!REAL :: Coef_RX=0.92
REAL :: Coef_RY=0.92
REAL :: Coef_RZ=0.85
!REAL :: Coef_AlphaX=2.90
REAL :: Coef_AlphaY=45.0
REAL :: Coef_AlphaZ=13.0
REAL :: Coef_1=1.0!3.28
REAL :: Coef_2=100.0!32.8
!REAL :: Coef_2=10.0!32.8

REAL :: TEMP_Y, TEMP_Z



V(:,:,:) = 0.0    ! initialize the velocity matrix

! ------------ arrays allocated -------------------------
Stat = 0.

!mlb IF ( .NOT. ALLOCATED( Dist ) ) ALLOCATE( Dist(NSize),      STAT=Stat )
!mlb IF ( Stat /= 0 ) THEN
!mlb    CALL TS_Abort('Error allocating '//TRIM( Int2LStr( ReKi*NSize/1024**2 ) )//' MB for the Dist coherence array.')
!mlb ENDIF
IF ( .NOT. ALLOCATED( Dist_Y ) ) ALLOCATE(Dist_Y(NSize),      STAT=Stat )
IF ( Stat /= 0 ) THEN
   CALL TS_Abort('Error allocating '//TRIM( Int2LStr( ReKi*NSize/1024**2 ) )//' MB for the Dist_Y coherence array.')
ENDIF
IF ( .NOT. ALLOCATED( Dist_Z ) ) ALLOCATE( Dist_Z(NSize),      STAT=Stat )
IF ( Stat /= 0 ) THEN
   CALL TS_Abort('Error allocating '//TRIM( Int2LStr( ReKi*NSize/1024**2 ) )//' MB for the Dist_Z coherence array.')
ENDIF
!mlb: IF ( .NOT. ALLOCATED( Dist_Z12 ) ) ALLOCATE( Dist_Z12(NSize),      STAT=Stat )
!mlb:
!mlb: IF ( Stat /= 0 ) THEN
!mlb:    CALL TS_Abort('Error allocating '//TRIM( Int2LStr( ReKi*NSize/1024**2 ) )//' MB for the Dist_Z12 coherence array.')
IF ( .NOT. ALLOCATED( Z1Z2 ) ) ALLOCATE( Z1Z2(NSize),      STAT=Stat )

IF ( Stat /= 0 ) THEN
   CALL TS_Abort('Error allocating '//TRIM( Int2LStr( ReKi*NSize/1024**2 ) )//' MB for the Z1Z2 coherence array.')
ENDIF



!! Initialize the velocity matrix  !V(:,:,:) = 0.0 - done outside the subroutine

!--------------------------------------------------------------------------------
! Calculate the distances and other parameters that don't change with frequency
!---------------------------------------------------------------------------------

IF ( SpecModel == SpecModel_IECKAI .or. SpecModel == SpecModel_IECVKM .OR. SpecModel == SpecModel_MODVKM .OR. SpecModel == SpecModel_API) THEN
   InCohB(:)    = 0.12/LC

ENDIF

         II = 0
POINT_I: DO IZ=1,p_grid%ZLim   !NumGrid_Z
            DO IY=1,p_grid%IYmax(IZ) !NumGrid_Y

               II = II + 1                            ! Index of point I: S(I)  !equivalent to II = ( IZ - 1 )*NumGrid_Y + IY
               IF (II > p_grid%NPoints) EXIT POINT_I            ! Don't go past the end of the array; this exits the IZ loop

               JJ = 0
POINT_J:       DO JZ=1,IZ
                  DO JY=1,p_grid%IYmax(JZ) !NumGrid_Y

                     JJ = JJ + 1                      ! Index of point J: S(J)  !equivalent to JJ = ( JZ - 1 )*NumGrid_Y + JY

                     IF ( JJ > II )  EXIT POINT_J     ! The coherence matrix is symmetric

                     IF ( IZ > p_grid%NumGrid_Z ) THEN       ! Get the correct location if we're using an extra point for the hub
                        I = p_grid%YLim
                        IF ( JZ > p_grid%NumGrid_Z ) THEN
                           J = p_grid%YLim
                        ELSE
                           J = JY
                        ENDIF
                     ELSE
                        I = IY
                        J = JY
                     ENDIF

                     JJ1       = JJ - 1
                     Indx      = p_grid%NPoints*JJ1 - JJ*JJ1/2 + II   !Index of matrix ExCoDW (now Matrix), coherence between points I & J

                         Dist_Y(Indx)= ABS( p_grid%Y(I) - p_grid%Y(J) )

                         Dist_Z(Indx)= ABS( p_grid%Z(IZ) - p_grid%Z(JZ) )
!mlb                         Dist_Z12(Indx)=ABS(Z(IZ)*Z(JZ))
                         Z1Z2(Indx) = p_grid%Z(IZ)*p_grid%Z(JZ)



                  ENDDO    ! JY
            ENDDO POINT_J  ! JZ

      ENDDO             ! IY
   ENDDO POINT_I        ! IZ

!IF ( COH_OUT ) THEN
IF ( .TRUE. ) THEN

      ! Write the coherence for three frequencies, for debugging purposes

      CALL OpenFOutFile( UC, TRIM(RootName)//'.coh' )
!      WRITE( UC, '(A4,X,3(A16),2(A11))' ) 'Comp','Distance', 'Avg Height', 'Avg Wind Speed', 'c(F1)','c(F2)'

      WRITE( UC, '(A4,X,A16,1X,'//INT2LSTR(NSize)//'(G10.4,1X))' ) 'Comp','Freq',(I,I=1,NSize)
!MLB: Dist is never set or, apparently, used.
!mlb      WRITE( UC,   '(5X,A16,1X,'//INT2LSTR(NSize)//'(G10.4,1X))' ) 'Distance',     Dist(:)
      WRITE( UC,   '(5X,A16,1X,'//INT2LSTR(NSize)//'(G10.4,1X))' ) 'Distance_Y',   Dist_Y(:)
      WRITE( UC,   '(5X,A16,1X,'//INT2LSTR(NSize)//'(G10.4,1X))' ) 'Distance_Z', Dist_Z(:)
!mlb      WRITE( UC,   '(5X,A16,1X,'//INT2LSTR(NSize)//'(G10.4,1X))' ) 'Distance_Z12', Dist_Z12(:)
      WRITE( UC,   '(5X,A16,1X,'//INT2LSTR(NSize)//'(G10.4,1X))' ) 'Z(IZ)*Z(JZ)', Z1Z2(:)
ENDIF

DO IVec = 1,3

   IF ( ( SpecModel == SpecModel_IECKAI .or. SpecModel == SpecModel_IECVKM .OR. SpecModel == SpecModel_MODVKM .OR. SpecModel == SpecModel_API ) .AND. ( IVec /= 1 ) ) THEN
         ! There is no coherence defined for the v or w component of the IEC spectral models
      IdentityCoh = .TRUE.
   ELSE
      IdentityCoh = .FALSE.
   ENDIF

   CALL WrScr ( '    '//Comp(IVec)//'-component matrices' )

   !--------------------------------------------------------------------------------
   ! Calculate the coherence, Veers' H matrix (CSDs), and the fourier coefficients
   !---------------------------------------------------------------------------------

   DO IFREQ = 1,p_grid%NumFreq

      !---------------------------------------------------
      ! Calculate the coherence and Veers' H matrix (CSDs)
      !---------------------------------------------------
      
      IF (.NOT. IdentityCoh) THEN

            ! -----------------------------------------------
            ! Create the coherence matrix for this frequency
            ! -----------------------------------------------

         DO II=1,p_grid%NPoints
            DO JJ=1,II

                  JJ1       = JJ - 1
                  Indx      = p_grid%NPoints*JJ1 - JJ*JJ1/2 + II   !Index of matrix ExCoDW (now Matrix), coherence between points I & J

!mlb: THis is where to look for the error.

!mlb                  TEMP_Y=Coef_AlphaY*p_grid%Freq(IFreq)**Coef_RY*(Dist_Y(Indx)/Coef_1)**Coef_QY*(Dist_Z12(Indx)/Coef_2)**(-0.5*Coef_PY)
!mlb                  TEMP_Z=Coef_AlphaZ*p_grid%Freq(IFreq)**Coef_RZ*(Dist_Z(Indx)/Coef_1)**Coef_QZ*(Dist_Z12(Indx)/Coef_2)**(-0.5*Coef_PZ)
                  TEMP_Y=Coef_AlphaY*p_grid%Freq(IFreq)**Coef_RY*(Dist_Y(Indx)/Coef_1)**Coef_QY*(Z1Z2(Indx)/Coef_2)**(-0.5*Coef_PY)
                  TEMP_Z=Coef_AlphaZ*p_grid%Freq(IFreq)**Coef_RZ*(Dist_Z(Indx)/Coef_1)**Coef_QZ*(Z1Z2(Indx)/Coef_2)**(-0.5*Coef_PZ)

!mlb                  TRH(Indx)=EXP(-Coef_1*SQRT(TEMP_Y**2+TEMP_Z**2)/U0_1HR)
                  TRH(Indx)=EXP(-Coef_1*SQRT(TEMP_Y**2+TEMP_Z**2)/URef)


            ENDDO ! JJ
         ENDDO ! II
      END IF
      
         
      CALL Coh2Coeffs( IdentityCoh, IVec, IFreq, TRH, S, PhaseAngles, V, NSize )
         

   ENDDO !IFreq

ENDDO !IVec

IF (COH_OUT)  CLOSE( UC )

!mlb IF ( ALLOCATED( Dist      ) ) DEALLOCATE( Dist      )
IF ( ALLOCATED( Dist_Y    ) ) DEALLOCATE( Dist_Y     )
IF ( ALLOCATED( Dist_Z ) ) DEALLOCATE( Dist_Z )
!mlb IF ( ALLOCATED( Dist_Z12 ) ) DEALLOCATE( Dist_Z12 )
IF ( ALLOCATED( Z1Z2 ) ) DEALLOCATE( Z1Z2 )
CALL WrScr( ' Two-dimensional API coherence matrix is generated!' )
RETURN

END SUBROUTINE CohSpecVMat_API
!=======================================================================
SUBROUTINE Coh2Coeffs( IdentityCoh, IVec, IFreq, TRH, S, PhaseAngles, V, NSize )

use NWTC_LAPACK
USE TSMods

LOGICAL,          INTENT(IN)     :: IdentityCoh
REAL(ReKi),       INTENT(INOUT)  :: TRH         (:)                          ! The transfer function  matrix (NumSteps).
REAL(ReKi),       INTENT(IN)     :: S           (:,:,:)                      ! The turbulence PSD array (NumFreq,NPoints,3).
REAL(ReKi),       INTENT(IN)     :: PhaseAngles (:,:,:)                      ! The array that holds the random phases [number of points, number of frequencies, number of wind components=3].
REAL(ReKi),       INTENT(INOUT)  :: V           (:,:,:)                      ! An array containing the summations of the rows of H (NumSteps,NPoints,3).
INTEGER(IntKi),   INTENT(IN)     :: IVec                                     ! loop counter (=number of wind components)
INTEGER(IntKi),   INTENT(IN)     :: IFreq                                    ! loop counter (=number of frequencies)
INTEGER(IntKi),   INTENT(IN)     :: NSize                                    ! Size of dimension 2 of Matrix


REAL(ReKi)                       :: CPh                                      ! Cosine of the random phase
REAL(ReKi)                       :: SPh                                      ! Sine of the random phase
INTEGER                          :: IF1                                      ! Index to real part of vector
INTEGER                          :: IF2                                      ! Index to complex part of vector

integer                          :: Indx, J, I, Stat
character(1024)                  :: ErrMsg

      IF (IdentityCoh) THEN

            ! -----------------------------------------------------------------------------------
            !  The coherence is the Identity (as is Cholesky); the Veers' H matrix is as follows:
            ! -----------------------------------------------------------------------------------

         Indx = 1
         DO J = 1,p_grid%NPoints ! The column number

               ! The diagonal entries of the matrix:

            TRH(Indx) = SQRT( ABS( S(IFreq,J,IVec) ) )

               ! The off-diagonal values:
            Indx = Indx + 1
            DO I = J+1,p_grid%NPoints ! The row number
               TRH(Indx) = 0.0
               Indx = Indx + 1
            ENDDO ! I
         ENDDO ! J

      ELSE

         IF (COH_OUT) THEN
!            IF (IFreq == 1 .OR. IFreq == p_grid%NumFreq) THEN
               WRITE( UC, '(I3,2X,F15.5,1X,'//INT2LSTR(NSize)//'(G10.4,1X))' ) IVec, p_grid%Freq(IFreq), TRH(1:NSize)
!            ENDIF
         ENDIF

            ! -------------------------------------------------------------
            ! Calculate the Cholesky factorization for the coherence matrix
            ! -------------------------------------------------------------

         CALL LAPACK_pptrf( 'L', p_grid%NPoints, TRH, Stat, ErrMsg )  ! 'L'ower triangular 'TRH' matrix (packed form), of order 'NPoints'; returns Stat

         IF ( Stat /= ErrID_None ) THEN
            CALL WrScr(ErrMsg)
            IF (Stat >= AbortErrLev) &
            CALL TS_Abort('Error '//TRIM(Int2LStr(Stat))//' in the Cholesky factorization occurred at frequency '//&
                           TRIM(Int2LStr(IFreq))//' ('//TRIM(Num2LStr(p_grid%Freq(IFreq)))//' Hz)'//&
                        '. The '//Comp(IVec)//'-component coherence matrix cannot be factored.  '//&
                        'Check the input file for invalid physical properties or modify the coherence exponent '//&
                        'or grid spacing.')
         ENDIF

            ! -------------------------------------------------------------
            ! Create the lower triangular matrix, H, from Veer's method
            ! -------------------------------------------------------------

         Indx = 1
         DO J = 1,p_grid%NPoints  ! Column
            DO I = J,p_grid%NPoints ! Row

                  ! S(IFreq,I,IVec) should never be less than zero, but the ABS makes sure...

               TRH(Indx) = TRH(Indx) * SQRT( ABS( S(IFreq,I,IVec) ) )

               Indx = Indx + 1

            ENDDO !I
         ENDDO !J

      ENDIF !IdentityCoh

      ! -------------------------------------------------------------
      ! Calculate the correlated fourier coefficients.
      ! -------------------------------------------------------------

      IF2      = IFreq*2
      IF1      = IF2 - 1

      DO J=1,p_grid%NPoints

            ! Apply a random phase to each of the columns of H to
            ! produce random phases in the wind component.
            ! Then sum each of the rows into the vector V.
         
         CPh   = COS( PhaseAngles(J,IFreq,IVec) )
         SPh   = SIN( PhaseAngles(J,IFreq,IVec) )

         Indx  = p_grid%NPoints*(J-1) - J*(J-1)/2 + J !Index of H(I,J)
         DO I=J,p_grid%NPoints

            V(IF1,I,IVec) = V(IF1,I,IVec) + TRH(Indx)*CPh  !Real part
            V(IF2,I,IVec) = V(IF2,I,IVec) + TRH(Indx)*SPh  !Imaginary part

            Indx = Indx + 1      !H(I,J)

         ENDDO ! I
      ENDDO ! J

END SUBROUTINE Coh2Coeffs
!=======================================================================

SUBROUTINE GetDefaultCoh(WS,Ht)
   ! These numbers come from Neil's analysis

   USE                     TSMods, ONLY : InCDec
   USE                     TSMods, ONLY : InCohB
   USE                     TSMods, ONLY : RICH_NO
   USE                     TSMods, ONLY : TurbModel, SpecModel
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

      SELECT CASE ( SpecModel )

         CASE ( SpecModel_GP_LLJ )
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


         CASE ( SpecModel_NWTCUP, SpecModel_USRVKM )
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

         CASE ( SpecModel_WF_UPW )
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

         CASE ( SpecModel_WF_07D, SpecModel_WF_14D )
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

         CASE ( SpecModel_USER )
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

   USE                     TSMods, ONLY : p_grid
   USE                     TSMods, ONLY : RefHt                           ! This is needed for the HYDRO spectral models.
   USE                     TSMods, ONLY : RICH_NO
   USE                     TSMods, ONLY : TurbModel, SpecModel
   USE                     TSMods, ONLY : UHub
   USE                     TSMods, ONLY : Ustar
   USE                     TSMods, ONLY : UVskip                                   ! Flag to determine if UV cross-feed term should be skipped or used
   USE                     TSMods, ONLY : UWskip                                   ! Flag to determine if UW cross-feed term should be skipped or used
   USE                     TSMods, ONLY : VWskip                                   ! Flag to determine if VW cross-feed term should be skipped or used
   USE                     TSMods, ONLY : WindProfileType
   USE                     TSMods, ONLY : zL
   
   use tsmods, only: p_RandNum
   use tsmods, only: OtherSt_RandNum
   
use TurbSim_Types
   
   REAL(ReKi), INTENT(INOUT)           :: UW  ! PC_UW
   REAL(ReKi), INTENT(OUT)             :: UV  ! PC_UV
   REAL(ReKi), INTENT(OUT)             :: VW  ! PC_VW
   REAL(ReKi)                          :: rndSgn
   REAL(ReKi)                          :: SignProb
   REAL(ReKi)                          :: Shr
   REAL(ReKi)                          :: Ustar2
   REAL(ReKi)                          :: V(2)
   REAL(ReKi)                          :: Z(2)
   REAL(ReKi)                          :: ZLtmp

   Z(2) = p_grid%HubHt + 0.5*p_grid%RotorDiameter    ! top of the grid
   Z(1) = Z(2) - p_grid%GridHeight     ! bottom of the grid
   V(:) = getWindSpeed(UHub, p_grid%HubHt, Z, p_grid%RotorDiameter, PROFILE_TYPE=WindProfileType)

   Shr = ( V(2)-V(1) ) / p_grid%GridHeight    ! dv/dz

!BJJ: check the ranges of our best-fit parameters, using domains of measured values

   SELECT CASE ( SpecModel )
      CASE ( SpecModel_GP_LLJ )
         ZLtmp  = MIN( MAX( ZL,    REAL(-1.00,ReKi) ), REAL(1.0,ReKi) )  !Limit the observed values of z/L
         UStar2 = MIN( MAX( UStar, REAL( 0.15,ReKi) ), REAL(1.0,ReKi) )  !Limit the observed values of u*
         Ustar2 = Ustar2*Ustar2
      CASE ( SpecModel_NWTCUP )
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

   CALL  RndUnif( p_RandNum, OtherSt_RandNum, rndSgn )
   SELECT CASE ( SpecModel )

      CASE ( SpecModel_GP_LLJ )

         IF (UW <= 0) THEN  !We don't have a local u* value to tie it to; otherwise, assume UW contains magnitude of value we want
            IF ( p_grid%HubHt >= 100.5 ) THEN     ! 116m
               UW =  0.0399 - 0.00371*Uhub - 0.00182*RICH_NO + 0.00251*ZLtmp - 0.402*Shr + 1.033*Ustar2
            ELSEIF ( p_grid%HubHt >= 76.0 ) THEN  ! 85 m
               UW = 0.00668 - 0.00184*Uhub + 0.000709*RICH_NO  + 0.264*Shr + 1.065*Ustar2  !magnitude
            ELSEIF ( p_grid%HubHt >= 60.5 ) THEN  ! 67 m
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

      CASE ( SpecModel_NWTCUP )

         IF ( p_grid%HubHt > 47.0 ) THEN      ! 58m data
            UW = 0.165 - 0.0232*UHub - 0.0129*RICH_NO + 1.337*Ustar2 - 0.758*SHR
         ELSEIF ( p_grid%HubHt >= 26.0 ) THEN ! 37m data
            UW = 0.00279 - 0.00139*UHub + 1.074*Ustar2 + 0.179*SHR
         ELSE                          ! 15m data
            UW = -0.1310 + 0.0239*UHub + 0.556*Ustar2
         ENDIF
         UW = MAX(UW,0.0)

         IF (UW > 0) THEN !i.e. not equal to zero
            SignProb = 0.765 + 0.57/PI * ATAN( 0.88356*LOG(UW)+2.47668)
            IF (rndSgn <= SignProb) UW = -UW
         ENDIF

      CASE ( SpecModel_WF_14D )

         UW = -Ustar2
         IF ( rndSgn > 0.9937 )  UW = -UW

      CASE ( SpecModel_USER )
         UW = 0.0
         UWskip = .true.

      CASE ( SpecModel_TIDAL, SpecModel_RIVER ) ! HYDROTURBSIM specific.
         UW = -Ustar2*(1-p_grid%HubHt/RefHt)
      CASE DEFAULT

         UW = -Ustar2

   END SELECT

   !-------------------------------------------------------------------------------------------------
   ! default UV Reynolds stress
   !-------------------------------------------------------------------------------------------------

   CALL  RndUnif( p_RandNum, OtherSt_RandNum, rndSgn )
   SELECT CASE ( SpecModel )

      CASE ( SpecModel_GP_LLJ )

         IF ( p_grid%HubHt >= 100.5 ) THEN     ! 116m
            UV = 0.199 - 0.0167*Uhub + 0.0115*ZLtmp + 1.143*Ustar2
            UV = MAX(UV,0.0)
            IF ( rndSgn < 0.6527 ) UV = -UV
         ELSEIF ( p_grid%HubHt >= 76.0 ) THEN  ! 85 m
            UV = 0.190 - 0.0156*Uhub + 0.00931*ZLtmp + 1.101*Ustar2
            UV = MAX(UV,0.0)
            IF ( rndSgn < 0.6394 ) UV = -UV
         ELSEIF ( p_grid%HubHt >= 60.5 ) THEN  ! 67 m
            UV = 0.178 - 0.0141*Uhub + 0.00709*ZLtmp + 1.072*Ustar2
            UV = MAX(UV,0.0)
            IF ( rndSgn < 0.6326 ) UV = -UV
         ELSE                           ! 54 m
            UV = 0.162 - 0.0123*Uhub + 0.00784*RICH_NO + 1.024*Ustar2
            UV = MAX(UV,0.0)
            IF ( rndSgn < 0.6191 ) UV = -UV
         ENDIF

      CASE ( SpecModel_NWTCUP )

            ! Get the magnitude and add the sign
         IF ( p_grid%HubHt > 47.0 ) THEN      ! 58m data
            UV = 0.669 - 0.0300*UHub - 0.0911*RICH_NO + 1.421*Ustar2 - 1.393*SHR
         ELSEIF ( p_grid%HubHt >= 26.0 ) THEN ! 37m data
            UV = 1.521 - 0.00635*UHub - 0.2200*RICH_NO + 3.214*Ustar2 - 3.858*SHR
         ELSE                          ! 15m data
            UV = 0.462 - 0.01400*UHub + 1.277*Ustar2
         ENDIF
         UV = MAX(UV,0.0)
         IF (UV > 0) THEN !i.e. not equal to zero
            SignProb = 0.33 + 0.64/PI * ATAN( -0.374775*LOG(UV)-0.205681)
            IF (rndSgn <= SignProb) UV = -UV
         ENDIF

      CASE ( SpecModel_WF_UPW )

         UV = 0.0202 + 0.890*Ustar2 - 2.461*Shr
         UV = MAX(UV,0.0)
         IF ( rndSgn < 0.7315 ) UV = -UV

      CASE ( SpecModel_WF_07D )

         UV = 0.5040 + 0.177*Ustar2
         UV = MAX(UV,0.0)
         IF ( rndSgn < 0.7355 ) UV = -UV

      CASE ( SpecModel_WF_14D )

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

   CALL  RndUnif( p_RandNum, OtherSt_RandNum, rndSgn )
   SELECT CASE ( SpecModel )

      CASE ( SpecModel_GP_LLJ )

         IF ( p_grid%HubHt >= 100.5 ) THEN     ! 116m
            VW =  0.0528  - 0.00210*Uhub - 0.00531*RICH_NO - 0.519*Shr + 0.283*Ustar2
            VW = MAX(VW,0.0)
            IF ( rndSgn < 0.2999 ) VW = -VW
         ELSEIF ( p_grid%HubHt >= 76.0 ) THEN  ! 85 m
            VW =  0.0482  - 0.00264*Uhub - 0.00391*RICH_NO - 0.240*Shr + 0.265*Ustar2
            VW = MAX(VW,0.0)
            IF ( rndSgn < 0.3061 ) VW = -VW
         ELSEIF ( p_grid%HubHt >= 60.5 ) THEN  ! 67 m
            VW =  0.0444  - 0.00249*Uhub - 0.00403*RICH_NO - 0.141*Shr + 0.250*Ustar2
            VW = MAX(VW,0.0)
            IF ( rndSgn < 0.3041 ) VW = -VW
         ELSE                           ! 54 m
            VW =  0.0443  - 0.00261*Uhub - 0.00371*RICH_NO - 0.107*Shr + 0.226*Ustar2
            VW = MAX(VW,0.0)
            IF ( rndSgn < 0.3111 ) VW = -VW
         ENDIF

      CASE ( SpecModel_NWTCUP )

         IF ( p_grid%HubHt > 47.0 ) THEN      ! 58m data
            VW = 0.174 + 0.00154*UHub - 0.0270*RICH_NO + 0.380*Ustar2 - 1.131*Shr - 0.00741*ZLtmp
         ELSEIF ( p_grid%HubHt >= 26.0 ) THEN ! 37m data
            VW = 0.120 + 0.00283*UHub - 0.0227*RICH_NO + 0.306*Ustar2 - 0.825*Shr
         ELSE                          ! 15m data
            VW = 0.0165 + 0.00833*UHub                 + 0.224*Ustar2
         ENDIF
         VW = MAX(VW,0.0)
         IF (VW > 0) THEN !i.e. not equal to zero
            SignProb = 0.725 + 0.65/PI * ATAN( 0.654886*LOG(VW)+1.777198)
            IF (rndSgn <= SignProb) VW = -VW
         ENDIF

      CASE ( SpecModel_WF_UPW )

         VW = 0.0263 + 0.273*Ustar2 - 0.684*Shr
         VW = MAX(VW,0.0)
         IF ( rndSgn < 0.3139 ) VW = -VW

      CASE ( SpecModel_WF_07D )

         VW = 0.241 + 0.118*Ustar2
         VW = MAX(VW,0.0)
         IF ( rndSgn < 0.0982 ) VW = -VW

      CASE ( SpecModel_WF_14D )

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
FUNCTION getUstarDiab(u_ref, z_ref)

   USE                                  TSMods, ONLY: z0             ! Surface roughness length -- It must be > 0 (which we've already checked for)
   USE                                  TSMods, ONLY: ZJetMax        ! Height of the low-level jet
   USE                                  TSMods, ONLY: ZL             ! M-O stability parameter

   IMPLICIT                              NONE

   REAL(ReKi),   INTENT(IN)           :: u_ref                       ! Wind speed at reference height
   REAL(ReKi),   INTENT(IN)           :: z_ref                       ! Reference height

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

   getUstarDiab = ( 0.4 * u_ref ) / ( LOG( z_ref / z0 ) - psiM )

END FUNCTION
!=======================================================================


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

SELECT CASE ( SpecModel )

   CASE (SpecModel_WF_UPW, SpecModel_NWTCUP)
      IF ( Ri_No > 0.0 ) THEN
         PowerLawExp = 0.14733
      ELSE
         PowerLawExp = 0.087687698 + 0.059641545*EXP(Ri_No/0.04717783)
      ENDIF

   CASE ( SpecModel_WF_07D, SpecModel_WF_14D )
      IF ( Ri_No > 0.04 ) THEN
         PowerLawExp = 0.17903
      ELSE
         PowerLawExp = 0.127704032 + 0.031228952*EXP(Ri_No/0.0805173)
      ENDIF

   CASE (SpecModel_SMOOTH, SpecModel_GP_LLJ, SpecModel_TIDAL, SpecModel_RIVER)
      ! A 1/7 power law seems to work ok for HYDRO spectral models also...
      PowerLawExp = 0.143

   CASE DEFAULT
      IF ( p_IEC%IEC_WindType == IEC_EWM1 .OR. p_IEC%IEC_WindType == IEC_EWM50 ) THEN
         PowerLawExp = 0.11         ! [IEC 61400-1 6.3.2.1 (14)]
      ELSEIF ( p_IEC%IECstandard == 3 ) THEN
         PowerLawExp = 0.14         ! [IEC 61400-3 Page 22 (3)]
      ELSE
         PowerLawExp = 0.2          ! [IEC 61400-1 6.3.1.2 (10)]
      ENDIF

END SELECT

RETURN
END FUNCTION PowerLawExp
!=======================================================================
SUBROUTINE PSDcal ( Ht, Ucmp, UShr, ZL_loc, Ustar_loc)

   ! This routine calls the appropriate spectral model.

USE                                TSMods

IMPLICIT                            NONE


REAL(ReKi), INTENT(IN)           :: Ht             ! Height
REAL(ReKi), INTENT(IN)           :: Ucmp           ! Velocity
REAL(ReKi), INTENT(IN), OPTIONAL :: UShr           ! Shear (du/dz)
REAL(ReKi), INTENT(IN), OPTIONAL :: Ustar_loc      ! Local ustar
REAL(ReKi), INTENT(IN), OPTIONAL :: ZL_loc         ! Local Z/L

SELECT CASE ( SpecModel )
   CASE ( SpecModel_GP_LLJ )
      IF ( PRESENT(Ustar_loc) .AND. PRESENT(ZL_loc) ) THEN
         CALL Spec_GPLLJ   ( Ht, Ucmp, ZL_loc, Ustar_loc, Work )
      ELSE
         CALL Spec_GPLLJ   ( Ht, Ucmp, ZL,     Ustar,     Work )
      ENDIF
   CASE ( SpecModel_IECKAI )
      CALL Spec_IECKAI  ( p_IEC%SigmaIEC, p_IEC%IntegralScale, Work )
   CASE ( SpecModel_IECVKM )
      CALL Spec_IECVKM  ( p_IEC%SigmaIEC(1), p_IEC%IntegralScale, Work )
   CASE ( SpecModel_API )
      CALL Spec_API ( Ht, Work )
   CASE (SpecModel_NWTCUP)
      CALL Spec_NWTCUP  ( Ht, Ucmp, Work )
   CASE ( SpecModel_SMOOTH )
      CALL Spec_SMOOTH   ( Ht, Ucmp, Work )
   CASE ( SpecModel_TIDAL, SpecModel_RIVER )
      CALL Spec_TIDAL  ( Ht, UShr, Work, SpecModel )
   CASE ( SpecModel_USER )
      CALL Spec_UserSpec   ( Work )
   CASE ( SpecModel_USRVKM )
      CALL Spec_vonKrmn   ( Ht, Ucmp, Work )
   CASE (SpecModel_WF_UPW)
      CALL Spec_WF_UPW  ( Ht, Ucmp, Work )
   CASE ( SpecModel_WF_07D, SpecModel_WF_14D )
      CALL Spec_WF_DW ( Ht, Ucmp, Work )
   CASE ( SpecModel_NONE )
      Work(:,:) = 0.0
!bjj TEST: CALL Spec_Test ( Ht, Ucmp, Work )

   CASE ( SpecModel_MODVKM )
      IF (MVK) THEN
 !        CALL Mod_vKrm( Ht, Ucmp, Work )
      ELSE
         CALL TS_Abort ( 'Specified turbulence PSD, "'//TRIM( TurbModel )//'", not availible.' )
      ENDIF

   CASE DEFAULT
      CALL TS_Abort ( 'Specified turbulence PSD, "'//TRIM( TurbModel )//'", not availible.' )
END SELECT



RETURN
END SUBROUTINE PSDcal
!=======================================================================




SUBROUTINE CreateGrid( p_grid, NumGrid_Y2, NumGrid_Z2, NSize, ErrStat, ErrMsg )

! Assumes that these variables are set:
!  GridHeight
!  GridWidth 
!  NumGrid_Y
!  NumGrid_Z

use TSMods
use FFT_Module, only: PSF

   TYPE(TurbSim_GridParameterType), INTENT(INOUT) :: p_grid
   INTEGER                        , INTENT(  OUT) :: NumGrid_Y2                      ! Y Index of the hub (or the nearest point left the hub if hub does not fall on the grid)
   INTEGER                        , INTENT(  OUT) :: NumGrid_Z2                      ! Z Index of the hub (or the nearest point below the hub if hub does not fall on the grid) 
   INTEGER                        , INTENT(  OUT) :: NSize                           ! Size of the spectral matrix at each frequency
   INTEGER(IntKi),                  intent(  out) :: ErrStat                         ! Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg                          ! Message describing error
   
   ! local variables:
   REAL(ReKi)                                     :: TmpReal                         ! A temporary variable holding a Z (height) value, and other misc. variables
   REAL(ReKi)                                     :: DelF                            ! Delta frequency
   INTEGER                                        :: IY, IZ, IFreq                   ! loop counters 
   INTEGER                                        :: FirstTwrPt                      ! Z index of first tower point
   INTEGER                                        :: NumSteps4                       ! one-fourth the number of steps
   
   INTEGER(IntKi)                                 :: ErrStat2                         ! Error level (local)
   CHARACTER(MaxMsgLen)                           :: ErrMsg2                          ! Message describing error (local)
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   ! First, let's deal with the frequencies:
   
      ! Calculate Total time and NumSteps.
      ! Find the product of small factors that is larger than NumSteps (prime #9 = 23).
      ! Make sure it is a multiple of 4 too.

   IF ( p_grid%Periodic ) THEN
      p_grid%NumSteps    = CEILING( p_grid%AnalysisTime / p_grid%TimeStep )
      NumSteps4   = ( p_grid%NumSteps - 1 )/4 + 1
      p_grid%NumSteps    = 4*PSF( NumSteps4 , 9 )  ! >= 4*NumSteps4 = NumOutSteps + 3 - MOD(NumOutSteps-1,4) >= NumOutSteps
      p_grid%NumOutSteps = p_grid%NumSteps
   ELSE
      p_grid%NumOutSteps = CEILING( ( p_grid%UsableTime + p_grid%GridWidth/UHub )/p_grid%TimeStep )
      p_grid%NumSteps    = MAX( CEILING( p_grid%AnalysisTime / p_grid%TimeStep ), p_grid%NumOutSteps )
      NumSteps4   = ( p_grid%NumSteps - 1 )/4 + 1
      p_grid%NumSteps    = 4*PSF( NumSteps4 , 9 )  ! >= 4*NumSteps4 = NumOutSteps + 3 - MOD(NumOutSteps-1,4) >= NumOutSteps
   END IF


   p_grid%NumFreq     = p_grid%NumSteps/2
   DelF        = 1.0/( p_grid%NumSteps*p_grid%TimeStep )
      
   
   CALL AllocAry( p_grid%Freq, p_grid%NumFreq, 'Freq (frequency array)', ErrStat2, ErrMsg2)
   CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'CreateGrid')
   IF (ErrStat >= AbortErrLev) RETURN

   DO IFreq=1,p_grid%NumFreq
      p_grid%Freq(IFreq) = IFreq*DelF
   ENDDO 
   
   
   ! Then, let's deal with the rectangular grid:
   
   
   p_grid%GridRes_Y = p_grid%GridWidth  / REAL( p_grid%NumGrid_Y - 1, ReKi )
   p_grid%GridRes_Z = p_grid%GridHeight / REAL( p_grid%NumGrid_Z - 1, ReKi )      

   p_grid%Zbottom = p_grid%HubHt + 0.5*p_grid%RotorDiameter                             ! height of the highest grid points
   p_grid%Zbottom = p_grid%Zbottom - p_grid%GridRes_Z * REAL(p_grid%NumGrid_Z - 1, ReKi)   ! height of the lowest grid points

   IF ( p_grid%Zbottom <= 0.0 ) THEN
      CALL TS_Abort ( 'The lowest grid point ('//TRIM(Num2LStr(p_grid%Zbottom))// ' m) must be above the ground. '//&
                      'Adjust the appropriate values in the input file.' )
   ENDIF

   NumGrid_Y2 = INT( ( p_grid%NumGrid_Y + 1 ) / 2 )                    ! These are the hub indicies, unless the hub is an extra point
   NumGrid_Z2 = INT( Tolerance + ( p_grid%HubHt - p_grid%Zbottom ) / p_grid%GridRes_Z ) + 1 
   
   p_grid%ExtraTwrPT = .FALSE.

   IF ( MOD(p_grid%NumGrid_Y, 2) == 0 ) THEN
      p_grid%ExtraTwrPT = .TRUE.
      p_grid%ExtraHubPT = .TRUE.
   ELSEIF ( ABS((NumGrid_Z2-1)*p_grid%GridRes_Z + p_grid%Zbottom - p_grid%HubHt) > Tolerance ) THEN
      p_grid%ExtraHubPT = .TRUE.
   ELSE
      p_grid%ExtraHubPT = .FALSE.
   ENDIF

   p_grid%NPoints    = p_grid%NumGrid_Y*p_grid%NumGrid_Z                ! Number of points in the regular grid

   
   ! Then, let's deal with the hub point:
   
   IF (p_grid%ExtraHubPT) THEN
      p_grid%NPoints = p_grid%NPoints + 1                               ! Add the hub point if necessary
      p_grid%ZLim = p_grid%NumGrid_Z+1
      p_grid%YLim = p_grid%NumGrid_Y+1
   ELSE
      p_grid%ZLim = p_grid%NumGrid_Z
      p_grid%YLim = p_grid%NumGrid_Y
   ENDIF

   
   ! Finally, let's deal with the tower "lollipop" points:
   
   IF ( WrFile(FileExt_TWR) ) THEN

         ! Compute the number of points between the bottom of the grid and the ground 
         ! ( but we don't want to be on the ground, just more than "Tolerance" from it )
 
      IZ = INT( ( p_grid%Zbottom - Tolerance ) / p_grid%GridRes_Z )

      IF ( p_grid%ExtraTwrPT ) THEN 
         IZ = IZ + 1  ! Let's add the point on the bottom of the grid so tower interpolation is easier in AeroDyn
      ENDIF

      IF ( IZ > 0 ) THEN

         p_grid%ZLim = p_grid%ZLim + IZ
         p_grid%NPoints = p_grid%NPoints + IZ                       ! Add the number of tower points

         IF ( .NOT. p_grid%ExtraHubPT ) THEN
            p_grid%YLim = p_grid%YLim + 1
         ENDIF

      ELSE

         CALL TS_Warn ( ' There are no extra tower data points below the grid. Tower output will be turned off.', US)  

         WrFile(FileExt_TWR) = .FALSE.

      ENDIF

   ENDIF

   NSize   = p_grid%NPoints*( p_grid%NPoints + 1 )/2   
   
   
   CALL AllocAry(p_grid%Z,     p_grid%ZLim, 'Z (vertical locations of the grid points)',   ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'CreateGrid')
   CALL AllocAry(p_grid%Y,     p_grid%YLim, 'Y (lateral locations of the grid points)',    ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'CreateGrid')
   CALL AllocAry(p_grid%IYmax, p_grid%ZLim, 'IYmax (horizontal locations at each height)', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'CreateGrid') !  Allocate the array of vertical locations of the grid points.

   IF (ErrStat >= AbortErrLev) RETURN
   
      ! Initialize cartesian Y,Z values of the grid.

   DO IY = 1,p_grid%NumGrid_Y
      p_grid%Y(IY)    = -0.5*p_grid%GridWidth  + p_grid%GridRes_Y*( IY - 1 )
   ENDDO

   DO IZ = 1,p_grid%NumGrid_Z
      p_grid%Z(IZ)     = p_grid%Zbottom + p_grid%GridRes_Z*( IZ - 1 )
      p_grid%IYmax(IZ) = p_grid%NumGrid_Y           ! Number of lateral points at this height
   ENDDO

   IF (p_grid%ExtraHubPT) THEN

      p_grid%Y(p_grid%NumGrid_Y+1)     = 0.0
      p_grid%Z(p_grid%NumGrid_Z+1)     = p_grid%HubHt
      p_grid%IYmax(p_grid%NumGrid_Z+1) = 1

      p_grid%HubIndx = p_grid%NumGrid_Y*p_grid%NumGrid_Z + 1

      FirstTwrPt = p_grid%NumGrid_Z + 2              ! The start of tower points, if they exist

   ELSE

      p_grid%HubIndx = p_grid%NumGrid_Y*(NumGrid_Z2-1) + NumGrid_Y2

      FirstTwrPt = p_grid%NumGrid_Z + 1              ! The start of tower points, if they exist

   ENDIF

   IF ( WrFile(FileExt_TWR) ) THEN

      p_grid%Y(p_grid%YLim) = 0.0
   
      TmpReal  = p_grid%Z(1)

      IF ( p_grid%ExtraTwrPT ) THEN 
         TmpReal = TmpReal + p_grid%GridRes_Z
      ENDIF

      DO IZ = FirstTwrPt,p_grid%ZLim
         p_grid%Z(IZ)     = TmpReal - p_grid%GridRes_Z
         TmpReal   = p_grid%Z(IZ)

         p_grid%IYmax(IZ) = 1                 ! The number of lateral points at this height
      ENDDO

   ENDIF   
         
END SUBROUTINE CreateGrid

!=======================================================================
!> This routine calculates the wind components in the Inertial reference
!!  frame.
SUBROUTINE CalculateWindComponents(v, ubar, HFlowAng, VFlowAng, V_Inertial, UH, UT)

   REAL(ReKi), INTENT(IN)           :: v(3)               !< u,v,w components (streamwise)
   REAL(ReKi), INTENT(IN)           :: ubar               !< mean streamwise component
   REAL(ReKi), INTENT(IN)           :: HFlowAng           !< horizontal flow angle
   REAL(ReKi), INTENT(IN)           :: VFlowAng           !< vertical flow angle
                                    
   REAL(ReKi), INTENT(OUT)          :: V_Inertial(3)      !< U,V,W components (inertial)
   REAL(ReKi), INTENT(OUT),OPTIONAL :: UH                 !< horizontal wind speed (U+V components)
   REAL(ReKi), INTENT(OUT),OPTIONAL :: UT                 !< total wind speed (U+V+W components)
   
   
   
   ! Local variables
   REAL(ReKi)              :: UTmp                            ! The instantaneous u-component wind speed at the hub
   REAL(ReKi)              :: UHTmp2                          ! The square of the instantaneous horizontal wind speed at the hub
   REAL(ReKi)              :: V_Inertial2(3)                  ! the U,V,W components (inertial) squared
   
   REAL(ReKi)              :: CVFA                            ! Cosine of the vertical flow angle
   REAL(ReKi)              :: SVFA                            ! Sine of the vertical flow angle
   REAL(ReKi)              :: CHFA                            ! Cosine of the horizontal flow angle
   REAL(ReKi)              :: SHFA                            ! Sine of the horizontal flow angle


   CHFA = COS( HFlowAng*D2R )
   SHFA = SIN( HFlowAng*D2R )

   CVFA = COS( VFlowAng*D2R )
   SVFA = SIN( VFlowAng*D2R )

   

      ! Calculate longitudinal (UTmp) value for point,
      ! as well as rotated (V_Inertial) 
      ! components applying specified flow angles.

      ! Add mean wind speed to the streamwise component
   UTmp = v(1) + ubar
   
      ! Rotate the wind components from streamwise orientation to the X-Y-Z grid       
   V_Inertial(1) = UTmp*CHFA*CVFA - v(2)*SHFA - v(3)*CHFA*SVFA
   V_Inertial(2) = UTmp*SHFA*CVFA + v(2)*CHFA - v(3)*SHFA*SVFA  
   V_Inertial(3) = UTmp*SVFA                  + v(3)*CVFA

   IF ( PRESENT( UH ) .OR. PRESENT( UT ) ) THEN
         ! Calculate hub horizontal wind speed (UHTmp) and Total wind speed (UTTmp)
   
      V_Inertial2 = V_Inertial*V_Inertial             !inertial frame coordinates   
      UHTmp2      = V_Inertial2(1) + V_Inertial2(2)   !inertial frame coordinates

      IF ( PRESENT( UH ) ) UH = SQRT( UHTmp2 )        !inertial frame coordinates
      IF ( PRESENT( UT ) ) UT = SQRT( UHTmp2 + V_Inertial2(3) )
   END IF

END SUBROUTINE CalculateWindComponents
!=======================================================================
! This routine calculates the instantaneous Reynolds stresses, including TKE and CTKE
SUBROUTINE CalculateStresses(v, uv, uw, vw, TKE, CTKE )
   REAL(ReKi), INTENT(IN)  :: v(3)               !< u,v,w components (streamwise, zero-mean)
   
   REAL(ReKi), INTENT(OUT) ::  uv                !< The instantaneous u'v' Reynolds stress at the hub
   REAL(ReKi), INTENT(OUT) ::  uw                !< The instantaneous u'w' Reynolds stress at the hub
   REAL(ReKi), INTENT(OUT) ::  vw                !< The instantaneous v'w' Reynolds stress at the hub
   REAL(ReKi), INTENT(OUT) ::  TKE               !< The instantaneous TKE at the hub
   REAL(ReKi), INTENT(OUT) ::  CTKE              !< The instantaneous CTKE the hub
   
   
   uv = v(1)*v(2)
   uw = v(1)*v(3)
   vw = v(2)*v(3)
   
   TKE  = 0.5*(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
   CTKE = 0.5*SQRT(uv*uv + uw*uw + vw*vw)
   
END SUBROUTINE
!=======================================================================
SUBROUTINE TimeSeriesScaling_IEC(p_grid, V, ScaleIEC, TargetSigma, ActualSigma, HubFactor)


   TYPE(TurbSim_GridParameterType), INTENT(IN)     ::  p_grid                          ! parameters defining TurbSim's grid (including time/frequency)
   REAL(ReKi),                      INTENT(INOUT)  ::  V(:,:,:)                        ! velocity 
   REAL(ReKi),                      INTENT(IN)     ::  TargetSigma(3)                  ! target standard deviations
   REAL(ReKi),                      INTENT(  OUT)  ::  ActualSigma(3)                  ! actual standard deviations
   REAL(ReKi),                      INTENT(  OUT)  ::  HubFactor(3)                    ! factor used to scale standard deviations at the hub point
   INTEGER(IntKi),                  INTENT(In   )  ::  ScaleIEC                        ! scaling method: 1=use only 1 factor; 2=scale each individual point
   
   REAL(DbKi)                                      ::  CGridSum                        ! The sums of the velocity components at the points surrounding the hub (or at the hub if it's on the grid)
   REAL(DbKi)                                      ::  CGridSum2                       ! The sums of the squared velocity components at the points surrouding the hub 
   REAL(ReKi)                                      ::  UGridMean                       ! Average wind speed at a point 
   REAL(ReKi)                                      ::  UGridSig                        ! Standard deviation of the wind speed at a point 
   INTEGER(IntKi)                                  ::  IT                              ! loop counter (time)
   INTEGER(IntKi)                                  ::  Indx                            ! loop counter (grid point)
   INTEGER(IntKi)                                  ::  IVec                            ! loop counter (wind component)
   
   
   DO IVec = 1,3
      CGridSum  = 0.0
      CGridSum2 = 0.0
                           
      DO IT=1,p_grid%NumSteps !BJJ: NumOutSteps  -- scale to the output value?
         CGridSum  = CGridSum  + V( IT, p_grid%HubIndx, IVec )
         CGridSum2 = CGridSum2 + V( IT, p_grid%HubIndx, IVec )* V( IT, p_grid%HubIndx, IVec )
      ENDDO ! IT
               
      UGridMean          = CGridSum/p_grid%NumSteps !BJJ: NumOutSteps  -- scale to the output value?
      ActualSigma(IVec)  = SQRT( ABS( (CGridSum2/p_grid%NumSteps) - UGridMean*UGridMean ) )
            

      HubFactor(IVec) = TargetSigma(IVec)/ActualSigma(IVec)
                  
      IF (ScaleIEC == 1 .OR. IVec > 1) THEN ! v and w have no coherence, thus all points have same std, so we'll save some calculations
               
         V(:,:,IVec) =     HubFactor(IVec) * V(:,:,IVec)
               
      ELSE  ! Scale each point individually
               
         DO Indx = 1,p_grid%NPoints             
            CGridSum  = 0.0
            CGridSum2 = 0.0
                                    
            DO IT=1,p_grid%NumSteps !BJJ: NumOutSteps  -- scale to the output value?
               CGridSum  = CGridSum  + V( IT, Indx, IVec )
               CGridSum2 = CGridSum2 + V( IT, Indx, IVec )* V( IT, Indx, IVec )
            ENDDO ! IT
                  
            UGridMean = CGridSum/p_grid%NumSteps !BJJ: NumOutSteps  -- scale to the output value?
            UGridSig  = SQRT( ABS( (CGridSum2/p_grid%NumSteps) - UGridMean*UGridMean ) )
            
            V(:,Indx,IVec) = (TargetSigma(IVec) / UGridSig) * V(:,Indx,IVec)   
         ENDDO ! Indx
               
      ENDIF            

   ENDDO !IVec   

END SUBROUTINE TimeSeriesScaling_IEC
!=======================================================================
SUBROUTINE CalcIECScalingParams( p_IEC )
! REQUires these be set prior to calling:NumTurbInp, IECedition, IECTurbC, IEC_WindType
! calculates SigmaIEC, Lambda, IntegralScale, Lc
use TurbSim_Types
use TSMods, only: UHub, p_grid, Fc, z0, SpecModel

   TYPE(IEC_ParameterType), INTENT(INOUT) :: p_IEC                       ! parameters for IEC models
      
   IF ( p_IEC%NumTurbInp )  THEN
   
         ! user specified a particular percent TI:
         
      p_IEC%TurbInt     = 0.01*p_IEC%PerTurbInt
      p_IEC%SigmaIEC(1) = p_IEC%TurbInt*UHub
      
      ! bjj: note Vave isn't set in this case, but we only print it to the summary file (and use it) if .not. NumTurbInp      

   ELSE

      
      SELECT CASE (p_IEC%IECedition)

         CASE ( 2 )

            IF ( p_IEC%IECTurbC  == 'A' ) THEN
               p_IEC%TurbInt15  = 0.18
               p_IEC%SigmaSlope = 2.0
            ELSEIF ( p_IEC%IECTurbC  == 'B' ) THEN
               p_IEC%TurbInt15  = 0.16
               p_IEC%SigmaSlope = 3.0
            ELSE   ! We should never get here, but just to be complete...
               CALL TS_Abort( ' Invalid IEC turbulence characteristic.' )
            ENDIF

            p_IEC%SigmaIEC(1) = p_IEC%TurbInt15*( ( 15.0 + p_IEC%SigmaSlope*UHub ) / ( p_IEC%SigmaSlope + 1.0 ) )
            p_IEC%TurbInt  = p_IEC%SigmaIEC(1)/UHub
         
         CASE ( 3 )

            IF ( p_IEC%IECTurbC == 'A' ) THEN
               p_IEC%TurbInt15  = 0.16
            ELSEIF ( p_IEC%IECTurbC == 'B' ) THEN
               p_IEC%TurbInt15  = 0.14
            ELSEIF ( p_IEC%IECTurbC == 'C' ) THEN
               p_IEC%TurbInt15  = 0.12
            ELSE   ! We should never get here, but just to be complete...
               CALL TS_Abort( ' Invalid IEC turbulence characteristic.' )
            ENDIF  

                   
            SELECT CASE ( p_IEC%IEC_WindType )
               CASE ( IEC_NTM )
                  p_IEC%SigmaIEC(1) = p_IEC%TurbInt15*( 0.75*UHub + 5.6 )                                      ! [IEC-1 Ed3 6.3.1.3 (11)]
               CASE ( IEC_ETM )
                  p_IEC%Vave        = 0.2*p_IEC%Vref                                                           ! [IEC-1 Ed3 6.3.1.1 ( 9)]
                  p_IEC%SigmaIEC(1) = p_IEC%ETMc * p_IEC%TurbInt15 * ( 0.072 * &
                                     ( p_IEC%Vave / p_IEC%ETMc + 3.0) * (Uhub / p_IEC%ETMc - 4.0)+10.0 )       ! [IEC-1 Ed3 6.3.2.3 (19)]
               CASE ( IEC_EWM1, IEC_EWM50 )
                  p_IEC%Vave        = 0.2*p_IEC%Vref                                                           ! [IEC-1 Ed3 6.3.1.1 ( 9)]
                  p_IEC%SigmaIEC(1) = 0.11*Uhub                                                                ! [IEC-1 Ed3 6.3.2.1 (16)]
               CASE DEFAULT 
                  CALL TS_Abort( 'Invalid IEC wind type.')
            END SELECT           
            p_IEC%TurbInt  = p_IEC%SigmaIEC(1)/UHub     
            
         CASE DEFAULT ! Likewise, this should never happen...

            CALL TS_Abort( 'Invalid IEC 61400-1 edition number.' )
            
         END SELECT                            
      

   ENDIF

   
      ! IEC turbulence scale parameter, Lambda(1), and IEC coherency scale parameter, LC

   IF ( p_IEC%IECedition == 2 ) THEN  
      
         ! section 6.3.1.3 Eq. 9
      IF ( p_grid%HubHt < 30.0 )  THEN
         p_IEC%Lambda(1) = 0.7*p_grid%HubHt
      ELSE
         p_IEC%Lambda(1) = 21.0
      ENDIF

      p_IEC%LC = 3.5*p_IEC%Lambda(1)

   ELSE !IF (p_IEC%IECedition == 3 )
      
         ! section 6.3.1.3 Eq. 9      
         
      IF ( p_grid%HubHt < 60.0 )  THEN
         p_IEC%Lambda(1) = 0.7*p_grid%HubHt
      ELSE
         p_IEC%Lambda(1) = 42.0
      ENDIF

      p_IEC%LC = 8.1*p_IEC%Lambda(1)

   ENDIF
   
      ! Set Lambda for Modified von Karman model:

   IF ( MVK .AND. SpecModel  == SpecModel_MODVKM ) THEN
      z0 = FindZ0(p_grid%HubHt, p_IEC%SigmaIEC(1), UHub, Fc)
      CALL ScaleMODVKM(p_grid%HubHt, UHub, p_IEC%Lambda(1), p_IEC%Lambda(2), p_IEC%Lambda(3))
   ENDIF   

   
      ! Sigma for v and w components and
      ! Integral scales (which depend on lambda)
   
   IF ( SpecModel == SpecModel_IECVKM ) THEN
      
      p_IEC%SigmaIEC(2)      =  1.0*p_IEC%SigmaIEC(1)
      p_IEC%SigmaIEC(3)      =  1.0*p_IEC%SigmaIEC(1)
      
      p_IEC%IntegralScale(:) =  3.5 *p_IEC%Lambda(1)   !L_k
            
   ELSE
      
      p_IEC%SigmaIEC(2)      =  0.8*p_IEC%SigmaIEC(1)
      p_IEC%SigmaIEC(3)      =  0.5*p_IEC%SigmaIEC(1)
            
      p_IEC%IntegralScale(1) =  8.1 *p_IEC%Lambda(1)   !L_k
      p_IEC%IntegralScale(2) =  2.7 *p_IEC%Lambda(1)   !L_k
      p_IEC%IntegralScale(3) =  0.66*p_IEC%Lambda(1)   !L_k
      
   END IF
   
   

END SUBROUTINE CalcIECScalingParams
!=======================================================================
SUBROUTINE TS_ValidateInput(ErrStat, ErrMsg)

use TSMods

   INTEGER(IntKi),                  intent(  out) :: ErrStat                         ! Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg                          ! Message describing error
   
   INTEGER(IntKi)                                 :: ErrStat2                        ! Error level (local)
   CHARACTER(MaxMsgLen)                           :: ErrMsg2                         ! Message describing error (local)
   



ErrStat = ErrID_None
ErrMsg  = ""


!BONNIE: UPPER LIMIT ON RICH_NO?
IF ( WrFile(FileExt_CTS) ) THEN
   
      ! models where coherent structures apply:
   IF ( SpecModel == SpecModel_GP_LLJ .OR. &
        SpecModel == SpecModel_NWTCUP .OR. &
        SpecModel == SpecModel_SMOOTH .OR. &
        SpecModel == SpecModel_WF_UPW .OR. &
        SpecModel == SpecModel_WF_07D .OR. &
        SpecModel == SpecModel_WF_14D ) THEN
      
      IF ( Rich_No <= -0.05 ) THEN
         CALL SetErrStat( ErrID_Info, 'A coherent turbulence time step file cannot be generated for RICH_NO <= -0.05.', ErrStat, ErrMsg, 'TS_ValidateInput')
         WrFile(FileExt_CTS) = .FALSE.      
      ELSEIF ( .NOT. ( WrFile(FileExt_BTS) .OR. WrFile(FileExt_WND) ) ) THEN
         CALL SetErrStat( ErrID_Info, 'AeroDyn Full-Field files(.bts) will be generated along with the coherent turbulence file.', ErrStat, ErrMsg, 'TS_ValidateInput')
         WrFile(FileExt_BTS) = .TRUE.
      ENDIF      
      
   ELSE
      CALL SetErrStat( ErrID_Info, 'A coherent turbulence time step file cannot be generated with the '//TRIM(TurbModel)//' model.', ErrStat, ErrMsg, 'TS_ValidateInput')
      WrFile(FileExt_CTS) = .FALSE.      
   END IF
      
ENDIF !WrAct


      ! Warn if EWM is used with incompatible times
      
IF ( ( p_IEC%IEC_WindType == IEC_EWM1 .OR. p_IEC%IEC_WindType == IEC_EWM50 ) .AND. & 
      ABS( 600.0_ReKi - MAX(p_grid%AnalysisTime,p_grid%UsableTime) ) > 90.0_ReKi ) THEN
   CALL SetErrStat( ErrID_Warn, 'The EWM parameters are valid for 10-min simulations only.', ErrStat, ErrMsg, 'TS_ValidateInput')   
ENDIF        


      ! Warn if Periodic is used with incompatible settings
      
IF ( p_grid%Periodic .AND. .NOT. EqualRealNos(p_grid%AnalysisTime, p_grid%UsableTime) ) THEN
   CALL SetErrStat( ErrID_Warn, 'Periodic output files will not be generated when AnalysisTime /= UsableTime. Setting Periodic = .FALSE.', ErrStat, ErrMsg, 'TS_ValidateInput')  
   p_grid%Periodic = .FALSE.
END IF


   ! Warn if tower points are output but grid is not:
IF ( WrFile(FileExt_TWR) .AND. .NOT. ( WrFile(FileExt_WND) .OR. WrFile(FileExt_BTS)) ) THEN 
   CALL SetErrStat( ErrID_Info, 'TurbSim .bts file will be generated to contain the tower points.', ErrStat, ErrMsg, 'TS_ValidateInput')  
   WrFile(FileExt_BTS) = .TRUE.
END IF



END SUBROUTINE TS_ValidateInput
!=======================================================================
SUBROUTINE TS_End()

   

END SUBROUTINE TS_END
!=======================================================================
END MODULE TSSubs
