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
MODULE TSSubs

   USE ModifiedvKrm_mod
   
   use TS_Profiles  
   use TS_RandNum
   use TS_VelocitySpectra
   USE NWTC_FFTPACK
   USE NWTC_LAPACK


   IMPLICIT NONE



CONTAINS

!=======================================================================
!> This subroutine returns the complex Fourier coefficients (packed in a
!! real array) of the simulated velocity (wind/water speed). It returns
!! values FOR ONLY the velocity components that use the IEC method for
!! computing spatial coherence; i.e., for i where SCMod(i) == CohMod_IEC
SUBROUTINE CalcFourierCoeffs_IEC( p, U, PhaseAngles, S, V, TRH, ErrStat, ErrMsg )

TYPE(TurbSim_ParameterType), INTENT(IN   )  :: p                            !< TurbSim parameters
REAL(ReKi),                  INTENT(IN)     :: U           (:)              !< The steady u-component wind speeds for the grid (NPoints).
REAL(ReKi),                  INTENT(IN)     :: PhaseAngles (:,:,:)          !< The array that holds the random phases [number of points, number of frequencies, number of wind components=3].
REAL(ReKi),                  INTENT(IN)     :: S           (:,:,:)          !< The turbulence PSD array (NumFreq,NPoints,3).
REAL(ReKi),                  INTENT(INOUT)  :: V           (:,:,:)          !< An array containing the summations of the rows of H (NumSteps,NPoints,3).
REAL(ReKi),                  INTENT(INOUT)  :: TRH (:)                      !< The transfer function matrix.  just used as a work array
INTEGER(IntKi),              INTENT(OUT)    :: ErrStat
CHARACTER(*),                INTENT(OUT)    :: ErrMsg

   
   ! Internal variables

REAL(ReKi), ALLOCATABLE       :: Dist(:)        ! The distance between points
REAL(ReKi), ALLOCATABLE       :: DistU(:)
   
INTEGER                       :: J
INTEGER                       :: I
INTEGER                       :: IFreq
INTEGER                       :: Indx
INTEGER                       :: IVec  ! wind component, 1=u, 2=v, 3=w
                              
INTEGER(IntKi)                :: ErrStat2
CHARACTER(MaxMsgLen)          :: ErrMsg2
   

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF (.NOT. ANY(p%met%SCMod == CohMod_IEC) ) RETURN

   !--------------------------------------------------------------------------------
   ! allocate arrays
   !--------------------------------------------------------------------------------
   CALL AllocAry( Dist,      p%grid%NPacked,      'Dist coherence array', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_IEC')
   CALL AllocAry( DistU,     p%grid%NPacked,     'DistU coherence array', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_IEC')
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF
   
   
   !--------------------------------------------------------------------------------
   ! Calculate the distances and other parameters that don't change with frequency
   !---------------------------------------------------------------------------------

   ! Calculate Dist array (distance between points I and J)
   ! and the DistU term, i.e., (r/u): u is uHub for IEC   
   Indx=0
   DO J=1,p%grid%NPoints
      DO I=J,p%grid%NPoints  ! The coherence matrix is symmetric so we're going to skip the other side 
         Indx = Indx + 1
         Dist(Indx)  = SQRT( ( p%grid%Y(I) - p%grid%Y(J) )**2  + ( p%grid%Z(I) - p%grid%Z(J) )**2 )
         DistU(Indx) = Dist(Indx)/p%UHub         
      END DO ! I
   END DO ! J
   
           
   !--------------------------------------------------------------------------------
   ! Calculate the fourier coefficients
   !---------------------------------------------------------------------------------
   
   DO IVec = 1,3
   
      IF (p%met%SCMod(IVec) /= CohMod_IEC) CYCLE ! Check the next component (this one doesn't use the IEC method)
   
      V(:,:,IVec) = 0.0_ReKi

      CALL WrScr ( '    '//Comp(IVec)//'-component matrices (IEC coherence method)' )

      !--------------------------------------------------------------------------------
      ! Calculate the coherence, Veers' H matrix (CSDs), and the fourier coefficients
      !---------------------------------------------------------------------------------

      DO IFREQ = 1,p%grid%NumFreq
         ! -----------------------------------------------
         ! Create the coherence matrix for this frequency
         ! -----------------------------------------------
         Indx = 1
         DO J = 1,p%usr%NPoints-1  ! start with the user-defined points (which don't get added coherence)
      
            TRH(Indx) =  1.0_ReKi
            Indx = Indx + 1
      
            DO I=J+1,p%grid%NPoints
               TRH(Indx) = 0.0_ReKi
               Indx = Indx + 1
            END DO !I
      
         END DO !J

         
         DO J=max(1, p%usr%NPoints),p%grid%NPoints
            DO I=J,p%grid%NPoints

                  TRH(Indx) = EXP( -1.0_ReKi * p%met%InCDec(IVec) * &
                              SQRT( (p%grid%Freq(IFreq)*DistU(Indx) )**2 + (p%met%InCohB(IVec)*Dist(Indx))**2 ) )
               
                  Indx = Indx  + 1

            ENDDO ! I
         ENDDO ! J
      
         ! -----------------------------------------------
         ! Now transform coherence to H matrix and then
         !   use H matrix to calculate coefficients
         ! -----------------------------------------------
         
         CALL Coh2H(    p, IVec, IFreq, TRH, S, ErrStat2, ErrMsg2 )       
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_IEC')
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF
         CALL H2Coeffs( IVec, IFreq, TRH, PhaseAngles, V, p%grid%NPoints )
      END DO !IFreq

   END DO !IVec   
               
   CALL Cleanup()
   RETURN
   
!............................................
CONTAINS
   SUBROUTINE Cleanup()

      IF ( ALLOCATED( Dist      ) ) DEALLOCATE( Dist      )
      IF ( ALLOCATED( DistU     ) ) DEALLOCATE( DistU     )
   END SUBROUTINE Cleanup
!............................................   
END SUBROUTINE CalcFourierCoeffs_IEC  
!=======================================================================
!> This subroutine returns the complex Fourier coefficients (packed in a
!! real array) of the simulated velocity (wind/water speed). It returns
!! values FOR ONLY the velocity components that use the general method for
!! computing spatial coherence; i.e., for i where SCMod(i) == CohMod_GENERAL
SUBROUTINE CalcFourierCoeffs_General( p, U, PhaseAngles, S, V, TRH, ErrStat, ErrMsg )

TYPE(TurbSim_ParameterType), INTENT(IN   )  :: p                            !< TurbSim parameters
REAL(ReKi),                  INTENT(IN)     :: U           (:)              !< The steady u-component wind speeds for the grid (NPoints).
REAL(ReKi),                  INTENT(IN)     :: PhaseAngles (:,:,:)          !< The array that holds the random phases [number of points, number of frequencies, number of wind components=3].
REAL(ReKi),                  INTENT(IN)     :: S           (:,:,:)          !< The turbulence PSD array (NumFreq,NPoints,3).
REAL(ReKi),                  INTENT(INOUT)  :: V           (:,:,:)          !< An array containing the summations of the rows of H (NumSteps,NPoints,3).
REAL(ReKi),                  INTENT(INOUT)  :: TRH (:)                      !< The transfer function matrix.  just used as a work array
INTEGER(IntKi),              INTENT(OUT)    :: ErrStat
CHARACTER(*),                INTENT(OUT)    :: ErrMsg

   
   ! Internal variables

INTEGER                      :: UC              ! I/O unit for Coherence debugging file.
LOGICAL,    PARAMETER        :: COH_OUT = .FALSE.                       ! This parameter has been added to replace the NON-STANDARD compiler directive previously used

REAL(ReKi), ALLOCATABLE       :: Dist(:)        ! The distance between points
REAL(ReKi), ALLOCATABLE       :: DistU(:)
REAL(ReKi), ALLOCATABLE       :: DistZMExp(:)
   
REAL(ReKi)                    :: dY             ! the lateral distance between two points
REAL(ReKi)                    :: UM             ! The mean wind speed of the two points
REAL(ReKi)                    :: ZM             ! The mean height of the two points

INTEGER                       :: J
INTEGER                       :: I
INTEGER                       :: IFreq
INTEGER                       :: Indx
INTEGER                       :: IVec  ! wind component, 1=u, 2=v, 3=w
                              
INTEGER(IntKi)                :: ErrStat2
CHARACTER(MaxMsgLen)          :: ErrMsg2
   

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   IF (.NOT. ANY(p%met%SCMod == CohMod_GENERAL) ) RETURN

   
   !--------------------------------------------------------------------------------
   ! allocate arrays
   !--------------------------------------------------------------------------------
   CALL AllocAry( Dist,      p%grid%NPacked,      'Dist coherence array', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_General')
   CALL AllocAry( DistU,     p%grid%NPacked,     'DistU coherence array', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_General')
   CALL AllocAry( DistZMExp, p%grid%NPacked, 'DistZMExp coherence array', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_General')
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF
   
      
   !--------------------------------------------------------------------------------
   ! Calculate the distances and other parameters that don't change with frequency
   !---------------------------------------------------------------------------------

      ! Calculate Dist array (distance between points I and J)
   IF ( .NOT. PeriodicY ) THEN
      Indx=0
      DO J=1,p%grid%NPoints
         DO I=J,p%grid%NPoints  ! The coherence matrix is symmetric so we're going to skip the other side 
            Indx = Indx + 1
            Dist(Indx)= SQRT( ( p%grid%Y(I) - p%grid%Y(J) )**2  + ( p%grid%Z(I) - p%grid%Z(J) )**2 )
         END DO ! I
      END DO ! J
   ELSE  
      ! bjj need to test ths more!!!
      Indx=0
      DO J=1,p%grid%NPoints
         DO I=J,p%grid%NPoints  ! The coherence matrix is symmetric so we're going to skip the other side 
      
            Indx = Indx + 1
            dY = p%grid%Y(I) - p%grid%Y(J)
            IF (dY > 0.5*p%grid%GridWidth ) THEN
               dY = dY - p%grid%GridWidth - p%grid%GridRes_Y
            ELSE IF (dY < -0.5*p%grid%GridWidth ) THEN
               dY = dY + p%grid%GridWidth + p%grid%GridRes_Y
            END IF

            Dist(Indx)= SQRT( ( dY )**2  + ( p%grid%Z(I) - p%grid%Z(J) )**2 )

         END DO
      END DO
   END IF


      ! Compute the DistZMExp term, i.e., -(r/z_m)^CohExp
   IF ( EqualRealNos( p%met%COHEXP, 0.0_ReKi ) ) THEN
      DistZMExp = -1.0_ReKi      ! value for entire array
   ELSE   
      Indx=0
      DO J=1,p%grid%NPoints
         DO I=J,p%grid%NPoints  ! The coherence matrix is symmetric so we're going to skip the other side 
      
            Indx            = Indx + 1
            ZM              = 0.5*( p%grid%Z(I) + p%grid%Z(J) )
            DistZMExp(Indx) = -1.0_ReKi*( Dist(Indx)/ZM )**p%met%COHEXP     ! Note: 0**0 = 1      
         END DO ! I  
      END DO ! J              
   END IF

      ! Compute the DistU term, i.e., (r/u): u is average u at points I and J
   Indx=0
   DO J=1,p%grid%NPoints
      DO I=J,p%grid%NPoints  ! The coherence matrix is symmetric so we're going to skip the other side      
         Indx        = Indx + 1
         UM          = p%UHub ! was: 0.5*( U(I) + U(J) )
         DistU(Indx) = Dist(Indx)/UM
      END DO ! I  
   END DO ! J 
   
   !.................
   ! DEBUGGING
   !.................
IF ( COH_OUT ) THEN !debugging info...

      ! Write the coherence for three frequencies, for debugging purposes
      CALL GetNewUnit( UC, ErrStat2, ErrMsg2 );  CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_General')
      
      CALL OpenFOutFile( UC, TRIM(p%RootName)//'.coh', ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_General')
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
      
      WRITE( UC, '(A4,X,A16,1X,'//Num2LSTR(p%grid%NPacked)//'(G10.4,1X))' ) 'Comp','Freq',(I,I=1,p%grid%NPacked)
      WRITE( UC,   '(5X,A16,1X,'//Num2LSTR(p%grid%NPacked)//'(G10.4,1X))' ) 'Distance',   Dist(:)
      WRITE( UC,   '(5X,A16,1X,'//Num2LSTR(p%grid%NPacked)//'(G10.4,1X))' ) '(r/u)', DistU(:)
      WRITE( UC,   '(5X,A16,1X,'//Num2LSTR(p%grid%NPacked)//'(G10.4,1X))' ) '(u)',  p%met%URef
      WRITE( UC,   '(5X,A16,1X,'//Num2LSTR(p%grid%NPacked)//'(G10.4,1X))' ) '-(r/z_m)^CohExp', DistZMExp(:)
ENDIF

   
   !--------------------------------------------------------------------------------
   ! Calculate the fourier coefficients
   !---------------------------------------------------------------------------------
   
   DO IVec = 1,3
   
      IF (p%met%SCMod(IVec) /= CohMod_GENERAL) CYCLE ! Check the next component (this one doesn't use the GENERAL method)
   
      V(:,:,IVec) = 0.0_ReKi

      CALL WrScr ( '    '//Comp(IVec)//'-component matrices (general coherence model)' )

      !--------------------------------------------------------------------------------
      ! Calculate the coherence, Veers' H matrix (CSDs), and the fourier coefficients
      !---------------------------------------------------------------------------------

      DO IFREQ = 1,p%grid%NumFreq
         ! -----------------------------------------------
         ! Create the coherence matrix for this frequency
         ! -----------------------------------------------
   
         Indx = 1
         DO J = 1,p%usr%NPoints-1 ! start with user-defined points (which don't get added coherence)
      
            TRH(Indx) =  1.0_ReKi
            Indx = Indx + 1
      
            DO I=J+1,p%grid%NPoints
               TRH(Indx) = 0.0_ReKi
               Indx = Indx + 1
            END DO !I
      
         END DO !J
         
         DO J=max(1, p%usr%NPoints),p%grid%NPoints
            DO I=J,p%grid%NPoints

                  TRH(Indx) = EXP( p%met%InCDec(IVec) * DistZMExp(Indx)* &
                              SQRT( (p%grid%Freq(IFreq)*DistU(Indx) )**2 + (p%met%InCohB(IVec)*Dist(Indx))**2 ) )
               
                  Indx = Indx  + 1

            ENDDO ! I
         ENDDO ! J
   !.................
   ! DEBUGGING
   !.................
         IF (COH_OUT) THEN
   !        IF (IFreq == 1 .OR. IFreq == p%grid%NumFreq) THEN
               WRITE( UC, '(I3,2X,F15.5,1X,'//Num2LSTR(p%grid%NPacked)//'(G10.4,1X))' ) IVec, p%grid%Freq(IFreq), TRH(1:p%grid%NPacked)
   !        ENDIF
         ENDIF 

         ! -----------------------------------------------
         ! Now transform coherence to H matrix and then
         !   use H matrix to calculate coefficients
         ! -----------------------------------------------
         
         CALL Coh2H(    p, IVec, IFreq, TRH, S, ErrStat2, ErrMsg2 )       
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_General')
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF
         CALL H2Coeffs( IVec, IFreq, TRH, PhaseAngles, V, p%grid%NPoints )
      END DO !IFreq

   END DO !IVec   
   
   CALL Cleanup()
   RETURN
!............................................
CONTAINS
   SUBROUTINE Cleanup()

      IF (COH_OUT .AND. UC > 0)  CLOSE( UC )

      IF ( ALLOCATED( Dist      ) ) DEALLOCATE( Dist      )
      IF ( ALLOCATED( DistU     ) ) DEALLOCATE( DistU     )
      IF ( ALLOCATED( DistZMExp ) ) DEALLOCATE( DistZMExp )
   END SUBROUTINE Cleanup
!............................................   
END SUBROUTINE CalcFourierCoeffs_General
!=======================================================================
!> This subroutine returns the complex Fourier coefficients (packed in a
!! real array) of the simulated velocity (wind/water speed).
!! It returns the values FOR ONLY the velocity components that use identity
!! spatial coherence; i.e., for i where SCMod(i) == CohMod_NONE
SUBROUTINE CalcFourierCoeffs_NONE( p, U, PhaseAngles, S, V, TRH, ErrStat, ErrMsg )

TYPE(TurbSim_ParameterType), INTENT(IN   )  :: p                            !< TurbSim parameters
REAL(ReKi),                  INTENT(IN)     :: U           (:)              !< The steady u-component wind speeds for the grid (NPoints).
REAL(ReKi),                  INTENT(IN)     :: PhaseAngles (:,:,:)          !< The array that holds the random phases [number of points, number of frequencies, number of wind components=3].
REAL(ReKi),                  INTENT(IN)     :: S           (:,:,:)          !< The turbulence PSD array (NumFreq,NPoints,3).
REAL(ReKi),                  INTENT(INOUT)  :: V           (:,:,:)          !< An array containing the summations of the rows of H (NumSteps,NPoints,3).
REAL(ReKi),                  INTENT(INOUT)  :: TRH (:)                      !< The transfer function matrix.  just used as a work array
INTEGER(IntKi),              INTENT(OUT)    :: ErrStat
CHARACTER(*),                INTENT(OUT)    :: ErrMsg

   
   ! Internal variables   
INTEGER                       :: IFreq
INTEGER                       :: IVec  ! wind component, 1=u, 2=v, 3=w
                                
   
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   !--------------------------------------------------------------------------------
   ! Calculate the fourier coefficients
   !---------------------------------------------------------------------------------
   
   DO IVec = 1,3
   
      IF (p%met%SCMod(IVec) /= CohMod_NONE) CYCLE ! Check the next component (this one doesn't use the identity coherence method)
   
      V(:,:,IVec) = 0.0_ReKi

      CALL WrScr ( '    '//Comp(IVec)//'-component matrices (identity coherence)' )
   
   
   ! now calculate coherence for compents that use this method
              
      ! -----------------------------------------------------------------------------------
      !  The coherence is the Identity (as is Cholesky Factorization); 
      !    the Veers' H matrix calculated in EyeCoh2H:
      ! -----------------------------------------------------------------------------------
      
      DO IFREQ = 1,p%grid%NumFreq
         CALL EyeCoh2H(  IVec, IFreq, TRH, S,              p%grid%NPoints )   
         CALL H2Coeffs(  IVec, IFreq, TRH, PhaseAngles, V, p%grid%NPoints )
      ENDDO !IFreq
      
   END DO ! IVec

   RETURN      
!............................................      
END SUBROUTINE CalcFourierCoeffs_NONE
!=======================================================================
!> This subroutine computes the coherence between two points on the grid,
!! forms the cross spectrum matrix, and returns the complex
!! Fourier coefficients of the simulated velocity (wind speed).
SUBROUTINE CalcFourierCoeffs( p, U, PhaseAngles, S, V, ErrStat, ErrMsg )

IMPLICIT                      NONE

   ! Passed variables

TYPE(TurbSim_ParameterType), INTENT(IN   )  :: p                            !< TurbSim parameters
REAL(ReKi),                  INTENT(in)     :: U           (:)              !< The steady u-component wind speeds for the grid (NPoints).
REAL(ReKi),                  INTENT(IN)     :: PhaseAngles (:,:,:)          !< The array that holds the random phases [number of points, number of frequencies, number of wind components=3].
REAL(ReKi),                  INTENT(IN)     :: S           (:,:,:)          !< The turbulence PSD array (NumFreq,NPoints,3).
REAL(ReKi),                  INTENT(  OUT)  :: V           (:,:,:)          !< An array containing the summations of the rows of H (NumSteps,NPoints,3).
INTEGER(IntKi),              INTENT(OUT)    :: ErrStat
CHARACTER(*),                INTENT(OUT)    :: ErrMsg

   ! Internal variables

REAL(ReKi),                   ALLOCATABLE   :: TRH (:)        ! The transfer function matrix.  
INTEGER(IntKi)                              :: ErrStat2
CHARACTER(MaxMsgLen)                        :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! with no turbulence, we return all zeros
   IF (p%met%TurbModel_ID == SpecModel_NONE) THEN
      V = 0.0_ReKi
      RETURN
   END IF
         
      ! otherwise, we use the coherence method specified by the user
      
   CALL AllocAry( TRH, p%grid%NPacked, 'TRH coherence array', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs')
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF

   CALL CalcFourierCoeffs_IEC(     p, U, PhaseAngles, S, V, TRH, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs')   
   CALL CalcFourierCoeffs_API(     p, U, PhaseAngles, S, V, TRH, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs')   
   CALL CalcFourierCoeffs_General( p, U, PhaseAngles, S, V, TRH, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs')   
   CALL CalcFourierCoeffs_NONE(    p, U, PhaseAngles, S, V, TRH, ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs')   
      
   CALL Cleanup()
   RETURN
!............................................
CONTAINS
   SUBROUTINE Cleanup()
      IF ( ALLOCATED( TRH ) ) DEALLOCATE( TRH )
   END SUBROUTINE Cleanup
!............................................
END SUBROUTINE CalcFourierCoeffs

!=======================================================================
!> This subroutine returns the complex Fourier coefficients (packed in a
!! real array) of the simulated velocity (wind/water speed). It returns
!! values FOR ONLY the velocity components that use the API method for
!! computing spatial coherence; i.e., for i where SCMod(i) == CohMod_API
SUBROUTINE CalcFourierCoeffs_API( p, U, PhaseAngles, S, V, TRH, ErrStat, ErrMsg )

IMPLICIT                      NONE

   ! Passed variables
TYPE(TurbSim_ParameterType), INTENT(IN   )  :: p                            !< TurbSim parameters
REAL(ReKi),                  INTENT(IN   )  :: U           (:)              !< The steady u-component wind speeds for the grid (NPoints).
REAL(ReKi),                  INTENT(IN   )  :: PhaseAngles (:,:,:)          !< The array that holds the phase angles [number of points, number of frequencies, number of wind components=3].
REAL(ReKi),                  INTENT(IN   )  :: S           (:,:,:)          !< The turbulence PSD array (NumFreq,NPoints,3).
REAL(ReKi),                  INTENT(INOUT)  :: V           (:,:,:)          !< An array containing the summations of the rows of H (NumSteps,NPoints,3).
REAL(ReKi),                  INTENT(INOUT)  :: TRH (:)                      !< The transfer function matrix.  just used as a work array
INTEGER(IntKi),              INTENT(  OUT)  :: ErrStat
CHARACTER(*),                INTENT(  OUT)  :: ErrMsg

   ! Internal variables

REAL(ReKi), ALLOCATABLE      :: Dist_Y(:)        ! The Y distance between points
REAL(ReKi), ALLOCATABLE      :: Dist_Z(:)        ! The Z distance between points
REAL(ReKi), ALLOCATABLE      :: z_g(:)           ! sqrt( Z(IZ)*Z(JZ) ) / H

INTEGER                      :: J
INTEGER                      :: I
INTEGER                      :: K
INTEGER                      :: IFreq
INTEGER                      :: Indx
INTEGER                      :: IVec            ! wind component, 1=u, 2=v, 3=w

INTEGER                      :: UC             ! I/O unit for Coherence debugging file.
LOGICAL,    PARAMETER        :: COH_OUT   = .FALSE.                       ! This parameter has been added to replace the NON-STANDARD compiler directive previously used

INTEGER(IntKi)                :: ErrStat2
CHARACTER(MaxMsgLen)          :: ErrMsg2

REAL, PARAMETER :: Qc(3)     = (/ 1.00, 1.00, 1.25 /)
REAL, PARAMETER :: Pc(3)     = (/ 0.40, 0.40, 0.50 /)
REAL, PARAMETER :: Rc(3)     = (/ 0.92, 0.92, 0.85 /)
REAL, PARAMETER :: Alpha(3) = (/ 2.9 ,45.0 ,13.0  /)

REAL, PARAMETER :: H             = 10. ! Reference height
REAL, PARAMETER :: Coef_1      =   1.0 !3.28

REAL :: A_Y, A_Z


   ! initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ""
   UC      = -1

   IF (.NOT. ANY(p%met%SCMod == CohMod_API) ) RETURN

   !--------------------------------------------------------------------------------
   ! allocate arrays
   !--------------------------------------------------------------------------------
   
   CALL AllocAry( Dist_Y,    p%grid%NPacked, 'Dist_Y coherence array', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_API')
   CALL AllocAry( Dist_Z,    p%grid%NPacked, 'Dist_Z coherence array', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_API')
   !CALL AllocAry( Dist_Z12, p%grid%NPacked, 'Dist_Z12 coherence array', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_API')
   CALL AllocAry( z_g,       p%grid%NPacked, 'z_g coherence array', ErrStat2, ErrMsg2 ); CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_API')
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF


   !--------------------------------------------------------------------------------
   ! Calculate the distances and other parameters that don't change with frequency
   !---------------------------------------------------------------------------------

   Indx=0
   DO J=1,p%grid%NPoints
      DO I=J,p%grid%NPoints  ! The coherence matrix is symmetric so we're going to skip the other side 
         Indx = Indx + 1      ! Indx = p%grid%NPoints*(J - 1) - J*(J - 1)/2 + I   !Index of packed V matrix, coherence between points I & J

         Dist_Y(Indx)= ABS( p%grid%Y(I) - p%grid%Y(J) )
         Dist_Z(Indx)= ABS( p%grid%Z(I) - p%grid%Z(J) )
         z_g(Indx)  = sqrt( p%grid%Z(I) * p%grid%Z(J) ) / H
                  
      END DO
   END DO
   

   !.................
   ! DEBUGGING
   !.................
IF ( COH_OUT ) THEN !debugging info...

      ! Write the coherence for three frequencies, for debugging purposes
      CALL GetNewUnit( UC, ErrStat2, ErrMsg2 );  CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_API')
      
      CALL OpenFOutFile( UC, TRIM(p%RootName)//'.coh', ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_API')
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         END IF
      
      WRITE( UC, '(A4,X,A16,1X,'//Num2LSTR(p%grid%NPacked)//'(G10.4,1X))' ) 'Comp','Freq',(I,I=1,p%grid%NPacked)
      WRITE( UC,   '(5X,A16,1X,'//Num2LSTR(p%grid%NPacked)//'(G10.4,1X))' ) 'Distance_Y',   Dist_Y(:)
      WRITE( UC,   '(5X,A16,1X,'//Num2LSTR(p%grid%NPacked)//'(G10.4,1X))' ) 'Distance_Z', Dist_Z(:)
      WRITE( UC,   '(5X,A16,1X,'//Num2LSTR(p%grid%NPacked)//'(G10.4,1X))' ) 'sqrt(Z(IZ)*Z(JZ))/H', z_g(:)
ENDIF

   !--------------------------------------------------------------------------------
   ! Calculate the fourier coefficients
   !---------------------------------------------------------------------------------
   
   DO IVec = 1,1 !BJJ: note that only the u component is defined, and I don't want to look at how to change the coherence in the other components....
   
      IF (p%met%SCMod(IVec) /= CohMod_API) CYCLE ! Check the next component (this one doesn't use the API method)

      CALL WrScr ( '    '//Comp(IVec)//'-component matrices (2-dimensional API coherence method)' )
      
   !--------------------------------------------------------------------------------
   ! Calculate the coherence, Veers' H matrix (CSDs), and the fourier coefficients
   !---------------------------------------------------------------------------------
      
      DO IFREQ = 1,p%grid%NumFreq
         ! -----------------------------------------------
         ! Create the coherence matrix for this frequency
         ! -----------------------------------------------
   
         Indx = 1
         DO J = 1,p%usr%NPoints-1 ! start with user-defined points (which don't get added coherence)
      
            TRH(Indx) =  1.0_ReKi
            Indx = Indx + 1
      
            DO I=J+1,p%grid%NPoints
               TRH(Indx) = 0.0_ReKi
               Indx = Indx + 1
            END DO !I
      
         END DO !J
         
         DO J=max(1, p%usr%NPoints),p%grid%NPoints
            DO I=J,p%grid%NPoints   
   
               A_Y = Alpha(2) * (p%grid%Freq(IFreq)**rc(2)) * ((Dist_Y(Indx)/Coef_1)**qc(2)) * (z_g(Indx)**(-pc(2)))
               A_Z = Alpha(3) * (p%grid%Freq(IFreq)**rc(3)) * ((Dist_Z(Indx)/Coef_1)**qc(3)) * (z_g(Indx)**(-pc(3)))

               TRH(Indx)=EXP(- Coef_1 * SQRT(A_Y**2 + A_Z**2) / p%met%URef)
               
               Indx = Indx  + 1                                    
               
            ENDDO ! I
         ENDDO ! J
      
      
   !.................
   ! DEBUGGING
   !.................
         IF (COH_OUT) THEN
   !        IF (IFreq == 1 .OR. IFreq == p%grid%NumFreq) THEN
               WRITE( UC, '(I3,2X,F15.5,1X,'//Num2LSTR(p%grid%NPacked)//'(G10.4,1X))' ) IVec, p%grid%Freq(IFreq), TRH(1:p%grid%NPacked)
   !        ENDIF
         ENDIF
      
         ! -----------------------------------------------
         ! Now transform coherence to H matrix and then
         !   use H matrix to calculate coefficients
         ! -----------------------------------------------
      
         CALL Coh2H(    p, IVec, IFreq, TRH, S, ErrStat2, ErrMsg2 )
            CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CalcFourierCoeffs_API')
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            END IF
      
         CALL H2Coeffs( IVec, IFreq, TRH, PhaseAngles, V, p%grid%NPoints )
      
      ENDDO !IFreq
   ENDDO !IVec

   CALL Cleanup()

RETURN
!............................................
CONTAINS
   SUBROUTINE Cleanup()

      IF (COH_OUT .AND. UC > 0)  CLOSE( UC )

      IF ( ALLOCATED( Dist_Y ) ) DEALLOCATE( Dist_Y )
      IF ( ALLOCATED( Dist_Z ) ) DEALLOCATE( Dist_Z )
      IF ( ALLOCATED( z_g    ) ) DEALLOCATE( z_g    )
   END SUBROUTINE Cleanup
!............................................
END SUBROUTINE CalcFourierCoeffs_API
!=======================================================================
SUBROUTINE EyeCoh2H( IVec, IFreq, TRH, S, NPoints )

REAL(ReKi),       INTENT(INOUT)  :: TRH         (:)                          ! The transfer function  matrix (length is >= p%grid%NPacked).
REAL(ReKi),       INTENT(IN)     :: S           (:,:,:)                      ! The turbulence PSD array (NumFreq,NPoints,3).
INTEGER(IntKi),   INTENT(IN)     :: IVec                                     ! loop counter (=number of wind components)
INTEGER(IntKi),   INTENT(IN)     :: IFreq                                    ! loop counter (=number of frequencies)
INTEGER(IntKi),   INTENT(IN)     :: NPoints                                  ! Size of dimension 2 of S

integer                          :: Indx, J, I

!NPoints = SIZE(S,2)

      ! -----------------------------------------------------------------------------------
      !  The coherence is the Identity (as is Cholesky); the Veers' H matrix is as follows:
      ! -----------------------------------------------------------------------------------

   Indx = 1
   DO J = 1,NPoints ! The column number

         ! The diagonal entries of the matrix:

      TRH(Indx) = SQRT( ABS( S(IFreq,J,IVec) ) )

         ! The off-diagonal values:
      Indx = Indx + 1
      DO I = J+1,NPoints ! The row number
         TRH(Indx) = 0.0
         Indx = Indx + 1
      ENDDO ! I
   ENDDO ! J

END SUBROUTINE EyeCoh2H
!=======================================================================
SUBROUTINE Coh2H( p, IVec, IFreq, TRH, S, ErrStat, ErrMsg )

!use NWTC_LAPACK

TYPE(TurbSim_ParameterType), INTENT(IN   )  :: p                                        ! TurbSim parameters
REAL(ReKi),                  INTENT(INOUT)  :: TRH         (:)                          ! The transfer function  matrix (size >= NumSteps).
REAL(ReKi),                  INTENT(IN)     :: S           (:,:,:)                      ! The turbulence PSD array (NumFreq,NPoints,3).
INTEGER(IntKi),              INTENT(IN)     :: IVec                                     ! loop counter (=number of wind components)
INTEGER(IntKi),              INTENT(IN)     :: IFreq                                    ! loop counter (=number of frequencies)

INTEGER(IntKi),              INTENT(OUT)    :: ErrStat
CHARACTER(*),                INTENT(OUT)    :: ErrMsg


integer                          :: Indx, J, I, NPts

         
      ! -------------------------------------------------------------
      ! Calculate the Cholesky factorization for the coherence matrix
      ! -------------------------------------------------------------      
   IF ( p%usr%NPoints > 0 ) THEN
      J = p%usr%NPoints
      Indx = p%grid%NPoints*(J-1) - J*(J-1)/2 + J !Index of H(J,J)
      NPts = p%grid%NPoints - p%usr%NPoints + 1
   ELSE
      Indx = 1
      NPts = p%grid%NPoints
   END IF
         
   CALL LAPACK_pptrf( 'L', NPts, TRH(Indx:), ErrStat, ErrMsg )  ! 'L'ower triangular 'TRH' matrix (packed form), of order 'NPoints'; returns Stat

   IF ( ErrStat /= ErrID_None ) THEN
      IF (ErrStat < AbortErrLev) then
         CALL WrScr(ErrMsg)
      ELSE
         ErrMsg = 'Error in Cholesky factorization: '//TRIM(ErrMsg)//newline//&
                  'The error occurred in the '//Comp(IVec)//'-component coherence matrix at frequency '//&
                   TRIM(Int2LStr(IFreq))//' ('//TRIM(Num2LStr(p%grid%Freq(IFreq)))//' Hz)'//&
                  ' Check the input file for invalid physical properties or modify the coherence.'
         RETURN
                     
      END IF      
   ENDIF

      ! -------------------------------------------------------------
      ! Create the lower triangular matrix, H, from Veer's method
      ! -------------------------------------------------------------
   
   Indx = 1
   DO J = 1,p%usr%NPoints-1  ! Column               
      ! use identity coherence for the user-input time series:
      
      !Indx = p%grid%NPoints*(J-1) - J*(J-1)/2 + J !Index of H(J,J)
      
      TRH(Indx) =  SQRT( ABS( S(IFreq,J,IVec) ) )
      Indx = Indx + 1
      
      DO I=J+1,p%grid%NPoints
         TRH(Indx) = 0.0_ReKi
         Indx = Indx + 1
      END DO
      
   END DO !J
         
   DO J = max(1,p%usr%NPoints),p%grid%NPoints  ! Column
      
      TRH(Indx) =  TRH(Indx) * SQRT( ABS( S(IFreq,J,IVec) ) )
      Indx = Indx + 1      
      
      DO I = J+1,p%grid%NPoints ! Row

            ! S(IFreq,I,IVec) should never be less than zero, but the ABS makes sure...
            !Indx = NPoints*(J-1) - J*(J-1)/2 + I !Index of H(I,J)
            
!         TRH(Indx) = TRH(Indx) * SQRT( ABS( S(IFreq,I,IVec) ) )
         TRH(Indx) = TRH(Indx) * SQRT( SQRT( ABS( S(IFreq,I,IVec) * S(IFreq,J,IVec) ) ) )

         Indx = Indx + 1

      ENDDO !I
   ENDDO !J
         

END SUBROUTINE Coh2H
!=======================================================================
SUBROUTINE H2Coeffs( IVec, IFreq, TRH, PhaseAngles, V, NPoints )


REAL(ReKi),       INTENT(IN)     :: TRH         (:)                          ! The transfer function  matrix (length is >= p%grid%NPacked).
REAL(ReKi),       INTENT(IN)     :: PhaseAngles (:,:,:)                      ! The array that holds the random phases [number of points, number of frequencies, number of wind components=3].
REAL(ReKi),       INTENT(INOUT)  :: V           (:,:,:)                      ! An array containing the summations of the rows of H (NumSteps,NPoints,3).
INTEGER(IntKi),   INTENT(IN)     :: IVec                                     ! loop counter (=number of wind components)
INTEGER(IntKi),   INTENT(IN)     :: IFreq                                    ! loop counter (=number of frequencies)
INTEGER(IntKi),   INTENT(IN)     :: NPoints                                  ! Size of dimension 2 of V


REAL(ReKi)                       :: CPh                                      ! Cosine of the random phase
REAL(ReKi)                       :: SPh                                      ! Sine of the random phase
INTEGER                          :: IF1                                      ! Index to real part of vector
INTEGER                          :: IF2                                      ! Index to complex part of vector

integer                          :: Indx, J, I


   ! -------------------------------------------------------------
   ! Calculate the correlated fourier coefficients.
   ! -------------------------------------------------------------

   IF2      = IFreq*2
   IF1      = IF2 - 1

   DO J=1,NPoints

         ! Apply a random phase to each of the columns of H to
         ! produce random phases in the wind component.
         ! Then sum each of the rows into the vector V.
         
      CPh   = COS( PhaseAngles(J,IFreq,IVec) )
      SPh   = SIN( PhaseAngles(J,IFreq,IVec) )

      Indx = NPoints*(J-1) - J*(J-1)/2 + J !Index of H(J,J)
      DO I=J,NPoints

         V(IF1,I,IVec) = V(IF1,I,IVec) + TRH(Indx)*CPh  !Real part
         V(IF2,I,IVec) = V(IF2,I,IVec) + TRH(Indx)*SPh  !Imaginary part

         Indx = Indx + 1      !H(I,J)

      ENDDO ! I
   ENDDO ! J

END SUBROUTINE H2Coeffs
!=======================================================================
!> This routine takes the Fourier coefficients and converts them to velocity
!! note that the resulting time series has zero mean.
SUBROUTINE Coeffs2TimeSeries( V, NumSteps, NPoints, NUsrPoints, ErrStat, ErrMsg )


   !USE NWTC_FFTPACK

   IMPLICIT NONE 
   

   ! passed variables
   INTEGER(IntKi),   INTENT(IN)     :: NumSteps                     !< Size of dimension 1 of V (number of time steps)
   INTEGER(IntKi),   INTENT(IN)     :: NPoints                      !< Size of dimension 2 of V (number of grid points)
   INTEGER(IntKi),   INTENT(IN)     :: NUsrPoints                   !< number of user-defined time series

   REAL(ReKi),       INTENT(INOUT)  :: V     (NumSteps,NPoints,3)   !< An array containing the summations of the rows of H (NumSteps,NPoints,3).

   INTEGER(IntKi),   intent(  out)  :: ErrStat                      !< Error level
   CHARACTER(*),     intent(  out)  :: ErrMsg                       !< Message describing error
   
   
   ! local variables
   TYPE(FFT_DataType)               :: FFT_Data                      ! data for applying FFT
   REAL(SiKi), ALLOCATABLE          :: Work ( : )                    ! working array to hold coefficients of fft  !bjj: made it allocatable so it doesn't take stack space

   
   INTEGER(IntKi)                   :: ITime                         ! loop counter for time step/frequency 
   INTEGER(IntKi)                   :: IVec                          ! loop counter for velocity components
   INTEGER(IntKi)                   :: IPoint                        ! loop counter for grid points
   
   INTEGER(IntKi)                   :: ErrStat2                      ! Error level (local)
  !CHARACTER(MaxMsgLen)             :: ErrMsg2                       ! Message describing error (local)
   

   ! initialize variables

 !ErrStat = ErrID_None
 !ErrMsg = ""
   
   CALL AllocAry(Work, NumSteps, 'Work',ErrStat,ErrMsg)
   if (ErrStat >= AbortErrLev) return
   
   !  Allocate the FFT working storage and initialize its variables

CALL InitFFT( NumSteps, FFT_Data, ErrStat=ErrStat2 )
   CALL SetErrStat(ErrStat2, 'Error in InitFFT', ErrStat, ErrMsg, 'Coeffs2TimeSeries' )
   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF


   ! Get the stationary-point time series.

CALL WrScr ( ' Generating time series for all points:' )

DO IVec=1,3

   CALL WrScr ( '    '//Comp(IVec)//'-component' )

   DO IPoint=1,NPoints    !NTotB

         ! Overwrite the first point with zero.  This sets the real (and 
         ! imaginary) part of the steady-state value to zero so that we 
         ! can add in the mean value later.

      Work(1) = 0.0_ReKi

!      DO ITime = 2,NumSteps-1
      DO ITime = 2,NumSteps
         Work(ITime) = V(ITime-1, IPoint, IVec)
      ENDDO ! ITime

      IF (iPoint > NUsrPoints) THEN
         ! BJJ: we can't override this for the user-input spectra or we don't get the correct time series out.
         ! Per JMJ, I will keep this here for the other points, but I personally think it could be skipped, too.
         
         ! Now, let's add a complex zero to the end to set the power in the Nyquist
         ! frequency to zero.
         
         Work(NumSteps) = 0.0
      END IF
      
      

         ! perform FFT

      CALL ApplyFFT( Work, FFT_Data, ErrStat2 )
         IF (ErrStat2 /= ErrID_None ) THEN
            CALL SetErrStat(ErrStat2, 'Error in ApplyFFT for point '//TRIM(Num2LStr(IPoint))//'.', ErrStat, ErrMsg, 'Coeffs2TimeSeries' )
            IF (ErrStat >= AbortErrLev) EXIT
         END IF
        
      V(:,IPoint,IVec) = Work

   ENDDO ! IPoint

ENDDO ! IVec 

CALL Cleanup()

RETURN
CONTAINS
!...........................................
SUBROUTINE Cleanup()
   
   CALL ExitFFT( FFT_Data, ErrStat2 )
   CALL SetErrStat(ErrStat2, 'Error in ExitFFT', ErrStat, ErrMsg, 'Coeffs2TimeSeries' )

   if (allocated(work)) deallocate(work)
   
   END SUBROUTINE Cleanup
END SUBROUTINE Coeffs2TimeSeries
!=======================================================================
!> This routine calculates the two-sided Fourier amplitudes of the frequencies
!! note that the resulting time series has zero mean.
SUBROUTINE CalcTargetPSD(p, S, U, ErrStat, ErrMsg)

   TYPE(TurbSim_ParameterType),  INTENT(in)     :: p                            !< TurbSim parameters
   REAL(ReKi),                   INTENT(in)     :: U           (:)              !< The steady u-component wind speeds for the grid (NPoints).
   REAL(ReKi),                   INTENT(  OUT)  :: S           (:,:,:)          !< The turbulence PSD array (NumFreq,NPoints,3).

   INTEGER(IntKi),               INTENT(  out)  :: ErrStat                      !< Error level
   CHARACTER(*),                 INTENT(  out)  :: ErrMsg                       !< Message describing error
   
   
   ! local variables
   
   INTEGER(IntKi)                   :: IFreq                            ! Index for frequency
   INTEGER(IntKi)                   :: LastIndex(2)                     ! Index for the last (Freq, Ht) used in models that interpolate/extrapolate user-input spectra or time series
   
   INTEGER(IntKi)                   :: IVec                             ! loop counter for velocity components
   INTEGER(IntKi)                   :: IPoint, iPointUsr                ! loop counter for grid points
   
   REAL(ReKi),   ALLOCATABLE        :: SSVS (:,:)                       ! A temporary work array (NumFreq,3) that holds a single-sided velocity spectrum.
   REAL(ReKi)                       :: DUDZ                             ! The steady u-component wind shear for the grid  [used in Hydro models only].
   REAL(ReKi)                       :: ZTmp, UTmp                       ! temporary height and velocity used for finite difference calculations
   
   REAL(ReKi)                       :: HalfDelF                         ! half of the delta frequency, used to discretize the continuous PSD at each point
      
   !INTEGER(IntKi)                   :: UP                               ! I/O unit for PSD debugging file.   
   !CHARACTER(200)                   :: FormStr                          ! String used to store format specifiers for PSD debugging.
   
   INTEGER(IntKi)                   :: ErrStat2                         ! Error level (local)
   CHARACTER(MaxMsgLen)             :: ErrMsg2                          ! Message describing error (local)
   CHARACTER(*), PARAMETER          :: RoutineName = 'CalcTargetPSD'

      ! initialize variables

   ErrStat = ErrID_None
   ErrMsg  = ""


      !  Allocate the array to hold the single-sided velocity spectrum.

   CALL AllocAry( SSVS, p%grid%NumFreq,3, 'SSVS', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF



      ! Calculate the single point Power Spectral Densities. 

   HalfDelF = 0.5*p%grid%Freq(1) 
   
   
   SELECT CASE ( p%met%TurbModel_ID )
      CASE ( SpecModel_IECKAI ) ! IECKAI has uniform spectra (does not vary with height or velocity)
         CALL Spec_IECKAI  ( p%UHub, p%IEC%SigmaIEC, p%IEC%IntegralScale, p%grid%Freq, p%grid%NumFreq, SSVS )
      
         DO IVec=1,3
            DO IFreq=1,p%grid%NumFreq
               S(IFreq,:,IVec) = SSVS(IFreq,IVec)*HalfDelF
            END DO ! IFreq
         END DO ! IVec

      
      CASE ( SpecModel_IECVKM )  ! IECVKM has uniform spectra (does not vary with height or velocity)
         CALL Spec_IECVKM  ( p%UHub, p%IEC%SigmaIEC(1), p%IEC%IntegralScale, p%grid%Freq, p%grid%NumFreq, SSVS )
   
         DO IVec=1,3
            DO IFreq=1,p%grid%NumFreq
               S(IFreq,:,IVec) = SSVS(IFreq,IVec)*HalfDelF
            END DO ! IFreq
         END DO ! IVec       
         
                        
      CASE ( SpecModel_API )
         DO IPoint=1,p%grid%NPoints
            CALL Spec_API ( p, p%grid%Z(IPoint), SSVS )
            S(:,IPoint,:) = SSVS*HalfDelF               
         ENDDO             
      
         
      CASE ( SpecModel_GP_LLJ )
         IF ( ALLOCATED( p%met%ZL_profile ) ) THEN !.AND. ALLOCATED( p%met%Ustar_profile ) )  THEN  
            DO IPoint=1,p%grid%NPoints
               CALL Spec_GPLLJ   ( p, p%grid%Z(IPoint), U(IPoint), p%met%ZL_profile(IPoint), p%met%Ustar_profile(IPoint), SSVS )
               S(:,IPoint,:) = SSVS*HalfDelF               
            ENDDO          
         ELSE
            DO IPoint=1,p%grid%NPoints
               CALL Spec_GPLLJ   ( p, p%grid%Z(IPoint), U(IPoint), p%met%ZL,                 p%met%Ustar,                 SSVS )
               S(:,IPoint,:) = SSVS*HalfDelF               
            ENDDO          
         ENDIF
                  
         
      CASE (SpecModel_NWTCUP)
         DO IPoint=1,p%grid%NPoints
            CALL Spec_NWTCUP  ( p, p%grid%Z(IPoint), U(IPoint), SSVS )
            S(:,IPoint,:) = SSVS*HalfDelF               
         ENDDO         
         
      
      CASE ( SpecModel_SMOOTH )
         DO IPoint=1,p%grid%NPoints
            CALL Spec_SMOOTH   ( P, p%grid%Z(IPoint), U(IPoint), SSVS )
            S(:,IPoint,:) = SSVS*HalfDelF               
         ENDDO         
         
      
      CASE ( SpecModel_TIDAL, SpecModel_RIVER )
         DO IPoint=1,p%grid%NPoints
            ZTmp = p%grid%Z(IPoint) + p%grid%GridRes_Z
            
            CALL getVelocity(p, p%UHub,p%grid%HubHt, ZTmp, UTmp, ErrStat2, ErrMsg2)   !get velocity Utmp at height ZTmp
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)            
            
            DUDZ = ( UTmp - U(IPoint) ) / p%grid%GridRes_Z
            CALL Spec_TIDAL  ( p, p%grid%Z(IPoint), DUDZ, SSVS, p%met%TurbModel_ID )
         
                  ! Discretize the continuous PSD and store it in matrix "S"
           
            S(:,IPoint,:) = SSVS*HalfDelF               
         ENDDO 
            
         
      CASE ( SpecModel_USER ) ! currently is uniform spectra
         CALL Spec_UserSpec   ( p, SSVS )
      
         DO IVec=1,3
            DO IFreq=1,p%grid%NumFreq
               S(IFreq,:,IVec) = SSVS(IFreq,IVec)*HalfDelF
            END DO ! IFreq
         END DO ! IVec         
         
         
      CASE ( SpecModel_TimeSer ) 
         
         DO iPointUsr = 1,p%usr%NPoints
      
            iPoint = iPointUsr 
            IF (iPointUsr == p%usr%RefPtID ) THEN
               iPoint = p%usr%NPoints !this is the point on the grid
            ELSEIF (iPointUsr  == p%usr%NPoints) THEN
               iPoint = p%usr%RefPtID
            END IF
       !bjj: make sure size(ssvs,1) = p%grid%NumFreq <= p%usr%nFreq = size(p%usr%S,1)
            ! initialize SSVS with extrapolated values if there are non-specified components or frequencies 
            ! i.e., fill the gaps where wind component or frequencies exceed what was specified in the time-series data with some numerical model
            ! (use zeros for known spectral values that will get overwritten later)
            CALL Spec_TimeSer_Extrap ( p, p%grid%Z(iPoint), U(iPoint), SSVS )
            
            ! overwrite the frequencies and wind components that were computed from measurements in the time-series file
            SSVS(1:p%usr%nFreq,1:p%usr%nComp) = p%usr%S(1:p%usr%nFreq,iPointUsr,1:p%usr%nComp)
                        
            S(:,iPoint,:) = SSVS*HalfDelF
         END DO
         
         
         DO iPoint=1+p%usr%NPoints,p%grid%NPoints
            CALL Spec_TimeSer( p, p%grid%Z(iPoint), U(iPoint), LastIndex, SSVS )
            S(:,iPoint,:) = SSVS*HalfDelF
         ENDDO                  
         
      CASE ( SpecModel_USRVKM )
         DO iPoint=1,p%grid%NPoints
            CALL Spec_vonKrmn   ( P, p%grid%Z(iPoint), U(iPoint), SSVS )
            S(:,iPoint,:) = SSVS*HalfDelF               
         ENDDO       
         
         
      CASE (SpecModel_WF_UPW)
         DO IPoint=1,p%grid%NPoints
            CALL Spec_WF_UPW  ( p, p%grid%Z(IPoint), U(IPoint), SSVS )
            S(:,IPoint,:) = SSVS*HalfDelF               
         ENDDO       
         
         
      CASE ( SpecModel_WF_07D, SpecModel_WF_14D )
         DO IPoint=1,p%grid%NPoints
            CALL Spec_WF_DW ( p, p%grid%Z(IPoint), U(IPoint), SSVS, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            S(:,IPoint,:) = SSVS*HalfDelF               
         ENDDO       
         
         
      CASE ( SpecModel_NONE )
         S = 0.0_ReKi ! whole matrix is zero
   !bjj TEST: CALL Spec_Test ( p%grid%Z(IPoint), U(IPoint), SSVS )

      CASE ( SpecModel_MODVKM )
         IF (MVK) THEN
         !  DO IPoint=1,p%grid%NPoints
         !     CALL Mod_vKrm( p%grid%Z(IPoint), U(IPoint), SSVS )
         !     S(:,IPoint,:) = SSVS*HalfDelF               
         !  ENDDO       
         ELSE
            CALL SetErrStat( ErrID_Fatal, 'Specified turbulence PSD, "'//TRIM( p%met%TurbModel )//'", not availible.', ErrStat, ErrMsg, RoutineName)
            CALL Cleanup()
            RETURN
         ENDIF

      CASE DEFAULT
         CALL SetErrStat( ErrID_Fatal, 'Specified turbulence PSD, "'//TRIM( p%met%TurbModel )//'", not availible.', ErrStat, ErrMsg, RoutineName)
         CALL Cleanup()
         RETURN
      END SELECT          
   
         
   !IF (PSD_OUT) THEN
   !   UP = -1
   !   CALL GetNewUnit( UP, ErrStat2, ErrMsg2 )
   !      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !   CALL OpenFOutFile ( UP, TRIM( p%RootName )//'.psd', ErrStat2, ErrMsg2)
   !      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !      IF (ErrStat >= AbortErrLev) THEN
   !         CALL Cleanup()
   !         RETURN
   !      END IF   
   !   
   !   WRITE (UP,"(A)")  'PSDs '
   !   WRITE (UP, "( A4,'"//TAB//"',A4,"//TRIM( Int2LStr( p%grid%NumFreq ) )//"('"//TAB//"',G10.4) )")  'Comp','Ht', p%grid%Freq(:)
   !   FormStr  = "( I4,"//TRIM( Int2LStr( p%grid%NumFreq+1 ) )//"('"//TAB//"',G10.4) )"
   !   
   !   DO IPoint=1,p%grid%NPoints           
   !                                                                         
   !      !IF ( ABS(Ht - p%grid%HubHt) < Tolerance ) THEN
   !         WRITE( UP, FormStr ) 1, p%grid%Z(IPoint), S(:,IPoint,1)/HalfDelF
   !         WRITE( UP, FormStr ) 2, p%grid%Z(IPoint), S(:,IPoint,2)/HalfDelF
   !         WRITE( UP, FormStr ) 3, p%grid%Z(IPoint), S(:,IPoint,3)/HalfDelF
   !      !ENDIF
   !   
   !   ENDDO ! IPoint
   !
   !   CLOSE( UP )
   !ENDIF
   

   CALL Cleanup()
   RETURN
   
CONTAINS 
!....................................
   SUBROUTINE Cleanup()

      !IF ( PSD_OUT .AND. UP > 0) CLOSE( UP )
         
      IF ( ALLOCATED( SSVS  ) )  DEALLOCATE( SSVS  )

   END SUBROUTINE Cleanup
END SUBROUTINE CalcTargetPSD
!=======================================================================
!> This routine creates the grid (cartesian + other points) that are
!!  to be simulated.
SUBROUTINE CreateGrid( p_grid, p_usr, UHub, AddTower, ErrStat, ErrMsg )

! Assumes that these variables are set:
!  GridHeight
!  GridWidth 
!  NumGrid_Y
!  NumGrid_Z

   TYPE(Grid_ParameterType),        INTENT(INOUT) :: p_grid
   TYPE(UserTSSpec_ParameterType),  INTENT(INOUT) :: p_usr
   
   REAL(ReKi)                     , INTENT(IN   ) :: UHub                            ! Mean wind speed at hub, used only when usable time is not "ALL" (i.e., periodic flag is false) 
   LOGICAL                        , INTENT(INOUT) :: AddTower                        ! Value of p%WrFile(FileExt_TWR) [determines if tower points should be generarated]

   INTEGER(IntKi),                  intent(  out) :: ErrStat                         ! Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg                          ! Message describing error
   
   ! local variables:
   REAL(DbKi)                                     :: DelF                            ! Delta frequency
   INTEGER(IntKi)                                 :: IY, IZ, IFreq                   ! loop counters 
   INTEGER(IntKi)                                 :: NTwrPts                         ! number of extra tower points
   INTEGER(IntKi)                                 :: NTwrIndx                        ! number of tower points to be placed in output file
   
   INTEGER(IntKi)                                 :: TmpIndex                        ! temporary index
      
   INTEGER(IntKi)                                 :: HubIndx_Y                       ! Index into Y dimension of grid for hub location
   INTEGER(IntKi)                                 :: HubIndx_Z                       ! Index into Z dimension of grid for hub location
   
   INTEGER(IntKi)                                 :: NumSteps2                       ! one-half the number of steps
   INTEGER(IntKi)                                 :: iPoint, iPointUsr               ! loop counter for points
   
   INTEGER(IntKi)                                 :: ErrStat2                        ! Error level (local)
   CHARACTER(MaxMsgLen)                           :: ErrMsg2                         ! Message describing error (local)
   CHARACTER(*), PARAMETER                        :: RoutineName = 'CreateGrid'
   LOGICAL                                        :: GenerateExtraHubPoint
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
   !.....................................................
   ! First, let's deal with time and frequencies:
   !.....................................................
   
      ! Calculate Total time and NumSteps.
      ! Find the product of small factors that is larger than NumSteps (prime #9 = 23).
!bjj: I have no idea why it is necessary to be a factor of 4, so I'm removing it for now:      ! Make sure it is a multiple of 2 too.

   IF ( p_grid%Periodic ) THEN
      p_grid%NumSteps    = CEILING( p_grid%AnalysisTime / p_grid%TimeStep )

         ! make sure NumSteps is an even number and a product of small primes
      NumSteps2          = ( p_grid%NumSteps - 1 )/2 + 1
      p_grid%NumSteps    = 2*PSF( NumSteps2 , 9 )  ! >= 2*NumSteps2 = NumSteps + 1 - MOD(NumSteps-1,2) >= NumSteps
      !p_grid%NumSteps    = PSF( p_grid%NumSteps , 9 )  
      
      p_grid%NumOutSteps = p_grid%NumSteps
   ELSE
      p_grid%NumOutSteps = CEILING( ( p_grid%UsableTime + p_grid%GridWidth / UHub )/p_grid%TimeStep )
      p_grid%NumSteps    = MAX( CEILING( p_grid%AnalysisTime / p_grid%TimeStep ), p_grid%NumOutSteps )
      
         ! make sure NumSteps is an even number and a product of small primes      
!      p_grid%NumSteps    = PSF( p_grid%NumSteps , 9 )  ! make sure it's a product of small primes
      NumSteps2          = ( p_grid%NumSteps - 1 )/2 + 1
      p_grid%NumSteps    = 2*PSF( NumSteps2 , 9 )  ! >= 2*NumSteps2 = NumOutSteps + 1 - MOD(NumOutSteps-1,2) >= NumOutSteps
      
   END IF

   !IF (p_grid%NumSteps < 2 )  THEN
   !   CALL SetErrStat( ErrID_Fatal, 'There must be at least 2 time steps. '//&
   !                    'Increase the usable length of the time series or decrease the time step.', ErrStat, ErrMsg, RoutineName )
   !   RETURN
   !END IF
   
   p_grid%NumFreq = p_grid%NumSteps / 2
   DelF           = 1.0_DbKi/( p_grid%NumSteps*p_grid%TimeStep )
      
      ! quick check that the frequency contents are the same as the user-input time series. (necessary because we want to keep the exact time series)  
   IF (p_usr%NPoints > 0) THEN
      ! IF ( .NOT. EqualRealNos( DelF, p_usr%f(1) ) .or. .not. EqualRealNos(p_grid%AnalysisTime,p_usr%f(1)*p_usr%NFreq ) .or. p_grid%NumFreq > size(p_usr%f) ) THEN
      IF ( .NOT. EqualRealNos( DelF, p_usr%DelF ) ) THEN
         CALL SetErrStat(ErrID_Fatal, 'Delta frequency in the user-input time series must be the same as the delta frequency in the simulated series. '//&
            'Change AnalysisTime or number of rows entered in user-defined time series file.', ErrStat, ErrMsg,RoutineName)
         RETURN
      END IF      
      
      if ( p_grid%NumFreq > p_usr%nFreq ) then
         CALL SetErrStat(ErrID_Fatal, 'Cannot output more frequency values than were entered in user-defined time series file.', ErrStat, ErrMsg,RoutineName)
         RETURN         
      end if
      
   END IF
   
            
   CALL AllocAry( p_grid%Freq, p_grid%NumFreq, 'Freq (frequency array)', ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      IF (ErrStat >= AbortErrLev) RETURN

   DO IFreq=1,p_grid%NumFreq
      p_grid%Freq(IFreq) = IFreq*DelF
   ENDDO 
   
   
   !.....................................................
   ! Now, figure out the points in space:
   !  1) user-specified time-series points
   !  2) regularly spaced y-z grid
   !  3) hub point
   !  4) lollipop-stick tower points
   !.....................................................
   
   ! start by determining how many points will be in the output files
   ! (1) the full-field grid:
   p_grid%GridRes_Y = p_grid%GridWidth  / REAL( p_grid%NumGrid_Y - 1, ReKi )
   p_grid%GridRes_Z = p_grid%GridHeight / REAL( p_grid%NumGrid_Z - 1, ReKi )      

   p_grid%Zbottom = p_grid%HubHt + 0.5*p_grid%RotorDiameter                                ! height of the highest grid points
   p_grid%Zbottom = p_grid%Zbottom - p_grid%GridRes_Z * REAL(p_grid%NumGrid_Z - 1, ReKi)   ! height of the lowest grid points

   IF ( p_grid%Zbottom <= 0.0_ReKi ) THEN
      CALL SetErrStat(ErrID_Fatal,'The lowest grid point ('//TRIM(Num2LStr(p_grid%Zbottom))// ' m) must be above the ground. '//&
                      'Adjust the appropriate values in the input file.',ErrStat,ErrMsg,RoutineName)
      RETURN
   ENDIF   
   

   ! (2) the tower points:   
   IF ( AddTower ) THEN

      IF ( MOD(p_grid%NumGrid_Y, 2) == 0 ) THEN
         p_grid%ExtraTwrPT = .TRUE.
      ELSE
         p_grid%ExtraTwrPT = .FALSE.
      END IF
         
         ! Compute the number of points between the bottom of the grid and the ground 
         ! ( but we don't want to be on the ground, just more than "Tolerance" from it )
 
      NTwrPts  = INT( ( p_grid%Zbottom - Tolerance ) / p_grid%GridRes_Z )
      NTwrIndx = NTwrPts + 1

      IF ( NTwrPts < 1 ) THEN        
         CALL SetErrStat(ErrID_Warn, ' There are no extra tower data points below the grid. Tower output will be turned off.',ErrStat,ErrMsg,RoutineName)
         AddTower = .FALSE.  ! bjj: change this so it doesn't actually modify this variable inside this routine???
         NTwrPts  = 0
         NTwrIndx = 0
      ENDIF
            
      IF ( p_grid%ExtraTwrPT ) THEN 
         NTwrPts = NTwrPts + 1  ! Let's add the point on the bottom of the grid so tower interpolation is easier in AeroDyn
      ENDIF
      
   ELSE
      NTwrPts  = 0
      NTwrIndx = 0
   ENDIF   
   
         ! we will set these index arrays to point to the grid/tower points or user-specified points
   CALL AllocAry(p_grid%GridPtIndx,p_grid%NumGrid_Y*p_grid%NumGrid_Z, 'GridPtIndx', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry(p_grid%TwrPtIndx, NTwrIndx,                          'TwrPtIndx',  ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
   IF (ErrStat >= AbortErrLev) RETURN

   
   !...............
   ! Now, let's see how many points we're going to simulate
   !..............
   
      ! here's our first estimate of how many points there will be. Later we will add a point for the hub if    
      ! necessary and subtract points from the grid or tower that are duplicates of the user-specified ones.   
   p_grid%NPoints =  p_usr%NPoints                       &   ! (1) the user-specified time-series points
                   + p_grid%NumGrid_Y*p_grid%NumGrid_Z   &   ! (2) the rectangular grid
                   + NTwrPts                                 ! (4) the tower points (the stick of the lollipop)
   
   ! Check if any of the user-specified time-series points are duplicated elsewhere:
   p_grid%GridPtIndx = 0
   p_grid%TwrPtIndx  = 0
   
   ! (2) the rectangular grid:
   DO iPointUsr = 1,p_usr%NPoints
      
      iPoint = iPointUsr
      IF (iPointUsr == p_usr%RefPtID ) THEN
         iPoint = p_usr%NPoints
      ELSEIF (iPointUsr == p_usr%NPoints) THEN
         iPoint = p_usr%RefPtID
      END IF
      
         ! Is this point on the regularly-spaced grid?   
      TmpIndex = IndexOnGrid( p_grid, p_usr%pointyi(iPointUsr), p_usr%pointzi(iPointUsr) )
      
      IF ( TmpIndex > 0 ) THEN
         p_grid%GridPtIndx( TmpIndex ) = iPoint
         ! it's a duplicate of a point on the rectangular grid, so subtract one from NPoints:
         p_grid%NPoints = p_grid%NPoints - 1
      ELSE
            ! Is this point on the tower?   
         IF ( NTwrPts > 0 ) THEN
            TmpIndex = IndexOnTower( p_grid, p_usr%pointyi(iPointUsr), p_usr%pointzi(iPointUsr) )
            
            IF ( TmpIndex > 0 ) THEN
               p_grid%TwrPtIndx( TmpIndex ) = iPoint    
               
               ! it's a duplicate of a tower point, so subtract one from NPoints:
               p_grid%NPoints = p_grid%NPoints - 1               
            END IF                
         END IF  ! NTwrPts > 0 
         
      END IF
      
   END DO
                  
                 
   ! (3) the hub point:
         
   IF ( MOD(p_grid%NumGrid_Y, 2) == 0 ) THEN
      
      p_grid%HubOnGrid  = .FALSE.
      
   ELSE
         ! This is the hub Z index if it falls on the grid
      HubIndx_Z = INT( Tolerance + ( p_grid%HubHt - p_grid%Zbottom ) / p_grid%GridRes_Z ) + 1 
      
      IF ( ABS((HubIndx_Z-1)*p_grid%GridRes_Z + p_grid%Zbottom - p_grid%HubHt) > Tolerance ) THEN
         p_grid%HubOnGrid = .FALSE.
      ELSE
         p_grid%HubOnGrid = .TRUE.
      END IF
      
   END IF

   p_grid%HubIndx = 0
   IF ( .NOT. p_grid%HubOnGrid ) THEN
      GenerateExtraHubPoint = .TRUE.      
      ! Is it a user-defined point?
      DO iPointUsr = 1,p_usr%NPoints
            
         IF ( EqualRealNos( p_usr%pointyi(iPointUsr), 0.0_ReKi ) .AND. EqualRealNos( p_usr%pointzi(iPointUsr), p_grid%HubHt ) ) THEN
                        
            IF (iPointUsr == p_usr%RefPtID ) THEN
               iPoint = p_usr%NPoints
            ELSEIF (iPointUsr == p_usr%NPoints) THEN
               iPoint = p_usr%RefPtID
            ELSE
               iPoint = iPointUsr
            END IF
                     
            p_grid%HubIndx = iPoint
            GenerateExtraHubPoint = .FALSE.
            EXIT  ! we found it
         END IF         
      END DO                  
   ELSE
      GenerateExtraHubPoint = .FALSE.
   END IF
   
   IF (GenerateExtraHubPoint) p_grid%NPoints = p_grid%NPoints + 1  

   p_grid%NPacked   = p_grid%NPoints*( p_grid%NPoints + 1 )/2    ! number of entries stored in the packed version of the symmetric matrix of size NPoints by NPoints

   
   ! we now know how many points there are going to be, so let's create the arrays that contains their locations and finish updating our index arrays
   
   CALL AllocAry(p_grid%Y, p_grid%NPoints, 'Y (lateral locations of the grid points)',  ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry(p_grid%Z, p_grid%NPoints, 'Z (vertical locations of the grid points)', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   IF (ErrStat >= AbortErrLev) RETURN
   
   ! (1) User-defined points 
   DO iPointUsr = 1,p_usr%NPoints
      
      iPoint = iPointUsr
      IF (iPointUsr == p_usr%RefPtID ) THEN
         iPoint = p_usr%NPoints
      ELSEIF (iPointUsr == p_usr%NPoints) THEN
         iPoint = p_usr%RefPtID
      END IF
      
      p_grid%y(iPoint) = p_usr%pointyi(iPointUsr)
      p_grid%z(iPoint) = p_usr%pointzi(iPointUsr)      
   END DO   

   
   ! (2) rectangular y-z grid:
   iPoint = p_usr%NPoints
   DO IZ = 1,p_grid%NumGrid_Z
      DO IY = 1,p_grid%NumGrid_Y
         
         TmpIndex = (IZ-1)*p_grid%NumGrid_Y + IY
         
         IF ( p_grid%GridPtIndx(TmpIndex) < 1 ) THEN ! we didn't find this grid point in the set of user-defined points, so create a new point
            iPoint = iPoint + 1
            
            p_grid%Y(iPoint)            = -0.5*p_grid%GridWidth  + p_grid%GridRes_Y*( IY - 1 )
            p_grid%Z(iPoint)            =      p_grid%Zbottom    + p_grid%GridRes_Z*( IZ - 1 )        
            p_grid%GridPtIndx(TmpIndex) = iPoint            
         END IF
         
      END DO
   END DO   
   
   ! note: GridPtIndx should be completely set now.
   
   ! (3) hub point:
         
   IF ( p_grid%HubOnGrid ) THEN
            
      HubIndx_Y = INT( ( p_grid%NumGrid_Y + 1 ) / 2 ) ! the center point      
      p_grid%HubIndx   = p_grid%GridPtIndx( p_grid%NumGrid_Y*( HubIndx_Z - 1 ) + HubIndx_Y )           
      
   ELSEIF ( GenerateExtraHubPoint ) THEN
      iPoint = iPoint + 1  
      
      p_grid%Y(iPoint) = 0.0_ReKi
      p_grid%Z(iPoint) = p_grid%HubHt        
      p_grid%HubIndx   = iPoint
      
   ! ELSE -> HubIndx is set already
   ENDIF

   
   ! (4) Finally, let's deal with the tower "lollipop" points:  
      
   IF ( AddTower ) THEN !p%WrFile(FileExt_TWR)
      
      IF ( .NOT. p_grid%ExtraTwrPT ) THEN
         p_grid%TwrPtIndx(1)  = p_grid%GridPtIndx( INT(p_grid%NumGrid_Y / 2) + 1 ) ! center y location on bottom height
      END IF
      
      
      DO IZ = 1,NTwrIndx
         IF ( p_grid%TwrPtIndx(IZ) < 1 ) THEN
            iPoint = iPoint + 1
            
            p_grid%Y(iPoint)     = 0.0_ReKi
            p_grid%Z(iPoint)     = p_grid%ZBottom - (IZ-1)*p_grid%GridRes_Z        
            p_grid%TwrPtIndx(IZ) = iPoint
            
         END IF
      END DO     

   ENDIF   
         
END SUBROUTINE CreateGrid
!=======================================================================
!> This routine determines if a point at location (y,z) is
!! on the regularly-spaced y-z grid. If it does, it returns the
!! index of the point on the grid. If it does not, it returns -1.
FUNCTION IndexOnGrid( p_grid, y, z )

   TYPE(Grid_ParameterType),        INTENT(IN) :: p_grid          !< grid parameters
   REAL(ReKi),                      INTENT(IN) :: y               !< y position of point we're querying
   REAL(ReKi),                      INTENT(IN) :: z               !< z position of point we're querying
   
   INTEGER(IntKi)                              :: IndexOnGrid     !< Index on regularly spaced grid
   
      ! local variables
   INTEGER(IntKi)                              :: YIndx           !  Index on regularly spaced grid
   INTEGER(IntKi)                              :: ZIndx           !  Index on regularly spaced grid
   REAL(ReKi)                                  :: y1              !  left-most location on grid

      
   y1 = -0.5_ReKi * p_grid%GridWidth
      
   ZIndx = INT( Tolerance + ( z - p_grid%Zbottom ) / p_grid%GridRes_Z ) + 1 
            
   IF ( .NOT. EqualRealNos( p_grid%Zbottom + (ZIndx-1)*p_grid%GridRes_Z , z) ) THEN
      IndexOnGrid = -1
      RETURN
   END IF
   
   
   
   YIndx = INT( Tolerance + ( y - y1             ) / p_grid%GridRes_Y ) + 1 
   
   IF ( .NOT. EqualRealNos( y1 + (YIndx-1)*p_grid%GridRes_Y , y) ) THEN
      IndexOnGrid = -1
      RETURN
   END IF
   
   IF ( YIndx < 1 .OR. YIndx > p_grid%NumGrid_Y .OR. &
        ZIndx < 1 .OR. ZIndx > p_grid%NumGrid_Z ) THEN
      IndexOnGrid = -1
      RETURN
   END IF
      
   
   IndexOnGrid = (ZIndx-1)*p_grid%NumGrid_Y + YIndx

END FUNCTION IndexOnGrid
!=======================================================================
!> This routine determines if a point at location (y,z) is
!! on the regularly-spaced y-z grid. If it does, it returns the
!! index of the point on the grid. If it does not, it returns -1.
FUNCTION IndexOnTower( p_grid, y, z )

   TYPE(Grid_ParameterType),        INTENT(IN) :: p_grid          !< grid parameters
   REAL(ReKi),                      INTENT(IN) :: y               !< y position of point we're querying
   REAL(ReKi),                      INTENT(IN) :: z               !< z position of point we're querying
   
   INTEGER(IntKi)                              :: IndexOnTower    !< Index on regularly spaced tower points
   
      ! local variables
   INTEGER(IntKi)                              :: ZIndx           !  Index on regularly spaced grid

      
   IF ( .NOT. EqualRealNos( 0.0_ReKi , y) ) THEN
      IndexOnTower = -1
      RETURN
   END IF   
      
   ZIndx = INT( Tolerance + ( p_grid%Zbottom - z ) / p_grid%GridRes_Z ) + 1 
            
   IF ( zIndx < 0 .OR. .NOT. EqualRealNos( p_grid%Zbottom - (ZIndx-1)*p_grid%GridRes_Z , z) ) THEN
      IndexOnTower = -1
      RETURN
   END IF           
   
   IndexOnTower = ZIndx


END FUNCTION IndexOnTower
!=======================================================================
!> This routine calculates the wind components in the Inertial reference
!!  frame.
SUBROUTINE SetPhaseAngles( p, OtherSt_RandNum, PhaseAngles, ErrStat, ErrMsg )


   TYPE(TurbSim_ParameterType),  INTENT(IN   ) :: p                                              !< parameters for TurbSim
   TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt_RandNum                                !< other states for random number generation
   INTEGER(IntKi)  ,             INTENT(  OUT) :: ErrStat                                        !< error level/status
   CHARACTER(*) ,                INTENT(  OUT) :: ErrMsg                                         !< error message
                                                                                              
   REAL(ReKi)                  , INTENT(  OUT) :: PhaseAngles(p%grid%NPoints,p%grid%NumFreq,3)   !< phases

      ! local variables
   INTEGER(IntKi)                              :: iPoint, iPointUsr                              ! points that have phases defined already 
   

      ! generate random phases for all the points
      
      ! bjj: todo: don't generate the angles for user-specified time-series points, which have phases already
   CALL RndPhases(p%RNG, OtherSt_RandNum, PhaseAngles, p%grid%NPoints, p%grid%NumFreq, p%US, ErrStat, ErrMsg)

   
   IF (p%met%TurbModel_ID == SpecModel_TimeSer) THEN
      
         ! note: setting the phase angles this way assumes that p%usr%f(1:p%usr%nFreq) = p%grid%f(1:p%usr%nFreq) [i.e., TMax, AnalysisTime are equal]; 
         ! however, the simulated time series may have more frequencies and/or smaller time step than the user time-series input file.
      DO iPointUsr = 1,p%usr%NPoints
      
         iPoint = iPointUsr
         IF (iPointUsr == p%usr%RefPtID ) THEN
            iPoint = p%usr%NPoints
         ELSEIF (iPointUsr == p%usr%NPoints) THEN
            iPoint = p%usr%RefPtID
         END IF
         PhaseAngles(iPoint,1:p%usr%nFreq,1:p%usr%nComp) = p%usr%phaseAngles(:,iPointUsr,:)
      END DO
      
   END IF   
   
   ! nyquist frequency must be real, thus phase angle must be 0:
   PhaseAngles(:,p%grid%NumFreq,:) = 0.0_ReKi

END SUBROUTINE SetPhaseAngles
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
   
END SUBROUTINE CalculateStresses
!=======================================================================
!> Scale the velocity aligned along the sreamwise direction.
SUBROUTINE ScaleTimeSeries(p, V, ErrStat, ErrMsg)


   TYPE(TurbSim_ParameterType),     INTENT(IN)     ::  p                               !< TurbSim's parameters
   REAL(ReKi),                      INTENT(INOUT)  ::  V(:,:,:)                        !< velocity, aligned along the streamwise direction without mean values added 
   INTEGER(IntKi),                  intent(  out)  :: ErrStat                     !< Error level
   CHARACTER(*),                    intent(  out)  :: ErrMsg                      !< Message describing error


   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! Crossfeed cross-axis components to u', v', w' components and scale IEC models if necessary

   SELECT CASE ( p%met%TurbModel_ID )
   !MLB: There does not seem to be a CASE for TurbModel=="API".
      
   CASE (SpecModel_GP_LLJ, &
         SpecModel_NWTCUP, &
         SpecModel_SMOOTH, &
         SpecModel_WF_UPW, &
         SpecModel_WF_07D, &
         SpecModel_WF_14D, &
         SpecModel_USRVKM, &
         SpecModel_TIDAL,  &
         SpecModel_RIVER,  &
         SpecModel_USER    )
               
   
      CALL TimeSeriesScaling_ReynoldsStress(p, V, ErrStat, ErrMsg)


   CASE ( SpecModel_IECKAI , SpecModel_IECVKM ) ! API is considered an IEC model in this code, so it should fall here, if we wanted it

      CALL TimeSeriesScaling_IEC(p, V)
                       
   END SELECT   
   
END SUBROUTINE ScaleTimeSeries
!=======================================================================
!> This routine scales the time series so that the output has the exact
!! statistics desired. This scaling has the effect of changing the amplitude
!! of the target spectra to account for discretizing the spectra over a 
!! finite length of time.
SUBROUTINE TimeSeriesScaling_IEC(p, V)


   TYPE(TurbSim_ParameterType),     INTENT(IN)     ::  p                               !< TurbSim's parameters
   REAL(ReKi),                      INTENT(INOUT)  ::  V(:,:,:)                        !< velocity, aligned along the streamwise direction without mean values added 
   
   
   REAL(DbKi)                                      ::  CGridSum                        ! The sums of the velocity components at the points surrounding the hub (or at the hub if it's on the grid)
   REAL(DbKi)                                      ::  CGridSum2                       ! The sums of the squared velocity components at the points surrouding the hub 
   REAL(ReKi)                                      ::  UGridMean                       ! Average wind speed at a point 
   REAL(ReKi)                                      ::  UGridSig                        ! Standard deviation of the wind speed at a point 
   INTEGER(IntKi)                                  ::  IT                              ! loop counter (time)
   INTEGER(IntKi)                                  ::  Indx                            ! loop counter (grid point)
   INTEGER(IntKi)                                  ::  IVec                            ! loop counter (wind component)
   
   
   REAL(ReKi)                                      ::  ActualSigma(3)                  ! actual standard deviations
   REAL(ReKi)                                      ::  HubFactor(3)                    ! factor used to scale standard deviations at the hub point
   
   IF (p%IEC%ScaleIEC < 1) RETURN
   
   DO IVec = 1,3
      CGridSum  = 0.0
      CGridSum2 = 0.0
                           
      DO IT=1,p%grid%NumSteps !BJJ: NumOutSteps  -- scale to the output value?
         CGridSum  = CGridSum  + V( IT, p%grid%HubIndx, IVec )
         CGridSum2 = CGridSum2 + V( IT, p%grid%HubIndx, IVec )* V( IT, p%grid%HubIndx, IVec )
      ENDDO ! IT
               
      UGridMean          = CGridSum/p%grid%NumSteps !BJJ: NumOutSteps  -- scale to the output value?
      ActualSigma(IVec)  = SQRT( ABS( (CGridSum2/p%grid%NumSteps) - UGridMean*UGridMean ) )
            

      HubFactor(IVec) = p%IEC%SigmaIEC(IVec)/ActualSigma(IVec)  ! factor = Target / actual
                  
      IF (p%IEC%ScaleIEC == 1 .OR. p%met%SCMod(IVec) == CohMod_None) THEN ! with no coherence, all points have same std, so we'll save some calculations
               
         V(:,:,IVec) =     HubFactor(IVec) * V(:,:,IVec)
               
      ELSE  ! Scale each point individually
               
         DO Indx = 1,p%grid%NPoints
            CGridSum  = 0.0
            CGridSum2 = 0.0
                                    
            DO IT=1,p%grid%NumSteps !BJJ: NumOutSteps  -- scale to the output value?
               CGridSum  = CGridSum  + V( IT, Indx, IVec )
               CGridSum2 = CGridSum2 + V( IT, Indx, IVec )* V( IT, Indx, IVec )
            ENDDO ! IT
                  
            UGridMean = CGridSum/p%grid%NumSteps !BJJ: NumOutSteps  -- scale to the output value?
            UGridSig  = SQRT( ABS( (CGridSum2/p%grid%NumSteps) - UGridMean*UGridMean ) )
            
            V(:,Indx,IVec) = (p%IEC%SigmaIEC(IVec) / UGridSig) * V(:,Indx,IVec)   
         ENDDO ! Indx
               
      ENDIF            

   ENDDO !IVec   
   
   IF (p%US > 0 ) THEN
      WRITE( p%US, "(//,'Scaling statistics from the hub grid point:',/)" )                  
      WRITE( p%US, "(2X,'Component  Target Sigma (m/s)  Simulated Sigma (m/s)  Scaling Factor')" )
      WRITE( p%US, "(2X,'---------  ------------------  ---------------------  --------------')" )
         
      DO IVec = 1,3
         WRITE( p%US, "(5X,A,7x,f11.3,9x,f12.3,11x,f10.3)") Comp(IVec)//"'", p%IEC%SigmaIEC(IVec), &
                                                           ActualSigma(IVec), HubFactor(IVec) 
      END DO
   END IF
   

END SUBROUTINE TimeSeriesScaling_IEC
!=======================================================================
!> This routine performs a linear combination of the uncorrelated zero-mean
!! velocity aligned along the streamwise direction to obtain the desired 
!! Reynolds Stress values at the hub.
SUBROUTINE TimeSeriesScaling_ReynoldsStress(p, V, ErrStat, ErrMsg)

      ! passed variables
   TYPE(TurbSim_ParameterType), INTENT(IN)     ::  p                 !< parameters 
   REAL(ReKi),                  INTENT(INOUT)  ::  V(:,:,:)          !< velocity, aligned along the streamwise direction without mean values added
   INTEGER(IntKi),              intent(  out)  :: ErrStat            !< Error level
   CHARACTER(*),                intent(  out)  :: ErrMsg             !< Message describing error

   
      ! local variables
   REAL(DbKi)                                  ::  UVsum             ! The sum of the u'v' Reynolds stress component at the hub
   REAL(DbKi)                                  ::  UWsum             ! The sum of the u'w' Reynolds stress component at the hub
   REAL(DbKi)                                  ::  VWsum             ! The sum of the v'w' Reynolds stress component at the hub
   REAL(DbKi)                                  ::  UUsum             ! The sum of the u'u' Reynolds stress component at the hub
   REAL(DbKi)                                  ::  VVsum             ! The sum of the v'v' Reynolds stress component at the hub
   REAL(DbKi)                                  ::  WWsum             ! The sum of the w'w' Reynolds stress component at the hub
                                               
   REAL(ReKi)                                  ::  UVmean            ! The mean u'v' Reynolds stress component at the hub
   REAL(ReKi)                                  ::  UWmean            ! The mean u'w' Reynolds stress component at the hub
   REAL(ReKi)                                  ::  VWmean            ! The mean v'w' Reynolds stress component at the hub
   REAL(ReKi)                                  ::  UUmean            ! The mean u'u' Reynolds stress component at the hub
  !REAL(ReKi)                                  ::  VVmean            ! The mean v'v' Reynolds stress component at the hub
   REAL(ReKi)                                  ::  WWmean            ! The mean w'w' Reynolds stress component at the hub
                                               
   REAL(ReKi)                                  ::  alpha_uw          ! The coefficient of the u component added to the w component for correlation
   REAL(ReKi)                                  ::  alpha_uv          ! The coefficient of the u component added to the v component for correlation
   REAL(ReKi)                                  ::  alpha_wv          ! The coefficient of the w component added to the v component for correlation
               
   REAL(ReKi)                                  ::  u_indept          ! temporary copy of the uncorrelated u component of the velocity
   REAL(ReKi)                                  ::  v_indept          ! temporary copy of the uncorrelated v component of the velocity
   REAL(ReKi)                                  ::  w_indept          ! temporary copy of the uncorrelated w component of the velocity
   
   
   REAL(ReKi)                                  ::  INumSteps         ! Multiplicative Inverse of the Number of time Steps
   
   INTEGER(IntKi)                              ::  ITime             ! loop counter for time step/frequency 
   INTEGER(IntKi)                              ::  IPoint            ! loop counter for grid points
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
      !...................
      ! Calculate coefficients for obtaining "correct" Reynold's stresses at the hub
      !...................
            
      ! compute mean values:            
   UWsum = 0.0_DbKi
   UVsum = 0.0_DbKi
   VWsum = 0.0_DbKi
   UUSum = 0.0_DbKi
   VVSum = 0.0_DbKi
   WWSum = 0.0_DbKi
                           
   DO ITime = 1,p%grid%NumSteps
      UWsum = UWsum + V(ITime,p%grid%HubIndx,1) * V(ITime,p%grid%HubIndx,3)
      UVsum = UVsum + V(ITime,p%grid%HubIndx,1) * V(ITime,p%grid%HubIndx,2)
      VWsum = VWsum + V(ITime,p%grid%HubIndx,2) * V(ITime,p%grid%HubIndx,3)
      UUSum = UUSum + V(ITime,p%grid%HubIndx,1) * V(ITime,p%grid%HubIndx,1)
      !VVSum = VVSum + V(ITime,p%grid%HubIndx,2) * V(ITime,p%grid%HubIndx,2)
      WWSum = WWSum + V(ITime,p%grid%HubIndx,3) * V(ITime,p%grid%HubIndx,3)
   ENDDO
         
   INumSteps   = 1.0/p%grid%NumSteps
         
   UWmean = UWsum * INumSteps  
   UVmean = UVsum * INumSteps  
   VWmean = VWsum * INumSteps  
   UUmean = UUSum * INumSteps  
   !VVmean = VVSum * INumSteps  
   WWmean = WWSum * INumSteps  
            
      !BJJ: this is for v=alpha1, w=alpha2, u=alpha3 using derivation equations
   alpha_uw = ( p%met%PC_UW - UWmean ) / UUmean                                                              !alpha23
   alpha_wv = ( UUmean*(p%met%PC_VW - VWmean - alpha_uw*UVmean) - p%met%PC_UW*(p%met%PC_UV - UVmean) ) / &   !alpha12
               ( UUmean*(WWmean + alpha_uw*UWmean) - UWmean*p%met%PC_UW )
   alpha_uv = ( p%met%PC_UV - UVmean - alpha_wv*UWmean) / UUmean                                             !alpha13         
         
                  
      ! if we enter "none" for any of the Reynolds-stress terms, don't scale that component:
   IF (p%met%UWskip) alpha_uw = 0.0_ReKi
   IF (p%met%UVskip) alpha_uv = 0.0_ReKi
   IF (p%met%VWskip) alpha_wv = 0.0_ReKi
         
      !bjj: I'm implementing limits on the range of values here so that the spectra don't get too
      !     out of whack.  We'll display a warning in this case.
            
   IF ( ABS(alpha_uw) > 1.0 .OR. ABS(alpha_uv) > 1.0 .OR. ABS(alpha_wv) > 1.0 ) THEN
      ErrStat = ErrID_Info
      ErrMsg  = "Scaling terms exceed 1.0.  Reynolds stresses may be affected."
            
      alpha_uw = MAX( MIN( alpha_uw, 1.0_ReKi ), -1.0_ReKi )
      alpha_uv = MAX( MIN( alpha_uv, 1.0_ReKi ), -1.0_ReKi )
      alpha_wv = MAX( MIN( alpha_wv, 1.0_ReKi ), -1.0_ReKi )
                        
   ENDIF
                                    
      ! calculate the correlated time series:
            
   DO IPoint = 1,p%grid%NPoints
      DO ITime = 1, p%grid%NumSteps
         u_indept = V(ITime,IPoint,1)
         v_indept = V(ITime,IPoint,2)
         w_indept = V(ITime,IPoint,3)
                  
            ! equation 16 [PC_UW] in TurbSim user's guide v1.50
         V(ITime,IPoint,2) = alpha_uv*u_indept + v_indept + alpha_wv*w_indept 
         V(ITime,IPoint,3) = alpha_uw*u_indept +                     w_indept
               
      ENDDO          
   ENDDO   
   
   IF ( p%US > 0 ) THEN
         
      WRITE( p%US, "(//,'Scaling statistics from the hub grid point:',/)" )
      WRITE( p%US, "(3X,               'Cross-Component  Scaling Factor')" )
      WRITE( p%US, "(3X,               '---------------  --------------')" )
      WRITE( p%US, "(3X,A,2X,E14.5)" ) "u'w'           ", alpha_uw
      WRITE( p%US, "(3X,A,2X,E14.5)" ) "u'v'           ", alpha_uv
      WRITE( p%US, "(3X,A,2X,E14.5)" ) "v'w'           ", alpha_wv
         
   END IF    
         
END SUBROUTINE TimeSeriesScaling_ReynoldsStress
!=======================================================================
SUBROUTINE AddMeanAndRotate(p, V, U, HWindDir, VWindDir)

      ! passed variables
   TYPE(TurbSim_ParameterType), INTENT(IN)     ::  p                 !< parameters 
   REAL(ReKi),                  INTENT(INOUT)  ::  V(:,:,:)          !< velocity, on input: aligned along with the mean velocity without mean values added
                                                                     !! on output, aligned in the inertial reference frame with mean velocities added
   REAL(ReKi),                  INTENT(IN)     ::  U       (:)       !< profile of steady wind speed
   REAL(ReKi),                  INTENT(IN)     ::  HWindDir(:)       !< profile of horizontal wind direction
   REAL(ReKi),                  INTENT(IN)     ::  VWindDir(:)       !< profile of vertical wind direction
                                                                                                                                          
      ! local variables                                                                     
   REAL(ReKi)                                  ::  v3(3)             ! temporary 3-component array containing velocity
   INTEGER(IntKi)                              ::  ITime             ! loop counter for time step 
   INTEGER(IntKi)                              ::  IPoint            ! loop counter for grid points
                                                                     
                                                                     
                                                                     
                                                                     
   !..............................................................................
   ! Add mean wind to u' components and rotate to inertial reference  
   !  frame coordinate system
   !..............................................................................        
   DO IPoint=1,p%grid%Npoints   
      DO ITime=1,p%grid%NumSteps

            ! Add mean wind speed to the streamwise component and
            ! Rotate the wind to the X-Y-Z (inertial) reference frame coordinates:
            
         v3 = V(ITime,IPoint,:) 
         CALL CalculateWindComponents( v3, U(IPoint), HWindDir(IPoint), VWindDir(IPoint), V(ITime,IPoint,:) )
                        
      ENDDO ! ITime
         
   ENDDO ! IPoint   
                                                                     
                                                                     
END SUBROUTINE AddMeanAndRotate
!=======================================================================
SUBROUTINE TS_ValidateInput(P, ErrStat, ErrMsg)

   TYPE(TurbSim_ParameterType), INTENT(INOUT)     ::  p                              !< parameters 

   INTEGER(IntKi),                  intent(  out) :: ErrStat                         ! Error level
   CHARACTER(*),                    intent(  out) :: ErrMsg                          ! Message describing error

      ! local variables
   INTEGER(IntKi)                                 :: UnOut                           ! unit for output files
   INTEGER(IntKi)                                 :: ErrStat2                        ! Error level (local)
   CHARACTER(MaxMsgLen)                           :: ErrMsg2                         ! Message describing error (local)
   CHARACTER(*), PARAMETER                        :: RoutineName = 'TS_ValidateInput'
   

ErrStat = ErrID_None
ErrMsg  = ""


!BONNIE: UPPER LIMIT ON RICH_NO?
IF ( p%WrFile(FileExt_CTS) ) THEN
   
      ! models where coherent structures apply:
   IF ( p%met%TurbModel_ID == SpecModel_GP_LLJ .OR. &
        p%met%TurbModel_ID == SpecModel_NWTCUP .OR. &
        p%met%TurbModel_ID == SpecModel_SMOOTH .OR. &
        p%met%TurbModel_ID == SpecModel_WF_UPW .OR. &
        p%met%TurbModel_ID == SpecModel_WF_07D .OR. &
        p%met%TurbModel_ID == SpecModel_WF_14D ) THEN
      
      IF ( p%met%RICH_NO <= -0.05 ) THEN
         CALL SetErrStat( ErrID_Info, 'A coherent turbulence time step file cannot be generated for RICH_NO <= -0.05.', ErrStat, ErrMsg, RoutineName)
         p%WrFile(FileExt_CTS) = .FALSE.      
      ELSEIF ( .NOT. ( p%WrFile(FileExt_BTS) .OR. p%WrFile(FileExt_WND) ) ) THEN
         CALL SetErrStat( ErrID_Info, 'AeroDyn Full-Field files(.bts) will be generated along with the coherent turbulence file.', ErrStat, ErrMsg, RoutineName)
         p%WrFile(FileExt_BTS) = .TRUE.
      ENDIF      
      
   ELSE
      CALL SetErrStat( ErrID_Info, 'A coherent turbulence time step file cannot be generated with the '//TRIM(p%met%TurbModel)//' model.', ErrStat, ErrMsg, RoutineName)
      p%WrFile(FileExt_CTS) = .FALSE.      
   END IF
      
ENDIF !WrAct

   ! Make sure inputs make sense for API model
IF ( p%met%TurbModel_ID == SpecModel_API ) THEN
   IF ( .not. EqualRealNos( p%grid%AnalysisTime, 3600.0_ReKi) ) then ! warn if AnalysisTime isn't 1 hour (which is what we're assuming the URef is defined for)
      CALL SetErrStat( ErrID_Warn, 'API model assumes URef is the 1-hour mean wind speed, even though AnalysisTime does not equal 3600 seconds.', ErrStat, ErrMsg, RoutineName )      
   END IF
   
   IF ( p%IEC%ScaleIEC /= 0 ) THEN
      CALL SetErrStat( ErrID_Fatal, 'API model does not provide time-series scaling to target standard deviation. Set ScaleIEC = 0.', ErrStat, ErrMsg, RoutineName )      
   END IF

   IF ( .NOT. EqualRealNos( p%met%RefHt, 10.0_ReKi ) ) THEN
      CALL SetErrStat( ErrID_Fatal, 'API model requires a 10-m reference height. Set RefHt = 10.', ErrStat, ErrMsg, RoutineName )      
   END IF
      
END IF



      ! Warn if EWM is used with incompatible times
      
IF ( ( p%IEC%IEC_WindType == IEC_EWM1 .OR. p%IEC%IEC_WindType == IEC_EWM50 .OR. p%IEC%IEC_WindType == IEC_EWM100) .AND. & 
      ABS( 600.0_ReKi - MAX(p%grid%AnalysisTime,p%grid%UsableTime) ) > 90.0_ReKi ) THEN
   CALL SetErrStat( ErrID_Warn, 'The EWM parameters are valid for 10-min simulations only.', ErrStat, ErrMsg, RoutineName)   
ENDIF        


      ! Warn if Periodic is used with incompatible settings      
IF ( p%grid%Periodic .AND. .NOT. EqualRealNos(p%grid%AnalysisTime, p%grid%UsableTime) ) THEN
   CALL SetErrStat( ErrID_Warn, 'Periodic output files will not be generated when AnalysisTime /= UsableTime. Setting Periodic = .FALSE.', ErrStat, ErrMsg, RoutineName)  
   p%grid%Periodic = .FALSE.
END IF


   ! Warn if tower points are output but grid is not:
IF ( p%WrFile(FileExt_TWR) .AND. .NOT. ( p%WrFile(FileExt_WND) .OR. p%WrFile(FileExt_BTS)) ) THEN 
   CALL SetErrStat( ErrID_Info, 'TurbSim .bts file will be generated to contain the tower points.', ErrStat, ErrMsg, RoutineName)  
   p%WrFile(FileExt_BTS) = .TRUE.
END IF



   ! Open appropriate output files.  We will open formatted FF files later, if requested.
   ! Mention the files in the summary file.

IF ( ANY (p%WrFile) )  THEN
   CALL GetNewUnit( UnOut )

   WRITE (p%US,"( // 'You have requested that the following file(s) be generated:' / )")

   IF ( p%WrFile(FileExt_BIN) )  THEN   

!      CALL OpenBOutFile ( UnOut, TRIM( p%RootName)//'.bin', ErrStat, ErrMsg )
      CALL OpenUOutfile ( UnOut , TRIM( p%RootName)//'.bin', ErrStat2, ErrMsg2 )  ! just making sure it can be opened (not locked elsewhere)
      CLOSE(UnOut)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)  
      
      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".bin (binary hub-height turbulence-parameter file)' )")  

   ENDIF

   IF ( p%WrFile(FileExt_DAT) )  THEN     

      CALL OpenFOutFile ( UnOut, TRIM( p%RootName)//'.dat', ErrStat2, ErrMsg2 ) ! just making sure it can be opened (not locked elsewhere)
      CLOSE( UnOut )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)  
      
      WRITE (p%US, "( 3X ,'"//TRIM( p%RootName)//".dat (formatted turbulence-parameter file)' )")  

   ENDIF

   IF ( p%WrFile(FileExt_HH) )  THEN     

      CALL OpenFOutFile ( UnOut, TRIM( p%RootName)//'.hh', ErrStat2, ErrMsg2 ) ! just making sure it can be opened (not locked elsewhere)
      CLOSE( UnOut )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)  
      
      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".hh  (AeroDyn hub-height file)' )")

   ENDIF

   IF ( p%WrFile(FileExt_BTS) )  THEN

      CALL OpenBOutFile ( UnOut, TRIM(p%RootName)//'.bts', ErrStat2, ErrMsg2 ) ! just making sure it can be opened (not locked elsewhere)
      CLOSE( UnOut )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)  

      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".bts (AeroDyn/TurbSim full-field wind file)' )")  

   ENDIF

   IF ( p%WrFile(FileExt_WND) )  THEN

      CALL OpenBOutFile ( UnOut, TRIM(p%RootName)//'.wnd', ErrStat2, ErrMsg2 ) ! just making sure it can be opened (not locked elsewhere)
      CLOSE(UnOut)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)  

      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".wnd (AeroDyn/BLADED full-field wnd file)' )")  

   ENDIF
   
   IF ( p%WrFile(FileExt_TWR) .AND. p%WrFile(FileExt_WND) )  THEN

      CALL OpenBOutFile ( UnOut, TRIM( p%RootName )//'.twr', ErrStat2, ErrMsg2 ) ! just making sure it can be opened (not locked elsewhere)
      CLOSE(UnOut)       
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)  

      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".twr (binary tower file)' )")  

   ENDIF

   IF ( p%WrFile(FileExt_CTS) ) THEN      
      CALL OpenBOutFile ( UnOut, TRIM( p%RootName )//'.cts', ErrStat2, ErrMsg2 ) ! just making sure it can be opened (not locked elsewhere)
      CLOSE(UnOut)       
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)  
      
      
      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".cts (coherent turbulence time step file)' )")  
   ENDIF

   IF ( p%WrFile(FileExt_HAWC) )  THEN
      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//"-u.bin (HAWC binary full-field u-component file)' )")  

      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//"-v.bin (HAWC binary full-field v-component file)' )")  

      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//"-w.bin (HAWC binary full-field w-component file)' )")  
   ENDIF

   IF ( p%WrFile(FileExt_UVW) )  THEN
      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".u (formatted full-field U-component file)' )")  

      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".v (formatted full-field V-component file)' )")  

      WRITE (p%US,"( 3X ,'"//TRIM( p%RootName)//".w (formatted full-field W-component file)' )")  
   ENDIF

ELSE
   CALL SetErrStat( ErrID_Fatal, 'You have requested no output.', ErrStat, ErrMsg, RoutineName)   
ENDIF

   ! WARN if using a large grid and not creating ff output files
IF ( p%grid%NumGrid_Y*p%grid%NumGrid_Z > 250 ) THEN 
   IF (.NOT. p%WrFile(FileExt_WND) .AND. .NOT. p%WrFile(FileExt_BTS) .AND. .NOT. p%WrFile(FileExt_HAWC) .AND. .NOT. p%WrFile(FileExt_UVW) ) THEN
   
      CALL SetErrStat( ErrID_Warn, 'You are using a large number of grid points but are not generating full-field output files.'//&
            ' The simulation will run faster if you reduce the number of points on the grid.', ErrStat, ErrMsg, RoutineName) 
   END IF   
END IF

END SUBROUTINE TS_ValidateInput
!=======================================================================
SUBROUTINE TimeSeriesToSpectra( p, ErrStat, ErrMsg )

  ! USE NWTC_FFTPACK
      
   ! passed variables
   TYPE(TurbSim_ParameterType), INTENT(INOUT)  :: p

   INTEGER(IntKi),              intent(  out)  :: ErrStat                      !< Error level
   CHARACTER(*),                intent(  out)  :: ErrMsg                       !< Message describing error
   
   
   ! local variables
   TYPE(FFT_DataType)                          :: FFT_Data                     ! data for applying FFT
   real(siki),         allocatable             :: work (:)                     ! working array for converting fourier coefficients to spectra
   real(reki)                                  :: Re, Im                       ! real and imaginary parts of complex variable returned from fft
   
   REAL(ReKi)                                  :: meanU                        ! mean value of the U component
   REAL(ReKi)                                  :: meanV                        ! mean value of the V component
   REAL(ReKi)                                  :: meanW                        ! mean value of the W component

   REAL(DbKi)                                  :: DelF                         ! delta frequency

   REAL(ReKi)                                  :: cosH, sinH                   ! cosine and sin of horizontal angles
   REAL(ReKi)                                  :: cosV, sinV                   ! cosine and sin of vertical angles
   REAL(ReKi)                                  :: rotateMatrix(3,3)            ! rotation matrix to align with direction of mean velocity
   REAL(ReKi)                                  :: uvw(3)                       ! temporary array to hold 3 velocity components in case < 3 were entered
   
   INTEGER(IntKi)                              :: Indx                         ! generic index
   INTEGER(IntKi)                              :: iFreq                        ! loop counter for frequency 
   INTEGER(IntKi)                              :: iTime                        ! loop counter for time 
   INTEGER(IntKi)                              :: iVec                         ! loop counter for velocity components
   INTEGER(IntKi)                              :: iPoint                       ! loop counter for grid points
   INTEGER(IntKi)                              :: NumSteps                     ! number of time steps
         
   INTEGER(IntKi)                              :: ErrStat2                     ! Error level (local)
   CHARACTER(MaxMsgLen)                        :: ErrMsg2                      ! error message (local)
   CHARACTER(*), parameter                     :: RoutineName = 'TimeSeriesToSpectra'
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
  
      ! BJJ: consider putting the call to PSF in the reading part
   p%usr%nFreq  = PSF ( p%usr%NTimes/2, 9, .TRUE.)
   NumSteps     = p%usr%nFreq*2
         
   CALL AllocAry(p%usr%meanU,                   p%usr%NPoints,p%usr%nComp,'meanU',      ErrStat2,ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry(p%usr%meanDir,                 p%usr%NPoints,            'meanDir',    ErrStat2,ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry(p%usr%meanVAng,                p%usr%NPoints,            'meanVAng',   ErrStat2,ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry(p%usr%S,          p%usr%nFreq ,p%usr%NPoints,p%usr%nComp,'S',          ErrStat2,ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry(p%usr%f,          p%usr%nFreq ,                          'f',          ErrStat2,ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry(p%usr%phaseAngles,p%usr%nFreq ,p%usr%NPoints,p%usr%nComp,'phaseAngles',ErrStat2,ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   CALL AllocAry(work,             NumSteps,                              'work',       ErrStat2,ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   IF (ErrStat >= AbortErrLev) THEN
      CALL Cleanup()
      RETURN
   END IF
   
   !.......................................................................   
   ! calculate wind direction in radians (before rotating or removing mean):
   !.......................................................................   
   
   meanV = 0.0_ReKi
   meanW = 0.0_ReKi
      
   DO iPoint = 1, p%usr%NPoints     
      meanU = sum(p%usr%v(:,iPoint,1))/p%usr%NTimes
      
      if (p%usr%nComp > 1 ) then
         meanV = sum(p%usr%v(:,iPoint,2))/p%usr%NTimes
         
         if ( p%usr%nComp > 2 ) then      
            meanW = sum(p%usr%v(:,iPoint,3))/p%usr%NTimes
         end if
         
      end if
      
      p%usr%meanDir( iPoint) = atan2( meanV, meanU )
      p%usr%meanVAng(iPoint) = atan2( meanW, sqrt( meanU**2 + meanV**2 ) )
   END DO
   
   !.......................................................................   
   ! rotate inputs based on angles at reference point:
   !.......................................................................   
   uvw = 0.0_ReKi
   DO iPoint = 1, p%usr%NPoints 
      
      cosH = cos(p%usr%meanDir(  iPoint ) )
      sinH = sin(p%usr%meanDir(  iPoint ) )
      cosV = cos(p%usr%meanVAng( iPoint ) )
      sinV = sin(p%usr%meanVAng( iPoint ) )
   
      rotateMatrix(1,1) =  cosH*cosV
      rotateMatrix(2,1) = -sinH
      rotateMatrix(3,1) = -cosH*sinV
   
      rotateMatrix(1,2) =  sinH*cosV
      rotateMatrix(2,2) =  cosH
      rotateMatrix(3,2) = -sinH*sinV
   
      rotateMatrix(1,3) =       sinV
      rotateMatrix(2,3) = 0.0_ReKi
      rotateMatrix(3,3) =       cosV
                                 
      DO iTime = 1, p%usr%nTimes
         uvw(1:p%usr%nComp) = p%usr%v(iTime,iPoint,:)
         uvw = MATMUL( rotateMatrix, uvw )
         p%usr%v(iTime,iPoint,:) = uvw(1:p%usr%nComp)
      END DO      
      
   END DO
   
   ! now convert angles to degrees:
   p%usr%meanDir  = p%usr%meanDir  * R2D
   p%usr%meanVAng = p%usr%meanVAng * R2D
         
   
   !.......................................................................   
   ! calculate and remove the mean wind components:         
   !.......................................................................   
   DO iVec = 1,p%usr%nComp
      DO iPoint = 1, p%usr%NPoints
         p%usr%meanU(iPoint,iVec) = SUM( p%usr%v(:,iPoint,iVec), 1 ) / p%usr%NTimes
         p%usr%v(:,  iPoint,iVec) =      p%usr%v(:,iPoint,iVec)    -   p%usr%meanU(iPoint,iVec)
      END DO
   END DO
               
   
   !.......................................................................   
   ! compute forward fft to get real and imaginary parts      
   ! S = Re^2 + Im^2 
   ! PhaseAngle = acos( Re / S )
   !.......................................................................   
   
   
   CALL InitFFT( NumSteps, FFT_Data, NormalizeIn=.TRUE., ErrStat=ErrStat2 )
      CALL SetErrStat(ErrStat2, 'Error in InitFFT', ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

      
   ! Get the stationary-point time series.
   DO iVec=1,p%usr%nComp
      DO iPoint=1,p%usr%NPoints    

         work = p%usr%v(:,iPoint,iVec)
         
            ! perform forward FFT

         CALL ApplyFFT_f( work, FFT_Data, ErrStat2 )
            IF (ErrStat2 /= ErrID_None ) THEN
               CALL SetErrStat(ErrStat2, 'Error in ApplyFFT_f for point '//TRIM(Num2LStr(iPoint))//'.', ErrStat, ErrMsg, RoutineName )
               IF (ErrStat >= AbortErrLev) EXIT
            END IF
            
            
         DO iFreq = 1,p%usr%nFreq-1
            Indx = iFreq*2
            Re = work(Indx)
            Im = work(Indx+1)
            
            p%usr%S(          iFreq,iPoint,iVec) = Re**2 + Im**2
            
            IF ( p%usr%S( iFreq,iPoint,iVec) /= 0.0_ReKi ) THEN
               p%usr%PhaseAngles(iFreq,iPoint,iVec) = atan2( Im, Re ) ! this gives us the angles in range -pi to pi
               if ( p%usr%PhaseAngles(iFreq,iPoint,iVec) < 0.0_ReKi  ) then ! we want it in the range 0 to 2pi
                  p%usr%PhaseAngles(iFreq,iPoint,iVec) = TwoPi + p%usr%PhaseAngles(iFreq,iPoint,iVec)
               end if                       
            ELSE
               p%usr%PhaseAngles(iFreq,iPoint,iVec) = 0.0_ReKi
            END IF                     
         END DO
                           
         p%usr%S(          p%usr%nFreq,iPoint,iVec) = work(NumSteps)**2 !0.0_ReKi !work(NumSteps)**2  ! this frequency doesn't seem to get used in the code, so I'm going to set it to zero. <<< bjj: will be used for extrapolation      
         p%usr%PhaseAngles(p%usr%nFreq,iPoint,iVec) = 0.0_ReKi         

      ENDDO ! IPoint
   ENDDO ! IVec 
   
   ! calculate associated frequencies:
   p%usr%DelF = 1.0_DbKi / ( NumSteps * ( p%usr%t(2) - p%usr%t(1) ) ) !store frequency in double precision for cases of very long time series
   do iFreq=1,p%usr%nFreq
      p%usr%f(iFreq) = p%usr%DelF * iFreq
   end do
   
   p%usr%S = p%usr%S*real(2.0_DbKi/p%usr%DelF,ReKi) ! p%usr%S*2.0_ReKi/p%usr%f(1)  ! make this the single-sided velocity spectra we're using in the rest of the code
      

   CALL Cleanup()

   RETURN
CONTAINS
!...........................................
   SUBROUTINE Cleanup()

   CALL ExitFFT( FFT_Data, ErrStat2 )         
   CALL SetErrStat(ErrStat2, 'Error in ExitFFT', ErrStat, ErrMsg, RoutineName )
   
   IF ( ALLOCATED(work) ) DEALLOCATE(work)
   
   ! we don't need these anymore:
   IF ( ALLOCATED( p%usr%t             ) )  DEALLOCATE( p%usr%t             )
   IF ( ALLOCATED( p%usr%v             ) )  DEALLOCATE( p%usr%v             )
   
   END SUBROUTINE Cleanup   
            
END SUBROUTINE TimeSeriesToSpectra
!=======================================================================
SUBROUTINE TS_End(p, OtherSt_RandNum)


   TYPE(TurbSim_ParameterType),  INTENT(INOUT) :: p                 !< parameters 
   TYPE(RandNum_OtherStateType), INTENT(INOUT) :: OtherSt_RandNum   !< other states for random numbers (next seed, etc)

   
   IF (p%US > 0) THEN
      CLOSE( p%US )
      p%US = -1
   END IF
   
      
!bjj: todo: add more; make sure everything is deallocated here; make sure files are closed, too.   
         
   IF ( ALLOCATED( p%grid%Y            ) )  DEALLOCATE( p%grid%Y            )
   IF ( ALLOCATED( p%grid%Z            ) )  DEALLOCATE( p%grid%Z            )
   IF ( ALLOCATED( p%grid%GridPtIndx   ) )  DEALLOCATE( p%grid%GridPtIndx   )
   IF ( ALLOCATED( p%grid%TwrPtIndx    ) )  DEALLOCATE( p%grid%TwrPtIndx    )
   IF ( ALLOCATED( p%grid%Freq         ) )  DEALLOCATE( p%grid%Freq         )
   
   IF ( ALLOCATED( p%met%ZL_profile    ) )  DEALLOCATE( p%met%ZL_profile    )
   IF ( ALLOCATED( p%met%Ustar_profile ) )  DEALLOCATE( p%met%Ustar_profile )
   
   IF ( ALLOCATED( p%met%USR_Z         ) )  DEALLOCATE( p%met%USR_Z         )
   IF ( ALLOCATED( p%met%USR_U         ) )  DEALLOCATE( p%met%USR_U         )
   IF ( ALLOCATED( p%met%USR_WindDir   ) )  DEALLOCATE( p%met%USR_WindDir   )
   IF ( ALLOCATED( p%met%USR_Sigma     ) )  DEALLOCATE( p%met%USR_Sigma     )
   IF ( ALLOCATED( p%met%USR_L         ) )  DEALLOCATE( p%met%USR_L         )
                                                                                           
   IF ( ALLOCATED( p%usr%pointyi       ) )  DEALLOCATE( p%usr%pointyi       )
   IF ( ALLOCATED( p%usr%pointzi       ) )  DEALLOCATE( p%usr%pointzi       )
   IF ( ALLOCATED( p%usr%t             ) )  DEALLOCATE( p%usr%t             )
   IF ( ALLOCATED( p%usr%v             ) )  DEALLOCATE( p%usr%v             )

   IF ( ALLOCATED( p%usr%meanU         ) )  DEALLOCATE( p%usr%meanU         )
   IF ( ALLOCATED( p%usr%meanDir       ) )  DEALLOCATE( p%usr%meanDir       )
   IF ( ALLOCATED( p%usr%meanVAng      ) )  DEALLOCATE( p%usr%meanVAng      )
   
   IF ( ALLOCATED( p%usr%f             ) )  DEALLOCATE( p%usr%f             )
   IF ( ALLOCATED( p%usr%S             ) )  DEALLOCATE( p%usr%S             )
   IF ( ALLOCATED( p%usr%phaseAngles   ) )  DEALLOCATE( p%usr%phaseAngles   )
   
   IF ( ALLOCATED( p%RNG%RandSeedAry   ) )  DEALLOCATE( p%RNG%RandSeedAry   )
   
   
   IF (ALLOCATED(OtherSt_RandNum%nextSeed) ) DEALLOCATE(OtherSt_RandNum%nextSeed)  
         

END SUBROUTINE TS_END
!=======================================================================
END MODULE TSSubs
