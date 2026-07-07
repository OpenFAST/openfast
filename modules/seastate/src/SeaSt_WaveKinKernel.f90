!**********************************************************************************************************************************
! The SeaSt_WaveKinKernel module holds the per-point wave-kinematics generation math shared by full-domain initialization
! (WvKinBlockMod=0) and on-demand block population (WvKinBlockMod=1).
!
! Numerical note: the per-frequency generation seeds (wave numbers, intrinsic frequencies, MacCamy-Fuchs coefficients) are
! captured inside VariousWaves_Init (Waves.f90) from the exact values computed there, NOT recomputed here — calling the
! dispersion solvers from a different compilation unit changes gfortran's inlining/FMA-contraction decisions and drifts
! results by 1 ulp, breaking bit-for-bit equivalence with full-domain precompute. For the same reason this module carries
! its own private copies of the hyperbolic kinematics helpers (verbatim from Waves.f90) so that the per-column kernel and
! its helpers are compiled together in one unit, mirroring how Waves.f90 compiles them today.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2026  National Renewable Energy Laboratory
!
!    This file is part of SeaState.
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
MODULE SeaSt_WaveKinKernel

   USE Waves_Types
   USE SeaSt_WaveField_Types
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

      ! ..... Public Subroutines ...................................................................................................
   PUBLIC :: WaveKinKernel_CaptureGridSeeds       ! Capture the grid z levels and steady current profile into a block store

CONTAINS

!----------------------------------------------------------------------------------------------------------------------------------
!> Capture the per-z-level generation seeds the on-demand block-population kernel needs: the grid z coordinates and the
!! steady current profile (pure copies — no floating-point arithmetic). The per-frequency seeds (WaveNmbrArr, OmegaIArr,
!! MCFCArr) are captured separately inside VariousWaves_Init; see the numerical note in the module header.
!! Must be called after Waves_Init has completed.
SUBROUTINE WaveKinKernel_CaptureGridSeeds( InitInp, Store, ErrStat, ErrMsg )

      TYPE(Waves_InitInputType),      INTENT(IN   ) :: InitInp     !< Waves initialization input (still allocated at SeaState init time)
      TYPE(SeaSt_WaveBlockStoreType), INTENT(INOUT) :: Store       !< Block store receiving the generation seeds
      INTEGER(IntKi),                 INTENT(  OUT) :: ErrStat     !< Error status of the operation
      CHARACTER(*),                   INTENT(  OUT) :: ErrMsg      !< Error message if ErrStat /= ErrID_None

         ! Local Variables:
      INTEGER(IntKi)               :: K                        ! Grid z-level index
      INTEGER(IntKi)               :: iFlat                    ! Flat grid-point index of point (1,1,K)
      INTEGER(IntKi)               :: NX, NY, NZ               ! Grid extents
      INTEGER(IntKi)               :: ErrStatTmp               ! Temporary error status
      CHARACTER(*),  PARAMETER     :: RoutineName = 'WaveKinKernel_CaptureGridSeeds'

      ErrStat = ErrID_None
      ErrMsg  = ""

      ! Grid z levels and the steady current profile (uniform over x and y; one value per z level).
      ! The flat index of grid point (i,j,k) is i + (j-1)*NX + (k-1)*NX*NY (x fastest, z slowest),
      ! matching the grid copy loops of VariousWaves_Init; point (1,1,K) is at 1 + (K-1)*NX*NY.
      NX = InitInp%NGrid(1)
      NY = InitInp%NGrid(2)
      NZ = InitInp%NGrid(3)
      ALLOCATE ( Store%zGrid(NZ), Store%CurrVxi(NZ), Store%CurrVyi(NZ), STAT=ErrStatTmp )
      IF ( ErrStatTmp /= 0 ) THEN
         CALL SetErrStat(ErrID_Fatal,'Error allocating the per-level seed arrays.',ErrStat,ErrMsg,RoutineName)
         RETURN
      END IF

      Store%HasCurr = ALLOCATED(InitInp%CurrVxi)
      DO K = 1, NZ
         iFlat = 1 + (K-1)*NX*NY
         Store%zGrid(K) = InitInp%WaveKinGridzi(iFlat)
         IF ( Store%HasCurr ) THEN
            Store%CurrVxi(K) = InitInp%CurrVxi(iFlat)
            Store%CurrVyi(K) = InitInp%CurrVyi(iFlat)
         ELSE
            Store%CurrVxi(K) = 0.0_SiKi
            Store%CurrVyi(K) = 0.0_SiKi
         END IF
      END DO

END SUBROUTINE WaveKinKernel_CaptureGridSeeds

!----------------------------------------------------------------------------------------------------------------------------------
! The hyperbolic kinematics helpers below are verbatim private copies of the ones in Waves.f90 (see the module header for
! why they are duplicated rather than shared). They are used by the per-column generation kernel added in later commits.
!----------------------------------------------------------------------------------------------------------------------------------
      !=======================================================================
      FUNCTION COSHNumOvrCOSHDen ( k, h, z )


         ! This FUNCTION computes the shallow water hyperbolic numerator
         ! over denominator term in the wave kinematics expressions:
         !
         !                    COSH( k*( z + h ) )/COSH( k*h )
         !
         ! given the wave number, k, water depth, h, and elevation z, as
         ! inputs.

      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(SiKi)                   :: COSHNumOvrCOSHDen                               ! This function = COSH( k*( z + h ) )/COSH( k*h ) (-)
      REAL(ReKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
      REAL(SiKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
      REAL(SiKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)



         ! Compute the hyperbolic numerator over denominator:

      IF ( k*h  > 89.4_SiKi )  THEN   ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, COSH( k*( z + h ) )/COSH( k*h ) = EXP( k*z ) + EXP( -k*( z + 2*h ) ) for large k*h.  This equals the deep water formulation, EXP( k*z ), except near z = -h, because h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

         COSHNumOvrCOSHDen = EXP( k*z ) + EXP( -k*( z + 2.0_SiKi*REAL(h,SiKi) ) )

      ELSE                       ! 0 < k*h <= 89.4; use the shallow water formulation.

         COSHNumOvrCOSHDen =REAL( COSH( k*( z + REAL(h,SiKi) ) ),R8Ki)/COSH( k*REAL(h,SiKi) )

      END IF



      RETURN
      END FUNCTION COSHNumOvrCOSHDen
!=======================================================================
      FUNCTION COSHNumOvrSINHDen ( k, h, z )


         ! This FUNCTION computes the shallow water hyperbolic numerator
         ! over denominator term in the wave kinematics expressions:
         !
         !                    COSH( k*( z + h ) )/SINH( k*h )
         !
         ! given the wave number, k, water depth, h, and elevation z, as
         ! inputs.



      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(SiKi)                   :: COSHNumOvrSINHDen                               ! This function = COSH( k*( z + h ) )/SINH( k*h ) (-)
      REAL(ReKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
      REAL(SiKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
      REAL(SiKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)



         ! Compute the hyperbolic numerator over denominator:


      IF (   k  < EPSILON(0.0_SiKi)  )  THEN  ! When .TRUE., the shallow water formulation is ill-conditioned; thus, HUGE(k) is returned to approximate the known value of infinity.

         COSHNumOvrSINHDen = HUGE( k )

      ELSEIF ( k*REAL(h,SiKi)  > 89.4_SiKi )  THEN  ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, COSH( k*( z + h ) )/SINH( k*h ) = EXP( k*z ) + EXP( -k*( z + 2*h ) ) for large k*h.  This equals the deep water formulation, EXP( k*z ), except near z = -h, because h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

         COSHNumOvrSINHDen = EXP( k*z ) + EXP( -k*( z + 2*REAL(h,SiKi) ) )

      ELSE                          ! 0 < k*h <= 89.4; use the shallow water formulation.

         COSHNumOvrSINHDen = COSH( k*( z + REAL(h,SiKi) ) )/SINH( k*REAL(h,SiKi) )

      END IF



      RETURN
      END FUNCTION COSHNumOvrSINHDen

      !=======================================================================
      FUNCTION SINHNumOvrSINHDen ( k, h, z )


         ! This FUNCTION computes the shallow water hyperbolic numerator
         ! over denominator term in the wave kinematics expressions:
         !
         !                    SINH( k*( z + h ) )/SINH( k*h )
         !
         ! given the wave number, k, water depth, h, and elevation z, as
         ! inputs.


      IMPLICIT                        NONE


         ! Passed Variables:

      REAL(SiKi)                   :: SINHNumOvrSINHDen                               ! This function = SINH( k*( z + h ) )/SINH( k*h ) (-)
      REAL(ReKi), INTENT(IN )      :: h                                               ! Water depth ( h      >  0 ) (meters)
      REAL(SiKi), INTENT(IN )      :: k                                               ! Wave number ( k      >= 0 ) (1/m)
      REAL(SiKi), INTENT(IN )      :: z                                               ! Elevation   (-h <= z <= 0 ) (meters)



         ! Compute the hyperbolic numerator over denominator:

      IF (     k   == 0.0_SiKi  )  THEN  ! When .TRUE., the shallow water formulation is ill-conditioned; thus, the known value of unity is returned.

         SINHNumOvrSINHDen = 1.0

      ELSEIF ( k*REAL(h,SiKi) >  89.4_SiKi )  THEN  ! When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, SINH( k*( z + h ) )/SINH( k*h ) = EXP( k*z ) - EXP( -k*( z + 2*h ) ) for large k*h.  This equals the deep water formulation, EXP( k*z ), except near z = -h, because h > 14.23*wavelength (since k = 2*Pi/wavelength) in this case.

         SINHNumOvrSINHDen = EXP( k*z ) - EXP( -k*( z + 2.0_SiKi*h ) )

      ELSE                          ! 0 < k*h <= 89.4; use the shallow water formulation.

         SINHNumOvrSINHDen = SINH( k*( z + REAL(h,SiKi) ) )/SINH( k*REAL(h,SiKi) )

      END IF



      RETURN
      END FUNCTION SINHNumOvrSINHDen

END MODULE SeaSt_WaveKinKernel
!**********************************************************************************************************************************
