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
MODULE TS_VelocitySpectra

   USE  TurbSim_Types
   
   IMPLICIT NONE

CONTAINS


!=======================================================================
!> This subroutine defines the Kaimal PSD model as specified by IEC 61400-1, 2nd Ed. & 3rd Ed.
!! the use of this subroutine requires that all variables have the units of meters and seconds.
SUBROUTINE Spec_IECKAI ( UHub, SigmaIEC, L_K, Freq, NumFreq, Spec )


   IMPLICIT                NONE

         ! Passed variables
   INTEGER(IntKi),             INTENT(IN   )  :: NumFreq                   !<  Input: Number of frequencies      
   REAL(ReKi),                 INTENT(IN   )  :: SigmaIEC (3)              !<  Input: sigma for 3 wind components specified by the IEC
   REAL(ReKi),                 INTENT(IN   )  :: L_k      (3)              !<  Input: L_k (integral scale parameter) for 3 wind components specified by the IEC
   REAL(ReKi),                 INTENT(IN   )  :: UHub                      !<  Input: mean wind speed at hub height
   REAL(ReKi),                 INTENT(IN   )  :: Freq     (NumFreq)        !<  Input: frequency array
   REAL(ReKi),                 INTENT(INOUT)  :: Spec     (NumFreq,3)      !<  Output: target spectrum

      ! Internal variables

   REAL(ReKi),PARAMETER  :: Exp1    = 5.0/3.0

   REAL(ReKi)            :: L_over_U      (3)
   REAL(ReKi)            :: SigmaLU (3)

   INTEGER               :: I
   INTEGER               :: IVec



   ! Create the spectrum.
L_over_U = L_k / UHub
SigmaLU  = 4.0 * SigmaIEC**2 * L_over_U     ! array operations

DO IVec = 1,3

   L_over_U(IVec) = 6.0*L_over_U(IVec)

   DO I = 1,NumFreq
      Spec(I,IVec) = SigmaLU(IVec) / ( 1.0 + L_over_U(IVec)*Freq(I) )**Exp1
   ENDDO !I

ENDDO !IVec


RETURN
END SUBROUTINE Spec_IECKAI
!=======================================================================
!> This subroutine defines the von Karman PSD model as specified by IEC 61400-1 (2nd Ed).
!! The use of this subroutine requires that all variables have the units of meters and seconds.
SUBROUTINE Spec_IECVKM ( UHub, SigmaIEC_u, IntegralScale, Freq, NumFreq, Spec )

   IMPLICIT                NONE

         ! Passed variables

   INTEGER(IntKi),               INTENT(IN   )  :: NumFreq                   !<  Input: Number of frequencies      
   REAL(ReKi),                   INTENT(IN   )  :: SigmaIEC_u                !<  Input: target standard deviation for u component
   REAL(ReKi),                   INTENT(IN   )  :: IntegralScale (3)         !<  Input: integral scale parameter, L (isotropic, so we only care about the 1st one)
   REAL(ReKi),                   INTENT(IN   )  :: UHub                      !<  Input: mean wind speed at hub height
   REAL(ReKi),                   INTENT(IN   )  :: Freq     (NumFreq)        !<  Input: frequency array
   REAL(ReKi),                   INTENT(  OUT)  :: Spec     (NumFreq,3)      !<  Output: target spectrum


         ! Internal variables

   REAL(ReKi),PARAMETER  :: Exp1 =  5.0/6.0
   REAL(ReKi),PARAMETER  :: Exp2 = 11.0/6.0
   REAL(ReKi)            :: FLU2
   REAL(ReKi)            :: L1_U
   REAL(ReKi)            :: SigmaL1_U
   REAL(ReKi)            :: Tmp

   INTEGER               :: I


   ! Set up scaling values.


   ! Define u-component integral scale.

L1_U   = IntegralScale(1)/UHub
SigmaL1_U = 2.0*SigmaIEC_u*SigmaIEC_u*L1_U

DO I=1,NumFreq

   FLU2      = ( Freq(I)*L1_U )**2
   Tmp       = 1.0 + 71.0*FLU2

   Spec(I,1) = 2.0*SigmaL1_U/Tmp**Exp1
   Spec(I,2) = SigmaL1_U*( 1.0 + 189.0*FLU2 )/Tmp**Exp2
   Spec(I,3) = Spec(I,2)

ENDDO ! I

RETURN
END SUBROUTINE Spec_IECVKM
!=======================================================================
!> This subroutine defines the API-BULLET-IN recommended extreme wind spectrum
!! The use of this subroutine requires that all variables have the units of meters and seconds.
!! See A.7.4 (Page 41) of API 2MET/ISO 19901-1:2005(E).
!! See https://rules.dnvgl.com/docs/pdf/DNV/codes/docs/2010-10/RP-C205.pdf (page 20 of 124), describing the
!! Froya model spectral density proposed by Andersen and Lovseth (1992, 2006) for wind over water.
SUBROUTINE Spec_API ( p, Ht, Spec )

   ! NOTE: This routine uses the Kaimal model to create the spectrum for all three components
   !       and then overwrites the u-component spectrum with the API model.


IMPLICIT                NONE

      ! Passed variables
   TYPE(TurbSim_ParameterType) , INTENT(IN   ) :: p                       !< Input: turbsim parameters
   REAL(ReKi),                   INTENT(IN)    :: Ht                      !<  Input: Height (Should be HubHt), value ignored !bjj: is this true????
   REAL(ReKi),                   INTENT(INOUT) :: Spec   (:,:)            !<  Output: target spectrum

!REAL(ReKi),INTENT(IN) :: URef ! Added by YG
!REAL(ReKi),INTENT(IN)           :: RefHt                       ! Reference height

      ! Internal variables

REAL(ReKi),PARAMETER  :: N    =  0.468 
REAL(ReKi),PARAMETER  :: Exp5 = 5.0/( 3.0*N )
!mlb REAL(ReKi),PARAMETER  :: Exp5 = 11.0/6.0
REAL(ReKi),PARAMETER  :: Ref_Ht = 10.0
REAL(ReKi),PARAMETER  :: Ref_WS = 10.0
REAL(ReKi)            :: Scale1
REAL(ReKi)            :: Scale2
REAL(ReKi)            :: Temp
!mlb REAL(ReKi)            :: FLU2
!mlb REAL(ReKi)            :: L1_U
!mlb REAL :: X0=10.0                                              ! Added by Y. Guo for calculating UHr_10
!mlb REAL :: X
INTEGER               :: I
!mlb REAL                  :: UHr_10

   ! Set up scaling values.

   ! calculate the spectra for the v and w components using IECKAI model
   ! because API doesn't specify a spectra for those components
CALL Spec_IECKAI ( p%UHub, p%IEC%SigmaIEC, p%IEC%IntegralScale, p%grid%Freq, p%grid%NumFreq, Spec )

   ! Define u-component integral scale.
!CALL WrScr ('Calling Froya/API wind spectrum.............')
!mlb L1_U   = 3.5*Lambda/p%UHub
!mlb SigmaL1_U = 2.0*SigmaIEC(1)*p%UHub*L1_U

!mlb CALL ROOT_SEARCHING(X0,X,p%UHub,Ht,Ht)
!mlb UHr_10=X;

   ! Compute some parameters that are independent of frequency.

Scale1 = 172.0*( Ht/Ref_Ht )**(2.0/3.0) * ( p%met%URef/Ref_WS )**(-0.75)
Scale2 = 320.0*( p%met%URef/Ref_WS )**2 * ( Ht/Ref_Ht )**0.45

DO I=1,p%grid%NumFreq


!mlb    Tmp1      = 172.0*p%grid%Freq(I)*(Ht/10.0)**Exp2*(UHr_10/10.0)**Exp3
!mlb    Tmp2       = (1.0+Tmp1**Exp1)**(5.0/3.0/Exp1)
!mlb 
!mlb    Spec(I,1) = 320.0*(UHr_10/10.0)**2*(Ht/10.0)**Exp4/Tmp2

   Temp      = Scale1*p%grid%Freq(I)
   Spec(I,1) = Scale2/( 1.0 + Temp**N )**Exp5

ENDDO ! I

!CALL WrScr ('Froya/API wind spectrum generated')

RETURN
END SUBROUTINE Spec_API
!=======================================================================
!> This subroutine defines the 3-D turbulence spectrum that can be expected over terrain
!! and heights similiar to the LLLJP project as developed by Neil Kelley & Bonnie Jonkman at NREL.
!! The use of this subroutine requires that variables have the units of meters and seconds.
SUBROUTINE Spec_GPLLJ ( p, Ht, Ucmp, ZL_tmp, UStar_tmp, Spec )

IMPLICIT                NONE

   ! Passed variables

   TYPE(TurbSim_ParameterType) , INTENT(IN   ) :: p                       !< Input: turbsim parameters
   REAL(ReKi),                   INTENT(IN)    :: Ht                      !< Height (local)
   REAL(ReKi),                   INTENT(IN)    :: Ucmp                    !< Longitudinal Velocity (local)
   REAL(ReKi),                   INTENT(IN)    :: Ustar_tmp               !< Local ustar
   REAL(ReKi),                   INTENT(IN)    :: ZL_tmp                  !< Local z/l
   REAL(ReKi),                   INTENT(  OUT) :: Spec   (:,:)            !<  Output: target spectrum


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
      DO I = 1,p%grid%NumFreq
         tmpX  = p%grid%Freq(I)*tmpF             ! reduced frequency divided by q (q = phiM here)
         X_l   = tmpX/fr_il(ic)
         X_h   = tmpX/fr_ih(ic)

         ps_l  = (Pr_il(ic)*scales(1,ic)*X_l*tmpPhi) / (1.0 + scales(2,ic)*X_l**Exp53);
         ps_h  = (Pr_ih(ic)*scales(1,ic)*X_h*tmpPhi) / (1.0 + scales(2,ic)*X_h**Exp53);

         Spec(I,IC) = (ps_l + ps_h)*uStar2/p%grid%Freq(I)
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

   tmpZIL = ( ABS(p%met%ZI / p%met%L) )**Exp23
   HtZI   = Ht / p%met%ZI

   tmpZIU = p%met%ZI / Ucmp
   tmpZU  = Ht / Ucmp
   HtZI2  = (1.0 - HtZI)**2
   UDen   = 1.0 + 15.0*HtZI
   VDen   = 1.0 +  2.8*HtZI

   DO I=1,p%grid%NumFreq
      fi      = p%grid%Freq(I)*tmpZIU
      tmpF    = p%grid%Freq(I)*tmpZU                ! reduced frequency
      Ustar2F = uStar2/p%grid%Freq(I)               ! Normalizing term

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
END SUBROUTINE Spec_GPLLJ
!=======================================================================
!> This subroutine defines the 3-D turbulence spectrum that can be expected
!! over terrain and heights similiar to the NWTC LIST project as developed
!! by Neil Kelley & Bonnie Jonkman at NREL. The use of this subroutine
!! requires that variables have the units of meters and seconds.
SUBROUTINE Spec_NWTCUP ( p, Ht, Ucmp, Spec )


IMPLICIT                NONE

   ! Passed variables

   TYPE(TurbSim_ParameterType) , INTENT(IN   ) :: p                       !< Input: turbsim parameters
   REAL(ReKi),                   INTENT(IN   ) :: Ht                      !< Height (local)
   REAL(ReKi),                   INTENT(IN   ) :: Ucmp                    !< Longitudinal Velocity (local)
   REAL(ReKi),                   INTENT(  OUT) :: Spec   (:,:)            !< Output: target spectrum



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

uStar2 = p%met%Ustar * p%met%Ustar

IF (p%met%zL >= 0) THEN

   zl_tmp = max( min(p%met%zl, 3.5_ReKi ), 0.005_ReKi )

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
      DO I = 1,p%grid%NumFreq
         tmpX  = p%grid%Freq(I)*tmpF             ! reduced frequency divided by q (q = phiM here)
         X_l   = tmpX/fr_il(ic)
         X_h   = tmpX/fr_ih(ic)

         ps_l  = (Pr_il(ic)*scales(1,ic)*X_l*tmpPhi) / (1.0 + scales(2,ic)*X_l**Exp53);
         ps_h  = (Pr_ih(ic)*scales(1,ic)*X_h*tmpPhi) / (1.0 + scales(2,ic)*X_h**Exp53);

         Spec(I,IC) = (ps_l + ps_h)*uStar2/p%grid%Freq(I)
      ENDDO
   ENDDO


ELSE
      ! Calculate UNSTABLE spectral estimates

   zl_tmp    = abs( min( max( p%met%zl  ,REAL(-0.5,ReKi) ),REAL( -0.025,ReKi) ) )
   ustar_tmp =      max( min(p%met%ustar,REAL( 1.4,ReKi) ),REAL(  0.2  ,ReKi) )

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

   tmpZIL = (-p%met%ZI / p%met%L)**Exp23
   HtZI   = Ht / p%met%ZI

   tmpZIU = p%met%ZI / Ucmp
   tmpZU  = Ht / Ucmp
   HtZI2  = (1.0 - HtZI)**2
   UDen   = 1.0 + 15.0*HtZI
   VDen   = 1.0 +  2.8*HtZI

   DO I=1,p%grid%NumFreq
      fi      = p%grid%Freq(I)*tmpZIU

      tmpF    = p%grid%Freq(I)*tmpZU                ! reduced frequency
      Ustar2F = uStar2/p%grid%Freq(I)               ! Normalizing term

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

END SUBROUTINE Spec_NWTCUP
!=======================================================================
!> This routine gets velocity spectra for each of 3 wind components (u,v,w)
!! by 2-D interpolation.
SUBROUTINE Spec_TimeSer ( p, Ht, Ucmp, LastIndex, Spec )


         ! Passed variables
   TYPE(TurbSim_ParameterType), INTENT(IN   ) :: p                       !< Input: turbsim parameters
   REAL(ReKi),                  INTENT(IN   ) :: Ht                      !< Input: height for which spectra are requested
   REAL(ReKi),                  INTENT(IN   ) :: Ucmp                    !< Input: wind speed for which spectra are requested (used for missing components)
   INTEGER(IntKi),              INTENT(INOUT) :: LastIndex(2)            !< Index for the last (Freq, Ht) used
   REAL(ReKi),                  INTENT(  OUT) :: Spec   (:,:)            !< Output: target spectrum (Frequency, component)
         ! Local variables
   REAL(ReKi)                                 :: InCoord(2)              ! Arranged as (Freq, Ht)
   INTEGER(IntKi)                             :: i                       ! loop counters
   
               

!bjj: fix me!!! (make use of nComp and height )                             

      ! initialize Spec with extrapolated values on non-specified components or frequencies
      ! i.e., fill the gaps where wind component or frequencies exceed what was specified in the time-series data with some numerical model
      ! (or use zeros for known spectral values that will get overwritten later)
   CALL Spec_TimeSer_Extrap ( p, Ht, Ucmp, Spec )


   InCoord(2) = Ht
   
      ! overwrite Spec at the frequencies and wind components by interpolating from known points
   DO I=1,p%usr%nFreq !p%grid%NumFreq ! note that this assumes TMax = AnalysisTime (i.e., we have the same delta frequencies)
      
      InCoord(1) = p%grid%Freq(i)      
      CALL UserSpec_Interp2D( InCoord, p%usr, LastIndex, Spec(I,:) )  ! sets only Spec(1:p%usr%nFreq, 1:p%usr%nComp) values
      
   ENDDO ! I
      
   RETURN      
      
                           
END SUBROUTINE Spec_TimeSer
!=======================================================================
!< This routine adds high-frequency content to user-supplied data, 
!! using the model specified in p%usr%TurbModel_ID
SUBROUTINE Spec_TimeSer_Extrap ( p, Ht, Ucmp, Spec )


         ! Passed variables
   TYPE(TurbSim_ParameterType), INTENT(IN   ) :: p                       !< Input: turbsim parameters
   REAL(ReKi),                  INTENT(IN   ) :: Ht                      !< Input: height for which spectra are requested
   REAL(ReKi),                  INTENT(IN   ) :: Ucmp                    !< Input: wind speed for which spectra are requested (used for missing components)
   REAL(ReKi),                  INTENT(  OUT) :: Spec   (:,:)            !< Output: target spectrum (Frequency, component)


   IF ( p%usr%nComp < 3 .OR. p%usr%nFreq < p%grid%NumFreq ) THEN
         
      SELECT CASE ( p%usr%TurbModel_ID )
         CASE ( SpecModel_IECKAI ) ! IECKAI has uniform spectra (does not vary with height or velocity)
            CALL Spec_IECKAI  ( p%UHub, p%IEC%SigmaIEC, p%IEC%IntegralScale, p%grid%Freq, p%grid%NumFreq, Spec )
            
         CASE ( SpecModel_IECVKM )  ! IECVKM has uniform spectra (does not vary with height or velocity)
            CALL Spec_IECVKM  ( p%UHub, p%IEC%SigmaIEC(1), p%IEC%IntegralScale, p%grid%Freq, p%grid%NumFreq, Spec )
                                    
         CASE ( SpecModel_API )
            CALL Spec_API ( p, Ht, Spec )
                              
         CASE ( SpecModel_SMOOTH )
               CALL Spec_SMOOTH   ( p, Ht, Ucmp, Spec )
                                                               
         CASE DEFAULT
            Spec = 0.0_ReKi ! whole matrix is zero
            
         END SELECT          
                                 
   ELSE
      Spec = 0.0_ReKi ! whole matrix is zero
   END IF

END SUBROUTINE Spec_TimeSer_Extrap
!=======================================================================
!< This routine linearly interpolates the p%usr spectral data. It is
!! set for a 2-d interpolation on frequency and height of the input point.
!! p%usr%f and p%usr%pointzi must be in increasing order. Each dimension
!! may contain only 1 value.
SUBROUTINE UserSpec_Interp2D( InCoord, p_usr, LastIndex, OutSpec )

      ! I/O variables

   REAL(ReKi),                     INTENT(IN   ) :: InCoord(2)                                   !< Arranged as (Freq, Ht)
   TYPE(UserTSSpec_ParameterType), INTENT(IN   ) :: p_usr                                        !<
   INTEGER(IntKi),                 INTENT(INOUT) :: LastIndex(2)                                 !< Index for the last (Freq, Ht) used
   REAL(ReKi),                     INTENT(INOUT) :: OutSpec(3)                                   !< The interpolated resulting PSD from each component of p%usr%S(:,:,1-3)


      ! Local variables
   
   INTEGER(IntKi)                                :: I                                            ! loop counter                                 
   INTEGER(IntKi)                                :: Indx_Lo(2)                                   ! index associated with lower bound of dimension 1,2 where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))
   INTEGER(IntKi)                                :: Indx_Hi(2)                                   ! index associated with upper bound of dimension 1,2 where val(Indx_lo(i)) <= InCoord(i) <= val(Indx_hi(i))
   REAL(ReKi)                                    :: Pos_Lo(2)                                    ! coordinate value with lower bound of dimension 1,2 
   REAL(ReKi)                                    :: Pos_Hi(2)                                    ! coordinate value with upper bound of dimension 1,2 
                                                 
   REAL(ReKi)                                    :: isopc(2)                                     ! isoparametric coordinates 
                                                                                                 
   REAL(ReKi)                                    :: N(4)                                         ! size 2^n
   REAL(ReKi)                                    :: u(4)                                         ! size 2^n
   
   
         
      ! find the indices into the arrays representing coordinates of each dimension:
      !  (by using LocateStp, we do not require equally spaced frequencies or points)
               
   CALL LocateStp( InCoord(1), p_usr%f,       LastIndex(1), p_usr%nFreq   )
   CALL LocateStp( InCoord(2), p_usr%pointzi, LastIndex(2), p_usr%nPoints )
   
   Indx_Lo = LastIndex  ! at this point, 0 <= Indx_Lo(i) <= n(i) for all i
   
   
   ! Frequency (indx 1)
   IF (Indx_Lo(1) == 0) THEN
      Indx_Lo(1) = 1
   ELSEIF (Indx_Lo(1) == p_usr%nFreq ) THEN
      Indx_Lo(1) = max( p_usr%nFreq - 1, 1 )                ! make sure it's a valid index
   END IF     
   Indx_Hi(1) = min( Indx_Lo(1) + 1 , p_usr%nFreq )         ! make sure it's a valid index

   ! Height (indx 2)
   IF (Indx_Lo(2) == 0) THEN
      Indx_Lo(2) = 1
   ELSEIF (Indx_Lo(2) == p_usr%nPoints ) THEN
      Indx_Lo(2) = max( p_usr%nPoints - 1, 1 )              ! make sure it's a valid index
   END IF     
   Indx_Hi(2) = min( Indx_Lo(2) + 1 , p_usr%nPoints )       ! make sure it's a valid index
      
         
      ! calculate the bounding box; the positions of all dimensions:
      
   pos_Lo(1) = p_usr%f( Indx_Lo(1) )
   pos_Hi(1) = p_usr%f( Indx_Hi(1) )
      
   pos_Lo(2) = p_usr%pointzi( Indx_Lo(2) )  ! note that this assumes z are in increasing order
   pos_Hi(2) = p_usr%pointzi (Indx_Hi(2) )
               
   
      ! 2-D linear interpolation:
      
   CALL IsoparametricCoords( InCoord, pos_Lo, pos_Hi, isopc )      ! Calculate iospc
   
   N(1)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi - isopc(2) )
   N(2)  = ( 1.0_ReKi + isopc(1) )*( 1.0_ReKi + isopc(2) )
   N(3)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi + isopc(2) )
   N(4)  = ( 1.0_ReKi - isopc(1) )*( 1.0_ReKi - isopc(2) )
   N     = N / REAL( SIZE(N), ReKi )  ! normalize
            
      
   do i = 1,p_usr%nComp
      u(1)  = p_usr%S( Indx_Hi(1), Indx_Lo(2), i )
      u(2)  = p_usr%S( Indx_Hi(1), Indx_Hi(2), i )
      u(3)  = p_usr%S( Indx_Lo(1), Indx_Hi(2), i )
      u(4)  = p_usr%S( Indx_Lo(1), Indx_Lo(2), i )
            
      OutSpec(i) = SUM ( N * u )                   
   end do
   
         
END SUBROUTINE UserSpec_Interp2D
!=======================================================================
!> This routine linearly interpolates data from an input file that 
!! specifies the velocity spectra for each of 3 wind components (u,v,w)
SUBROUTINE Spec_UserSpec ( p, Spec )

   IMPLICIT                NONE

         ! Passed variables
   type(TurbSim_ParameterType) , INTENT(IN   ) :: p                       !< Input: turbsim parameters
   REAL(ReKi),                   INTENT(  OUT) :: Spec   (:,:)            !< Output: target spectrum

      ! Internal variables

   REAL(ReKi)            :: Tmp


   INTEGER               :: I
   INTEGER               :: Indx
   INTEGER               :: J
   INTEGER,PARAMETER     :: iPoint = 1

      ! --------- Interpolate to the desired frequencies ---------------

   Indx = 1;

   DO I=1,p%grid%NumFreq

      IF ( p%grid%Freq(I) <= p%usr%f(1) ) THEN
         Spec(I,:) = p%usr%S(1,iPoint,:)
      ELSEIF ( p%grid%Freq(I) >= p%usr%f(p%usr%nFreq) ) THEN
         Spec(I,:) = p%usr%S(p%usr%nFreq,iPoint,:)
      ELSE

            ! Find the two points between which the frequency lies

         DO J=(Indx+1),p%usr%nFreq
            IF ( p%grid%Freq(I) <= p%usr%f(J) ) THEN
               Indx = J-1

                  ! Let's just do a linear interpolation for now

               Tmp  = (p%grid%Freq(I) - p%usr%f(Indx)) / ( p%usr%f(Indx) - p%usr%f(J) )

               Spec(I,:) = Tmp * ( p%usr%S(Indx,iPoint,:) - p%usr%S(J,iPoint,:) ) + p%usr%S(Indx,iPoint,:)

               EXIT
            ENDIF
         ENDDO ! J

      ENDIF

   ENDDO ! I

   RETURN


END SUBROUTINE Spec_UserSpec
!=======================================================================
!> This subroutine defines the 3-D turbulence spectrum that can be expected over flat,
!! homogeneous terrain as developed by RISO authors Hojstrup, Olesen, and Larsen.
!! The use of this subroutine requires that variables have the units of meters and seconds.
SUBROUTINE Spec_SMOOTH ( p, Ht, Ucmp, Spec )


IMPLICIT                NONE

   ! Passed variables

   TYPE(TurbSim_ParameterType) , INTENT(IN   ) :: p                       !< Input: turbsim parameters
   REAL(ReKi),                   INTENT(IN   ) :: Ht                      !< Height
   REAL(ReKi),                   INTENT(IN   ) :: Ucmp                    !< Longitudinal Velocity
   REAL(ReKi),                   INTENT(  OUT) :: Spec     (:,:)          !< output: target spectra

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

uStar2 = p%met%Ustar * p%met%Ustar

IF (p%met%zL >= 0) THEN

   ! Calculate NEUTRAL/STABLE spectral estimates

   phiM   =  1.0 + 4.7*(p%met%zL)            ! = q
   phiE   = (1.0 + 2.5*(p%met%zL)**0.6)**Exp3

   tmpPhi = uStar2 * ( (phiE / phiM)**Exp2 )
   tmpF   = Ht / (Ucmp * phiM)

   DO I = 1,p%grid%NumFreq
      tmpX  = p%grid%Freq(I)*tmpF             ! reduced frequency divided by q (q = phiM here)
      tmpXX = tmpX**Exp1
      tmpN  = tmpPhi / p%grid%Freq(I) * tmpX  ! normalization factor used to obtain power spectrum components

      Spec(I,1) = tmpN * (79.0) / (1.0 + 263.0*tmpXX)
      Spec(I,2) = tmpN * (13.0) / (1.0 +  32.0*tmpXX)
      Spec(I,3) = tmpN * ( 3.5) / (1.0 +   8.6*tmpXX)
   ENDDO

ELSE
   ! Calculate UNSTABLE spectral estimates
   tmpZIL = (- p%met%ZI / p%met%L)**Exp2
   HtZI   = Ht / p%met%ZI

   HtZI2  = (1.0 - HtZI)**2
   tmpZU  = Ht / Ucmp
   tmpZIU = p%met%ZI / Ucmp
   UDen   = 1.0 + 15.0*HtZI
   VDen   = 1.0 +  2.8*HtZI

   DO I = 1,p%grid%NumFreq

      Fi   = p%grid%Freq(I)*tmpZIU
      tmpF = p%grid%Freq(I)*tmpZU                ! reduced frequency

      Ustar2F = uStar2/p%grid%Freq(I)

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
END SUBROUTINE Spec_SMOOTH
!=======================================================================
!> This subroutine defines the 3-D turbulence expected in a tidal channel (HYDROTURBSIM specific).
!! It is similar to the 'smooth' spectral model (RISO; Hojstrup, Olesen and Larsen) for wind,
!! but is scaled by the TKE (SigmaU**2), and du/dz rather than Ustar and u/z.
!! The fit is based on data from Puget Sound, estimated by L. Kilcher.
!! The use of this subroutine requires that variables have the units of meters and seconds.
!! Note that this model does not require height.
SUBROUTINE Spec_TIDAL ( p, Ht, Shr_DuDz, Spec, SpecModel )

IMPLICIT                NONE

   ! Passed variables

   TYPE(TurbSim_ParameterType) , INTENT(IN   ) :: p                       !< Input: turbsim parameters
   REAL(ReKi),                   INTENT(IN   ) :: Ht                      !< Height (dz)
   REAL(ReKi),                   INTENT(IN   ) :: Shr_DuDz                !< Shear (du/dz)
   INTEGER(IntKi),               INTENT(IN   ) :: SpecModel               !< SpecModel (SpecModel_TIDAL .OR. SpecModel_RIVER)
   REAL(ReKi),                   INTENT(  OUT) :: Spec     (:,:)          !< output: target spectra

   ! Internal variables

   REAL(ReKi), PARAMETER :: Exp1  = 5.0 / 3.0
   REAL(ReKi)            :: Sigma_U2               ! Standard Deviation of U velocity, squared.
   REAL(ReKi)            :: Sigma_V2               ! Standard Deviation of V velocity, squared.
   REAL(ReKi)            :: Sigma_W2               ! Standard Deviation of W velocity, squared.

   REAL(ReKi)            :: tmpX                   ! Temporary variable for calculation of Spec
   REAL(ReKi)            :: tmpvec(3)              ! Temporary vector for calculation of Spec
   REAL(ReKi)            :: tmpa (3)               ! Spectra coefficients
   REAL(ReKi)            :: tmpb (3)               ! Spectra coefficients
   INTEGER               :: I                      ! DO LOOP counter



!print *, Ustar
!Sigma_U2=(TurbIntH20*U(IZ))**2 ! A fixed value of the turbulence intensity.  Do we want to implement this?
Sigma_U2=4.5*p%met%Ustar*p%met%Ustar*EXP(-2*Ht/p%met%RefHt)
Sigma_V2=0.5*Sigma_U2
Sigma_W2=0.2*Sigma_U2


SELECT CASE ( SpecModel )
   CASE ( SpecModel_TIDAL )
      tmpa = (/ 0.193, 0.053 , 0.0362 /)*TwoPi ! These coefficients were calculated using Shr_DuDz in units of 'radians', so we multiply these coefficients by 2*pi.
      tmpb = (/ 0.201, 0.0234, 0.0124 /)*(TwoPi**Exp1)
   CASE ( SpecModel_RIVER )
      ! THESE ARE NOT VERIFIED YET!!!, therefore they are undocumented.
      tmpa = (/ 0.081, 0.056 , 0.026 /)*TwoPi
      tmpb = (/ 0.16, 0.025, 0.020 /)*(TwoPi**Exp1)
END SELECT

tmpvec = tmpa*(/Sigma_U2, Sigma_V2, Sigma_W2/)/Shr_DuDz

DO I = 1,p%grid%NumFreq
   tmpX  = (p%grid%Freq(I)/Shr_DuDz)**Exp1
   Spec(I,1) = tmpvec(1) / (1.0 + tmpb(1)*tmpX)
   Spec(I,2) = tmpvec(2) / (1.0 + tmpb(2)*tmpX)
   Spec(I,3) = tmpvec(3) / (1.0 + tmpb(3)*tmpX)
ENDDO

RETURN
END SUBROUTINE Spec_TIDAL
!=======================================================================
!> This routine is just a test function to see if we get the requested
!! spectra from the TurbSim code.
SUBROUTINE Spec_Test ( Spec, Freq )

IMPLICIT                NONE

      ! Passed variables
   REAL(ReKi),                   intent(  out) :: Spec   (:,:)            !< Output: target spectrum
   REAL(ReKi),                   INTENT(IN   ) :: Freq(:)


INTEGER               :: I
INTEGER               :: IVec


   ! Create the spectrum.

DO IVec = 1,3

   DO I = 1,SIZE(Spec,1)
      Spec(I,IVec) = 0.0
   ENDDO !I
   !I = INT( NumFreq/2 )
   I = INT( 100 )
   Spec( I, IVec ) = 1/Freq(1)

   call WrScr( 'Test Spectra: sine wave with frequency '//trim(num2lstr(Freq(I)))//' Hz.' )

ENDDO !IVec


RETURN
END SUBROUTINE Spec_Test
!=======================================================================
!> This subroutine defines the von Karman PSD model.
!! The use of this subroutine requires that all variables have the units of meters and seconds.
SUBROUTINE Spec_vonKrmn ( p, Ht, Ucmp, Spec )


IMPLICIT                NONE

      ! Passed variables

   TYPE(TurbSim_ParameterType) , INTENT(IN   ) :: p                       !< Input: turbsim parameters
   REAL(ReKi),                   INTENT(IN   ) :: Ht                      !< local height
   REAL(ReKi),                   INTENT(IN   ) :: Ucmp                    !< local wind speed
   REAL(ReKi),                   INTENT(  OUT) :: Spec     (:,:)          !< Target spectra  

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
IF ( ALLOCATED( p%met%USR_L ) ) THEN
   IF ( Ht <= p%met%USR_Z(1) ) THEN
      Lvk = p%met%USR_L(1)   ! Extrapolation: nearest neighbor for heights below minimum height specified
   ELSEIF ( Ht >= p%met%USR_Z(p%met%NumUSRz) ) THEN
      Lvk = p%met%USR_L(p%met%NumUSRz)  ! Extrapolation: nearest neighbor for heights above maximum height specified
   ELSE !Interpolation: linear between user-defined height/integral scale curves
      DO I=2,p%met%NumUSRz
         IF ( Ht <= p%met%USR_Z(I) ) THEN
            Lvk = (Ht - p%met%USR_Z(I-1)) * ( p%met%USR_L(I-1) - p%met%USR_L(I) ) / ( p%met%USR_Z(I-1) - p%met%USR_Z(I) ) + p%met%USR_L(I-1)
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
IF ( ALLOCATED( p%met%USR_Sigma ) ) THEN
   IF ( Ht <= p%met%USR_Z(1) ) THEN
      Sigma = p%met%USR_Sigma(1)
   ELSEIF ( Ht >= p%met%USR_Z(p%met%NumUSRz) ) THEN
      Sigma = p%met%USR_Sigma(p%met%NumUSRz)
   ELSE
      DO I=2,p%met%NumUSRz
         IF ( Ht <= p%met%USR_Z(I) ) THEN
            Sigma = (Ht - p%met%USR_Z(I-1)) * ( p%met%USR_Sigma(I-1) - p%met%USR_Sigma(I) ) / ( p%met%USR_Z(I-1) - p%met%USR_Z(I) ) + p%met%USR_Sigma(I-1)
            EXIT
         ENDIF
      ENDDO
   ENDIF
ELSE
    Sigma = p%met%Ustar*2.15 !bjj: BONNIE, make sure this is defined, or else define ustar for this model...
ENDIF


L1_U   = Lvk/Ucmp
SigmaL1_U = 2.0*Sigma*Sigma*L1_U

DO I=1,p%grid%NumFreq

   FLU2      = ( p%grid%Freq(I)*L1_U )**2
   Tmp       = 1.0 + 71.0*FLU2

   Spec(I,1) = (p%met%USR_StdScale(1)**2)*2.0*SigmaL1_U/Tmp**Exp1
   Spec(I,2) = SigmaL1_U*( 1.0 + 189.0*FLU2 )/Tmp**Exp2
   Spec(I,3) = Spec(I,2)

   Spec(I,2) = (p%met%USR_StdScale(2)**2)*Spec(I,2)
   Spec(I,3) = (p%met%USR_StdScale(3)**2)*Spec(I,3)

ENDDO ! I

RETURN
END SUBROUTINE Spec_vonKrmn
!=======================================================================
!> This subroutine defines the 3-D turbulence spectrum that can be expected to exist upstream of a large, multi-row
!! wind park.  It is based on the smooth or homogeneous terrain models of Hojstrup, Olesen, and Larsen of RISO
!! National Laboratory in Denmark.  The RISO model has been adjusted to reflect the different spectral scaling present
!! in the flow upwind of a large wind park.  The scaling is based on measurements made by the National Renewable Energy
!! Laboratory (NREL) in San Gorgonio Pass, California.
SUBROUTINE Spec_WF_UPW ( p, Ht, Ucmp, Spec )

IMPLICIT                NONE

   ! Passed variables

   TYPE(TurbSim_ParameterType) , INTENT(IN   ) :: p                       !< Input: turbsim parameters
   REAL(ReKi),                   INTENT(IN   ) :: Ht                      !< Height   ( input )
   REAL(ReKi),                   INTENT(IN   ) :: Ucmp                    !< Velocity ( input )
   REAL(ReKi),                   INTENT(  out) :: Spec     (:,:)          !< Target velocity spectra ( output )

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


Ustar2 = p%met%Ustar*p%met%Ustar

IF ( p%met%ZL < 0 ) THEN
      ! BEGIN UNSTABLE FLOW LOOP

   ! Unstable high-frequency range scaling...

   X = - MAX( p%met%ZL, ZL_MinObs)

   Num    = 0.691114  + 0.0791666*X    ! was "Original" Puo_ih =
   Den    = 0.77991   + 0.1761624  / ( 1.0 + EXP( -(X - 0.0405364) / (-0.0184402) ) ) ! was "Measured" Pum_ih =
   Pur_ih = 0.10 * ( Num / Den )
   IF (p%met%ZI > ZI_UVlimit) Pur_ih = (p%met%ZI / ZI_UVlimit) * Pur_ih

   Num    = 0.421958 * EXP( 0.20739895*X )
   Den    = 0.5247865 + 0.0419204 / ( 1.0 + EXP( -(X - 0.0434172) / (-0.0179269) ) )
   Pvr_ih = Num / Den

   Num    = 0.222875  + 0.1347188*X
   Den    = 0.3542331 + 0.0168806 / ( 1.0 + EXP( -(X - 0.0388899) / (-0.0220998) ) )
   Pwr_ih = 0.80 * ( Num / Den )

   Num    = 0.047465  + 0.0132692*X
   Den    = 0.0599494 - 0.0139033*EXP(-X / 0.02603846)
   fur_ih = 1.75 * ( Num / Den )
   IF (p%met%ZI > ZI_UVlimit) fur_ih = (p%met%ZI / ZI_UVlimit)*fur_ih

   Num    = 0.18377384 * EXP( 0.2995136*X )
   Den    = 0.1581509  + 0.09501906*X
   fvr_ih = 1.50 * ( Num / Den )
   IF (p%met%ZI > ZI_UVlimit) fvr_ih = (p%met%ZI / ZI_UVlimit)*fvr_ih

   Num    = 0.3419874 + 0.24985029 * EXP(-X / 0.02619489)
   Den    = 0.451295  + 0.2355227*X
   fwr_ih = 2.0 * ( Num / Den )
   IF (p%met%ZI > ZI_Wlimit) fwr_ih = 0.35*(p%met%ZI / ZI_Wlimit)*fwr_ih


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
   IF (p%met%ZI > ZI_UVlimit) fvr_il = (p%met%ZI / ZI_UVlimit)*fvr_il

   Num    = 3.39482 * EXP( 0.279914*X )
   Den    = 4.59769 + 12.58881*EXP( -X / 0.03351852 )
   fwr_il = 2.25 * ( Num / Den )
   IF (p%met%ZI > ZI_Wlimit) fwr_il=4.0*(p%met%ZI / ZI_Wlimit)*fwr_il

   HtZI  = Ht / p%met%ZI
   HtZI2 = (1.0 - HtZI)**2
   ZInL  = ( p%met%ZI / ( -p%met%L ) )**Exp2
   HtU   = Ht / Ucmp
   ZIU   = p%met%ZI / Ucmp
   UDen  = 1.0 + 15.0*HtZI
   VDen  = 1.0 +  2.8*HtZI
   UDen2 = HtZI2 / UDen**Exp2
   VDen2 = HtZI2 / VDen**Exp2


   DO I = 1,p%grid%NumFreq

         F   = p%grid%Freq(I) * HtU
         Fi  = p%grid%Freq(I) * ZIU

         ! Bonnie: These () around 0.3 HtZI are incorrect as compared to the original SMOOTH model. (For now, leave as is since parameters were-supposedly-calculated with this formulation)
         Fw  = SQRT( (F**2 + (0.3*HtZI**2) ) / (F**2 + 0.0225) )

         Ustar2F = Ustar2 / p%grid%Freq(I)

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

      X = MIN(p%met%ZL, ZL_MaxObs)

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

      DO I = 1,p%grid%NumFreq

         ! CALCULATE Reduced Frequency, f

         f  = p%grid%Freq(I) * Ht / Ucmp
         fq = f / q     ! was XU = f/qu, XV = f/qv, XW = f/qw

         Ustar2F = Ustar2 / p%grid%Freq(I)

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
END SUBROUTINE Spec_WF_UPW
!=======================================================================
!> This subroutine defines the 3-D turbulence spectrum that can be expected to exist (7 to 14 rotor diameters) 
!! downstream of a large, multi-row wind park.  The scaling is based on measurements made by the National 
!! Renewable Energy Laboratory (NREL) in San Gorgonio Pass, California.
SUBROUTINE Spec_WF_DW ( p, Ht, Ucmp, Spec, ErrStat, ErrMsg )


IMPLICIT                NONE

   ! Passed variables

   TYPE(TurbSim_ParameterType) , INTENT(IN   ) :: p                       !< Input: turbsim parameters
   REAL(ReKi),                   INTENT(IN   ) :: Ht                      !< Height   ( input )
   REAL(ReKi),                   INTENT(IN   ) :: Ucmp                    !< Velocity ( input )
   REAL(ReKi),                   INTENT(  out) :: Spec     (:,:)          !< Target velocity spectra ( output )

   INTEGER(IntKi),              INTENT(OUT)    :: ErrStat
   CHARACTER(*),                INTENT(OUT)    :: ErrMsg
   
   
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

   
   ErrStat = ErrID_None
   ErrMsg  = ""

Ustar2 = p%met%Ustar*p%met%Ustar

IF (p%met%ZL < 0) THEN


      ! Get Unstable spectral peaks

   ! Unstable high-frequency range scaling...

   X = - MAX( p%met%ZL, ZL_MinObs )

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
   IF (p%met%ZI < 1350.0 ) fwr_oh = (p%met%ZI / 1350.0) * fwr_oh

      ! Unstable low-frequency range scaling ...

   Num    = 0.796264 + 0.316895 / (1.0 + EXP( -(X - 0.082483) / 0.027480 ))
   Den    = 0.07616  + EXP(0.303919 * X * 0.390906)   ! was Den = 0.07616 + EXP(0.303919*X)**0.390906
   Pur_ol = 4.0 * (Num / Den)
   IF (p%met%ZI < 1600.0) Pur_ol = (p%met%ZI / 1600.0) * Pur_ol

   Num    = 0.812483 + 0.1332134 * X
   Den    = 0.104132 + EXP(0.714674 * X * 0.495370)   ! was Den = 0.104132 + EXP(0.714674*X)**0.495370
   Pvr_ol = Num / Den
   Pvr_ol = (p%met%ZI / 1600.0)*Pvr_ol

   Num    = 0.371298  + 0.0425447 * X
   Den    = 0.0004375 + EXP(0.4145751 * X * 0.6091557)   ! was Den = 0.0004375 + EXP(0.4145751*X)**0.6091557
   Pwr_ol = 0.75 * (Num / Den)

   Num    = 0.859809 * EXP(0.157999 * X)
   Den    = 0.81459 + 0.021942 * X
   fur_ol = 1.5 * (Num / Den)
   IF (p%met%ZI > 1850.0) fur_ol = 2.6 * (p%met%ZI / 1850.0) * fur_ol

   !A0 =  0.8121775
   !A1 =  4.122E+15
   !A2 = -0.594909
   !A3 =  0.0559581
   Num    = 0.8121775 * EXP( -X / 4.122E+15 ) - 0.594909 * EXP( -X / 0.0559581 ) ! was Num = BETA4(X,A0,A1,A2,A3)
   Den    = 0.72535  - 0.0256291 * X
   fvr_ol = 3.0 * (Num / Den)
   fvr_ol = (p%met%ZI / 1600.0) * fvr_ol

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

   HtZI  = Ht / p%met%ZI
   HtZI2 = (1.0 - HtZI)**2
   ZInL  = ( p%met%ZI / (-p%met%L) )**Exp2
   UDen  = 1.0 + 15.0 * HtZI
   VDen  = 1.0 +  2.8 * HtZI
   UDen2 = HtZI2 / UDen**Exp2
   VDen2 = HtZI2 / VDen**Exp2

   DO I = 1,p%grid%NumFreq

      ! Calculate f,fi,fru,frv

      F   = p%grid%Freq(I)*Ht / Ucmp
      Fi  = p%grid%Freq(I)*p%met%ZI / Ucmp
      Fw  = SQRT( (F**2 + (0.3*HtZI**2)) / (F**2 + 0.0225) )

      Ustar2F = Ustar2 / p%grid%Freq(I)

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

   X      = MIN( p%met%ZL, ZL_MaxObs )

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

   DO I = 1,p%grid%NumFreq

      F       = p%grid%Freq(I) * Ht / Ucmp    ! Reduced frequency

      Fq      = F / q
      Ustar2F = Ustar2 / p%grid%Freq(I)

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


CONTAINS
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
         CALL SetErrStat( ErrID_Fatal, 'Invalid X or A inputs.', ErrStat, ErrMsg, 'Beta6' )
         RETURN
      ENDIF

      IF ( X < A + 1.0 ) THEN
         CALL Beta8( Beta6, A, X ) 
         ! Beta6 = GAMSER
      ELSE
         CALL Beta7( Beta6, A, X)  
         Beta6 = 1.0 - Beta6 
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
      CALL SetErrStat( ErrID_Fatal, 'Input variable X must be positive.', ErrStat, ErrMsg, 'Beta7' )
      RETURN
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
         CALL SetErrStat( ErrID_Fatal, 'Value of A is too large or ITMAX is too small.', ErrStat, ErrMsg, 'Beta7' )
         RETURN
      ENDIF

      GLN  = Beta9( A )
      IF (ErrStat >= AbortErrLev) RETURN

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
            CALL SetErrStat( ErrID_Fatal, 'Value of A is too large or ITMAX is too small.', ErrStat, ErrMsg, 'BETA8' )
            RETURN
         ENDIF

         GLN = Beta9( A )
         IF (ErrStat >= AbortErrLev) RETURN
         
         GAMSER = Sum * EXP( -X + A * LOG(X) - GLN)

      ELSEIF ( X == 0.0 ) THEN

            GAMSER = 0.0

      ELSE ! ( X < 0.0 )
         CALL SetErrStat( ErrID_Fatal, 'Invalid input.', ErrStat, ErrMsg, 'BETA8' )
         RETURN
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

   REAL(ReKi), PARAMETER  :: Cof(6) = (/ 76.18009173, -86.50532033, 24.01409822, -1.231739516, 0.120858003E-2, -0.536382E-5 /)
   REAL(ReKi), PARAMETER  :: STP  = 2.50662827465

   INTEGER                :: J       ! Loop counter

   IF ( XX <= -4.5 ) THEN
      CALL SetErrStat( ErrID_Fatal, 'Input variable XX must be larger than -4.5.', ErrStat, ErrMsg, 'Beta9' )
      RETURN      
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
      CALL SetErrStat( ErrID_Fatal, 'Variable SER must be larger than 0.0.', ErrStat, ErrMsg, 'Beta9' )
      RETURN      
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
END SUBROUTINE Spec_WF_DW

!=======================================================================
END MODULE TS_VelocitySpectra
