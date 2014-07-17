!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2014  National Renewable Energy Laboratory
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
MODULE TS_Profiles

   USE                     NWTC_Library   
   
use ts_errors

   IMPLICIT                NONE

      ! Create interface for a generic getWindSpeed that actually uses specific routines.

   INTERFACE getWindSpeed
      MODULE PROCEDURE getWindSpeedVal
      MODULE PROCEDURE getWindSpeedAry
   END INTERFACE


CONTAINS

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
!      CALL TS_Abort( 'The x and y vectors in ChebyshevVals() must be the same size.' )
      CALL ProgWarn( 'ChebyshevVals: The x and y vectors must be the same size.' )
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

FUNCTION getWindSpeedAry(URef, RefHt, Ht, RotorDiam, profile, UHangle)

   ! Determine the wind speed at a given height, with reference wind speed.
   ! Use power law if given height and reference height are within rotor disk.
   ! Use log profile if reference height is below rotor disk.

   USE                                  TSMods, ONLY: ChebyCoef_WD   ! Chebyshev coefficients (must be defined before calling this function!)
   USE                                  TSMods, ONLY: ChebyCoef_WS   ! Chebyshev coefficients (must be defined before calling this function!)
   USE                                  TSMods, ONLY: HFlowAng       ! The horizontal flow angle of the mean wind speed at hub height
   USE                                  TSMods, ONLY: HubHt          ! The hub height (used with HFlowAng)
   USE                                  TSMods, ONLY: H_Ref          ! The reference height (RefHt is being overwritten for some reason)
   USE                                  TSMods, ONLY: IEC_EWM1       ! IEC Extreme 1-yr wind speed model
   USE                                  TSMods, ONLY: IEC_EWM50      ! IEC Extreme 50-yr wind speed model
   USE                                  TSMods, ONLY: IEC_EWM100      ! IEC Extreme 100-yr wind speed model
   USE                                  TSMods, ONLY: IEC_WindType   ! Type of IEC wind
   USE                                  TSMods, ONLY: NumUSRz        ! Number of user-defined heights
   USE                                  TSMods, ONLY: PLExp          ! Power law exponent
   USE                                  TSMods, ONLY: U_Ref          ! The input wind speed at the reference height (URef is being overwritten for some reason)
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

   REAL                               :: C_factor
   REAL(ReKi)                         :: ZRef


   IF ( IEC_WindType == IEC_EWM50 ) THEN
      getWindSpeedAry(:) = VRef*( Ht(:)/HubHt )**0.11                ! [IEC 61400-1 6.3.2.1 (14)]  !bjj: this PLExp should be set to 0.11 already, why is this hard-coded?
      RETURN
   ELSEIF ( IEC_WindType == IEC_EWM1 ) THEN
      getWindSpeedAry(:) = 0.8*VRef*( Ht(:)/HubHt )**0.11            ! [IEC 61400-1 6.3.2.1 (14), (15)]
      RETURN
   ELSEIF ( IEC_WindType == IEC_EWM100 ) THEN
      getWindSpeedAry(:) = VRef*( Ht(:)/HubHt )**0.11               ! [API-IEC RECCOMENDATAION]  ADDED BY YGUO
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
     CASE ( 'API', 'A' ) !Panofsky, H.A.; Dutton, J.A. (1984). Atmospheric Turbulence: Models and Methods for Engineering Applications. New York: Wiley-Interscience; 397 pp.

     ! sample to write to screen.CALL WrScr ( '    A default value will be used for '//TRIM(VarName)//'.' )
!            CALL WrScr ('calling to API wind profile and write array')
                ! TO ADD THE FSOLVE PROGRAM TO CALCULATE C_factor
!            getWindSpeedAry(:)  = ZRef*(1+LOG( Ht(:) / 10) )*( 1-0.41*(0.06*(1+0.043*ZRef)*(Ht/10)**(-0.22)))*LOG(600.0/3600.0)
!            getWindSpeedAry(:)  = RefHt*(1+LOG( Ht(:) / 10) )*( 1-0.41*(0.06*(1+0.043*RefHt)*(Ht/10)**(-0.22)))*LOG(600.0/3600.0)
!            getWindSpeedAry(:)  = URef*(1+LOG( Ht(:) / RefHt) )*( 1-0.41*(0.06*(1+0.043*URef)*(Ht/RefHt)**(-0.22)))*LOG(600.0/3600.0)
!            getWindSpeedAry(:)  = URef*(1+LOG( Ht(:) / RefHt) )
!            getWindSpeedAry(:) = URef*( 1.0 + 0.0573*SQRT( 1.0 + 0.15*URef )*LOG( Ht(:)/RefHt) )
            getWindSpeedAry(:) = U_Ref*( 1.0 + 0.0573*SQRT( 1.0 + 0.15*U_Ref )*LOG( Ht(:)/H_Ref) )
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
   USE                                  TSMods, ONLY: H_Ref          ! The reference height (RefHt is being overwritten for some reason)
   USE                                  TSMods, ONLY: IEC_EWM1       ! IEC Extreme 1-yr wind speed model
   USE                                  TSMods, ONLY: IEC_EWM50      ! IEC Extreme 50-yr wind speed model
   USE                                  TSMods, ONLY: IEC_WindType   ! Type of IEC wind
   USE                                  TSMods, ONLY: NumUSRz        ! Number of user-defined heights
   USE                                  TSMods, ONLY: PLExp          ! Power law exponent
   USE                                  TSMods, ONLY: U_Ref          ! The input wind speed at the reference height (URef is being overwritten for some reason)
   USE                                  TSMods, ONLY: U_Usr          ! User-defined wind speeds
   USE                                  TSMods, ONLY: Ustar          ! Friction or shear velocity
   USE                                  TSMods, ONLY: Vref           ! IEC Extreme wind value
   USE                                  TSMods, ONLY: WindDir_USR    ! User-defined wind directions
   USE                                  TSMods, ONLY: Z_Usr          ! User-defined heights
   USE                                  TSMods, ONLY: z0             ! Surface roughness length -- It must be > 0 (which we've already checked for)
   USE                                  TSMods, ONLY: ZJetMax        ! Height of the low-level jet
   USE                                  TSMods, ONLY: ZL             ! M-O stability parameter
   USE                                  TSMods, ONLY: U0_1HR          !% ADDED BY Y. GUO
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

   !REAL(ReKi)                   :: U0_1HR                        ! Wind speed at reference height in 1hr time duration, added by Y.G. ON April 16 2013
   REAL                         :: C_factor                       ! Factor to convert wind speed from 10-min to 1 hr, added by Y.G. ON April 16 2013
   REAL :: X0                                              ! Added by Y. Guo for calculating C_factor
   REAL :: X, TEMP_1, TEMP_2                                                ! Added by Y. Guo for calculating C_factor
   !=======================================

      X0=URef
      !========================
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


      CASE ( 'PL', 'P' )  ! POWER LAW, commented by Y. Guo on April 16 2013

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

      CASE ( 'API', 'A' ) !Panofsky, H.A.; Dutton, J.A. (1984). Atmospheric Turbulence: Models and Methods for Engineering Applications. New York: Wiley-Interscience; 397 pp.

!MLB: We can exclude this logic by forcing the user to enter the 1-hour mean wind speed.
!     If we add the API stuff to the main version of TurbSim, we may want to eliminate that requirement, but we will have to
!     add a new input parameter saying how long a period was used to calculate the reference wind speed.

!           CALL ROOT_SEARCHING(X0,X,URef,RefHt,HubHt)  !URef
           !CALL Root_Searching(X0,X,42.5,10.0,10.0)  !URef, USED TO DEBUG THE CODE
!           U0_1HR=X  ! This is the wind speed at 10 m height within 1-hr window

!            CALL WrScr ('Calling to API wind profile')
!           TEMP_1=0.0573*(1.0+0.15*U0_1HR)**0.5
!           TEMP_2=0.06*(1+0.043*U0_1HR)*(Ht/10.0)**(-0.22)
!           getWindSpeedVal = U0_1HR*(1.0+TEMP_1*LOG( Ht / 10.0) )*( 1.0-0.41*TEMP_2*LOG(600.0/3600.0))
!           getWindSpeedVal = U0_1HR*( 1.0 + 0.0573*SQRT( 1.0 + 0.15*U0_1HR )*LOG( Ht/RefHt) )
!MLB: This assumes that the reference wind speed entered by the user is the 1-hour average wind speed at the input reference height.
           getWindSpeedVal = U_Ref*( 1.0 + 0.0573*SQRT( 1.0 + 0.15*U_Ref )*LOG( Ht/H_Ref) )

!            CALL WrScr ('API wind profile generated')

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
END MODULE TS_Profiles
