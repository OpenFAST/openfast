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
MODULE TS_Profiles

   USE                     NWTC_Library   
   USE                     TurbSim_Types

   IMPLICIT                NONE


CONTAINS

!=======================================================================
SUBROUTINE ChebyshevVals(coeffs,x,y,MinX,MaxX, ErrStat, ErrMsg)

   IMPLICIT NONE

      ! Passed variables

   REAL(ReKi), DIMENSION(:),   INTENT(IN   ) :: coeffs  ! Coefficients defined on [-1,1]
   REAL(ReKi), DIMENSION(:),   INTENT(IN   ) :: x       ! The x values where f(x)=y is desired
   REAL(ReKi), DIMENSION(:),   INTENT(  OUT) :: y       ! The desired function values
   REAL(SiKi),                 INTENT(IN   ) :: MinX    ! Min X of the interval the coeffs were originally calculated for
   REAL(SiKi),                 INTENT(IN   ) :: MaxX    ! Max X of the interval the coeffs were originally calculated for
   INTEGER(IntKi),             INTENT(  OUT) :: ErrStat ! Error level
   CHARACTER(*),               INTENT(  OUT) :: ErrMsg  ! Message describing error

   
   
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
      ErrStat = ErrID_Warn
      ErrMsg = 'ChebyshevVals:The x and y vectors must be the same size.' 
      SX = MIN(SX,SY)
      SY = SX
   ELSE
      ErrStat = ErrID_None
      ErrMsg  = ""
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
!> This subroutine determines what Chebyshev Coefficients will be used
!! for the jet wind speed and wind direction profiles
SUBROUTINE GetChebCoefs(p, UJetMax_IsKnown, ErrStat, ErrMsg)

   ! sets p%met%ChebyCoef_WS, p%met%ChebyCoef_WD, and 
   !    if .NOT. UJetMax_IsKnown, also sets p%met%UJetMax

   ! valid only for jet WindProfileType
   
IMPLICIT                   NONE

   TYPE(TurbSim_ParameterType),INTENT(INOUT) :: p                 ! TurbSim parameters
   LOGICAL,                    INTENT(IN)    :: UJetMax_IsKnown   
   INTEGER(IntKi),             intent(  out) :: ErrStat           !< Error level
   CHARACTER(*),               intent(  out) :: ErrMsg            !< Message describing error



   REAL(ReKi)              :: UH_coef(4,11)     ! The coefficients that Neil developed for calculating the Chebyshev coefficients
   REAL(ReKi)              :: WD_coef(4,11)     ! The coefficients that Neil developed for calculating the Chebyshev coefficients
   REAL(ReKi)              :: ChebyCoef_tmp(11)
   REAL(ReKi)              :: UTmp1             !
   REAL(ReKi)              :: UTmp2             !

   INTEGER                 :: I                 ! A loop counter


   ErrStat = ErrID_None
   ErrMsg  = ""
   
      ! Let's calculate the wind speed at the jet height

   CALL get_coefs(p%met%ZJetMax, UH_coef, WD_coef)


   IF ( UJetMax_IsKnown ) THEN

      DO I=1,11
         p%met%ChebyCoef_WS(I) = p%met%UJetMax*UH_coef(1,I) + p%met%Rich_No*UH_coef(2,I) &
                               + p%met%Ustar  *UH_coef(3,I) +               UH_coef(4,I)
      ENDDO

   ELSE

         ! Calculate the coefficients without UJetMax

      DO I=1,11
         p%met%ChebyCoef_WS(I) = p%met%Rich_No*UH_coef(2,I) + p%met%Ustar*UH_coef(3,I) + UH_coef(4,I) ! +UJetMax*UH_coef(1,I)
      ENDDO

      CALL getVelocity(p, p%met%URef, p%met%RefHt, p%met%RefHt, Utmp1, ErrStat, ErrMsg) ! reference p%met%URef, p%met%RefHt are unused; return velocity at RefHt with the WS coeffs missing the UJetMax term
         IF (ErrStat >= AbortErrLev) RETURN
         
         ! Now calculate the coefficients with just UJetMax term
         
      ChebyCoef_tmp(:)         = p%met%ChebyCoef_WS(:)
      p%met%ChebyCoef_WS(:)    = UH_coef(1,:)

      CALL getVelocity(p, p%met%URef, p%met%RefHt, p%met%RefHt, Utmp2, ErrStat, ErrMsg) ! reference p%met%URef, p%met%RefHt are unused
         IF (ErrStat >= AbortErrLev) RETURN
         
         ! this gives us UJetMax:
      p%met%UJetMax            = (p%met%URef - Utmp1)/Utmp2

         ! Get the final coefficients, using the computed UJetMax
      p%met%ChebyCoef_WS(:)    = p%met%UJetMax*p%met%ChebyCoef_WS(:) + ChebyCoef_tmp(:)

   ENDIF

   DO I=1,11
      p%met%ChebyCoef_WD(I)    = p%met%UJetMax*WD_coef(1,I) + p%met%Rich_No*WD_coef(2,I) &
                               + p%met%Ustar*WD_coef(3,I) +                 WD_coef(4,I)
   ENDDO


RETURN
END SUBROUTINE GetChebCoefs
!=======================================================================
!> This subroutine sets the array VelocityProfile, which contains the velocies
!! at each height specified by the input array Ht.
SUBROUTINE getVelocityProfile(p, U_Ref, z_Ref, Ht, VelocityProfile, ErrStat, ErrMsg )


   ! Determine the wind speed at a given height, with reference wind speed.
   
   IMPLICIT                              NONE

   TYPE(TurbSim_ParameterType),INTENT(IN)    :: p                           !< TurbSim  parameters
   REAL(ReKi),                 INTENT(IN)    :: U_Ref                       !< Velocity at reference height
   REAL(ReKi),                 INTENT(IN)    :: z_Ref                       !< Reference height
   REAL(ReKi),                 INTENT(IN)    :: Ht(:)                       !< Heights (array) in meters where wind/water velocity should be calculated
   REAL(ReKi),                 INTENT(  OUT) :: VelocityProfile(:)          !< Calculated velocity (wind/water speed) at Ht
   INTEGER(IntKi),             intent(  out) :: ErrStat                     !< Error level
   CHARACTER(*),               intent(  out) :: ErrMsg                      !< Message describing error

   
   
   REAL(SiKi),   PARAMETER                :: MinZ = 3.                   ! lower bound (height) for Cheby polynomial
   REAL(SiKi),   PARAMETER                :: MaxZ = 500.                 ! upper bound (height) for Cheby polynomial


   INTEGER                                :: I
   INTEGER                                :: Indx
   INTEGER                                :: J

!   REAL                                   :: C_factor
!   REAL(ReKi)                             :: ZRef

   ErrStat = ErrID_None
   ErrMsg  = ""

   IF ( p%IEC%IEC_WindType == IEC_EWM50 ) THEN
      VelocityProfile(:) = p%IEC%VRef*( Ht(:)/p%grid%HubHt )**p%met%PLExp                ! [IEC 61400-1 6.3.2.1 (14)]  
      RETURN
   ELSEIF ( p%IEC%IEC_WindType == IEC_EWM1 ) THEN
      VelocityProfile(:) = 0.8*p%IEC%VRef*( Ht(:)/p%grid%HubHt )**p%met%PLExp            ! [IEC 61400-1 6.3.2.1 (14), (15)]
      RETURN
   ELSEIF ( p%IEC%IEC_WindType == IEC_EWM100 ) THEN
      VelocityProfile(:) = p%IEC%VRef*( Ht(:)/p%grid%HubHt )**p%met%PLExp               ! [API-IEC RECCOMENDATAION]  ADDED BY YGUO  !bjj: this is the same as IEC_EWM50, but we should check that IEC_EWM100 is used in ALL the same places IEC_EWM50 is
      RETURN
   ENDIF


   SELECT CASE ( TRIM(p%met%WindProfileType) )

      CASE ( 'JET' )

         CALL ChebyshevVals( p%met%ChebyCoef_WS, Ht, VelocityProfile, MinZ, MaxZ, ErrStat, ErrMsg ) ! We originally calculated the coeffs for 3-500 m in height

      CASE ( 'LOG' )

            DO J = 1,SIZE(Ht)
               VelocityProfile(J) = getLogWindSpeed( Ht(J), z_Ref, U_Ref, p%met%ZL, p%met%Z0)
            END DO
            
            
            
      CASE ( 'H2L' )

               ! Calculate the windspeed.
               ! z_Ref and U_Ref both get modified consistently, therefore z_Ref is used instead of RefHt.
               VelocityProfile(:) = LOG( Ht(:)/z_Ref) * p%met%Ustar / 0.41_ReKi + U_Ref

      CASE ( 'PL' )

!         IF ( z_Ref > 0.0 .AND. Ht > 0.0 ) THEN
            VelocityProfile(:) = U_Ref*( Ht(:)/z_Ref )**p%met%PLExp

!         ENDIF

      CASE ( 'TS' )

            DO J = 1,SIZE(Ht)
               VelocityProfile(J) =  getTimeSeriesWindSpeed(p, Ht(J) )   
            END DO
         

     CASE ( 'API' )

     ! sample to write to screen.CALL WrScr ( '    A default value will be used for '//TRIM(VarName)//'.' )
!            CALL WrScr ('calling to API wind profile and write array')
                ! TO ADD THE FSOLVE PROGRAM TO CALCULATE C_factor
!            VelocityProfile(:)  = ZRef*(1+LOG( Ht(:) / 10) )*( 1-0.41*(0.06*(1+0.043*ZRef)*(Ht/10)**(-0.22)))*LOG(600.0/3600.0)
!            VelocityProfile(:)  = z_Ref*(1+LOG( Ht(:) / 10) )*( 1-0.41*(0.06*(1+0.043*z_Ref)*(Ht/10)**(-0.22)))*LOG(600.0/3600.0)
!            VelocityProfile(:)  = U_Ref*(1+LOG( Ht(:) / z_Ref) )*( 1-0.41*(0.06*(1+0.043*U_Ref)*(Ht/z_Ref)**(-0.22)))*LOG(600.0/3600.0)
!            VelocityProfile(:)  = U_Ref*(1+LOG( Ht(:) / z_Ref) )
!            VelocityProfile(:) = U_Ref*( 1.0 + 0.0573*SQRT( 1.0 + 0.15*p%met%URef )*LOG( Ht(:)/z_Ref) )
            VelocityProfile(:) = p%met%URef*( 1.0 + 0.0573*SQRT( 1.0 + 0.15*p%met%URef )*LOG( Ht(:)/p%met%RefHt) )
      CASE ( 'USR' )

         DO J = 1,SIZE(Ht)
            IF ( Ht(J) <= p%met%USR_Z(1) ) THEN
               VelocityProfile(J) = p%met%USR_U(1)
            ELSEIF ( Ht(J) >= p%met%USR_Z(p%met%NumUSRz) ) THEN
               VelocityProfile(J) = p%met%USR_U(p%met%NumUSRz)
            ELSE
               ! Find the two points between which the height lies

               DO I=2,p%met%NumUSRz
                  IF ( Ht(J) <= p%met%USR_Z(I) ) THEN
                     Indx = I-1

                     ! Let's just do a linear interpolation for now
                     VelocityProfile(J) = (Ht(J) - p%met%USR_Z(Indx)) * ( p%met%USR_U(Indx) - p%met%USR_U(I) ) / ( p%met%USR_Z(Indx) - p%met%USR_Z(I) ) &
                                        + p%met%USR_U(Indx)
                     EXIT
                  ENDIF
               ENDDO

            ENDIF

         ENDDO

      CASE DEFAULT   ! This is how it worked before

         DO I=1,SIZE(VelocityProfile)
            IF ( Ht(I) == z_Ref ) THEN
               VelocityProfile(I) = U_Ref
            ELSEIF ( ABS( Ht(I)-z_Ref ) <= 0.5*p%grid%RotorDiameter ) THEN
               VelocityProfile(I) = U_Ref*( Ht(I)/z_Ref )**p%met%PLExp
            ELSEIF ( Ht(I) > 0.0 .AND. z_Ref > 0.0 .AND. .NOT. EqualRealNos(z_Ref, p%met%Z0) ) THEN !Check that we don't have an invalid domain
               VelocityProfile(I) = U_Ref*LOG( Ht(I)/p%met%Z0 )/LOG( z_Ref/p%met%Z0 )
            ELSE
               VelocityProfile(I) = 0.0
            ENDIF
         ENDDO

   END SELECT

RETURN
END SUBROUTINE getVelocityProfile
!=======================================================================
!>  This subroutine sets the direction in degrees at each height in meters 
!!  specified by the input array Ht.
SUBROUTINE getDirectionProfile( p, Ht, DirectionProfile, VAngleProfile, ErrStat, ErrMsg )
  
   IMPLICIT                              NONE

   TYPE(TurbSim_ParameterType),INTENT(IN)    :: P                           !< TurbSim parameters  
   REAL(ReKi),                 INTENT(IN)    :: Ht(:)                       !< Array of heights (meters) where wind speed should be calculated
   REAL(ReKi)    ,             intent(  out) :: DirectionProfile(:)         !< Wind direction at Ht
   REAL(ReKi)    ,             intent(  out) :: VAngleProfile(:)            !< Vertical Wind angle at Ht
   INTEGER(IntKi),             intent(  out) :: ErrStat                     !< Error level
   CHARACTER(*),               intent(  out) :: ErrMsg                      !< Message describing error

   
   REAL(SiKi),   PARAMETER                :: MinZ = 3.                   ! lower bound (height) for Cheby polynomial
   REAL(SiKi),   PARAMETER                :: MaxZ = 500.                 ! upper bound (height) for Cheby polynomial


   REAL(ReKi)                             :: tmpHt(2)
   REAL(ReKi)                             :: tmpWD(2)
   REAL(ReKi)                             :: diff
                                          
   INTEGER                                :: IZ, IZm1
   INTEGER                                :: J
                                          

   ErrStat = ErrID_None
   ErrMsg  = ""

   SELECT CASE ( TRIM(p%met%WindProfileType) )

      CASE ( 'JET' )

            ! Calculate the wind direction at this height
         CALL ChebyshevVals( p%met%ChebyCoef_WD, Ht(:), DirectionProfile(:), MinZ, MaxZ, ErrStat, ErrMsg )
         IF (ErrStat >= AbortErrLev) RETURN
         
            ! Compute the wind direction at hub height & the jet height
         tmpHt(1) = p%grid%HubHt
         tmpHt(2) = p%met%ZJetMax
         CALL ChebyshevVals( p%met%ChebyCoef_WD, tmpHt, tmpWD(1:2), MinZ, MaxZ, ErrStat, ErrMsg )
         IF (ErrStat >= AbortErrLev) RETURN

            ! Make sure none of the directions are more than 45 degrees from the direction at the jet height
         IF ( ABS(tmpWD(1) - tmpWD(2) ) > 45. ) THEN  ! The direction at the hub height
            tmpWD(1) = tmpWD(2) + SIGN(REAL(45.,ReKi), tmpWD(1) - tmpWD(2))
         ENDIF

         DO J = 1,SIZE(DirectionProfile) ! The directions at all the heights
            IF ( ABS(DirectionProfile(J) - tmpWD(2) ) > 45. ) THEN
               DirectionProfile(J) = tmpWD(2) + SIGN(REAL(45.,ReKi), DirectionProfile(J) - tmpWD(2))
            ENDIF

            ! Remove the hub height direction so that we have a relative direction, then
            ! add the mean flow angle. (Note that the Chebyshev profile is cw looking upwind,
            ! but the horizontal angle is ccw looking upwind)

            DirectionProfile(J) = p%met%HFlowAng - (DirectionProfile(J) - tmpWD(1)) ! This is the counter-clockwise angle of the wind
         ENDDO
         VAngleProfile = p%met%VFlowAng

      CASE ( 'USR' )

         DO J = 1,SIZE(Ht)

               ! Calculate the wind direction at this height

            IF ( Ht(J) <= p%met%USR_Z(1) ) THEN
               DirectionProfile(J) = p%met%USR_WindDir(1)
            ELSEIF ( Ht(J) >= p%met%USR_Z(p%met%NumUSRz) ) THEN
               DirectionProfile(J) = p%met%USR_WindDir(p%met%NumUSRz)
            ELSE
               
               
               ! Find the two points between which the height lies

               DO IZ=2,p%met%NumUSRz
                  IF ( Ht(J) <= p%met%USR_Z(IZ) ) THEN
                     IZm1 = IZ-1

                        ! Let's just do a linear interpolation for now
                     !we need to check if the angle goes through 360, before we do the interpolation
                     diff = p%met%USR_WindDir(IZm1) - p%met%USR_WindDir(IZ)
                     IF ( diff > 180. ) THEN
                        tmpWD(1) = p%met%USR_WindDir(IZm1)
                        tmpWD(2) = p%met%USR_WindDir(IZ  ) + 360.
                     ELSEIF ( diff < -180. ) THEN
                        tmpWD(1) = p%met%USR_WindDir(IZm1) + 360.
                        tmpWD(2) = p%met%USR_WindDir(IZ  )
                     ELSE
                        tmpWD(1) = p%met%USR_WindDir(IZm1)
                        tmpWD(2) = p%met%USR_WindDir(IZ  )
                     ENDIF

                     DirectionProfile(J) = (Ht(J) - p%met%USR_Z(IZm1)) * ( tmpWD(1) - tmpWD(2) ) / ( p%met%USR_Z(IZm1) - p%met%USR_Z(IZ) ) + tmpWD(1)

                     EXIT
                  ENDIF
               ENDDO
               
               
            ENDIF

         ENDDO
         
!bjj: TODO: See if we can get this to have direction of HFlowAng at hub height.
         DirectionProfile = p%met%HFlowAng + DirectionProfile  ! This is the counter-clockwise angle of the wind
         VAngleProfile    = p%met%VFlowAng
         
      CASE ('TS')
         
         DO J = 1,SIZE(Ht)

               ! Calculate the wind direction at this height
         
            IF ( Ht(J) <= p%usr%pointzi(1) ) THEN
               DirectionProfile(J) = p%usr%meanDir(1)
               VAngleProfile(J)    = p%usr%meanVAng(1)
            ELSEIF ( Ht(J) >= p%usr%pointzi(p%usr%NPoints) ) THEN
               DirectionProfile(J) = p%usr%meanDir(p%usr%NPoints)
               VAngleProfile(J)    = p%usr%meanVAng(p%usr%NPoints)
            ELSE
               ! Find the two points between which the height lies

               DO IZ=2,p%usr%NPoints
                  IF ( Ht(J) <= p%usr%pointzi(IZ) ) THEN
                     IZm1 = IZ-1

                        ! Let's just do a linear interpolation for now
                     !we need to check if the angle goes through 360, before we do the interpolation
                     diff = p%usr%meanDir(IZm1) - p%usr%meanDir(IZ)
                     IF ( diff > 180. ) THEN
                        tmpWD(1) = p%usr%meanDir(IZm1)
                        tmpWD(2) = p%usr%meanDir(IZ  ) + 360.
                     ELSEIF ( diff < -180. ) THEN
                        tmpWD(1) = p%usr%meanDir(IZm1) + 360.
                        tmpWD(2) = p%usr%meanDir(IZ  )
                     ELSE
                        tmpWD(1) = p%usr%meanDir(IZm1)
                        tmpWD(2) = p%usr%meanDir(IZ  )
                     ENDIF

                     DirectionProfile(J) = (Ht(J) - p%usr%pointzi(IZm1)) * ( tmpWD(1) - tmpWD(2) ) / ( p%usr%pointzi(IZm1) - p%usr%pointzi(IZ) ) + tmpWD(1)    
                     
                     
                        ! Let's just do a linear interpolation for now
                     !we need to check if the angle goes through 360, before we do the interpolation
                     diff = p%usr%meanVAng(IZm1) - p%usr%meanVAng(IZ)
                     IF ( diff > 180. ) THEN
                        tmpWD(1) = p%usr%meanVAng(IZm1)
                        tmpWD(2) = p%usr%meanVAng(IZ  ) + 360.
                     ELSEIF ( diff < -180. ) THEN
                        tmpWD(1) = p%usr%meanVAng(IZm1) + 360.
                        tmpWD(2) = p%usr%meanVAng(IZ  )
                     ELSE
                        tmpWD(1) = p%usr%meanVAng(IZm1)
                        tmpWD(2) = p%usr%meanVAng(IZ  )
                     ENDIF
                                          
                     VAngleProfile(J)    = (Ht(J) - p%usr%pointzi(IZm1)) * ( tmpWD(1) - tmpWD(2) ) / ( p%usr%pointzi(IZm1) - p%usr%pointzi(IZ) ) + tmpWD(1)
                     
                     EXIT
                  ENDIF
               ENDDO

            ENDIF       
         
         END DO

         DirectionProfile = p%met%HFlowAng + DirectionProfile  ! This is the counter-clockwise angle of the wind         
         VAngleProfile    = p%met%VFlowAng + VAngleProfile
         
      CASE DEFAULT   

         DirectionProfile = p%met%HFlowAng
         VAngleProfile    = p%met%VFlowAng
         
   END SELECT

RETURN
END SUBROUTINE getDirectionProfile
!=======================================================================
!> This subroutine sets the scalar Velocity, which contains the velocity in m/s
!! at the height in mebers specified by the input value Ht.
SUBROUTINE getVelocity(p, U_Ref, z_Ref, Ht, Velocity, ErrStat, ErrMsg )
   
   IMPLICIT                              NONE

   TYPE(TurbSim_ParameterType), INTENT(IN)    :: P
   REAL(ReKi),                  INTENT(IN)    :: U_Ref                       ! Wind speed at reference height
   REAL(ReKi),                  INTENT(IN)    :: z_Ref                       ! Reference height
   REAL(ReKi),                  INTENT(IN)    :: Ht                          ! Height where wind speed should be calculated
   REAL(ReKi)    ,              intent(  out) :: Velocity                    ! This function, approximate wind/water speed at Ht
   INTEGER(IntKi),              intent(  out) :: ErrStat                     !< Error level
   CHARACTER(*),                intent(  out) :: ErrMsg                      !< Message describing error
                                           
   
   REAL(SiKi),   PARAMETER                 :: MinZ = 3.                   ! lower bound (height) for Cheby polynomial
   REAL(SiKi),   PARAMETER                 :: MaxZ = 500.                 ! upper bound (height) for Cheby polynomial
                                           
                                           
!   REAL(ReKi)                              :: psiM                        ! The diabatic term for the log wind profile
   REAL(ReKi)                              :: tmpHt(2)
   REAL(ReKi)                              :: tmpWS(2)
                                           
   INTEGER                                 :: I
   INTEGER                                 :: Indx

   !REAL(ReKi)                   :: U0_1HR                        ! Wind speed at reference height in 1hr time duration, added by Y.G. ON April 16 2013
!   REAL                         :: C_factor                       ! Factor to convert wind speed from 10-min to 1 hr, added by Y.G. ON April 16 2013
   REAL :: X0                                              ! Added by Y. Guo for calculating C_factor
   !REAL :: X, TEMP_1, TEMP_2                                                 ! Added by Y. Guo for calculating C_factor
   !=======================================

   ErrStat = ErrID_None
   ErrMsg  = ""
   
      X0=U_Ref
      !========================
   ! IF p%met%Z0 <= 0.0    CALL ProgAbort('The surface roughness must be a positive number')

   IF ( p%IEC%IEC_WindType == IEC_EWM50 ) THEN
      Velocity =    p%IEC%VRef*( Ht/p%grid%HubHt )**p%met%PLExp                   ! [IEC 61400-1 6.3.2.1 (14)]
      RETURN
   ELSEIF ( p%IEC%IEC_WindType == IEC_EWM1 ) THEN
      Velocity = 0.8*p%IEC%VRef*( Ht/p%grid%HubHt )**p%met%PLExp                  ! [IEC 61400-1 6.3.2.1 (14), (15)]
      RETURN
   ELSEIF ( p%IEC%IEC_WindType == IEC_EWM100 ) THEN
      Velocity =    p%IEC%VRef*( Ht/p%grid%HubHt )**p%met%PLExp                   ! [API-IEC RECCOMENDATAION]  ADDED BY YGUO  !bjj: this is the same as IEC_EWM50, but we should check that IEC_EWM100 is used in ALL the same places IEC_EWM50 is
      RETURN
   ENDIF


   SELECT CASE ( TRIM(p%met%WindProfileType) )

      CASE ( 'JET' )

         tmpHt(1) = Ht
         CALL ChebyshevVals( p%met%ChebyCoef_WS, tmpHt(1:1), tmpWS(1:1), MinZ, MaxZ, ErrStat, ErrMsg ) ! We originally calculated the coeffs for 3-500 m in height
         Velocity = tmpWS(1)

      CASE ( 'LOG' ) !Panofsky, H.A.; Dutton, J.A. (1984). Atmospheric Turbulence: Models and Methods for Engineering Applications. New York: Wiley-Interscience; 397 pp.

         Velocity = getLogWindSpeed(Ht, z_Ref, U_Ref, p%met%ZL, p%met%Z0)

      CASE ( 'H2L' )
               ! Calculate the windspeed.
               ! z_Ref and U_Ref both get modified consistently, therefore z_Ref is used instead of RefHt.
               Velocity = LOG( Ht/z_Ref ) * p%met%Ustar / 0.41_ReKi + U_Ref


      CASE ( 'PL' )  ! POWER LAW, commented by Y. Guo on April 16 2013

         IF ( z_Ref > 0.0 .AND. Ht > 0.0 ) THEN
            Velocity = U_Ref*( Ht/z_Ref )**p%met%PLExp      ! [IEC 61400-1 6.3.1.2 (10)]
         ELSE
            Velocity = 0.0
         ENDIF

      CASE ( 'USR' )

         IF ( Ht <= p%met%USR_Z(1) ) THEN
            Velocity = p%met%USR_U(1)
         ELSEIF ( Ht >= p%met%USR_Z(p%met%NumUSRz) ) THEN
            Velocity = p%met%USR_U(p%met%NumUSRz)
         ELSE
            ! Find the two points between which the height lies

            DO I=2,p%met%NumUSRz
               IF ( Ht <= p%met%USR_Z(I) ) THEN
                  Indx = I-1

                  ! Let's just do a linear interpolation for now
                  Velocity = (Ht - p%met%USR_Z(Indx)) * ( p%met%USR_U(Indx) - p%met%USR_U(I) ) / ( p%met%USR_Z(Indx) - p%met%USR_Z(I) ) + p%met%USR_U(Indx)
                  EXIT
               ENDIF
            ENDDO

         ENDIF

      CASE ( 'TS' )

         Velocity =  getTimeSeriesWindSpeed(p, Ht)
 
         
      CASE ( 'API' ) 

!MLB: We can exclude this logic by forcing the user to enter the 1-hour mean wind speed.
!     If we add the API stuff to the main version of TurbSim, we may want to eliminate that requirement, but we will have to
!     add a new input parameter saying how long a period was used to calculate the reference wind speed.

!           CALL ROOT_SEARCHING(X0,X,U_Ref,z_Ref,p%grid%HubHt)  !URef
           !CALL Root_Searching(X0,X,42.5,10.0,10.0)  !U_Ref, USED TO DEBUG THE CODE
!           U0_1HR=X  ! This is the wind speed at 10 m height within 1-hr window

!            CALL WrScr ('Calling to API wind profile')
!           TEMP_1=0.0573*(1.0+0.15*U0_1HR)**0.5
!           TEMP_2=0.06*(1+0.043*U0_1HR)*(Ht/10.0)**(-0.22)
!           Velocity = U0_1HR*(1.0+TEMP_1*LOG( Ht / 10.0) )*( 1.0-0.41*TEMP_2*LOG(600.0/3600.0))
!           Velocity = U0_1HR*( 1.0 + 0.0573*SQRT( 1.0 + 0.15*U0_1HR )*LOG( Ht/z_Ref) )
!MLB: This assumes that the reference wind speed entered by the user is the 1-hour average wind speed at the input reference height.
           Velocity = p%met%URef*( 1.0 + 0.0573*SQRT( 1.0 + 0.15*p%met%URef )*LOG( Ht/p%met%RefHt) )

!            CALL WrScr ('API wind profile generated')

      CASE DEFAULT   ! This is how it worked before

         IF ( Ht == z_Ref ) THEN
            Velocity = U_Ref
         ELSEIF ( ABS( Ht-z_Ref ) <= 0.5*p%grid%RotorDiameter ) THEN
            Velocity = U_Ref*( Ht/z_Ref )**p%met%PLExp                ! [IEC 61400-1 6.3.1.2 (10)]
         ELSEIF ( Ht > 0.0 .AND. z_Ref > 0.0 .AND. .NOT. EqualRealNos(z_Ref, p%met%Z0) ) THEN !Check that we don't have an invalid domain
            Velocity = U_Ref*LOG( Ht/p%met%Z0 )/LOG( z_Ref/p%met%Z0 )
         ELSE
            Velocity = 0.0
         ENDIF

   END SELECT


RETURN
END SUBROUTINE getVelocity

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the wind speed at Ht by linearly interpolating the mean wind speed at the points from the user-input
!! time series.
function getTimeSeriesWindSpeed(p, Ht)

   TYPE(TurbSim_ParameterType),INTENT(IN) :: p                         !< parameters
   REAL(ReKi),                 INTENT(IN) :: Ht                        !< height at which wind speed is requested [m]
   REAL(ReKi)                             :: getTimeSeriesWindSpeed    !< the calculated wind speed at Ht

   INTEGER(IntKi)                         :: IZ, IZm1
   
   
   
   IF ( Ht <= p%usr%pointzi(1) ) THEN
      getTimeSeriesWindSpeed = p%usr%meanU(1,1)
   ELSEIF ( Ht >= p%usr%pointzi(p%usr%NPoints) ) THEN
      getTimeSeriesWindSpeed = p%usr%meanU(p%usr%NPoints,1)
   ELSE
      ! Find the two points between which the height lies

      DO IZ=2,p%usr%NPoints
         IF ( Ht <= p%usr%pointzi(IZ) ) THEN
            IZm1 = IZ-1

            ! Let's just do a linear interpolation for now
            getTimeSeriesWindSpeed = (Ht - p%usr%pointzi(IZm1)) * ( p%usr%meanU(IZm1,1) - p%usr%meanU(IZ,1) ) / ( p%usr%pointzi(IZm1) - p%usr%pointzi(IZ) ) + p%usr%meanU(IZm1,1)
            EXIT
         ENDIF
      ENDDO

   ENDIF       

end function getTimeSeriesWindSpeed

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine calculates the wind speed at Ht assuming a logarithmic wind profile and inputs z_Ref, U_Ref, Z/L and Z0
!!
!!      U_{ref}*( LOG( Ht / Z0 ) - psiM )/( LOG( z_Ref / Z0 ) - psiM )
!! where
!!      psiM is a function of Z/L
!! In neutral conditions, psiM is 0 and we get the IEC log wind profile.
function getLogWindSpeed(Ht, z_Ref, U_Ref, ZL, Z0)

   REAL(ReKi),             INTENT(IN) :: Ht                 !< height at which wind speed is requested [m]
   REAL(ReKi),             INTENT(IN) :: z_Ref              !< height of the reference wind speed [m]
   REAL(ReKi),             INTENT(IN) :: U_Ref              !< reference wind speed [m/s]
   REAL(ReKi),             INTENT(IN) :: ZL                 !< a measure of stability [-]
   REAL(ReKi),             INTENT(IN) :: Z0                 !< surface roughness length [m]

   REAL(ReKi)                         :: getLogWindSpeed    !< the calculated wind speed at Ht
   
   
      ! local variables
   REAL(ReKi)                         :: psiM               ! The diabatic term for the log wind profile
   REAL(ReKi)                         :: tmp                ! A temporary variable for calculating psiM
   
   
      
      IF ( Ht > 0.0 .AND. z_Ref > 0.0 .AND. .NOT. EqualRealNos( z_Ref, Z0 ) ) THEN

         IF ( ZL >= 0 ) THEN !& ZL < 1
            psiM = -5.0*ZL
         ELSE
            tmp = (1.0 - 15.0*ZL)**0.25

            !psiM = -2.0*LOG( (1.0 + tmp)/2.0 ) - LOG( (1.0 + tmp*tmp)/2.0 ) + 2.0*ATAN( tmp ) - 0.5 * PI
            psiM = -LOG( 0.125 * ( (1.0 + tmp)**2 * (1.0 + tmp*tmp) ) ) + 2.0*ATAN( tmp ) - 0.5 * PI
            
            !bjj 11-may-2016: because of the negative sign in the equation below, I believe psiM needs to switch signs.
            ! if true, this has been implemented incorrectly for at least 15 years.
            psiM = -psiM
            
         ENDIF

!            IF ( p%met%Ustar > 0. ) THEN
!               getLogWindSpeed = ( p%met%UstarDiab / 0.4 ) * ( LOG( Ht / Z0 ) - psiM )
!            ELSE
            !In neutral conditions, psiM is 0 and we get the IEC log wind profile:
         getLogWindSpeed = U_Ref*( LOG( Ht / Z0 ) - psiM )/( LOG( z_Ref / Z0 ) - psiM )
!            ENDIF

      ELSE
         getLogWindSpeed = 0.0_ReKi
      ENDIF


end function getLogWindSpeed
!=======================================================================
SUBROUTINE get_coefs(JetHt,UH_coef,WD_coef)

      ! This subroutine just returns the coefficients that Neil calculated
      ! for getting the Chebyshev coefficients for jet wind profiles.

      ! The coefficients are
      ! Row 1 = Jet maximum wind speed coefficient
      ! Row 2 = Turbine layer Richardson number coefficient
      ! Row 3 = uStar over the rotor diameter coefficient
      ! Row 4 = constant coefficient
      ! Columns 1:11 = coefficients for 0-10th Chebyshev Basis Functions


   REAL(ReKi),INTENT(IN)     :: JetHt                                    ! The height of the jet
   REAL(ReKi),INTENT(OUT)    :: UH_coef(4,11)                            ! The coefficients for horizontal wind speed
   REAL(ReKi),INTENT(OUT)    :: WD_coef(4,11)                            ! The coefficients for (horizontal) wind direction

   INTEGER                   :: HtIndx

   HtIndx  = INT(JetHt - 50) / INT(20)
   HtIndx  = MIN( MAX( HtIndx, 1 ), 21 )

      ! The Horizontal Wind Speed coefficients
   SELECT CASE ( HtIndx )
      CASE ( 1 )     ! 70-90 m
         UH_coef(:, 1) = (/  0.856851,  7.51E-02,  1.39276,   0.894127 /)
         UH_coef(:, 2) = (/ -4.88E-02,  0.576344,  1.23582,   1.72687  /)
         UH_coef(:, 3) = (/  1.39E-02,  9.67E-02,  1.36737,  -0.723851 /)
         UH_coef(:, 4) = (/  0.100585,  0.234968, -1.06287,  -0.372353 /)
         UH_coef(:, 5) = (/ -7.69E-02, -0.154071, -0.301483,  0.150179 /)
         UH_coef(:, 6) = (/  8.53E-03,  0.104602, -0.382453,  0.520224 /)
         UH_coef(:, 7) = (/ -4.44E-03, -4.80E-02,  0.219135, -0.266775 /)
         UH_coef(:, 8) = (/  2.63E-02, -3.08E-02, -6.94E-02, -0.210521 /)
         UH_coef(:, 9) = (/ -2.01E-02, -5.61E-02,  0.220825,  0.179622 /)
         UH_coef(:,10) = (/  8.11E-03,  3.96E-02,  0.109793, -3.81E-02 /)
         UH_coef(:,11) = (/  4.99E-03,  5.00E-02, -0.124887, -0.11035  /)
      CASE ( 2 )     ! 90-110 m
         UH_coef(:, 1) = (/  0.741241, -0.122521,  0.875062,  1.43294  /)
         UH_coef(:, 2) = (/ -0.264131,  0.28827,   0.717571,  3.30541  /)
         UH_coef(:, 3) = (/ -5.92E-02,  3.86E-02,  1.09453,  -0.377399 /)
         UH_coef(:, 4) = (/  0.13792,   0.175628, -0.57163,  -0.539205 /)
         UH_coef(:, 5) = (/ -2.59E-02, -0.211126, -4.25E-02, -0.338308 /)
         UH_coef(:, 6) = (/ -1.02E-02,  0.153597, -0.197867,  0.570708 /)
         UH_coef(:, 7) = (/ -3.22E-02, -8.17E-02, -9.63E-02,  0.19095  /)
         UH_coef(:, 8) = (/  2.72E-02,  3.09E-02, -0.249399, -0.273684 /)
         UH_coef(:, 9) = (/ -1.60E-02,  8.88E-03,  0.132523,  9.58E-02 /)
         UH_coef(:,10) = (/ -5.29E-03,  2.98E-02,  0.205812,  9.27E-02 /)
         UH_coef(:,11) = (/  7.00E-03, -1.47E-02, -2.11E-02, -0.123083 /)
      CASE ( 3 )     ! 110-130 m
         UH_coef(:, 1) = (/  0.809492, -1.41752,  -0.817619,  1.64159  /)
         UH_coef(:, 2) = (/ -0.121866, -1.09012,  -2.60044,   3.63875  /)
         UH_coef(:, 3) = (/ -0.105142, -0.263657, -5.60E-02,  0.374811 /)
         UH_coef(:, 4) = (/  8.33E-02,  0.625103,  0.422112, -0.199598 /)
         UH_coef(:, 5) = (/ -1.69E-02, -7.09E-02,  1.76933,  -0.847721 /)
         UH_coef(:, 6) = (/  1.88E-02,  7.70E-02, -0.121062,  0.10533  /)
         UH_coef(:, 7) = (/ -3.15E-02,  2.50E-02, -7.39E-02,  0.299197 /)
         UH_coef(:, 8) = (/  3.48E-03,  4.25E-02, -6.52E-02, -4.29E-03 /)
         UH_coef(:, 9) = (/ -1.18E-02, -0.100754,  0.170602,  3.42E-02 /)
         UH_coef(:,10) = (/  2.09E-02,  3.36E-02, -0.104123, -8.49E-02 /)
         UH_coef(:,11) = (/ -2.91E-03, -3.52E-02, -0.258115,  4.81E-02 /)
      CASE ( 4 )     ! 130-150 m
         UH_coef(:, 1) = (/  0.694325, -0.463252,  2.11406,   1.28643  /)
         UH_coef(:, 2) = (/ -0.269118, -1.31381,   2.13374,   3.46187  /)
         UH_coef(:, 3) = (/ -8.40E-02, -5.97E-02,  2.09803,  -0.592335 /)
         UH_coef(:, 4) = (/  0.135657, -0.117732, -0.11134,  -0.28161  /)
         UH_coef(:, 5) = (/ -1.29E-02, -0.239685,  0.151264, -0.412806 /)
         UH_coef(:, 6) = (/  3.54E-02,  0.513824,  0.673662, -0.519536 /)
         UH_coef(:, 7) = (/ -1.55E-02,  7.49E-03,  0.393002,  2.07E-02 /)
         UH_coef(:, 8) = (/  2.37E-02,  0.225841,  3.84E-02, -0.202507 /)
         UH_coef(:, 9) = (/ -3.26E-02, -0.239615, -0.133893,  0.29135  /)
         UH_coef(:,10) = (/  1.52E-02,  7.15E-02,  0.25228,  -0.113016 /)
         UH_coef(:,11) = (/  7.19E-03,  9.79E-02,  0.252125, -0.173201 /)
      CASE ( 5 )     ! 150-170 m
         UH_coef(:, 1) = (/  0.909534,  0.581254, -2.90539,  -0.581377 /)
         UH_coef(:, 2) = (/  0.155834, -0.836954, -6.77075,   0.627044 /)
         UH_coef(:, 3) = (/ -8.99E-02, -5.28E-02, -2.0719,    2.44E-02 /)
         UH_coef(:, 4) = (/  7.01E-02, -0.152904, -0.348237,  0.460754 /)
         UH_coef(:, 5) = (/ -1.78E-02, -0.263166,  0.375798, -0.215738 /)
         UH_coef(:, 6) = (/  9.70E-03,  0.254932,  0.449286, -0.234    /)
         UH_coef(:, 7) = (/  7.46E-03, -0.304057, -0.122661, -7.14E-03 /)
         UH_coef(:, 8) = (/ -6.26E-03, -0.142341, -1.95E-02,  0.299841 /)
         UH_coef(:, 9) = (/ -2.59E-02,  0.174282,  0.193868, -5.81E-03 /)
         UH_coef(:,10) = (/  2.54E-03, -8.22E-02,  1.84E-02,  6.77E-02 /)
         UH_coef(:,11) = (/  5.77E-04, -5.43E-02, -7.69E-02,  2.96E-02 /)
      CASE ( 6 )     ! 170-190 m
         UH_coef(:, 1) = (/  0.885753, -1.15015,   0.155218, -0.707043 /)
         UH_coef(:, 2) = (/ -2.53E-02, -2.65126,   0.850151,  1.85279  /)
         UH_coef(:, 3) = (/ -7.23E-02, -0.399161,  0.142486, -0.917176 /)
         UH_coef(:, 4) = (/  3.78E-02,  0.178924,  0.227745,  0.528861 /)
         UH_coef(:, 5) = (/ -6.43E-03,  5.42E-02,  0.359052, -0.26111  /)
         UH_coef(:, 6) = (/  5.33E-02,  0.1546,   -0.335116, -0.602604 /)
         UH_coef(:, 7) = (/ -6.50E-03, -0.205907, -8.59E-02,  8.16E-02 /)
         UH_coef(:, 8) = (/  3.16E-02,  0.151199, -0.126411, -0.148609 /)
         UH_coef(:, 9) = (/ -3.95E-02,  0.127418,  0.158511,  0.20932  /)
         UH_coef(:,10) = (/ -2.53E-02, -5.32E-02,  0.36536,   0.214466 /)
         UH_coef(:,11) = (/  4.03E-03,  1.02E-02, -7.01E-03, -4.32E-02 /)
      CASE ( 7 )     ! 190-210 m
         UH_coef(:, 1) = (/  0.735269, -1.48574,   0.983734,  0.887351 /)
         UH_coef(:, 2) = (/  0.233065, -0.850536, -1.17754,  -0.880493 /)
         UH_coef(:, 3) = (/ -0.172346, -0.862128,  1.20075,   3.48E-02 /)
         UH_coef(:, 4) = (/  8.04E-02,  5.24E-02, -0.916548,  0.247144 /)
         UH_coef(:, 5) = (/  2.88E-02,  0.112064,  1.51E-04, -0.466186 /)
         UH_coef(:, 6) = (/ -2.75E-02, -9.01E-02, -0.321617,  0.379162 /)
         UH_coef(:, 7) = (/ -1.08E-02, -0.161368, -2.51E-04, -1.33E-02 /)
         UH_coef(:, 8) = (/  5.09E-02,  0.228507,  0.195942, -0.45807  /)
         UH_coef(:, 9) = (/ -1.98E-02, -7.23E-02,  6.66E-02,  0.133182 /)
         UH_coef(:,10) = (/ -5.57E-03, -5.31E-02,  2.44E-02,  5.60E-02 /)
         UH_coef(:,11) = (/  3.71E-03, -1.63E-02, -5.44E-02, -1.40E-02 /)
      CASE ( 8 )     ! 210-230 m
         UH_coef(:, 1) = (/  0.723721, -0.691359, -0.147971,  1.16041  /)
         UH_coef(:, 2) = (/  0.18799,   0.370199,  0.354538, -0.494962 /)
         UH_coef(:, 3) = (/ -0.204727, -0.166723,  0.682431,  0.367566 /)
         UH_coef(:, 4) = (/  1.40E-02,  0.334677,  0.169944,  0.494211 /)
         UH_coef(:, 5) = (/  3.84E-02,  0.258361,  0.389453, -0.625709 /)
         UH_coef(:, 6) = (/ -6.62E-03, -2.19E-02, -0.606278,  0.205521 /)
         UH_coef(:, 7) = (/ -2.54E-02, -0.17744,   7.49E-02,  7.61E-02 /)
         UH_coef(:, 8) = (/  5.03E-02,  7.97E-02, -9.98E-02, -0.312218 /)
         UH_coef(:, 9) = (/ -2.25E-02,  2.20E-02,  0.263227,  0.123311 /)
         UH_coef(:,10) = (/ -1.43E-02, -2.01E-02, -5.14E-02,  0.159391 /)
         UH_coef(:,11) = (/  2.64E-03,  3.46E-02, -0.12318,  -2.22E-02 /)
      CASE ( 9 )     ! 230-250 m
         UH_coef(:, 1) = (/  0.717665, -0.294178, -0.521541,  0.876418 /)
         UH_coef(:, 2) = (/  0.183182, -0.52658,  -1.34668,   0.414396 /)
         UH_coef(:, 3) = (/ -0.196162,  9.84E-02, -3.83E-02,  0.156018 /)
         UH_coef(:, 4) = (/  2.92E-02, -0.362193, -0.658593,  0.521854 /)
         UH_coef(:, 5) = (/  3.37E-02,  0.108203,  0.318667, -0.375309 /)
         UH_coef(:, 6) = (/ -8.24E-03,  0.128457, -0.149225,  0.1621   /)
         UH_coef(:, 7) = (/ -3.06E-02, -0.210106,  4.55E-02,  8.42E-02 /)
         UH_coef(:, 8) = (/  3.02E-02,  0.184626,  9.46E-02, -0.215191 /)
         UH_coef(:, 9) = (/  7.03E-03,  2.49E-02,  3.13E-02, -9.70E-02 /)
         UH_coef(:,10) = (/ -3.06E-03, -4.82E-02, -9.70E-02,  5.82E-02 /)
         UH_coef(:,11) = (/ -9.57E-03, -3.93E-02, -0.125623,  0.112639 /)
      CASE ( 10 )    ! 250-270 m
         UH_coef(:, 1) = (/  0.786229, -0.164848,  0.244948, -0.126263 /)
         UH_coef(:, 2) = (/  0.15218,  -0.153233, -0.558524,  0.84425  /)
         UH_coef(:, 3) = (/ -0.130716, -0.217411,  0.13439,  -0.536893 /)
         UH_coef(:, 4) = (/  1.70E-03,  5.49E-02,  0.551012,  0.335778 /)
         UH_coef(:, 5) = (/  2.47E-02,  2.82E-02,  0.290918, -0.223416 /)
         UH_coef(:, 6) = (/  1.48E-02,  5.94E-02, -0.277959,  3.91E-02 /)
         UH_coef(:, 7) = (/ -4.43E-02,  6.99E-03,  0.302386,  0.123719 /)
         UH_coef(:, 8) = (/  2.07E-02,  4.05E-02, -0.256155, -5.84E-02 /)
         UH_coef(:, 9) = (/  4.51E-03, -4.37E-02, -0.111911, -9.20E-03 /)
         UH_coef(:,10) = (/  4.05E-03, -6.90E-03,  0.14697,  -7.03E-02 /)
         UH_coef(:,11) = (/ -6.68E-03,  1.53E-02, -2.55E-02,  4.97E-02 /)
      CASE ( 11 )    ! 270-290 m
         UH_coef(:, 1) = (/  0.715734, -0.772062, -0.556396,  1.02929  /)
         UH_coef(:, 2) = (/  0.322509, -0.465616, -0.671711, -1.2413   /)
         UH_coef(:, 3) = (/ -0.166728, -0.281268,  0.924893, -0.282907 /)
         UH_coef(:, 4) = (/  1.27E-02, -0.342767, -1.10823,   0.516431 /)
         UH_coef(:, 5) = (/  3.80E-02,  5.35E-03,  0.833719, -0.510102 /)
         UH_coef(:, 6) = (/  1.97E-02,  0.279705, -0.179026, -4.36E-02 /)
         UH_coef(:, 7) = (/ -4.74E-02, -0.227673,  9.00E-02,  0.341958 /)
         UH_coef(:, 8) = (/  8.99E-03, -1.92E-02, -0.433969,  5.90E-02 /)
         UH_coef(:, 9) = (/  4.34E-03,  8.12E-02,  0.25764,  -0.148492 /)
         UH_coef(:,10) = (/  1.03E-02,  3.24E-02,  0.141971, -0.105207 /)
         UH_coef(:,11) = (/ -4.84E-03, -1.99E-02,  7.33E-02,  2.84E-02 /)
      CASE ( 12 )    ! 290-310 m
         UH_coef(:, 1) = (/  0.723348, -0.289581, -1.10618,   0.970713 /)
         UH_coef(:, 2) = (/  0.283383,  1.12986,  -0.152861, -0.653269 /)
         UH_coef(:, 3) = (/ -0.16513,   0.295047,  0.245326, -7.06E-02 /)
         UH_coef(:, 4) = (/  8.55E-03,  9.38E-02, -0.826824,  0.283436 /)
         UH_coef(:, 5) = (/  3.45E-02,  0.364581,  0.566317, -0.521081 /)
         UH_coef(:, 6) = (/  2.83E-02,  0.107252, -0.124867, -4.80E-02 /)
         UH_coef(:, 7) = (/ -3.57E-02, -0.230151, -6.88E-02,  0.231208 /)
         UH_coef(:, 8) = (/  5.62E-04,  1.40E-02, -0.334942,  0.121313 /)
         UH_coef(:, 9) = (/ -6.35E-03, -6.19E-02,  0.139396,  2.77E-02 /)
         UH_coef(:,10) = (/  1.14E-02, -2.67E-02,  0.24201,  -0.127337 /)
         UH_coef(:,11) = (/  1.71E-04, -6.37E-04,  4.39E-02, -5.61E-03 /)
      CASE ( 13 )    ! 310-330 m
         UH_coef(:, 1) = (/  0.736987, -0.103727,  9.95E-02,  0.343208 /)
         UH_coef(:, 2) = (/  0.28285,   0.370583,  1.17749,  -0.490259 /)
         UH_coef(:, 3) = (/ -0.130451, -0.557928, -0.272771, -0.230816 /)
         UH_coef(:, 4) = (/ -1.83E-02,  1.00E-01, -0.367321,  0.486971 /)
         UH_coef(:, 5) = (/  2.66E-02, -0.149206,  0.365342, -0.318809 /)
         UH_coef(:, 6) = (/  4.16E-02,  3.60E-02, -0.801161,  6.00E-06 /)
         UH_coef(:, 7) = (/ -2.36E-02,  1.96E-04,  0.340449,  2.72E-02 /)
         UH_coef(:, 8) = (/  1.30E-03,  0.214384,  0.125371, -8.47E-02 /)
         UH_coef(:, 9) = (/ -1.23E-02,  4.75E-02,  0.182118,  1.78E-02 /)
         UH_coef(:,10) = (/  4.63E-03, -0.1309,   -0.130584,  2.35E-02 /)
         UH_coef(:,11) = (/  9.03E-04, -6.18E-02, -7.85E-03,  1.17E-02 /)
      CASE ( 14 )    ! 330-350 m
         UH_coef(:, 1) = (/  0.706488, -1.21766,   1.08617,   0.674247 /)
         UH_coef(:, 2) = (/  0.341777,  2.27476,   3.81434,  -2.32363  /)
         UH_coef(:, 3) = (/ -0.112822,  7.53E-02,  0.221349, -0.700428 /)
         UH_coef(:, 4) = (/ -1.99E-02, -1.95E-02,  0.947788,  4.68E-02 /)
         UH_coef(:, 5) = (/  3.08E-02,  0.334947,  0.10847,  -0.534662 /)
         UH_coef(:, 6) = (/  5.21E-02,  0.349056, -1.14517,  -0.147474 /)
         UH_coef(:, 7) = (/ -1.67E-02, -0.143994, -0.409398,  0.228081 /)
         UH_coef(:, 8) = (/ -1.75E-03, -0.115198,  3.23E-03,  0.100094 /)
         UH_coef(:, 9) = (/ -2.30E-02, -5.63E-02,  0.168561,  0.159537 /)
         UH_coef(:,10) = (/ -6.41E-03, -8.48E-02,  0.135087,  8.81E-02 /)
         UH_coef(:,11) = (/  1.13E-03,  2.07E-02,  9.18E-02, -3.77E-02 /)
      CASE ( 15 )    ! 350-370 m
         UH_coef(:, 1) = (/  0.721629, -0.941544,  0.923908,  0.543678 /)
         UH_coef(:, 2) = (/  0.346956, -0.281582, -2.32358,  -0.244435 /)
         UH_coef(:, 3) = (/ -0.109484,  0.275053,  0.86928,  -0.771081 /)
         UH_coef(:, 4) = (/ -3.96E-02, -0.790621, -8.84E-02,  0.723378 /)
         UH_coef(:, 5) = (/  1.59E-02, -0.394222, -0.479505, -8.67E-02 /)
         UH_coef(:, 6) = (/  2.68E-02,  0.466895,  0.522378, -0.263669 /)
         UH_coef(:, 7) = (/ -9.57E-03, -8.52E-02,  1.11E-02,  3.20E-02 /)
         UH_coef(:, 8) = (/  3.46E-04, -5.34E-02,  0.15998,   0.108225 /)
         UH_coef(:, 9) = (/ -1.10E-02, -0.116864, -6.06E-02,  6.09E-02 /)
         UH_coef(:,10) = (/ -2.93E-03,  2.72E-02,  5.08E-02,  7.50E-03 /)
         UH_coef(:,11) = (/ -2.04E-03, -2.07E-02, -3.07E-02,  3.58E-02 /)
      CASE ( 16 )    ! 370-390 m
         UH_coef(:, 1) = (/  0.732127, -2.66819,  -7.94E-02,  0.676096 /)
         UH_coef(:, 2) = (/  0.285167,  3.89442,  -0.917426,  0.104248 /)
         UH_coef(:, 3) = (/ -8.38E-02,  0.235268, -2.19E-03, -0.914663 /)
         UH_coef(:, 4) = (/ -3.98E-02, -0.858603, -0.538194,  0.843739 /)
         UH_coef(:, 5) = (/ -1.64E-02,  0.287007, -5.39E-02,  0.108834 /)
         UH_coef(:, 6) = (/  3.31E-02,  0.218726,  0.175636, -0.329844 /)
         UH_coef(:, 7) = (/  3.10E-05, -6.89E-02,  3.76E-02, -4.73E-02 /)
         UH_coef(:, 8) = (/  1.06E-02, -5.03E-02,  1.99E-02,  3.74E-02 /)
         UH_coef(:, 9) = (/ -1.05E-02,  9.92E-02,  0.11293,   2.26E-02 /)
         UH_coef(:,10) = (/ -2.99E-03, -0.106831,  0.122628,  1.83E-02 /)
         UH_coef(:,11) = (/ -7.32E-03,  3.52E-02, -3.36E-02,  8.59E-02 /)
      CASE ( 17 )    ! 390-410 m
         UH_coef(:, 1) = (/  0.707698,  0.119876,  0.427545,  0.2468   /)
         UH_coef(:, 2) = (/  0.307273,  0.428003, -3.09224,   1.01117  /)
         UH_coef(:, 3) = (/ -7.33E-02,  0.51572,  -0.229086, -0.792402 /)
         UH_coef(:, 4) = (/ -4.73E-02,  8.49E-02, -0.52415,   0.571084 /)
         UH_coef(:, 5) = (/ -2.83E-02,  0.165455, -0.691726,  0.349932 /)
         UH_coef(:, 6) = (/  2.17E-02,  0.258434,  0.170597, -0.236707 /)
         UH_coef(:, 7) = (/ -4.59E-03, -0.130722,  0.182955, -3.40E-02 /)
         UH_coef(:, 8) = (/  1.82E-02,  9.79E-02,  0.189511, -0.158597 /)
         UH_coef(:, 9) = (/ -7.84E-04, -2.50E-02,  0.137171, -5.77E-02 /)
         UH_coef(:,10) = (/ -2.91E-03, -4.84E-02,  0.168698,  8.22E-03 /)
         UH_coef(:,11) = (/ -4.67E-03,  1.75E-03,  1.80E-02,  4.41E-02 /)
      CASE ( 18 )    ! 410-430 m
         UH_coef(:, 1) = (/  0.688761, -0.7286,   -1.55711,   1.27145  /)
         UH_coef(:, 2) = (/  0.300421,  0.633115,  0.881706, -8.38E-03 /)
         UH_coef(:, 3) = (/ -6.81E-02,  0.210301,  0.610772, -0.714435 /)
         UH_coef(:, 4) = (/ -5.93E-02, -0.373997, -0.593894,  1.01556  /)
         UH_coef(:, 5) = (/ -4.26E-02, -2.45E-02, -0.400705,  0.399717 /)
         UH_coef(:, 6) = (/  1.39E-02,  6.09E-02, -0.161239, -3.06E-02 /)
         UH_coef(:, 7) = (/ -4.41E-03, -1.98E-02,  0.293288, -0.110401 /)
         UH_coef(:, 8) = (/  1.42E-02,  8.22E-02, -1.50E-02, -1.54E-02 /)
         UH_coef(:, 9) = (/  6.30E-03, -1.50E-02, -7.57E-02, -7.10E-02 /)
         UH_coef(:,10) = (/  2.19E-03, -2.59E-02,  8.53E-02, -2.29E-02 /)
         UH_coef(:,11) = (/ -2.76E-03,  1.68E-02, -8.77E-02,  3.27E-02 /)
      CASE ( 19 )    ! 430-450 m
         UH_coef(:, 1) = (/  0.659495, -0.22327,  -1.75403,   1.65777  /)
         UH_coef(:, 2) = (/  0.384097,  1.06351,   2.53779,  -1.63428  /)
         UH_coef(:, 3) = (/ -2.42E-02,  0.113735, -1.42805,  -0.690773 /)
         UH_coef(:, 4) = (/ -3.30E-02,  8.60E-02, -1.00836,   0.764307 /)
         UH_coef(:, 5) = (/ -2.76E-02,  0.297567,  0.697445, -0.187071 /)
         UH_coef(:, 6) = (/  1.21E-02,  0.212621, -0.570822,  1.23E-02 /)
         UH_coef(:, 7) = (/ -2.22E-02,  0.166286,  0.50751,   1.87E-02 /)
         UH_coef(:, 8) = (/  1.52E-02,  5.81E-02, -0.256912, -5.10E-02 /)
         UH_coef(:, 9) = (/  2.11E-03, -1.45E-02, -8.94E-02, -2.00E-02 /)
         UH_coef(:,10) = (/  3.06E-03,  1.60E-02,  7.45E-02, -3.77E-02 /)
         UH_coef(:,11) = (/ -1.84E-04, -1.56E-02, -6.25E-02,  1.57E-02 /)
      CASE ( 20 )    ! 450-470 m
         UH_coef(:, 1) = (/  0.64099,  -2.02496,   0.427597,  1.52166  /)
         UH_coef(:, 2) = (/  0.391609,  2.03441,  -0.122486, -1.03579  /)
         UH_coef(:, 3) = (/  8.28E-03,  0.5942,   -0.42469,  -1.35655  /)
         UH_coef(:, 4) = (/ -2.54E-02, -0.826812, -0.812187,  0.911776 /)
         UH_coef(:, 5) = (/ -2.77E-02, -9.73E-03,  0.315974,  2.34E-02 /)
         UH_coef(:, 6) = (/  1.37E-02,  0.365984,  0.141952, -0.299349 /)
         UH_coef(:, 7) = (/ -1.95E-02, -0.406182,  2.32E-02,  0.184752 /)
         UH_coef(:, 8) = (/  7.34E-03,  8.54E-02, -0.255458,  7.08E-02 /)
         UH_coef(:, 9) = (/  1.54E-03,  5.82E-02, -5.72E-02, -6.37E-02 /)
         UH_coef(:,10) = (/  5.11E-03, -6.11E-02, -7.04E-03, -3.64E-02 /)
         UH_coef(:,11) = (/  1.97E-03, -1.09E-02, -8.18E-02, -6.03E-03 /)
      CASE ( 21 )    ! 470-490 m
         UH_coef(:, 1) = (/  0.547127, -0.327778,  2.00666,   2.67869  /)
         UH_coef(:, 2) = (/  0.427112,  8.56E-02, -1.61197,  -1.17989  /)
         UH_coef(:, 3) = (/  6.23E-02,  0.760714, -0.659927, -2.30882  /)
         UH_coef(:, 4) = (/ -4.04E-02, -0.873328, -0.118326,  1.19626  /)
         UH_coef(:, 5) = (/ -4.85E-03,  0.130813, -0.169613, -0.181674 /)
         UH_coef(:, 6) = (/  4.82E-03,  0.289038,  7.34E-02,  6.45E-03 /)
         UH_coef(:, 7) = (/ -2.49E-02, -0.375342,  0.15139,   0.208253 /)
         UH_coef(:, 8) = (/  9.48E-04,  5.23E-02, -0.213227,  0.137941 /)
         UH_coef(:, 9) = (/ -9.18E-03,  3.91E-02,  7.26E-02,  4.73E-02 /)
         UH_coef(:,10) = (/ -6.00E-05,  1.03E-02,  7.46E-03,  1.86E-02 /)
         UH_coef(:,11) = (/ -2.21E-03, -9.70E-05, -7.13E-02,  4.29E-02 /)
      CASE DEFAULT
         CALL ProgAbort ('Error getting UH coefficients' )
   END SELECT

   SELECT CASE ( HtIndx )
      CASE ( 1 )     ! 70-90 m
         WD_coef(:, 1) = (/  5.07735,   96.4785,   18.8465,   110.986    /)
         WD_coef(:, 2) = (/  0.75209,  -16.5103,  -25.9592,     9.05636  /)
         WD_coef(:, 3) = (/ -1.50806,    1.69319,  -7.7859,    13.3041   /)
         WD_coef(:, 4) = (/  1.11287,    3.711,    13.1084,   -11.9491   /)
         WD_coef(:, 5) = (/ -0.987363,  -2.93059,  -4.75454,    9.04282  /)
         WD_coef(:, 6) = (/  0.65727,    0.560223, -0.541911,  -5.33397  /)
         WD_coef(:, 7) = (/ -0.493572,  -0.455574,  2.03972,    3.53745  /)
         WD_coef(:, 8) = (/  0.244207,   0.390402,  1.5338,    -1.9793   /)
         WD_coef(:, 9) = (/ -1.26E-02,   0.19732,  -2.70454,    0.179412 /)
         WD_coef(:,10) = (/  9.13E-04,   9.65E-02,  0.304467,   4.79E-02 /)
         WD_coef(:,11) = (/ -7.71E-02,  -0.11096,   0.51028,    0.585717 /)
      CASE ( 2 )     ! 90-110 m
         WD_coef(:, 1) = (/  2.98622,   87.1045,   41.7453,   124.301    /)
         WD_coef(:, 2) = (/  0.241282, -10.9238,  -31.5696,    11.0764   /)
         WD_coef(:, 3) = (/ -0.380786,  -1.71395,  -8.35561,    3.68007  /)
         WD_coef(:, 4) = (/  0.287014,   6.76407,  17.1736,    -7.4345   /)
         WD_coef(:, 5) = (/ -0.682991,  -5.48805, -12.7947,    10.9313   /)
         WD_coef(:, 6) = (/  0.415999,   2.36938,   4.47285,   -5.47595  /)
         WD_coef(:, 7) = (/ -0.184533,  -7.04E-02,  0.81309,    1.06891  /)
         WD_coef(:, 8) = (/  0.152381,  -0.344921,  3.40496,   -1.81465  /)
         WD_coef(:, 9) = (/ -0.113556,  -1.02575,  -5.54619,    2.51668  /)
         WD_coef(:,10) = (/  3.87E-02,   1.0794,    0.98668,   -0.942351 /)
         WD_coef(:,11) = (/  7.37E-02,  -0.284347,  1.12315,   -1.04163  /)
      CASE ( 3 )     ! 110-130 m
         WD_coef(:, 1) = (/ -10.8064,   63.1523,   18.7751,   255.252    /)
         WD_coef(:, 2) = (/  1.89875,  -15.7662,  -27.2545,    -5.90699  /)
         WD_coef(:, 3) = (/ -1.81141,   -7.58E-03,  4.49E-02,  19.4007   /)
         WD_coef(:, 4) = (/ -0.420216,   4.54261,  16.6642,    -1.5632   /)
         WD_coef(:, 5) = (/  3.09E-02,   0.162346, -5.68196,    1.70168  /)
         WD_coef(:, 6) = (/  0.372585,  -0.888944, -0.400871,  -3.98736  /)
         WD_coef(:, 7) = (/  0.137532,  -1.86E-02, -1.97659,   -1.07897  /)
         WD_coef(:, 8) = (/  7.11E-02,   0.275322,  2.06716,   -0.99703  /)
         WD_coef(:, 9) = (/ -0.142081,   0.690143,  1.74256,    0.963168 /)
         WD_coef(:,10) = (/ -0.225792,  -0.215169,  0.660299,   1.89319  /)
         WD_coef(:,11) = (/  1.91E-02,  -0.23,     -1.69222,    0.190668 /)
      CASE ( 4 )     ! 130-150 m
         WD_coef(:, 1) = (/  0.270461, 107.786,   140.705,    143.549    /)
         WD_coef(:, 2) = (/  2.46519,   25.9261,   54.6629,   -43.2182   /)
         WD_coef(:, 3) = (/ -1.11746,   -4.09287,  -5.71316,   16.4144   /)
         WD_coef(:, 4) = (/ -0.104557,   2.88836,  14.657,     -5.58632  /)
         WD_coef(:, 5) = (/  1.4104,    -0.862421,  1.88282,  -13.3856   /)
         WD_coef(:, 6) = (/ -0.994103,   6.07897,   6.16378,    6.53327  /)
         WD_coef(:, 7) = (/  0.440338,  -7.14173, -12.2957,     0.653282 /)
         WD_coef(:, 8) = (/ -0.705677,   2.13336,   2.39331,    5.62277  /)
         WD_coef(:, 9) = (/  0.398742,  -3.5049,   -3.97854,   -1.68531  /)
         WD_coef(:,10) = (/ -7.72E-02,   2.14124,   3.42657,   -0.982025 /)
         WD_coef(:,11) = (/  0.120525,  -1.80518,  -3.44124,    0.391772 /)
      CASE ( 5 )     ! 150-170 m
         WD_coef(:, 1) = (/  10.3894,  203.711,    87.9736,     0.818669 /)
         WD_coef(:, 2) = (/  4.15105,   37.734,    56.1061,   -72.0928   /)
         WD_coef(:, 3) = (/ -1.60031,   -6.42686,   2.99983,   21.7355   /)
         WD_coef(:, 4) = (/  0.162421, -22.7335,    4.23498,    0.433394 /)
         WD_coef(:, 5) = (/ -1.00817,   -1.82237, -17.2291,    18.8346   /)
         WD_coef(:, 6) = (/  0.591051,   5.30019,  22.1782,   -15.2786   /)
         WD_coef(:, 7) = (/ -0.350898,  -1.35238, -14.9057,     9.09022  /)
         WD_coef(:, 8) = (/  0.512704,   5.33682,  12.0501,   -11.3284   /)
         WD_coef(:, 9) = (/ -0.294613,  -6.61282, -13.756,      9.48747  /)
         WD_coef(:,10) = (/  0.180824,   6.67558,   8.1748,    -6.39538  /)
         WD_coef(:,11) = (/ -0.168678,  -3.5973,   -2.92266,    3.62255  /)
      CASE ( 6 )     ! 170-190 m
         WD_coef(:, 1) = (/ -3.05838,   92.242,    -6.17694,  218.678    /)
         WD_coef(:, 2) = (/ -1.19176,   10.9436,    5.33317,   23.6574   /)
         WD_coef(:, 3) = (/  0.396791,   5.36609,  14.86,     -12.1807   /)
         WD_coef(:, 4) = (/ -0.260044,  -3.3155,   -1.83325,    3.07872  /)
         WD_coef(:, 5) = (/  0.147588,   3.54423,   2.61624,   -2.87076  /)
         WD_coef(:, 6) = (/ -3.09E-02,  -0.298005, -3.99378,    2.512    /)
         WD_coef(:, 7) = (/  3.52E-02,   0.476622,  0.917889,  -1.19482  /)
         WD_coef(:, 8) = (/ -0.10397,   -3.13393,  -1.34654,    2.38467  /)
         WD_coef(:, 9) = (/  0.111959,   0.768005,  1.09164,   -1.84864  /)
         WD_coef(:,10) = (/ -5.32E-02,  -0.753046,  0.517477,   0.77376  /)
         WD_coef(:,11) = (/  2.36E-02,  -0.255733, -0.765475,  -0.183366 /)
      CASE ( 7 )     ! 190-210 m
         WD_coef(:, 1) = (/  2.63747,   48.8574, -148.839,    198.635    /)
         WD_coef(:, 2) = (/  0.276349,   8.15568,  11.5466,     4.89475  /)
         WD_coef(:, 3) = (/ -0.161153,  -3.92434,  15.2465,    -2.75263  /)
         WD_coef(:, 4) = (/ -0.215546,  -6.05707,  -0.221136,   2.96778  /)
         WD_coef(:, 5) = (/ -0.174687,   0.722833,  2.58751,    1.43519  /)
         WD_coef(:, 6) = (/ -3.24E-03,   0.841219,  2.36677,   -0.541046 /)
         WD_coef(:, 7) = (/ -0.14379,   -0.422125,  6.03272,   -3.55E-02 /)
         WD_coef(:, 8) = (/  4.94E-02,  -0.165447, -1.64947,   -0.118004 /)
         WD_coef(:, 9) = (/  6.88E-03,   0.618011,  0.600728,  -0.312735 /)
         WD_coef(:,10) = (/ -2.96E-02,  -0.102388, -0.423526,   0.526055 /)
         WD_coef(:,11) = (/  3.77E-03,  -0.79762,  -1.48591,    0.487559 /)
      CASE ( 8 )     ! 210-230 m
         WD_coef(:, 1) = (/  1.25931,   81.7121,  -72.2497,   192.288    /)
         WD_coef(:, 2) = (/ -0.421425,   0.812039, 26.4136,    12.7087   /)
         WD_coef(:, 3) = (/ -0.477334,  -0.804493, 10.2938,     2.63738  /)
         WD_coef(:, 4) = (/  0.27025,   -1.48414,   6.44E-02,  -3.62925  /)
         WD_coef(:, 5) = (/ -0.206555,   2.60212,   4.78E-03,   1.41829  /)
         WD_coef(:, 6) = (/  0.199714,  -0.145286, -1.43609,   -1.0421   /)
         WD_coef(:, 7) = (/ -8.81E-02,  -1.11826,   0.562309,   0.568182 /)
         WD_coef(:, 8) = (/  4.38E-02,  -0.94946,  -1.20199,    0.184361 /)
         WD_coef(:, 9) = (/ -5.13E-02,  -0.157795, -0.596316,   0.747777 /)
         WD_coef(:,10) = (/  5.03E-02,   6.23E-02, -0.821348,  -0.411198 /)
         WD_coef(:,11) = (/ -2.45E-02,   3.66E-03,  0.61934,    0.147334 /)
      CASE ( 9 )     ! 230-250 m
         WD_coef(:, 1) = (/  4.99773,   45.439,   -22.9981,   142.166    /)
         WD_coef(:, 2) = (/  1.34923,   -0.690733,  1.11037,   -7.00256  /)
         WD_coef(:, 3) = (/ -4.58E-02,  -1.48399,   3.15438,   -1.20619  /)
         WD_coef(:, 4) = (/ -5.86E-02,  -0.324401, -0.520264,   0.827308 /)
         WD_coef(:, 5) = (/  6.67E-02,   1.95293,  -1.46579,   -1.66186  /)
         WD_coef(:, 6) = (/  2.23E-02,   1.10257,   1.61038,   -0.14154  /)
         WD_coef(:, 7) = (/  4.83E-02,  -0.46633,   0.318096,  -1.22718  /)
         WD_coef(:, 8) = (/ -3.56E-02,  -0.905797, -0.659337,   1.10221  /)
         WD_coef(:, 9) = (/ -6.54E-04,   0.514329,  0.38488,   -0.221416 /)
         WD_coef(:,10) = (/  2.40E-03,  -0.307029, -0.455799,   0.167602 /)
         WD_coef(:,11) = (/  5.79E-03,  -0.3575,   -6.82E-02,  -1.79E-02 /)
      CASE ( 10 )    ! 250-270 m
         WD_coef(:, 1) = (/  2.87491,   81.7603,  -14.221,    143.973    /)
         WD_coef(:, 2) = (/  0.176626,   0.711168, 14.3778,     3.41781  /)
         WD_coef(:, 3) = (/ -0.112353,  -4.44334,   5.01439,   -0.539061 /)
         WD_coef(:, 4) = (/  0.135496,   0.868787, -2.54952,   -1.4882   /)
         WD_coef(:, 5) = (/ -5.87E-02,   7.34E-02,  0.618705,   0.341871 /)
         WD_coef(:, 6) = (/  4.36E-02,   1.16076,  -2.2411,     0.371484 /)
         WD_coef(:, 7) = (/ -4.21E-03,  -0.219162,  3.07613,   -1.48294  /)
         WD_coef(:, 8) = (/  2.91E-02,  -7.90E-02, -2.06058,    0.637811 /)
         WD_coef(:, 9) = (/  6.84E-04,   0.398542, -0.227958,  -0.195655 /)
         WD_coef(:,10) = (/ -1.33E-02,  -0.148014,  0.112677,   0.28039  /)
         WD_coef(:,11) = (/  4.56E-02,  -0.4372,   -1.05259,   -0.39506  /)
      CASE ( 11 )    ! 270-290 m
         WD_coef(:, 1) = (/ -3.74E-02,   5.72313, -25.8459,   204.708    /)
         WD_coef(:, 2) = (/  0.387587,   5.70337,  37.0722,    -5.10619  /)
         WD_coef(:, 3) = (/  0.130067,   8.86213,   7.6219,    -6.77984  /)
         WD_coef(:, 4) = (/ -1.83E-02,  -4.80402,   1.26728,    1.1988   /)
         WD_coef(:, 5) = (/ -0.125984,   5.69111,  -2.4798,     0.370193 /)
         WD_coef(:, 6) = (/  7.02E-02,  -4.02809,   0.545202,   0.396538 /)
         WD_coef(:, 7) = (/ -4.89E-02,   1.99119,  -7.47E-02,  -0.617665 /)
         WD_coef(:, 8) = (/  7.28E-02,  -1.94844,  -0.9012,     0.174322 /)
         WD_coef(:, 9) = (/ -2.75E-02,   0.875895,  8.29E-02,   1.47E-02 /)
         WD_coef(:,10) = (/ -4.90E-03,  -0.26505,   0.684299,  -0.101304 /)
         WD_coef(:,11) = (/ -2.46E-03,  -9.03E-02, -0.25124,    0.130552 /)
      CASE ( 12 )    ! 290-310 m
         WD_coef(:, 1) = (/  4.48806,  101.681,   -24.2152,   108.849    /)
         WD_coef(:, 2) = (/  1.12228,  -11.8153,   -5.83094,   -3.59506  /)
         WD_coef(:, 3) = (/  0.152934,   0.610899, 10.1148,    -6.59595  /)
         WD_coef(:, 4) = (/  6.76E-02,   1.44362,  -8.36227,    1.70741  /)
         WD_coef(:, 5) = (/ -8.86E-02,   1.22016,   4.89384,   -1.422    /)
         WD_coef(:, 6) = (/  1.14E-02,  -0.801065, -4.6529,     2.29577  /)
         WD_coef(:, 7) = (/ -5.68E-03,  -0.156515,  3.48364,   -1.85745  /)
         WD_coef(:, 8) = (/  3.21E-02,   0.643855, -1.80571,    0.499593 /)
         WD_coef(:, 9) = (/ -5.96E-03,  -0.645,     1.0105,    -0.256849 /)
         WD_coef(:,10) = (/ -1.79E-02,   0.137457, -7.45E-03,   0.232805 /)
         WD_coef(:,11) = (/ -5.07E-04,  -1.20E-03, -0.280138,   9.13E-02 /)
      CASE ( 13 )    ! 310-330 m
         WD_coef(:, 1) = (/  0.253568,  43.3822,   42.3741,   166.917    /)
         WD_coef(:, 2) = (/ -0.210713,  14.3161,   12.187,      9.66539  /)
         WD_coef(:, 3) = (/  0.176871,  -3.28688,  -2.78059,   -1.64384  /)
         WD_coef(:, 4) = (/  0.30952,    2.34743,  -5.8261,    -3.72051  /)
         WD_coef(:, 5) = (/ -0.211586,  -1.38792,  -0.891686,   3.26282  /)
         WD_coef(:, 6) = (/  0.114874,  -1.0177,   -2.95833,   -0.285227 /)
         WD_coef(:, 7) = (/ -0.168163,   1.33608,   5.32715,    0.270668 /)
         WD_coef(:, 8) = (/  0.106821,   0.746965, -1.28128,   -1.11127  /)
         WD_coef(:, 9) = (/ -2.17E-02,   0.198171,  0.911532,   2.31E-02 /)
         WD_coef(:,10) = (/ -5.64E-03,   0.278658,  0.250055,  -9.16E-02 /)
         WD_coef(:,11) = (/  7.21E-03,   2.24E-02,  6.76E-02,  -0.1011   /)
      CASE ( 14 )    ! 330-350 m
         WD_coef(:, 1) = (/  1.4365,   104.113,    86.7884,   138.082    /)
         WD_coef(:, 2) = (/  1.01951,  -22.4231,    8.14651,   -3.0374   /)
         WD_coef(:, 3) = (/ -0.14238,    5.5217,   -8.37098,    1.9052   /)
         WD_coef(:, 4) = (/ -8.04E-02,   2.56411,   8.01756,    0.450076 /)
         WD_coef(:, 5) = (/  7.34E-03,  -3.31792,  -10.0037,    1.66433  /)
         WD_coef(:, 6) = (/ -3.82E-02,   3.00083,   6.14358,   -0.656165 /)
         WD_coef(:, 7) = (/  0.113861,  -4.41267,  -2.98194,   -1.24882  /)
         WD_coef(:, 8) = (/ -0.154066,   4.29174,   3.74587,    1.4816   /)
         WD_coef(:, 9) = (/  0.127996,  -2.88696,  -2.49795,   -1.24336  /)
         WD_coef(:,10) = (/ -6.71E-02,   1.70388,   0.935254,   0.748082 /)
         WD_coef(:,11) = (/  8.19E-03,  -4.50E-02, -0.263839,  -5.18E-02 /)
      CASE ( 15 )    ! 350-370 m
         WD_coef(:, 1) = (/ -0.675054, 121.016,     0.173435, 199.751    /)
         WD_coef(:, 2) = (/ -0.52795,   26.7663,   36.6465,     8.14164  /)
         WD_coef(:, 3) = (/  0.686068,  -2.58652,   1.37125,  -12.8021   /)
         WD_coef(:, 4) = (/ -0.115391,  -0.715049,  0.225913,   2.68255  /)
         WD_coef(:, 5) = (/  0.127924,   1.18619,  -3.81934,   -2.40047  /)
         WD_coef(:, 6) = (/ -0.201212,  -1.51136,   4.51548,    3.23679  /)
         WD_coef(:, 7) = (/  0.175571,  -0.664591, -5.74074,   -2.24143  /)
         WD_coef(:, 8) = (/ -0.107098,   0.889236,  3.25149,    1.18349  /)
         WD_coef(:, 9) = (/  3.15E-02,  -6.48E-02, -0.882842,  -0.404645 /)
         WD_coef(:,10) = (/ -9.69E-03,  -0.486174, -0.284323,   0.336898 /)
         WD_coef(:,11) = (/  1.04E-03,  -0.144399, -6.10E-02,   6.62E-02 /)
      CASE ( 16 )    ! 370-390 m
         WD_coef(:, 1) = (/  0.610558, -90.3161,  -86.1311,   221.346    /)
         WD_coef(:, 2) = (/ -0.878196,   0.234356, -1.96802,   30.3835   /)
         WD_coef(:, 3) = (/  0.536954,   2.31986,   0.611791, -11.624    /)
         WD_coef(:, 4) = (/ -0.203843,  -2.10521,  -1.77538,    5.20693  /)
         WD_coef(:, 5) = (/ -6.04E-02,  -1.53784,   0.391834,   1.09004  /)
         WD_coef(:, 6) = (/ -3.32E-02,   1.08307,   0.756223,   0.579045 /)
         WD_coef(:, 7) = (/  2.20E-03,   1.00851,   0.872176,  -1.24302  /)
         WD_coef(:, 8) = (/ -4.70E-02,   0.313443, -5.20E-02,   1.24129  /)
         WD_coef(:, 9) = (/  0.105906,   2.60251,  -0.805126,  -2.35033  /)
         WD_coef(:,10) = (/ -3.95E-02,  -0.866726,  0.244709,   0.996069 /)
         WD_coef(:,11) = (/  5.34E-02,   0.423689, -0.910358,  -0.888237 /)
      CASE ( 17 )    ! 390-410 m
         WD_coef(:, 1) = (/ -0.256694, -53.0924,  -28.899,    212.286    /)
         WD_coef(:, 2) = (/  0.368178,   0.200188,-15.1321,     9.40209  /)
         WD_coef(:, 3) = (/ -0.102825,  -4.83546,   9.24228,   -0.64019  /)
         WD_coef(:, 4) = (/  0.191961,   2.99238,  -4.8869,    -2.80575  /)
         WD_coef(:, 5) = (/ -9.33E-02,   0.237869,  3.72573,   -8.03E-02 /)
         WD_coef(:, 6) = (/  1.70E-02,   2.22246,  -0.874,      0.324301 /)
         WD_coef(:, 7) = (/ -4.39E-02,  -1.22545,   1.03253,   -7.41E-02 /)
         WD_coef(:, 8) = (/  9.07E-03,  -0.438369, -1.85468,    0.746178 /)
         WD_coef(:, 9) = (/ -2.97E-02,  -0.626331,  1.32958,    0.161941 /)
         WD_coef(:,10) = (/ -4.73E-03,  -0.639604, -0.50062,    0.398523 /)
         WD_coef(:,11) = (/  7.78E-04,   0.203885,  0.111938,  -9.66E-02 /)
      CASE ( 18 )    ! 410-430 m
         WD_coef(:, 1) = (/ -1.05454,   19.3432,   14.3866,   209.914    /)
         WD_coef(:, 2) = (/ -5.37E-02,  -6.69143,  -5.48868,   13.8188   /)
         WD_coef(:, 3) = (/  0.130461,   1.84379,  10.2975,    -6.85151  /)
         WD_coef(:, 4) = (/  0.120135,   3.25255,  -4.64527,   -0.957415 /)
         WD_coef(:, 5) = (/ -0.157071,  -1.87681,   4.37492,    1.52585  /)
         WD_coef(:, 6) = (/  0.220174,   1.14707,  -5.27774,   -2.10403  /)
         WD_coef(:, 7) = (/ -0.185849,  -8.73E-02,  4.5702,     1.45097  /)
         WD_coef(:, 8) = (/  5.77E-02,  -0.265271, -2.17262,    1.19E-02 /)
         WD_coef(:, 9) = (/ -3.19E-02,   0.159054,  1.11463,    9.91E-02 /)
         WD_coef(:,10) = (/ -9.31E-03,  -0.514427, -0.486658,   0.472324 /)
         WD_coef(:,11) = (/  5.84E-03,  -6.98E-02, -6.53E-02,  -7.68E-02 /)
      CASE ( 19 )    ! 430-450 m
         WD_coef(:, 1) = (/  0.624689,  63.9533, -115.139,    203.718    /)
         WD_coef(:, 2) = (/ -0.249911,   8.56489,  12.0426,    11.2274   /)
         WD_coef(:, 3) = (/  0.208499,  -2.38494,   8.76157,   -7.17681  /)
         WD_coef(:, 4) = (/ -0.205812,   3.60713,   5.60652,    2.51439  /)
         WD_coef(:, 5) = (/  0.320606,  -7.16713, -10.6408,    -3.32927  /)
         WD_coef(:, 6) = (/ -0.178674,   5.15743,   3.70481,    2.92097  /)
         WD_coef(:, 7) = (/  0.101549,  -5.22916,  -1.89887,   -1.64557  /)
         WD_coef(:, 8) = (/ -9.30E-02,   2.8729,    1.14221,    1.4604   /)
         WD_coef(:, 9) = (/  1.45E-02,  -1.29998,  -0.491218,  -6.91E-02 /)
         WD_coef(:,10) = (/ -6.95E-04,   0.830442,  1.25591,   -0.451134 /)
         WD_coef(:,11) = (/ -6.90E-04,   1.30E-02, -0.16423,    7.65E-02 /)
      CASE ( 20 )    ! 450-470 m
         WD_coef(:, 1) = (/  4.30205,   83.823,   -77.8869,   120.115    /)
         WD_coef(:, 2) = (/  0.11147,   -2.13123, -13.0305,    11.4506   /)
         WD_coef(:, 3) = (/  5.36E-02,  -9.82942,   3.21203,   -2.14437  /)
         WD_coef(:, 4) = (/  3.12E-02,  -0.694,    -2.56494,    0.846492 /)
         WD_coef(:, 5) = (/ -3.97E-02,   0.628515,  0.898384,  -0.403596 /)
         WD_coef(:, 6) = (/  0.187725,  -1.32489,  -3.10108,   -1.64756  /)
         WD_coef(:, 7) = (/ -8.75E-02,  -0.750003,  1.2358,     0.95118  /)
         WD_coef(:, 8) = (/  4.29E-02,   0.206995, -0.591777,  -0.495133 /)
         WD_coef(:, 9) = (/ -3.25E-02,   0.187007,  0.351131,   0.374602 /)
         WD_coef(:,10) = (/ -1.79E-02,  -0.651232, -0.437205,   0.653204 /)
         WD_coef(:,11) = (/  5.74E-03,   0.210108, -0.185616,  -8.91E-02 /)
      CASE ( 21 )    ! 470-490 m
         WD_coef(:, 1) = (/  0.685959,  76.5757,  -26.8137,   187.31     /)
         WD_coef(:, 2) = (/ -0.229648,   3.36903, -12.3466,    19.5787   /)
         WD_coef(:, 3) = (/  5.56E-02,  -6.33886,   2.64958,   -2.35925  /)
         WD_coef(:, 4) = (/ -3.42E-02,  -1.78314,   1.51304,    0.43034  /)
         WD_coef(:, 5) = (/  5.81E-02,   4.2818,   -1.08668,   -2.13185  /)
         WD_coef(:, 6) = (/ -1.94E-02,  -2.76039,  -0.573698,   1.97694  /)
         WD_coef(:, 7) = (/  1.26E-02,   0.932315,  0.974862,  -1.5273   /)
         WD_coef(:, 8) = (/  1.04E-02,  -0.143063, -0.728002,   0.464589 /)
         WD_coef(:, 9) = (/  1.21E-03,   0.262702, -0.133363,  -0.236706 /)
         WD_coef(:,10) = (/ -2.29E-04,  -0.162697, -0.138587,   0.17236  /)
         WD_coef(:,11) = (/  6.61E-03,  -5.47E-02, -0.104054,  -9.64E-02 /)
      CASE DEFAULT
         CALL ProgAbort ('Error getting WD coefficients' )
   END SELECT

   RETURN
END SUBROUTINE get_coefs
!=======================================================================
FUNCTION getUStarProfile(p, WS, Ht, UStarOffset, UstarSlope)

   IMPLICIT                              NONE

   TYPE(TurbSim_ParameterType), INTENT(IN)     ::  p                 !< parameters 
   REAL(ReKi),                  INTENT(IN)     :: Ht(:)                       ! Height at which ustar is defined
   REAL(ReKi),                  INTENT(IN)     :: WS(:)                       ! Wind speed(s) at heights, Ht
   REAL(ReKi),                  INTENT(IN)     :: UStarOffset                 ! A scaling/offset value used with the Ustar_profile to ensure that the mean hub u'w' and ustar inputs agree with the profile values
   REAL(ReKi),                  INTENT(IN)     :: UstarSlope                  ! A scaling/slope value used with the Ustar_profile to ensure that the mean hub u'w' and ustar inputs agree with the profile values

   REAL(ReKi)                                  :: tmpZ                        ! a temporary value
   REAL(ReKi)                                  :: getUStarProfile(SIZE(Ht))       ! the array of ustar values
                                               
   INTEGER(IntKi)                              :: IZ
   INTEGER(IntKi)                              :: Zindx
   INTEGER(IntKi)                              :: Zindx_mn (1)
   INTEGER(IntKi)                              :: Zindx_mx (1)
                                               
   LOGICAL                                     :: mask(SIZE(Ht))

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

            getUStarProfile(  IZ) = ( 0.045355367 +  4.47275E-8*tmpZ**3)                                                      &
                                  + ( 0.511491978 -  0.09691157*LOG(tmpZ) - 199.226951/tmpZ**2           ) * WS(Zindx)        &
                                  + (-0.00396447  - 55.7818832/tmpZ**2                                   ) * p%met%RICH_NO    &
                                  + (-5.35764429  +  0.102002162*tmpZ/LOG(tmpZ) + 25.30585136/SQRT(tmpZ) ) * p%met%UstarDiab
         ENDDO

      ELSE ! All are above the max height so we'll use the old relationship at all heights
         getUStarProfile(:) = 0.17454 + 0.72045*p%met%UstarDiab**1.36242
      ENDIF

   ELSE ! All are below the min height so we'll use the diabatic Ustar value
      getUStarProfile(:) = p%met%UstarDiab
   ENDIF

   getUStarProfile = UstarSlope * getUStarProfile(:) + UstarOffset  ! These terms are used to make the ustar profile match the rotor-disk averaged value and input hub u'w'

END FUNCTION
!=======================================================================
FUNCTION getZLProfile(WS, Ht, RichNo, ZL, L, ZLOffset, WindProfileType)

   IMPLICIT                              NONE

   
   REAL(ReKi),   INTENT(IN)           :: Ht(:)                       ! Height at which local z/L is defined
   REAL(ReKi),   INTENT(IN)           :: WS(:)                       ! Wind speed(s) at heights, Ht
   REAL(ReKi),   INTENT(IN)           :: RichNo                      ! Richardson Number
   REAL(ReKi),   INTENT(IN)           :: ZL                          ! z/L, an alternate measure of stability (M-O) for RichNo
   REAL(ReKi),   INTENT(IN)           :: L                           ! L, M-O length
   REAL(ReKi),   INTENT(IN)           :: ZLOffset                    ! Offset to align profile with rotor-disk averaged z/L 

   CHARACTER(*), INTENT(IN)           :: WindProfileType
   
   REAL(ReKi)                         :: tmpZ                        ! a temporary value
   REAL(ReKi)                         :: getZLProfile(SIZE(Ht))          ! the array of z/L values

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
               IF ( RichNo >= 0 ) THEN
                  getZLProfile( IZ) =                     - 0.352464*RichNo + 0.005272*WS(Zindx) + 0.465838
               ELSE
                  getZLProfile( IZ) =  0.004034*Ht(Zindx) + 0.809494*RichNo - 0.008298*WS(Zindx) - 0.386632
               ENDIF !RichNo
            ELSE
               IF ( RichNo >= 0 ) THEN
                  getZLProfile( IZ) =  0.003068*Ht(Zindx) + 1.140264*RichNo + 0.036726*WS(Zindx) - 0.407269
               ELSE
                  getZLProfile( IZ) =  0.003010*Ht(Zindx) + 0.942617*RichNo                      - 0.221886
               ENDIF 
            ENDIF
            getZLProfile( IZ) = MIN( getZLProfile( IZ), 1.0_ReKi )
            getZLProfile( IZ) = getZLProfile(IZ) * tmpZ

         ENDDO

      ELSE ! All are above the max height so instead of extrapolating, we'll use ZL at all heights
         getZLProfile(:) = ZL
      ENDIF

   ELSE ! All are below the min height so we'll keep L constant (as is the case in the surface layer)
      getZLProfile(:) = Ht(:) / L
   ENDIF

   getZLProfile = getZLProfile(:) + ZLOffset  ! This offset term is used to make the zl profile match the rotor-disk averaged value


END FUNCTION getZLProfile
!=======================================================================
END MODULE TS_Profiles
