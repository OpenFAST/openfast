!**********************************************************************************************************************************
! The Morison and Morison_Types modules make up a template for creating user-defined calculations in the FAST Modularization 
! Framework. Morison_Types will be auto-generated based on a description of the variables for the module.
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012-2015  National Renewable Energy Laboratory
!
!    This file is part of Morison.
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
!    
!**********************************************************************************************************************************
MODULE Morison
   USE Waves
   USE Morison_Types  
   USE Morison_Output
  ! USE HydroDyn_Output_Types
   USE NWTC_Library

   
   IMPLICIT NONE
   
   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: Morison_ProgDesc = ProgDesc( 'Morison', '', '' )

   
      ! ..... Public Subroutines ...................................................................................................
   PUBLIC:: Morison_ProcessMorisonGeometry
   
   PUBLIC :: Morison_Init                           ! Initialization routine
   PUBLIC :: Morison_End                            ! Ending routine (includes clean up)
   
   PUBLIC :: Morison_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating 
                                                    !   continuous states, and updating discrete states
   PUBLIC :: Morison_CalcOutput                     ! Routine for computing outputs
   
   PUBLIC :: Morison_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: Morison_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: Morison_UpdateDiscState                ! Tight coupling routine for updating discrete states
      
   
   
CONTAINS

SUBROUTINE Morison_DirCosMtrx( pos0, pos1, DirCos )

! Compute the direction cosine matrix given two points along the axis of a cylinder

   REAL(ReKi), INTENT( IN    )  ::   pos0(3), pos1(3)
   Real(ReKi), INTENT(   OUT )  ::   DirCos(3,3)
   Real(DbKi)                   ::   xz, xyz
   Real(DbKi)                   ::   x0, y0, z0
   Real(DbKi)                   ::   x1, y1, z1
!   Real(DbKi)                   ::   temp

   x0 = pos0(1)
   y0 = pos0(2)
   z0 = pos0(3)
   x1 = pos1(1)
   y1 = pos1(2)
   z1 = pos1(3)
   
      ! Need to verify that z0 <= z1, but this was already handled in the element construction process!!! GJH 9/24/13 
   !IF ( z0 > z1 ) THEN
   !   temp = x0
   !   x0   = x1
   !   x1   = temp
   !   temp = y0
   !   y0   = y1
   !   y1   = temp
   !   temp = z0
   !   z0   = z1
   !   z1   = temp
   !END IF
   
   xz  = sqrt((x0-x1)*(x0-x1)+(z0-z1)*(z0-z1))
   xyz = sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1)+(z0-z1)*(z0-z1))
   
   IF ( xz==0 ) THEN
      
      IF (y1<y0) THEN
         
         DirCos = transpose(reshape((/ 1, 0, 0, 0, 0, -1, 0, 1, 0 /), shape(DirCos)))
          
      ELSE
         
         DirCos = transpose(reshape((/ 1, 0, 0, 0, 0, 1, 0, -1, 0 /), shape(DirCos)))
         
      END IF
      
   ELSE
      
      DirCos(1, 1) = -(z0-z1)/xz
      DirCos(1, 2) = -(x0-x1)*(y0-y1)/(xz*xyz)
      DirCos(1, 3) = (x1-x0)/xyz
      
      DirCos(2, 1) = 0.0
      DirCos(2, 2) = xz/xyz
      DirCos(2, 3) = (y1-y0)/xyz
      
      DirCos(3, 1) = -(x1-x0)/xz
      DirCos(3, 2) = -(y0-y1)*(z0-z1)/(xz*xyz)
      DirCos(3, 3) = (z1-z0)/xyz
      
      ! DEBUG:  TODO : Remove
      !PRINT*, sqrt(DirCos(1,1)*DirCos(1,1) + DirCos(1,2)*DirCos(1,2)+DirCos(1,3)*DirCos(1,3))
      !PRINT*, sqrt(DirCos(2,1)*DirCos(2,1) + DirCos(2,2)*DirCos(2,2)+DirCos(2,3)*DirCos(2,3))
      !PRINT*, sqrt(DirCos(3,1)*DirCos(3,1) + DirCos(3,2)*DirCos(3,2)+DirCos(3,3)*DirCos(3,3))
   END IF    
   
END SUBROUTINE Morison_DirCosMtrx

!====================================================================================================
SUBROUTINE GetDistance ( a, b, l )
!    This private subroutine computes the distance between points a and b.
!---------------------------------------------------------------------------------------------------- 

   REAL(ReKi), INTENT ( IN    )  :: a(3)     ! the position of point a
   REAL(ReKi), INTENT ( IN    )  :: b(3)     ! the position of point b
   REAL(ReKi), INTENT (   OUT )  :: l        ! the distance between point a and b
   
   l = sqrt( ( a(1) - b(1) ) * ( a(1) - b(1) ) + ( a(2) - b(2) ) * ( a(2) - b(2) ) + ( a(3) - b(3) ) * ( a(3) - b(3) ) )
   
END SUBROUTINE GetDistance

!====================================================================================================
SUBROUTINE ElementCentroid ( Rs, Re, p1, h, DCM, centroid )
!    This private subroutine computes the centroid of a tapered right cylinder element.
!---------------------------------------------------------------------------------------------------- 

   REAL(ReKi), INTENT ( IN    )  :: Rs          ! starting radius
   REAL(ReKi), INTENT ( IN    )  :: Re          ! ending radius
   REAL(ReKi), INTENT ( IN    )  :: p1(3)       ! starting point of the element in global coordinates
   REAL(ReKi), INTENT ( IN    )  :: h           ! height of the element
   REAL(ReKi), INTENT ( IN    )  :: DCM(3,3)    ! direction cosine matrix to transform local element coordinates to global coordinates
   REAL(ReKi), INTENT (   OUT )  :: centroid(3) ! centroid of the element in local coordinates
   
   centroid(1) = 0.0
   centroid(2) = 0.0
   centroid(3) = h * (Rs*Rs + 2.0*Rs*Re +  3.0*Re*Re) / (4.0*( Rs*Rs + Rs*Re +  Re*Re  ) )                    !( 2.0*Re + Rs ) / ( 3.0 * ( Rs + Re ) )
   centroid    = matmul( DCM, centroid ) + p1
   
END SUBROUTINE ElementCentroid

!====================================================================================================
REAL(ReKi) FUNCTION ElementVolume ( Rs, Re, h )
!    This private function computes the volume of a tapered right cylinder element.
!---------------------------------------------------------------------------------------------------- 

   REAL(ReKi), INTENT ( IN    )  :: Rs          ! starting radius
   REAL(ReKi), INTENT ( IN    )  :: Re          ! ending radius
   REAL(ReKi), INTENT ( IN    )  :: h           ! height of the element
   
   ElementVolume = Pi*h*( Rs*Rs + Re*Re + Rs*Re  ) / 3.0
   
END FUNCTION ElementVolume

!====================================================================================================
SUBROUTINE    FindInterpFactor( p, p1, p2, s )

   REAL(ReKi),  INTENT ( IN    )  :: p, p1, p2
   REAL(ReKi),  INTENT (   OUT )  :: s
   
   REAL(ReKi)                     :: dp
! find normalized interpolation factor, s, such:
! p = p1*(1-s) + p2*s
!  *--------------*--------------------------------*
!  p1             p                                p2
!
!  0-----------------------------------------------1
!  <------- s ---->
   
   dp = p2 - p1
   IF ( EqualRealNos(dp, 0.0_ReKi) ) THEN
      s = 0
   ELSE
      s = ( p - p1  ) / dp 
   END IF
         
END SUBROUTINE FindInterpFactor
!=======================================================================
FUNCTION InterpWrappedStpInt( XValIn, XAry, YAry, Ind, AryLen )


      ! This funtion returns a y-value that corresponds to an input x-value which is wrapped back
      ! into the range [1-XAry(AryLen).  It finds a x-value which corresponds to a value in the XAry where XAry(Ind-1) < MOD(XValIn, XAry(AryLen)) <= XAry(Ind)
      ! It is assumed that XAry is sorted in ascending order.
      ! It uses the passed index as the starting point and does a stepwise interpolation from there.  This is
      ! especially useful when the calling routines save the value from the last time this routine was called
      ! for a given case where XVal does not change much from call to call.  .
      ! 
      ! This routine assumes YAry is INTEGER.


      ! Function declaration.

   INTEGER                  :: InterpWrappedStpInt                                  ! This function.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the arrays.
   INTEGER, INTENT(INOUT)       :: Ind                                             ! Initial and final index into the arrays.

   REAL(SiKi), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(SiKi), INTENT(IN)       :: XValIn                                          ! X value to be interpolated.
   INTEGER, INTENT(IN)          :: YAry    (AryLen)                                ! Array of Y values to be interpolated.

   REAL(SiKi)                   :: XVal                                            ! X value to be interpolated.
   
   
   
      ! Wrap XValIn into the range XAry(1) to XAry(AryLen)
   XVal = MOD(XValIn, XAry(AryLen))

      ! Set the Ind to the first index if we are at the beginning of XAry
   IF ( XVal <= XAry(2) )  THEN  
      Ind           = 1
   END IF
   
   
        ! Let's check the limits first.

   IF ( XVal <= XAry(1) )  THEN
      InterpWrappedStpInt = YAry(1)
      Ind           = 1
      RETURN
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      InterpWrappedStpInt = YAry(AryLen)
      Ind           = MAX(AryLen - 1, 1)
      RETURN
   END IF


     ! Let's interpolate!

   Ind = MAX( MIN( Ind, AryLen-1 ), 1 )

   DO

      IF ( XVal < XAry(Ind) )  THEN

         Ind = Ind - 1

      ELSE IF ( XVal >= XAry(Ind+1) )  THEN

         Ind = Ind + 1

      ELSE

         InterpWrappedStpInt = YAry(Ind) 
         RETURN

      END IF

   END DO

   RETURN
END FUNCTION InterpWrappedStpInt ! ( XVal, XAry, YAry, Ind, AryLen )
   
   
!=======================================================================
FUNCTION InterpWrappedStpLogical( XValIn, XAry, YAry, Ind, AryLen )


      ! This funtion returns a y-value that corresponds to an input x-value which is wrapped back
      ! into the range [0-XAry(AryLen) by interpolating into the arrays.  
      ! It is assumed that XAry is sorted in ascending order.
      ! It uses the passed index as the starting point and does a stepwise interpolation from there.  This is
      ! especially useful when the calling routines save the value from the last time this routine was called
      ! for a given case where XVal does not change much from call to call.  When there is no correlation
      ! from one interpolation to another, InterpBin() may be a better choice.
      ! It returns the first or last YAry() value if XVal is outside the limits of XAry().
      ! This routine assumes YAry is REAL.


      ! Function declaration.

   LOGICAL                  :: InterpWrappedStpLogical                                  ! This function.


      ! Argument declarations.

   INTEGER, INTENT(IN)          :: AryLen                                          ! Length of the arrays.
   INTEGER, INTENT(INOUT)       :: Ind                                             ! Initial and final index into the arrays.

   REAL(SiKi), INTENT(IN)       :: XAry    (AryLen)                                ! Array of X values to be interpolated.
   REAL(SiKi), INTENT(IN)       :: XValIn                                           ! X value to be interpolated.
   LOGICAL, INTENT(IN)          :: YAry    (AryLen)                                ! Array of Y values to be interpolated.

   REAL(SiKi)                   :: XVal                                           ! X value to be interpolated.
   
   
   
      ! Wrap XValIn into the range XAry(1) to XAry(AryLen)
   XVal = MOD(XValIn, XAry(AryLen))

      ! Set the Ind to the first index if we are at the beginning of XAry
   IF ( XVal <= XAry(2) )  THEN  
      Ind           = 1
   END IF
   
   
        ! Let's check the limits first.

   IF ( XVal <= XAry(1) )  THEN
      InterpWrappedStpLogical = YAry(1)
      Ind           = 1
      RETURN
   ELSE IF ( XVal >= XAry(AryLen) )  THEN
      InterpWrappedStpLogical = YAry(AryLen)
      Ind           = MAX(AryLen - 1, 1)
      RETURN
   END IF


     ! Let's interpolate!

   Ind = MAX( MIN( Ind, AryLen-1 ), 1 )

   DO

      IF ( XVal < XAry(Ind) )  THEN

         Ind = Ind - 1

      ELSE IF ( XVal >= XAry(Ind+1) )  THEN

         Ind = Ind + 1

      ELSE

         InterpWrappedStpLogical = YAry(Ind) 
         RETURN

      END IF

   END DO

   RETURN
END FUNCTION InterpWrappedStpLogical ! ( XVal, XAry, YAry, Ind, AryLen )

SUBROUTINE DistrBuoyancy( densWater, R, tMG, dRdz, Z, C, g, F_B  ) 

   REAL(ReKi),         INTENT ( IN    )  :: densWater
   REAL(ReKi),         INTENT ( IN    )  :: R
   REAL(ReKi),         INTENT ( IN    )  :: tMG
   REAL(ReKi),         INTENT ( IN    )  :: dRdz
   REAL(ReKi),         INTENT ( IN    )  :: Z
   REAL(ReKi),         INTENT ( IN    )  :: C(3,3)
   REAL(ReKi),         INTENT ( IN    )  :: g
   REAL(ReKi),         INTENT (   OUT )  :: F_B(6)
   
   REAL(DbKi)                           :: Reff,ReffSq,ReffCub,f1,f2,f3
   
   REAL(DbKi) :: CC(3,3)
   
   CC    = REAL(C,DbKi)
   Reff  = REAL(R + tMG,DbKi)
  
   
   
   ReffSq  = Reff*Reff 
   ReffCub = ReffSq*Reff
   f1      = REAL(denswater,DbKi)*REAL(g,DbKi)*Pi_D
   f2      = f1*ReffCub*REAL(dRdz,DbKi)
   f3      = Reff*REAL(dRdz,DbKi)*REAL(Z,DbKi)
  
   
   F_B(1) = f1*( (CC(1,1)*CC(3,1) + CC(1,2)*CC(3,2))*ReffSq - 2.0*CC(1,3)*f3 )
   F_B(2) = f1*( (CC(2,1)*CC(3,1) + CC(2,2)*CC(3,2))*ReffSq - 2.0*CC(2,3)*f3 )
   F_B(3) = f1*( (CC(3,1)*CC(3,1) + CC(3,2)*CC(3,2))*ReffSq - 2.0*CC(3,3)*f3 )
   F_B(4) = -f2*( CC(1,1)*CC(3,2) - CC(1,2)*CC(3,1) )
   F_B(5) = -f2*( CC(2,1)*CC(3,2) - CC(2,2)*CC(3,1) )
   F_B(6) = -f2*( CC(3,1)*CC(3,2) - CC(3,2)*CC(3,1) )
   
  
   
   
END SUBROUTINE DistrBuoyancy


SUBROUTINE DistrInertialLoads( nodeIndx, densWater, Ca, Cp, AxCa, AxCp, R, tMG, dRdZ, k, NStepWave, WaveAcc0, WaveDynP0, F_I, ErrStat, ErrMsg  )

   INTEGER,            INTENT ( IN    )  :: nodeIndx
   REAL(ReKi),         INTENT ( IN    )  :: densWater
   REAL(ReKi),         INTENT ( IN    )  :: Ca
   REAL(ReKi),         INTENT ( IN    )  :: Cp
   REAL(ReKi),         INTENT ( IN    )  :: AxCa
   REAL(ReKi),         INTENT ( IN    )  :: AxCp
   REAL(ReKi),         INTENT ( IN    )  :: R
   REAL(ReKi),         INTENT ( IN    )  :: tMG
   REAL(ReKi),         INTENT ( IN    )  :: dRdZ
   REAL(ReKi),         INTENT ( IN    )  :: k(3)
   INTEGER,            INTENT ( IN    )  :: NStepWave
   REAL(SiKi),         INTENT ( IN    )  :: WaveAcc0(0:,:,:)
   REAL(SiKi),         INTENT ( IN    )  :: WaveDynP0(0:,:)
   REAL(ReKi),ALLOCATABLE,  INTENT (   OUT )  :: F_I(:,:)
   INTEGER,            INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),       INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None

   INTEGER                               :: I
   REAL(ReKi)                            :: f, f1, f2, f3, adotk !, v_len
   REAL(ReKi)                            :: v(3), af(3) !p0(3), m(3), 
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = "" 
      
      ! Allocate F_I
   ALLOCATE ( F_I(0:NStepWave, 6), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating distributed inertial loads array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF  
   
   f  = (Ca + Cp)*densWater*Pi*(R+tMG)*(R+tMG) 
   f2 = AxCa*densWater*2.0*Pi*(R+tMG)*(R+tMG)*abs(dRdZ)       
   f1 = AxCp*2.0*Pi*(R+tMG)*dRdz
   
   DO I=0,NStepWave
      
      af    =  WaveAcc0(I,nodeIndx,:)       
      adotk = af(1)*k(1) + af(2)*k(2) + af(3)*k(3)   
      v     =  af - adotk*k
    
      ! NOTE: (k cross l) x k = l - (l dot k)k
      
      f3 = f1*WaveDynP0(I,nodeIndx)
      
      !CALL GetDistance( p0, v, v_len )  
      !TODO What about multiplying by the magnitude?
      
      
      
      F_I(I,1) = f*v(1) + (f3 + f2*adotk)*k(1)
      F_I(I,2) = f*v(2) + (f3 + f2*adotk)*k(2) 
      F_I(I,3) = f*v(3) + (f3 + f2*adotk)*k(3)
      F_I(I,4) = 0.0
      F_I(I,5) = 0.0
      F_I(I,6) = 0.0
      
   END DO
   
END SUBROUTINE DistrInertialLoads



SUBROUTINE DistrInertialLoads2( densWater, Ca, Cp, AxCa, AxCp, R, tMG, dRdZ, k, WaveAcc, WaveDynP, F_I, ErrStat, ErrMsg  )
                  
   REAL(ReKi),         INTENT ( IN    )  :: densWater
   REAL(ReKi),         INTENT ( IN    )  :: Ca
   REAL(ReKi),         INTENT ( IN    )  :: Cp
   REAL(ReKi),         INTENT ( IN    )  :: AxCa
   REAL(ReKi),         INTENT ( IN    )  :: AxCp
   REAL(ReKi),         INTENT ( IN    )  :: R
   REAL(ReKi),         INTENT ( IN    )  :: tMG
   REAL(ReKi),         INTENT ( IN    )  :: dRdZ
   REAL(ReKi),         INTENT ( IN    )  :: k(3)
   REAL(ReKi),         INTENT ( IN    )  :: WaveAcc(3)
   REAL(ReKi),         INTENT ( IN    )  :: WaveDynP
   REAL(ReKi),         INTENT (   OUT )  :: F_I(3)
   !REAL(ReKi),         INTENT (   OUT )  :: F_I(3)
   INTEGER,            INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),       INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None

   INTEGER                               :: I
   REAL(ReKi)                            :: f, f1, f2, f3, v_len, adotk
   REAL(ReKi)                            :: p0(3), m(3), v(3), af(3)
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = "" 
      

   
   f  = (Ca + Cp)*densWater*Pi*(R+tMG)*(R+tMG) 
   f2 = AxCa*densWater*2.0*Pi*(R+tMG)*(R+tMG)*abs(dRdZ)       
   f1 = AxCp*2.0*Pi*(R+tMG)*dRdz
   
   adotk = WaveAcc(1)*k(1) + WaveAcc(2)*k(2) + WaveAcc(3)*k(3)   
   v     =  WaveAcc - adotk*k
    
   ! NOTE: (k cross l) x k = l - (l dot k)k
      
   f3 = f1*WaveDynP
      
   !CALL GetDistance( p0, v, v_len )  
   !TODO What about multiplying by the magnitude?
      
      
      
   F_I(1) = f*v(1) + (f3 + f2*adotk)*k(1)
   F_I(2) = f*v(2) + (f3 + f2*adotk)*k(2) 
   F_I(3) = f*v(3) + (f3 + f2*adotk)*k(3)
   !F_I(4) = 0.0_ReKi
   !F_I(5) = 0.0_ReKi
   !F_I(6) = 0.0_ReKi
      
  
   
END SUBROUTINE DistrInertialLoads2


SUBROUTINE DistrMGLoads(MGdens, g, R, tMG, F_MG )  
   REAL(ReKi),         INTENT ( IN    )  :: MGdens
   REAL(ReKi),         INTENT ( IN    )  :: g
   REAL(ReKi),         INTENT ( IN    )  :: R
   REAL(ReKi),         INTENT ( IN    )  :: tMG
   REAL(ReKi),         INTENT (   OUT )  :: F_MG(6)
   
   F_MG(:) = 0.0
   F_MG(3) = -MGdens*g*Pi* ( (R + tMG ) * ( R + tMG ) - R*R )
   
END SUBROUTINE DistrMGLoads
          
SUBROUTINE DistrDragConst( densWater, Cd, R, tMG, DragConst  ) 

   ! This is used to minimize the computations which occur at each timestep
   
   REAL(ReKi),         INTENT ( IN    )  :: Cd
   REAL(ReKi),         INTENT ( IN    )  :: densWater
   REAL(ReKi),         INTENT ( IN    )  :: R
   REAL(ReKi),         INTENT ( IN    )  :: tMG
   REAL(ReKi),         INTENT (   OUT )  :: DragConst
   
   DragConst = Cd*densWater*(R+tMG)
   
END SUBROUTINE DistrDragConst


SUBROUTINE DistrFloodedBuoyancy( densFluid, Z_f, R, t, dRdz, Z, C, g, F_B  ) 

   REAL(ReKi),         INTENT ( IN    )  :: densFluid
   REAL(ReKi),         INTENT ( IN    )  :: Z_f
   REAL(ReKi),         INTENT ( IN    )  :: R
   REAL(ReKi),         INTENT ( IN    )  :: t
   REAL(ReKi),         INTENT ( IN    )  :: dRdz
   REAL(ReKi),         INTENT ( IN    )  :: Z
   REAL(ReKi),         INTENT ( IN    )  :: C(3,3)
   REAL(ReKi),         INTENT ( IN    )  :: g
   REAL(ReKi),         INTENT (   OUT )  :: F_B(6)
   
   REAL(DbKi)                           :: Zeff,Reff,ReffSq,ReffCub,f1,f2,f3
    
  
   REAL(DbKi) :: CC(3,3)
   CC = REAL(C,DbKi)
   
  
   Reff  =  REAL(R - t,DbKi)
   Zeff  = REAL(Z - Z_f,DbKi)
   
   ReffSq  = Reff*Reff 
   ReffCub = ReffSq*Reff
   f1      = -REAL(densFluid,DbKi)*REAL(g,DbKi)*Pi_D
   f2      = f1*ReffCub*REAL(dRdz,DbKi)
   f3      = Reff*REAL(dRdz,DbKi)*Zeff
   
  
   
   F_B(1) = f1*( (CC(1,1)*CC(3,1) + CC(1,2)*CC(3,2))*ReffSq - 2.0*CC(1,3)*f3 )
   F_B(2) = f1*( (CC(2,1)*CC(3,1) + CC(2,2)*CC(3,2))*ReffSq - 2.0*CC(2,3)*f3 )
   F_B(3) = f1*( (CC(3,1)*CC(3,1) + CC(3,2)*CC(3,2))*ReffSq - 2.0*CC(3,3)*f3 )
   F_B(4) = -f2*( CC(1,1)*CC(3,2) - CC(1,2)*CC(3,1) )
   F_B(5) = -f2*( CC(2,1)*CC(3,2) - CC(2,2)*CC(3,1) )
   F_B(6) = -f2*( CC(3,1)*CC(3,2) - CC(3,2)*CC(3,1) )
   
END SUBROUTINE DistrFloodedBuoyancy

SUBROUTINE DistrAddedMass( densWater, Ca, AxCa, C, R, tMG, dRdZ, AM_M)

   REAL(ReKi),         INTENT ( IN    )  :: densWater
   REAL(ReKi),         INTENT ( IN    )  :: Ca
   REAL(ReKi),         INTENT ( IN    )  :: AxCa
   REAL(ReKi),         INTENT ( IN    )  :: C(3,3)
   REAL(ReKi),         INTENT ( IN    )  :: R
   REAL(ReKi),         INTENT ( IN    )  :: tMG
   REAL(ReKi),         INTENT ( IN    )  :: dRdZ
   REAL(ReKi),         INTENT (   OUT )  :: AM_M(3,3)
   
   REAL(ReKi)                            :: f,f2
   
   f         = Ca*densWater*Pi*(R+tMG)*(R+tMG)
   f2        = AxCa*2.0*densWater*Pi*abs(dRdZ)*(R+tMG)*(R+tMG)
   !AM_M      = 0.0
   AM_M(1,1) = f*(  C(2,3)*C(2,3) + C(3,3)*C(3,3) ) -f2*C(1,3)*C(1,3)
   AM_M(1,2) = f*( -C(1,3)*C(2,3)                 ) -f2*C(1,3)*C(2,3)
   AM_M(1,3) = f*( -C(1,3)*C(3,3)                 ) -f2*C(1,3)*C(3,3)
   
   AM_M(2,1) = f*( -C(1,3)*C(2,3)                 ) -f2*C(2,3)*C(1,3)
   AM_M(2,2) = f*(  C(1,3)*C(1,3) + C(3,3)*C(3,3) ) -f2*C(2,3)*C(2,3)
   AM_M(2,3) = f*( -C(2,3)*C(3,3)                 ) -f2*C(2,3)*C(3,3)
   
   AM_M(3,1) = f*( -C(1,3)*C(3,3)                 ) -f2*C(3,3)*C(1,3)
   AM_M(3,2) = f*( -C(2,3)*C(3,3)                 ) -f2*C(3,3)*C(2,3)
   AM_M(3,3) = f*(  C(1,3)*C(1,3) + C(2,3)*C(2,3) ) -f2*C(3,3)*C(3,3)


END SUBROUTINE DistrAddedMass


SUBROUTINE DistrAddedMassMG( densMG, R, tMG, AM_MG)

   REAL(ReKi),         INTENT ( IN    )  :: densMG
   REAL(ReKi),         INTENT ( IN    )  :: R
   REAL(ReKi),         INTENT ( IN    )  :: tMG
   REAL(ReKi),         INTENT (   OUT )  :: AM_MG

   AM_MG = densMG*Pi*((R+tMG)*(R+tMG) - R*R)
  
   
END SUBROUTINE DistrAddedMassMG


SUBROUTINE DistrAddedMassFlood( densFluid, R, t, AM_F)

   REAL(ReKi),         INTENT ( IN    )  :: densFluid
   REAL(ReKi),         INTENT ( IN    )  :: R
   REAL(ReKi),         INTENT ( IN    )  :: t
   REAL(ReKi),         INTENT (   OUT )  :: AM_F
   
   
   AM_F = densFluid*Pi*(R-t)*(R-t)
  
END SUBROUTINE DistrAddedMassFlood
         
         

SUBROUTINE LumpDragConst( densWater, Cd, R, tMG, DragConst  ) 

   ! This is used to minimize the computations which occur at each timestep
   
   REAL(ReKi),         INTENT ( IN    )  :: Cd
   REAL(ReKi),         INTENT ( IN    )  :: densWater
   REAL(ReKi),         INTENT ( IN    )  :: R
   REAL(ReKi),         INTENT ( IN    )  :: tMG
   REAL(ReKi),         INTENT (   OUT )  :: DragConst
   
   DragConst = 0.5*Cd*densWater*(R+tMG)*(R+tMG)
   
END SUBROUTINE LumpDragConst

         
SUBROUTINE LumpDynPressure( nodeIndx, Cp, k, R, tMG, NStepWave, WaveDynP, F_DP, ErrStat, ErrMsg )


   INTEGER,            INTENT ( IN    )  :: nodeIndx
   REAL(ReKi),         INTENT ( IN    )  :: Cp
   REAL(ReKi),         INTENT ( IN    )  :: k(3)
   REAL(ReKi),         INTENT ( IN    )  :: R
   REAL(ReKi),         INTENT ( IN    )  :: tMG
   INTEGER,            INTENT ( IN    )  :: NStepWave
   REAL(SiKi),         INTENT ( IN    )  :: WaveDynP(0:,:)
   REAL(ReKi),ALLOCATABLE,         INTENT (   OUT )  :: F_DP(:,:)
   INTEGER,            INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),       INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None

   INTEGER                               :: I
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = "" 
   
   
      ! Allocate F_DP
      
   ALLOCATE ( F_DP(0:NStepWave,6), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating distributed dynamic pressure loads array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF  
   
   DO I=0,NStepWave
      F_DP(I,1) = Cp*k(1)*Pi*(R+tMG)*(R+tMG)*WaveDynP(I,nodeIndx) 
      F_DP(I,2) = Cp*k(2)*Pi*(R+tMG)*(R+tMG)*WaveDynP(I,nodeIndx) 
      F_DP(I,3) = Cp*k(3)*Pi*(R+tMG)*(R+tMG)*WaveDynP(I,nodeIndx) 
      F_DP(I,4) = 0.0
      F_DP(I,5) = 0.0
      F_DP(I,6) = 0.0
   END DO
   
   
END SUBROUTINE LumpDynPressure



SUBROUTINE LumpBuoyancy( sgn, densWater, R, tMG, Z, C, g, F_B  ) 

   REAL(ReKi),         INTENT ( IN    )  :: sgn
   REAL(ReKi),         INTENT ( IN    )  :: densWater
   REAL(ReKi),         INTENT ( IN    )  :: R
   REAL(ReKi),         INTENT ( IN    )  :: tMG
   REAL(ReKi),         INTENT ( IN    )  :: Z
   REAL(ReKi),         INTENT ( IN    )  :: C(3,3)
   REAL(ReKi),         INTENT ( IN    )  :: g
   REAL(ReKi),         INTENT (   OUT )  :: F_B(6)


   REAL(DbKi)                            :: f, f1, f2, f3, Reff, Rsq,R_4
   Reff = REAL(R+tMG,DbKi)
   Rsq  = Reff**2
   R_4  = Rsq**2
   f  = REAL(g,DbKi)*REAL(densWater,DbKi)*REAL(sgn,DbKi)
   f1 = -REAL(Z,DbKi)*Pi_D*Rsq
   f2 = f*Pi_D*R_4
   !f3 =  C(3,1)*R_4

   F_B(1) = C(1,3)*f1*f
   F_B(2) = C(2,3)*f1*f
   F_B(3) = C(3,3)*f1*f
   F_B(4) =  0.25*( -C(3,2)*C(1,1) + C(1,2)*C(3,1) )*f2   ! TODO: We flipped the signs of the moments because 1 member tapered integrated moments were not zero.  GJH 10/1/13  Jason is verifying.
   F_B(5) =  0.25*( -C(3,2)*C(2,1) + C(2,2)*C(3,1) )*f2
   F_B(6) =  0.25*( -C(3,2)*C(3,1) + C(3,2)*C(3,1) )*f2
   
   
END SUBROUTINE LumpBuoyancy



SUBROUTINE LumpFloodedBuoyancy( sgn, densFill, R, tM, FillFS, Z, C, g, F_BF  ) 

   REAL(ReKi),         INTENT ( IN    )  :: sgn
   REAL(ReKi),         INTENT ( IN    )  :: densFill
   REAL(ReKi),         INTENT ( IN    )  :: R
   REAL(ReKi),         INTENT ( IN    )  :: tM
   REAL(ReKi),         INTENT ( IN    )  :: FillFS
   REAL(ReKi),         INTENT ( IN    )  :: Z
   REAL(ReKi),         INTENT ( IN    )  :: C(3,3)
   REAL(ReKi),         INTENT ( IN    )  :: g
   REAL(ReKi),         INTENT (   OUT )  :: F_BF(6)

   
   REAL(ReKi)                            :: f, f1, f2 !, f3
   
   f  = -densFill*g*sgn
   
   f1 = -(Z - FillFS)*Pi*       (R-tM)*(R-tM)
   f2 = 0.25*Pi*(R-tM)*(R-tM)*(R-tM)*(R-tM)
   

   F_BF(1) = C(1,3)*f1*f
   F_BF(2) = C(2,3)*f1*f
   F_BF(3) = C(3,3)*f1*f
   F_BF(4) =  (-C(1,1)*C(3,2) + C(1,2)*C(3,1))*f*f2   ! TODO: We flipped the signs of the moments because 1 member tapered integrated moments were not zero.  GJH 10/1/13  Jason is verifying.
   F_BF(5) =  (-C(2,1)*C(3,2) + C(2,2)*C(3,1))*f*f2
   F_BF(6) =  (-C(3,1)*C(3,2) + C(3,2)*C(3,1))*f*f2
   
   
END SUBROUTINE LumpFloodedBuoyancy

LOGICAL FUNCTION IsThisSplitValueUsed(nSplits, splits, checkVal)

   INTEGER,                        INTENT    ( IN    )  :: nSplits
   REAL(ReKi),                     INTENT    ( IN    )  :: splits(:)
   REAL(ReKi),                     INTENT    ( IN    )  :: checkVal
   
   INTEGER             :: I
   
   DO I=1,nSplits
      IF ( EqualRealNos(splits(I), checkVal ) ) THEN
         IsThisSplitValueUsed = .TRUE.
         RETURN
      END IF
END DO

   IsThisSplitValueUsed = .FALSE.
   
END FUNCTION IsThisSplitValueUsed
!====================================================================================================
SUBROUTINE GetMaxSimQuantities( numMGDepths, MGTop, MGBottom, MSL2SWL, Zseabed, filledGroups, numJoints, joints, numMembers, members, maxNodes, maxElements, maxSuperMembers )
!     This private subroutine determines the maximum nodes, elements, and super members which may appear
!     in the final simulation mesh.  This is based on the following:
!     1) Member splitting at the marine growth boundaries ( new nodes and members )
!     2) Member splitting due to internal subdivision ( new nodes and members )
!     3) New nodes and super members if a joint marked with JointOvrlp = 1 (and additional conditions are satisfied)
!     
!---------------------------------------------------------------------------------------------------- 
   INTEGER,                        INTENT    ( IN    )  :: numMGDepths              ! number of MGDepths specified in the input table
   REAL(ReKi),                     INTENT    ( IN    )  :: MGTop                    ! Global Z-value of the upper marine growth boundary
   REAL(ReKi),                     INTENT    ( IN    )  :: MGBottom                 ! Global Z-value of the lower marine growth boundary
   REAL(ReKi),                     INTENT    ( IN    )  :: MSL2SWL                  ! Global Z-value of mean sea level
   REAL(ReKi),                     INTENT    ( IN    )  :: Zseabed                  ! Global Z-value of the top of the seabed
   TYPE(Morison_FilledGroupType),  INTENT    ( IN    )  :: filledGroups(:)
   INTEGER,                        INTENT    ( IN    )  :: numJoints                ! number of joints specified in the inputs
   TYPE(Morison_JointType),        INTENT    ( IN    )  :: joints(:)                ! array of input joint data structures
   INTEGER,                        INTENT    ( IN    )  :: numMembers               ! number of members specified in the inputs
   TYPE(Morison_MemberInputType),  INTENT    ( INOUT )  :: members(:)               ! array of input member data structures
   INTEGER,                        INTENT    (   OUT )  :: maxNodes                 ! maximum number of nodes which may appear in the final simulation mesh
   INTEGER,                        INTENT    (   OUT )  :: maxElements              ! maximum number of elements which may appear in the final simulation mesh
   INTEGER,                        INTENT    (   OUT )  :: maxSuperMembers          ! maximum number of super members which may appear in the final simulation mesh
   
      ! Local variables
   INTEGER                                              :: WtrSplitCount = 0         ! number of possible new members due to splitting at water boundaries   
   INTEGER                                              :: MGsplitCount = 0         ! number of possible new members due to splitting at marine growth boundaries
   INTEGER                                              :: maxSubMembers = 0        ! maximum added nodes and members due to member subdivision
   INTEGER                                              :: maxSuperMemNodes = 0     ! maximum number of new nodes due to super member generation
   INTEGER                                              :: I, J, j1, j2             ! generic integer for counting
   TYPE(Morison_JointType)                           :: joint1, joint2           ! joint data structures                               
   Real(ReKi)                                           :: z1, z2                   ! z values of the first and second joints
   Real(ReKi)                                           :: temp                     ! temporary variable
   REAL(ReKi)                                           :: memLen                   ! member length
   INTEGER                                              :: nSplits, totalSplits, nodeSplits
   REAL(ReKi)                                           :: possibleSplits(5)
      ! Initialize quantities
   maxNodes         = numJoints
   maxElements      = numMembers
   maxSuperMembers  = 0
   maxSuperMemNodes = 0
   maxSubMembers    = 0
   MGsplitCount     = 0
   WtrSplitCount    = 0
   nodeSplits       = 0
   totalSplits      = 0 
       
      ! Determine new members and nodes due to internal member subdivision
   DO I = 1,numMembers
       
              
      z1 = joints( members(I)%MJointID1Indx )%JointPos(3)
      z2 = joints( members(I)%MJointID2Indx )%JointPos(3)
      IF ( z1 > z2) THEN
         temp = z1
         z1   = z2
         z2   = temp
      END IF
      
      
      
         ! For this member, determine possible split conditions due to crossing through:
         ! MSL, seabed, marine growth boundaries, filled free surface location.
         !
         
      nSplits = 0  
      possibleSplits = -9999999.0  ! Initialize possibleSplit values to a number that never appears in the geometry.
      
         ! Is the member filled?
      IF ( members(I)%MmbrFilledIDIndx /= -1 ) THEN
         nSplits =  1
            ! The free surface is specified relative to the MSL.
         possibleSplits(1) = filledGroups(members(I)%MmbrFilledIDIndx)%FillFSLoc 
      END IF
      
      
         ! Check if MSL is equal to Zfs, if it is, then don't add MSL2SWL as an additional possible split, otherwise do add it.
     
         IF ( .NOT. IsThisSplitValueUsed(nSplits, possibleSplits, MSL2SWL )) THEN
            nSplits = nSplits + 1
            possibleSplits(nSplits) = MSL2SWL
         END IF  
      
      
        ! Is there a marine growth region?
        
      IF ( numMGDepths > 0 ) THEN   
         
            ! Recursively check to see if this
            IF ( .NOT. IsThisSplitValueUsed(nSplits, possibleSplits, MGTop) ) THEN
               nSplits = nSplits + 1
               possibleSplits(nSplits) = MGTop
            END IF
            IF ( .NOT. IsThisSplitValueUsed(nSplits, possibleSplits, MGBottom) ) THEN
               nSplits = nSplits + 1
               possibleSplits(nSplits) = MGBottom
            END IF
         
      END IF
      
        ! Check if seabed is equal to other possibleSplits
      
         IF ( .NOT. IsThisSplitValueUsed(nSplits, possibleSplits, Zseabed) ) THEN
            nSplits = nSplits + 1
            possibleSplits(nSplits) = Zseabed
         END IF  
     
         
       ! Now determine which possible splits this member actually crosses
       
      DO J=1,nSplits
         
         IF ( z1 < possibleSplits(J) .AND. z2 > possibleSplits(J) ) THEN
            members(I)%NumSplits = members(I)%NumSplits + 1
            members(I)%Splits(members(I)%NumSplits) = possibleSplits(J)
         END IF
      
      END DO
         ! Sort the splits from smallest Z value to largest Z value
      CALL BSortReal ( members(I)%Splits, members(I)%NumSplits )
      totalSplits = totalSplits + members(I)%NumSplits
      
      !   ! Determine new members due to elements crossing the MSL or the seabed
      !IF ( z2 > MSL2SWL ) THEN
      !   IF ( z1 < MSL2SWL .AND. z1 >= Zseabed ) THEN
      !      ! Split this member
      !      WtrSplitCount = WtrSplitCount + 1
      !      members(I).WtrSplitState = 1
      !   END IF
      !   IF ( z1 < Zseabed ) THEN
      !      ! Split this member twice because it crosses both boundaries
      !      WtrSplitCount = WtrSplitCount + 2
      !      members(I).WtrSplitState = 3
      !   END IF  
      !END IF
      !IF ( z2 < MSL2SWL .AND. z2 >= Zseabed ) THEN
      !   IF ( z1 < MGBottom ) THEN
      !      ! Split this member
      !      WtrSplitCount = WtrSplitCount + 1
      !      members(I).WtrSplitState = 2
      !   END IF
      !         
      !END IF
      !      
      !   ! Determine new members and nodes due to marine growth boundary splitting
      !   members(I).MGSplitState = 0
      !IF ( numMGDepths > 0 ) THEN
      !   
      !   IF ( z2 > MGTop ) THEN
      !      IF ( z1 < MGTop .AND. z1 >= MGBottom ) THEN
      !         ! Split this member
      !         MGsplitCount = MGsplitCount + 1
      !         members(I).MGSplitState = 1
      !      END IF
      !      IF ( z1 < MGBottom ) THEN
      !         ! Split this member twice because it crosses both boundaries
      !         MGsplitCount = MGsplitCount + 2
      !         members(I).MGSplitState = 3
      !      END IF  
      !   END IF
      !   IF ( z2 < MGTop .AND. z2 >= MGBottom ) THEN
      !      IF ( z1 < MGBottom ) THEN
      !         ! Split this member
      !         MGsplitCount = MGsplitCount + 1
      !         members(I).MGSplitState = 2
      !      END IF
      !         
      !   END IF
      !           
      !END IF
      
      j1 = members(I)%MJointID1Indx
      j2 = members(I)%MJointID2Indx
      joint1 = joints(j1)   ! note Inspector complains of uninitialized variables; this is due to copying types here (and some fields haven't been initialized)
      joint2 = joints(j2)   ! note Inspector complains of uninitialized variables; this is due to copying types here (and some fields haven't been initialized)
      CALL GetDistance(joint1%JointPos, joint2%JointPos, memLen)
      maxSubMembers = maxSubMembers + CEILING( memLen / members(I)%MDivSize  ) - 1
      
   END DO
   
      ! Look for all possible super member creation
   DO I = 1,numJoints
            
         ! Check #1 are there more than 2 members connected to the joint?
      IF ( joints(I)%JointOvrlp == 1 .AND. joints(I)%NConnections > 2) THEN
            maxSuperMemNodes = maxSuperMemNodes + ( joints(I)%NConnections - 1 )
            maxSuperMembers  = maxSuperMembers  + 1  
      ELSE
         nodeSplits = nodeSplits + joints(I)%NConnections - 1
      END IF
            
            
   END DO
   
   maxNodes        = maxNodes    + totalSplits*2 +  nodeSplits + maxSubMembers + maxSuperMemNodes
   maxElements     = maxElements + totalSplits + maxSubMembers
   
   
END SUBROUTINE GetMaxSimQuantities

SUBROUTINE WriteSummaryFile( UnSum, MSL2SWL, WtrDpth, numNodes, nodes, numElements, elements, NOutputs, OutParam, NMOutputs, MOutLst, distribToNodeIndx, NJOutputs, JOutLst, inLumpedMesh, outLumpedMesh, inDistribMesh, outDistribMesh, L_F_B, L_F_BF, D_F_B, D_F_BF, D_F_MG, g, ErrStat, ErrMsg )  !, numDistribMarkers, distribMarkers, numLumpedMarkers, lumpedMarkers

   REAL(ReKi),               INTENT ( IN    )  :: MSL2SWL
   REAL(ReKi),               INTENT ( IN    )  :: WtrDpth
   INTEGER,                  INTENT ( IN    )  :: UnSum
   INTEGER,                  INTENT ( IN    )  :: numNodes
   TYPE(Morison_NodeType),   INTENT ( IN    )  :: nodes(:)  
   INTEGER,                  INTENT ( IN    )  :: numElements
   TYPE(Morison_MemberType), INTENT ( IN    )  :: elements(:)
   INTEGER,                  INTENT ( IN    )  :: NOutputs
   TYPE(OutParmType),        INTENT ( IN    )  :: OutParam(:)
   INTEGER,                  INTENT ( IN    )  :: NMOutputs
   TYPE(Morison_MOutput),    INTENT ( IN    )  :: MOutLst(:)
   INTEGER,                  INTENT ( IN    )  :: distribToNodeIndx(:)
   INTEGER,                  INTENT ( IN    )  :: NJOutputs
   TYPE(Morison_JOutput),    INTENT ( IN    )  :: JOutLst(:)
   TYPE(MeshType),           INTENT ( INOUT )  :: inLumpedMesh
   TYPE(MeshType),           INTENT ( INOUT )  :: outLumpedMesh
   TYPE(MeshType),           INTENT ( INOUT )  :: inDistribMesh
   TYPE(MeshType),           INTENT ( INOUT )  :: outDistribMesh
   REAL(ReKi),               INTENT ( IN    )  :: L_F_B(:,:)           ! Lumped buoyancy force associated with the member
   REAL(ReKi),               INTENT ( IN    )  :: L_F_BF(:,:)          ! Lumped buoyancy force associated flooded/filled fluid within the member
   REAL(ReKi),               INTENT ( IN    )  :: D_F_B(:,:)           ! Lumped buoyancy force associated with the member
   REAL(ReKi),               INTENT ( IN    )  :: D_F_BF(:,:)          ! Lumped buoyancy force associated flooded/filled fluid within the member
   REAL(ReKi),               INTENT ( IN    )  :: D_F_MG(:,:)
   REAL(ReKi),               INTENT ( IN    )  :: g                    ! gravity
   !INTEGER,                  INTENT ( IN    )  :: numDistribMarkers
   !TYPE(Morison_NodeType),   INTENT ( IN    )  :: distribMarkers(:)
   !INTEGER,                  INTENT ( IN    )  :: numLumpedMarkers
   !TYPE(Morison_NodeType),   INTENT ( IN    )  :: lumpedMarkers(:)
   INTEGER,                  INTENT (   OUT )  :: ErrStat             ! returns a non-zero value when an error occurs  
   CHARACTER(*),             INTENT (   OUT )  :: ErrMsg              ! Error message if ErrStat /= ErrID_None

   INTEGER                                     :: I, J
   REAL(ReKi)                                  :: l                   ! length of an element
   LOGICAL                                     :: filledFlag          ! flag indicating if element is filled/flooded
   CHARACTER(2)                                :: strFmt
   CHARACTER(10)                               :: strNodeType         ! string indicating type of node: End, Interior, Super
   REAL(ReKi)                                  :: ident(3,3)          ! identity matrix
   REAL(ReKi)                                  :: ExtBuoyancy(6)      ! sum of all external buoyancy forces lumped at (0,0,0)
   REAL(ReKi)                                  :: IntBuoyancy(6)      ! sum of all internal buoyancy forces lumped at (0,0,0)
   REAL(ReKi)                                  :: MG_Wt(6)            ! weight of the marine growth as applied to (0,0,0)
   TYPE(MeshType)                              :: WRP_Mesh            ! mesh representing the WAMIT reference point (0,0,0)
   TYPE(MeshType)                              :: WRP_Mesh_position   ! mesh representing the WAMIT reference point (0,0,0)   (with no displaced position)
   TYPE(MeshMapType)                           :: M_L_2_P             ! Map  Morison Line2 to  WRP_Mesh point
   TYPE(MeshMapType)                           :: M_P_2_P             ! Map  Morison Point to  WRP_Mesh point
   REAL(ReKi)                                  :: elementVol            ! displaced volume of an element
   REAL(ReKi)                                  :: totalDisplVol       ! total displaced volume of the structure
   REAL(ReKi)                                  :: totalVol            ! total volume of structure
   REAL(ReKi)                                  :: MGvolume            ! volume of the marine growth material
   REAL(ReKi)                                  :: totalMGVol          !
   REAL(ReKi)                                  :: totalFillVol        !
   REAL(ReKi)                                  :: elemCentroid(3)     ! location of the element centroid
   REAL(ReKi)                                  :: COB(3)              ! center of buoyancy location in global coordinates
   INTEGER                                     :: m1, m2              ! Indices of the markers which surround the requested output location
   REAL(ReKi)                                  :: s                   ! The linear interpolation factor for the requested location
   REAL(ReKi)                                  :: outloc(3)           ! Position of the requested member output
   INTEGER                                     :: mbrIndx, nodeIndx
   CHARACTER(10)                               :: tmpName
   REAL(ReKi)                                  :: totalFillMass, mass_fill, fillVol
   REAL(ReKi)                                  :: totalMGMass, mass_MG
   TYPE(Morison_NodeType)                      ::  node1, node2
   REAL(ReKi)                                  :: Cd1, Cd2, Ca1, Ca2, Cp1, Cp2, AxCa1, AxCa2, AxCp1, AxCp2, JAxCd1, JAxCd2, JAxCa1, JAxCa2, JAxCp1, JAxCp2 ! tmp coefs
   
      ! Initialize data
   ErrStat       = ErrID_None
   ErrMsg        = ""
   ExtBuoyancy   = 0.0
   totalFillMass = 0.0
   totalDisplVol = 0.0
   totalVol      = 0.0
   totalMGVol    = 0.0
   totalFillVol  = 0.0
   totalMGMass   = 0.0
   COB           = 0.0
   
      ! Create identity matrix
   CALL EYE(ident,ErrStat,ErrMsg)
   
   IF ( UnSum > 0 ) THEN
      
      

      !WRITE( UnSum,  '(//)' ) 
      !WRITE( UnSum,  '(A8)' ) 'Elements'
      !WRITE( UnSum,  '(/)' ) 
      !WRITE( UnSum, '(1X,A5,2X,A5,2X,A5,5(2X,A12),2X,A12,17(2X,A12))' ) '  i  ', 'node1','node2','  Length  ', '  MGVolume  ', '  MGDensity ', 'PropPot ', 'FilledFlag', 'FillDensity', '  FillFSLoc ', '  FillMass  ', '     Cd1    ', '   CdMG1  ', '     Ca1    ', '    CaMG1   ', '      R1    ', '     t1     ','     Cd2    ', '    CdMG2   ', '     Ca2    ', '    CaMG2   ', '      R2    ', '     t2     '
      !WRITE( UnSum, '(1X,A5,2X,A5,2X,A5,5(2X,A12),2X,A12,17(2X,A12))' ) ' (-) ', ' (-) ',' (-) ','   (m)    ', '   (m^3)    ', '  (kg/m^3)  ', '   (-)    ', '   (-)    ', ' (kg/m^3)  ', '     (-)    ', '    (kg)    ', '     (-)    ', '    (-)   ', '     (-)    ', '     (-)    ', '     (m)    ', '     (m)    ','     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (m)    ', '     (m)    '
      !
      
      DO I = 1,numElements 
         
         node1   = nodes(elements(I)%Node1Indx)
         node2   = nodes(elements(I)%Node2Indx)
         IF ( ( (node1%tMG > 0.0_ReKi ) .AND. EqualRealNos(node2%tMG,0.0_ReKi) ) .OR. ( (node2%tMG > 0.0_ReKi ) .AND. EqualRealNos(node1%tMG,0.0_ReKi) ) ) THEN
            ErrStat = ErrID_Fatal
            ErrMsg  = 'If one node of an element has MG, then both must.  This is an internal code problem within HydroDyn.'
            RETURN
         END IF
         CALL GetDistance( nodes(elements(I)%Node1Indx)%JointPos, nodes(elements(I)%Node2Indx)%JointPos, l )
         
         
         IF (elements(I)%PropPot) THEN
            MGvolume    = 0.0
            elementVol  = 0.0
         ELSE
            elementVol  = ElementVolume(elements(I)%R1 + node1%tMG, elements(I)%R2 + node2%tMG, l)
            MGvolume    = elementVol  - ElementVolume(elements(I)%R1, elements(I)%R2, l)
         END IF
         totalMGVol  = totalMGVol  + MGvolume
         mass_MG     = MGvolume*elements(I)%FillDens
         totalMGMass = totalMGMass + mass_MG
         CALL ElementCentroid(elements(I)%R1 + node1%tMG, elements(I)%R2 + node2%tMG, node1%JointPos, l, elements(I)%R_LToG, elemCentroid)
         
         COB         = COB         + elementVol*elemCentroid
         
         totalVol    = totalVol    + elementVol
         
         IF ( node2%JointPos(3) <= MSL2SWL .AND. node1%JointPos(3) >= -WtrDpth) totalDisplVol = totalDisplVol + elementVol
         
         IF ( elements(I)%MmbrFilledIDIndx > 0 ) THEN          
           ! filledFlag = .TRUE.
            !IF ( ( node2%JointPos(3) <= elements(I)%FillFSLoc ) .AND. ( node1%JointPos(3) <= elements(I)%FillFSLoc ) ) THEN
               fillVol       = ElementVolume(elements(I)%R1 - elements(I)%t1, elements(I)%R2 - elements(I)%t2, l)
               totalFillVol  = totalFillVol  + fillVol
               mass_fill     = elements(I)%FillDens*fillVol
               totalFillMass = totalFillMass + mass_fill
            !END IF
         ELSE
           ! mass_fill  = 0.0
           ! filledFlag = .FALSE.
         END IF
         
         !WRITE( UnSum, '(1X,I5,2X,I5,2X,I5,3(2X,ES12.5),2(2X,L12),2X,ES12.5,17(2X,ES12.5))' ) I, elements(I)%Node1Indx, elements(I)%Node2Indx, l, MGvolume, node1%MGdensity, elements(I)%PropPot, filledFlag, elements(I)%FillDens, elements(I)%FillFSLoc, mass_fill, elements(I)%Cd1, elements(I)%CdMG1, elements(I)%Ca1, elements(I)%CaMG1, elements(I)%R1, elements(I)%t1, elements(I)%Cd2, elements(I)%CdMG2, elements(I)%Ca2, elements(I)%CaMG2, elements(I)%R2, elements(I)%t2

      END DO   ! I = 1,numElements 
               
      
     
      
      WRITE( UnSum,  '(//)' ) 
      WRITE( UnSum, '(A24)' )        'Volume Calculations(m^3)'
      WRITE( UnSum, '(A24)' )        '------------------------'
      WRITE( UnSum, '(A27,ES12.5)' ) '  Structure Volume     :   ', totalVol
      WRITE( UnSum, '(A27,ES12.5)' ) '  Submerged Volume     :   ', totalDisplVol
      WRITE( UnSum, '(A27,ES12.5)' ) '  Marine Growth Volume :   ', totalMGVol   
      WRITE( UnSum, '(A27,ES12.5)' ) '  Ballasted Volume     :   ', totalFillVol
      WRITE( UnSum, '(A111)') '              NOTE: Structure, Submerged and Marine Growth volumes are based on members not modelled with WAMIT'
      WRITE( UnSum, '(A149)') '                    Ballasted volume is computed from all members which are marked as filled in the HydroDyn input file, regardless of PropPot flag'
           
      
      
         ! Sum all buoyancy loads to the COB
         ! Do this by creating a temporary mesh which is for (0,0,0)
         
      !COB = COB / totalVol   
      
         ! Write out the Center of Buoyancy (geometric center of the displaced volume)
      !WRITE( UnSum,  '(//)' ) 
      !WRITE( UnSum, '(A18)' )        'Center of Buoyancy'
      !WRITE( UnSum, '(3(2X,A10  ))' ) ' COBXi ', ' COBYi ', ' COBZi '
      !WRITE( UnSum, '(3(2X,A10  ))' ) '  (m)  ', '  (m)  ', '  (m)  '
      !WRITE( UnSum, '(3(2X,F10.3))' ) COB(1)   , COB(2)   , COB(3)
      
      CALL MeshCreate( BlankMesh        = WRP_Mesh          &
                     ,IOS               = COMPONENT_INPUT   &
                     ,Nnodes            = 1                 &
                     ,ErrStat           = ErrStat           &
                     ,ErrMess           = ErrMsg            &
                     ,Force             = .TRUE.            &
                     ,Moment            = .TRUE.            &
                     )
         ! Create the node on the mesh
            
      CALL MeshPositionNode (WRP_Mesh                              &
                              , 1                                  &
                              , (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/)   &  
                              , ErrStat                            &
                              , ErrMsg                             &
                              )
      
      IF ( ErrStat /= 0 ) RETURN
       
      
         ! Create the mesh element
      CALL MeshConstructElement (  WRP_Mesh            &
                                  , ELEMENT_POINT      &                         
                                  , ErrStat            &
                                  , ErrMsg             &
                                  , 1                  &
                                )
      CALL MeshCommit ( WRP_Mesh           &
                      , ErrStat            &
                      , ErrMsg             )
   
      IF ( ErrStat /= ErrID_None ) RETURN
            
         ! we need the translation displacement mesh for loads transfer:
      CALL MeshCopy ( SrcMesh  = WRP_Mesh            &
                    , DestMesh = WRP_Mesh_position   &
                    , CtrlCode = MESH_SIBLING        &
                    , IOS      = COMPONENT_INPUT     &
                    , TranslationDisp = .TRUE.       &
                    , ErrStat  = ErrStat             &
                    , ErrMess  = ErrMsg              )  ! automatically sets    DestMesh%RemapFlag = .TRUE.
                    
      IF ( ErrStat /= ErrID_None ) RETURN
      WRP_Mesh_position%TranslationDisp = 0.0  ! bjj: this is actually initialized in the ModMesh module, but I'll do it here anyway.
      
      WRP_Mesh%RemapFlag  = .TRUE.
      
      
         ! Attach the external distributed buoyancy loads to the distributed mesh so they can be transferred to the WRP
      
         ! Because of wave stretching and user-supplied waves, we may have loads above the still water line (SWL) which will be used
         ! in the hydrodynamics for conditions where the wave height is > SWL.  So we now need to check that the vertical position
         ! is <= SWL for this summary file calculation.
      
      DO J = 1, outDistribMesh%Nnodes
         if ( outDistribMesh%Position(3,J) <= MSL2SWL ) then
            DO I=1,6
            
               IF (I < 4 ) THEN
               
                  outDistribMesh%Force(I   ,J) = D_F_B(I,J)
            
               ELSE
               
                  outDistribMesh%Moment(I-3,J) = D_F_B(I,J)
               
               END IF
            
            END DO  ! DO I
         end if              ! <= MSL2SWL check
      END DO ! DO J
      
 
         ! Transfer the loads from the distributed mesh to the (0,0,0) point mesh
         
      CALL MeshMapCreate           ( outDistribMesh, WRP_Mesh, M_L_2_P, ErrStat, ErrMsg                )
        !CALL CheckError( ErrStat, 'Message from MeshMapCreate HD_M_L_2_ED_P: '//NewLine//ErrMsg )
      CALL Transfer_Line2_to_Point( outDistribMesh, WRP_Mesh, M_L_2_P, ErrStat, ErrMsg, inDistribMesh, WRP_Mesh_position )
      
      ExtBuoyancy(1:3) = WRP_Mesh%Force (:,1)
      ExtBuoyancy(4:6) = WRP_Mesh%Moment(:,1)
    
      
      
         ! Transfer the loads from the lumped mesh to the (0,0,0) point mesh

         ! Because of wave stretching and user-supplied waves, we may have loads above the still water line (SWL) which will be used
         ! in the hydrodynamics for conditions where the wave height is > SWL.  So we now need to check that the vertical position
         ! is <= SWL for this summary file calculation.

      DO J = 1, outLumpedMesh%Nnodes
         if ( outLumpedMesh%Position(3,J) <= MSL2SWL ) then 
            DO I=1,6
            
               IF (I < 4 ) THEN           
               
                  outLumpedMesh%Force(I   ,J) = L_F_B(I,J)
            
               ELSE
               
                  outLumpedMesh%Moment(I-3,J) = L_F_B(I,J)
               
               END IF
            
            END DO  ! DO I
         end if              ! <= MSL2SWL check
      END DO ! DO J
      
         ! Remap for the lumped to WRP mesh transfer       
      WRP_Mesh%RemapFlag  = .TRUE.
      
      CALL MeshMapCreate           ( outLumpedMesh, WRP_Mesh, M_P_2_P, ErrStat, ErrMsg               )
      CALL Transfer_Point_to_Point( outLumpedMesh, WRP_Mesh, M_P_2_P, ErrStat, ErrMsg, inLumpedMesh, WRP_Mesh_position )
      
      ExtBuoyancy(1:3) = ExtBuoyancy(1:3) + WRP_Mesh%Force (:,1)
      ExtBuoyancy(4:6) = ExtBuoyancy(4:6) + WRP_Mesh%Moment(:,1)
      
      
      
         ! Write the buoyancy table headers and the external results

      WRITE( UnSum,  '(//)' ) 
      WRITE( UnSum, '(A45)' ) 'Buoyancy loads summed about ( 0.0, 0.0, 0.0 )'
      WRITE( UnSum, '(A45)' ) '---------------------------------------------'
      WRITE( UnSum, '(18x,6(2X,A20))' ) ' BuoyFxi ', ' BuoyFyi ', ' BuoyFzi ', ' BuoyMxi ', ' BuoyMyi ', ' BuoyMzi '
      WRITE( UnSum, '(18x,6(2X,A20))' ) '   (N)   ', '   (N)   ', '   (N)   ', '  (N-m)  ', '  (N-m)  ', '  (N-m)  '
      WRITE( UnSum, '(A18,6(2X,ES20.6))') 'External:        ', ExtBuoyancy(1), ExtBuoyancy(2), ExtBuoyancy(3), ExtBuoyancy(4), ExtBuoyancy(5), ExtBuoyancy(6)
      
      
         ! Now compute internal Buoyancy
         
      DO J = 1, outDistribMesh%Nnodes
         
         DO I=1,6
            
            IF (I < 4 ) THEN
               
               outDistribMesh%Force(I,J   ) = D_F_BF(I,J)
               
            ELSE
               
               outDistribMesh%Moment(I-3,J) = D_F_BF(I,J)
               
            END IF
            
         END DO  ! DO I
         
      END DO ! DO J
       
      IntBuoyancy = 0.0
      CALL Transfer_Line2_to_Point( outDistribMesh, WRP_Mesh, M_L_2_P, ErrStat, ErrMsg, inDistribMesh, WRP_Mesh_position )
      IntBuoyancy(1:3) = WRP_Mesh%Force(:,1)
      IntBuoyancy(4:6) = WRP_Mesh%Moment(:,1)
      
      
      DO J = 1, outLumpedMesh%Nnodes
         
         DO I=1,6
            
            IF (I < 4 ) THEN
               
               outLumpedMesh%Force(I,J) = L_F_BF(I,J)
            
            ELSE
               
               outLumpedMesh%Moment(I-3,J) = L_F_BF(I,J)
               
            END IF
            
         END DO  ! DO I
         
      END DO ! DO J 
      
      CALL Transfer_Point_to_Point( outLumpedMesh, WRP_Mesh, M_P_2_P, ErrStat, ErrMsg, inLumpedMesh, WRP_Mesh_position )
      IntBuoyancy(1:3) = IntBuoyancy(1:3) + WRP_Mesh%Force(:,1)
      IntBuoyancy(4:6) = IntBuoyancy(4:6) + WRP_Mesh%Moment(:,1)
      
      
         ! clean up
      
      CALL MeshMapDestroy( M_P_2_P, ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) CALL WrScr(TRIM(ErrMsg))
     
      WRITE( UnSum, '(A18,6(2X,ES20.6))') 'Internal:        ', IntBuoyancy(1), IntBuoyancy(2), IntBuoyancy(3), IntBuoyancy(4), IntBuoyancy(5), IntBuoyancy(6)
      IntBuoyancy = IntBuoyancy + ExtBuoyancy
      WRITE( UnSum, '(A18,6(2X,ES20.6))') 'Total   :        ', IntBuoyancy(1), IntBuoyancy(2), IntBuoyancy(3), IntBuoyancy(4), IntBuoyancy(5), IntBuoyancy(6)
      !WRITE( UnSum,  '(/)' ) 
      WRITE( UnSum, '(A81)') '              NOTE: External buoyancy is based on members not modelled with WAMIT'
      WRITE( UnSum, '(A150)') '                    Internal buoyancy is computed from all members which are marked as filled in the HydroDyn input file, regardless of PropPot flag'
      WRITE( UnSum, '(A88)') '                    Total buoyancy does not include WAMIT-modelled buoyancy contribution'
      
      
      
         ! Now compute marine growth weight at the WRP
         
      DO J = 1, outDistribMesh%Nnodes
         
         DO I=1,6
            
            IF (I < 4 ) THEN           
               
               outDistribMesh%Force(I,J) = D_F_MG(I,J) 
            
            ELSE
               
               outDistribMesh%Moment(I-3,J) = D_F_MG(I,J)
               
            END IF
            
         END DO  ! DO I
         
       END DO ! DO J
       
         
      MG_Wt = 0.0
      CALL Transfer_Line2_to_Point( outDistribMesh, WRP_Mesh, M_L_2_P, ErrStat, ErrMsg, inDistribMesh, WRP_Mesh_position )
      MG_Wt(1:3) = WRP_Mesh%Force(:,1)
      MG_Wt(4:6) = WRP_Mesh%Moment(:,1)
      
           
      WRITE( UnSum,  '(//)' ) 
      WRITE( UnSum, '(A36)' ) 'Weight loads about ( 0.0, 0.0, 0.0 )'
      WRITE( UnSum, '(A36)' ) '------------------------------------'
      WRITE( UnSum, '(18x,6(2X,A20))' ) '  MGFxi  ', '  MGFyi  ', '  MGFzi  ', '  MGMxi  ', '  MGMyi  ', '  MGMzi  '
      WRITE( UnSum, '(18x,6(2X,A20))' ) '   (N)   ', '   (N)   ', '   (N)   ', '  (N-m)  ', '  (N-m)  ', '  (N-m)  '
      !WRITE( UnSum, '(A16,6(2X,E20.6))') 'Structure    :  ',  M_Wt(1),  M_Wt(2),  M_Wt(3),  M_Wt(4),  M_Wt(5),  M_Wt(6)
      WRITE( UnSum, '(A18,6(2X,ES20.6))') 'Marine Growth:   ', MG_Wt(1), MG_Wt(2), MG_Wt(3), MG_Wt(4), MG_Wt(5), MG_Wt(6)
      !WRITE( UnSum, '(A16,6(2X,E20.6))') 'Filled Fluid :  ',  F_Wt(1),  F_Wt(2),  F_Wt(3),  F_Wt(4),  F_Wt(5),  F_Wt(6)
      !M_Wt = M_Wt + MG_Wt + F_Wt
      !WRITE( UnSum, '(A16,6(2X,E20.6))') 'Total        :  ',  M_Wt(1),  M_Wt(2),  M_Wt(3),  M_Wt(4),  M_Wt(5),  M_Wt(6)
      
      
      CALL MeshMapDestroy( M_L_2_P, ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) CALL WrScr(TRIM(ErrMsg))
      CALL MeshDestroy(WRP_Mesh, ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) CALL WrScr(TRIM(ErrMsg))
      CALL MeshDestroy(WRP_Mesh_position, ErrStat, ErrMsg ); IF ( ErrStat /= ErrID_None ) CALL WrScr(TRIM(ErrMsg))
      
         ! Write the header for this section
      WRITE( UnSum,  '(//)' ) 
      WRITE( UnSum,  '(A5)' ) 'Nodes'
      WRITE( UnSum,  '(/)' ) 
      WRITE( UnSum, '(1X,A5,24(2X,A10),2X,A5,2X,A15)' ) '  i  ', 'JointIndx ', 'JointOvrlp', 'InpMbrIndx', '   Nxi    ', '   Nyi    ', '   Nzi    ', 'InpMbrDist', '     R    ', '   dRdZ   ', '    t     ', '   tMG    ', '  MGDens  ', ' PropPot ', 'FilledFlag', ' FillDens ', 'FillFSLoc ', '    Cd    ', '    Ca    ', '    Cp    ', '   AxCa   ', '   AxCp   ', '   JAxCd  ', '   JAxCa  ', '   JAxCp  ', 'NConn ', 'Connection List'
      WRITE( UnSum, '(1X,A5,24(2X,A10),2X,A5,2X,A15)' ) ' (-) ', '   (-)    ', '   (-)    ', '   (-)    ', '   (m)    ', '   (m)    ', '   (m)    ', '    (-)   ', '    (m)   ', '    (-)   ', '   (m)    ', '   (m)    ', ' (kg/m^3) ', '   (-)    ', '   (-)    ', ' (kg/m^3) ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', ' (-)  ', '               '

         ! Write the data
      DO I = 1,numNodes   
         WRITE(strFmt,'(I2)') nodes(I)%NConnections
         IF ( nodes(I)%NodeType == 1 ) THEN
            strNodeType = 'End       '
         ELSE IF ( nodes(I)%NodeType == 2 ) THEN
            strNodeType = 'Interior  '
         ELSE IF ( nodes(I)%NodeType == 3 ) THEN
            strNodeType = 'Super     '
         ELSE
            strNodeType = 'ERROR     '
         END IF
         
         WRITE( UnSum, '(1X,I5,3(2X,I10),4(2X,F10.4),5(2X,ES10.3),2(2X,L10),10(2X,ES10.3),2X,I5,' // strFmt // '(2X,I4))' ) I, nodes(I)%JointIndx, nodes(I)%JointOvrlp, nodes(I)%InpMbrIndx, nodes(I)%JointPos, nodes(I)%InpMbrDist, nodes(I)%R, nodes(I)%DRDZ, nodes(I)%t, nodes(I)%tMG, nodes(I)%MGdensity, nodes(I)%PropPot, nodes(I)%FillFlag, nodes(I)%FillDensity, nodes(I)%FillFSLoc, nodes(I)%Cd, nodes(I)%Ca, nodes(I)%Cp, nodes(I)%AxCa, nodes(I)%AxCp, nodes(I)%JAxCd, nodes(I)%JAxCa, nodes(I)%JAxCp, nodes(I)%NConnections, nodes(I)%ConnectionList(1:nodes(I)%NConnections)
      END DO
      
       WRITE( UnSum,  '(//)' ) 
      WRITE( UnSum,  '(A8)' ) 'Elements'
      WRITE( UnSum,  '(/)' ) 
      WRITE( UnSum, '(1X,A5,2X,A5,2X,A5,13(2X,A12),2X,A12,21(2X,A12))' ) '  i  ', 'node1','node2','  Length  ', '   Volume   ', '  MGVolume  ', '      R1    ', '    tMG1    ', '     t1     ', '      R2    ', '    tMG2    ', '     t2     ', '   MGDens1  ', '   MGDens2  ', ' PropPot ', 'FilledFlag', 'FillDensity', '  FillFSLoc ', '  FillMass  ', '     Cd1    ', '    Ca1   ', '     Cp1    ', '    AxCa1   ', '    AxCp1   ', '   JAxCd1   ', '   JAxCa1   ', '  JAxCp1   ', '     Cd2    ', '     Ca2    ', '     Cp2    ', '    AxCa2   ', '    AxCp2   ', '   JAxCd2   ', '   JAxCa2   ', '   JAxCp2   '
      WRITE( UnSum, '(1X,A5,2X,A5,2X,A5,13(2X,A12),2X,A12,21(2X,A12))' ) ' (-) ', ' (-) ',' (-) ','   (m)    ', '   (m^3)    ', '   (m^3)    ', '     (m)    ', '     (m)    ', '     (m)    ', '     (m)    ', '     (m)    ', '     (m)    ', '  (kg/m^3)  ', '  (kg/m^3)  ', '   (-)    ', '   (-)    ', ' (kg/m^3)  ', '     (-)    ', '    (kg)    ', '     (-)    ', '    (-)   ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)   ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    '
      
      
   
   
      DO I = 1,numElements 
         
         node1   = nodes(elements(I)%Node1Indx)
         node2   = nodes(elements(I)%Node2Indx)
         IF ( ( (node1%tMG > 0.0_ReKi ) .AND. EqualRealNos(node2%tMG,0.0_ReKi) ) .OR. ( (node2%tMG > 0.0_ReKi ) .AND. EqualRealNos(node1%tMG,0.0_ReKi) ) ) THEN
            ErrStat = ErrID_Fatal
            ErrMsg  = 'If one node of an element has MG, then both must.  This is an internal code problem within HydroDyn.'
            RETURN
         END IF
         CALL GetDistance( nodes(elements(I)%Node1Indx)%JointPos, nodes(elements(I)%Node2Indx)%JointPos, l )
         
         
         IF (elements(I)%PropPot) THEN
            MGvolume    = 0.0
            elementVol  = 0.0
         ELSE
            elementVol  = ElementVolume(elements(I)%R1 + node1%tMG, elements(I)%R2 + node2%tMG, l)
            MGvolume    = elementVol  - ElementVolume(elements(I)%R1, elements(I)%R2, l)
         END IF
        ! totalMGVol  = totalMGVol  + MGvolume
        ! mass_MG     = MGvolume*elements(I)%FillDens
        ! totalMGMass = totalMGMass + mass_MG
         CALL ElementCentroid(elements(I)%R1 + node1%tMG, elements(I)%R2 + node2%tMG, node1%JointPos, l, elements(I)%R_LToG, elemCentroid)
         
        ! COB         = COB         + elementVol*elemCentroid
         
        ! totalVol    = totalVol    + elementVol
         
         !IF ( node2%JointPos(3) <= MSL2SWL .AND. node1%JointPos(3) >= ) totalDisplVol = totalDisplVol + elementVol
         
         IF ( elements(I)%MmbrFilledIDIndx > 0 ) THEN          
            filledFlag = .TRUE.
            !IF ( ( node2%JointPos(3) <= elements(I)%FillFSLoc ) .AND. ( node1%JointPos(3) <= elements(I)%FillFSLoc ) ) THEN
               fillVol       = ElementVolume(elements(I)%R1 - elements(I)%t1, elements(I)%R2 - elements(I)%t2, l)
              ! totalFillVol  = totalFillVol  + fillVol
               mass_fill     = elements(I)%FillDens*fillVol
              ! totalFillMass = totalFillMass + mass_fill
            !END IF
         ELSE
            mass_fill  = 0.0
            filledFlag = .FALSE.
         END IF
         
         IF (EqualRealNos(node1%tMG,0.0_ReKi)) THEN
            Cd1   = elements(I)%Cd1
            Cd2   = elements(I)%Cd2
            Ca1   = elements(I)%Ca1
            Ca2   = elements(I)%Ca2
            Cp1   = elements(I)%Cp1
            Cp2   = elements(I)%Cp2
            AxCa1 = elements(I)%AxCa1
            AxCa2 = elements(I)%AxCa2
            AxCp1 = elements(I)%AxCp1
            AxCp2 = elements(I)%AxCp2
         ELSE
            Cd1   = elements(I)%CdMG1
            Cd2   = elements(I)%CdMG2
            Ca1   = elements(I)%CaMG1
            Ca2   = elements(I)%CaMG2
            Cp1   = elements(I)%CpMG1
            Cp2   = elements(I)%CpMG2
            AxCa1 = elements(I)%AxCaMG1
            AxCa2 = elements(I)%AxCaMG2
            AxCp1 = elements(I)%AxCpMG1
            AxCp2 = elements(I)%AxCpMG2
         END IF
         
         JAxCd1 = node1%JAxCd
         JAxCa1 = node1%JAxCa
         JAxCp1 = node1%JAxCp
         JAxCd2 = node2%JAxCd
         JAxCa2 = node2%JAxCa
         JAxCp2 = node2%JAxCp
         
         WRITE( UnSum, '(1X,I5,2X,I5,2X,I5,11(2X,ES12.5),2(2X,L12),2X,ES12.5,21(2X,ES12.5))' ) I, &
                       elements(I)%Node1Indx, elements(I)%Node2Indx, l, elementVol, MGvolume, elements(I)%R1, &
                       node1%tMG, elements(I)%t1, elements(I)%R2, node2%tMG, elements(I)%t2, node1%MGdensity, node2%MGdensity, &
                       elements(I)%PropPot, filledFlag, elements(I)%FillDens, elements(I)%FillFSLoc, &
                       mass_fill, Cd1, Ca1, Cp1, AxCa1, AxCp1, JAxCd1, JAxCa1, JAxCp1, &
                       Cd2, Ca2, Cp2, AxCa2, AxCp2, JAxCd2, JAxCa2, JAxCp2

      END DO   ! I = 1,numElements 
               
      
      WRITE( UnSum,  '(//)' ) 
      WRITE( UnSum,  '(A24)' ) 'Requested Member Outputs'
      WRITE( UnSum,  '(/)' ) 
      WRITE( UnSum, '(1X,A10,11(2X,A10))' ) '  Label   ', '    Xi    ',  '    Yi    ', '    Zi    ', 'InpMbrIndx', ' StartXi  ',  ' StartYi  ', ' StartZi  ', '  EndXi   ', '  EndYi   ', '  EndZi   ', '   Loc    '
      WRITE( UnSum, '(1X,A10,11(2X,A10))' ) '   (-)    ', '    (m)   ',  '    (m)   ', '    (m)   ', '   (-)    ', '   (m)    ',  '   (m)    ', '   (m)    ', '   (m)    ', '   (m)    ', '   (m)    ', '   (-)    '
      
      
      DO I = 1,NOutputs
     ! DO J=1, NMOutputs     
         !DO I=1, MOutLst(J)%NOutLoc   
         
           
               ! Get the member index and node index for this output label.  If this is not a member output the indices will return 0 with no errcode.
           ! CALL MrsnOut_GetMemberOutputInfo(WriteOutputHdr(I), NMOutputs, MOutLst, mbrIndx, nodeIndx, ErrStat, ErrMsg )
          !  IF (ErrStat >= AbortErrLev ) RETURN
           ! IF ( mbrIndx > 0 ) THEN
         tmpName =  OutParam(I)%Name
         IF (OutParam(I)%SignM == -1 ) tmpName = tmpName(2:10)
               
         IF ( ( INDEX( 'mM', tmpName(1:1) ) > 0 ) .AND. (OutParam(I)%Units /= 'INVALID' ) ) THEN
               !Get Member index and Node index
            read (tmpName(2:2),*) mbrIndx
            read (tmpName(4:4),*) nodeIndx
            
            ! These indices are in the DistribMesh index system, not the overal nodes index system, so distribToNodeIndx() mapping needs to be performed if you want 
            !   to index into the nodes array or wave kinematics arrays
            
            m1 = MOutLst(mbrIndx)%Marker1(nodeIndx)
            m2 = MOutLst(mbrIndx)%Marker2(nodeIndx)
            s  = MOutLst(mbrIndx)%s      (nodeIndx)
         
               ! The member output is computed as a linear interpolation of the nearest two markers
            node1 = nodes(distribToNodeIndx((m1)))
            node2 = nodes(distribToNodeIndx((m2)))
            
            outLoc    = node1%JointPos*(1-s) + node2%JointPos*s
            WRITE( UnSum, '(1X,A10,3(2x,F10.4),2x,I10,7(2x,F10.4))' ) OutParam(I)%Name, outLoc, node1%InpMbrIndx, node1%JointPos, node2%JointPos, s
         END IF
         
          !  END IF 
           !WRITE( UnSum, '(1X,A10,11(2X,ES10.3))' ) WriteOutputHdr(I)
        ! END DO      
      END DO
      
      
      WRITE( UnSum,  '(//)' ) 
      WRITE( UnSum,  '(A24)' ) 'Requested Joint Outputs'
      WRITE( UnSum,  '(/)' ) 
      WRITE( UnSum, '(1X,A10,5(2X,A10))' ) '  Label   ', '    Xi    ',  '    Yi    ', '    Zi    ', 'InpJointID'
      WRITE( UnSum, '(1X,A10,5(2X,A10))' ) '   (-)    ', '    (m)   ',  '    (m)   ', '    (m)   ', '   (-)    '
      
      
      DO I = 1,NOutputs
     ! DO J=1, NMOutputs     
         !DO I=1, MOutLst(J)%NOutLoc   
         
           
               ! Get the member index and node index for this output label.  If this is not a member output the indices will return 0 with no errcode.
           ! CALL MrsnOut_GetMemberOutputInfo(WriteOutputHdr(I), NMOutputs, MOutLst, mbrIndx, nodeIndx, ErrStat, ErrMsg )
          !  IF (ErrStat >= AbortErrLev ) RETURN
           ! IF ( mbrIndx > 0 ) THEN
         tmpName =  OutParam(I)%Name
         IF (OutParam(I)%SignM == -1 ) tmpName = tmpName(2:10)
               
         IF ( ( INDEX( 'jJ', tmpName(1:1) ) > 0 ) .AND. (OutParam(I)%Units /= 'INVALID') ) THEN
            
               !Get Member index and Node index
            read (tmpName(2:2),*) nodeIndx
            m1 = JOutLst(nodeIndx)%Markers(1)     
            WRITE( UnSum, '(1X,A10,3(2x,F10.4),2x,I10)' ) OutParam(I)%Name, nodes(m1)%JointPos, JOutLst(nodeIndx)%JointID
            
         END IF
         
          !  END IF 
           !WRITE( UnSum, '(1X,A10,11(2X,ES10.3))' ) WriteOutputHdr(I)
        ! END DO      
      END DO
      
   END IF

END SUBROUTINE WriteSummaryFile

!====================================================================================================
SUBROUTINE SplitElementOnZBoundary( axis, boundary, iCurrentElement, numNodes, numElements, node1, node2, originalElement, newNode, newElement, ErrStat, ErrMsg )

   INTEGER,                  INTENT ( IN    )  :: axis
   REAL(ReKi),               INTENT ( IN    )  :: boundary
   INTEGER,                  INTENT ( IN    )  :: iCurrentElement
   INTEGER,                  INTENT ( INOUT )  :: numNodes
   TYPE(Morison_NodeType),   INTENT ( INOUT )  :: node1, node2   
   INTEGER,                  INTENT ( INOUT )  :: numElements
   TYPE(Morison_MemberType), INTENT ( INOUT )  :: originalElement
   TYPE(Morison_NodeType),   INTENT (   OUT )  :: newNode
   TYPE(Morison_MemberType), INTENT (   OUT )  :: newElement
   INTEGER,                  INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),             INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None

   INTEGER                                     :: I, J
   REAL(ReKi)                                  :: s
   INTEGER                                     :: newNodeIndx, newElementIndx
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
      ! Create new node and element indices
   newNodeIndx    = numNodes + 1
   newElementIndx = numElements + 1
   
      ! find normalized distance from 1nd node to the boundary
   CALL FindInterpFactor( boundary, node1%JointPos(axis), node2%JointPos(axis), s )
   newNode = node1 ! copy all base node properties
   DO I=axis,axis+1
      J = MOD(I,3) + 1
      newNode%JointPos(J) =  node1%JointPos(J)*(1-s) + node2%JointPos(J)*s
   END DO
   newNode%JointPos(axis) =  boundary
   newNode%R_LToG         =  node1%R_LToG   
      ! Create the new  node information.  
      ! Note that the caller will determine if this is an interior node (subdivide) or an endnode (split due to MG, MSL, seabed, filled free surface)
   newNode%JointOvrlp = 0
   newNode%NConnections = 2
   newNode%ConnectionList(1) = iCurrentElement
   newNode%ConnectionList(2) = newElementIndx
   
   
      
      ! Update node2 connectivity
   DO I = 1,10  ! 10 is the maximum number of connecting elements for a node, this should probably be a parameter !! TODO
      IF ( node2%ConnectionList(I) == iCurrentElement ) THEN
         node2%ConnectionList(I) = newElementIndx
         EXIT
      END IF
   END DO
   
   
      ! Create the new element properties by first copying all the properties from the existing element
   newElement = originalElement
      ! Linearly interpolate the coef values based on depth
   originalElement%R2          = originalElement%R1 * (1-s) + originalElement%R2*s 
   newElement%R1               = originalElement%R2 
   originalElement%t2          = originalElement%t1 * (1-s) + originalElement%t2*s 
   originalElement%InpMbrDist2 = originalElement%InpMbrDist1 * (1-s) + originalElement%InpMbrDist2*s 
   newElement%t1               = originalElement%t2 
   newElement%InpMbrDist1      = originalElement%InpMbrDist2 
   
      ! The end point of the new element is set to the original end point of the existing element, then
      ! the starting point of the new element and the ending point of the existing element are set to the 
      ! newly created node
   newElement%Node2Indx      = originalElement%Node2Indx
   originalElement%Node2Indx = newNodeIndx        
   newElement%Node1Indx      = newNodeIndx
   
END SUBROUTINE SplitElementOnZBoundary


!====================================================================================================
!SUBROUTINE SplitElementsForMG(MGTop, MGBottom, numNodes, nodes, numElements, elements, ErrStat, ErrMsg)   
!
!   REAL(ReKi),               INTENT ( IN    )  :: MGTop
!   REAL(ReKi),               INTENT ( IN    )  :: MGBottom
!   INTEGER,                  INTENT ( INOUT )  :: numNodes
!   TYPE(Morison_NodeType),   INTENT ( INOUT )  :: nodes(:)   
!   INTEGER,                  INTENT ( INOUT )  :: numElements
!   TYPE(Morison_MemberType), INTENT ( INOUT )  :: elements(:)
!   INTEGER,                  INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs  
!   CHARACTER(*),             INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None
!   
!   INTEGER                                     :: I, J, K
!   INTEGER                                     :: node1Indx, node2Indx
!   TYPE(Morison_NodeType)                      :: node1, node2, newNode, newNode2
!   TYPE(Morison_MemberType)                    :: element, newElement, newElement2
!   REAL(ReKi)                                  :: zBoundary
!   INTEGER                                     :: origNumElements
!   
!   origNumElements = numElements
!   
!   DO I=1,origNumElements
!     
!      IF ( elements(I)%MGSplitState > 0 ) THEN
!         
!         node1Indx =  elements(I)%Node1Indx
!         node1 = nodes(node1Indx)
!         node2Indx =  elements(I)%Node2Indx
!         node2 = nodes(node2Indx)
!         element = elements(I)
!         
!            ! Intersects top boundary
!         IF ( elements(I)%MGSplitState == 1 ) THEN        
!            zBoundary = MGTop          
!         ELSE  ! Intersects the bottom boundary           
!            zBoundary = MGBottom          
!         END IF
!         
!         
!         CALL SplitElementOnZBoundary( zBoundary, I, numNodes, numElements, node1, node2, element, newNode, newElement, ErrStat, ErrMsg )
!         newNode%NodeType = 1 ! end node
!            ! Update the number of nodes and elements by one each
!         numNodes    = numNodes + 1
!         newNode%JointIndx = numNodes
!         numElements = numElements + 1
!         
!            ! Copy the altered nodes and elements into the master arrays
!         nodes(node1Indx)      = node1
!         nodes(node2Indx)      = node2
!         nodes(numNodes)       = newNode
!         elements(I)           = element
!         elements(numElements) = newElement
!         
!            ! If the original element spanned both marine growth boundaries, then we need to make an additional split
!         IF ( elements(I)%MGSplitState == 3 ) THEN 
!            
!            CALL SplitElementOnZBoundary( MGTop, numElements, numNodes, numElements, newNode, node2, newElement, newNode2, newElement2, ErrStat, ErrMsg )
!            newNode2%NodeType = 1 ! end node
!            newNode2%JointIndx = numNodes + 1
!               ! Copy the altered nodes and elements into the master arrays
!            nodes(numNodes)         = newNode
!            nodes(node2Indx)        = node2
!            nodes(numNodes+1)       = newNode2
!            elements(numElements)   = newElement
!            elements(numElements+1) = newElement2
!            
!               ! Update the number of nodes and elements by one each
!            numNodes    = numNodes + 1
!            
!            numElements = numElements + 1
!            
!         END IF
!      END IF 
!   END DO     
!END SUBROUTINE SplitElementsForMG

SUBROUTINE SplitElements(numNodes, nodes, numElements, elements, ErrStat, ErrMsg)   

   
   INTEGER,                  INTENT ( INOUT )  :: numNodes
   TYPE(Morison_NodeType),   INTENT ( INOUT )  :: nodes(:)   
   INTEGER,                  INTENT ( INOUT )  :: numElements
   TYPE(Morison_MemberType), INTENT ( INOUT )  :: elements(:)
   INTEGER,                  INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),             INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
   INTEGER                                     :: I, J, iCurrent, nSplits !, K 
   REAL(ReKi)                                  :: splits(5)
   INTEGER                                     :: node1Indx, node2Indx
   TYPE(Morison_NodeType)                      :: node1, node2, newNode !, newNode2
   TYPE(Morison_MemberType)                    :: element, newElement !, newElement2
   REAL(ReKi)                                  :: zBoundary
   INTEGER                                     :: origNumElements
   
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = "" 
   
   
   origNumElements = numElements
   
   DO I=1,origNumElements
     
      IF ( elements(I)%NumSplits > 0 ) THEN
         
         
            ! The splits are presorted from smallest Z to largest Z
         nSplits  = elements(I)%NumSplits
         splits   = elements(I)%Splits
         iCurrent = I
         
         DO J=1,nSplits
            
            node1Indx =  elements(iCurrent)%Node1Indx
            node1     = nodes(node1Indx)
            node2Indx =  elements(iCurrent)%Node2Indx
            node2     = nodes(node2Indx)
            element   = elements(iCurrent)
            
            CALL SplitElementOnZBoundary( 3, splits(J), iCurrent, numNodes, numElements, node1, node2, element, newNode, newElement, ErrStat, ErrMsg )
            
               ! Was this split due to the location of an elements free surface location crossing through the element?  
            IF ( element%MmbrFilledIDIndx  /= -1 ) THEN
               IF ( EqualRealNos(element%FillFSLoc, splits(J)) )  THEN
                  newElement%MmbrFilledIDIndx = -1
               END IF
            END IF
            
            newNode%NodeType = 1 ! end node
               ! Update the number of nodes and elements by one each
            numNodes    = numNodes + 1
            newNode%JointIndx = numNodes
            numElements = numElements + 1
         
               ! Copy the altered nodes and elements into the master arrays
            !nodes(node1Indx)      = node1
            nodes(node2Indx)      = node2
            nodes(numNodes)       = newNode
            elements(iCurrent)    = element
            elements(numElements) = newElement
            
            
            
               ! now make element = newElement by setting iCurrent to numElements
            iCurrent = numElements
            
            
         END DO           
    
      END IF 
   END DO     
END SUBROUTINE SplitElements

!====================================================================================================
!SUBROUTINE SplitElementsForWtr(MSL2SWL, Zseabed, numNodes, nodes, numElements, elements, ErrStat, ErrMsg)   
!
!   REAL(ReKi),               INTENT ( IN    )  :: MSL2SWL
!   REAL(ReKi),               INTENT ( IN    )  :: Zseabed
!   INTEGER,                  INTENT ( INOUT )  :: numNodes
!   TYPE(Morison_NodeType),   INTENT ( INOUT )  :: nodes(:)   
!   INTEGER,                  INTENT ( INOUT )  :: numElements
!   TYPE(Morison_MemberType), INTENT ( INOUT )  :: elements(:)
!   INTEGER,                  INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs  
!   CHARACTER(*),             INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None
!   
!   INTEGER                                     :: I, J, K
!   INTEGER                                     :: node1Indx, node2Indx
!   TYPE(Morison_NodeType)                      :: node1, node2, newNode, newNode2
!   TYPE(Morison_MemberType)                    :: element, newElement, newElement2
!   REAL(ReKi)                                  :: zBoundary
!   INTEGER                                     :: origNumElements
!   
!   origNumElements = numElements
!   
!   DO I=1,origNumElements
!     
!      IF ( elements(I)%WtrSplitState > 0 ) THEN
!         
!         node1Indx =  elements(I)%Node1Indx
!         node1 = nodes(node1Indx)
!         node2Indx =  elements(I)%Node2Indx
!         node2 = nodes(node2Indx)
!         element = elements(I)
!         
!            ! Intersects top boundary
!         IF ( elements(I)%WtrSplitState == 1 ) THEN        
!            zBoundary = MSL2SWL          
!         ELSE  ! Intersects the bottom boundary           
!            zBoundary = Zseabed          
!         END IF
!         
!         
!         CALL SplitElementOnZBoundary( zBoundary, I, numNodes, numElements, node1, node2, element, newNode, newElement, ErrStat, ErrMsg )
!         newNode%NodeType = 1 ! end node
!            ! Update the number of nodes and elements by one each
!         numNodes    = numNodes + 1
!         newNode%JointIndx = numNodes
!         numElements = numElements + 1
!         
!            ! Copy the altered nodes and elements into the master arrays
!         nodes(node1Indx)      = node1
!         nodes(node2Indx)      = node2
!         nodes(numNodes)       = newNode
!         elements(I)           = element
!         elements(numElements) = newElement
!         
!            ! If the original element spanned both marine growth boundaries, then we need to make an additional split
!         IF ( elements(I)%WtrSplitState == 3 ) THEN 
!            
!            CALL SplitElementOnZBoundary( MSL2SWL, numElements, numNodes, numElements, newNode, node2, newElement, newNode2, newElement2, ErrStat, ErrMsg )
!            newNode2%NodeType = 1 ! end node
!            newNode2%JointIndx = numNodes + 1
!               ! Copy the altered nodes and elements into the master arrays
!            nodes(numNodes)         = newNode
!            nodes(node2Indx)        = node2
!            nodes(numNodes+1)       = newNode2
!            elements(numElements)   = newElement
!            elements(numElements+1) = newElement2
!            
!               ! Update the number of nodes and elements by one each
!            numNodes    = numNodes + 1
!            
!            numElements = numElements + 1
!            
!         END IF
!      END IF 
!   END DO     
!END SUBROUTINE SplitElementsForWtr
!====================================================================================================
SUBROUTINE SubdivideMembers( numNodes, nodes, numElements, elements, ErrStat, ErrMsg )

   INTEGER,                  INTENT ( INOUT )  :: numNodes
   TYPE(Morison_NodeType),   INTENT ( INOUT )  :: nodes(:)   
   INTEGER,                  INTENT ( INOUT )  :: numElements
   TYPE(Morison_MemberType), INTENT ( INOUT )  :: elements(:)   
   INTEGER,                  INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),             INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
   TYPE(Morison_NodeType)                      :: node1, node2, newNode
   TYPE(Morison_MemberType)                    :: element, newElement
   INTEGER                                     :: numDiv
   REAL(ReKi)                                  :: divSize(3)
   INTEGER                                     :: I, J, K
   REAL(ReKi)                                  :: memLen
   INTEGER                                     :: origNumElements
   INTEGER                                     :: node1Indx, node2Indx, elementIndx, axis
   REAL(ReKi)                                  :: start, Loc
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = "" 
   
   
   origNumElements = numElements
   
   DO I=1,origNumElements
      
      CALL Morison_CopyMemberType( elements(I), element, MESH_NEWCOPY, ErrStat, ErrMsg )    !   element = elements(I); bjj: I'm replacing the "equals" here with the copy routine, 
                                                                                            !   though at this point there are no allocatable/pointer fields in this type so no harm done.
                                                                                            !   mostly for me to remember that when Inspector complains about uninitialized values, it's 
                                                                                            !   because not all *fields* have been initialized.
      node1Indx  = element%Node1Indx
      CALL Morison_CopyNodeType( nodes(node1Indx), node1, MESH_NEWCOPY, ErrStat, ErrMsg )   !   node1 = nodes(node1Indx); bjj: note that not all fields have been initialized, but maybe that's okay.

      node2Indx   = element%Node2Indx          ! We need this index for the last sub-element
      CALL Morison_CopyNodeType( nodes(node2Indx), node2, MESH_NEWCOPY, ErrStat, ErrMsg )   !   node2 = nodes(node2Indx); bjj: note that not all fields have been initialized, but maybe that's okay.
      
      elementIndx = I
            
      
      CALL GetDistance(node1%JointPos, node2%JointPos, memLen)
      
      
         ! If the requested division size is less then the member length, we will subdivide the member
         
      IF ( element%MDivSize < memLen ) THEN
         IF ( .NOT. ( EqualRealNos( node2%JointPos(3) , node1%JointPos(3) ) ) ) THEN
            axis  = 3
         ELSE IF ( .NOT. ( EqualRealNos( node2%JointPos(2) , node1%JointPos(2)  ) ) ) THEN       
            axis  = 2
         ELSE IF ( .NOT. ( EqualRealNos( node2%JointPos(1) , node1%JointPos(1)  ) ) ) THEN
            axis  = 1
         ELSE
            ! ERROR
         END IF
         
         start = node1%JointPos(axis)
         numDiv = CEILING( memLen / element%MDivSize  )
      
         DO K=1,3
            divSize(K) = (node2%JointPos(K) - node1%JointPos(K)) / numDiv
         END DO
      
         DO J=1,numDiv - 1
            
            loc = start + divSize(axis)*J
            CALL SplitElementOnZBoundary( axis, loc, elementIndx, numNodes, numElements, node1, node2, element, newNode, newElement, ErrStat, ErrMsg )
            newNode%NodeType = 2 ! interior node
            newNode%JointIndx = -1
               ! Copy updated node and element information to the nodes and elements arrays
            nodes(node1Indx)       = node1                  ! this is copying all the fields in the type; type contains no allocatable arrays or pointers
            nodes(node2Indx)       = node2                  ! this is copying all the fields in the type; type contains no allocatable arrays or pointers
            numNodes               = numNodes + 1
            numElements            = numElements + 1
            nodes(numNodes)        = newNode               ! this is copying all the fields in the type; type contains no allocatable arrays or pointers
            elements(elementIndx)  = element               ! this is copying all the fields in the type; type contains no allocatable arrays or pointers
            elements(numElements)  = newElement            ! this is copying all the fields in the type; type contains no allocatable arrays or pointers
            
            node1                  = newNode               ! this is copying all the fields in the type; type contains no allocatable arrays or pointers
            element                = newElement            ! this is copying all the fields in the type; type contains no allocatable arrays or pointers
            node1Indx              = numNodes
            elementIndx            = numElements 
            
            
            
            
            
            ! Create a new node
            !newNode = node2
            !
            !DO K=1,3
            !   newNode%JointPos = node1%JointPos(K) + divSize(K)*J
            !END DO
            !
            !numNodes = numNodes + 1
            !element%Node2Indx = numNodes
            !nodes(numNodes) = newNode
            !
            !   ! Create a new element
            !newElement = element
            !newElement%Node1Indx = numNodes
            !element = newElement
            !numElements = numElements + 1
         END DO
         
         !element%Node2Indx = node2Indx
         
      END IF  
      
   END DO
         
         
END SUBROUTINE SubdivideMembers         
      
SUBROUTINE CreateSuperMembers( )

 
      
         ! Determine if any of the joints flagged for overlap elimination (super member creation) satisfy the necessary requirements
      !IF ( InitInp%NJoints  > 0 ) THEN
      !
      !   DO I = 1,InitInp%NJoints
      !      
      !         ! Check #1 are there more than 2 members connected to the joint?
      !      IF ( InitInp%InpJoints(I)%JointOvrlp == 1 .AND. InitInp%InpJoints(I)%NConnections > 2) THEN
      !            ! Check #2 are there two members whose local z-axis are the same?
      !         CALL Get180Members(joint, InitInp%Joints, InitInp%Members, member1, member2)
      !         
      !         !zVect1
      !         !zVect2
      !         !dot(zVect1, zVect2)
      !         
      !      END IF
      !      
      !      
      !   END DO
      !   
      !
      !END IF
      
      
END SUBROUTINE CreateSuperMembers


!====================================================================================================
SUBROUTINE SetDepthBasedCoefs( z, NCoefDpth, CoefDpths, Cd, CdMG, Ca, CaMG, Cp, CpMG, AxCa, AxCaMG, AxCp, AxCpMG )
   
   REAL(ReKi), INTENT ( IN )              :: z
   INTEGER, INTENT (IN   ) :: NCoefDpth
   TYPE(Morison_CoefDpths), INTENT (IN   )  :: CoefDpths(:)
   REAL(ReKi), INTENT (  OUT)             :: Cd
   REAL(ReKi), INTENT (  OUT)             :: CdMG
   REAL(ReKi), INTENT (  OUT)             :: Ca
   REAL(ReKi), INTENT (  OUT)             :: CaMG
   REAL(ReKi), INTENT (  OUT)             :: Cp
   REAL(ReKi), INTENT (  OUT)             :: CpMG
   REAL(ReKi), INTENT (  OUT)             :: AxCa
   REAL(ReKi), INTENT (  OUT)             :: AxCaMG
   REAL(ReKi), INTENT (  OUT)             :: AxCp
   REAL(ReKi), INTENT (  OUT)             :: AxCpMG
   
   INTEGER                 :: I, indx1, indx2
   REAL(ReKi)              :: dd, s
   LOGICAL                 :: foundLess 
   
   
      ! Find the table entry(ies) which match the node's depth value
      ! The assumption here is that the depth table is stored from largest
      ! to smallest in depth
   
   foundLess = .FALSE.
   indx1     = 1
   indx2     = 1 
   
   DO I = 1, NCoefDpth
      IF ( CoefDpths(I)%Dpth <= z .AND. .NOT. foundLess ) THEN
         indx1 = I
         foundLess = .TRUE.
      END IF
      IF ( CoefDpths(I)%Dpth >= z ) THEN
         indx2 = I
      END IF
      
   END DO
   
      ! Linearly interpolate the coef values based on depth
   !CALL FindInterpFactor( z, CoefDpths(indx1)%Dpth, CoefDpths(indx2)%Dpth, s )
      
   dd = CoefDpths(indx1)%Dpth - CoefDpths(indx2)%Dpth
   IF ( EqualRealNos(dd, 0.0_ReKi) ) THEN
      s = 0
   ELSE
      s = ( CoefDpths(indx1)%Dpth - z ) / dd
   END IF
   
   Cd     = CoefDpths(indx1)%DpthCd*(1-s)     + CoefDpths(indx2)%DpthCd*s
   Ca     = CoefDpths(indx1)%DpthCa*(1-s)     + CoefDpths(indx2)%DpthCa*s
   Cp     = CoefDpths(indx1)%DpthCp*(1-s)     + CoefDpths(indx2)%DpthCp*s
   AxCa   = CoefDpths(indx1)%DpthCa*(1-s)     + CoefDpths(indx2)%DpthAxCa*s
   AxCp   = CoefDpths(indx1)%DpthCp*(1-s)     + CoefDpths(indx2)%DpthAxCp*s
   CdMG   = CoefDpths(indx1)%DpthCdMG*(1-s)   + CoefDpths(indx2)%DpthCdMG*s
   CaMG   = CoefDpths(indx1)%DpthCaMG*(1-s)   + CoefDpths(indx2)%DpthCaMG*s
   CpMG   = CoefDpths(indx1)%DpthCpMG*(1-s)   + CoefDpths(indx2)%DpthCpMG*s
   AxCaMG = CoefDpths(indx1)%DpthAxCaMG*(1-s) + CoefDpths(indx2)%DpthAxCaMG*s
   AxCpMG = CoefDpths(indx1)%DpthAxCpMG*(1-s) + CoefDpths(indx2)%DpthAxCpMG*s

END SUBROUTINE SetDepthBasedCoefs


!====================================================================================================
SUBROUTINE SetSplitNodeProperties( numNodes, nodes, numElements, elements, ErrStat, ErrMsg )   
!     This private subroutine generates the properties of nodes after the mesh has been split
!     the input data.  
!---------------------------------------------------------------------------------------------------- 

   INTEGER,                  INTENT ( IN    )  :: numNodes
   INTEGER,                  INTENT ( IN    )  :: numElements
   TYPE(Morison_MemberType), INTENT ( INOUT )  :: elements(:)
   TYPE(Morison_NodeType),   INTENT ( INOUT )  :: nodes(:)
   INTEGER,                  INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),             INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
   
   INTEGER                                     :: I
   TYPE(Morison_MemberType)                    :: element
   REAL(ReKi)                                  :: dR, dz
!   REAL(ReKi)                                  :: DirCos(3,3)
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   
   DO I=1,numNodes
      
      IF ( nodes(I)%NodeType /= 3 ) THEN
         
            ! End point or internal member node
            ! Super member nodes already have their properties set
            
            
         !element = elements(nodes(I)%ConnectionList(1))
         
            ! Calculate the element-level direction cosine matrix and attach it to the entry in the elements array
            
        ! CALL Morison_DirCosMtrx( nodes(element%Node1Indx)%JointPos, nodes(element%Node2Indx)%JointPos, elements(nodes(I)%ConnectionList(1))%R_LToG )
         
         element = elements(nodes(I)%ConnectionList(1))
         
         nodes(I)%R_LToG     = element%R_LToG
         
         nodes(I)%InpMbrIndx = element%InpMbrIndx
         IF ( .NOT. ( ( nodes(element%Node1Indx)%tMG > 0 ) .AND. ( nodes(element%Node2Indx)%tMG > 0 ) .AND. (.NOT. element%PropPot) ) )  THEN
            nodes(element%Node1Indx)%tMG       = 0.0
            nodes(element%Node2Indx)%tMG       = 0.0
            nodes(element%Node1Indx)%MGdensity = 0.0
            nodes(element%Node2Indx)%MGdensity = 0.0
         END IF
         IF ( element%Node1Indx == I ) THEN
            
            IF ( nodes(I)%tMG > 0 ) THEN
               nodes(I)%Cd   = element%CdMG1
               nodes(I)%Ca   = element%CaMG1
               nodes(I)%Cp   = element%CpMG1
               nodes(I)%AxCa   = element%AxCaMG1
               nodes(I)%AxCp   = element%AxCpMG1
            ELSE
               nodes(I)%Cd   = element%Cd1
               nodes(I)%Ca   = element%Ca1
               nodes(I)%Cp   = element%Cp1
               nodes(I)%AxCa   = element%AxCa1
               nodes(I)%AxCp   = element%AxCp1
            END IF
            
            nodes(I)%R    = element%R1
            nodes(I)%t    = element%t1
            nodes(I)%InpMbrDist = element%InpMbrDist1
         ELSE
            
            IF ( nodes(I)%tMG > 0 ) THEN
               nodes(I)%Cd   = element%CdMG2
               nodes(I)%Ca   = element%CaMG2
               nodes(I)%Cp   = element%CpMG2
               nodes(I)%AxCa   = element%AxCaMG2
               nodes(I)%AxCp   = element%AxCpMG2
            ELSE
               nodes(I)%Cd   = element%Cd2
               nodes(I)%Ca   = element%Ca2
               nodes(I)%Cp   = element%Cp2
               nodes(I)%AxCa   = element%AxCa2
               nodes(I)%AxCp   = element%AxCp2
            END IF
            
            nodes(I)%R    = element%R2
            nodes(I)%t    = element%t2
            nodes(I)%InpMbrDist = element%InpMbrDist2
         END IF
         
         CALL GetDistance( nodes(element%Node1Indx)%JointPos, nodes(element%Node2Indx)%JointPos, dz )
         dR = ( element%R2 + nodes(element%Node2Indx)%tMG ) - ( element%R1 + nodes(element%Node1Indx)%tMG )
         IF ( EqualRealNos(dR, 0.0_ReKi) ) dR = 0.0
         IF ( EqualRealNos(dz, 0.0_ReKi) ) THEN
            nodes(I)%dRdz = 0.0
         ELSE   
            nodes(I)%dRdz = dR / dz
         END IF
         
         nodes(I)%PropPot = element%PropPot
         
         IF ( element%MmbrFilledIDIndx /= -1 ) THEN
            nodes(I)%FillFlag  = .TRUE.
            nodes(I)%FillFSLoc    = element%FillFSLoc  ! This is relative to the MSL.
            nodes(I)%FillDensity  = element%FillDens
            
         ELSE
            nodes(I)%FillFSLoc    = 0.0  ! This is the MSL.
            nodes(I)%FillDensity  = 0.0
            elements(nodes(I)%ConnectionList(1))%FillDens  = 0.0
            elements(nodes(I)%ConnectionList(1))%FillFSLoc = 0.0
         END IF
            
         
     
         
      END IF
      
END DO

END SUBROUTINE SetSplitNodeProperties


!====================================================================================================
!SUBROUTINE SetMemberCoefs( SimplCd, SimplCdMG, SimplCa, SimplCaMG, CoefMembers, NCoefDpth, CoefDpths, element, node1, node2 )
SUBROUTINE SetElementCoefs( SimplCd, SimplCdMG, SimplCa, SimplCaMG, SimplCp, SimplCpMG, SimplAxCa, SimplAxCaMG, SimplAxCp, SimplAxCpMG,CoefMembers, NCoefDpth, CoefDpths, numNodes, nodes, numElements, elements )   
!     This private subroutine generates the Cd, Ca, Cp, CdMG, CaMG and CpMG coefs for the member based on
!     the input data.  
!---------------------------------------------------------------------------------------------------- 

   REAL(ReKi),                INTENT( IN    )  :: SimplCd 
   REAL(ReKi),                INTENT( IN    )  :: SimplCdMG
   REAL(ReKi),                INTENT( IN    )  :: SimplCa
   REAL(ReKi),                INTENT( IN    )  :: SimplCaMG 
   REAL(ReKi),                INTENT( IN    )  :: SimplCp
   REAL(ReKi),                INTENT( IN    )  :: SimplCpMG 
   REAL(ReKi),                INTENT( IN    )  :: SimplAxCa
   REAL(ReKi),                INTENT( IN    )  :: SimplAxCaMG 
   REAL(ReKi),                INTENT( IN    )  :: SimplAxCp
   REAL(ReKi),                INTENT( IN    )  :: SimplAxCpMG 
   TYPE(Morison_CoefMembers), INTENT( IN    )  :: CoefMembers(:)
   INTEGER,                   INTENT( IN    )  :: NCoefDpth
   TYPE(Morison_CoefDpths),   INTENT( IN    )  :: CoefDpths(:)
   INTEGER,                   INTENT( IN    )  :: numNodes
   INTEGER,                   INTENT( IN    )  :: numElements
   TYPE(Morison_MemberType),  INTENT( INOUT )  :: elements(:)
   TYPE(Morison_NodeType),    INTENT( IN    )  :: nodes(:)
   
   TYPE(Morison_NodeType)                      :: node1, node2
   
   INTEGER                                     :: MCoefMod
   INTEGER                                     :: I, J
   REAL(ReKi)                                  :: Cd, CdMG, Ca, CaMG, Cp, CpMG, AxCa, AxCp, AxCaMG, AxCpMG
   DO I=1,numElements
      
      
      MCoefMod = elements(I)%MCoefMod
      node1    = nodes(elements(I)%Node1Indx)
      node2    = nodes(elements(I)%Node2Indx)
      
      SELECT CASE ( MCoefMod )
      
      CASE (1)
      
         elements(I)%Cd1   = SimplCd
         elements(I)%Cd2   = SimplCd
         elements(I)%Ca1   = SimplCa
         elements(I)%Ca2   = SimplCa
         elements(I)%Cp1   = SimplCp
         elements(I)%Cp2   = SimplCp
         elements(I)%AxCa1   = SimplAxCa
         elements(I)%AxCa2   = SimplAxCa
         elements(I)%AxCp1   = SimplAxCp
         elements(I)%AxCp2   = SimplAxCp
         elements(I)%CdMG1 = SimplCdMG
         elements(I)%CdMG2 = SimplCdMG
         elements(I)%CaMG1 = SimplCaMG
         elements(I)%CaMG2 = SimplCaMG
         elements(I)%CpMG1 = SimplCpMG
         elements(I)%CpMG2 = SimplCpMG
         elements(I)%AxCaMG1 = SimplAxCaMG
         elements(I)%AxCaMG2 = SimplAxCaMG
         elements(I)%AxCpMG1 = SimplAxCpMG
         elements(I)%AxCpMG2 = SimplAxCpMG
      
      CASE (2)
       
         CALL SetDepthBasedCoefs( node1%JointPos(3), NCoefDpth, CoefDpths, Cd, CdMG, Ca, CaMG, Cp, CpMG, AxCa, AxCaMG, AxCp, AxCpMG )
         elements(I)%Cd1     = Cd
         elements(I)%Ca1     = Ca
         elements(I)%Cp1     = Cp
         elements(I)%AxCa1     = AxCa
         elements(I)%AxCp1     = AxCp
         elements(I)%CdMG1   = CdMG
         elements(I)%CaMG1   = CaMG  
         elements(I)%CpMG1   = CpMG
         elements(I)%AxCaMG1   = AxCaMG  
         elements(I)%AxCpMG1   = AxCpMG
         
         CALL SetDepthBasedCoefs( node2%JointPos(3), NCoefDpth, CoefDpths, Cd, CdMG, Ca, CaMG, Cp, CpMG, AxCa, AxCaMG, AxCp, AxCpMG )
         elements(I)%Cd2     = Cd
         elements(I)%Ca2     = Ca
         elements(I)%Cp2     = Cp
         elements(I)%AxCa2     = Ca
         elements(I)%AxCp2     = Cp
         elements(I)%CdMG2   = CdMG
         elements(I)%CaMG2   = CaMG
         elements(I)%CpMG2   = CpMG
         elements(I)%AxCaMG2   = AxCaMG
         elements(I)%AxCpMG2   = AxCpMG
         
      CASE (3)
      
         J          = elements(I)%MmbrCoefIDIndx
         elements(I)%Cd1   = CoefMembers(J)%MemberCd1
         elements(I)%Cd2   = CoefMembers(J)%MemberCd2
         elements(I)%Ca1   = CoefMembers(J)%MemberCa1
         elements(I)%Ca2   = CoefMembers(J)%MemberCa2
         elements(I)%Cp1   = CoefMembers(J)%MemberCp1
         elements(I)%Cp2   = CoefMembers(J)%MemberCp2
         elements(I)%AxCa1   = CoefMembers(J)%MemberAxCa1
         elements(I)%AxCa2   = CoefMembers(J)%MemberAxCa2
         elements(I)%AxCp1   = CoefMembers(J)%MemberAxCp1
         elements(I)%AxCp2   = CoefMembers(J)%MemberAxCp2
         elements(I)%CdMG1 = CoefMembers(J)%MemberCdMG1
         elements(I)%CdMG2 = CoefMembers(J)%MemberCdMG2
         elements(I)%CaMG1 = CoefMembers(J)%MemberCaMG1
         elements(I)%CaMG2 = CoefMembers(J)%MemberCaMG2
         elements(I)%CpMG1 = CoefMembers(J)%MemberCpMG1
         elements(I)%CpMG2 = CoefMembers(J)%MemberCpMG2
         elements(I)%AxCaMG1 = CoefMembers(J)%MemberAxCaMG1
         elements(I)%AxCaMG2 = CoefMembers(J)%MemberAxCaMG2
         elements(I)%AxCpMG1 = CoefMembers(J)%MemberAxCpMG1
         elements(I)%AxCpMG2 = CoefMembers(J)%MemberAxCpMG2
         
      END SELECT
   
      
   END DO
   
END SUBROUTINE SetElementCoefs


SUBROUTINE SetAxialCoefs( NJoints, NAxCoefs, AxialCoefs, numNodes, nodes, numElements, elements )

   INTEGER,                    INTENT( IN    )  :: NJoints
   INTEGER,                    INTENT( IN    )  :: NAxCoefs
   TYPE(Morison_AxialCoefType),INTENT( IN    )  :: AxialCoefs(:)
   INTEGER,                    INTENT( IN    )  :: numNodes
   INTEGER,                    INTENT( IN    )  :: numElements
   TYPE(Morison_MemberType),   INTENT( INOUT )  :: elements(:)
   TYPE(Morison_NodeType),     INTENT( INOUT )  :: nodes(:)
   
 !  TYPE(Morison_NodeType)                       :: node1, node2 
   
   INTEGER                                     :: I !, J
   
   DO I=1,numNodes
      
      IF ( nodes(I)%JointAxIDIndx > 0 .AND. nodes(I)%JointIndx > 0 .AND. nodes(I)%JointIndx <= NJoints) THEN
         nodes(I)%JAxCd = AxialCoefs(nodes(I)%JointAxIDIndx)%AxCd
         nodes(I)%JAxCa = AxialCoefs(nodes(I)%JointAxIDIndx)%AxCa
         nodes(I)%JAxCp = AxialCoefs(nodes(I)%JointAxIDIndx)%AxCp
      ELSE ! These are end nodes that were generated by the software, and hence do not have lumped axial loads, or they are interior nodes.
         nodes(I)%JAxCd = 0.0
         nodes(I)%JAxCa = 0.0
         nodes(I)%JAxCp = 0.0
      END IF
      
      !node1    = nodes(elements(I)%Node1Indx)
      !node2    = nodes(elements(I)%Node2Indx)
      
   END DO
   
END SUBROUTINE SetAxialCoefs


SUBROUTINE SetNodeMG( numMGDepths, MGDepths, numNodes, nodes )

   INTEGER,                      INTENT( IN    )  :: numMGDepths
   TYPE(Morison_MGDepthsType),   INTENT( IN    )  :: MGDepths(:)
   INTEGER,                      INTENT( IN    )  :: numNodes
   TYPE(Morison_NodeType),       INTENT( INOUT )  :: nodes(:)

   INTEGER                                     :: I, J
   REAL(ReKi)              :: z
   INTEGER                 :: indx1, indx2
   REAL(ReKi)              :: dd, s
   LOGICAL                 :: foundLess = .FALSE.
   
   DO I=1,numNodes
      
         !Find the table entry(ies) which match the node's depth value
      ! The assumption here is that the depth table is stored from largest
      ! to smallest in depth
      z = nodes(I)%JointPos(3)
      foundLess = .FALSE.
      indx1 = 0
      indx2 = 0
      DO J = 1, numMGDepths
         IF ( MGDepths(J)%MGDpth <= z .AND. .NOT. foundLess ) THEN
            indx1 = J
            
            foundLess = .TRUE.
         END IF
         IF ( MGDepths(J)%MGDpth >= z ) THEN
            indx2 = J
         END IF
      
      END DO
      IF ( indx2 == 0 .OR. .NOT. foundLess ) THEN
         !Not at a marine growth depth
         nodes(I)%tMG       = 0.0
         nodes(I)%MGdensity = 0.0
      ELSE
         ! Linearly interpolate the coef values based on depth
         !CALL FindInterpFactor( z, CoefDpths(indx1)%Dpth, CoefDpths(indx2)%Dpth, s )
      
         dd = MGDepths(indx1)%MGDpth - MGDepths(indx2)%MGDpth
         IF ( EqualRealNos(dd, 0.0_ReKi) ) THEN
            s = 0.0_ReKi
         ELSE
            s = ( MGDepths(indx1)%MGDpth - z ) / dd
         END IF
         nodes(I)%tMG       = MGDepths(indx1)%MGThck*(1-s) + MGDepths(indx2)%MGThck*s
         nodes(I)%MGdensity = MGDepths(indx1)%MGDens*(1-s) + MGDepths(indx2)%MGDens*s
      END IF
      
   END DO
   

END SUBROUTINE SetNodeMG



SUBROUTINE SetElementFillProps( numFillGroups, filledGroups, numElements, elements )

   INTEGER,                       INTENT( IN    )  :: numFillGroups
   TYPE(Morison_FilledGroupType), INTENT( IN    )  :: filledGroups(:)
   INTEGER,                       INTENT( IN    )  :: numElements
   TYPE(Morison_MemberType),      INTENT( INOUT )  :: elements(:)  
   
   INTEGER                                         :: I !, J
   
   DO I=1,numElements
      
      
      IF ( elements(I)%MmbrFilledIDIndx > 0 ) THEN
         
         elements(I)%FillDens     =   filledGroups(elements(I)%MmbrFilledIDIndx)%FillDens
         elements(I)%FillFSLoc    =   filledGroups(elements(I)%MmbrFilledIDIndx)%FillFSLoc
      ELSE
         elements(I)%FillDens     =   0.0
         elements(I)%FillFSLoc    =   0.0
      END IF
      
     
      
   END DO
   
   
END SUBROUTINE SetElementFillProps

!SUBROUTINE CreateLumpedMarkers( numNodes, nodes, numElements, elements, numLumpedMarkers, lumpedMarkers, ErrStat, ErrMsg )
!
!   INTEGER,                   INTENT( IN    )  :: numNodes
!   INTEGER,                   INTENT( IN    )  :: numElements
!   TYPE(Morison_MemberType),  INTENT( IN    )  :: elements(:)  
!   TYPE(Morison_NodeType),    INTENT( IN    )  :: nodes(:)
!   INTEGER,                   INTENT(   OUT )  :: numLumpedMarkers
!   TYPE(Morison_NodeType), ALLOCATABLE,    INTENT(   OUT )  :: lumpedMarkers(:)
!   INTEGER,                  INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs  
!   CHARACTER(*),             INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None
!   
!   INTEGER                                     :: I, J, count
!   TYPE(Morison_MemberType)                    :: element
!   
!   numLumpedMarkers = 0
!   
!      ! Count how many distributed markers we need to create by looping over the nodes
!   DO I=1,numNodes
!      IF ( nodes(I)%NodeType == 1 .AND. nodes(I)%JointOvrlp == 0 ) THEN
!            ! end of a member that was not a part of super member creation
!         numLumpedMarkers = numLumpedMarkers + 1
!      END IF
!   END DO
!   
!      ! Allocate the array for the distributed markers
!   ALLOCATE ( lumpedMarkers(numLumpedMarkers), STAT = ErrStat )
!   IF ( ErrStat /= ErrID_None ) THEN
!      ErrMsg  = ' Error allocating space for the lumped load markers array.'
!      ErrStat = ErrID_Fatal
!      RETURN
!   END IF  
!   count = 1   
!   DO I=1,numNodes
!      
!      IF ( nodes(I)%NodeType == 1 .AND. nodes(I)%JointOvrlp == 0) THEN
!         
!         element = elements(nodes(I)%ConnectionList(1))
!         
!         IF ( element%Node1Indx == I ) THEN
!            lumpedMarkers(count)%Cd   = element%Cd1
!            lumpedMarkers(count)%CdMG = element%CdMG1
!            lumpedMarkers(count)%Ca   = element%Ca1
!            lumpedMarkers(count)%CaMG = element%CaMG1
!            lumpedMarkers(count)%R    = element%R1
!            lumpedMarkers(count)%t    = element%t1
!         ELSE
!            lumpedMarkers(count)%Cd   = element%Cd2
!            lumpedMarkers(count)%CdMG = element%CdMG2
!            lumpedMarkers(count)%Ca   = element%Ca2
!            lumpedMarkers(count)%CaMG = element%CaMG2
!            lumpedMarkers(count)%R    = element%R2
!            lumpedMarkers(count)%t    = element%t2
!         END IF
!         
!         lumpedMarkers(count)%PropPot = element%PropPot
!         lumpedMarkers(count)%tMG       = nodes(I)%tMG
!         lumpedMarkers(count)%MGdensity = nodes(I)%MGdensity
!         
!         
!            ! Compute all initialization forces now so we have access to the element information
!            
!         !IF ( element%PropPot == .FALSE. ) THEN
!         !   
!         !      ! Member is not modeled with WAMIT
!         !   CALL LumpedBuoyancy( )             
!         !   CALL LumpedMGLoads( )       
!         !   CALL LumpedDynPressure( )
!         !   CALL LumpedAddedMass( )
!         !   CALL LumpedAddedMassMG( )
!         !   CALL LumpedAddedMassFlood( )  ! Do we actually compute this??? TODO
!         !   
!         !END IF
!         !   
!         !   ! These are the only two loads we compute at initialization if the member is modeled with WAMIT
!         !CALL LumpedDragConst( ) 
!         !CALL LumpedFloodedBuoyancy( )  
!         
!         
!         count = count + 1        
!      
!      END IF
!      
!   END DO
!   
!END SUBROUTINE CreateLumpedMarkers


SUBROUTINE SplitMeshNodes( numNodes, nodes, numElements, elements, numSplitNodes, splitNodes, ErrStat, ErrMsg )

   INTEGER,                   INTENT( IN    )  :: numNodes
   INTEGER,                   INTENT( IN    )  :: numElements
   TYPE(Morison_MemberType),  INTENT( INOUT )  :: elements(:)
   TYPE(Morison_NodeType),    INTENT( IN    )  :: nodes(:)
   INTEGER,                   INTENT(   OUT )  :: numSplitNodes
   TYPE(Morison_NodeType), ALLOCATABLE,    INTENT(   OUT )  :: splitNodes(:)
   INTEGER,                  INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),             INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
   INTEGER                                     :: I, J, splitNodeIndx
!   TYPE(Morison_MemberType)                    :: element
   TYPE(Morison_NodeType)                      :: node1, node2, newNode
   
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = "" 
   
   numSplitNodes = 0
   
      ! Count how many distributed markers we need to create by looping over the nodes
   DO I=1,numNodes
      IF ( nodes(I)%NodeType == 1 ) THEN
            ! Nodes at the end of members get one node for each connecting member
         numSplitNodes = numSplitNodes + nodes(I)%NConnections     
      ELSE      
          ! Internal nodes and super member nodes only get one node
         numSplitNodes = numSplitNodes + 1
      END IF
   END DO
   
      ! Allocate the array for the distributed markers
   ALLOCATE ( splitNodes(numSplitNodes), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for split nodes array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF  
   
   splitNodes(1:numNodes)    = nodes(1:numNodes)
   
   IF ( numSplitNodes > numNodes ) THEN
      
      splitNodeIndx = numNodes + 1 
   
      DO I=1,numElements
         ! Loop over elements in the processed mesh and create additional nodes/markers at the end of elements if that node connects to other elements
         node1 = splitnodes(elements(I)%Node1Indx)
         node2 = splitnodes(elements(I)%Node2Indx)
      
         IF (node1%NodeType == 1 ) THEN ! end node
            IF ( node1%NConnections > 1 ) THEN
                  !create new node by copying the old one
               newNode = node1
               newNode%NConnections = 1
               splitnodes(splitNodeIndx) = newNode
               splitnodes(elements(I)%Node1Indx)%NConnections = node1%NConnections - 1
               !set the new node as the first node of this element
               elements(I)%Node1Indx = splitNodeIndx
               splitNodeIndx = splitNodeIndx + 1
               !NOTE: the node connection list entries are now bogus!!!!
            END IF
      
         END IF
      
         IF (node2%NodeType == 1 ) THEN ! end node
            IF ( node2%NConnections > 1 ) THEN
                  !create new node by copying the old one
               newNode = node2
               newNode%NConnections = 1
               splitnodes(splitNodeIndx) = newNode
               splitnodes(elements(I)%Node2Indx)%NConnections = node2%NConnections - 1
               !set the new node as the first node of this element
               elements(I)%Node2Indx = splitNodeIndx
               splitNodeIndx = splitNodeIndx + 1
               !NOTE: the node connection list entries are now bogus!!!!
            END IF
      
         END IF
      
      END DO
      
      ! Fix connections
      DO J = 1,numSplitNodes
         splitnodes(J)%NConnections = 0
      END DO
      
      DO I = 1,numElements
        
            
         DO J = 1,numSplitNodes
            IF ( elements(I)%Node1Indx == J ) THEN
               splitnodes(J)%NConnections = splitnodes(J)%NConnections + 1
               splitnodes(J)%ConnectionList(splitnodes(J)%NConnections) = I
            END IF 
            IF ( elements(I)%Node2Indx == J ) THEN
               splitnodes(J)%NConnections = splitnodes(J)%NConnections + 1
               splitnodes(J)%ConnectionList(splitnodes(J)%NConnections) = I
            END IF 
         END DO
      END DO
      
  END IF 
   
END SUBROUTINE SplitMeshNodes



SUBROUTINE GenerateLumpedLoads( nodeIndx, sgn, node, gravity, MSL2SWL, densWater, NStepWave, WaveDynP, dragConst, F_DP, F_B,  ErrStat, ErrMsg )

   INTEGER,                 INTENT( IN    )     ::  nodeIndx
   REAL(ReKi),              INTENT( IN    )     ::  sgn
   TYPE(Morison_NodeType),  INTENT( IN    )     ::  node
   REAL(ReKi),              INTENT( IN    )     ::  gravity
   REAL(ReKi),              INTENT( IN    )     ::  MSL2SWL
   REAL(ReKi),              INTENT( IN    )     ::  densWater
   INTEGER,                 INTENT( IN    )     ::  NStepWave
   REAL(SiKi),              INTENT( IN    )     ::  WaveDynP(:,:) ! TODO:  Verify it is ok to use (:,:) for the  zero-based  first array index GJH 2/5/14
   REAL(ReKi),ALLOCATABLE,  INTENT(   OUT )     ::  F_DP(:,:)
   REAL(ReKi),              INTENT(   OUT )     ::  F_B(6)
   REAL(ReKi),              INTENT(   OUT )     ::  dragConst
   INTEGER,                 INTENT(   OUT )     ::  ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),            INTENT(   OUT )     ::  ErrMsg               ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                                   ::  k(3)
 
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = "" 
   
   
   IF (.NOT. node%PropPot ) THEN
   
      k =  sgn * node%R_LToG(:,3)
      
      CALL LumpDynPressure( nodeIndx, node%JAxCp, k, node%R, node%tMG, NStepWave, WaveDynP, F_DP, ErrStat, ErrMsg)
      
         ! For buoyancy calculations we need to adjust the Z-location based on MSL2SWL. If MSL2SWL > 0 then SWL above MSL, and so we need to place the Z value at a deeper position.  
         !   SWL is at Z=0 for buoyancy calcs, but geometry was specified relative to MSL (MSL2SWL = 0) 
      CALL LumpBuoyancy( sgn, densWater, node%R, node%tMG, node%JointPos(3) - MSL2SWL, node%R_LToG, gravity, F_B  ) 
                    
       
      ! This one is tricky because we need to calculate a signed volume which is the signed sum of all connecting elements and then split the 
      ! result across all the connecting nodes.
      !CALL LumpAddedMass() 
   
   ELSE
      
      ALLOCATE ( F_DP(0:NStepWave, 6), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating distributed dynamic pressure loads array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF  
      F_DP = 0.0
      F_B  = 0.0
   END IF
   
   
   CALL LumpDragConst( densWater, node%Cd, node%R, node%tMG, dragConst ) 
   
   
     

END SUBROUTINE GenerateLumpedLoads



SUBROUTINE CreateLumpedMesh( densWater, gravity, MSL2SWL, wtrDpth, NStepWave, WaveDynP, WaveAcc, numNodes, nodes, numElements, elements, &
                                  numLumpedMarkers, lumpedMeshIn, lumpedMeshOut, lumpedToNodeIndx, L_An,        &
                                  L_F_B, L_F_I, L_F_BF, L_AM_M, L_dragConst, &
                                  ErrStat, ErrMsg )

   REAL(ReKi),                             INTENT( IN    )  ::  densWater
   REAL(ReKi),                             INTENT( IN    )  ::  gravity
   REAL(ReKi),                             INTENT( IN    )  ::  MSL2SWL
   REAL(ReKi),                             INTENT( IN    )  ::  wtrDpth
   INTEGER,                                INTENT( IN    )  ::  NStepWave
   REAL(SiKi),                             INTENT( IN    )  ::  WaveDynP(0:,:)
   REAL(SiKi),                             INTENT( IN    )  ::  WaveAcc(0:,:,:)
   INTEGER,                                INTENT( IN    )  ::  numNodes
   INTEGER,                                INTENT( IN    )  ::  numElements
   TYPE(Morison_MemberType),               INTENT( IN    )  ::  elements(:)
   TYPE(Morison_NodeType),                 INTENT( INOUT )  ::  nodes(:)
   INTEGER,                                INTENT(   OUT )  ::  numLumpedMarkers
   !TYPE(Morison_NodeType), ALLOCATABLE,    INTENT(   OUT )  ::  lumpedMarkers(:)
   TYPE(MeshType),                         INTENT(   OUT )  ::  lumpedMeshIn
   TYPE(MeshType),                         INTENT(   OUT )  ::  lumpedMeshOut 
   INTEGER, ALLOCATABLE,                   INTENT(   OUT )  ::  lumpedToNodeIndx(:)
   REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  L_An(:,:)                         ! The signed/summed end cap Area x k of all connected members at a common joint
   REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  L_F_B(:,:)                      ! Buoyancy force associated with the member
   REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  L_F_I(:,:,:)                     ! Inertial force.  TODO:  Eventually the dynamic pressure will be included in this force! GJH 4/15/14
   !REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  L_F_DP(:,:,:)                     ! Dynamic pressure force
   REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  L_F_BF(:,:)                     ! Flooded buoyancy force
   REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  L_AM_M(:,:,:)                   ! Added mass of member
   REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  L_dragConst(:)                   ! 
   INTEGER,                                INTENT(   OUT )  ::  ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                           INTENT(   OUT )  ::  ErrMsg               ! Error message if ErrStat /= ErrID_None
   
             
   INTEGER                    ::  I, J, M, count
   TYPE(Morison_MemberType)   ::  element
   TYPE(Morison_NodeType)     ::  node, node1, node2
   REAL(ReKi)                 ::  L, sgn
!   REAL(ReKi)                 ::  k(3)
   REAL(ReKi)                 ::  z0
   REAL(ReKi),ALLOCATABLE     ::  F_DP(:,:)
   REAL(ReKi)                 ::  F_B(6)
   REAL(ReKi)                 ::  F_BF(6)
   REAL(ReKi)                 ::  Vmat(3,1), F_I(6), AM_M(6,6) !AM(6,6), 
   REAL(ReKi)                 ::  dragConst
   
   
   INTEGER, ALLOCATABLE       :: nodeToLumpedIndx(:)
   INTEGER, ALLOCATABLE       :: commonNodeLst(:)
   LOGICAL, ALLOCATABLE       :: usedJointList(:)
   INTEGER                    :: nCommon
!   REAL(ReKi)                 :: CA
   REAL(ReKi)                 :: AMfactor
   REAL(ReKi)                 :: An(3), Vn(3), af(3)
!   REAL(ReKi)                 :: AM11, AM22, AM33
   REAL(ReKi)                 :: f1, VnDotAf, Vmag !f2, 
   
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = "" 
   
   
   numLumpedMarkers = 0
   z0                = -(wtrDpth) ! The total sea depth is the still water depth of the seabed 
   
   
   
      ! Count how many lumped markers we need to create by looping over the nodes
      
   DO I=1,numNodes
      
      IF ( (nodes(I)%NodeType == 3) .OR. ( nodes(I)%NodeType == 1 .AND. nodes(I)%JointOvrlp == 0 )  ) THEN
      
            numLumpedMarkers = numLumpedMarkers + 1
      
      END IF
      
   END DO
   
   
      ! Create the input and output meshes associated with lumped loads
      
   CALL MeshCreate( BlankMesh      = lumpedMeshIn           &
                     ,IOS          = COMPONENT_INPUT        &
                     ,Nnodes       = numLumpedMarkers       &
                     ,ErrStat      = ErrStat                &
                     ,ErrMess      = ErrMsg                 &
                     ,TranslationDisp = .TRUE.             &
                     ,Orientation     = .TRUE.                 &
                     ,TranslationVel  = .TRUE.              &
                     ,RotationVel     = .TRUE.              &
                     ,TranslationAcc  = .TRUE.              &
                     ,RotationAcc     = .TRUE.     )

    
   
   
   !   ! Allocate the array for the lumped markers
   !   
   !ALLOCATE ( lumpedMarkers(numLumpedMarkers), STAT = ErrStat )
   !IF ( ErrStat /= ErrID_None ) THEN
   !   ErrMsg  = ' Error allocating space for the lumped load markers array.'
   !   ErrStat = ErrID_Fatal
   !   RETURN
   !END IF  
   
   ALLOCATE ( commonNodeLst(10), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the commonNodeLst array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF 
   commonNodeLst = -1
   
   ALLOCATE ( usedJointList(numNodes), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the UsedJointList array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF  
   
   usedJointList = .FALSE.
   
   ALLOCATE ( lumpedToNodeIndx(numLumpedMarkers), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the lumped index array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF  
   
   ALLOCATE ( nodeToLumpedIndx(numNodes), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the lumped index array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF  
   
   
      
   ALLOCATE ( L_F_B( 6, numLumpedMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the lumped buoyancy forces/moments array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   L_F_B = 0.0
   
   ALLOCATE ( L_An( 3, numLumpedMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the L_An array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   L_An = 0.0
   
   ! This is 
   ALLOCATE ( L_F_I( 0:NStepWave, 6, numLumpedMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the lumped inertial forces/moments array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
    L_F_I = 0.0
    
   !ALLOCATE ( L_F_DP( 0:NStepWave, 6, numLumpedMarkers ), STAT = ErrStat )
   !IF ( ErrStat /= ErrID_None ) THEN
   !   ErrMsg  = ' Error allocating space for the lumped dynamic pressure forces/moments array.'
   !   ErrStat = ErrID_Fatal
   !   RETURN
   !END IF
   ! L_F_DP = 0.0
   
   ALLOCATE ( L_F_BF( 6, numLumpedMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the lumped buoyancy due to flooding forces/moments array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   L_F_BF = 0.0
   
   ALLOCATE ( L_AM_M( 6, 6, numLumpedMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the lumped member added mass.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   L_AM_M = 0.0
      
   ALLOCATE ( L_dragConst( numLumpedMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the lumped drag constants.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   L_dragConst = 0.0
   
   ! Loop over nodes to create all loads on the resulting markers except for the buoyancy loads
   ! For the buoyancy loads, loop over the elements and then apply one half of the resulting value
   ! to each of the interior element nodes but the full value to an end node.  This means that an internal member node will receive 1/2 of its
   ! load from element A and 1/2 from element B.  If it is the end of a member it will simply receive
   ! the element A load.
   
   count = 1 
   
   DO I=1,numNodes
      
          ! exclude internal member nodes and end nodes which were connected to a joint made into a super member
          
      IF ( (nodes(I)%NodeType == 3) .OR. ( nodes(I)%NodeType == 1 .AND. nodes(I)%JointOvrlp == 0 )  ) THEN
         
            lumpedToNodeIndx(count) = I
            nodeToLumpedIndx(I) = count
               
               ! If this is a super member node, then generate the lumped loads now, otherwise save it for the loop over elements
               
            IF ( nodes(I)%NodeType == 3 ) THEN
               
            END IF
            

            
            
               ! Create the node on the mesh
            
            CALL MeshPositionNode (lumpedMeshIn          &
                              , count                    &
                              , nodes(I)%JointPos        &  ! this info comes from FAST
                              , ErrStat                  &
                              , ErrMsg                   &
                              ) !, transpose(nodes(I)%R_LToG)          )
            IF ( ErrStat /= 0 ) THEN
               RETURN
            END IF 
         
               ! Create the mesh element
         
            CALL MeshConstructElement (  lumpedMeshIn   &
                                  , ELEMENT_POINT      &                                  
                                  , ErrStat            &
                                  , ErrMsg  &
                                  , count                  &
                                              )
            count = count + 1    
         
         
      END IF
            
   END DO
   

   
   
   ! CA is the added mass coefficient for three dimensional bodies in infinite fluid (far from boundaries) The default value is 2/Pi
     
   
   AMfactor = 2.0 * densWater * Pi / 3.0
 
       ! Loop over nodes again in order to create lumped axial drag. 
       
   usedJointList = .FALSE.   
   commonNodeLst = -1
   
   DO I=1,numNodes
      
         
      
            ! Determine bounds checking based on what load we are calculating, 
            ! This is for L_An
         IF ( nodes(I)%JointPos(3) >= z0 ) THEN
         
               ! exclude internal member nodes and end nodes which were connected to a joint made into a super member
      
         
               
                  ! If this is a super member node, then generate the lumped loads now, otherwise save it for the loop over elements
               
               IF ( nodes(I)%NodeType == 3 ) THEN
               
               ELSE
            
                  IF ( nodes(I)%JointIndx /= -1 ) THEN  ! TODO: MAYBE THIS SHOULD CHECK JointOvrlp value instead!!
               
                     ! Have we already set the added mass for this node?
                  IF ( .NOT. usedJointList(nodes(I)%JointIndx) ) THEN
                  
                     nCommon   = 0
                     An        = 0.0
                     Vn        = 0.0
                     
                     DO J=1,numNodes
                     
                           ! must match joint index but also cannot be modeled using WAMIT
                        IF  ( nodes(I)%JointIndx == nodes(J)%JointIndx ) THEN
                           
                           nCommon = nCommon + 1
                           commonNodeLst(nCommon) = J
                        
                              ! Compute the signed area*outward facing normal of this member
                           sgn = 1.0
                           
                           element = elements(nodes(J)%ConnectionList(1))
                           
                           IF ( element%Node1Indx == J ) THEN
                              sgn = -1.0                                ! Local coord sys points into element at starting node, so flip sign of local z vector
                           ELSE IF ( element%Node2Indx == J ) THEN
                              sgn = 1.0                                 ! Local coord sys points out of element at ending node, so leave sign of local z vector
                           ELSE
                              ErrMsg  = 'Internal Error in CreateLumpedMesh: could not find element node index match.'
                              ErrStat = ErrID_FATAL
                              RETURN
                           END IF
                              ! Compute the signed volume of this member
                           f1 = (nodes(J)%R+nodes(J)%tMG)*(nodes(J)%R+nodes(J)%tMG)*(nodes(J)%R+nodes(J)%tMG) 
                           Vn = Vn + sgn*f1*nodes(J)%R_LToG(:,3)
                           An = An + sgn*nodes(J)%R_LtoG(:,3)*Pi*(nodes(J)%R+nodes(J)%tMG)**2
                                         
                        END IF
                     
                     END DO
                  
                     nodes(I)%NConnectPreSplit = nCommon
                     
                        ! Divide the directed area equally across all connected markers 
                     Vmag = sqrt(Dot_Product(Vn,Vn))
                     
                     
                     
                     AM_M = 0.0
                     IF ( (Vmag > 0.0) .AND. (.NOT. nodes(I)%PropPot) ) THEN
                        Vmat = RESHAPE(Vn,(/3,1/))
                        AM_M(1:3,1:3) = (nodes(I)%JAxCa*AMfactor/(REAL( nCommon, ReKi)*Vmag) )*MatMul(Vmat,TRANSPOSE(Vmat))
                     END IF
                     
                     DO J=1,nCommon
                     
                        IF ( nodes(I)%JointPos(3) >= z0 ) THEN
                        
                           L_An  (:,  nodeToLumpedIndx(commonNodeLst(J)))    =  An / nCommon   
                           L_AM_M(:,:,nodeToLumpedIndx(commonNodeLst(J)))    =  AM_M
                           
                           DO M=0,NStepWave
                              ! The WaveAcc array has indices of (timeIndx, nodeIndx, vectorIndx), the nodeIndx needs to correspond to the total list of nodes for which
                              ! the wave kinematics were generated.  We can use the nodeToLumpedIndx however for L_F_I and it's indices are (timeIndx, nodeIndx, vectorIndx)
                              
                              
                              F_I      = 0.0
                              IF ( (Vmag > 0.0) .AND. (.NOT. nodes(I)%PropPot) ) THEN
                                 af =  WaveAcc(M,commonNodeLst(J),:)
                                 VnDotAf = Dot_Product(Vn,af)
                                 F_I(1:3) = ( nodes(I)%JAxCa*AMfactor*VnDotAf / ( REAL( nCommon, ReKi ) * Vmag ) ) * Vn
                              END IF
                              L_F_I(M, :,nodeToLumpedIndx(commonNodeLst(J)))   =  F_I
                           END DO
                           
                        ELSE
                           ! Should we ever land in here?
                           L_An(:,nodeToLumpedIndx(commonNodeLst(J))) = 0.0
                           L_AM_M(:,:,nodeToLumpedIndx(commonNodeLst(J))) = 0.0
                           L_F_I(:,:,nodeToLumpedIndx(commonNodeLst(J)))   = 0.0
                        END IF
                     
                     END DO
                  
                     usedJointList(nodes(I)%JointIndx) = .TRUE.
                  END IF   !IF ( .NOT. usedJointList(nodes(I)%JointIndx) )
               
                  END IF  !  IF ( nodes(I)%JointIndx /= -1 )
                  
              END IF ! IF ( nodes(I)%NodeType == 3 ) THEN
              
   
         
         END IF   ! ( nodes(I)%JointPos(3) >= z0 )
      
     
            
   END DO   ! I=1,numNodes
   
   
      ! Loop over elements and identify those end nodes which have a JointOvrlp option of 0.
      
   DO I=1,numElements
      
      element = elements(I)
      node1   = nodes(element%Node1Indx)
      node2   = nodes(element%Node2Indx)
      
      CALL GetDistance( node1%JointPos, node2%JointPos, L )
          
      IF ( node1%NodeType == 1 .AND.  node1%JointOvrlp == 0 ) THEN
         
            !Process Lumped loads for this node
         node = node1
         sgn = 1.0
         IF (  ( node%JointPos(3) >= z0 ) .AND. (.NOT. node%PropPot) )THEN
            CALL GenerateLumpedLoads( element%Node1Indx, sgn, node, gravity, MSL2SWL, densWater, NStepWave, WaveDynP, dragConst, F_DP, F_B,  ErrStat, ErrMsg )
            L_F_I(:, :, nodeToLumpedIndx(element%Node1Indx))  = L_F_I(:, :, nodeToLumpedIndx(element%Node1Indx)) + F_DP
           ! L_F_DP(:, :, nodeToLumpedIndx(element%Node1Indx)) = F_DP
            
            L_dragConst(nodeToLumpedIndx(element%Node1Indx))  = dragConst
         IF ( ( node%JointPos(3) >= z0 ) .AND. (.NOT. node%PropPot) )THEN
            L_F_B (:, nodeToLumpedIndx(element%Node1Indx))    = F_B
         END IF
            
         ELSE
            F_BF                                              = 0.0
            !L_F_DP(:, :, nodeToLumpedIndx(element%Node1Indx)) = 0.0
            L_F_B (:, nodeToLumpedIndx(element%Node1Indx))    = 0.0
            L_dragConst(nodeToLumpedIndx(element%Node1Indx))  = 0.0
            
         END IF
         IF ( node%FillFlag ) THEN
            IF ( (node%JointPos(3) <= (node%FillFSLoc))  .AND. (node%JointPos(3) >= z0) ) THEN
               CALL LumpFloodedBuoyancy( sgn, node%FillDensity, node%R, node%t, node%FillFSLoc, node%JointPos(3) , node%R_LToG, gravity, F_BF )      
      
               L_F_BF(:, nodeToLumpedIndx(element%Node1Indx))    = F_BF
            ELSE
               L_F_BF(:, nodeToLumpedIndx(element%Node1Indx))    = 0.0
            END IF
         ELSE
            L_F_BF(:, nodeToLumpedIndx(element%Node1Indx))       = 0.0
         END IF
         
         
      ENDIF
      
      
      IF ( node2%NodeType == 1 .AND.  node2%JointOvrlp == 0 ) THEN
         
            !Process Lumped loads for this node
         node = node2
         sgn = -1.0
         
            ! Generate the loads regardless of node location, and then make the bounds check per load type because the range is different
         CALL GenerateLumpedLoads( element%Node2Indx, sgn, node, gravity, MSL2SWL, densWater, NStepWave, WaveDynP, dragConst, F_DP, F_B, ErrStat, ErrMsg )
         IF ( ( node%JointPos(3) >= z0 ) .AND. (.NOT. node%PropPot) ) THEN
            L_F_I(:, :, nodeToLumpedIndx(element%Node2Indx))  = L_F_I(:, :, nodeToLumpedIndx(element%Node2Indx)) + F_DP
            !L_F_DP(:, :, nodeToLumpedIndx(element%Node2Indx)) = F_DP
            
            L_dragConst(nodeToLumpedIndx(element%Node2Indx))  = dragConst
            
           IF ( ( node%JointPos(3) >= z0 ) .AND. (.NOT. node%PropPot) ) THEN
              L_F_B (:, nodeToLumpedIndx(element%Node2Indx))    = F_B
           END IF
           
         ELSE
            F_BF                                              = 0.0
            !L_F_DP(:, :, nodeToLumpedIndx(element%Node2Indx)) = 0.0
            L_F_B (:, nodeToLumpedIndx(element%Node2Indx))    = 0.0
            L_dragConst(nodeToLumpedIndx(element%Node2Indx))  = 0.0
            
         END IF
         IF ( node%FillFlag ) THEN
            IF ( (node%JointPos(3) <= (node%FillFSLoc))   .AND. (node%JointPos(3) >= z0) ) THEN
               CALL LumpFloodedBuoyancy( sgn, node%FillDensity, node%R, node%t, node%FillFSLoc, node%JointPos(3), node%R_LToG, gravity, F_BF )      
               L_F_BF(:, nodeToLumpedIndx(element%Node2Indx))    = F_BF
            ELSE
               L_F_BF(:, nodeToLumpedIndx(element%Node2Indx))    = 0.0
            END IF
         ELSE
            L_F_BF(:, nodeToLumpedIndx(element%Node2Indx))       = 0.0
         END IF
      ENDIF
      
      
       
      
      IF ( ErrStat /= 0 ) THEN
            RETURN
      END IF 
      
      
   END DO    
      
   
   CALL MeshCommit ( lumpedMeshIn   &
                      , ErrStat            &
                      , ErrMsg             )
   
   IF ( ErrStat /= 0 ) THEN
         RETURN
   END IF 
   
      ! Initialize the inputs
   DO I=1,lumpedMeshIn%Nnodes
      lumpedMeshIn%Orientation(:,:,I) = lumpedMeshIn%RefOrientation(:,:,I)
   END DO
   lumpedMeshIn%TranslationDisp = 0.0
   lumpedMeshIn%TranslationVel  = 0.0
   lumpedMeshIn%RotationVel     = 0.0
   lumpedMeshIn%TranslationAcc  = 0.0
   lumpedMeshIn%RotationAcc     = 0.0
   
   CALL MeshCopy (   SrcMesh      = lumpedMeshIn            &
                     ,DestMesh     = lumpedMeshOut          &
                     ,CtrlCode     = MESH_SIBLING           &
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat                &
                     ,ErrMess      = ErrMsg                 &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 )
   
   lumpedMeshIn%RemapFlag  = .TRUE.
   lumpedMeshOut%RemapFlag = .TRUE.
   
   
END SUBROUTINE CreateLumpedMesh                                 

!subroutine ComputeDistributedLoadsAtNode( elementWaterState, densWater, JointZPos, &
!                                          PropPot, R, dRdz, t, tMG, MGdensity, &
!                                          R_LToG, Ca, Cp, AxCa, AxCp, Cd, WaveAcc, WaveDynP, D_dragConst_in, &
!                                          D_AM_M, D_dragConst, D_F_I, ErrStat, ErrMsg )   
!
!   INTEGER,                                INTENT( IN    )  ::  elementWaterState
!   REAL(ReKi),                             INTENT( IN    )  ::  densWater
!   REAL(ReKi),                             INTENT( IN    )  ::  JointZPos
!   LOGICAL,                                INTENT( IN    )  ::  PropPot
!   REAL(ReKi),                             INTENT( IN    )  ::  R
!   REAL(ReKi),                             INTENT( IN    )  ::  dRdz
!   REAL(ReKi),                             INTENT( IN    )  ::  t
!   REAL(ReKi),                             INTENT( IN    )  ::  tMG
!   REAL(ReKi),                             INTENT( IN    )  ::  MGdensity 
!   REAL(ReKi),                             INTENT( IN    )  ::  R_LToG(3,3)
!   REAL(ReKi),                             INTENT( IN    )  ::  Ca
!   REAL(ReKi),                             INTENT( IN    )  ::  Cp
!   REAL(ReKi),                             INTENT( IN    )  ::  AxCa
!   REAL(ReKi),                             INTENT( IN    )  ::  AxCp
!   REAL(ReKi),                             INTENT( IN    )  ::  Cd
!   REAL(ReKi),                             INTENT( IN    )  ::  WaveAcc(3)
!   REAL(ReKi),                             INTENT( IN    )  ::  WaveDynP
!   REAL(ReKi),                             INTENT( IN    )  ::  D_dragConst_in                   ! 
!   REAL(ReKi),                             INTENT(   OUT )  ::  D_AM_M(6,6)                   ! Added mass of member
!   
!   REAL(ReKi),                             INTENT(   OUT )  ::  D_dragConst                   ! 
!   REAL(ReKi),                             INTENT(   OUT )  ::  D_F_I(3)                      ! Inertial force associated with the member
!   INTEGER,                                INTENT(   OUT )  ::  ErrStat              ! returns a non-zero value when an error occurs  
!   CHARACTER(*),                           INTENT(   OUT )  ::  ErrMsg               ! Error message if ErrStat /= ErrID_None
!   REAL(ReKi)   ::  k(3)
!   
!   
!   IF ( .NOT. PropPot ) THEN        ! Member is not modeled with WAMIT      
!                   
!            
!      ! node is in the water, what about the entire element?
!      IF ( elementWaterState == 1 ) THEN
!                  
!         ! Element is in the water
!   
!         
!            ! For buoyancy calculations we need to adjust the Z-location based on MSL2SWL. If MSL2SWL > 0 then SWL above MSL, and so we need to place the Z value at a deeper position.  
!            !   SWL is at Z=0 for buoyancy calcs, but geometry was specified relative to MSL (MSL2SWL = 0) 
!         k = R_LToG(:,3)
!         CALL DistrInertialLoads( densWater, Ca, Cp, AxCa, AxCp, R, tMG, dRdZ, k, WaveAcc, WaveDynP, D_F_I, ErrStat, ErrMsg  )                
!         CALL DistrAddedMass( densWater, Ca, AxCa, R_LToG, R, tMG, dRdZ, D_AM_M )  
!                 
!      ELSE
!            ! Element is out of the water
!         D_F_I (:)   = 0.0
!         D_AM_M(:,:) = 0.0          ! This is not time-dependent
!      END IF
!                         
!            
!   END IF      ! IF ( .NOT. nodes(I)%PropPot )
!            
!      ! These are the only two loads we compute at initialization if the member is modeled with WAMIT, they are also computed when Morison is used.
!   IF  ( elementWaterState == 1 )THEN 
!         ! element is in the water
!      D_dragConst = D_dragConst_in
!   ELSE
!      D_dragConst = 0.0
!   END IF
!                       
!end subroutine ComputeDistributedLoadsAtNode
                                          
                                  
SUBROUTINE CreateDistributedMesh( densWater, gravity, MSL2SWL, wtrDpth, NStepWave, WaveAcc, WaveDynP, numNodes, nodes, nodeInWater, numElements, elements, &
                                  numDistribMarkers,  distribMeshIn, distribMeshOut, distribToNodeIndx,  D_F_I,      &
                                  D_F_B, D_F_DP, D_F_MG, D_F_BF, D_AM_M, D_AM_MG, D_AM_F, D_dragConst, elementWaterStateArr, &
                                  Morison_Rad, ErrStat, ErrMsg )

   REAL(ReKi),                             INTENT( IN    )  ::  densWater
   REAL(ReKi),                             INTENT( IN    )  ::  gravity
   REAL(ReKi),                             INTENT( IN    )  ::  MSL2SWL
   REAL(ReKi),                             INTENT( IN    )  ::  wtrDpth
   INTEGER,                                INTENT( IN    )  ::  NStepWave
   REAL(SiKi),                             INTENT( IN    )  ::  WaveAcc(0:,:,:)
   REAL(SiKi),                             INTENT( IN    )  ::  WaveDynP(0:,:)
   INTEGER,                                INTENT( IN    )  ::  numNodes
   INTEGER,                                INTENT( IN    )  ::  numElements
   TYPE(Morison_MemberType),               INTENT( IN    )  ::  elements(:)
   TYPE(Morison_NodeType),                 INTENT( IN    )  ::  nodes(:)
   INTEGER(IntKi),                         INTENT( IN    )  ::  nodeInWater(0:,:)   ! Flag indicating whether or not a node is in the water at a given wave time
   INTEGER,                                INTENT(   OUT )  ::  numDistribMarkers
   !TYPE(Morison_NodeType), ALLOCATABLE,    INTENT(   OUT )  ::  distribMarkers(:)
   TYPE(MeshType),                         INTENT(   OUT )  ::  distribMeshIn
   TYPE(MeshType),                         INTENT(   OUT )  ::  distribMeshOut 
   INTEGER, ALLOCATABLE,                   INTENT(   OUT )  ::  distribToNodeIndx(:)
   REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  D_F_I(:,:,:) 
   REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  D_F_B(:,:)                      ! Buoyancy force associated with the member
   REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  D_F_DP(:,:,:)                   ! Dynamic pressure force
   REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  D_F_MG(:,:)                     ! Marine growth weight
   REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  D_F_BF(:,:)                     ! Flooded buoyancy force
   REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  D_AM_M(:,:,:)                   ! Added mass of member
   REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  D_AM_MG(:)                      ! Added mass of marine growth   
   REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  D_AM_F(:)                   ! Added mass of flooded fluid
   REAL(ReKi),ALLOCATABLE,                 INTENT(   OUT)   ::  D_dragConst(:)                   ! 
   INTEGER,ALLOCATABLE,                    INTENT(   OUT)   ::  elementWaterStateArr(:,:)
   REAL(SiKi), ALLOCATABLE,                INTENT(   OUT )  ::  Morison_Rad(:)       ! radius of each node (for FAST visualization)
   INTEGER,                                INTENT(   OUT )  ::  ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                           INTENT(   OUT )  ::  ErrMsg               ! Error message if ErrStat /= ErrID_None
   
            
   
   INTEGER                    ::  I, J, count, node2Indx 
   INTEGER                    ::  secondNodeWaterState
   TYPE(Morison_MemberType)   ::  element
   TYPE(Morison_NodeType)     ::  node1, node2
   REAL(ReKi)                 ::  L
   REAL(ReKi)                 ::  k(3)
   REAL(R8Ki)                 :: orientation(3,3)
  
  ! REAL(ReKi),ALLOCATABLE     ::  F_DP(:,:)
   REAL(ReKi)                 ::  F_B(6)
   REAL(ReKi)                 ::  F_BF(6)
   REAL(ReKi),ALLOCATABLE     :: F_I(:,:)
   REAL(ReKi)                 ::  z0
   INTEGER, ALLOCATABLE       :: nodeToDistribIndx(:)
   
   numDistribMarkers = 0
   z0                = -(wtrDpth) ! The total sea depth is the still water depth of the seabed 
   
      ! Count how many distributed markers we need to create by looping over the nodes
      
   DO I=1,numNodes
      
      IF ( nodes(I)%NodeType /= 3 ) THEN ! exclude super member nodes
            
         numDistribMarkers = numDistribMarkers + 1
      
      END IF
      
   END DO
   
   
      ! Create the input and output meshes associated with distributed loads
      
   CALL MeshCreate( BlankMesh      = distribMeshIn          &
                     ,IOS          = COMPONENT_INPUT        &
                     ,Nnodes       = numDistribMarkers      &
                     ,ErrStat      = ErrStat                &
                     ,ErrMess      = ErrMsg                 &
                     ,TranslationDisp = .TRUE.              &
                     ,Orientation     = .TRUE.              &
                     ,TranslationVel  = .TRUE.              &
                     ,RotationVel     = .TRUE.              &
                     ,TranslationAcc  = .TRUE.              &
                     ,RotationAcc     = .TRUE.               )

   IF ( ErrStat >= AbortErrLev ) RETURN
    
   CALL AllocAry( Morison_Rad, numDistribMarkers, 'Morison_Rad', ErrStat, ErrMsg)
   IF ( ErrStat >= AbortErrLev ) RETURN
   
   
      ! Allocate the array for the distributed markers
      
   !ALLOCATE ( distribMarkers(numDistribMarkers), STAT = ErrStat )
   !IF ( ErrStat /= ErrID_None ) THEN
   !   ErrMsg  = ' Error allocating space for the distributed load markers array.'
   !   ErrStat = ErrID_Fatal
   !   RETURN
   !END IF  
   
   ALLOCATE ( distribToNodeIndx(numDistribMarkers), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the distributed index array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF  
   
   ALLOCATE ( nodeToDistribIndx(numNodes), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the distributed index array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF  
   
  
   
   ALLOCATE ( elementWaterStateArr( 0:NStepWave, numDistribMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the elementWaterStateArr array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   elementWaterStateArr = 0 ! out of the water
   
   ALLOCATE ( D_F_I( 0:NStepWave, 6, numDistribMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the distributed inertial forces/moments array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   D_F_I = 0.0
   
   ALLOCATE ( D_F_B( 6, numDistribMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the distributed buoyancy forces/moments array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   D_F_B = 0.0
   
   ALLOCATE ( D_F_DP( 0:NStepWave, 6, numDistribMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the distributed dynamic pressure forces/moments array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   D_F_DP = 0.0
   
   ALLOCATE ( D_F_MG( 6, numDistribMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the distributed marine growth weight forces/moments array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   D_F_MG = 0.0
   
   ALLOCATE ( D_F_BF( 6, numDistribMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the distributed buoyancy due to flooding forces/moments array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   D_F_BF = 0.0
   
  
    ALLOCATE ( D_AM_M( 3, 3, numDistribMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the distributed added mass of flooded fluid.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   D_AM_M = 0.0
   
   ALLOCATE ( D_AM_MG( numDistribMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the distributed added mass of marine growth.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   D_AM_MG = 0.0
   
   ALLOCATE ( D_AM_F( numDistribMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the distributed added mass of flooded fluid.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   D_AM_F = 0.0
   
   ALLOCATE ( D_dragConst( numDistribMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the distributed drag constants.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   D_dragConst = 0.0
   
   ! Loop over nodes to create all loads on the resulting markers except for the buoyancy loads
   ! For the buoyancy loads, loop over the elements and then apply one half of the resulting value
   ! to each of the interior element nodes but the full value to an end node.  This means that an internal member node will receive 1/2 of its
   ! load from element A and 1/2 from element B.  If it is the end of a member it will simply receive
   ! the element A load.
   
   count = 1 
   
   DO I=1,numNodes
      
      IF ( nodes(I)%NodeType /= 3 ) THEN
         
            ! End point or internal member node
            
            ! Find the node index for the other end of this element
         !IF ( nodes(I)%NodeType == 1 ) THEN
         !   element = elements(nodes(I)%ConnectionList(1))
         !   IF ( element%Node1Indx == I ) THEN
         !      node2Indx = element%Node2Indx
         !   ELSE
         !      node2Indx = element%Node1Indx
         !   END IF
         !ELSE
         !   node2Indx    = -1
         !END IF
         !
         !   ! Need to see if this node is connected to an element which goes above MSL2SWL or below Seabed.
         !IF ( node2Indx > 0 ) THEN
         !   IF ( nodes(node2Indx)%JointPos(3) > MSL2SWL ) THEN
         !      secondNodeWaterState = 1
         !   ELSE IF  ( nodes(node2Indx)%JointPos(3) < z0 ) THEN
         !      secondNodeWaterState = 2
         !   ELSE
         !      secondNodeWaterState = 0
         !   END IF
         !ELSE
         !   secondNodeWaterState = 0
         !END IF
         
         !CALL GetDistance( element
         !   ! Compute all initialization forces now so we have access to the element information
         !   
          IF ( .NOT. nodes(I)%PropPot .AND. nodes(I)%JointPos(3) >= z0 ) THEN
            
                ! Member is not modeled with WAMIT
            
            k =  nodes(I)%R_LToG(:,3)
            
          !  IF ( nodes(I)%JointPos(3) <= MSL2SWL .AND. nodes(I)%JointPos(3) >= z0 ) THEN
               
               
               CALL DistrAddedMass( densWater, nodes(I)%Ca, nodes(I)%AxCa, nodes(I)%R_LToG, nodes(I)%R, nodes(I)%tMG, nodes(I)%dRdZ, D_AM_M(:,:,count) )  
             !  IF ( secondNodeWaterState == 0 ) THEN
                     ! Element is in the water
                     
                  
                  !CALL DistrDynPressure( I, nodes(I)%AxCa,nodes(I)%AxCp, nodes(I)%R_LToG, nodes(I)%R, nodes(I)%tMG, nodes(I)%dRdz, NStepWave, WaveDynP, WaveAcc, F_DP, ErrStat, ErrMsg)
             !     D_F_DP(:,:,count) = 0.0 !F_DP
                     ! For buoyancy calculations we need to adjust the Z-location based on MSL2SWL. If MSL2SWL > 0 then SWL above MSL, and so we need to place the Z value at a deeper position.  
                     !   SWL is at Z=0 for buoyancy calcs, but geometry was specified relative to MSL (MSL2SWL = 0) 
                  CALL DistrBuoyancy( densWater, nodes(I)%R, nodes(I)%tMG, nodes(I)%dRdz, nodes(I)%JointPos(3) - MSL2SWL, nodes(I)%R_LToG, gravity, F_B  ) 
                  D_F_B(:,count)    = F_B
              !   ! Compute all initialization forces now so we have access to the element information
         !   
       ! IF ( ( .NOT. nodes(J)%PropPot ) .AND. ( nodes(J)%JointPos(3) >= z0 ) ) THEN
       !    k =  nodes(J)%R_LToG(:,3)
            CALL DistrInertialLoads( I, densWater, nodes(I)%Ca, nodes(I)%Cp, nodes(I)%AxCa, nodes(I)%AxCp, nodes(I)%R, nodes(I)%tMG, nodes(I)%dRdZ, k, NStepWave, WaveAcc, WaveDynP, F_I, ErrStat, ErrMsg  )              
            D_F_I(:,:,count)  = F_I     
            
       ! END IF 
                  
           !    ELSE
                     ! Element is out of the water
          !        D_F_DP(:,:,count) = 0.0
           !       D_F_B(:,count)    = 0.0
                 
           !    END IF
               
          !  ELSE 
               ! NOTE: Everything was initialized to zero so this isn't really necessary. GJH 9/24/13
             
        !       D_F_DP(:,:,count) = 0.0
               
         !      D_F_B(:,count)    = 0.0
         !   END IF
            
        !    IF ( ( nodes(I)%JointPos(3) >= z0 ) .AND. (secondNodeWaterState /= 2 ) ) THEN
                  ! if the node is at or above the seabed then the element is in the water
               CALL DistrMGLoads( nodes(I)%MGdensity, gravity, nodes(I)%R, nodes(I)%tMG, D_F_MG(:,count) )            
               CALL DistrAddedMassMG( nodes(I)%MGdensity, nodes(I)%R, nodes(I)%tMG, D_AM_MG(count) )
       !     ELSE
       !        D_F_MG(:,count)   = 0.0
       !        D_AM_MG(count)= 0.0
       !     END IF
            
          END IF      ! IF ( .NOT. nodes(I)%PropPot )
            
          ! This is always computed, but may be zereod out for any given timestep during the CalcOutput work
         CALL DistrDragConst( densWater, nodes(I)%Cd, nodes(I)%R, nodes(I)%tMG, D_dragConst(count) ) 
         
         IF ( nodes(I)%FillFlag ) THEN
            IF ( nodes(I)%JointPos(3) <= nodes(I)%FillFSLoc   .AND. nodes(I)%JointPos(3) >= z0 ) THEN
               
                           ! Find the node index for the other end of this element
               IF ( nodes(I)%NodeType == 1 ) THEN
                  element = elements(nodes(I)%ConnectionList(1))
                  IF ( element%Node1Indx == I ) THEN
                     node2Indx = element%Node2Indx
                  ELSE
                     node2Indx = element%Node1Indx
                  END IF
               ELSE
                  node2Indx    = -1
               END IF
               
                  ! different check for filled element, based on free-surface location
               IF ( node2Indx > 0 ) THEN
                  IF ( nodes(node2Indx)%JointPos(3) > nodes(I)%FillFSLoc ) THEN
                     secondNodeWaterState = 0
                  ELSE IF  ( nodes(node2Indx)%JointPos(3) < z0 ) THEN
                     secondNodeWaterState = 2
                  ELSE
                     secondNodeWaterState = 1
                  END IF
               ELSE
                  secondNodeWaterState = 1
               END IF
               
               IF (secondNodeWaterState == 1 ) THEN
                  CALL DistrAddedMassFlood( nodes(I)%FillDensity, nodes(I)%R, nodes(I)%t, D_AM_F(count) )
                     ! For buoyancy calculations we need to adjust the Z-location based on MSL2SWL. If MSL2SWL > 0 then SWL above MSL, and so we need to place the Z value at a deeper position.  
                     !   SWL is at Z=0 for buoyancy calcs, but geometry was specified relative to MSL (MSL2SWL = 0) 
                  CALL DistrFloodedBuoyancy( nodes(I)%FillDensity, nodes(I)%FillFSLoc, nodes(I)%R, nodes(I)%t, nodes(I)%dRdZ, nodes(I)%JointPos(3) - MSL2SWL, nodes(I)%R_LToG, gravity, F_BF )
                  D_F_BF(:,count  ) = F_BF
               ELSE
                  D_AM_F(count) = 0.0
                  D_F_BF(:,count  ) = 0.0
               END IF
               
            ELSE
               ! NOTE: Everything was initialized to zero so this isn't really necessary. GJH 9/24/13
               D_AM_F(count) = 0.0
               D_F_BF(:,count  ) = 0.0
            END IF      
         ELSE
               D_AM_F(count) = 0.0
               D_F_BF(:,count  ) = 0.0
         END IF
         
         
         
            ! Create the node on the mesh
            
         orientation = transpose(nodes(I)%R_LToG )
         CALL MeshPositionNode (distribMeshIn           &
                              , count                   &
                              , nodes(I)%JointPos       &  ! this info comes from FAST
                              , ErrStat                 &
                              , ErrMsg                  & ! , orient = orientation  & ! bjj: I need this orientation set for visualization. but, will it mess up calculations (because loads are in different positions?)
                              ) !, transpose(nodes(I)%R_LToG )     )
         IF ( ErrStat >= AbortErrLev ) RETURN
         
         Morison_Rad(count) = nodes(I)%R   ! set this for FAST visualization
         
         distribToNodeIndx(count) = I
         nodeToDistribIndx(I) = count
         count = count + 1    
         
      END IF
      
   END DO
   
   ! Now for time-varying values
    
   DO count=1,numDistribMarkers
     J = distribToNodeIndx(count)
     DO I=0,NStepWave    
      IF ( nodes(J)%NodeType /= 3 ) THEN
         
            ! End point or internal member node
            
            ! Find the node index for the other end of this element
         IF ( nodes(J)%NodeType == 1 ) THEN
            element = elements(nodes(J)%ConnectionList(1))
            IF ( element%Node1Indx == J ) THEN
               node2Indx = element%Node2Indx
            ELSE
               node2Indx = element%Node1Indx
            END IF
         ELSE
            node2Indx    = -1
         END IF
      
            ! Need to see if this node is connected to an element which goes above MSL2SWL or below Seabed.
         IF ( node2Indx > 0 ) THEN
            IF ( ( nodeInWater(I,node2Indx) == 0 ) .AND. ( nodes(node2Indx)%JointPos(3) >= z0 ) ) THEN
               secondNodeWaterState = 0
            ELSE IF  ( nodes(node2Indx)%JointPos(3) < z0 ) THEN
               secondNodeWaterState = 2
            ELSE
               secondNodeWaterState = 1
            END IF
         ELSE
            secondNodeWaterState = 1
         END IF
         
         
         
     
         IF ( (nodeInWater(I,J) == 1) .AND. nodes(J)%JointPos(3) >= z0 ) THEN
            IF ( secondNodeWaterState == 1 ) THEN
                  ! Element is in the water                   
               elementWaterStateArr(I,count) = 1
            END IF            
         END IF             
         
      END IF
      
     END DO   ! DO I=0,NStepWave  
   END DO     ! DO count=1,numDistribMarkers
                                  
   
  ! End of time-varying values
   
   
   DO I=1,numElements
   
       
         ! Create the mesh element
         
      CALL MeshConstructElement (  distribMeshIn   &
                                  , ELEMENT_LINE2       &                                 
                                  , ErrStat            &
                                  , ErrMsg  &
                                  , nodeToDistribIndx(elements(I)%Node1Indx)                  &
                                  , nodeToDistribIndx(elements(I)%Node2Indx)                )
      
      IF ( ErrStat /= 0 )    RETURN
      
 
   END DO
   
   
   CALL MeshCommit ( distribMeshIn   &
                      , ErrStat            &
                      , ErrMsg             )
   
   IF ( ErrStat /= 0 ) THEN
         RETURN
   END IF 
   
      ! Initialize the inputs
   DO I=1,distribMeshIn%Nnodes
      distribMeshIn%Orientation(:,:,I) = distribMeshIn%RefOrientation(:,:,I)
   END DO
   distribMeshIn%TranslationDisp = 0.0
   distribMeshIn%TranslationVel  = 0.0
   distribMeshIn%RotationVel     = 0.0
   distribMeshIn%TranslationAcc  = 0.0
   distribMeshIn%RotationAcc     = 0.0
   
   CALL MeshCopy (    SrcMesh      = distribMeshIn &
                     ,DestMesh     = distribMeshOut         &
                     ,CtrlCode     = MESH_SIBLING           &
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat                &
                     ,ErrMess      = ErrMsg                 &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 )

   distribMeshIn%RemapFlag  = .TRUE.
   distribMeshOut%RemapFlag = .TRUE.
   
END SUBROUTINE CreateDistributedMesh
                                  
  

!====================================================================================================
SUBROUTINE Morison_ProcessMorisonGeometry( InitInp, ErrStat, ErrMsg )
!     This public subroutine process the input geometry and parameters and eliminates joint overlaps,  
!     sub-divides members, sets joint-level properties, etc.
!----------------------------------------------------------------------------------------------------  

      ! Passed variables
   
   TYPE(Morison_InitInputType),   INTENT( INOUT )   :: InitInp              ! the Morison initialization data 
   !TYPE(Morison_ParameterType),   INTENT( INOUT )   :: p                    ! tge Morison parameter data
   INTEGER,                       INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                  INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None

  
      ! Local variables
         
   INTEGER                                      :: I    !, J, j1, j2, tempINT                ! generic integer for counting
!   TYPE(Morison_JointType)                      :: joint1, joint2                                   
!   Real(ReKi)                                   :: z1
!   Real(ReKi)                                   :: z2
   Real(ReKi)                                   :: d
   INTEGER                                      :: temp
   INTEGER                                      :: prop1Indx, prop2Indx, node1Indx, node2Indx
   INTEGER                                      :: maxNodes        = 0
   INTEGER                                      :: maxElements     = 0
   INTEGER                                      :: maxSuperMembers = 0
 !  TYPE(Morison_NodeType)                       :: node1, node2, tempNode
   TYPE(Morison_MemberPropType)                 :: propSet
   INTEGER                                      :: numSplitNodes
   TYPE(Morison_NodeType),ALLOCATABLE           :: splitNodes(:)
   LOGICAL                                      :: doSwap
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = "" 
   
   
   IF ( InitInp%NMembers > 0 ) THEN
      
      
         ! Determine the maximum number of nodes,  elements, and super members which might be generated for the simulation mesh
      CALL GetMaxSimQuantities( InitInp%NMGDepths, InitInp%MGTop, InitInp%MGBottom, InitInp%MSL2SWL, -InitInp%WtrDpth, InitInp%FilledGroups, InitInp%NJoints, InitInp%InpJoints, InitInp%NMembers, InitInp%InpMembers, maxNodes, maxElements, maxSuperMembers )
  
  
      ! Create a worse case size for the number of nodes and number of elements that will be generated for the simulation
      ! marine growth split + super member split + member subdivision all creates new nodes
         
      ! marine growth split + member subdivision creates new elements
      
      ! Create a worse case size for the number of super members
      
         ! 1) Let's start by generating a mirror of the input mesh (joints and members) as the initial version of the simulation mesh
         ! In doing so, create the initial mapping between the input mesh and this current version of the simulation mesh
         
         
         ! Allocate memory for Joint-related arrays
         
      InitInp%NNodes = InitInp%NJoints
      
      ALLOCATE ( InitInp%Nodes(maxNodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for Nodes array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF    
          
      
      DO I = 1,InitInp%NNodes
         ! Copy all necessary data from the input joints to these node data structures
         InitInp%Nodes(I)%JointPos       = InitInp%InpJoints(I)%JointPos
         InitInp%Nodes(I)%JointAxIDIndx  = InitInp%InpJoints(I)%JointAxIDIndx
         InitInp%Nodes(I)%JointOvrlp     = InitInp%InpJoints(I)%JointOvrlp
         InitInp%Nodes(I)%NConnections   = InitInp%InpJoints(I)%NConnections
         InitInp%Nodes(I)%ConnectionList = InitInp%InpJoints(I)%ConnectionList
         InitInp%Nodes(I)%JointIndx      = I
         InitInp%Nodes(I)%NodeType       = 1  ! 1 = end of a member, 2 = interior of a member, 3 = super member node
         InitInp%Nodes(I)%FillFSLoc      = InitInp%MSL2SWL  
         InitInp%Nodes(I)%FillFlag       = .FALSE.
         InitInp%Nodes(I)%FillDensity    = 0.0
         
         
         
      END DO
      
      
          ! Allocate memory for Members arrays
          
      InitInp%NElements = InitInp%NMembers  
      
      ALLOCATE ( InitInp%Elements(maxElements), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for Elements array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
          
      
      DO I = 1,InitInp%NMembers  
         
         InitInp%Elements(I)%Node1Indx = InitInp%InpMembers(I)%MJointID1Indx              ! Index of  the first node in the Morison_NodeType array
         InitInp%Elements(I)%Node2Indx = InitInp%InpMembers(I)%MJointID2Indx              ! Index of  the second node in the Morison_NodeType array
         node1Indx                     = InitInp%Elements(I)%Node1Indx
         node2Indx                     = InitInp%Elements(I)%Node2Indx
         prop1Indx = InitInp%InpMembers(I)%MPropSetID1Indx
         prop2Indx = InitInp%InpMembers(I)%MPropSetID2Indx
         
            ! Make sure that Node1 has the lower Z value, re-order if necessary
            ! We need to do this because the local element coordinate system is defined such that the first node is located with a smaller global Z value
            ! than the second node.
            ! The local element coordinate system requires that Z1 <= Z2, and if Z1=Z2 then X1 <= X2, and if Z1=Z2, X1=X2 then Y1<Y2
   
         InitInp%Elements(I)%InpMbrDist1         = 0.0
         InitInp%Elements(I)%InpMbrDist2         = 1.0
         doSwap = .FALSE.
                          
         IF ( EqualRealNos(InitInp%Nodes(node1Indx)%JointPos(3), InitInp%Nodes(node2Indx)%JointPos(3) ) ) THEN         ! Z1 = Z2          
            IF ( EqualRealNos(InitInp%Nodes(node1Indx)%JointPos(1), InitInp%Nodes(node2Indx)%JointPos(1) ) ) THEN      ! X1 = X2
               IF   ( InitInp%Nodes(node1Indx)%JointPos(2) > InitInp%Nodes(node2Indx)%JointPos(2) ) THEN
                  doSwap = .TRUE.  ! Y1 > Y2
               END IF
            ELSE IF ( InitInp%Nodes(node1Indx)%JointPos(1) > InitInp%Nodes(node2Indx)%JointPos(1) ) THEN
               doSwap = .TRUE.  ! X1 > X2
            END IF
         ELSE IF    ( InitInp%Nodes(node1Indx)%JointPos(3) > InitInp%Nodes(node2Indx)%JointPos(3) ) THEN
            doSwap = .TRUE.                                ! Z1 > Z2  
         END IF
         
         IF ( doSwap ) THEN
            
               ! Swap node indices to satisfy orientation rules for element nodes
            
            InitInp%Elements(I)%Node1Indx = InitInp%InpMembers(I)%MJointID2Indx              
            InitInp%Elements(I)%Node2Indx = InitInp%InpMembers(I)%MJointID1Indx  
            node1Indx                     = InitInp%Elements(I)%Node1Indx
            node2Indx                     = InitInp%Elements(I)%Node2Indx
            temp = prop1Indx
            prop1Indx = prop2Indx
            prop2Indx = temp
            InitInp%Elements(I)%InpMbrDist1         = 1.0
            InitInp%Elements(I)%InpMbrDist2         = 0.0
            
         END IF
         
         propSet = InitInp%MPropSets(prop1Indx)
         InitInp%Elements(I)%R1               = propSet%PropD / 2.0
         InitInp%Elements(I)%t1               = propSet%PropThck
         
         propSet = InitInp%MPropSets(prop2Indx)
         InitInp%Elements(I)%R2               = propSet%PropD / 2.0
         InitInp%Elements(I)%t2               = propSet%PropThck 
         
         InitInp%Elements(I)%NumSplits        = InitInp%InpMembers(I)%NumSplits
         InitInp%Elements(I)%Splits        = InitInp%InpMembers(I)%Splits
         !InitInp%Elements(I)%MGSplitState     = InitInp%InpMembers(I)%MGSplitState
         !InitInp%Elements(I)%WtrSplitState     = InitInp%InpMembers(I)%WtrSplitState
         InitInp%Elements(I)%MDivSize         = InitInp%InpMembers(I)%MDivSize
         InitInp%Elements(I)%MCoefMod         = InitInp%InpMembers(I)%MCoefMod
         InitInp%Elements(I)%MmbrCoefIDIndx   = InitInp%InpMembers(I)%MmbrCoefIDIndx
         InitInp%Elements(I)%MmbrFilledIDIndx = InitInp%InpMembers(I)%MmbrFilledIDIndx
      
         CALL GetDistance( InitInp%Nodes(node1Indx)%JointPos, InitInp%Nodes(node2Indx)%JointPos, d)
         
         InitInp%Elements(I)%InpMbrLen           = d
         InitInp%Elements(I)%InpMbrIndx          = I
         
            ! Direction cosines matrix which transforms a point in member coordinates to the global inertial system
         !CALL Morison_DirCosMtrx( node1%JointPos, node2%JointPos, InitInp%Elements(I)%R_LToG  )    
         
        
         InitInp%Elements(I)%PropPot  =  InitInp%InpMembers(I)%PropPot                  ! Flag specifying whether member is modelled in WAMIT [true = modelled in WAMIT, false = not modelled in WAMIT]
         
         
        
         
            ! Calculate the element-level direction cosine matrix and attach it to the entry in the elements array
            
         CALL Morison_DirCosMtrx( InitInp%Nodes(node1Indx)%JointPos, InitInp%Nodes(node2Indx)%JointPos, InitInp%Elements(I)%R_LToG )
        ! InitInp%Nodes(node1Indx)%R_LToG = InitInp%Elements(I)%R_LToG
        ! InitInp%Nodes(node2Indx)%R_LToG = InitInp%Elements(I)%R_LToG
      END DO
      
      
      
         ! Set the fill properties onto the elements
         
      CALL SetElementFillProps( InitInp%NFillGroups, InitInp%FilledGroups, InitInp%NElements, InitInp%Elements )
    
         ! Split elements
      CALL SplitElements(InitInp%NNodes, InitInp%Nodes, InitInp%NElements, InitInp%Elements, ErrStat, ErrMsg)
      
         ! Split element due to MSL2SWL location and seabed location
      !CALL SplitElementsForWtr(InitInp%MSL2SWL, -InitInp%WtrDpth, InitInp%NNodes, InitInp%Nodes, InitInp%NElements, InitInp%Elements, ErrStat, ErrMsg)
      
         ! Split elements if they cross the marine growth boundary. 
         
      !CALL SplitElementsForMG(InitInp%MGTop, InitInp%MGBottom, InitInp%NNodes, InitInp%Nodes, InitInp%NElements, InitInp%Elements, ErrStat, ErrMsg)
      
      
         ! Create any Super Members
      !CALL CreateSuperMembers( InitInp%NNodes, InitInp%Nodes, InitInp%NElements, InitInp%Elements, ErrStat, ErrMsg )
      
         ! Subdivide the members based on user-requested maximum division sizes (MDivSize)
         
      CALL SubdivideMembers( InitInp%NNodes, InitInp%Nodes, InitInp%NElements, InitInp%Elements, ErrStat, ErrMsg )  
      
    
         
         ! Set the element Cd, Ca, and Cp coefs
         
      CALL SetElementCoefs( InitInp%SimplCd, InitInp%SimplCdMG, InitInp%SimplCa, InitInp%SimplCaMG, InitInp%SimplCp, InitInp%SimplCpMG, InitInp%SimplAxCa, InitInp%SimplAxCaMG, InitInp%SimplAxCp, InitInp%SimplAxCpMG,InitInp%CoefMembers, InitInp%NCoefDpth, InitInp%CoefDpths, InitInp%NNodes, InitInp%Nodes, InitInp%NElements, InitInp%Elements )   
      
      
         ! Set the axial coefs AxCd and AxCa
     CALL SetAxialCoefs( InitInp%NJoints, InitInp%NAxCoefs, InitInp%AxialCoefs, InitInp%NNodes, InitInp%Nodes, InitInp%NElements, InitInp%Elements )      
      
         ! Set the marine growth thickness and density information onto the nodes (this is not a per-element quantity, but a per-node quantity
         
      CALL SetNodeMG( InitInp%NMGDepths, InitInp%MGDepths, InitInp%NNodes, InitInp%Nodes )
      
      
         ! Create duplicate nodes at the ends of elements so that only one element is connected to any given end node
         
      CALL SplitMeshNodes( InitInp%NNodes, InitInp%Nodes, InitInp%NElements, InitInp%Elements, numSplitNodes, splitNodes, ErrStat, ErrMsg )
      
      IF (numSplitNodes > InitInp%NNodes ) THEN
         
         InitInp%NNodes = numSplitNodes
         !Reallocate the Nodes array
         DEALLOCATE ( InitInp%Nodes )
         ALLOCATE ( InitInp%Nodes(numSplitNodes), STAT = ErrStat )
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Error allocating space for Nodes array.'
            ErrStat = ErrID_Fatal
            RETURN
         END IF
         InitInp%Nodes = splitNodes
         DEALLOCATE ( splitNodes )
         
      END IF
      
      
         ! Now that the nodes are split, we can push the element properties down to the individual nodes without an issue
         
      CALL SetSplitNodeProperties( InitInp%NNodes, InitInp%Nodes, InitInp%NElements, InitInp%Elements, ErrStat, ErrMsg ) 
      
      
      
         
         ! 6) Store information necessary to compute the user-requested member outputs and joint outputs.  The requested output locations
         !    may be located in between two simulation nodes, so quantities will need to be interpolated. qOutput = q1*s + q2*(1-s), where 0<= s <= 1.
         
         ! NOTE: since we need to mantain the input geometry, the altered members are now part of the simulation mesh and 
         !       we will generate a mapping between the input and simulation meshes which is needed to generate user-requested outputs.
   
    
         
   ELSE  
      
      
         ! No Morison elements, so no processing is necessary, but set nodes and elements to 0.
         
     ! p%NMorisonNodes    = 0  
    !  p%NMorisonElements = 0
      
   END IF
   
END SUBROUTINE Morison_ProcessMorisonGeometry

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps. 
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE Morison_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(Morison_InitInputType),       INTENT(INOUT)  :: InitInp     !< Input data for initialization routine !intent out because of MOVE_ALLOC
      TYPE(Morison_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
      TYPE(Morison_ParameterType),       INTENT(  OUT)  :: p           !< Parameters      
      TYPE(Morison_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
      TYPE(Morison_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
      TYPE(Morison_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
      TYPE(Morison_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states            
      TYPE(Morison_OutputType),          INTENT(  OUT)  :: y           !< Initial system outputs (outputs are not calculated; 
                                                                       !!   only the output mesh is initialized)
      TYPE(Morison_MiscVarType),         INTENT(  OUT)  :: m           !< Initial misc/optimization variables            
      REAL(DbKi),                        INTENT(INOUT)  :: Interval    !< Coupling interval in seconds: the rate that 
                                                                       !!   (1) Morison_UpdateStates() is called in loose coupling &
                                                                       !!   (2) Morison_UpdateDiscState() is called in tight coupling.
                                                                       !!   Input is the suggested time from the glue code; 
                                                                       !!   Output is the actual coupling interval that will be used 
                                                                       !!   by the glue code.
      TYPE(Morison_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      
     ! TYPE(Morison_InitInputType)                       :: InitLocal   ! Local version of the input data for the geometry processing routine
!      INTEGER, ALLOCATABLE                                          :: distribToNodeIndx(:)
!      INTEGER, ALLOCATABLE                                          :: lumpedToNodeIndx(:)
      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Initialize the NWTC Subroutine Library
         
      CALL NWTC_Init(  )


     ! InitLocal = InitInp
      p%WtrDens    = InitInp%WtrDens
      p%NumOuts    = InitInp%NumOuts
      p%NMOutputs  = InitInp%NMOutputs                       ! Number of members to output [ >=0 and <10]
      p%OutSwtch   = InitInp%OutSwtch
      ALLOCATE ( p%MOutLst(p%NMOutputs), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for MOutLst array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
IF (ALLOCATED(InitInp%MOutLst) ) &
      p%MOutLst =    InitInp%MOutLst           ! Member output data
      
      p%NJOutputs = InitInp%NJOutputs                        ! Number of joints to output [ >=0 and <10]
      
      ALLOCATE ( p%JOutLst(p%NJOutputs), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for JOutLst array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
IF (ALLOCATED(InitInp%JOutLst) ) &
      p%JOutLst =    InitInp%JOutLst            ! Joint output data
      
     
      
      
       
      p%NNodes   = InitInp%NNodes
      
      ALLOCATE ( p%Nodes(p%NNodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for Nodes array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      
      p%Nodes    = InitInp%Nodes
      
      
      p%NStepWave= InitInp%NStepWave
      
      ALLOCATE ( p%WaveVel(0:p%NStepWave, p%NNodes, 3), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for wave velocities array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      p%WaveVel = InitInp%WaveVel      
      
      ALLOCATE ( p%WaveAcc(0:p%NStepWave, p%NNodes, 3), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for wave accelerations array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      p%WaveAcc = InitInp%WaveAcc
      
       ALLOCATE ( p%WaveDynP(0:p%NStepWave, p%NNodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for wave dynamic pressure array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      p%WaveDynP = InitInp%WaveDynP
      
      
      
      ALLOCATE ( p%WaveTime(0:p%NStepWave), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for wave time array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF     
      p%WaveTime     = InitInp%WaveTime

      
      CALL MOVE_ALLOC( InitInp%nodeInWater, p%nodeInWater )   
      
      
      
         ! Use the processed geometry information to create the distributed load mesh and associated marker parameters
         
         ! We are storing the parameters in the DistribMarkers data structure instead of trying to hold this information within the DistribMesh.  But these two data structures
         ! must always be in sync.  For example, the 5th element of the DistribMarkers array must correspond to the 5th node in the DistribMesh data structure.
       
      CALL CreateDistributedMesh( InitInp%WtrDens, InitInp%Gravity, InitInp%MSL2SWL, InitInp%WtrDpth, InitInp%NStepWave, InitInp%WaveAcc, InitInp%WaveDynP, &
                                  p%NNodes, p%Nodes, p%nodeInWater, InitInp%NElements, InitInp%Elements, &
                                  p%NDistribMarkers, u%DistribMesh, y%DistribMesh, p%distribToNodeIndx, p%D_F_I, &
                                  p%D_F_B, p%D_F_DP, p%D_F_MG, p%D_F_BF, p%D_AM_M, p%D_AM_MG, p%D_AM_F, p%D_dragConst, p%elementWaterState, &                 ! 
                                  InitOut%Morison_Rad,  ErrStat, ErrMsg )
                                    
                                 

      
                                 
     IF ( ErrStat > ErrID_None ) RETURN
     
         
     CALL CreateLumpedMesh( InitInp%WtrDens, InitInp%Gravity, InitInp%MSL2SWL, InitInp%WtrDpth, InitInp%NStepWave, InitInp%WaveDynP, InitInp%WaveAcc,p%NNodes, p%Nodes, InitInp%NElements, InitInp%Elements, &
                                  p%NLumpedMarkers,  u%LumpedMesh, y%LumpedMesh, p%lumpedToNodeIndx,   p%L_An,     &
                                  p%L_F_B, p%L_F_I, p%L_F_BF, p%L_AM_M, p%L_dragConst, &
                                  ErrStat, ErrMsg )
     IF ( ErrStat > ErrID_None ) RETURN
     
      
      
     ! CALL CreateSuperMesh( InitInp%NNodes, InitInp%Nodes, InitInp%NElements, InitInp%Elements, p%NSuperMarkers, p%SuperMarkers, InitOut%LumpedMesh, ErrStat, ErrMsg )
      
     
     
      
      
         ! Define parameters here:
       
     
      p%DT  = Interval


         ! Define initial system states here:

      x%DummyContState           = 0
      xd%DummyDiscState          = 0
      z%DummyConstrState         = 0
      OtherState%DummyOtherState = 0
      m%LastIndWave              = 1

   IF ( p%OutSwtch > 0 ) THEN
      ALLOCATE ( m%D_F_D(3,y%DistribMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for D_F_D array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%D_F_D = 0.0_ReKi
      ALLOCATE ( m%D_F_I(3,y%DistribMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for D_F_I array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%D_F_I = 0.0_ReKi
      ALLOCATE ( m%D_F_B(6,y%DistribMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for D_F_B array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF      
      m%D_F_B = 0.0_ReKi
      ALLOCATE ( m%D_F_AM(6,y%DistribMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for D_F_AM array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%D_F_AM = 0.0_ReKi
      ALLOCATE ( m%D_F_AM_M(6,y%DistribMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for D_F_AM_M array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%D_F_AM_M = 0.0_ReKi
      ALLOCATE ( m%D_F_AM_MG(6,y%DistribMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for D_F_AM_MG array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%D_F_AM_MG = 0.0_ReKi
      ALLOCATE ( m%D_F_AM_F(6,y%DistribMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for D_F_AM_F array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%D_F_AM_F = 0.0_ReKi
      ALLOCATE ( m%D_FV(3,y%DistribMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for D_FV array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%D_FV = 0.0_ReKi
      ALLOCATE ( m%D_FA(3,y%DistribMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for D_FA array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%D_FA = 0.0_ReKi
      ALLOCATE ( m%D_FDynP(y%DistribMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for D_FDynP array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%D_FDynP = 0.0_ReKi
      
      ALLOCATE ( m%L_F_B(6,y%LumpedMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for L_F_B array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%L_F_B = 0.0_ReKi
      ALLOCATE ( m%L_F_D(3,y%LumpedMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for L_F_D array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%L_F_D = 0.0_ReKi
      ALLOCATE ( m%L_F_I(6,y%LumpedMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for L_F_I array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%L_F_I = 0.0_ReKi
      !ALLOCATE ( m%L_F_DP(6,y%LumpedMesh%Nnodes), STAT = ErrStat )
      !IF ( ErrStat /= ErrID_None ) THEN
      !   ErrMsg  = ' Error allocating space for L_F_DP array.'
      !   ErrStat = ErrID_Fatal
      !   RETURN
      !END IF
      ALLOCATE ( m%L_FV(3,y%LumpedMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for L_FV array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%L_FV = 0.0_ReKi
      ALLOCATE ( m%L_FA(3,y%LumpedMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for L_FA array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%L_FA = 0.0_ReKi
      ALLOCATE ( m%L_FDynP(y%LumpedMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for L_FDynP array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%L_FDynP = 0.0_ReKi
      ALLOCATE ( m%L_F_AM(6,y%LumpedMesh%Nnodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for L_F_AM array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      m%L_F_AM = 0.0_ReKi
      
         ! Define initial guess for the system inputs here:

  !    u%DummyInput = 0


         ! Define system output initializations (set up mesh) here:
  
         
         ! Define initialization-routine output here:
         
         ! Initialize the outputs
         
      CALL MrsnOUT_Init( InitInp, y, p, InitOut, ErrStat, ErrMsg )
      IF ( ErrStat > ErrID_None ) RETURN
      
         ! Determine if we need to perform output file handling
      
      IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN  
         CALL MrsnOUT_OpenOutput( Morison_ProgDesc%Name, TRIM(InitInp%OutRootName)//'.HD', p, InitOut, ErrStat, ErrMsg )
         IF ( ErrStat > ErrID_None ) RETURN
      END IF
      
   END IF  
   
   
      ! Write Summary information now that everything has been initialized.
      
   CALL WriteSummaryFile( InitInp%UnSum, InitInp%MSL2SWL, InitInp%WtrDpth, InitInp%NNodes, InitInp%Nodes, InitInp%NElements, InitInp%Elements, p%NumOuts, p%OutParam, p%NMOutputs, p%MOutLst, p%distribToNodeIndx, p%NJOutputs, p%JOutLst, u%LumpedMesh, y%LumpedMesh,u%DistribMesh, y%DistribMesh, p%L_F_B, p%L_F_BF, p%D_F_B, p%D_F_BF, p%D_F_MG, InitInp%Gravity, ErrStat, ErrMsg ) !p%NDistribMarkers, distribMarkers, p%NLumpedMarkers, lumpedMarkers,
   IF ( ErrStat > ErrID_None ) RETURN  
      
         ! If you want to choose your own rate instead of using what the glue code suggests, tell the glue code the rate at which
         !   this module must be called here:
         
       !Interval = p%DT                                               
   !Contains:
   !   SUBROUTINE CleanUpInitOnErr
   !   IF (ALLOCATED(sw(1)%array))  DEALLOCATE(sw(1)%array, STAT=aviFail)
   !   END SUBROUTINE

END SUBROUTINE Morison_Init


!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE Morison_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(Morison_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(Morison_ParameterType),       INTENT(INOUT)  :: p           !< Parameters     
      TYPE(Morison_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(Morison_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(Morison_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(Morison_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states            
      TYPE(Morison_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(Morison_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables            
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Place any last minute operations or calculations here:


         ! Close files here:     
                  
 

         ! Destroy the input data:
         
      CALL Morison_DestroyInput( u, ErrStat, ErrMsg )


         ! Determine if we need to close the output file
         
      IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN   
         CALL MrsnOut_CloseOutput( p, ErrStat, ErrMsg )         
      END IF 
         
         ! Destroy the parameter data:
         
      
      CALL Morison_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:
         
      CALL Morison_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL Morison_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL Morison_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL Morison_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
         
      CALL Morison_DestroyMisc( m, ErrStat, ErrMsg )

         ! Destroy the output data:
         
      CALL Morison_DestroyOutput( y, ErrStat, ErrMsg )


      

END SUBROUTINE Morison_End
!----------------------------------------------------------------------------------------------------------------------------------
!> This is a loose coupling routine for solving constraint states, integrating continuous states, and updating discrete and other 
!! states. Continuous, constraint, discrete, and other states are updated to values at t + Interval.
SUBROUTINE Morison_UpdateStates( Time, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................
   
      REAL(DbKi),                         INTENT(IN   ) :: Time        !< Current simulation time in seconds
      TYPE(Morison_InputType),            INTENT(IN   ) :: u           !< Inputs at Time                    
      TYPE(Morison_ParameterType),        INTENT(IN   ) :: p           !< Parameters                              
      TYPE(Morison_ContinuousStateType),  INTENT(INOUT) :: x           !< Input: Continuous states at Time; 
                                                                       !!   Output: Continuous states at Time + Interval
      TYPE(Morison_DiscreteStateType),    INTENT(INOUT) :: xd          !< Input: Discrete states at Time; 
                                                                       !!   Output: Discrete states at Time + Interval
      TYPE(Morison_ConstraintStateType),  INTENT(INOUT) :: z           !< Input: Constraint states at Time;
                                                                       !!   Output: Constraint states at Time + Interval
      TYPE(Morison_OtherStateType),       INTENT(INOUT) :: OtherState  !< Input: Other states at Time;
                                                                       !!   Output: Other states at Time + Interval
      TYPE(Morison_MiscVarType),          INTENT(INOUT) :: m           !< Misc/optimization variables            
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat     !< Error status of the operation     
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg      !< Error message if ErrStat /= ErrID_None

         ! Local variables
                  
      INTEGER(IntKi)                                    :: ErrStat2    ! Error status of the operation (occurs after initial error)
      CHARACTER(ErrMsgLen)                              :: ErrMsg2     ! Error message if ErrStat2 /= ErrID_None
                        
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
           
      

      
END SUBROUTINE Morison_UpdateStates

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE Morison_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )   
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(Morison_InputType),           INTENT(IN   )  :: u           !< Inputs at Time
      TYPE(Morison_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(Morison_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(Morison_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(Morison_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(Morison_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time
      TYPE(Morison_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at Time (Input only so that mesh con-
                                                                       !!   nectivity information does not have to be recalculated)
      TYPE(Morison_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables            
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      REAL(ReKi)                                        :: F_D(6), F_DP(6), D_F_I(3), kvec(3), v(3),  vf(3), vrel(3), vmag
      INTEGER                                           :: I, J, K, nodeIndx
      REAL(ReKi)                                        :: elementWaterState
      REAL(ReKi)                                        :: AllOuts(MaxMrsnOutputs)  ! TODO: think about adding to OtherState
      REAL(ReKi)                                        :: qdotdot(6) ,qdotdot2(3)     ! The structural acceleration of a mesh node
      !REAL(ReKi)                                        :: accel_fluid(6) ! Acceleration of fluid at the mesh node
      REAL(ReKi)                                        :: dragFactor     ! The lumped drag factor
      REAL(ReKi)                                        :: AnProd         ! Dot product of the directional area of the joint
      REAL(ReKi)                                        :: F_B(6)
      REAL(ReKi)                                        :: C(3,3)
      REAL(ReKi)                                        :: sgn
      REAL(ReKi)                                        :: D_AM_M(6,6)
      REAL(ReKi)                                        :: nodeInWater
      REAL(ReKi)                                        :: D_dragConst     ! The distributed drag factor
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Compute outputs here:
      
      ! We need to attach the distributed drag force (D_F_D), distributed inertial force (D_F_I), and distributed dynamic pressure force (D_F_DP) to the Misc type so that we don't need to
      ! allocate their data storage at each time step!  If we could make them static local variables (like in C) then we could avoid adding them to the OtherState datatype.  
      ! The same is true for the lumped drag (L_F_D) and the lumped dynamic pressure (L_F_DP)
         
      DO J = 1, y%DistribMesh%Nnodes
         
            ! Obtain the node index because WaveVel, WaveAcc, and WaveDynP are defined in the node indexing scheme, not the markers
         nodeIndx = p%distribToNodeIndx(J)
          
            ! Determine in or out of water status for the element which this node is a part of.        
            ! NOTE: This will find the closest WaveTime index (wvIndx) which is has waveTime(wvIndx) > = Time.  If WaveDT = DT then waveTime(wvIndx) will equal Time
            ! For WaveMod = 6 or WaveMod = 5 WaveDT must equal DT for the returned value of elementWaterState to be meaningful, for other WaveMod, 
            ! elementWaterState is the same for all time for a given node, J.
        elementWaterState = REAL( InterpWrappedStpInt( REAL(Time, SiKi), p%WaveTime(:), p%elementWaterState(:,J), m%LastIndWave, p%NStepWave + 1 ), ReKi )
       
         
         ! Determine the dynamic pressure at the marker
         m%D_FDynP(J) = InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%WaveDynP(:,nodeIndx), &
                                    m%LastIndWave, p%NStepWave + 1 )
         
            
         DO I=1,3
               ! Determine the fluid acceleration and velocity at the marker
            m%D_FA(I,J) = InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%WaveAcc(:,nodeIndx,I), &
                                    m%LastIndWave, p%NStepWave + 1       )
            m%D_FV(I,J) = InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%WaveVel(:,nodeIndx,I), &
                                    m%LastIndWave, p%NStepWave + 1       )
            
            vrel(I) =  m%D_FV(I,J) - u%DistribMesh%TranslationVel(I,J)
            
            m%D_F_I(I,J) = elementWaterState * InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%D_F_I(:,I,J), &
                                    m%LastIndWave, p%NStepWave + 1       )
         END DO
         
            ! (k x vrel x k)
         kvec =  p%Nodes(nodeIndx)%R_LToG(:,3)
         v = vrel - Dot_Product(kvec,vrel)*kvec
         vmag = sqrt( v(1)*v(1) + v(2)*v(2) + v(3)*v(3)  )
         
                   
            ! Distributed added mass loads
            ! need to multiply by elementInWater value to zero out loads when out of water 
         qdotdot2(1)    =       elementWaterState *u%DistribMesh%TranslationAcc(1,J)
         qdotdot2(2)    =       elementWaterState *u%DistribMesh%TranslationAcc(2,J)
         qdotdot2(3)    =       elementWaterState *u%DistribMesh%TranslationAcc(3,J)
            ! calculated the added mass forces (moments are zero)
         m%D_F_AM_M(1:3,J)  = -matmul( p%D_AM_M (:,:,J) , qdotdot2 )  !bjj: these lines take up a lot of time. are the matrices sparse?

         DO I=1,6
            IF (I < 4 ) THEN
                  ! We are now combining the dynamic pressure term into the inertia term  
               m%D_F_AM_MG(I,J) = -p%D_AM_MG(J)*u%DistribMesh%TranslationAcc(I,J)
               m%D_F_AM_F(:,J)  = -p%D_AM_F(J)*u%DistribMesh%TranslationAcc(I,J)
               m%D_F_AM(I,J)    = m%D_F_AM_M(I,J) + m%D_F_AM_MG(I,J) + m%D_F_AM_F(I,J)           
               m%D_F_D(I,J) = elementWaterState * vmag*v(I) * p%D_dragConst(J)      
               m%D_F_B(I,J) = elementWaterState * p%D_F_B(I,J)
               y%DistribMesh%Force(I,J) = m%D_F_AM(I,J) + m%D_F_D(I,J)  + m%D_F_I(I,J) + m%D_F_B(I,J) +  p%D_F_MG(I,J) + p%D_F_BF(I,J)
            ELSE
               m%D_F_B(I,J) = elementWaterState * p%D_F_B(I,J)
               y%DistribMesh%Moment(I-3,J) =   m%D_F_B(I,J) + p%D_F_BF(I,J)     
            END IF
         END DO  ! DO I
         
         
         
      ENDDO

      
      ! NOTE:  All wave kinematics have already been zeroed out above the SWL or instantaneous wave height (for WaveStMod > 0), so loads derived from the kinematics will be correct
      !        without the use of a nodeInWater value, but other loads need to be multiplied by nodeInWater to zero them out above the SWL or instantaneous wave height.
      
      DO J = 1, y%LumpedMesh%Nnodes
         
            ! Obtain the node index because WaveVel, WaveAcc, and WaveDynP are defined in the node indexing scheme, not the markers

         nodeIndx = p%lumpedToNodeIndx(J)
         nodeInWater = REAL( InterpWrappedStpInt( REAL(Time, SiKi), p%WaveTime(:), p%nodeInWater(:,nodeIndx), m%LastIndWave, p%NStepWave + 1 ), ReKi )
            ! Determine the dynamic pressure at the marker
         m%L_FDynP(J) = InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%WaveDynP(:,nodeIndx), &
                                    m%LastIndWave, p%NStepWave + 1       )
         
         
         DO I=1,3
               ! Determine the fluid acceleration and velocity at the marker
            m%L_FA(I,J) = InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%WaveAcc(:,nodeIndx,I), &
                                    m%LastIndWave, p%NStepWave + 1       )
               
            m%L_FV(I,J) = InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%WaveVel(:,nodeIndx,I), &
                                    m%LastIndWave, p%NStepWave + 1       )
            vrel(I)     = m%L_FV(I,J) - u%LumpedMesh%TranslationVel(I,J)
         END DO
         
         
        
            ! Compute the dot product of the relative velocity vector with the directional Area of the Joint
         vmag =  nodeInWater * ( vrel(1)*p%L_An(1,J) + vrel(2)*p%L_An(2,J) + vrel(3)*p%L_An(3,J) )
         AnProd = p%L_An(1,J)**2 + p%L_An(2,J)**2 + p%L_An(3,J)**2
         IF (EqualRealNos(AnProd, 0.0_ReKi)) THEN
            dragFactor = 0.0
         ELSE
            dragFactor = p%Nodes(nodeIndx)%JAxCd*p%WtrDens*abs(vmag)*vmag / ( 4.0_ReKi * AnProd )
         END IF
         
 
            ! Lumped added mass loads
         qdotdot                 = reshape((/u%LumpedMesh%TranslationAcc(:,J),u%LumpedMesh%RotationAcc(:,J)/),(/6/))   
         m%L_F_AM(:,J)           = matmul( p%L_AM_M(:,:,J) , ( - qdotdot) )
         DO I=1,3
            m%L_F_AM(I,J) = nodeInWater * m%L_F_AM(I,J)  ! Note that the rotational components are zero because L_AM_M is populated with only the upper-left 3x3
         END DO
         
         DO I=1,6
                        
            ! We are now combining the dynamic pressure term into the inertia term
            m%L_F_I(I,J) = InterpWrappedStpReal ( REAL(Time, SiKi), p%WaveTime(:), p%L_F_I(:,I,J), &
                                    m%LastIndWave, p%NStepWave + 1       ) 
            
            IF (I < 4 ) THEN
      
               m%L_F_D(I,J) =  p%L_An(I,J)*dragFactor
               m%L_F_B(I,J) =  nodeInWater*p%L_F_B(I,J)
               y%LumpedMesh%Force(I,J) = m%L_F_AM(I,J) + m%L_F_D(I,J) +  m%L_F_B(I,J) + m%L_F_I(I,J)  +  p%L_F_BF(I,J)
            ELSE
               m%L_F_B(I,J) =  nodeInWater*p%L_F_B(I,J)
               y%LumpedMesh%Moment(I-3,J) =   m%L_F_AM(I,J) + m%L_F_B(I,J) +   p%L_F_BF(I,J)
            END IF
            
            
         END DO      
      ENDDO
     
         ! OutSwtch determines whether or not to actually output results via the WriteOutput array
         ! 1 = Morison will generate an output file of its own.  2 = the caller will handle the outputs, but
         ! Morison needs to provide them.  3 = Both 1 and 2, 0 = No one needs the Morison outputs provided
         ! via the WriteOutput array.
         
      IF ( p%OutSwtch > 0 ) THEN
     
            ! Map calculated results into the AllOuts Array
         CALL MrsnOut_MapOutputs(Time, y, p, u, m, AllOuts, ErrStat, ErrMsg)
               
      
            ! Put the output data in the WriteOutput array
   
         DO I = 1,p%NumOuts

            y%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      
         END DO
         
         
            ! Generate output into the output file
            
         IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN
            CALL MrsnOut_WriteOutputs( p%UnOutFile, Time, y, p, ErrStat, ErrMsg )         
         END IF
      END IF
      
   
END SUBROUTINE Morison_CalcOutput


!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
SUBROUTINE Morison_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )  
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(Morison_InputType),           INTENT(IN   )  :: u           !< Inputs at Time                    
      TYPE(Morison_ParameterType),       INTENT(IN   )  :: p           !< Parameters                             
      TYPE(Morison_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(Morison_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(Morison_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(Morison_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time                   
      TYPE(Morison_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables            
      TYPE(Morison_ContinuousStateType), INTENT(  OUT)  :: dxdt        !< Continuous state derivatives at Time
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation     
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Compute the first time derivatives of the continuous states here:
      
      dxdt%DummyContState = 0.0
         

END SUBROUTINE Morison_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for updating discrete states
SUBROUTINE Morison_UpdateDiscState( Time, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )   
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        !< Current simulation time in seconds   
      TYPE(Morison_InputType),           INTENT(IN   )  :: u           !< Inputs at Time                       
      TYPE(Morison_ParameterType),       INTENT(IN   )  :: p           !< Parameters                                 
      TYPE(Morison_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(Morison_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Input: Discrete states at Time; 
                                                                       !!   Output: Discrete states at Time + Interval
      TYPE(Morison_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(Morison_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time          
      TYPE(Morison_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables            
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Update discrete states here:
      
      ! StateData%DiscState = 

END SUBROUTINE Morison_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
SUBROUTINE Morison_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )   
!..................................................................................................................................
   
      REAL(DbKi),                        INTENT(IN   )  :: Time        !< Current simulation time in seconds   
      TYPE(Morison_InputType),           INTENT(IN   )  :: u           !< Inputs at Time                       
      TYPE(Morison_ParameterType),       INTENT(IN   )  :: p           !< Parameters                           
      TYPE(Morison_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(Morison_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(Morison_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(Morison_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time       
      TYPE(Morison_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables            
      TYPE(Morison_ConstraintStateType), INTENT(  OUT)  :: z_residual  !< Residual of the constraint state equations using  
                                                                       !!     the input values described above      
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

               
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Solve for the constraint states here:
      
      z_residual%DummyConstrState = 0.0_ReKi

END SUBROUTINE Morison_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
   
END MODULE Morison
!**********************************************************************************************************************************
