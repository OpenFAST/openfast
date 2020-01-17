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
   ! This calculates the distributed buoyancy forces and moments on a given node

   REAL(ReKi),         INTENT ( IN    )  :: densWater
   REAL(ReKi),         INTENT ( IN    )  :: R                      ! Radius of node [m]
   REAL(ReKi),         INTENT ( IN    )  :: tMG                    ! Thickness of marine growth (adds to radius) [m]
   REAL(ReKi),         INTENT ( IN    )  :: dRdz                   ! Rate of change in radius with length at node [-]
   REAL(ReKi),         INTENT ( IN    )  :: Z                      ! z elevation of node [m] (not currently used)
   REAL(ReKi),         INTENT ( IN    )  :: C(3,3)
   REAL(ReKi),         INTENT ( IN    )  :: g
   REAL(ReKi),         INTENT (   OUT )  :: F_B(6)                 ! Distributed force and moment vector [N/m and N-m/m]
   
   REAL(DbKi)                           :: Reff,ReffSq,ReffCub,f1,f2,f3
   
   REAL(DbKi) :: CC(3,3)
   
   CC    = REAL(C,DbKi)
   Reff  = REAL(R + tMG,DbKi)                                      ! Effective radius after adding marine growth
  
   
   
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
          
SUBROUTINE DistrDragConst( densWater, Cd, R, tMG, DragConst  )   !@mhall: is there any reason to have this function?
                                                         ! It's only called once and it's a simple multiplication.

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
   ! This calculates the distributed hydrodynamic added mass matrix for a given node.

   REAL(ReKi),         INTENT ( IN    )  :: densWater
   REAL(ReKi),         INTENT ( IN    )  :: Ca                         ! Transverse added mass coefficient
   REAL(ReKi),         INTENT ( IN    )  :: AxCa                       ! Axial added mass coefficient (applied to tapered portions)
   REAL(ReKi),         INTENT ( IN    )  :: C(3,3)                     
   REAL(ReKi),         INTENT ( IN    )  :: R                          ! Radius at node [m]
   REAL(ReKi),         INTENT ( IN    )  :: tMG                        ! Thickness of marine growth (adds to radius) [m]
   REAL(ReKi),         INTENT ( IN    )  :: dRdZ                       ! Rate of change in radius with length at node [-]
   REAL(ReKi),         INTENT (   OUT )  :: AM_M(3,3)                  ! Distributed added mass matrix to be calculated [kg/m]
   
   REAL(ReKi)                            :: f,f2
   
   f         = Ca*densWater*Pi*(R+tMG)*(R+tMG)                         ! Transverse added mass scaler, applied to I - k*k^T
   f2        = AxCa*2.0*densWater*Pi*abs(dRdZ)*(R+tMG)*(R+tMG)         ! Axial added mass scaler, applied to k k^T 
   !AM_M      = 0.0
   AM_M(1,1) = f*(  C(2,3)*C(2,3) + C(3,3)*C(3,3) ) -f2*C(1,3)*C(1,3)  !<----@mhall: why is the f2 term being subtracted rather than added?
   AM_M(1,2) = f*( -C(1,3)*C(2,3)                 ) -f2*C(1,3)*C(2,3)
   AM_M(1,3) = f*( -C(1,3)*C(3,3)                 ) -f2*C(1,3)*C(3,3)
   
   AM_M(2,1) = f*( -C(1,3)*C(2,3)                 ) -f2*C(2,3)*C(1,3)
   AM_M(2,2) = f*(  C(1,3)*C(1,3) + C(3,3)*C(3,3) ) -f2*C(2,3)*C(2,3)  !<----@mhall: would it be cleaner to just use the k unit vector?  (also, diagonal terms can be shortened (1-k*kT))
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
   ! This calculates lumped buoyancy forces/moments on a member end.

   REAL(ReKi),         INTENT ( IN    )  :: sgn                  ! @mhall: this indicates if this is the start or end node so that direction is correct?
   REAL(ReKi),         INTENT ( IN    )  :: densWater
   REAL(ReKi),         INTENT ( IN    )  :: R                    ! Radius of end [m]
   REAL(ReKi),         INTENT ( IN    )  :: tMG                  ! Thickness of marine growth (adds to radius) [m]
   REAL(ReKi),         INTENT ( IN    )  :: Z                    ! z elevation of end [m] 

   REAL(ReKi),         INTENT ( IN    )  :: C(3,3)
   REAL(ReKi),         INTENT ( IN    )  :: g
   REAL(ReKi),         INTENT (   OUT )  :: F_B(6)               ! Lumped force and moment vector [N and N-m]


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
      
      ! Start with no splits.	  
      nSplits = 0  
      possibleSplits = -9999999.0  ! Initialize possibleSplit values to a number that never appears in the geometry.
      
      ! If the member is filled, add a possible split location at the global Z location of the filled free surface.
      IF ( members(I)%MmbrFilledIDIndx /= -1 ) THEN
         nSplits =  1
            ! The free surface is specified relative to the MSL.
         possibleSplits(1) = filledGroups(members(I)%MmbrFilledIDIndx)%FillFSLoc 
      END IF
      
      
      ! If the filled fluid hasn't already caused a split in the member at the still water line, then add a possible split there (at the water free surface, MSL2SWL).
      IF ( .NOT. IsThisSplitValueUsed(nSplits, possibleSplits, MSL2SWL )) THEN
            nSplits = nSplits + 1
            possibleSplits(nSplits) = MSL2SWL
      END IF  
      
      ! If there are one more depth-defined marine growth regions, add a possible split at each boundary if one doesn't already exist.        
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
      
      ! Add a possible split a the seabed if there isn't already one there.
       ! Check if seabed is equal to other possibleSplits      
      IF ( .NOT. IsThisSplitValueUsed(nSplits, possibleSplits, Zseabed) ) THEN
            nSplits = nSplits + 1
            possibleSplits(nSplits) = Zseabed
      END IF  
     
         
       ! Now determine which possible splits this member actually crosses and record them in the member's data structure.
       
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
   REAL(ReKi),               INTENT ( IN    )  :: D_F_B(:,:)           ! Distributed buoyancy force associated with the member
   REAL(ReKi),               INTENT ( IN    )  :: D_F_BF(:,:)          ! Distributed buoyancy force associated flooded/filled fluid within the member
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
SUBROUTINE SplitElementOnZBoundary( axis, boundary, numNodes, node1, node2, newNode, ErrStat, ErrMsg )
   !@mhall: This splits an element at a given global x, y, or z location?

   INTEGER,                  INTENT ( IN    )  :: axis                 ! which axis to work with in calculating positions along element (global x,y, or z)?
   REAL(ReKi),               INTENT ( IN    )  :: boundary             ! [axis] coordinate of boundary?
   INTEGER,                  INTENT ( INOUT )  :: numNodes
   TYPE(Morison_NodeType),   INTENT ( INOUT )  :: node1, node2  
   TYPE(Morison_NodeType),   INTENT (   OUT )  :: newNode
   INTEGER,                  INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),             INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None

   INTEGER                                     :: I, J
   REAL(ReKi)                                  :: s  ! interpolation factor
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   
      ! find normalized distance from 1nd node to the boundary
   CALL FindInterpFactor( boundary, node1%JointPos(axis), node2%JointPos(axis), s )
   newNode = node1 ! copy all base node properties
   
   newNode%JointPos(axis) =  boundary                    !@mhall: set [axis] coordinate of new node based on provided input [boundary]
   
   !@mthall: set other two coordinates of new node based on interpolation of original end node coordinates
   DO I=axis,axis+1
      J = MOD(I,3) + 1
      newNode%JointPos(J) =  node1%JointPos(J)*(1-s) + node2%JointPos(J)*s
   END DO
   
   newNode%R_LToG         =  node1%R_LToG   
   
END SUBROUTINE SplitElementOnZBoundary




SUBROUTINE DiscretizeMember( numNodes, nodes, member, propSet1, propSet2, numMGDepths, MGDepths, ErrStat, ErrMsg )
   ! This subdivides a Morison member according to its maximum desired 
   ! element length (MDivSize), allocating the member's arrays, and
   ! adding resuling new nodes to the master node array.
   
   INTEGER,                  INTENT ( INOUT )  :: numNodes             ! the current number of nodes being used in the model
   TYPE(Morison_NodeType),   INTENT ( INOUT )  :: nodes(:)             ! the array of nodes (maximum size for now)
   TYPE(Morison_MemberType), INTENT ( INOUT )  :: member               ! the morison member to discretize
   TYPE(Morison_MemberPropType) INTENT( IN  )  :: propSet1             ! property set of node 1
   TYPE(Morison_MemberPropType) INTENT( IN  )  :: propSet2             ! property set of node N+1
   INTEGER,                     INTENT( IN  )  :: numMGDepths
   TYPE(Morison_MGDepthsType),  INTENT( IN  )  :: MGDepths(:)
   INTEGER,                  INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),             INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
   TYPE(Morison_NodeType)                      :: node1, node2, newNode
   !TYPE(Morison_MemberType)                    :: element, newElement
   INTEGER                                     :: numDiv
   REAL(ReKi)                                  :: divSize(3)
   INTEGER                                     :: I, J, K
   INTEGER                                     :: origNumElements
   INTEGER                                     :: node1Indx, node2Indx, axis
   REAL(ReKi)                                  :: start, Loc
   REAL(ReKi)                                  :: s  ! interpolation factor
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = "" 
   
     
      ! some copying for convenience
      node1Indx  = member%Node1Indx
      CALL Morison_CopyNodeType( nodes(node1Indx), node1, MESH_NEWCOPY, ErrStat, ErrMsg )   !   node1 = nodes(node1Indx); bjj: note that not all fields have been initialized, but maybe that's okay.

      node2Indx   = member%Node2Indx          ! We need this index for the last sub-element
      CALL Morison_CopyNodeType( nodes(node2Indx), node2, MESH_NEWCOPY, ErrStat, ErrMsg )   !   node2 = nodes(node2Indx); bjj: note that not all fields have been initialized, but maybe that's okay.
      
      
      ! calculate and save member length
      CALL GetDistance(node1%JointPos, node2%JointPos, member%Len)                           ! Calculate member length.
      
      
      ! If the requested division size is less then the member length, create nodes along the member         
      IF ( member%MDivSize < member%Len ) THEN
	  
         ! Ensure a safe choice of x/y/z axis to use for splitting.
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
         numDiv = CEILING( memLen / member%MDivSize  )
         
         ! set number of elements in member and element size
         member%NElements = numDiv
         member%dl        = member%Len/numDiv
         
         ! allocate member arrays
         CALL AllocateMemberVariables(member, ErrStat2, ErrMsg2)
         
         ! now that the arrays are allocated, fill in values for the ends of the member
         member%R(  1)   = propSet1%PropD / 2.0            
         member%RMG(1)   = propSet1%PropD / 2.0 + node1%tMG 
         member%Rin(1)   = propSet1%PropD / 2.0 - propSet1%PropThck  
         member%R(  N+1) = propSet2%PropD / 2.0            
         member%RMG(N+1) = propSet2%PropD / 2.0 + node2%tMG 
         member%Rin(N+1) = propSet2%PropD / 2.0 - propSet2%PropThck 
         
         member%NodeIndx(  1) = node1Indx
         member%NodeIndx(N+1) = node2Indx
         
      
         ! get the division size along each axis
         DO K=1,3
            divSize(K) = (node2%JointPos(K) - node1%JointPos(K)) / numDiv
         END DO
         
      
         ! loop through the new node locations along the member and add a node at each
         DO J=1,numDiv - 1
            
            loc = start + divSize(axis)*J
            
            ! below was previously CALL SplitElementOnZBoundary( axis, loc, numNodes, node1, node2, newNode, ErrStat, ErrMsg )
            
            ! find normalized distance from 1nd node to the boundary
            CALL FindInterpFactor( boundary, node1%JointPos(axis), node2%JointPos(axis), s )
            newNode = node1 ! copy all base node properties
            
            newNode%JointPos(axis) =  boundary                    !@mhall: set [axis] coordinate of new node based on provided input [boundary]
            
            ! Set properties of new node based on interpolation of original end node coordinates
            
            !@mthall: set other two coordinates of new node based on interpolation of original end node coordinates
            DO I=axis,axis+1
               K = MOD(I,3) + 1
               newNode%JointPos(K) =  node1%JointPos(K)*(1-s) + node2%JointPos(K)*s
            END DO
            
            newNode%R_LToG         =  node1%R_LToG   
            
            newNode%NodeType = 2 ! interior node
            newNode%JointIndx = -1
            
            ! Set the marine growth thickness and density information for the new node
            CALL SetNodeMG( numMGDepths, MGDepths, newNode )
            
            ! Copy new node into the nodes array (this now just appends the new node to the list, the variables node1 and node2 don't need to change)
            numNodes               = numNodes + 1
            nodes(numNodes)        = newNode               ! this is copying all the fields in the type; type contains no allocatable arrays or pointers
            
            ! Now fill in properties of the member at the new node location
            I = J+1 ! the node index along the member, from 1 to N+1
            member%R(  I) =  member%R(  1)*(1-s) + member%R(  N+1)*s
            member%Rin(I) =  member%Rin(1)*(1-s) + member%Rin(N+1)*s
            member%RMG(I) =  member%R(I) + newNode%tMG
            
            ! record the node index in the member's node index list
            member%NodeIndx(I) = node1Indx + J
            
         END DO
         

      ELSE  ! if the member will only be one element long, record that and allocate its arrays
            
         ! set number of elements in member and element size
         member%NElements = 1
         member%dl        = member%Len
         
         ! allocate member arrays
         CALL AllocateMemberVariables(member, ErrStat2, ErrMsg2)
         
         ! now that the arrays are allocated, fill in values for the ends of the member
         member%R(  1)   = propSet1%PropD / 2.0            
         member%RMG(1)   = propSet1%PropD / 2.0 + node1%tMG 
         member%Rin(1)   = propSet1%PropD / 2.0 - propSet1%PropThck  
         member%R(  N+1) = propSet2%PropD / 2.0            
         member%RMG(N+1) = propSet2%PropD / 2.0 + node2%tMG 
         member%Rin(N+1) = propSet2%PropD / 2.0 - propSet2%PropThck 
         
         member%NodeIndx(  1) = node1Indx
         member%NodeIndx(N+1) = node2Indx
      
      END IF        
         
END SUBROUTINE DiscretizeMember         
      

SUBROUTINE AllocateMemberVariables(member, ErrStat, ErrMsg)
  
   TYPE(Morison_MemberType), INTENT ( INOUT )  :: member               ! the morison member to discretize
   INTEGER,                  INTENT (   OUT )  :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),             INTENT (   OUT )  :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
   INTEGER                 :: N
   
   
   N = member%NElements
   
   !@mhall: I'm hoping the arrays for each member can be allocated like this
   
   ALLOCATE ( member%NodeIndx(N+1), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the member node index array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF 
      
   ALLOCATE ( member%R(N+1), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the member node radius array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF 
   
   !@mhall: etc.
   
   
END SUBROUTINE AllocateMemberVariables


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
   
   if (NCoefDpth == 0) return
   
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


SUBROUTINE SetNodeMG( numMGDepths, MGDepths, node )
   ! sets the margine growth thickness of a single node (previously all nodes)
   INTEGER,                      INTENT( IN    )  :: numMGDepths
   TYPE(Morison_MGDepthsType),   INTENT( IN    )  :: MGDepths(:)
   TYPE(Morison_NodeType),       INTENT( INOUT )  :: node

   INTEGER                                     :: I, J
   REAL(ReKi)              :: z
   INTEGER                 :: indx1, indx2
   REAL(ReKi)              :: dd, s
   LOGICAL                 :: foundLess = .FALSE.
   
      
         !Find the table entry(ies) which match the node's depth value
      ! The assumption here is that the depth table is stored from largest
      ! to smallest in depth
      z = node%JointPos(3)
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
         node%tMG       = 0.0
         node%MGdensity = 0.0
      ELSE
         ! Linearly interpolate the coef values based on depth
         !CALL FindInterpFactor( z, CoefDpths(indx1)%Dpth, CoefDpths(indx2)%Dpth, s )
      
         dd = MGDepths(indx1)%MGDpth - MGDepths(indx2)%MGDpth
         IF ( EqualRealNos(dd, 0.0_ReKi) ) THEN
            s = 0.0_ReKi
         ELSE
            s = ( MGDepths(indx1)%MGDpth - z ) / dd
         END IF
         node%tMG       = MGDepths(indx1)%MGThck*(1-s) + MGDepths(indx2)%MGThck*s
         node%MGdensity = MGDepths(indx1)%MGDens*(1-s) + MGDepths(indx2)%MGDens*s
      END IF
      

END SUBROUTINE SetNodeMG


  

!====================================================================================================
SUBROUTINE Morison_ProcessMorisonGeometry( InitInp, ErrStat, ErrMsg )
!     This public subroutine process the input geometry and parameters and eliminates joint overlaps,  
!     sub-divides members, sets joint-level properties, etc. It is called by HydroDyn before Morison_Init
!     so that the node positions are available beforehand for other parts of HydroDyn.
!----------------------------------------------------------------------------------------------------  

      ! Passed variables
   
   TYPE(Morison_InitInputType),   INTENT( INOUT )   :: InitInp              ! the Morison initialization data 
   !TYPE(Morison_ParameterType),   INTENT( INOUT )   :: p                    ! tge Morison parameter data
   INTEGER,                       INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                  INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None

  
      ! Local variables
         
   INTEGER                                      :: I, J  !, j1, j2, tempINT                ! generic integer for counting
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
   INTEGER                                      :: numSplitNodes
   INTEGER                                      :: Node1Indx, Node2Indx   ! indices of joint nodes at member ends
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
         
         ! Set the marine growth thickness and density information for each joint node
         CALL SetNodeMG( InitInp%NMGDepths, InitInp%MGDepths, InitInp%Nodes(I) )
         
      END DO
      
      
      
      
          ! Allocate memory for Members arrays
          
      InitInp%NElements = InitInp%NMembers  
      
      ALLOCATE ( InitInp%Members(maxElements), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for Members array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
          
      
      ! loop through members, assign each its nodes, etc.
      DO I = 1,InitInp%NMembers  
         
         InitInp%Members(I)%Node1Indx = InitInp%InpMembers(I)%MJointID1Indx              ! Index of  the first node in the Morison_NodeType array
         InitInp%Members(I)%Node2Indx = InitInp%InpMembers(I)%MJointID2Indx              ! Index of  the second node in the Morison_NodeType array
         node1Indx                     = InitInp%Members(I)%Node1Indx
         node2Indx                     = InitInp%Members(I)%Node2Indx
         prop1Indx = InitInp%InpMembers(I)%MPropSetID1Indx
         prop2Indx = InitInp%InpMembers(I)%MPropSetID2Indx
         
            ! Make sure that Node1 has the lower Z value, re-order if necessary
            ! We need to do this because the local element coordinate system is defined such that the first node is located with a smaller global Z value
            ! than the second node.
            ! The local element coordinate system requires that Z1 <= Z2, and if Z1=Z2 then X1 <= X2, and if Z1=Z2, X1=X2 then Y1<Y2
   
         InitInp%Members(I)%InpMbrDist1         = 0.0
         InitInp%Members(I)%InpMbrDist2         = 1.0
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
            
            InitInp%Members(I)%Node1Indx = InitInp%InpMembers(I)%MJointID2Indx              
            InitInp%Members(I)%Node2Indx = InitInp%InpMembers(I)%MJointID1Indx  
            node1Indx                     = InitInp%Members(I)%Node1Indx
            node2Indx                     = InitInp%Members(I)%Node2Indx
            temp = prop1Indx
            prop1Indx = prop2Indx
            prop2Indx = temp
            InitInp%Members(I)%InpMbrDist1         = 1.0
            InitInp%Members(I)%InpMbrDist2         = 0.0
            
            ! --- Swap member coeffs if needed. 
            ! Fine in this loop since there is a unique CoefMember per Member (otherwise we could swap them several times).
            J = InitInp%InpMembers(I)%MmbrCoefIDIndx ! Index in CoefMembers table
            IF (J>0) THEN 
                ! NOTE: SWAP defined at the end of the current subroutine
                CALL SWAP(InitInp%CoefMembers(J)%MemberCd1    , InitInp%CoefMembers(J)%MemberCd2)
                CALL SWAP(InitInp%CoefMembers(J)%MemberCa1    , InitInp%CoefMembers(J)%MemberCa2)
                CALL SWAP(InitInp%CoefMembers(J)%MemberCp1    , InitInp%CoefMembers(J)%MemberCp2)
                CALL SWAP(InitInp%CoefMembers(J)%MemberAxCa1  , InitInp%CoefMembers(J)%MemberAxCa2)
                CALL SWAP(InitInp%CoefMembers(J)%MemberAxCp1  , InitInp%CoefMembers(J)%MemberAxCp2)
                CALL SWAP(InitInp%CoefMembers(J)%MemberCdMG1  , InitInp%CoefMembers(J)%MemberCdMG2)
                CALL SWAP(InitInp%CoefMembers(J)%MemberCaMG1  , InitInp%CoefMembers(J)%MemberCaMG2)
                CALL SWAP(InitInp%CoefMembers(J)%MemberCpMG1  , InitInp%CoefMembers(J)%MemberCpMG2)
                CALL SWAP(InitInp%CoefMembers(J)%MemberAxCaMG1, InitInp%CoefMembers(J)%MemberAxCaMG2)
                CALL SWAP(InitInp%CoefMembers(J)%MemberAxCpMG1, InitInp%CoefMembers(J)%MemberAxCpMG2)
            END IF 
            
         END IF
         
         
         ! assign member to its joints' lists of connected members
         Node1Indx = InitInp%Members(I)%Node1Indx
         Node2Indx = InitInp%Members(I)%Node2Indx
         InitInp%Nodes(Node1Indx)%Nconnections = InitInp%Nodes(Node1Indx)%Nconnections + 1   ! increment the joint's number of connections
         InitInp%Nodes(Node2Indx)%Nconnections = InitInp%Nodes(Node2Indx)%Nconnections + 1   ! increment the joint's number of connections
         InitInp%Nodes(Node1Indx)%ConnectionList(InitInp%Nodes(Node1Indx)%Nconnections) = I  ! assign the member ID to the joint's connection list (positive since end 1)
         InitInp%Nodes(Node2Indx)%ConnectionList(InitInp%Nodes(Node2Indx)%Nconnections) =-I  ! assign the member ID to the joint's connection list (negative since end 2)
         
         ! assign member some input quantities
         InitInp%Members(I)%InpMbrIndx       = I
         InitInp%Members(I)%MDivSize         = InitInp%InpMembers(I)%MDivSize
         InitInp%Members(I)%MCoefMod         = InitInp%InpMembers(I)%MCoefMod
         InitInp%Members(I)%MmbrCoefIDIndx   = InitInp%InpMembers(I)%MmbrCoefIDIndx
         InitInp%Members(I)%MmbrFilledIDIndx = InitInp%InpMembers(I)%MmbrFilledIDIndx
      
         ! assign member length
         CALL GetDistance( InitInp%Nodes(node1Indx)%JointPos, InitInp%Nodes(node2Indx)%JointPos, d)
         InitInp%Members(I)%InpMbrLen           = d
         
         ! Calculate the element-level direction cosine matrix and attach it to the entry in the elements array
         CALL Morison_DirCosMtrx( InitInp%Nodes(node1Indx)%JointPos, InitInp%Nodes(node2Indx)%JointPos, InitInp%Members(I)%R_LToG )
        
         InitInp%Members(I)%PropPot  =  InitInp%InpMembers(I)%PropPot                  ! Flag specifying whether member is modelled in WAMIT [true = modelled in WAMIT, false = not modelled in WAMIT]
         
         ! InitInp%Nodes(node1Indx)%R_LToG = InitInp%Members(I)%R_LToG
        ! InitInp%Nodes(node2Indx)%R_LToG = InitInp%Members(I)%R_LToG
        
        !@mhall: no longer any need to split or subdivide members. 
        ! Instead, we need to discretize a member and create node points along it.
        
        DO I = 1, InitInp%NMembers
            ! discretize a member and create node points along it. This will be done in Morison_Init
            CALL DiscretizeMember( InitInp%NNodes, InitInp%Nodes, InitInp%Members(I), InitInp%MPropSets(prop1Indx), InitInp%MPropSets(prop2Indx), InitInp%NMGDepths, InitInp%MGDepths, ErrStat, ErrMsg )  
           !@mhall: hoping passing one entry of Members is okay - otherwise can pass all and put for loop inside DiscretizeMember
         
        END DO 
         
      
         ! Set the fill properties onto the elements
      !@mthall: This is now done in member setup loop in Morison_Init.  CALL SetElementFillProps( InitInp%NFillGroups, InitInp%FilledGroups, InitInp%NElements, InitInp%Members )
     
     
      !@mhall: The below should be changed to operate on each member m%Members(I) and the arrays within it, or each joint.
      !        The hydro coefficients should be set for each node along a member.

         ! Set the element Cd, Ca, and Cp coefs
      CALL SetElementCoefs( InitInp%SimplCd, InitInp%SimplCdMG, InitInp%SimplCa, InitInp%SimplCaMG, InitInp%SimplCp, InitInp%SimplCpMG, InitInp%SimplAxCa, InitInp%SimplAxCaMG, InitInp%SimplAxCp, InitInp%SimplAxCpMG,InitInp%CoefMembers, InitInp%NCoefDpth, InitInp%CoefDpths, InitInp%NNodes, InitInp%Nodes, InitInp%NElements, InitInp%Members )   
      
         ! Set the axial coefs AxCd and AxCa
     CALL SetAxialCoefs( InitInp%NJoints, InitInp%NAxCoefs, InitInp%AxialCoefs, InitInp%NNodes, InitInp%Nodes, InitInp%NElements, InitInp%Members )      
      
      
      !@mhall: duplicate nodes at member ends no longer needed
      
      
         
         ! 6) Store information necessary to compute the user-requested member outputs and joint outputs.  The requested output locations
         !    may be located in between two simulation nodes, so quantities will need to be interpolated. qOutput = q1*s + q2*(1-s), where 0<= s <= 1.
         
         ! NOTE: since we need to mantain the input geometry, the altered members are now part of the simulation mesh and 
         !       we will generate a mapping between the input and simulation meshes which is needed to generate user-requested outputs.
   
        
        
      END DO  !I = 1,InitInp%NMembers  
      
    
         
   ELSE  
      
      
         ! No Morison elements, so no processing is necessary, but set nodes and elements to 0.
         
     ! p%NMorisonNodes    = 0  
    !  p%NMorisonElements = 0
      
   END IF
      CONTAINS

        SUBROUTINE SWAP(x1,x2)
           Real(Reki),intent(inout) :: x1,x2
           Real(Reki) :: tmp
           tmp = x1
           x1  = x2
           x2  = tmp
        END SUBROUTINE SWAP
END SUBROUTINE Morison_ProcessMorisonGeometry

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps. 
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
!! A lot of the model setup has been done previously in Morison_ProcessMorisonGeometry, and stored in InitInp.
SUBROUTINE Morison_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(Morison_InitInputType),       INTENT(INOUT)  :: InitInp     !< Input data for initialization routine !intent out because of MOVE_ALLOC
      TYPE(Morison_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
      TYPE(Morison_ParameterType),       INTENT(  OUT)  :: p           !< Parameters      
      TYPE(Morison_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
      TYPE(Morison_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
      TYPE(Morison_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
      TYPE(Morison_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states (this contains the Members array) 
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

      TYPE(Morison_MemberType) :: mem    ! the current member
      INTEGER                  :: N
      REAL(ReKi)               :: dl
      REAL(ReKi)               :: vec(3)
      REAL(ReKi)               :: phi    ! member tilt angle
      REAL(ReKi)               :: beta   ! member tilt heading
      REAL(ReKi)               :: cosPhi
      REAL(ReKi)               :: sinPhi
      REAL(ReKi)               :: tanPhi
      REAL(ReKi)               :: sinBeta
      REAL(ReKi)               :: cosBeta
      REAL(ReKi)               :: Za
      REAL(ReKi)               :: Zb
      
      
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
      
     
     
      ! ----------------------- set up the members -----------------------
      
      ! allocate and copy in the Members array
      
      ALLOCATE ( m%Members(p%NMembers), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for the members array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF     
          
      m%Members = InitInp%Members
      
     
     ! ----------------------- set up the nodes -----------------------
      
      ! allocate and copy in the nodes list
       
      p%NNodes   = InitInp%NNodes
      
      ALLOCATE ( p%Nodes(p%NNodes), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for Nodes array.'
         ErrStat = ErrID_Fatal
         RETURN
      END IF
      
      p%Nodes    = InitInp%Nodes   
      
      
      ! allocate and copy in hydrodynamic arrays for the nodes
      
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
      
      
  ! allocate and initialize some wave-related arrays   
   
   ALLOCATE ( elementWaterStateArr( 0:NStepWave, p%NNodes ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the elementWaterStateArr array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   elementWaterStateArr = 0 ! out of the water
   
   
   
   ALLOCATE ( m%F_I( 0:NStepWave, 6, NNodes ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the inertial forces/moments array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%F_I = 0.0
   
   ALLOCATE ( m%F_DP( 0:NStepWave, 6, NNodes ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the dynamic pressure forces/moments array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%F_DP = 0.0
   
     
   
   ! allocate and initialize joint-specific arrays   
      
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
  
   ALLOCATE ( L_An( 3, numLumpedMarkers ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the joint directional area array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   L_An = 0.0
   
   
   
   ! initialize load arrays for all nodes
   
   CALL AllocateNodeLoadVariables(m, p%NNodes, ErrStat, ErrMsg)
   
   ! a few additional loads that
   
  
   
   
   ! Create the input and output meshes associated with loads at the nodes
      
   CALL MeshCreate( BlankMesh      = u%Mesh          &
                     ,IOS          = COMPONENT_INPUT        &
                     ,Nnodes       = p%NNodes      &
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
   
   
   DO I=1,p%NNodes
         
      IF ( p%Nodes(I)%NodeType == 3 ) THEN
         
      END IF
      
         ! Create the node on the mesh
      
      !orientation = transpose(p%Nodes(I)%R_LToG )
      
      CALL MeshPositionNode (u%Mesh          &
                        , count                    &
                        , p%Nodes(I)%JointPos      &  ! this info comes from FAST
                        , ErrStat                  &
                        , ErrMsg                   &
                        ) !, transpose(p%Nodes(I)%R_LToG)          )
      IF ( ErrStat /= 0 ) RETURN
   
   
      Morison_Rad(count) = p%Nodes(I)%R   ! set this for FAST visualization
      
      !@mhall: what is happening in these lines?
      distribToNodeIndx(count) = I
      nodeToDistribIndx(I) = count
   
         ! Create the mesh element
   
      CALL MeshConstructElement (u%Mesh   &
                            , ELEMENT_POINT      &                                  
                            , ErrStat            &
                            , ErrMsg  &
                            , count                  &
                                        )
      count = count + 1    
   
   END DO

   
   
   CALL MeshCommit ( u%Mesh   &
                      , ErrStat            &
                      , ErrMsg             )
   
   IF ( ErrStat /= 0 ) THEN
         RETURN
   END IF 
   
      ! Initialize the inputs
   DO I=1,u%Mesh%Nnodes
      u%Mesh%Orientation(:,:,I) = u%Mesh%RefOrientation(:,:,I)
   END DO
   
   u%Mesh%TranslationDisp = 0.0
   u%Mesh%TranslationVel  = 0.0
   u%Mesh%RotationVel     = 0.0
   u%Mesh%TranslationAcc  = 0.0
   u%Mesh%RotationAcc     = 0.0
   
   
   ! Duplicate the input mesh to create the output mesh
   
   CALL MeshCopy (    SrcMesh      = u%Mesh &
                     ,DestMesh     = y%Mesh         &
                     ,CtrlCode     = MESH_SIBLING           &
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat                &
                     ,ErrMess      = ErrMsg                 &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 )

   u%Mesh%RemapFlag = .TRUE.
   y%Mesh%RemapFlag = .TRUE.
   
      
      
     !--------------- 
      
      ! Define parameters here:
       
     
      p%DT  = Interval


         ! Define initial system states here:

      x%DummyContState           = 0
      xd%DummyDiscState          = 0
      z%DummyConstrState         = 0
      OtherState%DummyOtherState = 0
      m%LastIndWave              = 1

   ! IF ( p%OutSwtch > 0 ) THEN  @mhall: I think the below need to be allocated in all cases



      
      
      
 ! loop through elements of each member and precalculate the required quantities
 DO im = 1, p%NMembers
 
   
   mem = m%Members(im)

   N = mem%NElements
   dl = mem%dl

   ! find fill location of member (previously in SetElementFillProps)
   IF ( mem%MmbrFilledIDIndx > 0 ) THEN
      
      mem%FillDens     =  InitInp%FilledGroups(elements(I)%MmbrFilledIDIndx)%FillDens
      mem%FillFSLoc    =  InitInp%FilledGroups(elements(I)%MmbrFilledIDIndx)%FillFSLoc
   ELSE
      mem%FillDens     =   0.0
      mem%FillFSLoc    =   0.0
   END IF


   ! calculate instantaneous incline angle and heading, and related trig values
   vec = p%Nodes(mem%NodeIndx(N+1))%JointPos - p%Nodes(mem%NodeIndx(1))%JointPos 
   
   phi = arccos(vec(3)/SQRT(Dot_Product(vec,vec)))  ! incline angle   
!   beta = arctan2(vec(2), vec(1))                   ! heading of incline
   
   sinPhi = sin(phi)
   cosPhi = cos(phi)  
!   tanPhi = tan(phi)
   
!   sinBeta = sin(beta)
!   cosBeta = cos(beta) 


 
   ! calculate l_fill or z_overfill for the member
   mem%l_fill = (mem%FillFSLoc - )/cosPhi   ! fill length along cylinder axis

   if mem%l_fill > mem%Len then
      mem%z_overfill = mem%FillFSLoc - Zb
   else 
      mem%z_overfill = 0.0
   end if
         

   ! calculate h_floor if seabed-piercing
   if (Za < -WtrDepth) then
      do i=2,N+1
         if (Za > -WtrDepth) then            ! find the lowest node above the seabed
            mem%h_floor = (-WtrDepth-Za)/cosPhi  ! get the distance from the node to the seabed along the member axis (negative value)
            mem%i_floor = i-1                    ! record the number of the element that pierces the seabed
            break
         end if
      end do
   end if


   ! calculate element-level values
   DO i = 1,N
      mem%m_mg(i) = (mem%RMG(i+1) - mem%RMG(i))/dl
      mem%m_in(i) = (mem%Rin(i+1) - mem%Rin(i))/dl
      
      mem%alpha(   i) = GetAlpha(mem%RMG(i), mem%RMG(i+1))
      mem%alpha_fb(i) = GetAlpha(mem%Rin(i), mem%Rin(i+1))
      
   END DO


   ! force-related constants for each element
   DO i = 1,N

      Za = p%Nodes(mem%NodeIndx(  i))%JointPos(3)   ! z location of node i
      Zb = p%Nodes(mem%NodeIndx(i+1))%JointPos(3)   ! z location of node i+1

      ! ------------------ marine growth weight and inertia, and flooded ballast inertia--------------------
         
      if (i > mem%i_floor) then         
         ! full marine growth: get the properties for each half-element lumped to the appropriate node
                  
         Rmid   = 0.5*(mem%R(  i)+mem%R(  i+1))  ! radius at middle of segment, where division occurs
         RmidMG = 0.5*(mem%RMG(i)+mem%RMG(i+1))  ! radius with marine growth at middle of segment, where division occurs
         Rmidin = 0.5*(mem%Rin(i)+mem%Rin(i+1))  ! radius of member interior at middle of segment, where division occurs
         Lmid   = 0.5*dl   ! = 0.5*(R2-R1)/m  half-length of segment

         CALL MarineGrowthPartSegment(mem%R(i  ), Rmid, mem%RMG(i  ),RmidMG, Lmid, rhoMG,  mem%m_mg_l(i), mem%h_cmg_l(i), mem%I_lmg_l(i), mem%I_rmg_l(i))   ! get precomupted quantities for lower half-segment
         CALL MarineGrowthPartSegment(mem%R(i+1), Rmid, mem%RMG(i+1),RmidMG,-Lmid, rhoMG,  mem%m_mg_u(i), mem%h_cmg_u(i), mem%I_lmg_u(i), mem%I_rmg_u(i))   ! get precomupted quantities for upper half-segment
         CALL FloodedBallastPartSegment(mem%Rin(i  ), mem%Rmidin,  Lmid, FillDens,  mem%m_fb_l(i), mem%h_cfb_l(i), mem%I_lfb_l(i), mem%I_rfb_l(i))   ! get precomupted quantities for lower half-segment
         CALL FloodedBallastPartSegment(mem%Rin(i+1), mem%Rmidin, -Lmid, FillDens,  mem%m_fb_u(i), mem%h_cfb_u(i), mem%I_lfb_u(i), mem%I_rfb_u(i))   ! get precomupted quantities for upper half-segment

      else if (i == mem%ifloor) then         
         ! crossing seabed: get the properties for part-element above the seabed and lump to the upper node      

         Rmid   = (-mem%h_floor*mem%R(  i) +(dl+mem%h_floor)*mem%R(  i+1))/dl
         RmidMG = (-mem%h_floor*mem%RMG(i) +(dl+mem%h_floor)*mem%RMG(i+1))/dl
         Rmidin = (-mem%h_floor*mem%Rin(i) +(dl+mem%h_floor)*mem%Rin(i+1))/dl
         Lmid   = -mem%h_floor
         
         CALL MarineGrowthPartSegment(mem%R(i+1), Rmid, mem%RMG(i+1),RmidMG, -Lmid, rhoMG,  mem%m_mg_u(i), mem%h_cmg_u(i), mem%I_lmg_u(i), mem%I_rmg_u(i))   ! get precomupted quantities for upper half-segment
         CALL FloodedBallastPartSegment(mem%Rin(i+1), Rmidin, -Lmid, FillDens,  mem%m_fb_u(i), mem%h_cfb_u(i), mem%I_lfb_u(i), mem%I_rfb_u(i))   ! get precomupted quantities for upper half-segment
      
      end if
         
      
      ! ------------------ flooded ballast weight (done) --------------------

      ! fully buried element
      if (Zb < -WtrDepth) then
         mem%floodstatus(i) = 0
      
      ! fully filled elements 
      if (mem%FillFSLoc > Zb) then  
         mem%floodstatus(i) = 1
      
         ! depth-adjusted force distribution constant
         mem%alpha_fb_star(i) = mem%alpha_fb(i)*( Zb - mem%FillFSLoc )^3 / ( ( (1-mem%alpha_fb(i))*(Za - mem%FillFSLoc))^3 + mem%alpha_fb(i)*(Zb - mem%FillFSLoc)^3 )
         
         ! force and moment magnitude constants
         mem%Cfl_fb(i) = TwoPi * m * mem%FillDens * g * dl *( (li - mem%lfill)*Rin(i) + 0.5*((li - mem%lfill)*m_in + mem%Rin(i))*dl + 1/3*m_in*dl^2 )
         mem%Cfr_fb(i) =    Pi *     mem%FillDens * g * dl *( mem%Rin(i)^2 + m_in_in*Rin(i)*dl +1/3 m_in^2 *dl^2 )
         mem%CM0_fb(i) = TwoPi *     mem%FillDens * g * dl *( 0.25*dl^3*m_in^4 + 0.25*dl^3*m_in^2 + dl^2*m_in^3*Rin(i) + 2/3*dl^2*m_in*mem%Rin(i) + 1.5*dl*m_in^2*mem%Rin(i)^2 + 0.5*dl*mem%Rin(i)^2 + m_in*mem%Rin(i)^3 )
         
         
      ! partially filled element
      else if ((mem%FillFSLoc > Za) .AND. (mem%FillFSLoc < Zb)) then
         mem%floodstatus(i) = 2
      
         ! length along axis from node i to fill level
         mem%h_fill = mem%l_fill - (i-1)*dl
      
         ! depth-adjusted force distribution constant
         mem%alpha_fb_star(i) = (1 - mem%alpha_fb(i))*( Za - mem%FillFSLoc )^3 / ( ( (1-mem%alpha_fb(i))*(Za - mem%FillFSLoc))^3 - mem%alpha_fb(i)*(Zb - mem%FillFSLoc)^3 )
         
         ! force and moment magnitude constants
         mem%Cfl_fb(i) = TwoPi * m * mem%FillDens * g * mem%h_fill *( (li - mem%lfill)*mem%Rin(i) + 0.5*((li - mem%lfill)*m_in + mem%Rin(i))*mem%h_fill + 1/3*m_in*mem%h_fill^2 )
         mem%Cfr_fb(i) =    Pi *     mem%FillDens * g * mem%h_fill *( mem%Rin(i)^2 + m_in_in*mem%Rin(i)*mem%h_fill +1/3 m_in^2 *mem%h_fill^2 )
         mem%CM0_fb(i) = TwoPi *     mem%FillDens * g * mem%h_fill *( 0.25*mem%h_fill^3*m_in^4 + 0.25*mem%h_fill^3*m_in^2 + mem%h_fill^2*m_in^3*mem%Rin(i) + 2/3*mem%h_fill^2*m_in*mem%Rin(i) + 1.5*mem%h_fill*m_in^2*mem%Rin(i)^2 + 0.5*mem%h_fill*mem%Rin(i)^2 + m_in*mem%Rin(i)^3 )
      
      ! unflooded element
      else
         mem%floodstatus(i) = 0
      
      end if
      

   end do ! end looping through elements   
  
end do ! looping through members  
      

      
! loop through joints to calculate joint quantities (the joints are the first NJoints nodes)


   ! CA is the added mass coefficient for three dimensional bodies in infinite fluid (far from boundaries) The default value is 2/Pi
   AMfactor = 2.0 * densWater * Pi / 3.0
     
   usedJointList = .FALSE.   
   commonNodeLst = -1

DO I = 1,p%NJoints      

   An        = 0.0
   Vn        = 0.0
   I_n       = 0.0
   
   IF ( p%Nodes(I)%JointPos(3) >= z0 ) THEN
   
      ! loop through each member attached to the joint, getting the radius of its appropriate end
      DO J = 1, p%Nodes(I)%NConnections
      
         ! identify attached member and which end to us
         IF (p%Nodes(I)%ConnectionList(J) > 0) THEN         ! set up for end node 1
            member = p%Members(p%Nodes(I)%ConnectionList(J))
            MemberEndIndx = 1
         ELSE                                               ! set up for end node N+1
            member = p%Members(-p%Nodes(I)%ConnectionList(J))
            MemberEndIndx = member%NElements + 1
         END IF
      
         ! attached member cannot be modeled using WAMIT if we're to count it
         IF  (.NOT. member%PropPot) THEN
            
            ! Compute the signed area*outward facing normal of this member
            sgn = 1.0
            
            IF ( MemberEndIndx == 1 ) THEN
               sgn = -1.0                                ! Local coord sys points into member at starting node, so flip sign of local z vector
            ELSE
               sgn = 1.0                                 ! Local coord sys points out of member at ending node, so leave sign of local z vector
            END IF
            
            ! Compute the signed quantities for this member end, and add them to the joint values
            An = An + sgn* member%k*Pi*(member%RMG(MemberEndIndx))**2     ! area-weighted normal vector
            Vn = Vn + sgn* member%k*   (member%RMG(MemberEndIndx))**3     ! r^3-weighted normal vector used for mass
            I_n=I_n + sgn* member%k*Pi*(member%RMG(MemberEndIndx))**4     ! r^4-weighted normal vector used for moments of inertia
            
         END IF
         
      END DO   !J = 1, p%Nodes(I)%NConnections

      ! magnitudes of normal-weighted values
      Amag = sqrt(Dot_Product(An ,An))
      Vmag = sqrt(Dot_Product(Vn ,Vn))
      Imag = sqrt(Dot_Product(I_n,I_n))
      
          
      ! marine growth mass/inertia magnitudes
      p%Nodes(I)%m_MG = p%Nodes(I)%MGdensity * p%Nodes(I)%tMG * Amag         ! marine growth mass at joint
      Ir_MG_end = 0.25* p%Nodes(I)%MGdensity * p%Nodes(I)%tMG * Imag  ! radial moment of inertia magnitude
      Il_MG_end =  0.5* p%Nodes(I)%MGdensity * p%Nodes(I)%tMG * Imag  ! axial moment of inertia magnitude
      
      ! get rotation matrix for moment of inertia orientations
      RodrigMat(I_n, R_I, ErrStat, ErrMsg)

      ! globally-oreinted moment of inertia matrix for joint
      Irl_mat = 0
      Irl_mat(1,1) = Ir_MG_end
      Irl_mat(2,2) = Ir_MG_end
      Irl_mat(3,3) = Il_MG_end
      
      p%Nodes(I)%I_MG = MatMul( MatMul(R_I, Irl_mat), Transpose(R_I) ) ! final moment of inertia matrix for node
      
      
      
      ! pre-compute wave inertia loads on joint
      DO M=0,NStepWave
         ! The WaveAcc array has indices of (timeIndx, nodeIndx, vectorIndx), the nodeIndx needs to correspond to the total list of nodes for which
         ! the wave kinematics were generated.  We can use the nodeToLumpedIndx however for L_F_I and it's indices are (timeIndx, nodeIndx, vectorIndx)
         
         F_I      = 0.0
         IF ( (Vmag > 0.0) .AND. (.NOT. p%Nodes(I)%PropPot) ) THEN
            af =  p%WaveAcc(M, i,:)
            VnDotAf = Dot_Product(Vn,af)
            F_I(1:3) = ( p%Nodes(I)%JAxCa*AMfactor*VnDotAf / ( REAL( nCommon, ReKi ) * Vmag ) ) * Vn
         END IF
         m%L_F_I(M, :, i)   =  F_I
      END DO
      
   ELSE    ! if joint is below seabed
   
      p%Nodes(I)%m_MG = 0.0
      p%Nodes(I)%I_MG = 0.0
         m%L_F_I(:,:, i)   = 0.0
      
   END IF  ! nodes(I)%JointPos(3) >= z0
   
END DO ! looping through nodes that are joints
     
      
      
         ! Define initial guess for the system inputs here:

  !    u%DummyInput = 0


         ! Define system output initializations (set up mesh) here:
  
         
         ! Define initialization-routine output here:
         
         ! Initialize the outputs
         
   IF ( p%OutSwtch > 0) then  !@mhall: moved this "if" to after allocations
   
      CALL MrsnOUT_Init( InitInp, y, p, InitOut, ErrStat, ErrMsg )
      IF ( ErrStat > ErrID_None ) RETURN
      
         ! Determine if we need to perform output file handling
      
      IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN  
         CALL MrsnOUT_OpenOutput( Morison_ProgDesc%Name, TRIM(InitInp%OutRootName)//'.HD', p, InitOut, ErrStat, ErrMsg )
         IF ( ErrStat > ErrID_None ) RETURN
      END IF
      
   END IF  
   
   
      ! Write Summary information now that everything has been initialized.
      
   CALL WriteSummaryFile( InitInp%UnSum, InitInp%MSL2SWL, InitInp%WtrDpth, InitInp%NNodes, InitInp%Nodes, InitInp%NElements, InitInp%Members, p%NumOuts, p%OutParam, p%NMOutputs, p%MOutLst, p%distribToNodeIndx, p%NJOutputs, p%JOutLst, u%LumpedMesh, y%LumpedMesh,u%DistribMesh, y%DistribMesh, p%L_F_B, p%L_F_BF, p%D_F_B, p%D_F_BF, p%D_F_MG, InitInp%Gravity, ErrStat, ErrMsg ) !p%NDistribMarkers, distribMarkers, p%NLumpedMarkers, lumpedMarkers,
   IF ( ErrStat > ErrID_None ) RETURN  
      
         ! If you want to choose your own rate instead of using what the glue code suggests, tell the glue code the rate at which
         !   this module must be called here:
         
       !Interval = p%DT                                               
   !Contains:
   !   SUBROUTINE CleanUpInitOnErr
   !   IF (ALLOCATED(sw(1)%array))  DEALLOCATE(sw(1)%array, STAT=aviFail)
   !   END SUBROUTINE

END SUBROUTINE Morison_Init


SUBROUTINE RodrigMat(a, R, ErrStat, ErrMsg)
   ! calculates rotation matrix R to rotate unit vertical vector to direction of input vector a
   
   REAL(ReKi),      INTENT ( IN    )  :: a(3)    ! input vector
   REAL(ReKi),      INTENT ( INOUT )  :: R(3,3)  ! rotation matrix from Rodrigues's rotation formula
   INTEGER(IntKi),  INTENT(  OUT)     :: ErrStat     ! Error status of the operation
   CHARACTER(*),    INTENT(  OUT)     :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   REAL(ReKi)                         :: vec(3)  ! scaled and adjusted input vector
   REAL(ReKi)                         :: factor  ! denomenator used for scaling


   IF ((a(1) == 0) .AND. (a(2)==0)) THEN    ! return identity if vertical
         CALL EYE(R, ErrStat,ErrMsg)
      IF (a(3) < 0) THEN
         R = -R
      END IF
   
   ELSE   
      vec = a/SQRT(Dot_Product(a,a))
      vec(3) = vec(3) + 1
   
      factor = 2.0/Dot_Product(vec, vec)
      
      R(1,1) = factor*vec(1)*vec(1) - 1
      R(1,2) = factor*vec(1)*vec(2)
      R(1,3) = factor*vec(1)*vec(3)
      R(2,1) = factor*vec(2)*vec(1)
      R(2,2) = factor*vec(2)*vec(2) - 1
      R(2,3) = factor*vec(2)*vec(3)
      R(3,1) = factor*vec(3)*vec(1)
      R(3,2) = factor*vec(3)*vec(2)
      R(3,3) = factor*vec(3)*vec(3) - 1
   END IF
   
END SUBROUTINE RodrigMat


FUNCTION GetAlpha(R1,R2)
   ! calculates relative center of volume location for a (tapered) cylindrical element
   
   REAL(ReKi),                     INTENT    ( IN    )  :: R1  ! interior radius of element at node point
   REAL(ReKi),                     INTENT    ( IN    )  :: R2  ! interior radius of other end of part-element
   
   
   REAL(ReKi)  :: alpha  ! relative location of volumentric centroid between radius 1 and 2 (0= at radius 1, 1= at radius 2)
   
   alpha = (R1*R1 + 2*R1*R2 + 3*R2*R2)/4/(R1*R1 + R1*R2 + R2*R2)

   return alpha

END FUNCTION GetAlpha


SUBROUTINE AllocateNodeLoadVariables(m, NNodes, ErrStat, ErrMsg )
      TYPE(Morison_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables            
      INTEGER(IntKi),                    INTENT(IN   )  :: NNodes      ! number of nodes in node list
      INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      

   ALLOCATE ( m%F_D(3,NNodes), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for D_F_D array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%F_D = 0.0_ReKi
      
   ALLOCATE ( m%F_B( 6, NNodes ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the buoyancy forces/moments array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%F_B = 0.0
   
   
   ALLOCATE ( m%F_MG( 6, NNodes ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the marine growth weight forces/moments array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%F_MG = 0.0
   
   ALLOCATE ( m%F_BF( 6, NNodes ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the buoyancy due to flooding forces/moments array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%F_BF = 0.0
   
  
    ALLOCATE ( m%AM_M( 3, 3, NNodes ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the added mass of flooded fluid.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%AM_M = 0.0
   
   ALLOCATE ( m%AM_MG( NNodes ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the added mass of marine growth.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%AM_MG = 0.0
   
   ALLOCATE ( m%AM_F( NNodes ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the added mass of flooded fluid.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%AM_F = 0.0
   
      ALLOCATE ( m%F_AM(6,NNodes), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for D_F_AM array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%F_AM = 0.0_ReKi
   
   ALLOCATE ( m%F_AM_M(6,NNodes), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for D_F_AM_M array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%F_AM_M = 0.0_ReKi
   
   ALLOCATE ( m%F_AM_MG(6,NNodes), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for D_F_AM_MG array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%F_AM_MG = 0.0_ReKi
   
   ALLOCATE ( m%dragConst( NNodes ), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for the drag constants.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%dragConst = 0.0
   
   
   ALLOCATE ( m%FV(3,NNodes), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for D_FV array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%FV = 0.0_ReKi
   
   ALLOCATE ( m%FA(3,NNodes), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for D_FA array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%FA = 0.0_ReKi
   
   ALLOCATE ( m%FDynP(NNodes), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for D_FDynP array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   m%FDynP = 0.0_ReKi
   

END SUBROUTINE AllocateNodeLoadVariables



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
      
      
      TYPE(Morison_MemberType) :: mem    ! the current member
      INTEGER                  :: N
      REAL(ReKi)               :: dl
      REAL(ReKi)               :: vec(3)
      REAL(ReKi)               :: phi    ! member tilt angle
      REAL(ReKi)               :: beta   ! member tilt heading
      REAL(ReKi)               :: cosPhi
      REAL(ReKi)               :: sinPhi
      REAL(ReKi)               :: tanPhi
      REAL(ReKi)               :: sinBeta
      REAL(ReKi)               :: cosBeta
      REAL(ReKi)               :: z1
      REAL(ReKi)               :: z2
      REAL(ReKi)               :: r1
      REAL(ReKi)               :: r2
      REAL(ReKi)               :: m     ! shorthand for taper including marine growth of segment i
      REAL(ReKi)               :: Rmid  
      REAL(ReKi)               :: RmidMG
      REAL(ReKi)               :: Rmidin
      REAL(ReKi)               :: Lmid  
      REAL(ReKi)               :: Imat
      REAL(ReKi)               :: Fl    ! various element axial force 
      REAL(ReKi)               :: Fr    ! various element radial force
      REAL(ReKi)               :: M     ! various element radial moment
      REAL(ReKi)               :: h0    ! distances along cylinder centerline from point 1 to the waterplane
      REAL(ReKi)               :: rh    ! radius of cylinder at point where its centerline crosses the waterplane
      REAL(ReKi)               :: l1    ! distance from cone end to bottom node
      REAL(ReKi)               :: Vs    ! segment submerged volume
      REAL(ReKi)               :: a0    ! waterplane ellipse shape
      REAL(ReKi)               :: b0    
      REAL(ReKi)               :: cr    ! centroid of segment submerged volume relative to its lower node
      REAL(ReKi)               :: cl 
      REAL(ReKi)               :: cx 
      REAL(ReKi)               :: cz 
      REAL(ReKi)               :: pwr   ! exponent for buoyancy node distribution smoothing
      REAL(ReKi)               :: alpha ! final load distribution factor for element
      REAL(ReKi)               :: Fb    !buoyant force
      REAL(ReKi)               :: Fr    !radial component of buoyant force
      REAL(ReKi)               :: Fl    !axial component of buoyant force
      REAL(ReKi)               :: M     !moment induced about the center of the cylinder's bottom face
      REAL(ReKi)               :: BuoyF(3) ! buoyancy force vector aligned with an element
      REAL(ReKi)               :: BuoyM(3) ! buoyancy moment vector aligned with an element
      
      
      
      
      
         ! Initialize ErrStat
         
      ErrStat = ErrID_None         
      ErrMsg  = ""               
      
      
         ! Compute outputs here:
      
      ! We need to attach the distributed drag force (D_F_D), distributed inertial force (D_F_I), and distributed dynamic pressure force (D_F_DP) to the Misc type so that we don't need to
      ! allocate their data storage at each time step!  If we could make them static local variables (like in C) then we could avoid adding them to the OtherState datatype.  
      ! The same is true for the lumped drag (L_F_D) and the lumped dynamic pressure (L_F_DP)
         
      DO J = 1, y%Mesh%Nnodes
         
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
               m%D_F_AM_F(:,J)  = -p%D_AM_F(J)*u%DistribMesh%TranslationAcc(I,J)  !@mhall: should D_F_AM_F(:,J) be D_F_AM_F(I,J) ?
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
      
      
      
      
   !!! below from GenerateLumpedLoads - need to make sure nothing is forgotten <<<<<<
   IF (.NOT. node%PropPot ) THEN
      k =  sgn * node%R_LToG(:,3)
      CALL LumpDynPressure( nodeIndx, node%JAxCp, k, node%R, node%tMG, NStepWave, WaveDynP, F_DP, ErrStat, ErrMsg)
      
         ! For buoyancy calculations we need to adjust the Z-location based on MSL2SWL. If MSL2SWL > 0 then SWL above MSL, and so we need to place the Z value at a deeper position.  
         !   SWL is at Z=0 for buoyancy calcs, but geometry was specified relative to MSL (MSL2SWL = 0) 
     ! CALL LumpBuoyancy( sgn, densWater, node%R, node%tMG, node%JointPos(3) - MSL2SWL, node%R_LToG, gravity, F_B  ) 
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
    
      
<<<<< most of above likely not used anymore!!!  
    

   ! Clear saved loads on each node since these are accumulated
   !@mhall: Faster way to do this? Is this list right?
   DO i = 1, p%NNodes
    DO J=1,6
      F_D    (i,J) = 0.0
      F_I    (i,J) = 0.0
      F_B    (i,J) = 0.0
      F_FB   (i,J) = 0.0
      F_FBI  (i,J) = 0.0
      F_MG   (i,J) = 0.0
      F_MGI  (i,J) = 0.0
      F_AM   (i,J) = 0.0
      F_AM_M (i,J) = 0.0
      F_AM_MG(i,J) = 0.0
      F_AM_F (i,J) = 0.0
    END DO
   END DO
   
   
   
   
   
   ! ==============================================================================================
   ! Calculate instantaneous loads on each member except for the hydrodynamic loads on member ends.
   ! This covers aspects of the load calculations previously in CreateDistributedMesh.  


   ! Loop through each member
   DO im = 1, p%NMembers
      
      N = m%Members(im)%NElements
      
      mem = m%Members(im)   !@mhall: does this have much overhead?
      
      
      ! calculate isntantaneous incline angle and heading, and related trig values
      vec = u%Mesh%TranslationDisp(:, mem%NodeIndx(N+1)) - u%Mesh%TranslationDisp(:, mem%NodeIndx(1  )) 
      
      phi = arccos(vec(3)/SQRT(Dot_Product(vec,vec)))  ! incline angle   
      beta = arctan2(vec(2), vec(1))                   ! heading of incline
      
      sinPhi = sin(phi)
      cosPhi = cos(phi)  
      tanPhi = tan(phi)
      
      sinBeta = sin(beta)
      cosBeta = cos(beta) 
      
      


      DO i =1,N    ! loop through member elements

        ! save some commonly used variables
     
        dl = mem%dl

        z1 = u%Mesh%TranslationDisp(3, mem%NodeIndx(i  ))   ! get node z locations from input mesh
        z2 = u%Mesh%TranslationDisp(3, mem%NodeIndx(i+1))
        r1 = m%Members(im)%RMG(i  )
        r2 = m%Members(im)%RMG(i+1)
        
        m = mem%m_mg(i)
            
        ! should i_floor theshold be applied to below calculations to avoid wasting time on computing zero-valued things? <<<<<
        ! should lumped half-element coefficients get combined at initialization? <<<
              
        ! ------------------ marine growth --------------------  

        ! lower node
        m%F_MG(mem%NodeIndx(i  ), 3) = m%F_MG(mem%NodeIndx(i  ), 3) + mem%m_mg_l(i)*g ! weight force
        m%F_MG(mem%NodeIndx(i  ), 4) = m%F_MG(mem%NodeIndx(i  ), 4) - mem%m_mg_l(i)*g * mem%h_cmg_l(i)* sinPhi * sinBeta! weight force
        m%F_MG(mem%NodeIndx(i  ), 5) = m%F_MG(mem%NodeIndx(i  ), 5) + mem%m_mg_l(i)*g * mem%h_cmg_l(i)* sinPhi * cosBeta! weight force
        ! upper node
        m%F_MG(mem%NodeIndx(i+1), 3) = m%F_MG(mem%NodeIndx(i+1), 3) + mem%m_mg_u(i)*g ! weight force
        m%F_MG(mem%NodeIndx(i+1), 4) = m%F_MG(mem%NodeIndx(i+1), 4) - mem%m_mg_u(i)*g * mem%h_cmg_u(i)* sinPhi * sinBeta! weight force
        m%F_MG(mem%NodeIndx(i+1), 5) = m%F_MG(mem%NodeIndx(i+1), 5) + mem%m_mg_u(i)*g * mem%h_cmg_u(i)* sinPhi * cosBeta! weight force

        ! lower node
        Imat = diag(mem%I_rmg_l(i), mem%I_rmg_l(i), mem%I_lmg_l(i))
        m%F_MGI(mem%NodeIndx(i  ), 1:3) = m%F_MGI(mem%NodeIndx(i  ), 1:3) - a_struct * mem%m_mg_l(i)
        m%F_MGI(mem%NodeIndx(i  ), 4:6) = m%F_MGI(mem%NodeIndx(i  ), 4:6) - cross (a_struct * mem%m_mg_l(i), mem%h_cmg_l(i) * mem%k) + C * Imat * Ctrans  
        ! upper node
        Imat = diag(mem%I_rmg_u(i), mem%I_rmg_u(i), mem%I_lmg_u(i))
        m%F_MGI(mem%NodeIndx(i+1), 1:3) = m%F_MGI(mem%NodeIndx(i+1), 1:3) - a_struct * mem%m_mg_u(i)
        m%F_MGI(mem%NodeIndx(i+1), 4:6) = m%F_MGI(mem%NodeIndx(i+1), 4:6) - cross (a_struct * mem%m_mg_u(i), mem%h_cmg_u(i) * mem%k) + C * Imat * Ctrans  
           

        ! ------------------ flooded ballast inertia ---------------------

        ! lower node
        Imat = diag(mem%I_rfb_l(i), mem%I_rfb_l(i), mem%I_lfb_l(i))
        m%F_FBI(mem%NodeIndx(i  ), 1:3) = m%F_FBI(mem%NodeIndx(i  ), 1:3) - a_struct * mem%m_fb_l(i)
        m%F_FBI(mem%NodeIndx(i  ), 4:6) = m%F_FBI(mem%NodeIndx(i  ), 4:6) - cross (a_struct * mem%m_fb_l(i), mem%h_cfb_l(i) * mem%k) + C * Imat * Ctrans  
        ! upper node
        Imat = diag(mem%I_rfb_u(i), mem%I_rfb_u(i), mem%I_lfb_u(i))
        m%F_FBI(mem%NodeIndx(i+1), 1:3) = m%F_FBI(mem%NodeIndx(i+1), 1:3) - a_struct * mem%m_fb_u(i)
        m%F_FBI(mem%NodeIndx(i+1), 4:6) = m%F_FBI(mem%NodeIndx(i+1), 4:6) - cross (a_struct * mem%m_fb_u(i), mem%h_cfb_u(i) * mem%k) + C * Imat * Ctrans  
        
           
        ! ------------------ flooded ballast weight ---------------------
        
        ! fully filled elements
        if (floodstatus(i) == 1) then
        
           ! forces and moment in tilted coordinates about node i
           Fl = mem%Cfl_fb(i)*cosPhi     
           Fr = mem%Cfr_fb(i)*sinPhi     
           M  = mem%CM0_fb(i)*sinPhi - Fr(i)*alpha_fb_star(i)*dl
           
           ! calculate full vector and distribute to nodes
           DistributeElementLoads(Fl, Fr, M, sinPhi, cosPhi, SinBeta, cosBeta, (1-mem%alpha_fb_star(i)), m%F_FB(mem%NodeIndx(i), :), m%F_FB(mem%NodeIndx(i+1), :))
           
           
        ! partially filled element
        else if (floodstatus(i) == 2) then
           
           ! forces and moment in tilted coordinates about node i
           Fl = mem%Cfl_fb(i)*cosPhi     
           Fr = mem%Cfr_fb(i)*sinPhi     
           M  = mem%CM0_fb(i)*sinPhi + Fr(i)*(1 - alpha_fb_star(i))*dl
        
           ! calculate full vector and distribute to nodes
           DistributeElementLoads(Fl, Fr, M, sinPhi, cosPhi, SinBeta, cosBeta, mem%alpha_fb_star(i), m%F_FB(mem%NodeIndx(i), :), m%F_FB(mem%NodeIndx(i-1), :))
                
        
        ! no load for unflooded element or element fully below seabed
        
        end if
        
           
        ! ------------------- buoyancy loads... ------------------------

        if (z1 <= 0) then   ! if segment is at least partially submerged ...
              
              
           if (z1*z2 <= 0) then ! special calculation if the slice is partially submerged
        
              h0 = -z1/cosPhi             ! distances along element centerline from point 1 to the waterplane
              rh = r1 + h0*m    ! radius of element at point where its centerline crosses the waterplane
              
              if (abs(m) < 0.0001) then      ! untapered cylinder case

                 Vs =    Pi*r1*r1*h0   ! volume of total submerged portion
                 
                 cr = 0.25*r1*r1*tanPhi/h0
                 cl = 0.5*h0 + 0.125*r1*r1*tanPhi*tanPhi/h0

                 cx = cr*cosPhi + cl*sinPhi
                 cz = cl*cosPhi - cr*sinPhi
                    
                 !alpha0 = 0.5*h0/dl            ! force distribution between end nodes
                 
              else       ! inclined tapered cylinder case (note I've renamed r0 to rh here!!)
              
                 l1 = r1/m  ! distance from cone end to bottom node
                                
                 ! waterplane ellipse shape
                 b0 = rh/sqrt(1 - m^2 * tanPhi**2)
                 a0 = rh/((1 - m^2*tanPhi**2)*cosPhi)             ! simplified from what's in ConicalCalcs.ipynb
              
                 ! segment submerged volume
                 Vs = pi*(a0*b0*rh*cosPhi - l1^3*m^3)/(3*m)
                 
                 ! centroid of segment submerged volume (relative to bottom node)
                 cx = -0.25*(3*a0*b0*rh*rh*(m**2 + 1)*cosPhi + 3.0*l1**4*m**4*(m**2*tanPhi**2 - 1) + 4*l1*m*(m**2*tanPhi**2 - 1)*(a0*b0*rh*cosPhi - 1.0*l1**3*m**3))*sin(phi)/(m*(m**2*tanPhi**2 - 1)*(a0*b0*rh*cosPhi - l1**3*m**3))
                 cz = 0.25*(-4.0*a0*b0*l1*m*rh*cosPhi + 3.0*a0*b0*rh*rh*cosPhi + 1.0*l1**4*m**4)*cosPhi/(m*(a0*b0*rh*cosPhi - l1**3*m**3))
                             
                 !alpha0 = (r1*r1 + 2*r1*r2 + 3*r2**2)/4/(r1*r1 + r1*r2 + r2**2)  ! this can be precomputed
              
              end if

              pwr = 1
              alpha    = (1-mem%alpha(i))*z1**pwr/(-mem%alpha(i)*z2**pwr + (1-mem%alpha(i))*z1**pwr)

              Fb  = Vs*rho*g       !buoyant force
              Fr  = -Fb*sinPhi     !radial component of buoyant force
              Fl  = Fb*cosPhi      !axial component of buoyant force
              M = -Fb*cx           !moment induced about the center of the cylinder's bottom face

              ! calculate (imaginary) bottom plate forces/moment to subtract from displacement-based values
              Fl  = Fl  + rho*g*z1* Pi *r1*r1        
              M0  = M   + rho*g* sinPhi * Pi/4*r1**4       


              ! reduce taper-based moment to remove (not double count) radial force distribution to each node 
              M  = M0 - Fr*alpha*dl + Fr*dl

              BuoyF(1) =  Fr*cosBeta
              BuoyF(2) =  Fr*sinBeta
              BuoyF(3) =  Fl
              BuoyM(1) = -M*sinBeta
              BuoyM(2) =  M*cosBeta
              BuoyM(3) =  0

              ! distribute force and moment to each adjacent node of the element BELOW THIS ONE
              m%F_B(mem%NodeIndx(i-1), 1:3) = MatMul(C, BuoyF )*(1-alpha)
              m%F_B(mem%NodeIndx(i-1), 4:6) = Matmul(C, BuoyMT)*(1-alpha)
              
              m%F_B(mem%NodeIndx(i  ), 1:3) = MatMul(C, BuoyF )*alpha
              m%F_B(mem%NodeIndx(i  ), 4:6) = Matmul(C, BuoyMT)*alpha


           else ! normal, fully submerged case
              
              Fl = -2.0*Pi*m*rho*g*dl*( z1*r1 + 0.5*(z1*m + r1*cosPhi)*dl + 1/3*m*cosPhi*dl*dl )   ! from CylinderCalculationsR1.ipynb
              
              Fr = -Pi*rho*g*dl*(r1*r1 + m*r1*dl + 1/3*m**2*dl**2)*sinPhi                          ! from CylinderCalculationsR1.ipynb
                 
              M0	= -Pi*dl*g*rho*(3*dl**3*m**4 + 3*dl**3*m**2 + 12*dl**2*m**3*r1 + 8*dl**2*m*r1 + 18*dl*m**2*r1*r1 + 6*dl*r1*r1 + 12*m*r1**3)*sin(phi)/12   ! latest from CylinderCalculationsR1.ipynb
                 
              ! precomputed as mem%alpha(i) ... alpha0 = (r1*r1 + 2*r1*r2 + 3*r2**2)/4/(r1*r1 + r1*r2 + r2**2)
              
              z1d = -min(0,z1)
              z2d = -min(0,z2)
                    
              alpha = mem%alpha(i)*z2d/(mem%alpha(i)*z2d+(1-mem%alpha(i))*z1d)
                             
              
              ! reduce moment to remove (not double count) radial force distribution to each node
              M = M0 - Fr*alpha*dl
              
              BuoyF(1) =  Fr*cosBeta
              BuoyF(2) =  Fr*sinBeta
              BuoyF(3) =  Fl
              BuoyM(1) = -M*sinBeta
              BuoyM(2) =  M*cosBeta
              BuoyM(3) =  0
              
              ! distribute force and moment to each adjacent node
              m%F_B(mem%NodeIndx(i  ), 1:3) = MatMul(C, BuoyF)*(1-alpha)
              m%F_B(mem%NodeIndx(i  ), 4:6) = Matmul(C, BuoyM)*(1-alpha)
              
              m%F_B(mem%NodeIndx(i+1), 1:3) = MatMul(C, BuoyF)*alpha
              m%F_B(mem%NodeIndx(i+1), 4:6) = Matmul(C, BuoyM)*alpha

           end if  ! submergence cases
             
        end if ! element at least partially submerged

      END DO ! i =1,N    ! loop through member elements       

       
      ! Any end plate loads that are modeled on a per-member basis
      
      ! reassign convenience variables to correspond to member ends
      z1 = u%Mesh%TranslationDisp(3, m%Members(im)%NodeIndx(1  )) 
      z2 = u%Mesh%TranslationDisp(3, m%Members(im)%NodeIndx(N+1))

      ! Water ballast buoyancy 
      if (z1 >= -wtrDpth) then   ! end load only if end is above seabed
         
         ! if member is fully flooded
         if (mem%z_overfill > 0) then 
            Fl = -mem%FillDens * g * pi *mem%Rin(  1)**2* (mem%z_overfill + max(z2-z1, 0))
            M  = -mem%FillDens * g * pi *0.25*mem%Rin(  1)**4*sinPhi
            AddEndLoad(Fl, M, sinPhi, cosPhi, SinBeta, cosBeta, m%F_FB(  1,:))
            
            Fl =  mem%FillDens * g * pi *mem%Rin(N+1)**2* (mem%z_overfill + max(z1-z2, 0))
            M  =  mem%FillDens * g * pi *0.25*mem%Rin(N+1)**4*sinPhi            
            AddEndLoad(Fl, M, sinPhi, cosPhi, SinBeta, cosBeta, m%F_FB(N+1,:))
            
         ! if member is partially flooded
         else if (l_fill > 0) then 
            Fl = -mem%FillDens * g * pi *Rin(1)**2*l_fill*cosPhi
            M  = -mem%FillDens * g * pi *0.25*Rin(1)**4*sinPhi
            AddEndLoad(Fl, M, sinPhi, cosPhi, SinBeta, cosBeta, m%F_FB(1,:))

         ! no load if member is not flooded at all
         end if
      ! no load if end is below seabed
      end if

      ! --- no inertia loads from water ballast modeled on ends

   end do ! im - looping through members
      
      
   !@mhall: hydrodynamic loads not yet implemented, below not yet adjusted   

      
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



! Takes loads on node i in element tilted frame and converts to 6DOF loads at node i and adjacent node
SUBROUTINE DistributeElementLoads(Fl, Fr, M, sinPhi, cosPhi, SinBeta, cosBeta, alpha, Fi, F2)
   
   REAL(ReKi),                     INTENT    ( IN    )  :: Fl        ! (N)   axial load about node i
   REAL(ReKi),                     INTENT    ( IN    )  :: Fr        ! (N)   radial load about node i in direction of tilt
   REAL(ReKi),                     INTENT    ( IN    )  :: M         ! (N-m) radial moment about node i, positive in direction of tilt angle
   REAL(ReKi),                     INTENT    ( IN    )  :: sinPhi    ! trig functions of  tilt angle 
   REAL(ReKi),                     INTENT    ( IN    )  :: cosPhi   
   REAL(ReKi),                     INTENT    ( IN    )  :: sinBeta   ! trig functions of heading of tilt
   REAL(ReKi),                     INTENT    ( IN    )  :: cosBeta  
   REAL(ReKi),                     INTENT    ( IN    )  :: alpha     ! fraction of load staying with node i (1-alpha goes to other node)  
   
   REAL(ReKi),                     INTENT    ( OUT   )  :: Fi(6)   ! (N, Nm) force/moment vector for node i
   REAL(ReKi),                     INTENT    ( OUT   )  :: F2(6)   ! (N, Nm) force/moment vector for the other node (whether i+1, or i-1)
         

   Fi(1) =  cosBeta*(Fl*sinPhi + Fr*cosPhi)*alpha
   Fi(2) = -sinBeta*(Fl*sinPhi + Fr*cosPhi)*alpha
   Fi(3) =          (Fl*cosPhi + Fr*sinPhi)*alpha
   Fi(4) =  sinBeta * M                    *alpha
   Fi(5) =  cosBeta * M                    *alpha
   Fi(6) = 0.0
      
   Fi(1) =  cosBeta*(Fl*sinPhi + Fr*cosPhi)*(1-alpha)
   Fi(2) = -sinBeta*(Fl*sinPhi + Fr*cosPhi)*(1-alpha)
   Fi(3) =          (Fl*cosPhi + Fr*sinPhi)*(1-alpha)
   Fi(4) =  sinBeta * M                    *(1-alpha)
   Fi(5) =  cosBeta * M                    *(1-alpha)
   Fi(6) = 0.0

END SUBROUTINE DistributeElementLoads


! Takes loads on end node i and converts to 6DOF loads, adding to the nodes existing loads
SUBROUTINE AddEndLoad(Fl, M, sinPhi, cosPhi, SinBeta, cosBeta, Fi)
   
   REAL(ReKi),                     INTENT    ( IN    )  :: Fl        ! (N)   axial load about node i
   REAL(ReKi),                     INTENT    ( IN    )  :: M         ! (N-m) radial moment about node i, positive in direction of tilt angle
   REAL(ReKi),                     INTENT    ( IN    )  :: sinPhi    ! trig functions of  tilt angle 
   REAL(ReKi),                     INTENT    ( IN    )  :: cosPhi   
   REAL(ReKi),                     INTENT    ( IN    )  :: sinBeta   ! trig functions of heading of tilt
   REAL(ReKi),                     INTENT    ( IN    )  :: cosBeta  
   REAL(ReKi),                     INTENT    ( INOUT )  :: Fi(6)   ! (N, Nm) force/moment vector for end node i
   
   Fi(1) = Fi(1) + Fl*sinPhi*cosBeta
   Fi(2) = Fi(2) + Fl*sinPhi*sinBeta
   Fi(3) = Fi(3) + Fl*cosPhi
   Fi(4) = Fi(4) + M*sinBeta
   Fi(5) = Fi(5) - M*cosBeta
   
END SUBROUTINE AddEndLoad


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
