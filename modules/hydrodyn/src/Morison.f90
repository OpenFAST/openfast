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
   USE SeaSt_WaveField
  ! USE HydroDyn_Output_Types
   USE NWTC_Library
   USE YawOffset

   
   IMPLICIT NONE
   
   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: Morison_ProgDesc = ProgDesc( 'Morison', '', '' )

   INTERFACE Morison_DirCosMtrx
      MODULE PROCEDURE Morison_DirCosMtrx_Spin
      MODULE PROCEDURE Morison_DirCosMtrx_noSpin
   END INTERFACE

      ! ..... Public Subroutines ...................................................................................................
   PUBLIC:: Morison_GenerateSimulationNodes
   
   PUBLIC :: Morison_Init                           ! Initialization routine
   PUBLIC :: Morison_CalcOutput                     ! Routine for computing outputs
   PUBLIC :: Morison_UpdateDiscState                ! Routine for updating discrete states
   
CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Morison_DirCosMtrx_Spin( pos0, pos1, spin, DirCos )
! Compute the direction cosine matrix given two end points and a spin angle for rectangular members
! Left multiplying DirCos with a vector in element local sys returns vector in global sys
! Updated to match the convention in SubDyn for consistency

   REAL(ReKi), INTENT( IN    )  ::   pos0(3), pos1(3), spin
   Real(ReKi), INTENT(   OUT )  ::   DirCos(3,3)
   Real(DbKi)                   ::   Le, Lexy
   Real(DbKi)                   ::   dx, dy, dz
   Real(DbKi)                   ::   Rspin(3,3)

   DirCos = 0.0
   Rspin  = 0.0

   dx = pos1(1) - pos0(1)
   dy = pos1(2) - pos0(2)
   dz = pos1(3) - pos0(3)

   Lexy = sqrt( dx*dx + dy*dy )

   IF ( EqualRealNos(Lexy, 0.0_DbKi) ) THEN
      IF (dz > 0) THEN
         DirCos(1,1) =  1.0
         DirCos(2,2) =  1.0
         DirCos(3,3) =  1.0
      ELSE
         DirCos(1,1) =  1.0
         DirCos(2,2) = -1.0
         DirCos(3,3) = -1.0
      END IF
   ELSE
      Le = sqrt( dx*dx + dy*dy + dz*dz )

      DirCos(1, 1) =  dy/Lexy
      DirCos(1, 2) =  dx*dz/(Lexy*Le)
      DirCos(1, 3) =  dx/Le

      DirCos(2, 1) = -dx/Lexy
      DirCos(2, 2) =  dy*dz/(Lexy*Le)
      DirCos(2, 3) =  dy/Le

      DirCos(3, 1) =  0.0
      DirCos(3, 2) = -Lexy/Le
      DirCos(3, 3) =  dz/Le
   END IF

   IF ( .not. EqualRealNos(spin, 0.0) ) THEN
      ! Spin the member about its axis first
      Rspin(1,1) =  cos(spin)
      Rspin(2,2) =  Rspin(1,1)
      Rspin(1,2) = -sin(spin)
      Rspin(2,1) = -Rspin(1,2)
      Rspin(3,3) =  1.0
      DirCos  =  matmul(DirCos,Rspin)
   END IF

END SUBROUTINE Morison_DirCosMtrx_Spin

SUBROUTINE Morison_DirCosMtrx_noSpin( pos0, pos1, DirCos )
! Compute the direction cosine matrix given two end points without spin for cylindrical members
! Left multiplying DirCos with a vector in element local sys returns vector in global sys
! Updated to match the convention in SubDyn for consistency

   REAL(ReKi), INTENT( IN    )  ::   pos0(3), pos1(3)
   Real(ReKi), INTENT(   OUT )  ::   DirCos(3,3)
   Real(DbKi)                   ::   Le, Lexy
   Real(DbKi)                   ::   dx, dy, dz

   DirCos = 0.0

   dx = pos1(1) - pos0(1)
   dy = pos1(2) - pos0(2)
   dz = pos1(3) - pos0(3)

   Lexy = sqrt( dx*dx + dy*dy )

   IF ( EqualRealNos(Lexy, 0.0_DbKi) ) THEN
      IF (dz > 0) THEN
         DirCos(1,1) =  1.0
         DirCos(2,2) =  1.0
         DirCos(3,3) =  1.0
      ELSE
         DirCos(1,1) =  1.0
         DirCos(2,2) = -1.0
         DirCos(3,3) = -1.0
      END IF
   ELSE
      Le = sqrt( dx*dx + dy*dy + dz*dz )

      DirCos(1, 1) =  dy/Lexy
      DirCos(1, 2) =  dx*dz/(Lexy*Le)
      DirCos(1, 3) =  dx/Le

      DirCos(2, 1) = -dx/Lexy
      DirCos(2, 2) =  dy*dz/(Lexy*Le)
      DirCos(2, 3) =  dy/Le

      DirCos(3, 1) =  0.0
      DirCos(3, 2) = -Lexy/Le
      DirCos(3, 3) =  dz/Le
   END IF

END SUBROUTINE Morison_DirCosMtrx_noSpin

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
!----------------------------------------------------------------------------------------------------------------------------------
subroutine GetOrientationAngles(p1, p2, phi, sinPhi, cosPhi, tanPhi, sinBeta, cosBeta, k_hat, errStat, errMsg)
   real(ReKi),   intent(in   ) :: p1(3),p2(3)
   real(ReKi),   intent(  out) :: phi, sinPhi, cosPhi, tanPhi, sinBeta, cosBeta, k_hat(3)
   integer,      intent(  out) :: errStat              ! returns a non-zero value when an error occurs  
   character(*), intent(  out) :: errMsg               ! Error message if errStat /= ErrID_None
   character(*), parameter     :: RoutineName = 'GetOrientationAngles'
   
   real(ReKi) :: vec(3), vecLen, vecLen2D, beta
   
      ! Initialize errStat
         
   errStat = ErrID_None         
   errMsg  = "" 
   
            ! calculate isntantaneous incline angle and heading, and related trig values
         ! the first and last NodeIndx values point to the corresponding Joint nodes idices which are at the start of the Mesh
         vec      = p2 - p1   
         vecLen   = SQRT(Dot_Product(vec,vec))
         vecLen2D = SQRT(vec(1)**2+vec(2)**2)
         if ( vecLen < 0.000001 ) then
            call SeterrStat(ErrID_Fatal, 'An element of the Morison structure has co-located endpoints!  This should never occur.  Please review your model.', errStat, errMsg, RoutineName )
            return
         else
            k_hat = vec / vecLen 
            phi   = atan2(vecLen2D, vec(3))  ! incline angle   
         end if
         if ( EqualRealNos(phi, 0.0_ReKi) ) then
            beta = 0.0_ReKi
         else
            beta = atan2(vec(2), vec(1))                    ! heading of incline     
         endif
         sinPhi  = sin(phi)
         cosPhi  = cos(phi)  
         tanPhi  = tan(phi)     
         sinBeta = sin(beta)
         cosBeta = cos(beta)
         
end subroutine GetOrientationAngles
!----------------------------------------------------------------------------------------------------------------------------------
!function to return conical taper geometry calculations (volume and center of volume)
SUBROUTINE CylTaperCalc(R1, R2, H, taperV, h_c)
   REAL(ReKi),                     INTENT    ( IN    )  :: R1
   REAL(ReKi),                     INTENT    ( IN    )  :: R2
   REAL(ReKi),                     INTENT    ( IN    )  :: H
   REAL(ReKi),                     INTENT    ( OUT   )  :: taperV   ! volume of tapered section
   REAL(ReKi),                     INTENT    ( OUT   )  :: h_c    ! center of mass offset from first node
   
   real(ReKi) :: m
   
   m = (R2-R1)/H

   if ( EqualRealNos(R1, R2) ) then             ! if just a cylinder
      taperV = abs(pi*R1*R1*H)
      h_c = H/2.0
   elseif ( EqualRealNos(R1,0.0_ReKi) ) then             ! seperate this case out because it gives a divide by zero in general formula
      taperV = abs(1.0/3.0*pi*R2*R2*H)                                            ! cone volume
      h_c = 3.0/4.0*H                                                        ! from base           
   else
     taperV = abs(pi/3.0/m*(R2**3 - R1**3))
     h_c = H*(R1**2 + 2*R1*R2 + 3*R2**2)/4.0/(R1**2 + R1*R2 + R2**2) !( coneV*1./4.*coneH - coneVtip*(1./4.*(coneH-H) + H) )/ taperV ! from base
   end if
   
END SUBROUTINE CylTaperCalc
!----------------------------------------------------------------------------------------------------------------------------------
!function to return pyramidal taper geometry calculations (volume and center of volume)
SUBROUTINE RecTaperCalc(a0, a1, b0, b1, H, taperV, h_c)
   REAL(ReKi),                     INTENT    ( IN    )  :: a0
   REAL(ReKi),                     INTENT    ( IN    )  :: a1
   REAL(ReKi),                     INTENT    ( IN    )  :: b0
   REAL(ReKi),                     INTENT    ( IN    )  :: b1
   REAL(ReKi),                     INTENT    ( IN    )  :: H
   REAL(ReKi),                     INTENT    ( OUT   )  :: taperV   ! volume of tapered section
   REAL(ReKi),                     INTENT    ( OUT   )  :: h_c      ! center of mass offset from first node

   REAL(ReKi) :: a0b0, a0b1, a1b0, a1b1, tmp   

   a0b0 = a0*b0
   a0b1 = a0*b1
   a1b0 = a1*b0
   a1b1 = a1*b1

   tmp    = 2.0*a0b0 + a0b1 + a1b0 + 2.0*a1b1
   if ( EqualRealNos(tmp, 0.0_ReKi) ) then
      taperV = 0.0_ReKi
      h_c    = 0.0_ReKi
   else
      taperV = abs(H/6.0*tmp)
      h_c    = 0.5*H*(a0b0 + a0b1 + a1b0 + 3.0*a1b1)/tmp ! from base
   end if
   
END SUBROUTINE RecTaperCalc
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE CylInertia(R1, R2, H, rho, Il, Ir)
   REAL(ReKi),                     INTENT    ( IN    )  :: R1
   REAL(ReKi),                     INTENT    ( IN    )  :: R2
   REAL(ReKi),                     INTENT    ( IN    )  :: H
   REAL(ReKi),                     INTENT    ( IN    )  :: rho ! density of material
   REAL(ReKi),                     INTENT    ( OUT   )  :: Il
   REAL(ReKi),                     INTENT    ( OUT   )  :: Ir
   
   real(ReKi) :: m, Ir_tip, h_c
   
   m = (R2-R1)/H

   if ( EqualRealNos(R1, R2) ) then             ! if just a cylinder
      Ir = abs(1.0/12.0* rho*pi*R1*R1*H *(3.0*R1*R1 + 4.0*H*H)) ! radial inertia about node 1 
      Il = abs(0.5* rho*pi*R1*R1*H *R1*R1)    
   ELSEIF ( EqualRealNos(R1,0.0_ReKi) ) then        ! seperate this case out because it gives a divide by zero in general formula
      Ir = abs(rho*pi*(1.0/20.0/m + 1.0/5.0/m**3) * R2**5)      
      Il = abs(1.0/10.0*rho*pi/m*R2**5)            
   ELSE 
     h_c = H*(R1**2 + 2*R1*R2 + 3*R2**2)/4.0/(R1**2 + R1*R2 + R2**2) 
     !l_c = R1/M + (R2-R1)/m *(R1**2 + 2*R1*R2 + 3*R2**2)/4/(R1**2 + R1*R2 + R2**2) 
     Ir_tip = abs(pi/20.0 *rho/m*(1.0 + 4.0/m**2) * (R2**5 - R1**5))                    ! radial moment of inertia about tip of cone
     Ir = abs(Ir_tip - rho/3.0/m*pi*(R2**3-R1**3) * (R1/m + 2.0*h_c)*R1/m )  ! radial moment of inertia about node 1
     Il = abs(1.0/10.0/m*rho*pi*(R2**5 - R1**5))  
   END IF
   
END SUBROUTINE CylInertia
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE RecInertia(a0, a1, b0, b1, H, rho, Ize, Ixe, Iye)
   REAL(ReKi),                     INTENT    ( IN    )  :: a0  ! Length of side A at node 0
   REAL(ReKi),                     INTENT    ( IN    )  :: a1  ! Length of side A at node 1
   REAL(ReKi),                     INTENT    ( IN    )  :: b0  ! Length of side B at node 0
   REAL(ReKi),                     INTENT    ( IN    )  :: b1  ! Length of side B at node 1
   REAL(ReKi),                     INTENT    ( IN    )  :: H   ! Element height/distance from node 0 to node 1
   REAL(ReKi),                     INTENT    ( IN    )  :: rho ! density of material
   REAL(ReKi),                     INTENT    ( OUT   )  :: Ize ! Moment of inertia about element local z-axis (along member axis)
   REAL(ReKi),                     INTENT    ( OUT   )  :: Ixe ! Moment of inertia about element local x-axis (aligned with sides A)
   REAL(ReKi),                     INTENT    ( OUT   )  :: Iye ! Moment of inertia about element local y-axis (aligned with sides B)
   real(ReKi) :: I1, I2, I3, da, da2, da3, db, db2, db3, a02, a03, b02, b03
   ! All moment of inertia computed about node 0
   da  = a1 - a0
   da2 = da * da
   da3 = da * da2
   db  = b1 - b0
   db2 = db * db
   db3 = db * db2
   a02 = a0 * a0
   a03 = a0 * a02
   b02 = b0 * b0
   b03 = b0 * b02
   I1 = rho*H/12.0 * ( db3*(0.2*a1+0.05*a0) + db2*b0*(0.75*a1+0.25*a0) + db*b02*(a1+0.5*a0) + 0.5*(a1+a0)*b03 )
   I2 = rho*H/12.0 * ( da3*(0.2*b1+0.05*b0) + da2*a0*(0.75*b1+0.25*b0) + da*a02*(b1+0.5*b0) + 0.5*(b1+b0)*a03 )
   I3 = rho*H**3   * ( 0.2*a1*b1 + 0.05*a1*b0 + 0.05*a0*b1 + a0*b0/30.0 )
   Ixe = abs(I1 + I3)
   Iye = abs(I2 + I3)
   Ize = abs(I1 + I2)
END SUBROUTINE RecInertia
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE MarineGrowthPartSegmentCyl(R1, R2, Rmg1, Rmg2, L, rho,  Vinner, Vouter, m_mg, h_c, Ilmg, Irmg)
   
   REAL(ReKi),                     INTENT    ( IN    )  :: R1
   REAL(ReKi),                     INTENT    ( IN    )  :: R2
   REAL(ReKi),                     INTENT    ( IN    )  :: Rmg1
   REAL(ReKi),                     INTENT    ( IN    )  :: Rmg2
   REAL(ReKi),                     INTENT    ( IN    )  :: L
   REAL(ReKi),                     INTENT    ( IN    )  :: rho ! density of material
   REAL(ReKi),                     INTENT    ( OUT   )  :: Vinner   ! volume from inner radius
   REAL(ReKi),                     INTENT    ( OUT   )  :: Vouter   ! volume from outer radius
   REAL(ReKi),                     INTENT    ( OUT   )  :: m_mg   ! mass of marine growth
   REAL(ReKi),                     INTENT    ( OUT   )  :: h_c    ! center of mass offset from first node
   REAL(ReKi),                     INTENT    ( OUT   )  :: Ilmg   ! moment of inertia about axis
   REAL(ReKi),                     INTENT    ( OUT   )  :: Irmg   ! moment of inertia about radial axis from first node
   
   ! Local variables
   
   REAL(ReKi)                         :: cVinner  ! center of volume from inner radius
   REAL(ReKi)                         :: cVouter  ! center of volume from outer radius
   REAL(ReKi)                         :: Ilinner
   REAL(ReKi)                         :: Irinner
   REAL(ReKi)                         :: Ilouter
   REAL(ReKi)                         :: Irouter
      
   ! get V and CV for element
   call CylTaperCalc(R1, R2, L, Vinner, cVinner) 

   ! get V and CV for marine growth displacement
   call CylTaperCalc(Rmg1, Rmg2, L, Vouter, cVouter) 
   
   ! get mass and CV specific to marine growth thickness
   m_mg = (Vouter - Vinner)*rho
   if ( EqualRealNos(m_mg, 0.0_ReKi) ) then
      h_c = 0.0
   else
      h_c = (cVouter*Vouter - Vinner*cVinner)/(Vouter - Vinner)
   end if
   
   ! get two moments of inertia for marine growth as if solid...
   call CylInertia(Rmg1, Rmg2, L, rho, Ilouter, Irouter)  ! inertias for marine growth if solid
   call CylInertia(R1  , R2  , L, rho, Ilinner, Irinner)  ! inertias for element if filled with marine growth

   ! subtract to get moments of inertia of marine growth shell
   Ilmg = Ilouter - Ilinner
   Irmg = Irouter - Irinner

END SUBROUTINE MarineGrowthPartSegmentCyl
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE MarineGrowthPartSegmentRec(a1, a2, b1, b2, amg1, amg2, bmg1, bmg2, L, rho,  Vinner, Vouter, m_mg, h_c, Izemg, Ixemg, Iyemg)
   
   REAL(ReKi),                     INTENT    ( IN    )  :: a1
   REAL(ReKi),                     INTENT    ( IN    )  :: a2
   REAL(ReKi),                     INTENT    ( IN    )  :: b1
   REAL(ReKi),                     INTENT    ( IN    )  :: b2
   REAL(ReKi),                     INTENT    ( IN    )  :: amg1
   REAL(ReKi),                     INTENT    ( IN    )  :: amg2
   REAL(ReKi),                     INTENT    ( IN    )  :: bmg1
   REAL(ReKi),                     INTENT    ( IN    )  :: bmg2
   REAL(ReKi),                     INTENT    ( IN    )  :: L
   REAL(ReKi),                     INTENT    ( IN    )  :: rho ! density of material
   REAL(ReKi),                     INTENT    ( OUT   )  :: Vinner   ! volume from inner radius
   REAL(ReKi),                     INTENT    ( OUT   )  :: Vouter   ! volume from outer radius
   REAL(ReKi),                     INTENT    ( OUT   )  :: m_mg   ! mass of marine growth
   REAL(ReKi),                     INTENT    ( OUT   )  :: h_c    ! center of mass offset from first node
   REAL(ReKi),                     INTENT    ( OUT   )  :: Izemg  ! moment of inertia about axis at first node
   REAL(ReKi),                     INTENT    ( OUT   )  :: Ixemg  ! moment of inertia about element local x-axis at first node
   REAL(ReKi),                     INTENT    ( OUT   )  :: Iyemg  ! moment of inertia about element local y-axis at first node
   
   ! Local variables
   
   REAL(ReKi)                         :: cVinner  ! center of volume from inner radius
   REAL(ReKi)                         :: cVouter  ! center of volume from outer radius
   REAL(ReKi)                         :: Izeinner
   REAL(ReKi)                         :: Ixeinner
   REAL(ReKi)                         :: Iyeinner
   REAL(ReKi)                         :: Izeouter
   REAL(ReKi)                         :: Ixeouter
   REAL(ReKi)                         :: Iyeouter
      
   ! get V and CV for element
   call RecTaperCalc(a1, a2, b1, b2, L, Vinner, cVinner) 

   ! get V and CV for marine growth displacement
   call RecTaperCalc(amg1, amg2, bmg1, bmg2, L, Vouter, cVouter) 
   
   ! get mass and CV specific to marine growth thickness
   m_mg = (Vouter - Vinner)*rho
   if ( EqualRealNos(m_mg, 0.0_ReKi) ) then
      h_c = 0.0
   else
      h_c = (cVouter*Vouter - cVinner*Vinner)/(Vouter - Vinner)
   end if
   
   ! get two moments of inertia for marine growth as if solid...
   call RecInertia(amg1, amg2, bmg1, bmg2, L, rho, Izeouter, Ixeouter, Iyeouter)  ! inertias for marine growth if solid
   call RecInertia(a1  , a2  , b1,   b2,   L, rho, Izeinner, Ixeinner, Iyeinner)  ! inertias for element if filled with marine growth

   ! subtract to get moments of inertia of marine growth shell
   Izemg = Izeouter - Izeinner
   Ixemg = Ixeouter - Ixeinner
   Iyemg = Iyeouter - Iyeinner

END SUBROUTINE MarineGrowthPartSegmentRec
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FloodedBallastPartSegmentCyl(R1, R2, L, rho, V, m, h_c, Il, Ir)
   
   REAL(ReKi),                     INTENT    ( IN    )  :: R1  ! interior radius of element at node point
   REAL(ReKi),                     INTENT    ( IN    )  :: R2  ! interior radius of other end of part-element
   REAL(ReKi),                     INTENT    ( IN    )  :: L   ! distance positive along axis to end of part-element
   REAL(ReKi),                     INTENT    ( IN    )  :: rho ! density of ballast
   REAL(ReKi),                     INTENT    ( OUT   )  :: V   ! volume from inner radius
   REAL(ReKi),                     INTENT    ( OUT   )  :: m   ! mass of material
   REAL(ReKi),                     INTENT    ( OUT   )  :: h_c  ! center of mass offset from first node
   REAL(ReKi),                     INTENT    ( OUT   )  :: Il   ! moment of inertia about axis
   REAL(ReKi),                     INTENT    ( OUT   )  :: Ir   ! moment of inertia about radial axis from first node
   
   ! get V and CV for flooded part of part-element
   call CylTaperCalc(R1, R2, L, V, h_c) 
   m = rho*V
   
   call CylInertia(R1, R2, L, rho, Il, Ir)  ! inertias for filled section

END SUBROUTINE FloodedBallastPartSegmentCyl
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FloodedBallastPartSegmentRec(a1, a2, b1, b2, L, rho, V, m, h_c, Ize, Ixe, Iye)
   
   REAL(ReKi),                     INTENT    ( IN    )  :: a1  ! interior length of side A at node 1
   REAL(ReKi),                     INTENT    ( IN    )  :: a2  ! interior length of side A at node 2
   REAL(ReKi),                     INTENT    ( IN    )  :: b1  ! interior length of side B at node 1
   REAL(ReKi),                     INTENT    ( IN    )  :: b2  ! interior length of side B at node 2
   REAL(ReKi),                     INTENT    ( IN    )  :: L   ! distance positive along axis to end of part-element
   REAL(ReKi),                     INTENT    ( IN    )  :: rho ! density of ballast
   REAL(ReKi),                     INTENT    ( OUT   )  :: V   ! volume from inner radius
   REAL(ReKi),                     INTENT    ( OUT   )  :: m   ! mass of material
   REAL(ReKi),                     INTENT    ( OUT   )  :: h_c    ! center of mass offset from first node
   REAL(ReKi),                     INTENT    ( OUT   )  :: Ize   ! moment of inertia about element local z-axis
   REAL(ReKi),                     INTENT    ( OUT   )  :: Ixe   ! moment of inertia about element local x-axis at first node
   REAL(ReKi),                     INTENT    ( OUT   )  :: Iye   ! moment of inertia about element local y-axis at first node   
   
   ! get V and CV for flooded part of part-element
   call RecTaperCalc(a1, a2, b1, b2, L, V, h_c) 
   m = rho*V
   
   call RecInertia(a1, a2, b1, b2, L, rho, Ize, Ixe, Iye)  ! inertias for filled section

END SUBROUTINE FloodedBallastPartSegmentRec
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE WriteSummaryFile( UnSum, numJoints, numNodes, nodes, numMembers, members, &
                             NOutputs, OutParam, MOutLst, JOutLst, uMesh, yMesh, p, m, errStat, errMsg ) 
                             
   INTEGER,                               INTENT ( IN    )  :: UnSum
   INTEGER,                               INTENT ( IN    )  :: numJoints
   INTEGER,                               INTENT ( IN    )  :: numNodes
   TYPE(Morison_NodeType),   ALLOCATABLE, INTENT ( IN    )  :: nodes(:)  
   INTEGER,                               INTENT ( IN    )  :: numMembers
   TYPE(Morison_MemberType), ALLOCATABLE, INTENT ( IN    )  :: members(:)
   INTEGER,                               INTENT ( IN    )  :: NOutputs
   TYPE(OutParmType),        ALLOCATABLE, INTENT ( IN    )  :: OutParam(:)
   TYPE(Morison_MOutput),    ALLOCATABLE, INTENT ( IN    )  :: MOutLst(:)
   TYPE(Morison_JOutput),    ALLOCATABLE, INTENT ( IN    )  :: JOutLst(:)
   TYPE(MeshType),                        INTENT ( INOUT )  :: uMesh
   TYPE(MeshType),                        INTENT ( INOUT )  :: yMesh
   TYPE(Morison_ParameterType),           INTENT ( IN    ) :: p
   TYPE(Morison_MiscVarType),             INTENT ( IN    ) :: m 
   INTEGER,                               INTENT (   OUT )  :: errStat             ! returns a non-zero value when an error occurs  
   CHARACTER(*),                          INTENT (   OUT )  :: errMsg              ! Error message if errStat /= ErrID_None

   INTEGER                                     :: I, J, II
   LOGICAL                                     :: filledFlag          ! flag indicating if element is filled/flooded
   REAL(ReKi)                                  :: ident(3,3)          ! identity matrix
   REAL(ReKi)                                  :: ExtBuoyancy(6)      ! sum of all external buoyancy forces lumped at (0,0,0)
   REAL(ReKi)                                  :: IntBuoyancy(6)      ! sum of all internal buoyancy forces lumped at (0,0,0)
   REAL(ReKi)                                  :: MG_Wt(6)            ! weight of the marine growth as applied to (0,0,0)
   TYPE(MeshType)                              :: WRP_Mesh            ! mesh representing the WAMIT reference point (0,0,0)
   TYPE(MeshType)                              :: WRP_Mesh_position   ! mesh representing the WAMIT reference point (0,0,0)   (with no displaced position)
   TYPE(MeshMapType)                           :: M_P_2_P             ! Map  Morison Point to  WRP_Mesh point
   REAL(ReKi)                                  :: totalDisplVol       ! total displaced volume of the structure
   REAL(ReKi)                                  :: totalVol            ! total volume of structure
   REAL(ReKi)                                  :: MGvolume            ! volume of the marine growth material
   REAL(ReKi)                                  :: totalMGVol          !
   REAL(ReKi)                                  :: totalFillVol        !
   REAL(ReKi)                                  :: COB(3)              ! center of buoyancy location in global coordinates
   INTEGER                                     :: m1                  ! Indices of the markers which surround the requested output location
   REAL(ReKi)                                  :: s                   ! The linear interpolation factor for the requested location
   REAL(ReKi)                                  :: outloc(3)           ! Position of the requested member output
   real(ReKi)                                  :: pos(3), pos2(3)     ! Position of a node or joint in the MSL inertial system
   INTEGER                                     :: mbrIndx, nodeIndx, c, N
   CHARACTER(ChanLen)                          :: tmpName
   REAL(ReKi)                                  :: totalFillMass, mass_fill, memberVol
   REAL(ReKi)                                  :: totalMGMass
   TYPE(Morison_NodeType)                      :: node1, node2
   real(ReKi)                                  :: ptLoad(6)
   logical                                     :: fillFlag
   type(Morison_MemberType)                    :: mem
   REAL(ReKi)                                  :: Cd1, CdA1, CdB1, Cd2, CdA2, CdB2, Ca1, CaA1, CaB1, Ca2, CaA2, CaB2, Cp1, Cp2, Cb1, Cb2
   REAL(ReKi)                                  :: AxCd1, AxCd2, AxCa1, AxCa2, AxCp1, AxCp2, JAxCd1, JAxCd2, JAxCa1, JAxCa2, JAxCp1, JAxCp2
   real(ReKi)                                  :: F_B(6, numNodes), F_BF(6, numNodes), F_WMG(6, numNodes)
   
   INTEGER                                     :: ErrStat2
   CHARACTER(ErrMsgLen)                        :: ErrMsg2
   CHARACTER(*), PARAMETER                     :: RoutineName = 'WriteSummaryFile'
   
      ! Initialize data
   errStat       = ErrID_None
   errMsg        = ""
   
   IF ( UnSum <= 0 ) RETURN ! can't write to the file (no summary file requested)

      
   ExtBuoyancy   = 0.0
   totalFillMass = 0.0
   totalDisplVol = 0.0
   totalVol      = 0.0
   totalMGVol    = 0.0
   totalFillVol  = 0.0
   totalMGMass   = 0.0
   COB           = 0.0
   F_B           = 0.0
   F_BF          = 0.0
   F_WMG         = 0.0
   
      ! Create identity matrix
   CALL EYE(ident,errStat,errMsg)
   
   do j = 1, numMembers   
      mem = members(j)
      totalVol      = totalVol      + mem%Vouter
      totalMGVol    = totalMGVol    + mem%Vouter - mem%Vinner
      totalDisplVol = totalDisplVol + mem%Vsubmerged
      totalFillVol  = totalFillVol  + mem%Vballast
      
      do i = 1, mem%NElements 
         totalMGMass = totalMGMass + mem%m_mg_l(i)
         totalMGMass = totalMGMass + mem%m_mg_u(i)
      end do
      do i = 1, mem%NElements+1 
         F_B  (:,mem%NodeIndx(i)) = F_B  (:,mem%NodeIndx(i)) + m%memberLoads(j)%F_B  (:,i)
         F_BF (:,mem%NodeIndx(i)) = F_BF (:,mem%NodeIndx(i)) + m%memberLoads(j)%F_BF (:,i)
         F_WMG(:,mem%NodeIndx(i)) = F_WMG(:,mem%NodeIndx(i)) + m%memberLoads(j)%F_WMG(:,i)
      end do
   end do
      
      
      WRITE( UnSum,  '(//)' ) 
      WRITE( UnSum, '(A37)' )        'Strip-Theory Volume Calculations(m^3)'
      WRITE( UnSum, '(A37)' )        '-------------------------------------'
      WRITE( UnSum, '(A27,ES12.5)' ) '  Structure Volume     :   ', totalVol
      WRITE( UnSum, '(A27,ES12.5)' ) '  Submerged Volume     :   ', totalDisplVol
      WRITE( UnSum, '(A27,ES12.5)' ) '  Marine Growth Volume :   ', totalMGVol   
      WRITE( UnSum, '(A27,ES12.5)' ) '  Ballasted Volume     :   ', totalFillVol
      WRITE( UnSum, '(A111)') '              NOTE: Structure, Submerged and Marine Growth volumes are based on members not modelled with WAMIT'
      WRITE( UnSum, '(A149)') '                    Ballasted volume is computed from all members which are marked as filled in the HydroDyn input file, regardless of PropPot flag'
           

   !  Create a point mesh at (0,0,0) so that we can integrate the Morison load contributions to a single point for reporting purposes
      
      CALL MeshCreate( BlankMesh        = WRP_Mesh          &
                     ,IOS               = COMPONENT_INPUT   &
                     ,Nnodes            = 1                 &
                     ,errStat           = errStat2          &
                     ,ErrMess           = errMsg2           &
                     ,Force             = .TRUE.            &
                     ,Moment            = .TRUE.            &
                     )
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         ! Create the node on the mesh
 
      CALL MeshPositionNode (WRP_Mesh                              &
                              , 1                                  &
                              , (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/)   &  
                              , errStat2                           &
                              , errMsg2                            &
                              )
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
      IF ( errStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if
       
      
         ! Create the mesh element
      CALL MeshConstructElement (  WRP_Mesh            &
                                  , ELEMENT_POINT      &                         
                                  , errStat2           &
                                  , errMsg2            &
                                  , 1                  &
                                )
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
      CALL MeshCommit ( WRP_Mesh           &
                      , errStat2           &
                      , errMsg2            )
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
            
         ! we need the translation displacement mesh for loads transfer:
      CALL MeshCopy ( SrcMesh  = WRP_Mesh            &
                    , DestMesh = WRP_Mesh_position   &
                    , CtrlCode = MESH_SIBLING        &
                    , IOS      = COMPONENT_INPUT     &
                    , TranslationDisp = .TRUE.       &
                    , errStat  = errStat2            &
                    , ErrMess  = errMsg2             )  ! automatically sets    DestMesh%RemapFlag = .TRUE.
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
      IF ( errStat >= AbortErrLev ) then
         call cleanup()
         RETURN
      end if
      
      WRP_Mesh_position%TranslationDisp = 0.0  ! bjj: this is actually initialized in the ModMesh module, but I'll do it here anyway.
      
      WRP_Mesh%RemapFlag  = .TRUE.
      
      
         ! Attach the external distributed buoyancy loads to the distributed mesh so they can be transferred to the WRP
      
         ! Because of wave stretching and user-supplied waves, we may have loads above the still water line (SWL) which will be used
         ! in the hydrodynamics for conditions where the wave height is > SWL.  So we now need to check that the vertical position
         ! is <= SWL for this summary file calculation.
      
      
      DO J = 1, yMesh%Nnodes
         
         if ( yMesh%Position(3,J) <= p%WaveField%MSL2SWL ) then  ! need to check relative to MSL2SWL offset because the Mesh Positons are relative to MSL
            
            if (J <= numJoints) then
               ptLoad = F_B(:,J) + m%F_B_end(:,J)
            else
               ptLoad = F_B(:,J)
            end if
            yMesh%Force(:,J)   = ptLoad(1:3)
            yMesh%Moment(:,J)  = ptLoad(4:6)
         else
            yMesh%Force(:,J)   = 0.0
            yMesh%Moment(:,J)  = 0.0
         end if              ! <= still water line check
      END DO ! DO J
      
   
         ! Transfer the loads from the distributed mesh to the (0,0,0) point mesh
         
      CALL MeshMapCreate           ( yMesh, WRP_Mesh, M_P_2_P, errStat2, errMsg2               )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF ( errStat >= AbortErrLev ) then
            call cleanup()
            RETURN
         end if
      
      CALL Transfer_Point_to_Point( yMesh, WRP_Mesh, M_P_2_P, errStat2, errMsg2, uMesh, WRP_Mesh_position );  call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      ExtBuoyancy(1:3) = WRP_Mesh%Force (:,1)
      ExtBuoyancy(4:6) = WRP_Mesh%Moment(:,1)
    
  
         ! Write the buoyancy table headers and the external results
   
      WRITE( UnSum,  '(//)' ) 
      WRITE( UnSum, '(A51)' ) 'Total Buoyancy loads summed about ( 0.0, 0.0, 0.0 )'
      WRITE( UnSum, '(A51)' ) '---------------------------------------------------'
      WRITE( UnSum, '(18x,6(2X,A20))' ) ' BuoyFxi ', ' BuoyFyi ', ' BuoyFzi ', ' BuoyMxi ', ' BuoyMyi ', ' BuoyMzi '
      WRITE( UnSum, '(18x,6(2X,A20))' ) '   (N)   ', '   (N)   ', '   (N)   ', '  (N-m)  ', '  (N-m)  ', '  (N-m)  '
      WRITE( UnSum, '(A18,6(2X,ES20.6))') 'External:        ', ExtBuoyancy(1), ExtBuoyancy(2), ExtBuoyancy(3), ExtBuoyancy(4), ExtBuoyancy(5), ExtBuoyancy(6)
      
      
         ! Now compute internal Buoyancy
         
      DO J = 1, yMesh%Nnodes
         
         if (J <= numJoints) then
            ptLoad = F_BF(:,J) + m%F_BF_end(:,J)
         else
            ptLoad = F_BF(:,J)
         end if
         yMesh%Force(:,J)   = ptLoad(1:3)
         yMesh%Moment(:,J)  = ptLoad(4:6)
         
      END DO ! DO J
       
      IntBuoyancy = 0.0
      CALL Transfer_Point_to_Point( yMesh, WRP_Mesh, M_P_2_P, errStat2, errMsg2, uMesh, WRP_Mesh_position );  call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IntBuoyancy(1:3) = WRP_Mesh%Force(:,1)
      IntBuoyancy(4:6) = WRP_Mesh%Moment(:,1)
      
     
      WRITE( UnSum, '(A18,6(2X,ES20.6))') 'Internal:        ', IntBuoyancy(1), IntBuoyancy(2), IntBuoyancy(3), IntBuoyancy(4), IntBuoyancy(5), IntBuoyancy(6)
      IntBuoyancy = IntBuoyancy + ExtBuoyancy
      WRITE( UnSum, '(A18,6(2X,ES20.6))') 'Total   :        ', IntBuoyancy(1), IntBuoyancy(2), IntBuoyancy(3), IntBuoyancy(4), IntBuoyancy(5), IntBuoyancy(6)
      WRITE( UnSum, '(A81)') '              NOTE: External buoyancy is based on members not modelled with WAMIT'
      WRITE( UnSum, '(A150)') '                    Internal buoyancy is computed from all members which are marked as filled in the HydroDyn input file, regardless of PropPot flag'
      WRITE( UnSum, '(A88)') '                    Total buoyancy does not include WAMIT-modelled buoyancy contribution'
      
      
      
      !   ! Now compute marine growth weight at the WRP
         
      DO J = 1, yMesh%Nnodes
         if (J <= numJoints) then
            yMesh%Force(:,J)   = F_WMG(1:3,J) + p%F_WMG_End(:,J)
         else
            yMesh%Force(:,J)   = F_WMG(1:3,J)
         end if
         
         yMesh%Moment(:,J)  = F_WMG(4:6,J)
         
      END DO ! DO J
         
      MG_Wt = 0.0
      CALL Transfer_Point_to_Point( yMesh, WRP_Mesh, M_P_2_P, errStat2, errMsg2, uMesh, WRP_Mesh_position );  call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      MG_Wt(1:3) = WRP_Mesh%Force(:,1)
      MG_Wt(4:6) = WRP_Mesh%Moment(:,1)
      !
           
      WRITE( UnSum,  '(//)' ) 
      WRITE( UnSum, '(A36)' ) 'Weight loads about ( 0.0, 0.0, 0.0 )'
      WRITE( UnSum, '(A36)' ) '------------------------------------'
      WRITE( UnSum, '(18x,6(2X,A20))' ) '  MGFxi  ', '  MGFyi  ', '  MGFzi  ', '  MGMxi  ', '  MGMyi  ', '  MGMzi  '
      WRITE( UnSum, '(18x,6(2X,A20))' ) '   (N)   ', '   (N)   ', '   (N)   ', '  (N-m)  ', '  (N-m)  ', '  (N-m)  '
      WRITE( UnSum, '(A18,6(2X,ES20.6))') 'Marine Growth:   ', MG_Wt(1), MG_Wt(2), MG_Wt(3), MG_Wt(4), MG_Wt(5), MG_Wt(6)

      
      !
      !   ! Write the header for this section
      WRITE( UnSum,  '(//)' ) 
      WRITE( UnSum,  '(A14,I4,A44)' ) 'Nodes (first [',numJoints,'] are joints, remainder are internal nodes)'
      WRITE( UnSum,  '(/)' ) 
      WRITE( UnSum, '(1X,A5,31(2X,A10))' ) &
         '  i  '     , '  MbrIndx ', '   Nxi    ', '   Nyi    ', '   Nzi    ', '  Shape   ', '     R    ', '    SA    ', '    SB    ', '    t     ', '   tMG    ', '  MGDens  ', &
         ' PropPot  ', 'FilledFlag', 'FilledMass', '    Cd    ', '    CdA   ', '    CdB   ', '    Ca    ', '    CaA   ', '    CaB   ', '    Cp    ', '    Cb    ', &
         '   AxCd   ', '   AxCa   ', '   AxCp   ', '   JAxCd  ', '   JAxCa  ', '   JAxCp  ', ' JAxFDMod ', 'JAxVnCOff ', 'JAxFDLoFSc'
      WRITE( UnSum, '(1X,A5,31(2X,A10))' ) &
         ' (-) '     , '    (-)   ', '    (m)   ', '    (m)   ', '    (m)   ', '    (-)   ', '    (m)   ', '    (m)   ', '    (m)   ', '   (m)    ', '   (m)    ', ' (kg/m^3) ', &
         '    (-)   ', '    (-)   ', '   (kg)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '   (-)    ', '   (-)    ', &
         '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   '
      
         ! Write the node data
      do I = 1,numJoints   
         ! need to add MSL2SWL offset from this because the Positons are relative to SWL, but we should report them relative to MSL here
         pos = nodes(i)%Position
         pos(3) = pos(3) + p%WaveField%MSL2SWL
         write( UnSum, '(1X,I5,(2X,A10),3(2X,F10.4),5(2X,A10),2(2X,ES10.3),14(2X,A10),3(2X,ES10.3),1X,I10,2(2X,ES10.3))' ) &
            i,         '     -    ', pos,                                      '     -    ', '     -    ', '     -    ', '     -    ', '     -    ', nodes(i)%tMG, nodes(i)%MGdensity, &
         '     -    ', '     -    ', '     -    ', '     -    ', '     -    ', '     -    ', '     -    ', '     -    ', '     -    ', '     -    ', '     -    ', &
         '     -    ', '     -    ', '     -    ',nodes(i)%JAxCd,nodes(i)%JAxCa,nodes(i)%JAxCp,nodes(i)%JAxFDMod,nodes(i)%JAxVnCOff,nodes(i)%JAxFDLoFSc
      end do
      c = numJoints
      do j= 1, numMembers
         do i = 2, members(j)%NElements
            c = c + 1
            if (members(j)%l_fill - members(j)%dl*(i-1) > 0.0) then
               fillFlag = .true.
            else
               fillFlag = .false.
            end if
            ! need to add MSL2SWL offset from this because the Positons are relative to SWL, but we should report them relative to MSL here
            pos = nodes(c)%Position
            pos(3) = pos(3) + p%WaveField%MSL2SWL
            if (members(j)%flipped) then
               II=members(j)%NElements+2-I
            else
               II=I
            endif
            if ( members(j)%MSecGeom == MSecGeom_Cyl ) then
               write( UnSum, '(1X,I5,(2X,I10),3(2X,F10.4),(2X,A10),(2X,ES10.3),2(2X,A10),3(2X,ES10.3),2(6X,L6),2(2X,ES10.3),2(2X,A10),(2X,ES10.3),2(2X,A10),5(2X,ES10.3),6(7x,A5))' ) &
                  c, members(j)%MemberID, pos, '    C     ', members(j)%R(ii),      '    -     ',      '    -     ',         members(j)%R(ii)-members(j)%Rin(ii), members(j)%tMG(ii), members(j)%MGdensity(ii), &
                  members(j)%PropPot,  fillFlag,  members(j)%m_fb_u(ii)+members(j)%m_fb_l(ii), members(j)%Cd(ii),       '    -     ',       '    -     ', members(j)%Ca(ii), &
                  '    -     ',       '    -     ',       members(j)%Cp(ii), members(j)%Cb(ii), members(j)%AxCd(ii), members(j)%AxCa(ii), members(j)%AxCp(ii), &
                  '  -  ', '  -  ', '  -  ', '  -  ', '  -  ', '  -  '
            else if ( members(j)%MSecGeom == MSecGeom_Rec ) then
               write( UnSum, '(1X,I5,(2X,I10),3(2X,F10.4),(2X,A10),(2X,A10),5(2X,ES10.3),2(6X,L6),(2X,ES10.3),(2X,A10),2(2X,ES10.3),(2X,A10),7(2X,ES10.3),6(7x,A5))' ) &
                  c, members(j)%MemberID, pos, '    R     ',     '    -     ', members(j)%Sa(ii), members(j)%Sb(ii), 0.5*(members(j)%Sa(ii)-members(j)%Sain(ii)), members(j)%tMG(ii), members(j)%MGdensity(ii), &
                  members(j)%PropPot,  fillFlag,  members(j)%m_fb_u(ii)+members(j)%m_fb_l(ii),      '    -     ', members(j)%CdA(ii), members(j)%CdB(ii),      '    -     ', &
                  members(j)%CaA(ii), members(j)%CaB(ii), members(j)%Cp(ii), members(j)%Cb(ii), members(j)%AxCd(ii), members(j)%AxCa(ii), members(j)%AxCp(ii), &
                  '  -  ', '  -  ', '  -  ', '  -  ', '  -  ', '  -  '
            end if
         end do
      end do
      
      write( UnSum,  '(//)' ) 
      write( UnSum,  '(A8)' ) 'Members'
      write( UnSum,  '(/)' ) 
      write( UnSum, '(1X,A8,2X,A6,2X,A6,45(2X,A12))' ) &
         'MemberID'    , 'joint1'      , 'joint2'      , '  Length  '  , '   NElem    ', '   Volume   ', '  MGVolume  ', '      R1    ', '     SA1    ', '     SB1    ', '     t1     ', &
         '      R2    ', '     SA2    ', '     SB2    ', '     t2     ', ' PropPot  '  , 'FilledFlag'  , 'FillDensity' , '  FillFSLoc ', '  FillMass  ', '     Cd1    ', '    CdA1    ', '    CdB1    ', &
         '     Ca1    ', '    CaA1    ', '    CaB1   ' , '     Cp1    ', '     Cb1    ', '    AxCd1   ', '    AxCa1   ', '    AxCp1   ', '   JAxCd1   ', '  JAxCa1    ', '   JAxCp1   ', &
         '     Cd2    ', '    CdA2    ', '    CdB2    ', '     Ca2    ', '    CaA2    ', '    CaB2    ', '     Cp2    ', '     Cb2    ', '    AxCd2   ', '    AxCa2   ', '    AxCp2   ', &
         '   JAxCd2   ', '   JAxCa2   ', '   JAxCp2   '
      write( UnSum, '(1X,A8,2X,A6,2X,A6,45(2X,A12))' ) &
         '  (-)   '    , ' (-)  '      , ' (-)  '      , '   (m)    '  , '    (-)     ', '   (m^3)    ', '   (m^3)    ', '      (m)   ', '     (m)    ', '     (m)    ', '     (m)    ', &
         '      (m)   ', '     (m)    ', '     (m)    ', '     (m)    ', '   (-)    '  , '   (-)    '  , ' (kg/m^3)  ' , '     (-)    ', '    (kg)    ', '     (-)    ', '     (-)    ', '     (-)    ', &
         '     (-)    ', '     (-)    ', '     (-)   ' , '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', &
         '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', &
         '     (-)    ', '     (-)    ', '     (-)    '
      
      do i = 1,numMembers
         N = members(i)%NElements

         IF (members(i)%PropPot) THEN
            MGvolume    = 0.0
            memberVol   = 0.0
         ELSE
            memberVol   = members(i)%Vouter
            MGvolume    = members(i)%Vouter - members(i)%Vinner
         END IF
   
         IF ( members(i)%l_fill > 0.0 ) THEN          
            filledFlag = .TRUE.         
            mass_fill     = members(i)%FillDens*members(i)%Vballast
         ELSE
            mass_fill  = 0.0
            filledFlag = .FALSE.
         END IF
     
         IF ( members(i)%MSecGeom == MSecGeom_Cyl ) THEN
            Cd1   = members(i)%Cd(1)
            Cd2   = members(i)%Cd(N+1)
            Ca1   = members(i)%Ca(1)
            Ca2   = members(i)%Ca(N+1)
         ELSE IF ( members(i)%MSecGeom == MSecGeom_Rec ) THEN
            CdA1   = members(i)%CdA(1)
            CdA2   = members(i)%CdA(N+1)
            CaA1   = members(i)%CaA(1)
            CaA2   = members(i)%CaA(N+1)
            CdB1   = members(i)%CdB(1)
            CdB2   = members(i)%CdB(N+1)
            CaB1   = members(i)%CaB(1)
            CaB2   = members(i)%CaB(N+1)
         END IF

         Cp1   = members(i)%Cp(1)
         Cp2   = members(i)%Cp(N+1)
         AxCd1 = members(i)%AxCd(1)
         AxCd2 = members(i)%AxCd(N+1)
         AxCa1 = members(i)%AxCa(1)
         AxCa2 = members(i)%AxCa(N+1)
         AxCp1 = members(i)%AxCp(1)
         AxCp2 = members(i)%AxCp(N+1)
         Cb1   = members(i)%Cb(1)
         Cb2   = members(i)%Cb(N+1)
     
         JAxCd1 = nodes(members(i)%NodeIndx(1  ))%JAxCd
         JAxCd2 = nodes(members(i)%NodeIndx(1+N))%JAxCd
         JAxCa1 = nodes(members(i)%NodeIndx(1  ))%JAxCa
         JAxCa2 = nodes(members(i)%NodeIndx(1+N))%JAxCa
         JAxCp1 = nodes(members(i)%NodeIndx(1  ))%JAxCp
         JAxCp2 = nodes(members(i)%NodeIndx(1+N))%JAxCp
       
         IF ( members(i)%MSecGeom == MSecGeom_Cyl ) THEN
            write( UnSum, '(1X,I8,2X,I6,2X,I6,2X,ES12.5,2X,I12, 3(2X,ES12.5),2(2X,A12),2(2X,ES12.5),2(2X,A12),(2X,ES12.5),2(2X,L12), 4(2X,ES12.5),2(2X,A12),(2X,ES12.5),2(2X,A12),9(2X,ES12.5),2(2X,A12),(2X,ES12.5),2(2X,A12),8(2X,ES12.5))' )  members(i)%MemberID, &
                          members(i)%NodeIndx(1), members(i)%NodeIndx(N+1), members(i)%RefLength, N, &
                          memberVol, MGvolume, members(i)%Rmg(1), '     -      ', '     -      ', members(i)%Rmg(1)-members(i)%Rin(1), &
                          members(i)%Rmg(N+1), '     -      ', '     -      ', members(i)%Rmg(N+1)-members(i)%Rin(N+1),  &
                          members(i)%PropPot, filledFlag, members(i)%FillDens, members(i)%FillFSLoc, &
                          mass_fill, Cd1, '     -      ', '     -      ', Ca1, '     -      ', '     -      ', Cp1, Cb1, AxCd1, AxCa1, AxCp1, JAxCd1, JAxCa1, JAxCp1, &
                          Cd2, '     -      ', '     -      ', Ca2, '     -      ', '     -      ', Cp2, Cb2, AxCd2, AxCa2, AxCp2, JAxCd2, JAxCa2, JAxCp2
         ELSE IF ( members(i)%MSecGeom == MSecGeom_Rec ) THEN
            write( UnSum, '(1X,I8,2X,I6,2X,I6,2X,ES12.5,2X,I12, 2(2X,ES12.5),(2X,A12),3(2X,ES12.5),(2X,A12),3(2X,ES12.5),2(2X,L12), 3(2X,ES12.5),(2X,A12),2(2X,ES12.5),(2X,A12),10(2X,ES12.5),(2X,A12),2(2X,ES12.5),(2X,A12),10(2X,ES12.5))' )  members(i)%MemberID, &
                          members(i)%NodeIndx(1), members(i)%NodeIndx(N+1), members(i)%RefLength, N, &
                          memberVol, MGvolume, '     -      ', members(i)%SaMG(1), members(i)%SbMG(1), 0.5*(members(i)%SaMG(1)-members(i)%SaIn(1)), &
                          '     -      ', members(i)%SaMG(N+1), members(i)%SbMG(N+1), 0.5*(members(i)%SaMG(N+1)-members(i)%SaIn(N+1)),  &
                          members(i)%PropPot, filledFlag, members(i)%FillDens, members(i)%FillFSLoc, &
                          mass_fill, '     -      ', CdA1, CdB1, '     -      ', CaA1, CaB1, Cp1, Cb1, AxCd1, AxCa1, AxCp1, JAxCd1, JAxCa1, JAxCp1, &
                          '     -      ', CdA2, CdB2, '     -      ', CaA2, CaB2, Cp2, Cb2, AxCd2, AxCa2, AxCp2, JAxCd2, JAxCa2, JAxCp2
         END IF
      
      end do   ! i = 1,numMembers
               
      
      WRITE( UnSum,  '(//)' ) 
      WRITE( UnSum,  '(A24)' ) 'Requested Member Outputs'
      WRITE( UnSum,  '(/)' ) 
      WRITE( UnSum, '(1X,A10,11(2X,A10))' ) '  Label   ', '    Xi    ',  '    Yi    ', '    Zi    ', ' MemberID ', ' StartXi  ',  ' StartYi  ', ' StartZi  ', '  EndXi   ', '  EndYi   ', '  EndZi   ', '   Loc    '
      WRITE( UnSum, '(1X,A10,11(2X,A10))' ) '   (-)    ', '    (m)   ',  '    (m)   ', '    (m)   ', '   (-)    ', '   (m)    ',  '   (m)    ', '   (m)    ', '   (m)    ', '   (m)    ', '   (m)    ', '   (-)    '
      
      
      DO I = 1,NOutputs

         tmpName =  OutParam(I)%Name
         IF (OutParam(I)%SignM == -1 ) tmpName = tmpName(2:10)
               
         IF ( ( INDEX( 'mM', tmpName(1:1) ) > 0 ) .AND. (OutParam(I)%Units /= 'INVALID' ) ) THEN
               !Get Member index and Node index
            read (tmpName(2:2),*) mbrIndx
            read (tmpName(4:4),*) nodeIndx
            
             
           
            s  = MOutLst(mbrIndx)%NodeLocs(nodeIndx)
            ! Find the member starting and ending node locations
               ! The member output is computed as a linear interpolation of the nearest two markers
            mem   = members(MOutLst(mbrIndx)%MemberIDIndx)
            node1 = nodes(mem%NodeIndx(1))
            node2 = nodes(mem%NodeIndx(mem%NElements+1))
            ! need to add MSL2SWL offset from this because the Positons are relative to SWL, but we should report them relative to MSL here
            pos = node1%Position
            pos(3) = pos(3) + p%WaveField%MSL2SWL
            pos2 = node2%Position
            pos2(3) = pos2(3) + p%WaveField%MSL2SWL
            outLoc    = pos*(1-s) + pos2*s
            WRITE( UnSum, '(1X,A10,3(2x,F10.4),2x,I10,7(2x,F10.4))' ) OutParam(I)%Name, outLoc,  MOutLst(mbrIndx)%MemberID, pos,pos2, s
         END IF
  
      END DO
      
      
      WRITE( UnSum,  '(//)' ) 
      WRITE( UnSum,  '(A24)' ) 'Requested Joint Outputs'
      WRITE( UnSum,  '(/)' ) 
      WRITE( UnSum, '(1X,A10,5(2X,A10))' ) '  Label   ', '    Xi    ',  '    Yi    ', '    Zi    ', 'InpJointID'
      WRITE( UnSum, '(1X,A10,5(2X,A10))' ) '   (-)    ', '    (m)   ',  '    (m)   ', '    (m)   ', '   (-)    '
      
      
      DO I = 1,NOutputs
         tmpName =  OutParam(I)%Name
         IF (OutParam(I)%SignM == -1 ) tmpName = tmpName(2:10)
               
         IF ( ( INDEX( 'jJ', tmpName(1:1) ) > 0 ) .AND. (OutParam(I)%Units /= 'INVALID') ) THEN
            
               !Get Member index and Node index
            read (tmpName(2:2),*) nodeIndx
            m1 = JOutLst(nodeIndx)%JointIDIndx 
            ! need to add MSL2SWL offset from this because the Positons are relative to SWL, but we should report them relative to MSL here
            pos = nodes(m1)%Position
            pos(3) = pos(3) + p%WaveField%MSL2SWL
            WRITE( UnSum, '(1X,A10,3(2x,F10.4),2x,I10)' ) OutParam(I)%Name, pos, JOutLst(nodeIndx)%JointID
            
         END IF
         
         
      END DO
      
   
   call cleanup()
   
contains
!...................................
   subroutine cleanup()
      call MeshDestroy(WRP_Mesh, ErrStat2, ErrMsg2)
      call MeshDestroy(WRP_Mesh_position, ErrStat2, ErrMsg2)
      call MeshMapDestroy(M_P_2_P, ErrStat2, ErrMsg2)
      
      call Morison_DestroyNodeType(node1, ErrStat2, ErrMsg2)
      call Morison_DestroyNodeType(node2, ErrStat2, ErrMsg2)
      call Morison_DestroyMemberType(mem, ErrStat2, ErrMsg2)
   end subroutine cleanup
END SUBROUTINE WriteSummaryFile

!----------------------------------------------------------------------------------------------------------------------------------
subroutine Morison_GenerateSimulationNodes( MSL2SWL, numJoints, inpJoints, numMembers, inpMembers, numNodes, nodes, errStat, errMsg )
   ! This subdivides a Morison member according to its maximum desired 
   ! element length (MDivSize), allocating the member's arrays, and
   ! adding resuling new nodes to the master node array.
   real(ReKi),                          intent (in   ) :: MSL2SWL              ! mean sea level To still water level offset value
   integer,                             intent (in   ) :: numJoints            ! number of joints in the input file
   type(Morison_JointType),             intent (in   ) :: inpJoints(:)         ! array of input joint data structures
   integer,                             intent (in   ) :: numMembers           ! number of members specified in the inputs
   type(Morison_MemberInputType),       intent (inout) :: inpMembers(:)        ! array of input member data structures                                    
   integer,                             intent (inout) :: numNodes             ! the total number of nodes being used for the simulation model
   type(Morison_NodeType), allocatable, intent (inout) :: nodes(:)             ! the array of simulation nodes
   integer,                             intent (  out) :: errStat              ! returns a non-zero value when an error occurs  
   character(*),                        intent (  out) :: errMsg               ! Error message if errStat /= ErrID_None
   

   integer                                     :: numDiv, maxNodes
   integer                                     :: i, j
   real(ReKi)                                  :: s  ! interpolation factor
   real(ReKi)                                  :: memLength  ! member length
   integer                                     :: j1, j2                ! generic integer for counting
   INTEGER(IntKi)                              :: errStat2    ! Error status of the operation (occurs after initial error)
   CHARACTER(errMsgLen)                        :: errMsg2     ! Error message if errStat2 /= ErrID_None

      ! Initialize errStat
         
   errStat = ErrID_None         
   errMsg  = "" 
   
      ! Initialize quantities
   maxNodes         = numJoints
   
   ! Determine maximum nodes in simulation mesh due to internal member subdivision
   do i = 1,numMembers

      j1 = inpMembers(I)%MJointID1Indx
      j2 = inpMembers(I)%MJointID2Indx
      call GetDistance(inpJoints(j1)%Position, inpJoints(j2)%Position, memLength)
      if ( EqualRealNos(memLength, 0.0_ReKi) ) then
         errMsg  = ' Input file member with ID: '//trim(num2lstr(inpMembers(i)%MemberID))//' must have length greater than zero.'
         errStat = ErrID_Fatal
         return
      end if
      numDiv = CEILING( memLength / inpMembers(i)%MDivSize  )
      ! set number of elements in member and element size
      inpMembers(i)%NElements = numDiv
      inpMembers(i)%dl        = memLength/numDiv
      inpMembers(i)%refLength = memLength
      maxNodes = maxNodes + numDiv - 1
      
   end do 
   
   ! Allocate nodes array
   allocate ( nodes(maxNodes), STAT = errStat2 )
      if ( errStat2 /= 0 ) then
         errMsg  = ' Error allocating space for Nodes array for Morison Module.'
         errStat = ErrID_Fatal
         return
      end if
    
   ! Loop over the input file joints and add there positions as the node positions at the beginning of the nodes array
   do i = 1, numJoints
      nodes(i)%Position(1:2) = inpJoints(i)%Position(1:2)
      nodes(i)%Position(3)   = inpJoints(i)%Position(3)   - MSL2SWL  ! Correct the Z-coordinate based on the mean sea level To still water level offset value
   end do
   
   numNodes = numJoints
   ! Now loop over the input file members and create necessary internal nodes and add them to the nodes array
   ! Also augment the input file members data with the new discretization information
   do i = 1,numMembers
      call AllocAry(inpMembers(i)%NodeIndx, inpMembers(i)%NElements+1, 'inpMembers(i)%NodeIndx', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'Morison_GenerateSimulationNodes')
      if (errStat >= AbortErrLev) return
      
      j1 = inpMembers(i)%MJointID1Indx
      j2 = inpMembers(i)%MJointID2Indx
      numDiv = inpMembers(i)%NElements 
      inpMembers(i)%NodeIndx(1) = j1
      inpMembers(i)%NodeIndx(1+numDiv) = j2
      ! If the requested division size is less then the member length, create nodes along the member         
      if (numDiv > 1 ) THEN
         ! loop through the new node locations along the member and add a node at each
         do j = 1, numDiv-1
            numNodes = numNodes + 1
            s = real(j, ReKi) / real(numDiv, ReKi)
            nodes(numNodes)%Position =  inpJoints(j1)%Position*(1-s) + inpJoints(j2)%Position*s
            nodes(numNodes)%Position(3) = nodes(numNodes)%Position(3) - MSL2SWL  ! Correct the Z-coordinate based on the mean sea level To still water level offset value
            inpMembers(i)%NodeIndx(j+1) = numNodes
         end do 
      end if        
    end do     
end subroutine Morison_GenerateSimulationNodes
   
   


!====================================================================================================
SUBROUTINE SetDepthBasedCoefs_Cyl( z, tMG, NCoefDpth, CoefDpths, Cd, Ca, Cp, AxCd, AxCa, AxCp, Cb )
   
   REAL(ReKi), INTENT (IN   )             :: z ! Z location relative to MSL inertial system
   REAL(ReKi), INTENT (IN   )             :: tMG
   INTEGER,    INTENT (IN   )             :: NCoefDpth
   TYPE(Morison_CoefDpthsCyl), INTENT (IN   ):: CoefDpths(:)
   REAL(ReKi), INTENT (  OUT)             :: Cd
   REAL(ReKi), INTENT (  OUT)             :: Ca
   REAL(ReKi), INTENT (  OUT)             :: Cp
   REAL(ReKi), INTENT (  OUT)             :: AxCd
   REAL(ReKi), INTENT (  OUT)             :: AxCa
   REAL(ReKi), INTENT (  OUT)             :: AxCp
   REAL(ReKi), INTENT (  OUT)             :: Cb
   
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
   if ( tMG > 0.0_ReKi ) then    
      Cd     = CoefDpths(indx1)%DpthCdMG*(1-s)   + CoefDpths(indx2)%DpthCdMG*s
      Ca     = CoefDpths(indx1)%DpthCaMG*(1-s)   + CoefDpths(indx2)%DpthCaMG*s
      Cp     = CoefDpths(indx1)%DpthCpMG*(1-s)   + CoefDpths(indx2)%DpthCpMG*s
      AxCd   = CoefDpths(indx1)%DpthAxCdMG*(1-s) + CoefDpths(indx2)%DpthAxCdMG*s
      AxCa   = CoefDpths(indx1)%DpthAxCaMG*(1-s) + CoefDpths(indx2)%DpthAxCaMG*s
      AxCp   = CoefDpths(indx1)%DpthAxCpMG*(1-s) + CoefDpths(indx2)%DpthAxCpMG*s
      Cb     = CoefDpths(indx1)%DpthCbMG*(1-s)   + CoefDpths(indx2)%DpthCbMG*s
   else
      Cd     = CoefDpths(indx1)%DpthCd*(1-s)     + CoefDpths(indx2)%DpthCd*s
      Ca     = CoefDpths(indx1)%DpthCa*(1-s)     + CoefDpths(indx2)%DpthCa*s
      Cp     = CoefDpths(indx1)%DpthCp*(1-s)     + CoefDpths(indx2)%DpthCp*s
      AxCd   = CoefDpths(indx1)%DpthAxCd*(1-s)   + CoefDpths(indx2)%DpthAxCd*s
      AxCa   = CoefDpths(indx1)%DpthAxCa*(1-s)   + CoefDpths(indx2)%DpthAxCa*s
      AxCp   = CoefDpths(indx1)%DpthAxCp*(1-s)   + CoefDpths(indx2)%DpthAxCp*s
      Cb     = CoefDpths(indx1)%DpthCb*(1-s)     + CoefDpths(indx2)%DpthCb*s
   end if
   

END SUBROUTINE SetDepthBasedCoefs_Cyl


!====================================================================================================
SUBROUTINE SetDepthBasedCoefs_Rec( z, tMG, NCoefDpth, CoefDpths, CdA, CdB, CaA, CaB, Cp, AxCd, AxCa, AxCp, Cb )
   
   REAL(ReKi), INTENT (IN   )             :: z ! Z location relative to MSL inertial system
   REAL(ReKi), INTENT (IN   )             :: tMG
   INTEGER,    INTENT (IN   )             :: NCoefDpth
   TYPE(Morison_CoefDpthsRec), INTENT (IN   ):: CoefDpths(:)
   REAL(ReKi), INTENT (  OUT)             :: CdA
   REAL(ReKi), INTENT (  OUT)             :: CdB
   REAL(ReKi), INTENT (  OUT)             :: CaA
   REAL(ReKi), INTENT (  OUT)             :: CaB
   REAL(ReKi), INTENT (  OUT)             :: Cp
   REAL(ReKi), INTENT (  OUT)             :: AxCd
   REAL(ReKi), INTENT (  OUT)             :: AxCa
   REAL(ReKi), INTENT (  OUT)             :: AxCp
   REAL(ReKi), INTENT (  OUT)             :: Cb
   
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
   if ( tMG > 0.0_ReKi ) then    
      CdA    = CoefDpths(indx1)%DpthCdAMG*(1-s)  + CoefDpths(indx2)%DpthCdAMG*s
      CdB    = CoefDpths(indx1)%DpthCdBMG*(1-s)  + CoefDpths(indx2)%DpthCdBMG*s
      CaA    = CoefDpths(indx1)%DpthCaAMG*(1-s)  + CoefDpths(indx2)%DpthCaAMG*s
      CaB    = CoefDpths(indx1)%DpthCaBMG*(1-s)  + CoefDpths(indx2)%DpthCaBMG*s
      Cp     = CoefDpths(indx1)%DpthCpMG*(1-s)   + CoefDpths(indx2)%DpthCpMG*s
      AxCd   = CoefDpths(indx1)%DpthAxCdMG*(1-s) + CoefDpths(indx2)%DpthAxCdMG*s
      AxCa   = CoefDpths(indx1)%DpthAxCaMG*(1-s) + CoefDpths(indx2)%DpthAxCaMG*s
      AxCp   = CoefDpths(indx1)%DpthAxCpMG*(1-s) + CoefDpths(indx2)%DpthAxCpMG*s
      Cb     = CoefDpths(indx1)%DpthCbMG*(1-s)   + CoefDpths(indx2)%DpthCbMG*s
   else
      CdA    = CoefDpths(indx1)%DpthCdA*(1-s)    + CoefDpths(indx2)%DpthCdA*s
      CdB    = CoefDpths(indx1)%DpthCdB*(1-s)    + CoefDpths(indx2)%DpthCdB*s
      CaA    = CoefDpths(indx1)%DpthCaA*(1-s)    + CoefDpths(indx2)%DpthCaA*s
      CaB    = CoefDpths(indx1)%DpthCaB*(1-s)    + CoefDpths(indx2)%DpthCaB*s
      Cp     = CoefDpths(indx1)%DpthCp*(1-s)     + CoefDpths(indx2)%DpthCp*s
      AxCd   = CoefDpths(indx1)%DpthAxCd*(1-s)   + CoefDpths(indx2)%DpthAxCd*s
      AxCa   = CoefDpths(indx1)%DpthAxCa*(1-s)   + CoefDpths(indx2)%DpthAxCa*s
      AxCp   = CoefDpths(indx1)%DpthAxCp*(1-s)   + CoefDpths(indx2)%DpthAxCp*s
      Cb     = CoefDpths(indx1)%DpthCb*(1-s)     + CoefDpths(indx2)%DpthCb*s
   end if
   

END SUBROUTINE SetDepthBasedCoefs_Rec



!====================================================================================================
SUBROUTINE SetExternalHydroCoefs_Cyl(  MSL2SWL, MCoefMod, MmbrCoefIDIndx, SimplCd, SimplCdMG, SimplCa, SimplCaMG, SimplCp, &
                                   SimplCpMG, SimplAxCd, SimplAxCdMG, SimplAxCa, SimplAxCaMG, SimplAxCp, SimplAxCpMG, SimplCb, SimplCbMG, SimplMCF, CoefMembers,    &
                                   NCoefDpth, CoefDpths, nodes, member )   
!     This private subroutine generates the Cd, Ca, Cp, Cb, CdMG, CaMG, CpMG, and CbMG coefs for the member based on
!     the input data.  
!---------------------------------------------------------------------------------------------------- 
   real(ReKi),                             intent(in   )  :: MSL2SWL
   integer(IntKi),                         intent(in   )  :: MCoefMod
   integer(IntKi),                         intent(in   )  :: MmbrCoefIDIndx
   real(ReKi),                             intent(in   )  :: SimplCd 
   real(ReKi),                             intent(in   )  :: SimplCdMG
   real(ReKi),                             intent(in   )  :: SimplCa
   real(ReKi),                             intent(in   )  :: SimplCaMG 
   real(ReKi),                             intent(in   )  :: SimplCp
   real(ReKi),                             intent(in   )  :: SimplCpMG 
   real(ReKi),                             intent(in   )  :: SimplAxCd
   real(ReKi),                             intent(in   )  :: SimplAxCdMG 
   real(ReKi),                             intent(in   )  :: SimplAxCa
   real(ReKi),                             intent(in   )  :: SimplAxCaMG 
   real(ReKi),                             intent(in   )  :: SimplAxCp
   real(ReKi),                             intent(in   )  :: SimplAxCpMG 
   real(ReKi),                             intent(in   )  :: SimplCb
   real(ReKi),                             intent(in   )  :: SimplCbMG
   logical,                                intent(in   )  :: SimplMCF
   type(Morison_CoefMembersCyl), allocatable, intent(in   )  :: CoefMembers(:)
   integer(IntKi),                         intent(in   )  :: NCoefDpth
   type(Morison_CoefDpthsCyl),allocatable, intent(in   )  :: CoefDpths(:)
   type(Morison_NodeType),    allocatable, intent(in   )  :: nodes(:)
   type(Morison_MemberType),               intent(inout)  :: member
   
   integer(IntKi)                              :: i
   real(ReKi)                                  :: s
  
   select case ( MCoefMod )
      
   case (1)  ! Simple model : all nodes receive the same coefficients
      do i = 1, member%NElements + 1
         if ( member%tMG(i) > 0.0_ReKi ) then
            member%Cd    (i) = SimplCdMG
            member%Ca    (i) = SimplCaMG
            member%Cp    (i) = SimplCpMG
            member%AxCd  (i) = SimplAxCdMG
            member%AxCa  (i) = SimplAxCaMG
            member%AxCp  (i) = SimplAxCpMG
            member%Cb    (i) = SimplCbMG
         else
            member%Cd    (i) = SimplCd
            member%Ca    (i) = SimplCa
            member%Cp    (i) = SimplCp
            member%AxCd  (i) = SimplAxCd
            member%AxCa  (i) = SimplAxCa
            member%AxCp  (i) = SimplAxCp
            member%Cb    (i) = SimplCb
         end if
      end do
      member%PropMCF = SimplMCF
   CASE (2) ! Depth-based model: coefficients are set using depth-based table data
      do i = 1, member%NElements + 1
         CALL SetDepthBasedCoefs_Cyl( nodes(member%NodeIndx(i))%Position(3)+MSL2SWL,  member%tMG(i), NCoefDpth, CoefDpths, member%Cd(i), member%Ca(i), &
                                    member%Cp(i), member%AxCd(i), member%AxCa(i), member%AxCp(i), member%Cb(i) )
      end do
      member%PropMCF = CoefDpths(1)%DpthMCF
   CASE (3) ! Member-based model: coefficients set using member-specific coefficient tables
       do i = 1, member%NElements + 1
         ! Pull member  end-node data from the tables and then linearly interpolate it onto the interior member nodes    
         s = (real(i,ReKi)-1.0) / real(member%NElements,ReKi)
         if ( member%tMG(i) > 0.0_ReKi ) then
            member%Cd    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCdMG1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCdMG2  *s
            member%Ca    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCaMG1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCaMG2  *s
            member%Cp    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCpMG1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCpMG2  *s
            member%Cb    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCbMG1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCbMG2  *s
            member%AxCd  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCaMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCdMG2*s
            member%AxCa  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCaMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCaMG2*s
            member%AxCp  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCpMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCpMG2*s
         else
            member%Cd    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCd1    *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCd2    *s
            member%Ca    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCa1    *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCa2    *s
            member%Cp    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCp1    *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCp2    *s
            member%Cb    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCb1    *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCb2    *s
            member%AxCd  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCd1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCd2  *s
            member%AxCa  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCa1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCa2  *s
            member%AxCp  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCp1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCp2  *s
         end if
      end do
      member%propMCF = CoefMembers(MmbrCoefIDIndx)%MemberMCF
   end select
  
end subroutine SetExternalHydroCoefs_Cyl

!====================================================================================================
SUBROUTINE SetExternalHydroCoefs_Rec(  MSL2SWL, MCoefMod, MmbrCoefIDIndx, SimplCdA, SimplCdAMG, SimplCdB, SimplCdBMG, SimplCaA, SimplCaAMG, SimplCaB, SimplCaBMG, SimplCp, &
                                   SimplCpMG, SimplAxCd, SimplAxCdMG, SimplAxCa, SimplAxCaMG, SimplAxCp, SimplAxCpMG, SimplCb, SimplCbMG, SimplMCF, CoefMembers,    &
                                   NCoefDpth, CoefDpths, nodes, member )   
!     This private subroutine generates the Cd, Ca, Cp, Cb, CdMG, CaMG, CpMG, and CbMG coefs for the member based on
!     the input data.  
!---------------------------------------------------------------------------------------------------- 
   real(ReKi),                                intent(in   )  :: MSL2SWL
   integer(IntKi),                            intent(in   )  :: MCoefMod
   integer(IntKi),                            intent(in   )  :: MmbrCoefIDIndx
   real(ReKi),                                intent(in   )  :: SimplCdA 
   real(ReKi),                                intent(in   )  :: SimplCdAMG
   real(ReKi),                                intent(in   )  :: SimplCdB 
   real(ReKi),                                intent(in   )  :: SimplCdBMG
   real(ReKi),                                intent(in   )  :: SimplCaA
   real(ReKi),                                intent(in   )  :: SimplCaAMG 
   real(ReKi),                                intent(in   )  :: SimplCaB
   real(ReKi),                                intent(in   )  :: SimplCaBMG 
   real(ReKi),                                intent(in   )  :: SimplCp
   real(ReKi),                                intent(in   )  :: SimplCpMG 
   real(ReKi),                                intent(in   )  :: SimplAxCd
   real(ReKi),                                intent(in   )  :: SimplAxCdMG 
   real(ReKi),                                intent(in   )  :: SimplAxCa
   real(ReKi),                                intent(in   )  :: SimplAxCaMG 
   real(ReKi),                                intent(in   )  :: SimplAxCp
   real(ReKi),                                intent(in   )  :: SimplAxCpMG 
   real(ReKi),                                intent(in   )  :: SimplCb
   real(ReKi),                                intent(in   )  :: SimplCbMG
   logical,                                   intent(in   )  :: SimplMCF
   type(Morison_CoefMembersRec), allocatable, intent(in   )  :: CoefMembers(:)
   integer(IntKi),                            intent(in   )  :: NCoefDpth
   type(Morison_CoefDpthsRec),   allocatable, intent(in   )  :: CoefDpths(:)
   type(Morison_NodeType),       allocatable, intent(in   )  :: nodes(:)
   type(Morison_MemberType),                  intent(inout)  :: member
   
   integer(IntKi)                              :: i
   real(ReKi)                                  :: s
  
   select case ( MCoefMod )
      
   case (1)  ! Simple model : all nodes receive the same coefficients
      do i = 1, member%NElements + 1
         if ( member%tMG(i) > 0.0_ReKi ) then
            member%CdA   (i) = SimplCdAMG
            member%CdB   (i) = SimplCdBMG
            member%CaA   (i) = SimplCaAMG
            member%CaB   (i) = SimplCaBMG
            member%Cp    (i) = SimplCpMG
            member%AxCd  (i) = SimplAxCdMG
            member%AxCa  (i) = SimplAxCaMG
            member%AxCp  (i) = SimplAxCpMG
            member%Cb    (i) = SimplCbMG
         else
            member%CdA   (i) = SimplCdA
            member%CdB   (i) = SimplCdB
            member%CaA   (i) = SimplCaA
            member%CaB   (i) = SimplCaB
            member%Cp    (i) = SimplCp
            member%AxCd  (i) = SimplAxCd
            member%AxCa  (i) = SimplAxCa
            member%AxCp  (i) = SimplAxCp
            member%Cb    (i) = SimplCb
         end if
      end do
      member%PropMCF = SimplMCF
   CASE (2) ! Depth-based model: coefficients are set using depth-based table data
      do i = 1, member%NElements + 1
         CALL SetDepthBasedCoefs_Rec( nodes(member%NodeIndx(i))%Position(3)+MSL2SWL,  member%tMG(i), NCoefDpth, CoefDpths, member%CdA(i), member%CdB(i), & 
                                     member%CaA(i), member%CaB(i), member%Cp(i), member%AxCd(i), member%AxCa(i), member%AxCp(i), member%Cb(i) )
      end do
      member%PropMCF = CoefDpths(1)%DpthMCF
   CASE (3) ! Member-based model: coefficients set using member-specific coefficient tables
       do i = 1, member%NElements + 1
         ! Pull member  end-node data from the tables and then linearly interpolate it onto the interior member nodes    
         s = (real(i,ReKi)-1.0) / real(member%NElements,ReKi)
         if ( member%tMG(i) > 0.0_ReKi ) then
            member%CdA   (i) = CoefMembers(MmbrCoefIDIndx)%MemberCdAMG1 *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCdAMG2 *s
            member%CdB   (i) = CoefMembers(MmbrCoefIDIndx)%MemberCdBMG1 *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCdBMG2 *s
            member%CaA   (i) = CoefMembers(MmbrCoefIDIndx)%MemberCaAMG1 *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCaAMG2 *s
            member%CaB   (i) = CoefMembers(MmbrCoefIDIndx)%MemberCaBMG1 *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCaBMG2 *s
            member%Cp    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCpMG1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCpMG2  *s
            member%Cb    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCbMG1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCbMG2  *s
            member%AxCd  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCaMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCdMG2*s
            member%AxCa  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCaMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCaMG2*s
            member%AxCp  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCpMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCpMG2*s
         else
            member%CdA   (i) = CoefMembers(MmbrCoefIDIndx)%MemberCdA1   *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCdA2   *s
            member%CdB   (i) = CoefMembers(MmbrCoefIDIndx)%MemberCdB1   *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCdB2   *s
            member%CaA   (i) = CoefMembers(MmbrCoefIDIndx)%MemberCaA1   *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCaA2   *s
            member%CaB   (i) = CoefMembers(MmbrCoefIDIndx)%MemberCaB1   *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCaB2   *s
            member%Cp    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCp1    *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCp2    *s
            member%Cb    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCb1    *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCb2    *s
            member%AxCd  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCd1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCd2  *s
            member%AxCa  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCa1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCa2  *s
            member%AxCp  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCp1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCp2  *s
         end if
      end do
      member%propMCF = CoefMembers(MmbrCoefIDIndx)%MemberMCF
   end select
  
end subroutine SetExternalHydroCoefs_Rec


!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetNodeMG( numMGDepths, MGDepths, node, MSL2SWL, tMG, MGdensity )
   ! sets the margine growth thickness of a single node (previously all nodes)
   INTEGER,                                  INTENT( IN    )  :: numMGDepths
   TYPE(Morison_MGDepthsType), ALLOCATABLE,  INTENT( IN    )  :: MGDepths(:)
   TYPE(Morison_NodeType),                   INTENT( IN    )  :: node
   real(ReKi),                               intent( in    )  :: MSL2SWL
   real(ReKi),                               intent( inout )  :: tMG
   real(ReKi),                               intent( inout )  :: MGdensity
   
   INTEGER                 :: J
   REAL(ReKi)              :: z
   INTEGER                 :: indx1, indx2
   REAL(ReKi)              :: dd, s
   LOGICAL                 :: foundLess = .FALSE.
   

         !Find the table entry(ies) which match the node's depth value
      ! The assumption here is that the depth table is stored from largest
      ! to smallest in depth
      z = node%Position(3) + MSL2SWL ! Place in MSL coordinate system
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
         tMG       = 0.0
         MGdensity = 0.0
      ELSE
         ! Linearly interpolate the marine growth thickness and density based on depth      
         dd = MGDepths(indx1)%MGDpth - MGDepths(indx2)%MGDpth
         IF ( EqualRealNos(dd, 0.0_ReKi) ) THEN
            s = 0.0_ReKi
         ELSE
            s = ( MGDepths(indx1)%MGDpth - z ) / dd
         END IF
         tMG       = MGDepths(indx1)%MGThck*(1-s) + MGDepths(indx2)%MGThck*s
         MGdensity = MGDepths(indx1)%MGDens*(1-s) + MGDepths(indx2)%MGDens*s
      END IF
      

END SUBROUTINE SetNodeMG


!----------------------------------------------------------------------------------------------------------------------------------
subroutine AllocateMemberDataArrays( member, memberLoads, errStat, errMsg )
   type(Morison_MemberType),     intent (inout)  :: member
   type(Morison_MemberLoads),    intent (inout)  :: memberLoads
   integer(IntKi),               intent (  out)  :: errStat              ! returns a non-zero value when an error occurs            
   character(*),                 intent (  out)  :: errMsg               ! Error message if errStat /= ErrID_None
   
   integer(IntKi) :: errStat2              ! returns a non-zero value when an error occurs            
   CHARACTER(errMsgLen)  :: errMsg2     ! Error message if errStat2 /= ErrID_None
   character(*), parameter :: routineName = 'AllocateMemberDataArrays'
   
   errStat = ErrID_None
   errMSg  = ''
   call AllocAry(member%NodeIndx     , member%NElements+1, 'member%NodeIndx'     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%floodstatus  , member%NElements,   'member%floodstatus'  , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%alpha        , member%NElements,   'member%alpha'        , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%alpha_fb     , member%NElements,   'member%alpha_fb'     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%alpha_fb_star, member%NElements,   'member%alpha_fb_star', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%m_fb_l       , member%NElements,   'member%m_fb_l       ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%m_fb_u       , member%NElements,   'member%m_fb_u       ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%h_cfb_l      , member%NElements,   'member%h_cfb_l      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%h_cfb_u      , member%NElements,   'member%h_cfb_u      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%I_lfb_l      , member%NElements,   'member%I_lfb_l      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%I_lfb_u      , member%NElements,   'member%I_lfb_u      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%m_mg_l       , member%NElements,   'member%m_mg_l       ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%m_mg_u       , member%NElements,   'member%m_mg_u       ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%h_cmg_l      , member%NElements,   'member%h_cmg_l      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%h_cmg_u      , member%NElements,   'member%h_cmg_u      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%I_lmg_l      , member%NElements,   'member%I_lmg_l      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%I_lmg_u      , member%NElements,   'member%I_lmg_u      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%tMG          , member%NElements+1, 'member%tMG          ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%MGdensity    , member%NElements+1, 'member%MGdensity    ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%Cp           , member%NElements+1, 'member%Cp           ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%AxCd         , member%NElements+1, 'member%AxCd         ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%AxCa         , member%NElements+1, 'member%AxCa         ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%AxCp         , member%NElements+1, 'member%AxCp         ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%Cb           , member%NElements+1, 'member%Cb           ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_D    , 6, member%NElements+1, 'memberLoads%F_D'  , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_A    , 6, member%NElements+1, 'memberLoads%F_A'  , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_B    , 6, member%NElements+1, 'memberLoads%F_B'  , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_BF   , 6, member%NElements+1, 'memberLoads%F_BF' , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_I    , 6, member%NElements+1, 'memberLoads%F_I'  , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_If   , 6, member%NElements+1, 'memberLoads%F_If' , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_WMG  , 6, member%NElements+1, 'memberLoads%F_WMG', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_IMG  , 6, member%NElements+1, 'memberLoads%F_IMG', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)

   ! Shape dependent variables
   if (member%MSecGeom == MSecGeom_Cyl) then
      call AllocAry(member%dRdl_mg      , member%NElements,   'member%dRdl_mg'      , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%dRdl_mg_b    , member%NElements,   'member%dRdl_mg_b'    , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%dRdl_in      , member%NElements,   'member%dRdl_in'      , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%I_rfb_l      , member%NElements,   'member%I_rfb_l      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%I_rfb_u      , member%NElements,   'member%I_rfb_u      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%I_rmg_l      , member%NElements,   'member%I_rmg_l      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%I_rmg_u      , member%NElements,   'member%I_rmg_u      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%R            , member%NElements+1, 'member%R            ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%RMG          , member%NElements+1, 'member%RMG          ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%RMGB         , member%NElements+1, 'member%RMGB         ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%Rin          , member%NElements+1, 'member%Rin          ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%Cd           , member%NElements+1, 'member%Cd           ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%Ca           , member%NElements+1, 'member%Ca           ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   else if (member%MSecGeom == MSecGeom_Rec) then
      call AllocAry(member%dSadl_mg     , member%NElements,   'member%dSadl_mg'     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%dSadl_mg_b   , member%NElements,   'member%dSadl_mg_b'   , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%dSadl_in     , member%NElements,   'member%dSadl_in'     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%dSbdl_mg     , member%NElements,   'member%dSbdl_mg'     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%dSbdl_mg_b   , member%NElements,   'member%dSbdl_mg_b'   , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%dSbdl_in     , member%NElements,   'member%dSbdl_in'     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%I_xfb_l      , member%NElements,   'member%I_xfb_l      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%I_xfb_u      , member%NElements,   'member%I_xfb_u      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%I_yfb_l      , member%NElements,   'member%I_yfb_l      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%I_yfb_u      , member%NElements,   'member%I_yfb_u      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%I_xmg_l      , member%NElements,   'member%I_xmg_l      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%I_xmg_u      , member%NElements,   'member%I_xmg_u      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%I_ymg_l      , member%NElements,   'member%I_ymg_l      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%I_ymg_u      , member%NElements,   'member%I_ymg_u      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%Sa           , member%NElements+1, 'member%Sa           ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%SaMG         , member%NElements+1, 'member%SaMG         ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%SaMGB        , member%NElements+1, 'member%SaMGB        ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%Sain         , member%NElements+1, 'member%Sain         ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%Sb           , member%NElements+1, 'member%Sb           ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%SbMG         , member%NElements+1, 'member%SbMG         ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%SbMGB        , member%NElements+1, 'member%SbMGB        ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%Sbin         , member%NElements+1, 'member%Sbin         ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%CdA          , member%NElements+1, 'member%CdA          ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%CdB          , member%NElements+1, 'member%CdB          ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%CaA          , member%NElements+1, 'member%CaA          ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
      call AllocAry(member%CaB          , member%NElements+1, 'member%CaB          ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   end if
   
   if (ErrStat >= AbortErrLev) return

   ! Initialize everything to zero
   member%NodeIndx      = 0.0_ReKi
   member%floodstatus   = 0.0_ReKi
   member%alpha         = 0.0_ReKi
   member%alpha_fb      = 0.0_ReKi
   member%alpha_fb_star = 0.0_ReKi
   member%m_fb_l        = 0.0_ReKi
   member%m_fb_u        = 0.0_ReKi
   member%h_cfb_l       = 0.0_ReKi
   member%h_cfb_u       = 0.0_ReKi
   member%I_lfb_l       = 0.0_ReKi
   member%I_lfb_u       = 0.0_ReKi
   member%m_mg_l        = 0.0_ReKi
   member%m_mg_u        = 0.0_ReKi
   member%h_cmg_l       = 0.0_ReKi
   member%h_cmg_u       = 0.0_ReKi
   member%I_lmg_l       = 0.0_ReKi
   member%I_lmg_u       = 0.0_ReKi

   member%tMG           = 0.0_ReKi
   member%MGdensity     = 0.0_ReKi
   member%Cp            = 0.0_ReKi
   member%AxCd          = 0.0_ReKi
   member%AxCa          = 0.0_ReKi
   member%AxCp          = 0.0_ReKi
   member%Cb            = 0.0_ReKi
   memberLoads%F_D      = 0.0_ReKi
   memberLoads%F_A      = 0.0_ReKi
   memberLoads%F_B      = 0.0_ReKi
   memberLoads%F_BF     = 0.0_ReKi
   memberLoads%F_I      = 0.0_ReKi
   memberLoads%F_If     = 0.0_ReKi
   memberLoads%F_WMG    = 0.0_ReKi
   memberLoads%F_IMG    = 0.0_ReKi

   if (member%MSecGeom == MSecGeom_Cyl) then
      member%dRdl_mg       = 0.0_ReKi
      member%dRdl_mg_b     = 0.0_ReKi
      member%dRdl_in       = 0.0_ReKi
      member%I_rfb_l       = 0.0_ReKi
      member%I_rfb_u       = 0.0_ReKi
      member%I_rmg_l       = 0.0_ReKi
      member%I_rmg_u       = 0.0_ReKi
      member%R             = 0.0_ReKi
      member%RMG           = 0.0_ReKi
      member%RMGB          = 0.0_ReKi
      member%Rin           = 0.0_ReKi
      member%Cd            = 0.0_ReKi
      member%Ca            = 0.0_ReKi
   else if (member%MSecGeom == MSecGeom_Rec) then
      member%dSadl_mg      = 0.0_ReKi
      member%dSadl_mg_b    = 0.0_ReKi
      member%dSadl_in      = 0.0_ReKi
      member%dSbdl_mg      = 0.0_ReKi
      member%dSbdl_mg_b    = 0.0_ReKi
      member%dSbdl_in      = 0.0_ReKi
      member%I_xfb_l       = 0.0_ReKi
      member%I_xfb_u       = 0.0_ReKi
      member%I_yfb_l       = 0.0_ReKi
      member%I_yfb_u       = 0.0_ReKi
      member%I_xmg_l       = 0.0_ReKi
      member%I_xmg_u       = 0.0_ReKi
      member%I_ymg_l       = 0.0_ReKi
      member%I_ymg_u       = 0.0_ReKi
      member%Sa            = 0.0_ReKi
      member%SaMG          = 0.0_ReKi
      member%SaMGB         = 0.0_ReKi
      member%Sain          = 0.0_ReKi
      member%Sb            = 0.0_ReKi
      member%SbMG          = 0.0_ReKi
      member%SbMGB         = 0.0_ReKi
      member%Sbin          = 0.0_ReKi
      member%CdA           = 0.0_ReKi
      member%CdB           = 0.0_ReKi
      member%CaA           = 0.0_ReKi
      member%CaB           = 0.0_ReKi
   end if

end subroutine AllocateMemberDataArrays
!----------------------------------------------------------------------------------------------------------------------------------
subroutine FlipMemberNodeData( member, nodes, doSwap)
   type(Morison_MemberType),     intent (inout)  :: member
   type(Morison_NodeType),       intent (in   )  :: nodes(:)
   logical,                      intent (  out)  :: doSwap
   
   integer(IntKi) :: i, j1, j2, numMemNodes, indx
   
   
   doSwap = .FALSE.
   numMemNodes = member%NElements + 1
   j1 = member%NodeIndx(1)
   j2 = member%NodeIndx(numMemNodes)
   IF ( EqualRealNos(nodes(j1)%Position(3), nodes(j2)%Position(3) ) ) THEN         ! Z1 = Z2          
      IF ( EqualRealNos(nodes(j1)%Position(1), nodes(j2)%Position(1) ) ) THEN      ! X1 = X2
         IF   ( nodes(j1)%Position(2) > nodes(j2)%Position(2) ) THEN
            doSwap = .TRUE.  ! Y1 > Y2
         END IF
      ELSE IF ( nodes(j1)%Position(1) > nodes(j2)%Position(1) ) THEN
         doSwap = .TRUE.  ! X1 > X2
      END IF
   ELSE IF    ( nodes(j1)%Position(3) > nodes(j2)%Position(3) ) THEN
      doSwap = .TRUE.                                ! Z1 > Z2  
   END IF
         
   ! If we swap the the nodes, we need know this later when calculating the normal vector to the ends
   member%Flipped = doSwap
   IF ( doSwap ) THEN
      member%NodeIndx(1) = j2
      member%NodeIndx(numMemNodes) = j1
      
      ! Loop over half the interior nodes and swap their indices
      do i = 1, ceiling( (numMemNodes-2.0_ReKi)/2.0_ReKi)
         indx = member%NodeIndx(1+i)
         member%NodeIndx(1+i) = member%NodeIndx(numMemNodes-i)
         member%NodeIndx(numMemNodes-i) = indx
      end do

      ! Flip the sign of the spin angle
      member%MSpinOrient = -member%MSpinOrient
      
   end if    
   
end subroutine FlipMemberNodeData
!----------------------------------------------------------------------------------------------------------------------------------
subroutine SetMemberProperties_Cyl( gravity, member, MCoefMod, MmbrCoefIDIndx, MmbrFilledIDIndx, propSet1, propSet2, InitInp, errStat, errMsg )
   real(ReKi),                   intent (in   )  :: gravity
   type(Morison_MemberType),     intent (inout)  :: member
   integer(IntKi),               intent (in   )  :: MCoefMod
   integer(IntKi),               intent (in   )  :: MmbrCoefIDIndx
   integer(IntKi),               intent (in   )  :: MmbrFilledIDIndx
   type(Morison_MemberPropTypeCyl), intent (in   )  :: propSet1             ! property set of node 1
   type(Morison_MemberPropTypeCyl), intent (in   )  :: propSet2             ! property set of node N+1
   type(Morison_InitInputType),  intent (in   )  :: InitInp
   integer(IntKi),               intent (  out)  :: errStat              ! returns a non-zero value when an error occurs            
   character(*),                 intent (  out)  :: errMsg               ! Error message if errStat /= ErrID_None

   character(*), parameter                       :: RoutineName = 'SetMemberProperties_Cyl'

   integer(IntKi) :: N, i
   real(ReKi)     :: s, dl
   real(ReKi)     :: vec(3)
   real(ReKi)     :: memLength 
   real(ReKi)     :: Za 
   real(ReKi)     :: Zb 
   real(ReKi)     :: phi 
   real(ReKi)     :: sinPhi
   real(ReKi)     :: cosPhi
   real(ReKi)     :: Rmid  
   real(ReKi)     :: RmidMG
   real(ReKi)     :: Rmidin
   real(ReKi)     :: Lmid
   real(ReKi)     :: li
   real(ReKi)     :: Vinner_l, Vinner_u, Vouter_l, Vouter_u, Vballast_l, Vballast_u
   real(ReKi)     :: tk(1,3), Imat(3,3)
   REAL(ReKi)     :: h_c    ! center of mass offset from first node
   
   errStat = ErrID_None
   errMSg  = ''
   
   N  = member%NElements
   dl = member%dl
   
   vec     = InitInp%Nodes(member%NodeIndx(N+1))%Position - InitInp%Nodes(member%NodeIndx(1))%Position   
   
   ! calculate reference orientation information.  Note: members are straight to start
   memLength = member%RefLength 
   member%k(1:3) = (vec/memLength)  ! vector along member from start to end point, length > 0 was already checked when the members were parsed and generated from the input file data
   tk(1,1) = member%k(1)
   tk(1,2) = member%k(2)
   tk(1,3) = member%k(3)
   member%kkt    = matmul(transpose(tk),tk)
   call Eye(Imat,errStat,errMsg)
   member%Ak     =  Imat - member%kkt
   phi = acos( max(-1.0_ReKi, min(1.0_ReKi, vec(3)/memLength) ) )  ! incline angle   
   sinPhi = sin(phi)
   cosPhi = cos(phi)  
   member%cosPhi_ref = cosPhi
   
   ! These are all per node and not done here, yet
   
   do i = 1, member%NElements+1
      call SetNodeMG( InitInp%NMGDepths, InitInp%MGDepths, InitInp%Nodes(member%NodeIndx(i)), InitInp%WaveField%MSL2SWL, member%tMG(i), member%MGDensity(i) )
   end do

   member%R(  1)   = propSet1%PropD / 2.0            
   member%RMG(1)   = propSet1%PropD / 2.0 + member%tMG(1) 
   member%Rin(1)   = propSet1%PropD / 2.0 - propSet1%PropThck  
   member%R(  N+1) = propSet2%PropD / 2.0            
   member%RMG(N+1) = propSet2%PropD / 2.0 + member%tMG(N+1)
   member%Rin(N+1) = propSet2%PropD / 2.0 - propSet2%PropThck 
   do i = 2,  member%NElements
      s = (real(i,ReKi)-1.0) / real(member%NElements,ReKi)
      member%R(  i) =  member%R(  1)*(1-s) + member%R(  N+1)*s
      member%Rin(i) =  member%Rin(1)*(1-s) + member%Rin(N+1)*s
      member%RMG(i) =  member%R(i) + member%tMG(i)
   end do

   call SetExternalHydroCoefs_Cyl( InitInp%WaveField%MSL2SWL, MCoefMod, MmbrCoefIDIndx, InitInp%SimplCd, InitInp%SimplCdMG, InitInp%SimplCa, InitInp%SimplCaMG, InitInp%SimplCp, &
                                   InitInp%SimplCpMG, InitInp%SimplAxCd, InitInp%SimplAxCdMG, InitInp%SimplAxCa, InitInp%SimplAxCaMG, InitInp%SimplAxCp, InitInp%SimplAxCpMG, &
                                   InitInp%SimplCb, InitInp%SimplCbMG, InitInp%SimplMCF, & 
                                   InitInp%CoefMembersCyl, InitInp%NCoefDpthCyl, InitInp%CoefDpthsCyl, InitInp%Nodes, member )
   
   ! calculate member radius with marine growth scaled by sqrt(Cb) for buoyancy/hydrostatic load calculation
   do i = 1, member%NElements+1
      member%RMGB(i) = member%RMG(i) * SQRT(member%Cb(i))
   end do

   ! calculate reference incline angle and heading, and related trig values.  Note: members are straight to start
   Za = InitInp%Nodes(member%NodeIndx(1  ))%Position(3) 
   Zb = InitInp%Nodes(member%NodeIndx(N+1))%Position(3)

   ! Check if members with the MacCamy-Fuchs diffraction model and not modeled by potential flow satisfy the necessary criteria.
   IF ( member%PropMCF .AND. ( .NOT. member%PropPot )) THEN
      ! Check if surface piercing
      IF ( Za*Zb > 0 ) THEN ! Two end joints of the member on the same side of the SWL
         CALL SetErrStat(ErrID_Fatal, 'MacCamy-Fuchs members must be surface piercing.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, RoutineName )   
         RETURN
      END IF
      ! Check inclination
      If ( ABS(phi) .GE. 0.174533 ) THEN ! If inclination from vertical is greater than 10 deg
         CALL SetErrStat(ErrID_Fatal, 'MacCamy-Fuchs members must be within 10 degrees from vertical.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, RoutineName )   
         RETURN
      END IF
      ! Check radius
      DO i = 1, member%NElements+1
         IF ( (member%RMG(i) .GT. 1.1_ReKi*REAL(0.5_SiKi*InitInp%WaveField%MCFD)) .OR. (member%RMG(i) .LT. 0.9_ReKi*REAL(0.5_SiKi*InitInp%WaveField%MCFD)) ) THEN
            ! Error because MacCamy-Fuchs members must have a diameter within +/-10% of MCFD specified in seastate.
            CALL SetErrStat(ErrID_Fatal, 'MacCamy-Fuchs members must have a diameter within +/-10% of MCFD specified in the SeaState input file.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, RoutineName )   
            RETURN
         END IF
      END DO
      ! Check draft-to-radius ratio
      IF ( (-InitInp%Nodes(member%NodeIndx(1))%Position(3)) < 0.5_SiKi*InitInp%WaveField%MCFD ) THEN
         CALL SetErrStat(ErrID_Fatal, 'Initial draft of MacCamy-Fuchs members should be at least as large as their radius.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, RoutineName )   
         RETURN
      END IF
   END IF

   ! find fill location of member (previously in SetElementFillProps)
   member%MmbrFilledIDIndx = MmbrFilledIDIndx ! Set this to the parameter version of this member data
   if ( MmbrFilledIDIndx > 0 ) then    
      member%FillDens     =  InitInp%FilledGroups(MmbrFilledIDIndx)%FillDens
      member%FillFSLoc    =  InitInp%FilledGroups(MmbrFilledIDIndx)%FillFSLoc - InitInp%WaveField%MSL2SWL
       if (member%FillFSLoc >= Zb) then
         member%z_overfill = member%FillFSLoc - Zb
         member%l_fill = member%RefLength
         member%memfloodstatus = 1  ! fully flooded   
       elseif (Za >= member%FillFSLoc) then
          ! No ballast
         member%memfloodstatus = 0  
         member%z_overfill = 0.0_ReKi
         member%l_fill = 0.0_ReKi
      else
         member%z_overfill =0
         if ( Zb <= -InitInp%WaveField%EffWtrDpth ) then
            member%memfloodstatus = 0  ! member fully buried in seabed
            member%l_fill = 0
         else
            member%memfloodstatus = 2  ! partially flooded member
            member%l_fill = (member%FillFSLoc - Za)/cosPhi
         end if
      
      end if
      
   else
      member%FillDens     =   0.0
      member%FillFSLoc    =   0.0  ! Future calculations for ballasting MUST verify that MbrFilledIDIndx > 0 for any ballasting calcs or this value will cause errors
      member%z_overfill =0
      member%l_fill = 0
      member%memfloodstatus = 0
   end if

    ! Check the member does not exhibit any of the following conditions
   if (.not. member%PropPot) then 
      if (member%MHstLMod == 1) then ! Only cylindrical members are allowed to use MHstLMod = 1
         if ( abs(Zb) < abs(member%Rmg(N+1)*sinPhi) ) then
            call SetErrStat(ErrID_Fatal, 'The upper end-plate of a member must not cross the water plane.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, RoutineName )   
         end if
         if ( abs(Za) < abs(member%Rmg(1)*sinPhi) ) then
            call SetErrStat(ErrID_Fatal, 'The lower end-plate of a member must not cross the water plane.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, RoutineName )   
         end if
      end if
      if ( ( Za < -InitInp%WaveField%EffWtrDpth .and. Zb >= -InitInp%WaveField%EffWtrDpth ) .and. ( phi > 10.0*d2r .or. abs((member%RMG(N+1) - member%RMG(1))/member%RefLength)>0.1 ) ) then
         call SetErrStat(ErrID_Fatal, 'A member which crosses the seabed must not be inclined more than 10 degrees from vertical or have a taper larger than 0.1.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, RoutineName )   
      end if
      
   end if
   

   ! calculate h_floor if seabed-piercing
   member%h_floor = 0.0_ReKi
   member%i_floor = member%NElements+1  ! Default to entire member is below the seabed
   member%doEndBuoyancy = .false.
   if (Za < -InitInp%WaveField%EffWtrDpth) then
      do i= 2, member%NElements+1
         Za = InitInp%Nodes(member%NodeIndx(i))%Position(3)
         if (Za > -InitInp%WaveField%EffWtrDpth) then            ! find the lowest node above the seabed
            
            if (cosPhi < 0.173648178 ) then ! phi > 80 degrees and member is seabed crossing
               call SetErrStat(ErrID_Fatal, 'A seabed crossing member must have an inclination angle of <= 80 degrees from vertical.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, RoutineName )
            end if
            
            member%h_floor = (-InitInp%WaveField%EffWtrDpth-Za)/cosPhi  ! get the distance from the node to the seabed along the member axis (negative value)
            member%i_floor = i-1                    ! record the number of the element that pierces the seabed
            member%doEndBuoyancy = .true.
            exit
         else if ( EqualRealNos(Za, -InitInp%WaveField%EffWtrDpth ) ) then
            member%doEndBuoyancy = .true.
         end if
      end do
   else
      member%i_floor = 0 ! lower end is at or above the seabed
   end if


  
   ! calculate element-level values
   
   do i = 1, member%NElements
      member%dRdl_mg(  i) = (member%RMG( i+1) - member%RMG( i))/dl
      member%dRdl_in(  i) = (member%Rin( i+1) - member%Rin( i))/dl
      member%dRdl_mg_b(i) = (member%RMGB(i+1) - member%RMGB(i))/dl
      
      member%alpha(   i) = GetAlphaCyl(member%RMGB(i), member%RMGB(i+1))   ! Only used to distribute external buoyancy load to nodes
      member%alpha_fb(i) = GetAlphaCyl(member%Rin( i), member%Rin( i+1))
      
   end do

   member%Vinner    = 0.0_ReKi  ! Total  volume of member without marine growth
   member%Vouter    = 0.0_ReKi  ! Total outer volume of member including marine growth
   member%Vballast  = 0.0_ReKi  ! Total ballasted volume of member
   member%elem_fill = 0         ! Last (partially) filled element of the member
   member%h_fill    = 0.0       ! Axial length of elem_fill occupied by water ballast
   
   ! force-related constants for each element
   do i = 1, member%NElements
   
      Za = InitInp%Nodes(member%NodeIndx(  i))%Position(3)   ! z location of node i
      Zb = InitInp%Nodes(member%NodeIndx(i+1))%Position(3)   ! z location of node i+1
      
      ! ------------------ marine growth weight and inertia ------------------------------------------------
      Vinner_l   = 0.0
      Vouter_l   = 0.0
      Vinner_U   = 0.0
      Vouter_U   = 0.0
      if (i > member%i_floor) then         
         ! full marine growth: get the properties for each half-element lumped to the appropriate node
                  
         Rmid   = 0.5*(member%R(  i)+member%R(  i+1))  ! radius at middle of segment, where division occurs
         RmidMG = 0.5*(member%RMG(i)+member%RMG(i+1))  ! radius with marine growth at middle of segment, where division occurs
         Lmid   = 0.5*dl   ! = 0.5*(R2-R1)/m  half-length of segment

         CALL MarineGrowthPartSegmentCyl(member%R(i  ), Rmid, member%RMG(i  ),RmidMG, Lmid, member%MGDensity(i),  Vinner_l, Vouter_l, member%m_mg_l(i), member%h_cmg_l(i), member%I_lmg_l(i), member%I_rmg_l(i))   ! get precomputed quantities for lower half-segment
         CALL MarineGrowthPartSegmentCyl(member%R(i+1), Rmid, member%RMG(i+1),RmidMG,-Lmid, member%MGDensity(i),  Vinner_u, Vouter_u, member%m_mg_u(i), member%h_cmg_u(i), member%I_lmg_u(i), member%I_rmg_u(i))   ! get precomputed quantities for upper half-segment
         
      else if (i == member%i_floor) then         
         ! crossing seabed: get the properties for part-element above the seabed and lump to the upper node      

         Rmid   = (-member%h_floor*member%R(  i) +(dl+member%h_floor)*member%R(  i+1))/dl
         RmidMG = (-member%h_floor*member%RMG(i) +(dl+member%h_floor)*member%RMG(i+1))/dl
         Lmid   = -member%h_floor

         CALL MarineGrowthPartSegmentCyl(member%R(i+1), Rmid, member%RMG(i+1),RmidMG, -Lmid, member%MGDensity(i),  Vinner_u, Vouter_u, member%m_mg_u(i), member%h_cmg_u(i), member%I_lmg_u(i), member%I_rmg_u(i))   ! get precomputed quantities for upper half-segment
         Vinner_l   = 0.0
         Vouter_l   = 0.0
      end if

      ! ------------------ flooded ballast inertia ---------------------------------------------------------
      Vballast_l = 0.0
      Vballast_U = 0.0
      if (member%memfloodstatus > 0 .and. (member%FillFSLoc > Za)) then
         ! Fully filled element, so split in middle
         if ((i > member%i_floor) .and. (member%FillFSLoc >= Zb)) then

            ! get the properties for each half-element lumped to the appropriate node
            Rmidin = 0.5*(member%Rin(i)+member%Rin(i+1))  ! radius of member interior at middle of segment, where division occurs
            Lmid   = 0.5*dl   ! = 0.5*(R2-R1)/m  half-length of segment
            CALL FloodedBallastPartSegmentCyl(member%Rin(i  ), Rmidin,  Lmid, member%FillDens, Vballast_l, member%m_fb_l(i), member%h_cfb_l(i), member%I_lfb_l(i), member%I_rfb_l(i))   ! get precomputed quantities for lower half-segment
            CALL FloodedBallastPartSegmentCyl(member%Rin(i+1), Rmidin, -Lmid, member%FillDens, Vballast_u, member%m_fb_u(i), member%h_cfb_u(i), member%I_lfb_u(i), member%I_rfb_u(i))   ! get precomputed quantities for upper half-segment
 
         ! partially filled element, so split at FillFSLoc
         else if ((i > member%i_floor)  .AND. (member%FillFSLoc < Zb)) then

            ! get the properties for each partial-element lumped to the appropriate node
            Lmid   = member%FillFSLoc - Za 
            Rmidin = member%Rin(i)+(Lmid/(Zb-Za))*(member%Rin(i+1)-member%Rin(i))  ! radius of member interior at middle of segment, where division occurs
            CALL FloodedBallastPartSegmentCyl(member%Rin(i  ), Rmidin,  Lmid, member%FillDens, Vballast_l, member%m_fb_l(i), member%h_cfb_l(i), member%I_lfb_l(i), member%I_rfb_l(i))   ! get precomputed quantities for lower half-segment
            CALL FloodedBallastPartSegmentCyl(member%Rin(i+1), Rmidin, -Lmid, 0.0_ReKi, Vballast_u, member%m_fb_u(i), member%h_cfb_u(i), member%I_lfb_u(i), member%I_rfb_u(i))   ! get precomputed quantities for upper half-segment
 
         else if (i == member%i_floor) then     ! Hopefully we don't have a partially filled element crossing the seabed.
 
            ! crossing seabed: get the properties for part-element above the seabed and lump to the upper node
            Rmidin = (-member%h_floor*member%Rin(i) +(dl+member%h_floor)*member%Rin(i+1))/dl
            Lmid   = -member%h_floor
            CALL FloodedBallastPartSegmentCyl(member%Rin(i+1), Rmidin, -Lmid, member%FillDens,  Vballast_u, member%m_fb_u(i), member%h_cfb_u(i), member%I_lfb_u(i), member%I_rfb_u(i))   ! get precomputed quantities for upper half-segment
            Vballast_l = 0.0
 
         end if
      else  ! Either no ballast flooding in member, or this particular element isn't flooded at all
         Vballast_u        = 0.0
         Vballast_l        = 0.0
         member%m_fb_u(i)  = 0.0
         member%h_cfb_u(i) = 0.0
         member%I_lfb_u(i) = 0.0
         member%I_rfb_u(i) = 0.0
      endif
         

      
      ! Determine volumes to add to Non-WAMIT modeled members, etc.
      if (.not. member%PropPot) then
         
         if (Zb < -InitInp%WaveField%EffWtrDpth) then
            ! fully buried element, do not add these volume contributions to totals
         else if (0.0 >= Zb) then   ! Bug fix per OpenFAST issue #844   GJH 2/3/2022
            ! fully submerged elements.  
            ! NOTE: For an element which is fractionaly in the seabed, the entire element volume is added to totals
            member%Vinner = member%Vinner + Vinner_l + Vinner_u
            member%Vouter = member%Vouter + Vouter_l + Vouter_u
            member%Vsubmerged = member%Vsubmerged + Vouter_l + Vouter_u
         else if ((0.0 > Za) .AND. (0.0 < Zb)) then ! Bug fix per OpenFAST issue #844   GJH 2/3/2022
            ! if (i == 1) then
            !    call SetErrStat(ErrID_Fatal, 'The lowest element of a member must not cross the free surface.  This is true for MemberID '//trim(num2lstr(member%MemberID)), errStat, errMsg, RoutineName)
            ! end if
            
            ! partially submerged element
            member%Vinner = member%Vinner + Vinner_l + Vinner_u
            member%Vouter = member%Vouter + Vouter_l + Vouter_u
            ! compute volume portion which is submerged
            Lmid = -Za/cosPhi 
            call CylTaperCalc( member%Rmg(i), member%Rmg(i)+Lmid*member%dRdl_mg(i), Lmid, Vouter_l, h_c)
            
            member%Vsubmerged = member%Vsubmerged + Vouter_l 
            
         else ! fully above the water
            member%Vinner = member%Vinner + Vinner_l + Vinner_u
            member%Vouter = member%Vouter + Vouter_l + Vouter_u
         end if 
      end if
      
      ! ------------------ flooded ballast weight (done) --------------------
      ! NOTE: this section of code is somewhat redundant with "flooded ballast inertia" section above

      li = dl*(i-1)

      if (Zb < -InitInp%WaveField%EffWtrDpth) then                                                                             ! Fully buried element

         member%floodstatus(i) = 0
      
      else if (member%memfloodstatus > 0 .and. member%FillFSLoc >= Zb) then                                                    ! Fully flooded elements

         member%floodstatus(i) = 1
         if ( EqualRealNos(member%FillFSLoc, Zb) .or. (i==member%NElements) ) then  ! No partially filled elements
            member%elem_fill = i
            member%h_fill    = dl
         end if
         member%Vballast = member%Vballast + Vballast_l + Vballast_u

      else if ((member%memfloodstatus > 0) .and. (member%FillFSLoc > Za) .AND. (member%FillFSLoc < Zb)) then                   ! Partially flooded element
         
         member%floodstatus(i) = 2
         member%elem_fill      = i
         member%h_fill         = member%l_fill - (i-1)*dl
         call CylTaperCalc( member%Rin(i), member%Rin(i)+member%h_fill*member%dRdl_in(i), member%h_fill, Vballast_l, h_c)
         member%Vballast = member%Vballast + Vballast_l ! Note: Vballast_l will match calculations above

      else                                                                                                                     ! Unflooded element

         member%floodstatus(i) = 0
      
      end if

   end do ! end looping through elements

end subroutine SetMemberProperties_Cyl

subroutine SetMemberProperties_Rec( gravity, member, MCoefMod, MmbrCoefIDIndx, MmbrFilledIDIndx, propSet1, propSet2, InitInp, errStat, errMsg )
   real(ReKi),                   intent (in   )  :: gravity
   type(Morison_MemberType),     intent (inout)  :: member
   integer(IntKi),               intent (in   )  :: MCoefMod
   integer(IntKi),               intent (in   )  :: MmbrCoefIDIndx
   integer(IntKi),               intent (in   )  :: MmbrFilledIDIndx
   type(Morison_MemberPropTypeRec), intent (in   )  :: propSet1             ! property set of node 1
   type(Morison_MemberPropTypeRec), intent (in   )  :: propSet2             ! property set of node N+1
   type(Morison_InitInputType),  intent (in   )  :: InitInp
   integer(IntKi),               intent (  out)  :: errStat              ! returns a non-zero value when an error occurs            
   character(*),                 intent (  out)  :: errMsg               ! Error message if errStat /= ErrID_None

   character(*), parameter                       :: RoutineName = 'SetMemberProperties_Rec'

   integer(IntKi) :: N, i
   real(ReKi)     :: s, dl
   real(ReKi)     :: vec(3)
   real(ReKi)     :: memLength 
   real(ReKi)     :: Za 
   real(ReKi)     :: Zb 
   real(ReKi)     :: phi 
   real(ReKi)     :: sinPhi
   real(ReKi)     :: cosPhi
   real(ReKi)     :: SaMid, SbMid  
   real(ReKi)     :: SaMidMG, SbMidMG
   real(ReKi)     :: SaMidIn, SbMidIn
   real(ReKi)     :: Lmid
   real(ReKi)     :: li
   real(ReKi)     :: Vinner_l, Vinner_u, Vouter_l, Vouter_u, Vballast_l, Vballast_u
   real(ReKi)     :: tk(1,3), Imat(3,3), CMatrix(3,3)
   REAL(ReKi)     :: h_c    ! center of mass offset from first node
   
   errStat = ErrID_None
   errMSg  = ''
   
   N  = member%NElements
   dl = member%dl
   
   vec     = InitInp%Nodes(member%NodeIndx(N+1))%Position - InitInp%Nodes(member%NodeIndx(1))%Position   
   
   ! calculate reference orientation information.  Note: members are straight to start
   memLength = member%RefLength 
   member%k(1:3) = (vec/memLength)  ! vector along member from start to end point, length > 0 was already checked when the members were parsed and generated from the input file data
   tk(1,1) = member%k(1)
   tk(1,2) = member%k(2)
   tk(1,3) = member%k(3)
   member%kkt    = matmul(transpose(tk),tk)
   call Eye(Imat,errStat,errMsg)
   member%Ak     =  Imat - member%kkt
   IF (member%MSecGeom == MSecGeom_Rec) THEN
      CALL Morison_DirCosMtrx( InitInp%Nodes(member%NodeIndx(1))%Position, InitInp%Nodes(member%NodeIndx(N+1))%Position, member%MSpinOrient, CMatrix )
      member%x_hat = CMatrix(1:3,1)
      member%y_hat = CMatrix(1:3,2)
   END IF
   phi = acos( max(-1.0_ReKi, min(1.0_ReKi, vec(3)/memLength) ) )  ! incline angle   
   sinPhi = sin(phi)
   cosPhi = cos(phi)  
   member%cosPhi_ref = cosPhi
   
   ! These are all per node and not done here, yet
   
   do i = 1, member%NElements+1
      call SetNodeMG( InitInp%NMGDepths, InitInp%MGDepths, InitInp%Nodes(member%NodeIndx(i)), InitInp%WaveField%MSL2SWL, member%tMG(i), member%MGDensity(i) )
   end do

   member%Sa(  1)   = propSet1%PropA
   member%SaMG(1)   = propSet1%PropA + 2.0 * member%tMG(1) 
   member%Sain(1)   = propSet1%PropA - 2.0 * propSet1%PropThck  
   member%Sb(  1)   = propSet1%PropB            
   member%SbMG(1)   = propSet1%PropB + 2.0 * member%tMG(1) 
   member%Sbin(1)   = propSet1%PropB - 2.0 * propSet1%PropThck  
   member%Sa(  N+1) = propSet2%PropA
   member%SaMG(N+1) = propSet2%PropA + 2.0 * member%tMG(N+1)
   member%Sain(N+1) = propSet2%PropA - 2.0 * propSet2%PropThck 
   member%Sb(  N+1) = propSet2%PropB
   member%SbMG(N+1) = propSet2%PropB + 2.0 * member%tMG(N+1)
   member%Sbin(N+1) = propSet2%PropB - 2.0 * propSet2%PropThck 
   do i = 2,  member%NElements
      s = (real(i,ReKi)-1.0) / real(member%NElements,ReKi)
      member%Sa(  i) =  member%Sa(  1)*(1-s) + member%Sa(  N+1)*s
      member%Sain(i) =  member%Sain(1)*(1-s) + member%Sain(N+1)*s
      member%SaMG(i) =  member%Sa(i) + 2.0 * member%tMG(i)
      member%Sb(  i) =  member%Sb(  1)*(1-s) + member%Sb(  N+1)*s
      member%Sbin(i) =  member%Sbin(1)*(1-s) + member%Sbin(N+1)*s
      member%SbMG(i) =  member%Sb(i) + 2.0 * member%tMG(i)
   end do

   call SetExternalHydroCoefs_Rec( InitInp%WaveField%MSL2SWL, MCoefMod, MmbrCoefIDIndx, InitInp%SimplRecCdA, InitInp%SimplRecCdAMG, InitInp%SimplRecCdB, InitInp%SimplRecCdBMG, &
                                   InitInp%SimplRecCaA, InitInp%SimplRecCaAMG, InitInp%SimplRecCaB, InitInp%SimplRecCaBMG, InitInp%SimplRecCp, &
                                   InitInp%SimplRecCpMG, InitInp%SimplRecAxCd, InitInp%SimplRecAxCdMG, InitInp%SimplRecAxCa, InitInp%SimplRecAxCaMG, InitInp%SimplRecAxCp, InitInp%SimplRecAxCpMG, &
                                   InitInp%SimplRecCb, InitInp%SimplRecCbMG, InitInp%SimplRecMCF, & 
                                   InitInp%CoefMembersRec, InitInp%NCoefDpthRec, InitInp%CoefDpthsRec, InitInp%Nodes, member )
   
   ! calculate member radius with marine growth scaled by sqrt(Cb) for buoyancy/hydrostatic load calculation
   do i = 1, member%NElements+1
      member%SaMGB(i) = member%SaMG(i) * SQRT(member%Cb(i))
      member%SbMGB(i) = member%SbMG(i) * SQRT(member%Cb(i))
   end do

   ! calculate reference incline angle and heading, and related trig values.  Note: members are straight to start
   Za = InitInp%Nodes(member%NodeIndx(1  ))%Position(3) 
   Zb = InitInp%Nodes(member%NodeIndx(N+1))%Position(3)

   ! Check if members with the MacCamy-Fuchs diffraction model and not modeled by potential flow satisfy the necessary criteria.
   IF ( member%PropMCF .AND. ( .NOT. member%PropPot )) THEN
      CALL SetErrStat(ErrID_Fatal, 'MacCamy-Fuchs diffraction correction cannot be applied to rectangular members. Check member with ID '//trim(num2lstr(member%MemberID))//'. ', errStat, errMsg, RoutineName )   
   END IF

   ! find fill location of member (previously in SetElementFillProps)
   member%MmbrFilledIDIndx = MmbrFilledIDIndx ! Set this to the parameter version of this member data
   if ( MmbrFilledIDIndx > 0 ) then    
      member%FillDens     =  InitInp%FilledGroups(MmbrFilledIDIndx)%FillDens
      member%FillFSLoc    =  InitInp%FilledGroups(MmbrFilledIDIndx)%FillFSLoc - InitInp%WaveField%MSL2SWL
       if (member%FillFSLoc >= Zb) then
         member%z_overfill = member%FillFSLoc - Zb
         member%l_fill = member%RefLength
         member%memfloodstatus = 1  ! fully flooded   
       elseif (Za >= member%FillFSLoc) then
          ! No ballast
         member%memfloodstatus = 0  
         member%z_overfill = 0.0_ReKi
         member%l_fill = 0.0_ReKi
      else
         member%z_overfill =0
         if ( Zb <= -InitInp%WaveField%EffWtrDpth ) then
            member%memfloodstatus = 0  ! member fully buried in seabed
            member%l_fill = 0
         else
            member%memfloodstatus = 2  ! partially flooded member
            member%l_fill = (member%FillFSLoc - Za)/cosPhi
         end if
      
      end if
      
   else
      member%FillDens     =   0.0
      member%FillFSLoc    =   0.0  ! Future calculations for ballasting MUST verify that MbrFilledIDIndx > 0 for any ballasting calcs or this value will cause errors
      member%z_overfill =0
      member%l_fill = 0
      member%memfloodstatus = 0
   end if

    ! Check the member does not exhibit any of the following conditions
   if (.not. member%PropPot) then 
      ! MHstLMod=1 is not allowed for rectangular members at the moment. Skip the following check.
      ! if (member%MHstLMod == 1) then
      !    if ( abs(Zb) < abs(member%Rmg(N+1)*sinPhi) ) then
      !       call SetErrStat(ErrID_Fatal, 'The upper end-plate of a member must not cross the water plane.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, RoutineName )   
      !    end if
      !    if ( abs(Za) < abs(member%Rmg(1)*sinPhi) ) then
      !       call SetErrStat(ErrID_Fatal, 'The lower end-plate of a member must not cross the water plane.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, RoutineName )   
      !    end if
      ! end if
      if ( ( Za < -InitInp%WaveField%EffWtrDpth .and. Zb >= -InitInp%WaveField%EffWtrDpth ) .and. ( phi > 10.0*d2r .or. abs((member%SaMG(N+1) - member%SaMG(1))/member%RefLength)>0.1 .or. abs((member%SbMG(N+1) - member%SbMG(1))/member%RefLength)>0.1 ) ) then
         call SetErrStat(ErrID_Fatal, 'A member which crosses the seabed must not be inclined more than 10 degrees from vertical or have a taper larger than 0.1.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, RoutineName )   
      end if      
   end if

   ! calculate h_floor if seabed-piercing
   member%h_floor = 0.0_ReKi
   member%i_floor = member%NElements+1  ! Default to entire member is below the seabed
   member%doEndBuoyancy = .false.
   if (Za < -InitInp%WaveField%EffWtrDpth) then
      do i= 2, member%NElements+1
         Za = InitInp%Nodes(member%NodeIndx(i))%Position(3)
         if (Za > -InitInp%WaveField%EffWtrDpth) then            ! find the lowest node above the seabed
            
            if (cosPhi < 0.173648178 ) then ! phi > 80 degrees and member is seabed crossing
               call SetErrStat(ErrID_Fatal, 'A seabed crossing member must have an inclination angle of <= 80 degrees from vertical.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, RoutineName )
            end if
            
            member%h_floor = (-InitInp%WaveField%EffWtrDpth-Za)/cosPhi  ! get the distance from the node to the seabed along the member axis (negative value)
            member%i_floor = i-1                    ! record the number of the element that pierces the seabed
            member%doEndBuoyancy = .true.
            exit
         else if ( EqualRealNos(Za, -InitInp%WaveField%EffWtrDpth ) ) then
            member%doEndBuoyancy = .true.
         end if
      end do
   else
      member%i_floor = 0 ! lower end is at or above the seabed
   end if
  
   ! calculate element-level values
   do i = 1, member%NElements
      member%dSadl_mg(  i) = (member%SaMG( i+1) - member%SaMG( i))/dl
      member%dSadl_in(  i) = (member%Sain( i+1) - member%Sain( i))/dl
      member%dSadl_mg_b(i) = (member%SaMGB(i+1) - member%SaMGB(i))/dl

      member%dSbdl_mg(  i) = (member%SbMG( i+1) - member%SbMG( i))/dl
      member%dSbdl_in(  i) = (member%Sbin( i+1) - member%Sbin( i))/dl
      member%dSbdl_mg_b(i) = (member%SbMGB(i+1) - member%SbMGB(i))/dl
      
      member%alpha(   i) = GetAlphaRec(member%SaMGB(i), member%SaMGB(i+1), member%SbMGB(i), member%SbMGB(i+1))   ! Only used to distribute external buoyancy load to nodes
      member%alpha_fb(i) = GetAlphaRec(member%Sain( i), member%Sain( i+1), member%Sbin( i), member%Sbin( i+1))
      
   end do

   member%Vinner    = 0.0_ReKi  ! Total volume of member without marine growth
   member%Vouter    = 0.0_ReKi  ! Total outer volume of member including marine growth
   member%Vballast  = 0.0_ReKi  ! Total ballasted volume of member
   member%elem_fill = 0         ! Last (partially) filled element of the member
   member%h_fill    = 0.0       ! Axial length of elem_fill occupied by water ballast
   
   ! force-related constants for each element
   do i = 1, member%NElements
   
      Za = InitInp%Nodes(member%NodeIndx(  i))%Position(3)   ! z location of node i
      Zb = InitInp%Nodes(member%NodeIndx(i+1))%Position(3)   ! z location of node i+1
      
      ! ------------------ marine growth weight and inertia ------------------------------------------------
      Vinner_l   = 0.0
      Vouter_l   = 0.0
      Vinner_u   = 0.0
      Vouter_u   = 0.0
      if (i > member%i_floor) then         
         ! full marine growth: get the properties for each half-element lumped to the appropriate node
                  
         SaMid   = 0.5*(member%Sa(  i)+member%Sa(  i+1))  ! length of Side A at middle of segment, where division occurs
         SaMidMG = 0.5*(member%SaMG(i)+member%SaMG(i+1))  ! length of Side A with marine growth at middle of segment, where division occurs
         SbMid   = 0.5*(member%Sb(  i)+member%Sb(  i+1))  ! length of Side B at middle of segment, where division occurs
         SbMidMG = 0.5*(member%SbMG(i)+member%SbMG(i+1))  ! length of Side B with marine growth at middle of segment, where division occurs
         Lmid    = 0.5*dl   ! half-length of segment

         CALL MarineGrowthPartSegmentRec(member%Sa(i  ), SaMid, member%Sb(i  ), SbMid, member%SaMG(i  ), SaMidMG, member%SbMG(i  ), SbMidMG, Lmid, member%MGDensity(i), Vinner_l, Vouter_l, member%m_mg_l(i), member%h_cmg_l(i), member%I_lmg_l(i), member%I_xmg_l(i), member%I_ymg_l(i))   ! get precomputed quantities for lower half-segment
         CALL MarineGrowthPartSegmentRec(member%Sa(i+1), SaMid, member%Sb(i+1), SbMid, member%SaMG(i+1), SaMidMG, member%SbMG(i+1), SbMidMG,-Lmid, member%MGDensity(i), Vinner_u, Vouter_u, member%m_mg_u(i), member%h_cmg_u(i), member%I_lmg_u(i), member%I_xmg_u(i), member%I_ymg_u(i))   ! get precomputed quantities for upper half-segment
         
      else if (i == member%i_floor) then         
         ! crossing seabed: get the properties for part-element above the seabed and lump to the upper node      

         SaMid   = (-member%h_floor*member%Sa(  i) +(dl+member%h_floor)*member%Sa(  i+1))/dl
         SaMidMG = (-member%h_floor*member%SaMG(i) +(dl+member%h_floor)*member%SaMG(i+1))/dl
         SbMid   = (-member%h_floor*member%Sb(  i) +(dl+member%h_floor)*member%Sb(  i+1))/dl
         SbMidMG = (-member%h_floor*member%SbMG(i) +(dl+member%h_floor)*member%SbMG(i+1))/dl
         Lmid    = -member%h_floor

         CALL MarineGrowthPartSegmentRec(member%Sa(i+1), SaMid, member%Sb(i+1), SbMid, member%SaMG(i+1), SaMidMG, member%SbMG(i+1), SbMidMG,-Lmid, member%MGDensity(i), Vinner_u, Vouter_u, member%m_mg_u(i), member%h_cmg_u(i), member%I_lmg_u(i), member%I_xmg_u(i), member%I_ymg_u(i))   ! get precomputed quantities for upper half-segment
         Vinner_l   = 0.0
         Vouter_l   = 0.0
      end if

      ! ------------------ flooded ballast inertia ---------------------------------------------------------
      Vballast_l = 0.0
      Vballast_U = 0.0
      if (member%memfloodstatus > 0 .and. (member%FillFSLoc > Za)) then
         ! Fully filled element, so split in middle
         if ((i > member%i_floor) .and. (member%FillFSLoc >= Zb)) then

            ! get the properties for each half-element lumped to the appropriate node
            SaMidIn = 0.5*(member%Sain(i)+member%Sain(i+1))  ! length of side A of member interior at middle of segment, where division occurs
            SbMidIn = 0.5*(member%Sbin(i)+member%Sbin(i+1))  ! length of side B of member interior at middle of segment, where division occurs
            Lmid   = 0.5*dl   ! half-length of segment
            CALL FloodedBallastPartSegmentRec(member%Sain(i  ), SaMidIn, member%Sbin(i  ), SbMidIn, Lmid, member%FillDens, Vballast_l, member%m_fb_l(i), member%h_cfb_l(i), member%I_lfb_l(i), member%I_xfb_l(i), member%I_yfb_l(i))   ! get precomputed quantities for lower half-segment
            CALL FloodedBallastPartSegmentRec(member%Sain(i+1), SaMidIn, member%Sbin(i+1), SbMidIn,-Lmid, member%FillDens, Vballast_u, member%m_fb_u(i), member%h_cfb_u(i), member%I_lfb_u(i), member%I_xfb_u(i), member%I_yfb_u(i))   ! get precomputed quantities for upper half-segment
 
         ! partially filled element, so split at FillFSLoc
         else if ((i > member%i_floor)  .AND. (member%FillFSLoc < Zb)) then

            ! get the properties for each partial-element lumped to the appropriate node
            Lmid    = member%FillFSLoc - Za 
            SaMidIn = member%Sain(i)+(Lmid/(Zb-Za))*(member%Sain(i+1)-member%Sain(i))  ! length of side A of member interior at middle of segment, where division occurs
            SbMidIn = member%Sbin(i)+(Lmid/(Zb-Za))*(member%Sbin(i+1)-member%Sbin(i))  ! length of side A of member interior at middle of segment, where division occurs
            CALL FloodedBallastPartSegmentRec(member%Sain(i  ), SaMidIn, member%Sbin(i  ), SbMidIn, Lmid, member%FillDens, Vballast_l, member%m_fb_l(i), member%h_cfb_l(i), member%I_lfb_l(i), member%I_xfb_l(i), member%I_yfb_l(i))   ! get precomputed quantities for lower half-segment
            CALL FloodedBallastPartSegmentRec(member%Sain(i+1), SaMidIn, member%Sbin(i+1), SbMidIn,-Lmid,        0.0_ReKi, Vballast_u, member%m_fb_u(i), member%h_cfb_u(i), member%I_lfb_u(i), member%I_xfb_u(i), member%I_yfb_u(i))   ! get precomputed quantities for upper half-segment
 
         else if (i == member%i_floor) then     ! Hopefully we don't have a partially filled element crossing the seabed.
 
            ! crossing seabed: get the properties for part-element above the seabed and lump to the upper node
            SaMidIn = (-member%h_floor*member%Sain(i) +(dl+member%h_floor)*member%Sain(i+1))/dl
            SbMidIn = (-member%h_floor*member%Sbin(i) +(dl+member%h_floor)*member%Sbin(i+1))/dl
            Lmid    = -member%h_floor
            CALL FloodedBallastPartSegmentRec(member%Sain(i+1), SaMidIn, member%Sbin(i+1), SbMidIn,-Lmid, member%FillDens, Vballast_u, member%m_fb_u(i), member%h_cfb_u(i), member%I_lfb_u(i), member%I_xfb_u(i), member%I_yfb_u(i))   ! get precomputed quantities for upper half-segment
            Vballast_l = 0.0
 
         end if
      else  ! Either no ballast flooding in member, or this particular element isn't flooded at all
         Vballast_u        = 0.0
         Vballast_l        = 0.0
         member%m_fb_u(i)  = 0.0
         member%h_cfb_u(i) = 0.0
         member%I_lfb_u(i) = 0.0
         member%I_xfb_u(i) = 0.0
         member%I_yfb_u(i) = 0.0
      endif
      
      ! Determine volumes to add to Non-WAMIT modeled members, etc.
      if (.not. member%PropPot) then
         
         if (Zb < -InitInp%WaveField%EffWtrDpth) then
            ! fully buried element, do not add these volume contributions to totals
         else if (0.0 >= Zb) then   ! Bug fix per OpenFAST issue #844   GJH 2/3/2022
            ! fully submerged elements.  
            ! NOTE: For an element which is fractionaly in the seabed, the entire element volume is added to totals
            member%Vinner = member%Vinner + Vinner_l + Vinner_u
            member%Vouter = member%Vouter + Vouter_l + Vouter_u
            member%Vsubmerged = member%Vsubmerged + Vouter_l + Vouter_u
         else if ((0.0 > Za) .AND. (0.0 < Zb)) then ! Bug fix per OpenFAST issue #844   GJH 2/3/2022
            ! if (i == 1) then
            !    call SetErrStat(ErrID_Fatal, 'The lowest element of a member must not cross the free surface.  This is true for MemberID '//trim(num2lstr(member%MemberID)), errStat, errMsg, RoutineName)
            ! end if
            
            ! partially submerged element
            member%Vinner = member%Vinner + Vinner_l + Vinner_u
            member%Vouter = member%Vouter + Vouter_l + Vouter_u
            ! compute volume portion which is submerged
            Lmid = -Za/cosPhi 
            call RecTaperCalc( member%SaMG(i), member%SaMG(i)+Lmid*member%dSadl_mg(i), member%SbMG(i), member%SbMG(i)+Lmid*member%dSbdl_mg(i), Lmid, Vouter_l, h_c)
            member%Vsubmerged = member%Vsubmerged + Vouter_l
            
         else ! fully above the water
            member%Vinner = member%Vinner + Vinner_l + Vinner_u
            member%Vouter = member%Vouter + Vouter_l + Vouter_u
         end if 
      end if
      
      ! ------------------ flooded ballast weight (done) --------------------
      ! NOTE: this section of code is somewhat redundant with "flooded ballast inertia" section above

      li = dl*(i-1)

      if (Zb < -InitInp%WaveField%EffWtrDpth) then                                                                    ! Fully buried element

         member%floodstatus(i) = 0

      else if (member%memfloodstatus > 0 .and. member%FillFSLoc >= Zb) then                                           ! Fully flooded elements

         member%floodstatus(i) = 1
         if ( EqualRealNos(member%FillFSLoc, Zb) .or. (i==member%NElements) ) then  ! No partially filled elements
            member%elem_fill = i
            member%h_fill    = dl
         end if
         member%Vballast = member%Vballast + Vballast_l + Vballast_u

      else if ((member%memfloodstatus > 0) .and. (member%FillFSLoc > Za) .AND. (member%FillFSLoc < Zb)) then          ! Partially flooded element
         
         member%floodstatus(i) = 2
         member%elem_fill      = i
         member%h_fill         = member%l_fill - (i-1)*dl
         call RecTaperCalc( member%Sain(i), member%Sain(i)+member%h_fill*member%dSadl_in(i), member%Sbin(i), member%Sbin(i)+member%h_fill*member%dSbdl_in(i), member%h_fill, Vballast_l, h_c)
         member%Vballast = member%Vballast + Vballast_l  ! Note: Vballast_l will match calculations above

      else                                                                                                            ! Unflooded element

         member%floodstatus(i) = 0

      end if
      
   end do ! end looping through elements   

end subroutine SetMemberProperties_Rec

!----------------------------------------------------------------------------------------------------------------------------------
subroutine SetupMembers( InitInp, p, m, errStat, errMsg )
   type(Morison_InitInputType),  intent (inout)  :: InitInp
   type(Morison_ParameterType),  intent (inout)  :: p
   type(Morison_MiscVarType),    intent (inout)  :: m
   integer(IntKi),               intent (  out)  :: errStat              ! returns a non-zero value when an error occurs            
   character(*),                 intent (  out)  :: errMsg               ! Error message if errStat /= ErrID_None

   integer(IntKi) :: i, prop1Indx, prop2Indx
   integer(IntKi) :: errStat2              ! returns a non-zero value when an error occurs            
   CHARACTER(errMsgLen)                        :: errMsg2     ! Error message if errStat2 /= ErrID_None
   logical       :: doSwap
   
   
   errStat = ErrID_None
   errMSg  = ''   
   
   ! allocate and copy in the InpMembers array
   p%NMembers = InitInp%NMembers
   ALLOCATE ( p%Members(p%NMembers), STAT = errStat2 )
   IF ( errStat2 /= 0 ) THEN
      errMsg  = ' Error allocating space for the members array.'
      errStat = ErrID_Fatal
      RETURN
   END IF   
   
   ALLOCATE ( m%MemberLoads(p%NMembers), STAT = errStat2 )
   IF ( errStat2 /= 0 ) THEN
      errMsg  = ' Error allocating space for the memberLoads array.'
      errStat = ErrID_Fatal
      RETURN
   END IF  
        
   do i = 1, p%NMembers
      p%Members(i)%MemberID    = InitInp%InpMembers(i)%MemberID
      p%Members(i)%RefLength   = InitInp%InpMembers(i)%RefLength
      p%Members(i)%dl          = InitInp%InpMembers(i)%dl
      p%Members(i)%NElements   = InitInp%InpMembers(i)%NElements
      p%Members(i)%PropPot     = InitInp%InpMembers(i)%PropPot
      p%Members(i)%MHstLMod    = InitInp%InpMembers(i)%MHstLMod
      p%Members(i)%MSecGeom    = InitInp%InpMembers(i)%MSecGeom
      p%Members(i)%MSpinOrient = InitInp%InpMembers(i)%MSpinOrient
      
      call AllocateMemberDataArrays(p%Members(i), m%MemberLoads(i), errStat2, errMsg2)
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'SetupMembers')
      if (ErrStat >= AbortErrLev) return
        
      p%Members(i)%NodeIndx  = InitInp%InpMembers(i)%NodeIndx ! now that the parameter version is allocated, copy the data from the InitInp version
      
      ! only reorder the nodes if the end nodes do not follow the necessary coordinate ordering rules
      call FlipMemberNodeData(p%Members(i), InitInp%nodes, doSwap)
      if (doSwap) then
            prop2Indx = InitInp%InpMembers(I)%MPropSetID1Indx
            prop1Indx = InitInp%InpMembers(I)%MPropSetID2Indx
      else
            prop1Indx = InitInp%InpMembers(I)%MPropSetID1Indx
            prop2Indx = InitInp%InpMembers(I)%MPropSetID2Indx
      end if
      ! Now populate the various member data arrays using the HydroDyn input file data
      if (p%Members(i)%MSecGeom == MSecGeom_Cyl) then
         call SetMemberProperties_Cyl( InitInp%Gravity, p%Members(i), InitInp%InpMembers(i)%MCoefMod, InitInp%InpMembers(i)%MmbrCoefIDIndx, InitInp%InpMembers(i)%MmbrFilledIDIndx, InitInp%MPropSetsCyl(prop1Indx), InitInp%MPropSetsCyl(prop2Indx), InitInp, errStat2, errMsg2 ) 
      else if (p%Members(i)%MSecGeom == MSecGeom_Rec) then
         call SetMemberProperties_Rec( InitInp%Gravity, p%Members(i), InitInp%InpMembers(i)%MCoefMod, InitInp%InpMembers(i)%MmbrCoefIDIndx, InitInp%InpMembers(i)%MmbrFilledIDIndx, InitInp%MPropSetsRec(prop1Indx), InitInp%MPropSetsRec(prop2Indx), InitInp, errStat2, errMsg2 ) 
      end if
      call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'SetupMembers')
      if (ErrStat >= AbortErrLev) return
   end do
      
end subroutine SetupMembers
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps. 
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
!! A lot of the model setup has been done previously in Morison_ProcessMorisonGeometry, and stored in InitInp.
SUBROUTINE Morison_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, errStat, errMsg )
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
   REAL(DbKi),                        INTENT(IN   )  :: Interval    !< Coupling interval in seconds: the rate that 
                                                                     !!   (1) Morison_UpdateStates() is called in loose coupling &
                                                                     !!   (2) Morison_UpdateDiscState() is called in tight coupling.
                                                                     !!   Input is the suggested time from the glue code; 
                                                                     !!   Output is the actual coupling interval that will be used 
                                                                     !!   by the glue code.
   TYPE(Morison_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
   INTEGER(IntKi),                    INTENT(  OUT)  :: errStat     !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: errMsg      !< Error message if errStat /= ErrID_None
   
   character(*), parameter                           :: RoutineName = 'Morison_Init'

   TYPE(Morison_MemberType) :: member      ! the current member
   INTEGER                  :: i, j, im
   REAL(ReKi)               :: v2D(3,1), pos(3)
   real(ReKi)               :: An(3), An_drag(3), Vn(3), I_n(3), sgn, Amag, Amag_drag, Vmag, Imag, Ir_MG_end, Il_MG_end, R_I(3,3), IRl_mat(3,3), tMG, MGdens
   integer(IntKi)           :: MemberEndIndx
   INTEGER, ALLOCATABLE     :: commonNodeLst(:)
   LOGICAL, ALLOCATABLE     :: usedJointList(:)
   integer(IntKi)           :: errStat2    ! returns a non-zero value when an error occurs            
   CHARACTER(errMsgLen)     :: errMsg2     ! Error message if errStat2 /= ErrID_None
   
   ! Initialize errStat        
   errStat = ErrID_None         
   errMsg  = ""               
   
   ! Define parameters here:  
   p%DT         = Interval
   p%Gravity    = InitInp%Gravity
   p%NNodes     = InitInp%NNodes
   p%NJoints    = InitInp%NJoints
   p%NumOuts    = InitInp%NumOuts
   p%NMOutputs  = InitInp%NMOutputs                       ! Number of members to output [ >=0 and <10]
   p%WaveDisp   = InitInp%WaveDisp
   p%AMMod      = InitInp%AMMod
   p%VisMeshes  = InitInp%VisMeshes                       ! visualization mesh for morison elements
   p%PtfmYMod   = InitInp%PtfmYMod
   p%NFillGroups = InitInp%NFillGroups

   ! Pointer to SeaState WaveField
   p%WaveField => InitInp%WaveField
   
   ! Only compute added-mass force up to the free surface if wave stretching is enabled
   IF ( p%WaveField%WaveStMod .EQ. 0_IntKi ) THEN
       ! Setting AMMod to zero just in case. Probably redundant.
       p%AMMod = 0_IntKi
   END IF


   ALLOCATE ( p%MOutLst(p%NMOutputs), STAT = errStat2 )
   IF ( errStat2 /= 0 ) THEN
      call SetErrStat(ErrID_Fatal,'Error allocating space for MOutLst array.', ErrStat, ErrMsg, RoutineName)
      RETURN
   END IF
   IF (ALLOCATED(InitInp%MOutLst) ) then
      do i=1,size(InitInp%MOutLst)
         call  Morison_CopyMOutput( InitInp%MOutLst(i), p%MOutLst(i), MESH_NEWCOPY, ErrStat2, ErrMsg2 )                 ! Member output data
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'Morison_Init' )
      end do
   end if
      
   p%NJOutputs = InitInp%NJOutputs                        ! Number of joints to output [ >=0 and <10]
      
   ALLOCATE ( p%JOutLst(p%NJOutputs), STAT = errStat2 )
   IF ( errStat2 /= 0 ) THEN
      call SetErrStat(ErrID_Fatal,'Error allocating space for JOutLst array.', ErrStat, ErrMsg, RoutineName)
      RETURN
   END IF
   IF (ALLOCATED(InitInp%JOutLst) ) &
      p%JOutLst =    InitInp%JOutLst            ! Joint output data
 
   ! ----------------------- set up the members -----------------------
   call SetupMembers( InitInp, p, m, errStat2, errMsg2 ) 
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if ( errStat >= AbortErrLev ) return
   
   !------------------------ set up joint (or joint-node) properties --
   do i = 1, InitInp%NJoints
      InitInp%Nodes(i)%JAxCd = InitInp%AxialCoefs(InitInp%InpJoints(i)%JointAxIDIndx)%AxCd
      InitInp%Nodes(i)%JAxCa = InitInp%AxialCoefs(InitInp%InpJoints(i)%JointAxIDIndx)%AxCa
      InitInp%Nodes(i)%JAxCp = InitInp%AxialCoefs(InitInp%InpJoints(i)%JointAxIDIndx)%AxCp     
      InitInp%Nodes(i)%JAxFDMod   = InitInp%AxialCoefs(InitInp%InpJoints(i)%JointAxIDIndx)%AxFDMod
      InitInp%Nodes(i)%JAxVnCOff  = InitInp%AxialCoefs(InitInp%InpJoints(i)%JointAxIDIndx)%AxVnCOff
      InitInp%Nodes(i)%JAxFDLoFSc = InitInp%AxialCoefs(InitInp%InpJoints(i)%JointAxIDIndx)%AxFDLoFSc
  
      ! Redundant work (these are already assigned to the member data arrays, 
      ! but is needed on the joint data because we report the tMG, and MGDensity at each Joint node in the Summary File
      call SetNodeMG( InitInp%NMGDepths, InitInp%MGDepths, InitInp%Nodes(i), p%WaveField%MSL2SWL, InitInp%Nodes(i)%tMG, InitInp%Nodes(i)%MGDensity )
   end do

   ! allocate and copy in node-based load and hydrodynamic arrays
   call AllocateNodeLoadVariables(InitInp, p, m, p%NNodes, errStat2, errMsg2 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if ( errStat >= AbortErrLev ) return
   
   ! Create the input and output meshes associated with loads at the nodes     
   CALL MeshCreate( BlankMesh      = u%Mesh          &
                     ,IOS          = COMPONENT_INPUT        &
                     ,Nnodes       = p%NNodes      &
                     ,errStat      = errStat                &
                     ,ErrMess      = errMsg2                &
                     ,TranslationDisp = .TRUE.              &
                     ,Orientation     = .TRUE.              &
                     ,TranslationVel  = .TRUE.              &
                     ,RotationVel     = .TRUE.              &
                     ,TranslationAcc  = .TRUE.              &
                     ,RotationAcc     = .TRUE.               )

   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if ( errStat >= AbortErrLev ) return

!TODO: Do we still need this for visualization?  How is it used? GJH 3/26/2020 Actually need a line mesh to properly visualize the members
   !CALL AllocAry( Morison_Rad, numDistribMarkers, 'Morison_Rad', errStat, errMsg)
   !IF ( errStat >= AbortErrLev ) RETURN
   
   
   DO I=1,p%NNodes
   ! This needs to change so that the Position is relative to MSL NOT SWL:
      pos = InitInp%Nodes(I)%Position
      pos(3) = pos(3) + p%WaveField%MSL2SWL
         ! Create the node on the mesh 
      CALL MeshPositionNode (u%Mesh                &
                        , i                        &      
                        , pos                      &  ! this info comes from HydroDyn input file and the subroutine: Morison_GenerateSimulationNodes
                        , errStat2                  &
                        , errMsg2                   &
                        ) !, transpose(p%Nodes(I)%R_LToG)          )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      if ( errStat >= AbortErrLev ) return

!TODO: Do we still need this for visualization?  How is it used? GJH 3/26/2020  Actually need a line mesh to properly visualize the members
     ! Morison_Rad(count) = p%Nodes(I)%R   ! set this for FAST visualization
      
     
   
         ! Create the mesh element
   
      CALL MeshConstructElement (u%Mesh   &
                            , ELEMENT_POINT      &                                  
                            , errStat2            &
                            , errMsg2  &
                            , i                  &
                                        )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         if ( errStat >= AbortErrLev ) return
                            

   END DO

   CALL MeshCommit ( u%Mesh   &
                      , errStat2            &
                      , errMsg2             )
   
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if ( errStat >= AbortErrLev ) return

   
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
                     ,errStat      = errStat2               &
                     ,ErrMess      = errMsg2                &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if ( errStat >= AbortErrLev ) return
   u%Mesh%RemapFlag = .TRUE.
   y%Mesh%RemapFlag = .TRUE.

         ! Define initial system states here:

   x%DummyContState           = 0
   !xd%DummyDiscState          = 0
   ALLOCATE ( xd%V_rel_n_FiltStat(p%NJoints), STAT = ErrStat )
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating space for V_rel_n_FiltStat array.'
      ErrStat = ErrID_Fatal
      RETURN
   END IF
   xd%V_rel_n_FiltStat = 0.0_ReKi

   z%DummyConstrState         = 0
   OtherState%DummyOtherState = 0

   
   ! allocate and initialize joint-specific arrays   
      
   ALLOCATE ( commonNodeLst(10), STAT = errStat2 )
   IF ( errStat2 /= 0 ) THEN
      call SetErrStat(ErrID_Fatal,'Error allocating space for commonNodeLst array.', ErrStat, ErrMsg, RoutineName)
      RETURN
   END IF 
   commonNodeLst = -1
   
   ALLOCATE ( usedJointList(p%NJoints), STAT = errStat2 )
   IF ( errStat2 /= 0 ) THEN
      call SetErrStat(ErrID_Fatal,'Error allocating space for UsedJointList array.', ErrStat, ErrMsg, RoutineName)
      RETURN
   END IF  
   usedJointList = .FALSE.
   
   ! loop through joints to calculate joint quantities (the joints are the first NJoints nodes)

   usedJointList = .FALSE.   
   commonNodeLst = -1
   !TODO: Error Handling
   
   
   do i = 1,p%NJoints      

      An        = 0.0
      Vn        = 0.0
      I_n       = 0.0
      MGdens    = 0.0
      tMG       = -999.0
      An_drag   = 0.0
      
      IF ( InitInp%InpJoints(i)%Position(3) >= -InitInp%WaveField%WtrDpth ) THEN
   
         ! loop through each member attached to the joint, getting the radius of its appropriate end
         DO J = 1, InitInp%InpJoints(I)%NConnections
      
            ! identify attached member and which end to use
            IF (InitInp%InpJoints(I)%ConnectionList(J) > 0) THEN         ! set up for end node 1
               !TODO: Should not perform a copy here?  A pointer to data would be better?
               member = p%Members(InitInp%InpJoints(I)%ConnectionList(J))   
               MemberEndIndx = 1
            ELSE     
               ! set up for end node N+1
               ! NOTE:  %ConnectionList(J) is negative valued if InitInp%Morison%InpMembers(I)%MJointID2 == InitInp%Morison%InpJoints(J)%JointID.  See HydroDynInput_ProcessInitData, members section
               member = p%Members(-InitInp%InpJoints(I)%ConnectionList(J))
               MemberEndIndx = member%NElements + 1
            END IF
      
            ! Compute the signed area*outward facing normal of this member
            sgn = 1.0
            
            IF ( MemberEndIndx == 1 ) THEN
               sgn = -1.0                                ! Local coord sys points into member at starting node, so flip sign of local z vector
            ELSE
               sgn =  1.0                                ! Local coord sys points out of member at ending node, so leave sign of local z vector
            END IF

            ! Account for reordering of what the original node for the end was -- This affects the sign of the An term which can pose a problem for members crossing the waterline
            if (member%Flipped)   sgn = -1.0 * sgn
               
            IF ( member%MSecGeom == MSecGeom_Cyl ) THEN
               ! Compute the signed quantities for this member end (for drag regardless of PropPot value), and add them to the joint values
               An_drag = An_drag + sgn* member%k*Pi*(member%RMG(MemberEndIndx))**2     ! area-weighted normal vector
               ! For the following quantities, the attached member cannot be modeled using WAMIT if we're to count it
               IF  (.NOT. member%PropPot) THEN
                  ! Compute the signed quantities for this member end, and add them to the joint values
                  An = An + sgn* member%k*Pi*(member%RMG(MemberEndIndx))**2                 ! area-weighted normal vector
                  Vn = Vn + sgn* member%k*TwoPi/3.0_ReKi*(member%RMG(MemberEndIndx))**3     ! r^3-weighted normal vector used for mass
                  I_n=I_n + sgn* member%k*Pi*(member%RMG(MemberEndIndx))**4                 ! r^4-weighted normal vector used for moments of inertia
               END IF
            ELSE IF ( member%MSecGeom == MSecGeom_Rec ) THEN
               ! Compute the signed quantities for this member end (for drag regardless of PropPot value), and add them to the joint values
               An_drag = An_drag + sgn* member%k*(member%SaMG(MemberEndIndx) * member%SbMG(MemberEndIndx))     ! area-weighted normal vector
               ! For the following quantities, the attached member cannot be modeled using WAMIT if we're to count it
               IF  (.NOT. member%PropPot) THEN
                  ! Compute the signed quantities for this member end, and add them to the joint values
                  An = An + sgn* member%k*(member%SaMG(MemberEndIndx) * member%SbMG(MemberEndIndx))                                       ! area-weighted normal vector
                  Vn = Vn + sgn* member%k*2.0_ReKi/3.0_ReKi*SQRT((member%SaMG(MemberEndIndx)*member%SbMG(MemberEndIndx)))**3/SQRT(Pi)     ! volume-weighted normal vector used for mass
                  ! Don't have a good way to handle marine growth moment of inertia for rectangular endplates. Use an approximation for now, which should be reasonble for nearly square members
                  I_n=I_n + sgn* member%k* ( (member%SaMG(MemberEndIndx)**3 * member%SbMG(MemberEndIndx)) + (member%SaMG(MemberEndIndx) * member%SbMG(MemberEndIndx)**3) ) /6.0_ReKi
               END IF
            ELSE
               call SetErrStat( ErrId_Fatal, " Unrecognized member cross section geometry ", errStat, errMsg, RoutineName ); return
            END IF

            IF (tMG == -999.0) THEN
               ! All member nodes at this joint will have the same MG thickness and density, so only do this once
               tMG = member%tMG(MemberEndIndx)
               MGdens = member%MGdensity(MemberEndIndx) 
            END IF

         END DO   !J = 1, InitInp%InpJoints(I)%NConnections
         
         p%An_End(:,i) = An_drag 
         Amag_drag = Dot_Product(An_drag ,An_drag)
         Amag = Dot_Product(An ,An)
         IF (EqualRealNos(Amag_drag, 0.0_ReKi)) THEN
            p%DragConst_End(i) =  0.0
         ELSE
            p%DragConst_End(i) = InitInp%Nodes(i)%JAxCd*p%WaveField%WtrDens / ( 4.0_ReKi * Amag_drag )
         END IF
         ! magnitudes of normal-weighted values
         Amag = sqrt(Amag)
         Vmag = sqrt(Dot_Product(Vn ,Vn))
         Imag = sqrt(Dot_Product(I_n,I_n))
      
         ! Constant part of the external hydrodynamic added mass term
         if ( Vmag > 0.0 ) then
            v2D(:,1) = Vn        
            p%AM_End(:,:,i) = (InitInp%Nodes(I)%JAxCa*p%WaveField%WtrDens/ Vmag)*matmul(v2D, transpose(v2D))
         end if
         
         ! Constant part of the external hydrodynamic dynamic pressure force
         if ( Amag > 0.0 ) then
            p%DP_Const_End(:,i) = -InitInp%Nodes(i)%JAxCp*An 
         endif
         
         ! marine growth mass/inertia magnitudes
         p%Mass_MG_End(i) =  MGdens * tMG * Amag
         p%F_WMG_End(3,i) = -MGdens * tMG * Amag * InitInp%Gravity  ! Z component of the directional force due to marine growth mass at joint
         Ir_MG_end   =  0.25 * MGdens * tMG * Imag  ! radial moment of inertia magnitude
         Il_MG_end   =  0.5  * MGdens * tMG * Imag  ! axial moment of inertia magnitude
      
         ! get rotation matrix for moment of inertia orientations
         call RodrigMat(I_n, R_I, errStat, errMsg)
         IF ( errStat >= AbortErrLev ) RETURN

         ! globally-oriented moment of inertia matrix for joint
         Irl_mat = 0.0
         Irl_mat(1,1) = Ir_MG_end
         Irl_mat(2,2) = Ir_MG_end
         Irl_mat(3,3) = Il_MG_end
      
         p%I_MG_End(:,:,i) = MatMul( MatMul(R_I, Irl_mat), Transpose(R_I) ) ! final moment of inertia matrix for node
         

      END IF  ! InitInp%InpJoints(i)%Position(3) >= -WtrDpth
   
      p%DragMod_End   (i) = InitInp%Nodes(i)%JAxFDMod
      IF ( InitInp%Nodes(i)%JAxVnCOff .LE. 0.0_ReKi) THEN
         p%VRelNFiltConst(i) = 1.0_ReKi
         p%DragLoFSc_End (i) = 1.0_ReKi
      ELSE
         p%VRelNFiltConst(i) = exp(-2.0*Pi*InitInp%Nodes(i)%JAxVnCOff * p%DT)
         p%DragLoFSc_End (i) = InitInp%Nodes(i)%JAxFDLoFSc
      END IF

   END DO ! looping through nodes that are joints, i
          
   ! Copy ballast group information to parameters
   ALLOCATE ( p%FilledGroups(p%NFillGroups), STAT = ErrStat2 )
      IF ( ErrStat2 /= 0 ) THEN
         call SetErrStat(ErrID_Fatal,'Error allocating space for FilledGroups array.',errStat,errMsg,RoutineName ); return
      END IF
   DO i=1,p%NFillGroups
      CALL Morison_CopyFilledGroupType(InitInp%FilledGroups(i),p%FilledGroups(i),0,ErrStat2,ErrMsg2)
      IF ( ErrStat2 /= 0 ) THEN
         call SetErrStat(ErrID_Fatal,'Error copying FilledGroups array.',errStat,errMsg,RoutineName ); return
      END IF
   END DO
   ! Determine if the filled group is open to the environment through the open end of members buried in the seabed
   DO i = 1,p%NFillGroups
      p%FilledGroups(i)%IsOpen = .false.
      DO j = 1,p%FilledGroups(i)%FillNumM
         im = p%FilledGroups(i)%FillMList(j)
         IF ( p%Members(im)%i_floor > 0 ) THEN
            p%FilledGroups(i)%IsOpen = .true.
            IF ( p%FilledGroups(i)%FillFSLoc > p%WaveField%MSL2SWL ) THEN
               call SetErrStat(ErrID_Fatal,' FillFSLoc cannot be higher than MSL2SWL if FillMList contains any member that is fully or partially buried in the seabed. ',errStat,errMsg,RoutineName ); return
            END IF
            IF ( .not. EqualRealNos( p%FilledGroups(i)%FillDens, p%WaveField%WtrDens) ) THEN
               call SetErrStat(ErrID_Fatal,' FillDens must be the same as the external water density if FillMList contains any member that is fully or partially buried in the seabed. ',errStat,errMsg,RoutineName ); return
            END IF
            EXIT
         END IF
      END DO
   END DO

         ! Define initial guess for the system inputs here:
         !    u%DummyInput = 0
         ! Define system output initializations (set up mesh) here:  
         ! Define initialization-routine output here:
         
   ! Initialize the outputs      
   CALL MrsnOUT_Init( InitInp, y, p, InitOut, errStat2, errMsg2 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if ( errStat >= AbortErrLev ) return
   
   ! visualization Line2 mesh
   if (p%VisMeshes) then
      call VisMeshSetup(u,p,y,m,InitOut,ErrStat2,ErrMsg2); call SetErrStat( ErrStat2, ErrMsg2, errStat, errMsg, 'Morison_Init' )
      if ( errStat >= AbortErrLev ) return
   endif

   ! We will call CalcOutput to compute the loads for the initial reference position
   ! Then we can use the computed load components in the Summary File
   ! NOTE: Morison module has no states, otherwise we could no do this. GJH
   
   call Morison_CalcOutput(0.0_DbKi, u, p, x, xd, z, OtherState, y, m, errStat2, errMsg2 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if ( errStat >= AbortErrLev ) return
   
      ! Write Summary information to *HydroDyn* summary file now that everything has been initialized. 
   CALL WriteSummaryFile( InitInp%UnSum, InitInp%NJoints, InitInp%NNodes, InitInp%Nodes, p%NMembers, p%Members, &
                          p%NumOuts, p%OutParam, p%MOutLst, p%JOutLst, u%Mesh, y%Mesh, p, m, errStat2, errMsg2 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if ( errStat >= AbortErrLev ) return
                                                       
   !Contains:
   !   SUBROUTINE CleanUpInitOnErr
   !   IF (ALLOCATED(sw(1)%array))  DEALLOCATE(sw(1)%array, STAT=aviFail)
   !   END SUBROUTINE

END SUBROUTINE Morison_Init


!----------------------------------------------------------------------------------------------------------------------------------
subroutine VisMeshSetup(u,p,y,m,InitOut,ErrStat,ErrMsg)
   type(Morison_InputType),      intent(inout)  :: u
   type(Morison_ParameterType),  intent(in   )  :: p
   type(Morison_OutputType),     intent(inout)  :: y
   type(Morison_MiscVarType),    intent(inout)  :: m
   type(Morison_InitOutputType), intent(inout)  :: InitOut
   integer(IntKi),               intent(  out)  :: ErrStat
   character(*),                 intent(  out)  :: ErrMsg

   integer(IntKi)          :: TotNodes                ! total nodes in all elements (may differ from p%NNodes due to overlaps)
   integer(IntKi)          :: TotElems                ! total number of elements
   integer(IntKi)          :: NdIdx, iMem, iNd, NdNum ! indexing
   real(ReKi)              :: NdPos(3),Pos1(3),Pos2(3)
   real(R8Ki)              :: MemberOrient(3,3)
   real(R8Ki)              :: Theta(3)                ! Euler rotations
   integer(IntKi)          :: ErrStat2
   character(ErrMsgLen)    :: ErrMsg2
   character(*), parameter :: RoutineName = 'VisMeshSetup'

   ErrStat = ErrID_None
   ErrMsg  = ""

   ! Total number of nodes = sum of all member nodes
   ! Total number of elements = sum of all member elements
   TotNodes=0
   TotElems=0
   do iMem=1,size(p%Members)
      TotElems = TotElems + p%Members(iMem)%NElements
      TotNodes = TotNodes + size(p%Members(iMem)%NodeIndx)
   enddo

   ! Storage for the radius associated with each node
   call AllocAry( InitOut%MorisonVisRad, TotNodes, 'MorisonVisRad', ErrStat2, ErrMsg2)
   if (Failed())  return

   call MeshCreate( BlankMesh    = y%VisMesh,         &
                    IOS          = COMPONENT_OUTPUT,  &
                    Nnodes       = TotNodes,          &
                    ErrStat      = ErrStat2,          &
                    ErrMess      = ErrMsg2,           &
                    TranslationDisp = .TRUE.,         &
                    Orientation     = .TRUE.          )
   if (Failed())  return

   ! Position the nodes
   NdNum=0   ! node number in y%VisMesh
   do iMem=1,size(p%Members)

!FIXME:MemberOrient This is not correct for non-circular or curved members
      ! calculate an orientation using yaw-pitch-roll sequence with roll defined as zero (insufficient info)
      Pos1=u%Mesh%Position(:,p%Members(iMem)%NodeIndx(1))                              ! start node position of member
      Pos2=u%Mesh%Position(:,p%Members(iMem)%NodeIndx(size(p%Members(iMem)%NodeIndx))) ! end   node position of member
      Theta(1) = 0.0_R8Ki                                                        ! roll (assumed since insufficient info)
      Theta(2) = acos(real((Pos2(3)-Pos1(3))/TwoNorm(Pos2-Pos1),R8Ki))           ! pitch
      Theta(3) = atan2(real(Pos2(2)-Pos1(2),R8Ki),real(Pos2(1)-Pos1(1),R8Ki))    ! yaw
      MemberOrient=EulerConstructZYX(Theta)  ! yaw-pitch-roll sequence

      ! Set mesh postion, orientation, and radius
      do iNd=1,size(p%Members(iMem)%NodeIndx)
         NdNum=NdNum+1                             ! node number in y%VisMesh
         NdIdx = p%Members(iMem)%NodeIndx(iNd)     ! node number in u%Mesh
         NdPos = u%Mesh%Position(:,NdIdx)          ! node position
         call MeshPositionNode (y%VisMesh, NdNum, u%Mesh%Position(:,NdIdx), ErrStat2,  ErrMsg2, Orient=MemberOrient)
         if (Failed())  return
         InitOut%MorisonVisRad(NdNum) = p%Members(iMem)%RMG(iNd)   ! radius (including marine growth) for visualization
      enddo
   enddo

   ! make elements (line nodes start at 0 index, so N+1 total nodes)
   NdNum=0   ! node number in y%VisMesh
   do iMem=1,size(p%Members)
      do iNd=1,size(p%Members(iMem)%NodeIndx)
         NdNum=NdNum+1                             ! node number in y%VisMesh
         if (iNd==1) cycle
         call MeshConstructElement ( Mesh      = y%VisMesh,    &
                                    Xelement = ELEMENT_LINE2,  &
                                    P1=NdNum-1, P2=NdNum,      &  ! nodes to connect
                                    errStat      = ErrStat2,   &
                                    ErrMess      = ErrMsg2     )
         if (Failed())  return
      enddo
   enddo

   ! commit the assembled mesh
   call MeshCommit ( y%VisMesh, ErrStat2, ErrMsg2)
   if (Failed())  return

   ! map the mesh to u%Mesh
   call MeshMapCreate( u%Mesh, y%VisMesh, m%VisMeshMap, ErrStat2, ErrMsg2 )
   if (Failed())  return

contains
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      !if (Failed) then
      !   call FailCleanup()
      !endif
   end function Failed
end subroutine VisMeshSetup


SUBROUTINE RodrigMat(a, R, errStat, errMsg)
   ! calculates rotation matrix R to rotate unit vertical vector to direction of input vector a
   
   REAL(ReKi),      INTENT ( IN    )  :: a(3)    ! input vector
   REAL(ReKi),      INTENT ( INOUT )  :: R(3,3)  ! rotation matrix from Rodrigues's rotation formula
   INTEGER(IntKi),  INTENT(  OUT)     :: errStat     ! Error status of the operation
   CHARACTER(*),    INTENT(  OUT)     :: errMsg      ! Error message if errStat /= ErrID_None

   REAL(ReKi)                         :: vec(3)  ! scaled and adjusted input vector
   REAL(ReKi)                         :: factor  ! denomenator used for scaling                     
   ErrStat  = ErrID_None
   ErrMsg   = ""
   factor = Dot_Product(a,a)
   ! Return the identity if the vector is zero.  We are defining it this way because of how this is used
   if ( EqualRealNos(factor, 0.0_ReKi) ) then
      CALL EYE(R, errStat,errMsg)
   else IF ( EqualRealNos(a(1), 0.0_ReKi) .AND. EqualRealNos(a(2), 0.0_ReKi) ) THEN    ! return identity if vertical
      CALL EYE(R, errStat,errMsg)
      IF (a(3) < 0) THEN
         R = -R
      END IF   
   else   
      vec = a/SQRT(factor) ! normalize a
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
   end if
   
END SUBROUTINE RodrigMat

!----------------------------------------------------------------------------------------------------------------------------------
FUNCTION GetAlphaCyl(R1,R2)
   ! calculates relative center of volume location for a (tapered) cylindrical element
   real(ReKi)    :: GetAlphaCyl
   REAL(ReKi),                     INTENT    ( IN    )  :: R1  ! interior radius of element at node point
   REAL(ReKi),                     INTENT    ( IN    )  :: R2  ! interior radius of other end of part-element
   
   IF ( EqualRealNos(R1, R2) ) THEN ! Also cover the case where R1=R2=0
      GetAlphaCyl = 0.5
   ELSE
      GetAlphaCyl = (R1*R1 + 2.0*R1*R2 + 3.0*R2*R2)/4.0/(R1*R1 + R1*R2 + R2*R2)
   END IF
   
END FUNCTION GetAlphaCyl

FUNCTION GetAlphaRec(a0,a1,b0,b1)
   ! calculates relative center of volume location for a (tapered) rectangular element
   real(ReKi)    :: GetAlphaRec
   REAL(ReKi),                     INTENT    ( IN    )  :: a0  ! Length of side A of element at node 1
   REAL(ReKi),                     INTENT    ( IN    )  :: a1  ! Length of side A of element at node 2
   REAL(ReKi),                     INTENT    ( IN    )  :: b0  ! Length of side B of element at node 1
   REAL(ReKi),                     INTENT    ( IN    )  :: b1  ! Length of side B of element at node 2
   
   REAL(ReKi) :: tmp
   tmp = 2.0*a0*b0+a0*b1+a1*b0+2.0*a1*b1

   IF ( EqualRealNos(tmp, 0.0_ReKi) ) THEN ! zero volume
      GetAlphaRec = 0.5
   ELSE
      GetAlphaRec = 0.5 * (a0*b0+a0*b1+a1*b0+3.0*a1*b1) / tmp
   END IF
   
END FUNCTION GetAlphaRec

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AllocateNodeLoadVariables(InitInp, p, m, NNodes, errStat, errMsg )
   TYPE(Morison_InitInputType),       INTENT(IN   )  :: InitInp     ! Initialization inputs
   TYPE(Morison_ParameterType),       INTENT(INOUT)  :: p           ! parameter variables
   TYPE(Morison_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables            
   INTEGER(IntKi),                    INTENT(IN   )  :: NNodes      ! number of nodes in node list
   INTEGER(IntKi),                    INTENT(  OUT)  :: errStat     ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: errMsg      ! Error message if errStat /= ErrID_None
   integer(IntKi)                                    :: errStat2    ! Returns a non-zero value when an error occurs            
   CHARACTER(errMsgLen)                              :: errMsg2     ! Error message if errStat2 /= ErrID_None
   character(*), parameter                           :: routineName = 'AllocateNodeLoadVariables'
   
   ! Initialize errStat        
   errStat = ErrID_None         
   errMsg  = ""               
   
   call AllocAry( m%DispNodePosHdn,   3, NNodes   , 'm%DispNodePosHdn', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName) 
   call AllocAry( m%DispNodePosHst,   3, NNodes   , 'm%DispNodePosHst', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName) 
   call AllocAry( m%nodeInWater  ,       NNodes   , 'm%nodeInWater'   , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%vrel         ,    3, NNodes   , 'm%vrel'          , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%FV           ,    3, NNodes   , 'm%FV'            , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%FA           ,    3, NNodes   , 'm%FA'            , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%FAMCF        ,    3, NNodes   , 'm%FAMCF'         , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%FDynP        ,       NNodes   , 'm%FDynP'         , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%WaveElev     ,       NNodes   , 'm%WaveElev'      , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%WaveElev1    ,       NNodes   , 'm%WaveElev1'     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%WaveElev2    ,       NNodes   , 'm%WaveElev2'     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%An_End       ,    3, p%NJoints, 'p%An_End'        , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%DragConst_End,       p%NJoints, 'p%DragConst_End' , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%F_I_End      ,    3, p%NJoints, 'm%F_I_End'       , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%F_BF_End     ,    6, p%NJoints, 'm%F_BF_End'      , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%F_A_End      ,    3, p%NJoints, 'm%F_A_End'       , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%F_D_End      ,    3, p%NJoints, 'm%F_D_End'       , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%F_B_End      ,    6, p%NJoints, 'm%F_B_End'       , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%F_IMG_End    ,    6, p%NJoints, 'm%F_IMG_End'     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%I_MG_End     , 3, 3, p%NJoints, 'p%I_MG_End'      , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%F_WMG_End    ,    3, p%NJoints, 'p%F_WMG_End'     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%Mass_MG_End  ,       p%NJoints, 'p%Mass_MG_End'   , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%AM_End       , 3, 3, p%NJoints, 'p%AM_End'        , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%DP_Const_End ,    3, p%NJoints, 'p%DP_Const_End'  , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%V_rel_n        ,     p%NJoints, 'm%V_rel_n'       , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%V_rel_n_HiPass ,     p%NJoints, 'm%V_rel_n_HiPass', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%zFillGroup   ,   p%NFillGroups, 'm%zFillGroup'    , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%DragMod_End    ,     p%NJoints, 'p%DragMod_End'   , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%DragLoFSc_End  ,     p%NJoints, 'p%DragLoFSc_End' , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%VRelNFiltConst ,     p%NJoints, 'p%VRelNFiltConst', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)

   if (errStat >= AbortErrLev) return
   
   m%DispNodePosHdn = 0.0_ReKi
   m%DispNodePosHst = 0.0_ReKi
   m%nodeInWater    = 0
   m%vrel           = 0.0_ReKi
   m%FV             = 0.0_ReKi
   m%FA             = 0.0_ReKi
   m%FDynP          = 0.0_ReKi
   p%An_End         = 0.0
   p%DragConst_End  = 0.0
   m%F_I_End        = 0.0
   m%F_BF_End       = 0.0
   m%F_A_End        = 0.0
   m%F_D_End        = 0.0
   m%F_B_End        = 0.0
   m%F_IMG_End      = 0.0
   p%DP_Const_End   = 0.0
   p%I_MG_End       = 0.0
   p%Mass_MG_End    = 0.0
   p%F_WMG_End      = 0.0
   p%AM_End         = 0.0
   m%V_rel_n        = 0.0_ReKi
   m%V_rel_n_HiPass = 0.0_ReKi
   
END SUBROUTINE AllocateNodeLoadVariables
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE Morison_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, errStat, errMsg )   
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
   INTEGER(IntKi),                    INTENT(  OUT)  :: errStat     !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: errMsg      !< Error message if errStat /= ErrID_None

   ! Local variables
   INTEGER(IntKi)                                    :: errStat2    ! Error status of the operation (occurs after initial error)
   CHARACTER(errMsgLen)                              :: errMsg2     ! Error message if errStat2 /= ErrID_None
   character(*), parameter                           :: RoutineName = 'Morison_CalcOutput'
      
   REAL(ReKi)                                        :: vmag, vmagf
   INTEGER                                           :: I, J
   REAL(ReKi)                                        :: qdotdot(6)      ! The structural acceleration of a mesh node
      
   TYPE(Morison_MemberType) :: mem     ! the current member
   INTEGER                  :: N       ! Number of elements within a given member
   REAL(ReKi)               :: dl      ! Element length within a given member, m
   REAL(ReKi)               :: vec(3)  ! Vector pointing from a member's 1st node to its last node
   REAL(ReKi)               :: phi, phi1, phi2     ! member tilt angle
   REAL(ReKi)               :: cosPhi, cosPhi1, cosPhi2
   REAL(ReKi)               :: sinPhi, sinPhi1, sinPhi2
   REAL(ReKi)               :: tanPhi
   REAL(ReKi)               :: sinBeta, sinBeta1, sinBeta2
   REAL(ReKi)               :: cosBeta, cosBeta1, cosBeta2
   REAL(ReKi)               :: CMatrix(3,3), CMatrix1(3,3), CMatrix2(3,3), CTrans(3,3) ! Direction cosine matrix for element, and its transpose
   REAL(ReKi)               :: l, z1, z2, zMid, r1, r2, r1b, r2b, r1In, r2In, rMidIn, rn, rn1, rn2, z_hi, zFillGroup
   REAL(ReKi)               :: Sa1, Sa2, Sa1b, Sa2b, SaMidb, Sa1In, Sa2In, SaMidIn
   REAL(ReKi)               :: Sb1, Sb2, Sb1b, Sb2b, SbMidb, Sb1In, Sb2In, SbMidIn
   REAL(ReKi)               :: dRdl_mg,   dSadl_mg,   dSbdl_mg    ! shorthand for taper including marine growth of element i
   REAL(ReKi)               :: dRdl_mg_b, dSadl_mg_b, dSbdl_mg_b  ! shorthand for taper including marine growth of element i with radius scaling by sqrt(Cb)
   REAL(ReKi)               :: RMGFSInt, SaMGFSInt, SbMGFSInt     ! Member radius with marine growth at the intersection with the instantaneous free surface
   REAL(ReKi)               :: g     ! gravity constant
   REAL(ReKi)               :: k_hat(3), k_hat1(3), k_hat2(3) ! Elemental unit vector pointing from 1st node to 2nd node of the element
   REAL(ReKi)               :: n_hat(3) ! Free surface unit normal vector pointing from water to air
   REAL(ReKi)               :: Fr    !radial component of buoyant force
   REAL(ReKi)               :: Fl    !axial component of buoyant force
   REAL(ReKi)               :: Moment     !moment induced about the center of the cylinder's bottom face
   INTEGER(IntKi)           :: im    ! counter   
   REAL(ReKi)               :: a_s1(3)       
   REAL(ReKi)               :: alpha_s1(3)
   REAL(ReKi)               :: omega_s1(3)
   REAL(ReKi)               :: a_s2(3)       
   REAL(ReKi)               :: alpha_s2(3)
   REAL(ReKi)               :: omega_s2(3)
   REAL(ReKi)               :: pos1(3), pos2(3)
   REAL(ReKi)               :: Imat(3,3)
   REAL(ReKi)               :: iArm(3), iTerm(3), h_c, dRdl_p, dRdl_pp, dSadl_p, dSadl_pp, dSbdl_p, dSbdl_pp, f_hydro(3), Am(3,3), lstar, deltal, deltalLeft, deltalRight
   REAL(ReKi)               :: h, h_c_AM, deltal_AM
   REAL(ReKi)               :: F_WMG(6), F_IMG(6), F_If(6), F_B0(6), F_B1(6), F_B2(6), F_B_End(6)
   REAL(ReKi)               :: AM_End(3,3), An_End(3), DP_Const_End(3), I_MG_End(3,3)

   ! Local variables needed for wave stretching and load smoothing/redistribution
   INTEGER(IntKi)           :: FSElem
   REAL(ReKi)               :: SubRatio
   REAL(ReKi)               :: Zeta1
   REAL(ReKi)               :: Zeta2
   REAL(ReKi)               :: FSInt(3)
   REAL(ReKi)               :: F_D0(3)
   REAL(ReKi)               :: F_A0(3)
   REAL(ReKi)               :: F_I0(3)
   REAL(ReKi)               :: F_0(3)
   REAL(ReKi)               :: F_DS(3)
   REAL(ReKi)               :: F_AS(3)
   REAL(ReKi)               :: F_IS(3)
   REAL(ReKi)               :: F_S(3)
   REAL(ReKi)               :: f_redist
   REAL(ReKi)               :: Df_hydro(3)
   REAL(ReKi)               :: DM_hydro(3)
   REAL(ReKi)               :: Df_hydro_lumped(6)
   REAL(ReKi)               :: FVFSInt(3)
   REAL(ReKi)               :: FAFSInt(3)
   REAL(ReKi)               :: SAFSInt(3)
   REAL(ReKi)               :: FDynPFSInt
   REAL(ReKi)               :: vrelFSInt(3)
   REAL(ReKi)               :: FAMCFFSInt(3)
   INTEGER(IntKi)           :: MemSubStat, NumFSX
   REAL(DbKi)               :: theta1, theta2
   REAL(ReKi)               :: x_hat(3), x_hat1(3), x_hat2(3), y_hat(3), y_hat1(3), y_hat2(3), z_hat(3), posMid(3), zetaMid, FSPt(3)
   INTEGER(IntKi)           :: secStat
   INTEGER(IntKi)           :: nodeInWater
   REAL(SiKi)               :: WaveElev1, WaveElev2, WaveElev, FDynP, FV(3), FA(3), FAMCF(3)
   LOGICAL                  :: Is1stElement, Is1stFloodedMember

   ! Initialize errStat
   errStat = ErrID_None
   errMsg  = ""
   Imat    = 0.0_ReKi
   g       = p%Gravity
   
   !===============================================================================================
   ! Get displaced positions of the hydrodynamic nodes   
   CALL GetDisplacedNodePosition( .FALSE., m%DispNodePosHdn ) ! For hydrodynamic loads; depends on WaveDisp and WaveStMod
   CALL GetDisplacedNodePosition( .TRUE. , m%DispNodePosHst ) ! For hydrostatic loads;  always use actual displaced position

   !===============================================================================================
   ! Calculate the fluid kinematics at all mesh nodes and store for use in the equations below
   CALL WaveField_GetWaveKin( p%WaveField, m%WaveField_m, Time, m%DispNodePosHdn, .FALSE., .TRUE., m%nodeInWater, m%WaveElev1, m%WaveElev2, m%WaveElev, m%FDynP, m%FV, m%FA, m%FAMCF, ErrStat2, ErrMsg2 )
     CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   ! Compute fluid velocity relative to the structure
   DO j = 1, p%NNodes
      m%vrel(:,j)  = ( m%FV(:,j) - u%Mesh%TranslationVel(:,j) ) * m%nodeInWater(j)
   END DO

   !===============================================================================================
   ! Get the instantaneous highest point of internal ballast for each filled group
   ! This is the elevation with zero internal hydrostatic pressure
   DO i = 1,p%NFillGroups
      IF ( p%FilledGroups(i)%IsOpen ) THEN
         m%zFillGroup(i) = 0.0  ! SWL because ballast group open to the environment follows the external hydrostatic pressure field
      ELSE
         Is1stFloodedMember = .true.
         DO j = 1,p%FilledGroups(i)%FillNumM
            im = p%FilledGroups(i)%FillMList(j)
            IF (p%Members(im)%memfloodstatus>0) THEN
               CALL getMemBallastHiPt(p%Members(im),z_hi,ErrStat2,ErrMsg2); if (Failed()) return
               IF ( Is1stFloodedMember ) THEN
                  m%zFillGroup(i) = z_hi
                  Is1stFloodedMember = .false.
               ELSE
                  m%zFillGroup(i) = MAX(m%zFillGroup(i), z_hi)
               END IF
            END IF
         END DO
      END IF
   END DO

   ! ==============================================================================================
   ! Calculate instantaneous loads on each member except for the hydrodynamic loads on member ends.
   ! This covers aspects of the load calculations previously in CreateDistributedMesh.  

   ! Zero out previous time-steps loads (these are loads which are computed at the member-level and summed onto a node, 
   !    so they need to be zeroed out before the summations happen)
   m%F_BF_End    = 0.0_ReKi
   m%F_B_End     = 0.0_ReKi
   y%Mesh%Force  = 0.0_ReKi
   y%Mesh%Moment = 0.0_ReKi
   
   ! Loop through each member
   DO im = 1, p%NMembers    
      mem = p%Members(im)
      N   = mem%NElements
      call YawMember(mem, u%PtfmRefY, ErrStat2, ErrMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      !zero member loads
      m%memberLoads(im)%F_B   = 0.0_ReKi
      m%memberLoads(im)%F_BF  = 0.0_ReKi
      m%memberLoads(im)%F_D   = 0.0_ReKi
      m%memberLoads(im)%F_A   = 0.0_ReKi
      m%memberLoads(im)%F_I   = 0.0_ReKi
      m%memberLoads(im)%F_WMG = 0.0_ReKi
      m%memberLoads(im)%F_IMG = 0.0_ReKi
      m%memberLoads(im)%F_If  = 0.0_ReKi

      ! Determine member submergence status
      IF ( p%WaveField%WaveStMod .EQ. 0_IntKi ) THEN ! No wave stretching - Only need to check the two ends
         IF ( m%nodeInWater(mem%NodeIndx(1)) .NE. m%nodeInWater(mem%NodeIndx(N+1)) ) THEN
            MemSubStat = 1_IntKi  ! Member centerline crosses the SWL once
         ELSE IF ( m%nodeInWater(mem%NodeIndx(1)) .EQ. 0_IntKi ) THEN
            MemSubStat = 3_IntKi  ! Member centerline completely above water
         ELSE
            MemSubStat = 0_IntKi  ! Member centerline fully submerged
         END IF 
      ELSE IF ( p%WaveField%WaveStMod > 0_IntKi ) THEN ! Has wave stretching - Need to check every node
         NumFSX = 0_IntKi ! Number of free-surface crossing
         DO i = 1, N ! loop through member elements
            IF ( m%nodeInWater(mem%NodeIndx(i)) .NE. m%nodeInWater(mem%NodeIndx(i+1)) ) THEN
               NumFSX = NumFSX + 1
            END IF
         END DO
         IF (NumFSX .EQ. 1_IntKi) THEN
            MemSubStat = 1_IntKi  ! Member centerline crosses the free surface once
         ELSE IF (NumFSX .GT. 1_IntKi) THEN
            MemSubStat = 2_IntKi  ! Member centerline crosses the free surface multiple time
         ELSE ! Member centerline does not cross the free surface
            IF ( m%nodeInWater(mem%NodeIndx(1)) .EQ. 0_IntKi ) THEN
               MemSubStat = 3_IntKi  ! Member centerline completely above water
            ELSE
               MemSubStat = 0_IntKi  ! Member centerline completely submerged
            END IF
         END IF
      END IF

      !---------------- Marine growth and Buoyancy: Sides: Only if member not modeled with potential flow theory ----------------  
      IF ( .NOT. mem%PropPot ) THEN ! Member is NOT modeled with Potential Flow Theory
         DO i = max(mem%i_floor,1), N    ! loop through member elements that are not completely buried in the seabed
         
            ! calculate instantaneous incline angle and heading, and related trig values
            ! the first and last NodeIndx values point to the corresponding Joint nodes indices which are at the start of the Mesh
            pos1    = m%DispNodePosHst(:, mem%NodeIndx(i  ))
            pos2    = m%DispNodePosHst(:, mem%NodeIndx(i+1))

            call GetOrientationAngles( pos1, pos2, phi, sinPhi, cosPhi, tanPhi, sinBeta, cosBeta, k_hat, errStat2, errMsg2 )
              call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            ! Compute element to global DirCos matrix for undisplaced structure first
            call Morison_DirCosMtrx( u%Mesh%Position(:,mem%NodeIndx(i  )), u%Mesh%Position(:,mem%NodeIndx(i+1)), mem%MSpinOrient, CMatrix )
            ! Prepend body motion - Assuming the rotation of the starting node is representative of the whole element
            CMatrix = matmul(transpose(u%Mesh%Orientation(:,:,mem%NodeIndx(i))),CMatrix)
            CTrans  = transpose(CMatrix)
            ! Note: CMatrix is element local to global displaced. CTrans is the opposite.
            ! save some commonly used variables   
            dl        = mem%dl
            z1        = pos1(3)          ! get node z locations from input mesh
            z2        = pos2(3)
            a_s1      = u%Mesh%TranslationAcc(:, mem%NodeIndx(i  ))
            alpha_s1  = u%Mesh%RotationAcc   (:, mem%NodeIndx(i  ))
            omega_s1  = u%Mesh%RotationVel   (:, mem%NodeIndx(i  ))
            a_s2      = u%Mesh%TranslationAcc(:, mem%NodeIndx(i+1))
            alpha_s2  = u%Mesh%RotationAcc   (:, mem%NodeIndx(i+1))
            omega_s2  = u%Mesh%RotationVel   (:, mem%NodeIndx(i+1))
            IF (mem%MSecGeom == MSecGeom_Cyl) THEN
               r1         = mem%RMG(i  )      ! outer radius at element nodes including marine growth
               r2         = mem%RMG(i+1)
               r1b        = mem%RMGB(i  )     ! outer radius at element nodes including marine growth scaled by sqrt(Cb)
               r2b        = mem%RMGB(i+1)
               dRdl_mg    = mem%dRdl_mg(i)    ! Taper of element including marine growth
               dRdl_mg_b  = mem%dRdl_mg_b(i)  ! Taper of element including marine growth with radius scaling by sqrt(Cb)
            ELSE IF (mem%MSecGeom == MSecGeom_Rec) THEN
               Sa1        = mem%SaMG(i  )     ! outer side A at element nodes including marine growth
               Sa2        = mem%SaMG(i+1)
               Sb1        = mem%SbMG(i  )     ! outer side B at element nodes including marine growth
               Sb2        = mem%SbMG(i+1)
               Sa1b       = mem%SaMGB(i  )    ! outer side A at element nodes including marine growth scaled by sqrt(Cb)
               Sa2b       = mem%SaMGB(i+1)
               Sb1b       = mem%SbMGB(i  )    ! outer side B at element nodes including marine growth scaled by sqrt(Cb)
               Sb2b       = mem%SbMGB(i+1)
               dSadl_mg   = mem%dSadl_mg(i)   ! Taper of element side A including marine growth
               dSadl_mg_b = mem%dSadl_mg_b(i) ! Taper of element side A including marine growth with radius scaling by sqrt(Cb)
               dSbdl_mg   = mem%dSbdl_mg(i)   ! Taper of element side B including marine growth
               dSbdl_mg_b = mem%dSbdl_mg_b(i) ! Taper of element side B including marine growth with radius scaling by sqrt(Cb)
            END IF

            ! ------------------ marine growth: Sides: Section 4.1.2 --------------------  
            ! ----- marine growth weight
            F_WMG = 0.0_ReKi

            ! lower node
            F_WMG(3) = - mem%m_mg_l(i)*g ! weight force  : Note: this is a constant
            F_WMG(4) = - mem%m_mg_l(i)*g * mem%h_cmg_l(i)* sinPhi * sinBeta! weight force
            F_WMG(5) =   mem%m_mg_l(i)*g * mem%h_cmg_l(i)* sinPhi * cosBeta! weight force
            m%memberLoads(im)%F_WMG(:,i) = m%memberLoads(im)%F_WMG(:,i) + F_WMG
            y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + F_WMG(1:3)
            y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + F_WMG(4:6)
            
            ! upper node
            F_WMG(3) = - mem%m_mg_u(i)*g ! weight force  : Note: this is a constant 
            F_WMG(4) = - mem%m_mg_u(i)*g * mem%h_cmg_u(i)* sinPhi * sinBeta! weight force
            F_WMG(5) =   mem%m_mg_u(i)*g * mem%h_cmg_u(i)* sinPhi * cosBeta! weight force
            m%memberLoads(im)%F_WMG(:,i+1) = m%memberLoads(im)%F_WMG(:,i+1) + F_WMG  
            y%Mesh%Force (:,mem%NodeIndx(i+1)) = y%Mesh%Force (:,mem%NodeIndx(i+1)) + F_WMG(1:3)
            y%Mesh%Moment(:,mem%NodeIndx(i+1)) = y%Mesh%Moment(:,mem%NodeIndx(i+1)) + F_WMG(4:6)
            
            ! ----- marine growth inertial load
            ! lower node
            Imat      = 0.0_ReKi
            IF (mem%MSecGeom == MSecGeom_Cyl) THEN
               Imat(1,1) = mem%I_rmg_l(i)
               Imat(2,2) = mem%I_rmg_l(i)
            ELSE IF (mem%MSecGeom == MSecGeom_Rec) THEN
               Imat(1,1) = mem%I_xmg_l(i)
               Imat(2,2) = mem%I_ymg_l(i)
            END IF
            Imat(3,3) = mem%I_lmg_l(i)
            Imat      =  matmul(matmul(CMatrix, Imat), CTrans)
            iArm = mem%h_cmg_l(i) * k_hat
            iTerm     = ( -a_s1 - cross_product(omega_s1, cross_product(omega_s1,iArm )) - cross_product(alpha_s1,iArm) ) * mem%m_mg_l(i)
            F_IMG(1:3) = iTerm
            F_IMG(4:6) = - matmul(Imat, alpha_s1) - cross_product(iArm,a_s1 * mem%m_mg_l(i)) &
                         - cross_product(omega_s1,matmul(Imat,omega_s1))
            m%memberLoads(im)%F_IMG(:,i) = m%memberLoads(im)%F_IMG(:,i) + F_IMG
            y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + F_IMG(1:3)
            y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + F_IMG(4:6)

            ! upper node
            Imat      = 0.0_ReKi
            IF (mem%MSecGeom == MSecGeom_Cyl) THEN
               Imat(1,1) = mem%I_rmg_u(i)
               Imat(2,2) = mem%I_rmg_u(i)
            ELSE IF (mem%MSecGeom == MSecGeom_Rec) THEN
               Imat(1,1) = mem%I_xmg_u(i)
               Imat(2,2) = mem%I_ymg_u(i)
            END IF
            Imat(3,3) = mem%I_lmg_u(i)
            Imat      =  matmul(matmul(CMatrix, Imat), CTrans)
            iArm = mem%h_cmg_u(i) * k_hat
            iTerm     = ( -a_s2 - cross_product(omega_s2, cross_product(omega_s2,iArm )) - cross_product(alpha_s2,iArm) ) * mem%m_mg_u(i)
            F_IMG(1:3) = iTerm
            F_IMG(4:6) = - matmul(Imat, alpha_s2) - cross_product(iArm,a_s2 * mem%m_mg_u(i)) &
                         - cross_product(omega_s2,matmul(Imat,omega_s2))
            m%memberLoads(im)%F_IMG(:,i+1) = m%memberLoads(im)%F_IMG(:,i+1) + F_IMG
            y%Mesh%Force (:,mem%NodeIndx(i+1)) = y%Mesh%Force (:,mem%NodeIndx(i+1)) + F_IMG(1:3)
            y%Mesh%Moment(:,mem%NodeIndx(i+1)) = y%Mesh%Moment(:,mem%NodeIndx(i+1)) + F_IMG(4:6)

            ! ------------------- buoyancy loads: sides: Sections 3.1 and 3.2 ------------------------
            IF (mem%MHstLMod == 1) THEN
               IF ( p%WaveField%WaveStMod > 0_IntKi ) THEN ! If wave stretching is enabled, compute buoyancy up to free surface
                  CALL GetTotalWaveElev( Time, pos1, Zeta1, ErrStat2, ErrMsg2 )
                  CALL GetTotalWaveElev( Time, pos2, Zeta2, ErrStat2, ErrMsg2 )
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               ELSE ! Without wave stretching, compute buoyancy based on SWL
                  Zeta1 = 0.0_ReKi
                  Zeta2 = 0.0_ReKi
               END IF
               Is1stElement = ( i .EQ. 1)
               CALL getElementHstLds_Mod1(mem, Time, pos1, pos2, Zeta1, Zeta2, k_hat, r1b, r2b, dl, mem%alpha(i), Is1stElement, F_B0, F_B1, F_B2, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               ! Add nodal loads to mesh
               IF ( .NOT. Is1stElement ) THEN
                  m%memberLoads(im)%F_B(:, i-1) = m%memberLoads(im)%F_B(:, i-1) + F_B0
                  y%Mesh%Force (:,mem%NodeIndx(i-1)) = y%Mesh%Force (:,mem%NodeIndx(i-1)) + F_B0(1:3)
                  y%Mesh%Moment(:,mem%NodeIndx(i-1)) = y%Mesh%Moment(:,mem%NodeIndx(i-1)) + F_B0(4:6)
               END IF
               m%memberLoads(im)%F_B(:, i  ) = m%memberLoads(im)%F_B(:, i  ) + F_B1
               m%memberLoads(im)%F_B(:, i+1) = m%memberLoads(im)%F_B(:, i+1) + F_B2
               y%Mesh%Force (:,mem%NodeIndx(i  )) = y%Mesh%Force (:,mem%NodeIndx(i  )) + F_B1(1:3)
               y%Mesh%Moment(:,mem%NodeIndx(i  )) = y%Mesh%Moment(:,mem%NodeIndx(i  )) + F_B1(4:6)
               y%Mesh%Force (:,mem%NodeIndx(i+1)) = y%Mesh%Force (:,mem%NodeIndx(i+1)) + F_B2(1:3)
               y%Mesh%Moment(:,mem%NodeIndx(i+1)) = y%Mesh%Moment(:,mem%NodeIndx(i+1)) + F_B2(4:6)
            ELSE IF (mem%MHstLMod == 2) THEN ! Alternative hydrostatic load calculation
               ! Get free surface elevation and normal at the element midpoint (both assumed constant over the element)
               posMid = 0.5 * (pos1+pos2)
               ! rn is only used to estimate free surface normal numerically
               IF (mem%MSecGeom == MSecGeom_Cyl) THEN
                  rn  = 0.5 * (r1b +r2b )
               ELSE IF (mem%MSecGeom == MSecGeom_Rec) THEN
                  rn  = MAX( 0.5*(Sa1b+Sa2b), 0.5*(Sb1b+Sb2b) )
               END IF
               IF (p%WaveField%WaveStMod > 0) THEN
                  CALL GetTotalWaveElev( Time, posMid, ZetaMid, ErrStat2, ErrMsg2 )
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  CALL GetFreeSurfaceNormal( Time, posMid, rn, n_hat, ErrStat2, ErrMsg2 )
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  FSPt = (/posMid(1),posMid(2),ZetaMid/) ! Reference point on the free surface
               ELSE
                  FSPt = (/posMid(1),posMid(2),0.0_ReKi/)
                  n_hat = (/0.0,0.0,1.0/)
               END IF

               IF (mem%MSecGeom == MSecGeom_Cyl) THEN
                  CALL GetSectionUnitVectors_Cyl( k_hat, y_hat, z_hat )
                  CALL getElementHstLds_Mod2_Cyl( pos1, pos2, FSPt, k_hat, y_hat, z_hat, n_hat, r1b, r2b, dl, F_B1, F_B2, ErrStat2, ErrMsg2)
               ELSE IF (mem%MSecGeom == MSecGeom_Rec) THEN
                  CALL GetSectionUnitVectors_Rec( CMatrix, x_hat, y_hat )
                  CALL getElementHstLds_Mod2_Rec( pos1, pos2, FSPt, k_hat, x_hat, y_hat, n_hat, Sa1b, Sa2b, Sb1b, Sb2b, dl, F_B1, F_B2, ErrStat2, ErrMsg2)
               END IF
               CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

               ! Add nodal loads to mesh
               m%memberLoads(im)%F_B(:,i  ) = m%memberLoads(im)%F_B(:,i  ) + F_B1
               m%memberLoads(im)%F_B(:,i+1) = m%memberLoads(im)%F_B(:,i+1) + F_B2
               y%Mesh%Force (:,mem%NodeIndx(i  )) = y%Mesh%Force (:,mem%NodeIndx(i  )) + F_B1(1:3)
               y%Mesh%Moment(:,mem%NodeIndx(i  )) = y%Mesh%Moment(:,mem%NodeIndx(i  )) + F_B1(4:6)
               y%Mesh%Force (:,mem%NodeIndx(i+1)) = y%Mesh%Force (:,mem%NodeIndx(i+1)) + F_B2(1:3)
               y%Mesh%Moment(:,mem%NodeIndx(i+1)) = y%Mesh%Moment(:,mem%NodeIndx(i+1)) + F_B2(4:6)
            END IF ! MHstLMod
         END DO ! i = max(mem%i_floor,1), N    ! loop through member elements that are not fully buried in the seabed
      END IF ! NOT Modeled with Potential flow theory

      ! --------------------------- flooded ballast: sides: Always compute regardless of PropPot setting ------------------------------
      ! NOTE: For memfloodstatus and floodstatus: 0 = fully buried or not ballasted, 1 = fully flooded, 2 = partially flooded
      IF ( mem%memfloodstatus > 0 ) THEN  ! Fully or partially flooded member
         zFillGroup = m%zFillGroup(mem%MmbrFilledIDIndx)
         DO i = max(mem%i_floor,1), N    ! loop through member elements that are not completely buried in the seabed
            IF (mem%floodstatus(i)>0) THEN
               ! calculate instantaneous incline angle and heading, and related trig values
               ! the first and last NodeIndx values point to the corresponding Joint nodes indices which are at the start of the Mesh
               pos1 = m%DispNodePosHst(:,mem%NodeIndx(i  ))
               pos2 = m%DispNodePosHst(:,mem%NodeIndx(i+1))

               call GetOrientationAngles( pos1, pos2, phi, sinPhi, cosPhi, tanPhi, sinBeta, cosBeta, k_hat, errStat2, errMsg2 ); if (Failed()) return
               ! Compute element to global DirCos matrix for undisplaced structure first
               call Morison_DirCosMtrx( u%Mesh%Position(:,mem%NodeIndx(i  )), u%Mesh%Position(:,mem%NodeIndx(i+1)), mem%MSpinOrient, CMatrix )
               ! Prepend body motion - Assuming the rotation of the starting node is representative of the whole element
               CMatrix = matmul(transpose(u%Mesh%Orientation(:,:,mem%NodeIndx(i))),CMatrix)
               CTrans  = transpose(CMatrix)
               ! Note: CMatrix is element local to global displaced. CTrans is the opposite.
               ! save some commonly used variables
               dl        = mem%dl
               z1        = pos1(3)          ! get displaced node z locations
               z2        = pos2(3)
               a_s1      = u%Mesh%TranslationAcc(:, mem%NodeIndx(i  ))
               alpha_s1  = u%Mesh%RotationAcc   (:, mem%NodeIndx(i  ))
               omega_s1  = u%Mesh%RotationVel   (:, mem%NodeIndx(i  ))
               a_s2      = u%Mesh%TranslationAcc(:, mem%NodeIndx(i+1))
               alpha_s2  = u%Mesh%RotationAcc   (:, mem%NodeIndx(i+1))
               omega_s2  = u%Mesh%RotationVel   (:, mem%NodeIndx(i+1))
               IF (mem%MSecGeom == MSecGeom_Cyl) THEN
                  r1In = mem%Rin(i  )      ! outer radius at element nodes including marine growth
                  r2In = mem%Rin(i+1)
                  IF ( mem%floodstatus(i) == 1 ) THEN    ! Fully flooded element
                     zMid   = 0.5 * (z1   + z2  )
                     rMidIn = 0.5 * (r1In + r2In)
                  ELSE                                   ! Partially flooded element
                     zMid   = z1 + mem%h_fill * k_hat(3)
                     l      = mem%h_fill/mem%dl
                     rMidIn = r1In * (1.0-l) + r2In * l
                  END IF
               ELSE IF (mem%MSecGeom == MSecGeom_Rec) THEN
                  Sa1In      = mem%Sain(i  )     ! outer side A at element nodes including marine growth
                  Sa2In      = mem%Sain(i+1)
                  Sb1In      = mem%Sbin(i  )     ! outer side B at element nodes including marine growth
                  Sb2In      = mem%Sbin(i+1)
                  IF ( mem%floodstatus(i) == 1 ) THEN    ! Fully flooded element
                     zMid    = 0.5 * (z1   + z2  )
                     SaMidIn = 0.5 * (Sa1In + Sa2In)
                     SbMidIn = 0.5 * (Sb1In + Sb2In)
                  ELSE                                   ! Partially flooded element
                     zMid    = z1 + mem%h_fill * k_hat(3)
                     l       = mem%h_fill/mem%dl
                     SaMidIn = Sa1In * (1.0-l) + Sa2In * l
                     SbMidIn = Sb1In * (1.0-l) + Sb2In * l
                  END IF
               END IF

               ! ------------------ flooded ballast inertia: sides: Section 6.1.1 : Always compute regardless of PropPot setting ---------------------
               ! lower node
               Imat      = 0.0_ReKi
               IF (mem%MSecGeom == MSecGeom_Cyl) THEN
                  Imat(1,1) = mem%I_rfb_l(i)
                  Imat(2,2) = mem%I_rfb_l(i)
               ELSE IF (mem%MSecGeom == MSecGeom_Rec) THEN
                  Imat(1,1) = mem%I_xfb_l(i)
                  Imat(2,2) = mem%I_yfb_l(i)
               END IF
               Imat(3,3) = mem%I_lfb_l(i)
               Imat      =  matmul(matmul(CMatrix, Imat), CTrans)
               iArm = mem%h_cfb_l(i) * k_hat
               iTerm     = ( -a_s1  - cross_product(omega_s1, cross_product(omega_s1,iArm ))  -  cross_product(alpha_s1,iArm) ) * mem%m_fb_l(i)
               F_If(1:3) =  iTerm
               F_If(4:6) =  - matmul(Imat, alpha_s1) - cross_product(iArm,a_s1 * mem%m_fb_l(i)) &
                            - cross_product(omega_s1,matmul(Imat,omega_s1))
               m%memberLoads(im)%F_If(:,i) = m%memberLoads(im)%F_If(:,i) + F_If
               y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + F_If(1:3)
               y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + F_If(4:6)

               ! upper node
               Imat      = 0.0_ReKi
               IF (mem%MSecGeom == MSecGeom_Cyl) THEN
                  Imat(1,1) = mem%I_rfb_u(i)
                  Imat(2,2) = mem%I_rfb_u(i)
               ELSE IF (mem%MSecGeom == MSecGeom_Rec) THEN
                  Imat(1,1) = mem%I_xfb_u(i)
                  Imat(2,2) = mem%I_yfb_u(i)
               END IF
               Imat(3,3) = mem%I_lfb_u(i)
               Imat      =  matmul(matmul(CMatrix, Imat), CTrans)
               iArm = mem%h_cfb_u(i) * k_hat
               iTerm     = ( -a_s2  - cross_product(omega_s2, cross_product(omega_s2,iArm ))  -  cross_product(alpha_s2,iArm) ) * mem%m_fb_u(i)
               F_If(1:3) = iTerm
               F_If(4:6) = - matmul(Imat, alpha_s2) - cross_product(iArm,a_s2 * mem%m_fb_u(i)) &
                           - cross_product(omega_s2,matmul(Imat,omega_s2))
               m%memberLoads(im)%F_If(:,i+1) = m%memberLoads(im)%F_If(:,i+1) + F_If
               y%Mesh%Force (:,mem%NodeIndx(i+1)) = y%Mesh%Force (:,mem%NodeIndx(i+1)) + F_If(1:3)
               y%Mesh%Moment(:,mem%NodeIndx(i+1)) = y%Mesh%Moment(:,mem%NodeIndx(i+1)) + F_If(4:6)

               ! ------------------ flooded ballast weight : sides : Section 5.1.2 & 5.2.2  : Always compute regardless of PropPot setting ---------------------
               F_B1 = 0.0
               F_B2 = 0.0
               IF (mem%MSecGeom == MSecGeom_Cyl) THEN
                  F_B1(3)   = - p%gravity * mem%m_fb_l(i)
                  F_B1(1:3) = F_B1(1:3) + mem%FillDens * p%gravity * pi * ( rMidIn*rMidIn*(zMid-zFillGroup) - r1In*r1In*(z1-zFillGroup) ) * k_hat
                  F_B1(4:6) = -( p%gravity * mem%m_fb_l(i) * mem%h_cfb_l(i) + mem%FillDens * p%gravity * 0.25*pi*(rMidIn**4-r1In**4) ) * Cross_Product(k_hat,(/0.0,0.0,1.0/))
                  IF ( mem%FloodStatus(i) == 1 ) THEN
                     F_B2(3)   = - p%gravity * mem%m_fb_u(i)
                     F_B2(1:3) = F_B2(1:3) + mem%FillDens * p%gravity * pi * ( r2In*r2In*(z2-zFillGroup) - rMidIn*rMidIn*(zMid-zFillGroup) ) * k_hat
                     F_B2(4:6) = -( p%gravity * mem%m_fb_u(i) * mem%h_cfb_u(i) + mem%FillDens * p%gravity * 0.25*pi*(r2In**4-rMidIn**4) ) * Cross_Product(k_hat,(/0.0,0.0,1.0/))
                  ELSE IF ( i == mem%elem_fill ) THEN ! Need to include end load here
                     F_B1(1:3) = F_B1(1:3) + mem%FillDens * p%gravity *        pi * rMidIn**2* (zFillGroup - zMid) * k_hat
                     F_B1(4:6) = F_B1(4:6) + mem%FillDens * p%gravity * 0.25 * pi * rMidIn**4* Cross_Product(k_hat,(/0.0,0.0,1.0/))
                  END IF
               ELSE IF (mem%MSecGeom == MSecGeom_Rec) THEN
                  CALL GetSectionUnitVectors_Rec( CMatrix, x_hat, y_hat )
                  F_B1(3)   = - p%gravity * mem%m_fb_l(i)
                  F_B1(1:3) = F_B1(1:3) + mem%FillDens * p%gravity * ( SaMidIn*SbMidIn*(zMid-zFillGroup) - Sa1In*Sb1In*(z1-zFillGroup) ) * k_hat
                  F_B1(4:6) = - p%gravity * mem%m_fb_l(i) * mem%h_cfb_l(i) * Cross_Product(k_hat,(/0.0,0.0,1.0/)) &
                              + mem%FillDens * p%gravity / 12.0 * ( (Sa1In**3*Sb1In*x_hat(3)*y_hat - Sa1In*Sb1In**3*y_hat(3)*x_hat ) - &
                                                                    (SaMidIn**3*SbMidIn*x_hat(3)*y_hat - SaMidIn*SbMidIn**3*y_hat(3)*x_hat ) )
                  IF ( mem%FloodStatus(i) == 1 ) THEN
                     F_B2(3)   = - p%gravity * mem%m_fb_u(i)
                     F_B2(1:3) = F_B2(1:3) + mem%FillDens * p%gravity * ( Sa2In*Sb2In*(z2-zFillGroup) - SaMidIn*SbMidIn*(zMid-zFillGroup) ) * k_hat
                     F_B2(4:6) = - p%gravity * mem%m_fb_u(i) * mem%h_cfb_u(i) * Cross_Product(k_hat,(/0.0,0.0,1.0/)) &
                                 + mem%FillDens * p%gravity / 12.0 * ( (SaMidIn**3*SbMidIn*x_hat(3)*y_hat - SaMidIn*SbMidIn**3*y_hat(3)*x_hat ) - &
                                                                       (Sa2In**3*Sb2In*x_hat(3)*y_hat - Sa2In*Sb2In**3*y_hat(3)*x_hat ) )
                  ELSE IF ( i == mem%elem_fill ) THEN ! Need to include end load here
                     F_B1(1:3) = F_B1(1:3) + mem%FillDens * p%gravity * SaMidIn*SbMidIn*(zFillGroup-zMid) * k_hat
                     F_B1(4:6) = F_B1(4:6) + mem%FillDens * p%gravity / 12.0 * (SaMidIn**3*SbMidIn*x_hat(3)*y_hat - SaMidIn*SbMidIn**3*y_hat(3)*x_hat)
                  END IF
               END IF

               m%memberLoads(im)%F_BF(:, i  ) = m%memberLoads(im)%F_BF(:, i  ) + F_B1
               m%memberLoads(im)%F_BF(:, i+1) = m%memberLoads(im)%F_BF(:, i+1) + F_B2
               y%Mesh%Force (:,mem%NodeIndx(i  )) = y%Mesh%Force (:,mem%NodeIndx(i  )) + F_B1(1:3)
               y%Mesh%Moment(:,mem%NodeIndx(i  )) = y%Mesh%Moment(:,mem%NodeIndx(i  )) + F_B1(4:6)
               y%Mesh%Force (:,mem%NodeIndx(i+1)) = y%Mesh%Force (:,mem%NodeIndx(i+1)) + F_B2(1:3)
               y%Mesh%Moment(:,mem%NodeIndx(i+1)) = y%Mesh%Moment(:,mem%NodeIndx(i+1)) + F_B2(4:6)

            END IF                                     ! mem%floodstatus(i) > 0
         END DO ! i = max(mem%i_floor,1), N        ! loop through member elements that are not fully buried in the seabed
      END IF                                    ! Fully or partially flooded member

      !-----------------------------------------------------------------------------------------------------!
      !                               External Hydrodynamic Side Loads - Start                              !
      !-----------------------------------------------------------------------------------------------------!
      IF ( p%WaveField%WaveStMod > 0 .AND. MemSubStat == 1 .AND. (m%NodeInWater(mem%NodeIndx(N+1)).EQ.0_IntKi) ) THEN 
      !----------------------------Apply load smoothing----------------------------!
      ! only when:
      ! 1. wave stretching is enabled
      ! 2. member centerline crosses the free surface exactly once
      ! 3. the last node is out of water, which implies the first node is in water
                
        FSElem = -1 ! Initialize the No. of the partially wetted element as -1
      
        DO i = mem%i_floor+1,N ! loop through member nodes starting from the first node above seabed, but skip the last node which should not be submerged anyways
           
           ! Get positions of node i and i+1
           pos1 = m%DispNodePosHdn(:,mem%NodeIndx(i  ))
           pos2 = m%DispNodePosHdn(:,mem%NodeIndx(i+1))

           ! Free surface elevation above or below node i and i+1
           Zeta1 = m%WaveElev(mem%NodeIndx(i))
           Zeta2 = m%WaveElev(mem%NodeIndx(i+1))

           ! Compute deltal and h_c
           IF ( i == 1 ) THEN ! First node
              deltal = mem%dl/2.0_ReKi
              h_c    = mem%dl/4.0_ReKi
           ELSE IF ( i == mem%i_floor + 1 ) THEN ! This node is the upper node of an element which crosses the seabed
              ! Superceded by i==1 above if mem%i_floor = 0
              deltal = mem%dl/2.0_ReKi - mem%h_floor
              h_c    = 0.5_ReKi*(mem%dl/2.0_ReKi + mem%h_floor)
           ELSE
              ! This node is an interior node. Note: Element crossing the free surface will be handled at the end in conjunction with wave stretching
              deltal = mem%dl
              h_c    = 0.0_ReKi           
           END IF ! Note: No need to consider i==N+1 because we do not allow the top node to become submerged. The loop also does not reach N+1.
           
           IF ( pos1(3) <= Zeta1 .AND. pos2(3) > Zeta2 ) THEN ! element is partially wetted
             ! Record the number of the partially wetted element
             FSElem = i
             ! Calculate submergence ratio
             SubRatio = ( Zeta1-pos1(3) ) / ( (Zeta1-pos1(3)) - (Zeta2-pos2(3)) )
             ! Calculate the position of the intersection between the free surface and the element
             FSInt = SubRatio * (pos2-pos1) + pos1
           END IF
         
           ! Compute the slope of member radius/side length
           IF (mem%MSecGeom==MSecGeom_Cyl) THEN
              IF (i == 1) THEN
                 dRdl_p  = abs(mem%dRdl_mg(i))
                 dRdl_pp = mem%dRdl_mg(i)
              ELSE IF ( i > 1 .AND. i < (N+1)) THEN
                 dRdl_p  = 0.5*( abs(mem%dRdl_mg(i-1)) + abs(mem%dRdl_mg(i)) )
                 dRdl_pp = 0.5*( mem%dRdl_mg(i-1) + mem%dRdl_mg(i) )
              ELSE
                 dRdl_p  = abs(mem%dRdl_mg(N))
                 dRdl_pp = mem%dRdl_mg(N)
              END IF
           ELSE IF (mem%MSecGeom==MSecGeom_Rec) THEN
              IF (i == 1) THEN
                 dSadl_p  = abs(mem%dSadl_mg(i))
                 dSadl_pp = mem%dSadl_mg(i)
                 dSbdl_p  = abs(mem%dSbdl_mg(i))
                 dSbdl_pp = mem%dSbdl_mg(i)
              ELSE IF ( i > 1 .AND. i < (N+1)) THEN
                 dSadl_p  = 0.5*( abs(mem%dSadl_mg(i-1)) + abs(mem%dSadl_mg(i)) )
                 dSadl_pp = 0.5*( mem%dSadl_mg(i-1) + mem%dSadl_mg(i) )
                 dSbdl_p  = 0.5*( abs(mem%dSbdl_mg(i-1)) + abs(mem%dSbdl_mg(i)) )
                 dSbdl_pp = 0.5*( mem%dSbdl_mg(i-1) + mem%dSbdl_mg(i) )
              ELSE
                 dSadl_p  = abs(mem%dSadl_mg(N))
                 dSadl_pp = mem%dSadl_mg(N)
                 dSbdl_p  = abs(mem%dSbdl_mg(N))
                 dSbdl_pp = mem%dSbdl_mg(N)
              END IF
           END IF

           !-------------------- hydrodynamic drag loads: sides: Section 7.1.2 ------------------------!
           vec = matmul( mem%Ak,m%vrel(:,mem%NodeIndx(i)) )
           IF (mem%MSecGeom==MSecGeom_Cyl) THEN
              f_hydro = mem%Cd(i)*p%WaveField%WtrDens*mem%RMG(i)*TwoNorm(vec)*vec  +  &                                              ! radial part
                        0.5*mem%AxCd(i)*p%WaveField%WtrDens * pi*mem%RMG(i)*dRdl_p * &                                               ! axial part
                        abs(dot_product( mem%k, m%vrel(:,mem%NodeIndx(i)) )) * matmul( mem%kkt, m%vrel(:,mem%NodeIndx(i)) )          ! axial part cont'd
           ELSE IF (mem%MSecGeom==MSecGeom_Rec) THEN
              f_hydro = 0.5*mem%CdB(i)*p%WaveField%WtrDens*mem%SbMG(i)*TwoNorm(vec)*Dot_Product(vec,mem%x_hat)*mem%x_hat  +  &       ! local x-direction
                        0.5*mem%CdA(i)*p%WaveField%WtrDens*mem%SaMG(i)*TwoNorm(vec)*Dot_Product(vec,mem%y_hat)*mem%y_hat  +  &       ! local z-direction
                        0.25*mem%AxCd(i)*p%WaveField%WtrDens * (dSadl_p*mem%SbMG(i) + dSbdl_p*mem%SaMG(i)) * &                       ! axial part
                        abs(dot_product( mem%k, m%vrel(:,mem%NodeIndx(i)) )) * matmul( mem%kkt, m%vrel(:,mem%NodeIndx(i)) )          ! axial part cont'd
           END IF
           CALL LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%memberLoads(im)%F_D(:, i) )
           y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + m%memberLoads(im)%F_D(1:3, i)
           y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + m%memberLoads(im)%F_D(4:6, i)
           IF (i == FSElem) THEN ! Save the distributed load at the first node below the free surface
             F_D0 = f_hydro
           END IF
           
           IF ( .NOT. mem%PropPot ) THEN
              !-------------------- hydrodynamic added mass loads: sides: Section 7.1.3 ------------------------!
              IF (mem%MSecGeom==MSecGeom_Cyl) THEN
                 Am = mem%Ca(i)*p%WaveField%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*mem%Ak + 2.0*mem%AxCa(i)*p%WaveField%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*dRdl_p*mem%kkt
                 f_hydro = -matmul( Am, u%Mesh%TranslationAcc(:,mem%NodeIndx(i)) )
              ELSE IF (mem%MSecGeom==MSecGeom_Rec) THEN
                 f_hydro = -p%WaveField%WtrDens*mem%CaB(i) * 0.25*pi*mem%SbMG(i)*mem%SbMG(i) * Dot_Product(u%Mesh%TranslationAcc(:,mem%NodeIndx(i)),mem%x_hat)*mem%x_hat &
                           -p%WaveField%WtrDens*mem%CaA(i) * 0.25*pi*mem%SaMG(i)*mem%SaMG(i) * Dot_Product(u%Mesh%TranslationAcc(:,mem%NodeIndx(i)),mem%y_hat)*mem%y_hat &
                       -0.5*p%WaveField%WtrDens*mem%AxCa(i) * (dSbdl_p*mem%SaMG(i)+dSadl_p*mem%SbMG(i))*SQRT(mem%SaMG(i)*mem%SbMG(i)) * Dot_Product(u%Mesh%TranslationAcc(:,mem%NodeIndx(i)),mem%k)*mem%k
              END IF
              IF ( p%AMMod .EQ. 0_IntKi ) THEN ! Compute added-mass force up to the SWL
                 z1 = u%Mesh%Position(3, mem%NodeIndx(i)) - p%WaveField%MSL2SWL ! Undisplaced z-position of the current node
                 IF ( z1 > 0.0_ReKi ) THEN ! Node is above SWL undisplaced; zero added-mass force
                    f_hydro = 0.0_ReKi
                    CALL LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%memberLoads(im)%F_A(:, i) )
                 ELSE
                    ! Need to compute deltal_AM and h_c_AM based on the formulation without wave stretching.
                    z2 = u%Mesh%Position(3, mem%NodeIndx(i+1)) - p%WaveField%MSL2SWL ! Undisplaced z-position of the next node
                    IF ( z2 > 0.0_ReKi ) THEN ! Element i crosses the SWL
                       h = -z1 / mem%cosPhi_ref ! Length of Element i between SWL and node i, h>=0
                       deltal_AM = mem%dl/2.0 + h
                       h_c_AM    = 0.5*(h-mem%dl/2.0)
                    ELSE
                       deltal_AM = deltal;
                       h_c_AM    = h_c
                    END IF
                    ! Note: Do not overwrite deltal and h_c here. Still need them for the fluid inertia and drag forces.
                    CALL LumpDistrHydroLoads( f_hydro, mem%k, deltal_AM, h_c_AM, m%memberLoads(im)%F_A(:, i) )
                 END IF
              ELSE ! Compute added-mass force up to the instantaneous free surface
                 f_hydro = f_hydro * m%nodeInWater(mem%NodeIndx(i)) ! Zero the force if node above free surface
                 CALL LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%memberLoads(im)%F_A(:, i) )
                 IF (i == FSElem) THEN ! Save the distributed load at the first node below the free surface
                    F_A0 = f_hydro
                 END IF
              END IF ! AMMod 0 or 1
              y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + m%memberLoads(im)%F_A(1:3, i)
              y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + m%memberLoads(im)%F_A(4:6, i)
              
              !--------------------- hydrodynamic inertia loads: sides: Section 7.1.4 --------------------------!
              IF (mem%MSecGeom==MSecGeom_Cyl) THEN
                 IF (mem%PropMCF) THEN
                    f_hydro=                     p%WaveField%WtrDens*pi*mem%RMG(i)*mem%RMG(i)       * matmul( mem%Ak,  m%FAMCF(:,mem%NodeIndx(i)) ) + &
                                 2.0*mem%AxCa(i)*p%WaveField%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*dRdl_p * matmul( mem%kkt, m%FA(:,mem%NodeIndx(i)) ) + &
                                 2.0*m%FDynP(mem%NodeIndx(i))*mem%AxCp(i)*pi*mem%RMG(i)*dRdl_pp*mem%k
                 ELSE
                    f_hydro=(mem%Ca(i)+mem%Cp(i))*p%WaveField%WtrDens*pi*mem%RMG(i)*mem%RMG(i)        * matmul( mem%Ak,  m%FA(:,mem%NodeIndx(i)) ) + &
                                 2.0*mem%AxCa(i) *p%WaveField%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*dRdl_p * matmul( mem%kkt, m%FA(:,mem%NodeIndx(i)) ) + &
                                 2.0*m%FDynP(mem%NodeIndx(i))*mem%AxCp(i)*pi*mem%RMG(i)*dRdl_pp*mem%k
                 END IF
              ELSE IF (mem%MSecGeom==MSecGeom_Rec) THEN
                 ! Note: MacCamy-Fuchs correction cannot be applied to rectangular members
                 f_hydro= mem%Cp(i)*p%WaveField%WtrDens* mem%SaMG(i)*mem%SbMG(i) * matmul( mem%Ak,  m%FA(:,mem%NodeIndx(i)) ) + &                            ! transver FK component
                          m%FDynP(mem%NodeIndx(i))*mem%AxCp(i)* (mem%SaMG(i)*dSbdl_pp+dSadl_pp*mem%SbMG(i)) *mem%k + &                                       ! axial FK component
                          p%WaveField%WtrDens*mem%CaB(i) * 0.25*pi*mem%SbMG(i)*mem%SbMG(i) * Dot_Product(m%FA(:,mem%NodeIndx(i)),mem%x_hat)*mem%x_hat + &    ! x-component of diffraction part
                          p%WaveField%WtrDens*mem%CaA(i) * 0.25*pi*mem%SaMG(i)*mem%SaMG(i) * Dot_Product(m%FA(:,mem%NodeIndx(i)),mem%y_hat)*mem%y_hat + &    ! y-component of diffraction part
                      0.5*p%WaveField%WtrDens*mem%AxCa(i) * (dSbdl_p*mem%SaMG(i)+dSadl_p*mem%SbMG(i))*SQRT(mem%SaMG(i)*mem%SbMG(i)) * &                      ! axial component of diffraction part
                          Dot_Product(m%FA(:,mem%NodeIndx(i)),mem%k)*mem%k                                                                                   ! axial component of diffraction part cont'd
              END IF
              CALL LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%memberLoads(im)%F_I(:, i) )
              y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + m%memberLoads(im)%F_I(1:3, i)
              y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + m%memberLoads(im)%F_I(4:6, i)
              IF (i == FSElem) THEN ! Save the distributed load at the first node below the free surface
                 F_I0 = f_hydro
              END IF
           END IF
           
        END DO ! i =1,N+1    ! loop through member nodes  

        !----------------------------------------------------------------------------------------------------!
        ! Compute the distributed loads at the point of intersection between the member and the free surface !
        !----------------------------------------------------------------------------------------------------!   
        ! Get wave kinematics at the free-surface intersection. Set forceNodeInWater=.TRUE. to guarantee the free-surface intersection is in water.
        CALL WaveField_GetNodeWaveKin( p%WaveField, m%WaveField_m, Time, FSInt, .TRUE., .TRUE., nodeInWater, WaveElev1, WaveElev2, WaveElev, FDynP, FV, FA, FAMCF, ErrStat2, ErrMsg2 )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        FDynPFSInt = REAL(FDynP,ReKi)
        FVFSInt    = REAL(FV,   ReKi)
        FAFSInt    = REAL(FA,   ReKi)
        IF ( mem%PropMCF .AND. ( .NOT. mem%PropPot ) ) THEN
           FAMCFFSInt = REAL(FAMCF,ReKi)
        END IF
        ! Structure translational acceleration at the free surface intersection
        SAFSInt = SubRatio  * u%Mesh%TranslationAcc(:,mem%NodeIndx(FSElem+1)) + &
             (1.0-SubRatio) * u%Mesh%TranslationAcc(:,mem%NodeIndx(FSElem  ))

        ! Viscous drag:
        ! Compute relative velocity at the free surface intersection. 
        ! Linear interpolation between the two nodes of the element is used to estimate velocity of the structure
        vrelFSInt = FVFSInt - ( & 
               SubRatio  * u%Mesh%TranslationVel(:,mem%NodeIndx(FSElem+1)) + &
          (1.0-SubRatio) * u%Mesh%TranslationVel(:,mem%NodeIndx(FSElem  ))   &
        )

        IF (mem%MSecGeom==MSecGeom_Cyl) THEN
           dRdl_p  = abs(mem%dRdl_mg(FSElem))
           dRdl_pp =     mem%dRdl_mg(FSElem)
           RMGFSInt = SubRatio * mem%RMG(FSElem+1) + (1.0-SubRatio) * mem%RMG(FSElem)

           vec = matmul( mem%Ak,vrelFSInt )
           F_DS = mem%Cd(FSElem)*p%WaveField%WtrDens*RMGFSInt*TwoNorm(vec)*vec  +  &
                     0.5*mem%AxCd(FSElem)*p%WaveField%WtrDens*pi*RMGFSInt*dRdl_p * &
                     abs(dot_product( mem%k, vrelFSInt )) * matmul( mem%kkt, vrelFSInt )

           ! Hydrodynamic added mass and inertia loads
           IF ( .NOT. mem%PropPot ) THEN

              ! ------------------- hydrodynamic added mass loads: sides: Section 7.1.3 ------------------------
              IF (p%AMMod > 0_IntKi) THEN
                 Am =      mem%Ca(FSElem)*p%WaveField%WtrDens*pi*RMGFSInt*RMGFSInt*mem%Ak + &
                     2.0*mem%AxCa(FSElem)*p%WaveField%WtrDens*pi*RMGFSInt*RMGFSInt*dRdl_p*mem%kkt
                 F_AS = -matmul( Am, &
                            SubRatio  * u%Mesh%TranslationAcc(:,mem%NodeIndx(FSElem+1)) + &
                       (1.0-SubRatio) * u%Mesh%TranslationAcc(:,mem%NodeIndx(FSElem  )) )
              END IF

              ! ------------------- hydrodynamic inertia loads: sides: Section 7.1.4 ------------------------
              IF ( mem%PropMCF) THEN
                 F_IS=                             p%WaveField%WtrDens*pi*RMGFSInt*RMGFSInt   * matmul( mem%Ak,  FAMCFFSInt ) + &
                              2.0*mem%AxCa(FSElem)*p%WaveField%WtrDens*pi*RMGFSInt*RMGFSInt*dRdl_p  * matmul( mem%kkt, FAFSInt ) + &
                              2.0*mem%AxCp(FSElem)          *pi*RMGFSInt                *dRdl_pp * FDynPFSInt*mem%k
              ELSE
                 F_IS=(mem%Ca(FSElem)+mem%Cp(FSElem))*p%WaveField%WtrDens*pi*RMGFSInt*RMGFSInt   * matmul( mem%Ak,  FAFSInt ) + &
                              2.0*mem%AxCa(FSElem)*p%WaveField%WtrDens*pi*RMGFSInt*RMGFSInt*dRdl_p  * matmul( mem%kkt, FAFSInt ) + &
                              2.0*mem%AxCp(FSElem)          *pi*RMGFSInt                *dRdl_pp * FDynPFSInt*mem%k
              END IF
           END IF
        ELSE IF (mem%MSecGeom==MSecGeom_Rec) THEN
           dSadl_p  = abs(mem%dSadl_mg(FSElem))
           dSadl_pp =     mem%dSadl_mg(FSElem)
           dSbdl_p  = abs(mem%dSbdl_mg(FSElem))
           dSbdl_pp =     mem%dSbdl_mg(FSElem)
           SaMGFSInt = SubRatio * mem%SaMG(FSElem+1) + (1.0-SubRatio) * mem%SaMG(FSElem)
           SbMGFSInt = SubRatio * mem%SbMG(FSElem+1) + (1.0-SubRatio) * mem%SbMG(FSElem)

           vec = matmul( mem%Ak,vrelFSInt )
           F_DS = 0.5*mem%CdB(FSElem)*p%WaveField%WtrDens*SbMGFSInt*TwoNorm(vec)*Dot_Product(vec,mem%x_hat)*mem%x_hat  +  &       ! local x-direction
                  0.5*mem%CdA(FSElem)*p%WaveField%WtrDens*SaMGFSInt*TwoNorm(vec)*Dot_Product(vec,mem%y_hat)*mem%y_hat  +  &       ! local z-direction
                  0.25*mem%AxCd(FSElem)*p%WaveField%WtrDens * (dSadl_p*SbMGFSInt + dSbdl_p*SaMGFSInt) * &                         ! axial part
                  abs(dot_product( mem%k, vrelFSInt )) * matmul( mem%kkt, vrelFSInt )                                             ! axial part cont'd

           ! Hydrodynamic added mass and inertia loads
           IF ( .NOT. mem%PropPot ) THEN

              ! ------------------- hydrodynamic added mass loads: sides: Section 7.1.3 ------------------------
              IF (p%AMMod > 0_IntKi) THEN
                 F_AS = -p%WaveField%WtrDens*mem%CaB(FSElem) * 0.25*pi*SbMGFSInt*SbMGFSInt * Dot_Product(SAFSInt,mem%x_hat)*mem%x_hat &
                        -p%WaveField%WtrDens*mem%CaA(FSElem) * 0.25*pi*SaMGFSInt*SaMGFSInt * Dot_Product(SAFSInt,mem%y_hat)*mem%y_hat &
                    -0.5*p%WaveField%WtrDens*mem%AxCa(FSElem) * (dSbdl_p*SaMGFSInt+dSadl_p*SbMGFSInt)*SQRT(SaMGFSInt*SbMGFSInt) * Dot_Product(SAFSInt,mem%k)*mem%k
              END IF
         
              ! ------------------- hydrodynamic inertia loads: sides: Section 7.1.4 ------------------------
              F_IS= mem%Cp(FSElem)*p%WaveField%WtrDens* SaMGFSInt*SbMGFSInt * matmul( mem%Ak,  FAFSInt ) + &                            ! transver FK component
                    FDynPFSInt*mem%AxCp(FSElem)* (SaMGFSInt*dSbdl_pp+dSadl_pp*SbMGFSInt) *mem%k + &                                     ! axial FK component
                    p%WaveField%WtrDens*mem%CaB(FSElem) * 0.25*pi*SbMGFSInt*SbMGFSInt * Dot_Product(FAFSInt,mem%x_hat)*mem%x_hat + &    ! x-component of diffraction part
                    p%WaveField%WtrDens*mem%CaA(FSElem) * 0.25*pi*SaMGFSInt*SaMGFSInt * Dot_Product(FAFSInt,mem%y_hat)*mem%y_hat + &    ! y-component of diffraction part
                0.5*p%WaveField%WtrDens*mem%AxCa(FSElem) * (dSbdl_p*SaMGFSInt+dSadl_p*SbMGFSInt)*SQRT(SaMGFSInt*SbMGFSInt) * &          ! axial component of diffraction part
                    Dot_Product(FAFSInt,mem%k)*mem%k                                                                                    ! axial component of diffraction part cont'd

           END IF
        END IF
        !----------------------------------------------------------------------------------------------------!
        !                         Perform the load redistribution for smooth time series                     !
        !----------------------------------------------------------------------------------------------------!
        ! Evaluate the load redistribution function
        f_redist = 0.0_ReKi
        IF (FSElem > 1_IntKi) THEN ! At least one fully submerged element
           f_redist = 2.0_ReKi * SubRatio**3 - 3.5_ReKi * SubRatio**2 + SubRatio + 0.5_ReKi
        END IF
        
        ! deltal = mem%dl and h_c = 0 should always be used here by design. Moment correction will be applied separately
        deltal = mem%dl
        h_c    = 0.0_ReKi
        
        ! Viscous drag
        ! Apply load redistribution to the first node below the free surface
        Df_hydro = ((SubRatio-1.0_ReKi)/(2.0_ReKi)-f_redist)*F_D0 + SubRatio/2.0_ReKi*F_DS
        CALL LumpDistrHydroLoads( Df_hydro, mem%k, deltal, h_c, Df_hydro_lumped)
        m%memberLoads(im)%F_D(:, FSElem) = m%memberLoads(im)%F_D(:, FSElem) + Df_hydro_lumped
        y%Mesh%Force (:,mem%NodeIndx(FSElem)) = y%Mesh%Force (:,mem%NodeIndx(FSElem)) + Df_hydro_lumped(1:3)
        y%Mesh%Moment(:,mem%NodeIndx(FSElem)) = y%Mesh%Moment(:,mem%NodeIndx(FSElem)) + Df_hydro_lumped(4:6)
        
        ! Apply load redistribution to the second node below the free surface
        IF (FSElem > 1_IntKi) THEN ! Note: Only need to modify the loads on the second node below the free surface when there is at least one fully submerged element.
           Df_hydro = f_redist * F_D0
           CALL LumpDistrHydroLoads( Df_hydro, mem%k, deltal, h_c, Df_hydro_lumped)
           m%memberLoads(im)%F_D(:, FSElem-1) = m%memberLoads(im)%F_D(:, FSElem-1) + Df_hydro_lumped
           y%Mesh%Force (:,mem%NodeIndx(FSElem-1)) = y%Mesh%Force (:,mem%NodeIndx(FSElem-1)) + Df_hydro_lumped(1:3)
           y%Mesh%Moment(:,mem%NodeIndx(FSElem-1)) = y%Mesh%Moment(:,mem%NodeIndx(FSElem-1)) + Df_hydro_lumped(4:6)
        END IF

        ! Hydrodynamic added mass and inertia loads
        IF ( .NOT. mem%PropPot ) THEN
           
           IF ( p%AMMod > 0_IntKi ) THEN
              !-------------------- hydrodynamic added mass loads: sides: Section 7.1.3 ------------------------!
              ! Apply load redistribution to the first node below the free surface
              Df_hydro = ((SubRatio-1.0_ReKi)/(2.0_ReKi)-f_redist)*F_A0 + SubRatio/2.0_ReKi*F_AS
              CALL LumpDistrHydroLoads( Df_hydro, mem%k, deltal, h_c, Df_hydro_lumped)
              m%memberLoads(im)%F_A(:, FSElem) = m%memberLoads(im)%F_A(:, FSElem) + Df_hydro_lumped
              y%Mesh%Force (:,mem%NodeIndx(FSElem)) = y%Mesh%Force (:,mem%NodeIndx(FSElem)) + Df_hydro_lumped(1:3)
              y%Mesh%Moment(:,mem%NodeIndx(FSElem)) = y%Mesh%Moment(:,mem%NodeIndx(FSElem)) + Df_hydro_lumped(4:6)
        
              ! Apply load redistribution to the second node below the free surface
              IF (FSElem > 1_IntKi) THEN
                  Df_hydro = f_redist * F_A0
                  CALL LumpDistrHydroLoads( Df_hydro, mem%k, deltal, h_c, Df_hydro_lumped)
                  m%memberLoads(im)%F_A(:, FSElem-1) = m%memberLoads(im)%F_A(:, FSElem-1) + Df_hydro_lumped
                  y%Mesh%Force (:,mem%NodeIndx(FSElem-1)) = y%Mesh%Force (:,mem%NodeIndx(FSElem-1)) + Df_hydro_lumped(1:3)
                  y%Mesh%Moment(:,mem%NodeIndx(FSElem-1)) = y%Mesh%Moment(:,mem%NodeIndx(FSElem-1)) + Df_hydro_lumped(4:6)
              END IF
           END IF
           
           !-------------------- hydrodynamic inertia loads: sides: Section 7.1.4 --------------------------!
           ! Apply load redistribution to the first node below the free surface
           Df_hydro = ((SubRatio-1.0_ReKi)/(2.0_ReKi)-f_redist)*F_I0 + SubRatio/2.0_ReKi*F_IS
           CALL LumpDistrHydroLoads( Df_hydro, mem%k, deltal, h_c, Df_hydro_lumped)
           m%memberLoads(im)%F_I(:, FSElem) = m%memberLoads(im)%F_I(:, FSElem) + Df_hydro_lumped
           y%Mesh%Force (:,mem%NodeIndx(FSElem)) = y%Mesh%Force (:,mem%NodeIndx(FSElem)) + Df_hydro_lumped(1:3)
           y%Mesh%Moment(:,mem%NodeIndx(FSElem)) = y%Mesh%Moment(:,mem%NodeIndx(FSElem)) + Df_hydro_lumped(4:6)
        
           ! Apply load redistribution to the second node below the free surface
           IF (FSElem > 1_IntKi) THEN
               Df_hydro = f_redist * F_I0
               CALL LumpDistrHydroLoads( Df_hydro, mem%k, deltal, h_c, Df_hydro_lumped)
               m%memberLoads(im)%F_I(:, FSElem-1) = m%memberLoads(im)%F_I(:, FSElem-1) + Df_hydro_lumped
               y%Mesh%Force (:,mem%NodeIndx(FSElem-1)) = y%Mesh%Force (:,mem%NodeIndx(FSElem-1)) + Df_hydro_lumped(1:3)
               y%Mesh%Moment(:,mem%NodeIndx(FSElem-1)) = y%Mesh%Moment(:,mem%NodeIndx(FSElem-1)) + Df_hydro_lumped(4:6)
           END IF
        END IF

        !----------------------------------------------------------------------------------------------------!
        !                     Perform moment correction to compensate for load redistribution                !
        !----------------------------------------------------------------------------------------------------!
        ! Moment correction to the first and second nodes below the free surface
        F_S = F_DS
        F_0 = F_D0
        IF ( .NOT. mem%PropPot) THEN
           F_S = F_S + F_IS
           F_0 = F_0 + F_I0
           IF ( p%AMMod > 0_IntKi) THEN
              F_S = F_S + F_AS
              F_0 = F_0 + F_A0
           END IF
        END IF
        ! First node below the free surface
        DM_hydro = 0.5_ReKi * SubRatio**2 * deltal * cross_product(mem%k, F_S)
        y%Mesh%Moment(:,mem%NodeIndx(FSElem))   = y%Mesh%Moment(:,mem%NodeIndx(FSElem))   + DM_hydro * deltal
        ! Second node below the free surface
        IF (FSElem > 1_IntKi) THEN
            DM_hydro =               f_redist * deltal * cross_product(mem%k, F_0)
            y%Mesh%Moment(:,mem%NodeIndx(FSElem-1)) = y%Mesh%Moment(:,mem%NodeIndx(FSElem-1)) + DM_hydro * deltal
        END IF

      ELSE IF ( MemSubStat .NE. 3_IntKi) THEN ! Skip members with centerline completely out of water
        !----------------------------No load smoothing----------------------------!
        DO i = mem%i_floor+1,N+1    ! loop through member nodes starting from the first node above seabed
           z1 = m%DispNodePosHdn(3, mem%NodeIndx(i))
           !---------------------------------------------Compute deltal and h_c------------------------------------------!
           ! Cannot make any assumption about WaveStMod and member orientation 
           IF ( m%NodeInWater(mem%NodeIndx(i)) .EQ. 0_IntKi ) THEN ! Node is out of water
              deltal = 0.0_ReKi
              h_c    = 0.0_ReKi
           ELSE ! Node in water
              ! Look to the "left" toward node 1
              IF ( i == 1 ) THEN ! First node. Note: Having i == 1 also implies mem%i_floor = 0.
                 deltalLeft = 0.0_ReKi
              ELSE IF ( i == mem%i_floor+1 ) THEN ! First node above seabed.
                 ! Note: This part is superceded by i==1 above when mem%i_floor = 0.
                 !       This is the correct behavior.
                 deltalLeft = -mem%h_floor
              ELSE ! Regular internal node
                 IF ( m%NodeInWater(mem%NodeIndx(i-1)) .EQ. 1_IntKi ) THEN ! Node to the left is submerged
                    deltalLeft = 0.5_ReKi * mem%dl
                 ELSE ! Element i-1 crosses the free surface
                    z2 = m%DispNodePosHdn(3, mem%NodeIndx(i-1))
                    IF ( p%WaveField%WaveStMod > 0_IntKi ) THEN ! Wave stretching enabled
                       zeta1 = m%WaveElev(mem%NodeIndx(i  ))
                       zeta2 = m%WaveElev(mem%NodeIndx(i-1))
                    ELSE
                       zeta1 = 0.0_ReKi
                       zeta2 = 0.0_ReKi
                    END IF
                    SubRatio = (zeta1-z1)/((zeta1-z1)-(zeta2-z2))
                    deltalLeft = SubRatio * mem%dl ! Portion of element i-1 in water
                 END IF
              END IF
              ! Look to the "right" toward node N+1
              IF ( i == N+1 ) THEN ! Last node
                 deltalRight = 0.0_ReKi
              ELSE ! Regular internal node
                 IF ( m%NodeInWater(mem%NodeIndx(i+1)) .EQ. 1_IntKi ) THEN ! Node to the right is submerged
                    deltalRight = 0.5_ReKi * mem%dl
                 ELSE ! Element i crosses the free surface
                    z2 = m%DispNodePosHdn(3, mem%NodeIndx(i+1))
                    IF ( p%WaveField%WaveStMod > 0_IntKi ) THEN ! Wave stretching enabled
                       zeta1 = m%WaveElev(mem%NodeIndx(i  ))
                       zeta2 = m%WaveElev(mem%NodeIndx(i+1))
                    ELSE
                       zeta1 = 0.0_ReKi
                       zeta2 = 0.0_ReKi
                    END IF
                    SubRatio = (zeta1-z1)/((zeta1-z1)-(zeta2-z2))
                    deltalRight = SubRatio * mem%dl ! Portion of element i in water
                 END IF
              END IF
              ! Combine left and right contributions
              deltal =              deltalRight + deltalLeft
              h_c    = 0.5_ReKi * ( deltalRight - deltalLeft )
           END IF

           ! Compute the slope of member radius/side length
           IF (mem%MSecGeom==MSecGeom_Cyl) THEN
              IF (i == 1) THEN
                 dRdl_p  = abs(mem%dRdl_mg(i))
                 dRdl_pp = mem%dRdl_mg(i)
              ELSE IF ( i > 1 .AND. i < (N+1)) THEN
                 dRdl_p  = 0.5*( abs(mem%dRdl_mg(i-1)) + abs(mem%dRdl_mg(i)) )
                 dRdl_pp = 0.5*( mem%dRdl_mg(i-1) + mem%dRdl_mg(i) )
              ELSE
                 dRdl_p  = abs(mem%dRdl_mg(N))
                 dRdl_pp = mem%dRdl_mg(N)
              END IF
           ELSE IF (mem%MSecGeom==MSecGeom_Rec) THEN
              IF (i == 1) THEN
                 dSadl_p  = abs(mem%dSadl_mg(i))
                 dSadl_pp = mem%dSadl_mg(i)
                 dSbdl_p  = abs(mem%dSbdl_mg(i))
                 dSbdl_pp = mem%dSbdl_mg(i)
              ELSE IF ( i > 1 .AND. i < (N+1)) THEN
                 dSadl_p  = 0.5*( abs(mem%dSadl_mg(i-1)) + abs(mem%dSadl_mg(i)) )
                 dSadl_pp = 0.5*( mem%dSadl_mg(i-1) + mem%dSadl_mg(i) )
                 dSbdl_p  = 0.5*( abs(mem%dSbdl_mg(i-1)) + abs(mem%dSbdl_mg(i)) )
                 dSbdl_pp = 0.5*( mem%dSbdl_mg(i-1) + mem%dSbdl_mg(i) )
              ELSE
                 dSadl_p  = abs(mem%dSadl_mg(N))
                 dSadl_pp = mem%dSadl_mg(N)
                 dSbdl_p  = abs(mem%dSbdl_mg(N))
                 dSbdl_pp = mem%dSbdl_mg(N)
              END IF
           END IF
         
           !--------------------- hydrodynamic drag loads: sides: Section 7.1.2 --------------------------------! 
           vec = matmul( mem%Ak,m%vrel(:,mem%NodeIndx(i)) )
           IF (mem%MSecGeom==MSecGeom_Cyl) THEN
              f_hydro = mem%Cd(i)*p%WaveField%WtrDens*mem%RMG(i)*TwoNorm(vec)*vec  +  &                                              ! radial part
                        0.5*mem%AxCd(i)*p%WaveField%WtrDens*pi*mem%RMG(i)*dRdl_p * &                                                 ! axial part
                        abs(dot_product( mem%k, m%vrel(:,mem%NodeIndx(i)) )) * matmul( mem%kkt, m%vrel(:,mem%NodeIndx(i)) )          ! axial part cont'd
           ELSE IF (mem%MSecGeom==MSecGeom_Rec) THEN
              f_hydro = 0.5*mem%CdB(i)*p%WaveField%WtrDens*mem%SbMG(i)*TwoNorm(vec)*Dot_Product(vec,mem%x_hat)*mem%x_hat  +  &       ! local x-direction
                        0.5*mem%CdA(i)*p%WaveField%WtrDens*mem%SaMG(i)*TwoNorm(vec)*Dot_Product(vec,mem%y_hat)*mem%y_hat  +  &       ! local z-direction
                        0.25*mem%AxCd(i)*p%WaveField%WtrDens * (dSadl_p*mem%SbMG(i) + dSbdl_p*mem%SaMG(i)) * &                       ! axial part
                        abs(dot_product( mem%k, m%vrel(:,mem%NodeIndx(i)) )) * matmul( mem%kkt, m%vrel(:,mem%NodeIndx(i)) )          ! axial part cont'd
           END IF
           CALL LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%memberLoads(im)%F_D(:, i) )
           y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + m%memberLoads(im)%F_D(1:3, i)
           y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + m%memberLoads(im)%F_D(4:6, i)
            
           IF ( .NOT. mem%PropPot ) THEN
              !-------------------- hydrodynamic added mass loads: sides: Section 7.1.3 ------------------------!
              IF (mem%MSecGeom==MSecGeom_Cyl) THEN
                 Am = mem%Ca(i)*p%WaveField%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*mem%Ak + 2.0*mem%AxCa(i)*p%WaveField%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*dRdl_p*mem%kkt
                 f_hydro = -matmul( Am, u%Mesh%TranslationAcc(:,mem%NodeIndx(i)) )
              ELSE IF (mem%MSecGeom==MSecGeom_Rec) THEN
                 f_hydro = -p%WaveField%WtrDens*mem%CaB(i) * 0.25*pi*mem%SbMG(i)*mem%SbMG(i) * Dot_Product(u%Mesh%TranslationAcc(:,mem%NodeIndx(i)),mem%x_hat)*mem%x_hat &
                           -p%WaveField%WtrDens*mem%CaA(i) * 0.25*pi*mem%SaMG(i)*mem%SaMG(i) * Dot_Product(u%Mesh%TranslationAcc(:,mem%NodeIndx(i)),mem%y_hat)*mem%y_hat &
                       -0.5*p%WaveField%WtrDens*mem%AxCa(i) * (dSbdl_p*mem%SaMG(i)+dSadl_p*mem%SbMG(i))*SQRT(mem%SaMG(i)*mem%SbMG(i)) * Dot_Product(u%Mesh%TranslationAcc(:,mem%NodeIndx(i)),mem%k)*mem%k
              END IF
              IF ( p%AMMod .EQ. 0_IntKi ) THEN ! Always compute added-mass force on nodes below SWL when undisplaced
                 z1 = u%Mesh%Position(3, mem%NodeIndx(i)) - p%WaveField%MSL2SWL ! Undisplaced z-position of the current node
                 IF ( z1 > 0.0_ReKi ) THEN ! Node is above SWL when undisplaced; zero added-mass force
                    f_hydro = 0.0_ReKi
                    CALL LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%memberLoads(im)%F_A(:, i) )
                 ELSE ! Node at or below SWL when undisplaced
                    IF ( i == 1 ) THEN
                       deltalLeft = 0.0_ReKi
                    ELSE IF ( i == mem%i_floor+1 ) THEN
                       deltalLeft = -mem%h_floor
                    ELSE
                       deltalLeft = 0.5_ReKi * mem%dl
                    END IF
                    IF ( i == N+1 ) THEN
                       deltalRight = 0.0_ReKi
                    ELSE
                       z2 = u%Mesh%Position(3, mem%NodeIndx(i+1)) - p%WaveField%MSL2SWL
                       IF ( z2 > 0.0_ReKi ) THEN ! Element i crosses the SWL
                          deltalRight = -z1 / mem%cosPhi_ref
                       ELSE
                          deltalRight = 0.5_ReKi * mem%dl
                       END IF
                    END IF
                    deltal_AM =              deltalRight + deltalLeft
                    h_c_AM    = 0.5_ReKi * ( deltalRight - deltalLeft )
                    CALL LumpDistrHydroLoads( f_hydro, mem%k, deltal_AM, h_c_AM, m%memberLoads(im)%F_A(:, i) )
                 END IF
              ELSE ! Compute added-mass force on the instantaneous wetted section of the member
                 f_hydro = f_hydro * m%nodeInWater(mem%NodeIndx(i)) ! Zero the force if node above free surface
                 CALL LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%memberLoads(im)%F_A(:, i) )
              END IF ! AMMod 0 or 1
              y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + m%memberLoads(im)%F_A(1:3, i)
              y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + m%memberLoads(im)%F_A(4:6, i)
              
              !-------------------- hydrodynamic inertia loads: sides: Section 7.1.4 ---------------------------!
              IF (mem%MSecGeom==MSecGeom_Cyl) THEN
                 IF ( mem%PropMCF ) THEN
                    f_hydro=                     p%WaveField%WtrDens*pi*mem%RMG(i)*mem%RMG(i)        * matmul( mem%Ak,  m%FAMCF(:,mem%NodeIndx(i)) ) + &
                                 2.0*mem%AxCa(i)*p%WaveField%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*dRdl_p * matmul( mem%kkt, m%FA(:,mem%NodeIndx(i)) ) + &
                                 2.0*m%FDynP(mem%NodeIndx(i))*mem%AxCp(i)*pi*mem%RMG(i)*dRdl_pp*mem%k
                 ELSE
                    f_hydro=(mem%Ca(i)+mem%Cp(i))*p%WaveField%WtrDens*pi*mem%RMG(i)*mem%RMG(i)       * matmul( mem%Ak,  m%FA(:,mem%NodeIndx(i)) ) + &
                                 2.0*mem%AxCa(i) *p%WaveField%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*dRdl_p * matmul( mem%kkt, m%FA(:,mem%NodeIndx(i)) ) + &
                                 2.0*m%FDynP(mem%NodeIndx(i))*mem%AxCp(i)*pi*mem%RMG(i)*dRdl_pp*mem%k
                 END IF
              ELSE IF (mem%MSecGeom==MSecGeom_Rec) THEN
                 ! Note: MacCamy-Fuchs correction cannot be applied to rectangular members
                 f_hydro= mem%Cp(i)*p%WaveField%WtrDens* mem%SaMG(i)*mem%SbMG(i) * matmul( mem%Ak,  m%FA(:,mem%NodeIndx(i)) ) + &                            ! transver FK component
                          m%FDynP(mem%NodeIndx(i))*mem%AxCp(i)* (mem%SaMG(i)*dSbdl_pp+dSadl_pp*mem%SbMG(i)) *mem%k + &                                       ! axial FK component
                          p%WaveField%WtrDens*mem%CaB(i) * 0.25*pi*mem%SbMG(i)*mem%SbMG(i) * Dot_Product(m%FA(:,mem%NodeIndx(i)),mem%x_hat)*mem%x_hat + &    ! x-component of diffraction part
                          p%WaveField%WtrDens*mem%CaA(i) * 0.25*pi*mem%SaMG(i)*mem%SaMG(i) * Dot_Product(m%FA(:,mem%NodeIndx(i)),mem%y_hat)*mem%y_hat + &    ! y-component of diffraction part
                      0.5*p%WaveField%WtrDens*mem%AxCa(i) * (dSbdl_p*mem%SaMG(i)+dSadl_p*mem%SbMG(i))*SQRT(mem%SaMG(i)*mem%SbMG(i)) * &                      ! axial component of diffraction part
                          Dot_Product(m%FA(:,mem%NodeIndx(i)),mem%k)*mem%k                                                                                   ! axial component of diffraction part cont'd
              END IF
              CALL LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%memberLoads(im)%F_I(:, i) )
              y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + m%memberLoads(im)%F_I(1:3, i)
              y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + m%memberLoads(im)%F_I(4:6, i)
           END IF

        END DO ! i = 1,N+1    ! loop through member nodes       
      
      END IF ! Check if the member is surface piercing
      !-----------------------------------------------------------------------------------------------------!
      !                                External Hydrodynamic Side Loads - End                               !
      !-----------------------------------------------------------------------------------------------------!

      !-----------------------------------------------------------------------------------------------------!
      !             Any end plate loads that are modeled on a per-member basis: F_B and F_BF                !
      !-----------------------------------------------------------------------------------------------------!
      ! reassign convenience variables to correspond to member ends
      ! We need to subtract the MSL2SWL offset to place this  in the SWL reference system
      pos1 = m%DispNodePosHst(:,mem%NodeIndx(1))
      pos2 = m%DispNodePosHst(:,mem%NodeIndx(2))
      call GetOrientationAngles( pos1, pos2, phi1, sinPhi1, cosPhi1, tanPhi, sinBeta1, cosBeta1, k_hat1, errStat2, errMsg2 )
        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if ( N == 1 ) then       ! Only one element in member
         sinPhi2  = sinPhi1
         cosPhi2  = cosPhi1
         sinBeta2 = sinBeta1
         cosBeta2 = cosBeta1
         k_hat2   = k_hat1
      else
         !  We need to subtract the MSL2SWL offset to place this  in the SWL reference system
         pos1 = m%DispNodePosHst(:, mem%NodeIndx(N  ))
         pos2 = m%DispNodePosHst(:, mem%NodeIndx(N+1))
         call GetOrientationAngles( pos1, pos2, phi2, sinPhi2, cosPhi2, tanPhi, sinBeta2, cosBeta2, k_hat2, errStat2, errMsg2 )
           call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      end if
      ! z-coordinates of the two ends of the member
      z1 = m%DispNodePosHst(3,mem%NodeIndx(  1))
      z2 = m%DispNodePosHst(3,mem%NodeIndx(N+1))
      
      if (mem%MSecGeom == MSecGeom_Rec) then
         ! Compute total orientation matrix of starting and ending joints
         call Morison_DirCosMtrx( u%Mesh%Position(:,mem%NodeIndx(1)), u%Mesh%Position(:,mem%NodeIndx(N+1)), mem%MSpinOrient, CMatrix )
         CMatrix1 = matmul(transpose(u%Mesh%Orientation(:,:,mem%NodeIndx(1  ))),CMatrix)
         CALL GetSectionUnitVectors_Rec( CMatrix1, x_hat1, y_hat1 )
         CMatrix2 = matmul(transpose(u%Mesh%Orientation(:,:,mem%NodeIndx(N+1))),CMatrix)
         CALL GetSectionUnitVectors_Rec( CMatrix2, x_hat2, y_hat2 )
      end if


      !----------------------------------- filled buoyancy loads: starts -----------------------------------!
      if ( mem%memfloodstatus > 0 ) then

         if ( mem%i_floor == 0 ) then                                                ! If the member is not buried in the seabed, compute the internal hydrostatic load on the starting endplate
            if ( mem%MSecGeom == MSecGeom_Cyl ) then
               m%F_BF_End(1:3, mem%NodeIndx(  1)) = m%F_BF_End(1:3, mem%NodeIndx(  1)) - mem%FillDens * g *        pi * mem%Rin(  1)**2* (zFillGroup - z1) * k_hat1
               m%F_BF_End(4:6, mem%NodeIndx(  1)) = m%F_BF_End(4:6, mem%NodeIndx(  1)) - mem%FillDens * g * 0.25 * pi * mem%Rin(  1)**4* Cross_Product(k_hat1,(/0.0,0.0,1.0/))
            else if ( mem%MSecGeom == MSecGeom_Rec ) then
               m%F_BF_End(1:3, mem%NodeIndx(  1)) = m%F_BF_End(1:3, mem%NodeIndx(  1)) - mem%FillDens * g *         mem%SaIn(  1)   *mem%SbIn(  1)* (zFillGroup - z1) * k_hat1
               m%F_BF_End(4:6, mem%NodeIndx(  1)) = m%F_BF_End(4:6, mem%NodeIndx(  1)) - mem%FillDens * g / 12.0 * (mem%SaIn(  1)**3*mem%SbIn(  1)*x_hat1(3)*y_hat1 - mem%SaIn(1)*mem%SbIn(1)**3*y_hat1(3)*x_hat1)
            end if
         end if

         if ( (mem%i_floor<mem%NElements+1) .and. (mem%memfloodstatus==1) ) then     ! If the member is not fully buried in the seabed and fully filled, compute the internal hydrostatic load on the ending endplate
            ! Note: If member is not fully filled, the endplate load is added to the appropriate member internal node under F_BF above
            if ( mem%MSecGeom == MSecGeom_Cyl ) then
               m%F_BF_End(1:3, mem%NodeIndx(N+1)) = m%F_BF_End(1:3, mem%NodeIndx(N+1)) + mem%FillDens * g *        pi * mem%Rin(N+1)**2* (zFillGroup - z2) * k_hat2
               m%F_BF_End(4:6, mem%NodeIndx(N+1)) = m%F_BF_End(4:6, mem%NodeIndx(N+1)) + mem%FillDens * g * 0.25 * pi * mem%Rin(N+1)**4* Cross_Product(k_hat2,(/0.0,0.0,1.0/))
            else if ( mem%MSecGeom == MSecGeom_Rec ) then
               m%F_BF_End(1:3, mem%NodeIndx(N+1)) = m%F_BF_End(1:3, mem%NodeIndx(N+1)) + mem%FillDens * g *         mem%SaIn(N+1)   *mem%SbIn(N+1)* (zFillGroup - z2) * k_hat2
               m%F_BF_End(4:6, mem%NodeIndx(N+1)) = m%F_BF_End(4:6, mem%NodeIndx(N+1)) + mem%FillDens * g / 12.0 * (mem%SaIn(N+1)**3*mem%SbIn(N+1)*x_hat2(3)*y_hat2 - mem%SaIn(N+1)*mem%SbIn(N+1)**3*y_hat2(3)*x_hat2)
            end if
         end if

       end if

      !------------------------------------ filled buoyancy loads: ends ------------------------------------!

      ! --- no inertia loads from water ballast modeled on ends

      !---------------------------------- external buoyancy loads: starts ----------------------------------!
      if ( (.not. mem%PropPot) .AND. (mem%MHstLMod /= 0) ) then
         ! Get positions and scaled radii of member end nodes
         pos1 = m%DispNodePosHst(:,mem%NodeIndx(  1))
         pos2 = m%DispNodePosHst(:,mem%NodeIndx(N+1))
         if (mem%MSecGeom==MSecGeom_Cyl) then
            r1      = mem%RMGB(  1)
            r2      = mem%RMGB(N+1)
            rn1     = r1
            rn2     = r2
         else if (mem%MSecGeom==MSecGeom_Rec) then
            Sa1     = mem%SaMGB(  1)
            Sa2     = mem%SaMGB(N+1)
            Sb1     = mem%SbMGB(  1)
            Sb2     = mem%SbMGB(N+1)
            rn1     = MAX(Sa1,Sb1)
            rn2     = MAX(Sa2,Sb2)
         end if
         if (mem%i_floor == 0) then  ! both ends above or at seabed
            ! Compute loads on the end plate of node 1
            IF (p%WaveField%WaveStMod > 0) THEN
               CALL GetTotalWaveElev( Time, pos1, Zeta1, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               CALL GetFreeSurfaceNormal( Time, pos1, rn1, n_hat, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FSPt = (/pos1(1),pos1(2),Zeta1/) ! Reference point on the free surface
            ELSE
               FSPt = (/pos1(1),pos1(2),0.0_ReKi/)
               n_hat = (/0.0,0.0,1.0/)
            END IF

            if (mem%MSecGeom==MSecGeom_Cyl) then
               CALL GetSectionUnitVectors_Cyl( k_hat1, y_hat, z_hat )
               CALL GetSectionFreeSurfaceIntersects_Cyl( REAL(pos1,DbKi), REAL(FSPt,DbKi), k_hat1, y_hat, z_hat, n_hat, REAL(r1,DbKi), theta1, theta2, secStat)
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               CALL GetEndPlateHstLds_Cyl(pos1, k_hat1, y_hat, z_hat, r1, theta1, theta2, F_B_End)
               IF (mem%MHstLMod == 1) THEN ! Check for partially wetted end plates
                  IF ( .NOT.( EqualRealNos((theta2-theta1),0.0_DbKi) .OR. EqualRealNos((theta2-theta1),2.0_DbKi*PI_D) ) ) THEN
                      CALL SetErrStat(ErrID_Warn, 'End plate is partially wetted with MHstLMod = 1. The buoyancy load and distribution potentially have large error. This has happened to the first node of Member ID ' //trim(num2lstr(mem%MemberID)), errStat, errMsg, RoutineName )
                  END IF
               END IF
            else if (mem%MSecGeom==MSecGeom_Rec) then
               CALL GetSectionUnitVectors_Rec( CMatrix1, x_hat, y_hat )
               CALL GetEndPlateHstLds_Rec(pos1, k_hat1, x_hat, y_hat, Sa1, Sb1, FSPt, n_hat, F_B_End)
            end if
            m%F_B_End(:, mem%NodeIndx(  1)) = m%F_B_End(:, mem%NodeIndx(  1)) + F_B_End

            ! Compute loads on the end plate of node N+1
            IF (p%WaveField%WaveStMod > 0) THEN
               CALL GetTotalWaveElev( Time, pos2, Zeta2, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               CALL GetFreeSurfaceNormal( Time, pos2, rn2, n_hat, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FSPt = (/pos2(1),pos2(2),Zeta2/) ! Reference point on the free surface
            ELSE
               FSPt = (/pos2(1),pos2(2),0.0_ReKi/)
               n_hat = (/0.0,0.0,1.0/)
            END IF

            if (mem%MSecGeom==MSecGeom_Cyl) then
               CALL GetSectionUnitVectors_Cyl( k_hat2, y_hat, z_hat )
               CALL GetSectionFreeSurfaceIntersects_Cyl( REAL(pos2,DbKi), REAL(FSPt,DbKi), k_hat2, y_hat, z_hat, n_hat, REAL(r2,DbKi), theta1, theta2, secStat)
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               CALL GetEndPlateHstLds_Cyl(pos2, k_hat2, y_hat, z_hat, r2, theta1, theta2, F_B_End)
               IF (mem%MHstLMod == 1) THEN ! Check for partially wetted end plates
                  IF ( .NOT.( EqualRealNos((theta2-theta1),0.0_DbKi) .OR. EqualRealNos((theta2-theta1),2.0_DbKi*PI_D) ) ) THEN
                      CALL SetErrStat(ErrID_Warn, 'End plate is partially wetted with MHstLMod = 1. The buoyancy load and distribution potentially have large error. This has happened to the last node of Member ID ' //trim(num2lstr(mem%MemberID)), errStat, errMsg, RoutineName )
                  END IF
               END IF
            else if (mem%MSecGeom==MSecGeom_Rec) then
               CALL GetSectionUnitVectors_Rec( CMatrix2, x_hat, y_hat )
               CALL GetEndPlateHstLds_Rec(pos2, k_hat2, x_hat, y_hat, Sa2, Sb2, FSPt, n_hat, F_B_End)
            end if
            m%F_B_End(:, mem%NodeIndx(N+1)) = m%F_B_End(:, mem%NodeIndx(N+1)) - F_B_End

         elseif ( mem%doEndBuoyancy ) then ! The member crosses the seabed line so only the upper end potentially have hydrostatic load
            ! Only compute the loads on the end plate of node N+1
            IF (p%WaveField%WaveStMod > 0) THEN
               CALL GetTotalWaveElev( Time, pos2, Zeta2, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               CALL GetFreeSurfaceNormal( Time, pos2, rn2, n_hat, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FSPt = (/pos2(1),pos2(2),Zeta2/) ! Reference point on the free surface
            ELSE
               FSPt = (/pos2(1),pos2(2),0.0_ReKi/)
               n_hat = (/0.0,0.0,1.0/)
            END IF

            if (mem%MSecGeom==MSecGeom_Cyl) then
               CALL GetSectionUnitVectors_Cyl( k_hat2, y_hat, z_hat )
               CALL GetSectionFreeSurfaceIntersects_Cyl( REAL(pos2,DbKi), REAL(FSPt,DbKi), k_hat2, y_hat, z_hat, n_hat, REAL(r2,DbKi), theta1, theta2, secStat)
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               CALL GetEndPlateHstLds_Cyl(pos2, k_hat2, y_hat, z_hat, r2, theta1, theta2, F_B_End)
               IF (mem%MHstLMod == 1) THEN ! Check for partially wetted end plates
                  IF ( .NOT.( EqualRealNos((theta2-theta1),0.0_DbKi) .OR. EqualRealNos((theta2-theta1),2.0_DbKi*PI_D) ) ) THEN
                      CALL SetErrStat(ErrID_Warn, 'End plate is partially wetted with MHstLMod = 1. The buoyancy load and distribution potentially have large error. This has happened to the last node of Member ID ' //trim(num2lstr(mem%MemberID)), errStat, errMsg, RoutineName )
                  END IF
               END IF
            else if (mem%MSecGeom==MSecGeom_Rec) then
               CALL GetSectionUnitVectors_Rec( CMatrix2, x_hat, y_hat )
               CALL GetEndPlateHstLds_Rec(pos2, k_hat2, x_hat, y_hat, Sa2, Sb2, FSPt, n_hat, F_B_End)
            end if
            m%F_B_End(:, mem%NodeIndx(N+1)) = m%F_B_End(:, mem%NodeIndx(N+1)) - F_B_End

         ! else
            ! entire member is buried below the seabed
         end if

      end if   ! PropPot
      !----------------------------------- external buoyancy loads: ends -----------------------------------!

   end do ! im - looping through members

   !---------------------------------------------------------------------------------------------------------------!
   !                                     External Hydrodynamic Joint Loads - Start                                 !
   !                                        F_D_End, F_I_End, F_A_End, F_IMG_End                                   !
   !---------------------------------------------------------------------------------------------------------------!      
   ! NOTE:  All wave kinematics have already been zeroed out above the SWL or instantaneous wave height (for WaveStMod > 0), 
   ! so loads derived from the kinematics will be correct without the use of a nodeInWater value, but other loads need to be 
   ! multiplied by nodeInWater to zero them out above the SWL or instantaneous wave height.

   DO J = 1, p%NJoints
      ! Obtain the node index because WaveVel, WaveAcc, and WaveDynP are defined in the node indexing scheme, not the markers (No longer relevant?)
      ! The first NJoints nodes are all the joints with the rest being the internal nodes. See Morison_GenerateSimulationNodes.
            
      ! NOTE: 
      ! The PropPot values are only for members, and when the p%AM_End, p%DP_Const_End, p%Mass_MG_End, and p%I_MG_End are computed at init,
      ! contributions to these values are added only if the member connecting to the joint is NOT modeled with potential flow theory
      ! However, the p%An_End term used data from ALL members attached to a node, regardless of the PropPot setting, because the drag force is alway on.
      ! Therefore, no need to check PropPot here.
      
      ! Effect of wave stretching already baked into m%FDynP, m%FA, and m%vrel. No additional modification needed.

      ! Joint yaw offset
      call YawJoint(J,u%PtfmRefY,AM_End,An_End,DP_Const_End,I_MG_End,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      
      ! Lumped added mass loads
      qdotdot                 = reshape((/u%Mesh%TranslationAcc(:,J),u%Mesh%RotationAcc(:,J)/),(/6/)) 
      m%F_A_End(:,J)          = m%nodeInWater(j) * matmul( AM_End, ( - qdotdot(1:3)) )
         
      ! TODO: The original code did not multiply by nodeInWater, but should we? GJH
      ! Should be ok because m%FDynP and m%FA are both zeroed above the SWL (when WaveStMod=0) or the instantaneous free surface (when WaveStMod>0)
      m%F_I_End(:,J) =   (DP_Const_End * m%FDynP(j) + matmul(AM_End,m%FA(:,j)))
         
      ! Marine growth inertia: ends: Section 4.2.2
      m%F_IMG_End(1:3,j) = -p%Mass_MG_End(j)*qdotdot(1:3)
      m%F_IMG_End(4:6,j) = -matmul(I_MG_End,qdotdot(4:6)) - cross_product(u%Mesh%RotationVel(:,J),matmul(I_MG_End,u%Mesh%RotationVel(:,J)))

      ! Compute the dot product of the relative velocity vector with the directional Area of the Joint
      ! m%nodeInWater(j) is probably not necessary because m%vrel is zeroed when the node is out of water
      vmag  = m%nodeInWater(j) * ( m%vrel(1,j)*An_End(1) + m%vrel(2,j)*An_End(2) + m%vrel(3,j)*An_End(3) )
      ! High-pass filtering
      vmagf = p%VRelNFiltConst(J) * (vmag + xd%v_rel_n_FiltStat(J))

      ! Record most up-to-date vmagf and vmag at join J
      m%v_rel_n(j) = vmag
      m%v_rel_n_HiPass(j) = vmagf

      ! Evaluate drag force and combine all per-joint loads
      DO I=1,6
         IF (I < 4 ) THEN ! Three force components
            IF ( p%DragMod_End(J) .EQ. 0_IntKi ) THEN
               ! Note: vmag is zero if node is not in the water
               m%F_D_End(i,j) = (1.0_ReKi - p%DragLoFSc_End(j)) * An_End(i) * p%DragConst_End(j) * abs(vmagf)*vmagf &   
                                          + p%DragLoFSc_End(j)  * An_End(i) * p%DragConst_End(j) * abs(vmag )*vmag  
            ELSE IF (p%DragMod_End(J) .EQ. 1_IntKi) THEN
               ! Note: vmag is zero if node is not in the water
               m%F_D_End(i,j) = (1.0_ReKi - p%DragLoFSc_End(j)) * An_End(i) * p%DragConst_End(j) * abs(vmagf)*max(vmagf,0.0_ReKi) &
                                          + p%DragLoFSc_End(j)  * An_End(i) * p%DragConst_End(j) * abs(vmag) *max(vmag, 0.0_ReKi)
               m%F_D_End(i,j) = 2.0_ReKi * m%F_D_End(i,j)
            END IF
            
            y%Mesh%Force(i,j)    = y%Mesh%Force(i,j)    + m%F_D_End(i,j) + m%F_I_End(i,j) + p%F_WMG_End(i,j) + m%F_B_End(i,j) + m%F_BF_End(i,j) + m%F_A_End(i,j) + m%F_IMG_End(i,j)
         ELSE ! Three moment components
            y%Mesh%Moment(i-3,j) = y%Mesh%Moment(i-3,j) + m%F_B_End(i,j) + m%F_BF_End(i,j)  + m%F_IMG_End(i,j)
         END IF
      END DO  ! I=1,6
         
   END DO  ! J = 1, p%NJoints

   !---------------------------------------------------------------------------------------------------------------!
   !                                      External Hydrodynamic Joint Loads - End                                  !
   !---------------------------------------------------------------------------------------------------------------!    
   ! Map calculated results into the y%WriteOutput Array
   CALL MrsnOut_MapOutputs(y, p, u, m)


   ! map the motion to the visulization mesh
   if (p%VisMeshes) then
      !FIXME: error handling is incorrect here (overwrites all previous errors/warnings)
      call Transfer_Point_to_Line2( u%Mesh, y%VisMesh, m%VisMeshMap, ErrStat, ErrMsg )
   endif


   CONTAINS

   SUBROUTINE GetDisplacedNodePosition( forceDisplaced, pos )
      LOGICAL,         INTENT( IN    ) :: forceDisplaced ! Set to true to return the exact displaced position no matter WaveDisp or WaveStMod
      REAL(ReKi),      INTENT(   OUT ) :: pos(:,:) ! Displaced node positions
      REAL(ReKi)                       :: Orient(3,3)

      INTEGER(IntKi)                   :: ErrStat2
      CHARACTER(ErrMsgLen)             :: ErrMsg2

      ! Undisplaced node position
      pos      = u%Mesh%Position
      pos(3,:) = pos(3,:) - p%WaveField%MSL2SWL ! Z position measured from the SWL
      IF ( (p%WaveDisp /= 0) .OR. forceDisplaced ) THEN 
         ! Use displaced X and Y position
         pos(1,:) = pos(1,:) + u%Mesh%TranslationDisp(1,:)
         pos(2,:) = pos(2,:) + u%Mesh%TranslationDisp(2,:)
         IF ( (p%WaveField%WaveStMod > 0) .OR. forceDisplaced ) THEN
            ! Use displaced Z position only when wave stretching is enabled
            pos(3,:) = pos(3,:) + u%Mesh%TranslationDisp(3,:)
         END IF
      ELSE ! p%WaveDisp=0 implies PtfmYMod=0
         ! Rotate the structure based on PtfmRefY (constant)
         call GetPtfmRefYOrient(u%PtfmRefY, Orient, ErrStat2, ErrMsg2)
         pos = matmul(transpose(Orient),pos)
      END IF

   END SUBROUTINE GetDisplacedNodePosition

   SUBROUTINE GetTotalWaveElev( Time, pos, Zeta, ErrStat, ErrMsg )
      REAL(DbKi),      INTENT( IN    ) :: Time
      REAL(ReKi),      INTENT( IN    ) :: pos(*)  ! Position at which free-surface elevation is to be calculated. Third entry ignored if present.
      REAL(ReKi),      INTENT(   OUT ) :: Zeta    ! Total free-surface elevation with first- and second-order contribution (if present)
      INTEGER(IntKi),  INTENT(   OUT ) :: ErrStat ! Error status of the operation
      CHARACTER(*),    INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None
      CHARACTER(*),    PARAMETER       :: RoutineName = 'GetTotalWaveElev'
      INTEGER(IntKi)                   :: errStat2
      CHARACTER(ErrMsgLen)             :: errMsg2
      ErrStat   = ErrID_None
      ErrMsg    = ""

      Zeta = WaveField_GetNodeTotalWaveElev( p%WaveField, m%WaveField_m, Time, pos, ErrStat2, ErrMsg2 )
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   END SUBROUTINE GetTotalWaveElev

   SUBROUTINE GetFreeSurfaceNormal( Time, pos, r, n, ErrStat, ErrMsg)
      REAL(DbKi),      INTENT( In    ) :: Time
      REAL(ReKi),      INTENT( In    ) :: pos(:)  ! Position at which free-surface normal is to be calculated. Third entry ignored if present.
      REAL(ReKi),      INTENT( In    ) :: r       ! Distance for central differencing
      REAL(ReKi),      INTENT(   OUT ) :: n(3)    ! Free-surface normal vector
      INTEGER(IntKi),  INTENT(   OUT ) :: ErrStat ! Error status of the operation
      CHARACTER(*),    INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None
      CHARACTER(*),    PARAMETER       :: RoutineName = 'GetFreeSurfaceNormal'
      INTEGER(IntKi)                   :: errStat2
      CHARACTER(ErrMsgLen)             :: errMsg2
      ErrStat   = ErrID_None
      ErrMsg    = ""

      CALL WaveField_GetNodeWaveNormal( p%WaveField, m%WaveField_m, Time, pos, r, n, ErrStat2, ErrMsg2 )
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   END SUBROUTINE GetFreeSurfaceNormal

   SUBROUTINE GetSectionUnitVectors_Cyl( k, y, z )
      REAL(ReKi),      INTENT( In    ) :: k(3) ! Member axial unit vector
      REAL(ReKi),      INTENT(   OUT ) :: y(3) ! Horizontal unit vector perpendicular to k
      REAL(ReKi),      INTENT(   OUT ) :: z(3) ! Unit vector perpendicular to k and y with positive vertical component
      IF ( ABS(k(3)) > 0.999999_ReKi ) THEN ! k is effectively vertical
         y = (/0.0,1.0,0.0/)
      ELSE
         y = (/-k(2),k(1),0.0_ReKi/)
         y = y / SQRT(Dot_Product(y,y))      
      ENDIF
      z = cross_product(k,y)
      IF ( z(3) < 0.0 ) THEN ! Flip y and z so z points upward
         y = -y;
         z = -z;
      END IF
   END SUBROUTINE GetSectionUnitVectors_Cyl

   SUBROUTINE GetSectionUnitVectors_Rec( DirCos, x, y )
      REAL(ReKi),      INTENT( In    ) :: DirCos(3,3) ! Element local to global rotation matrix
      REAL(ReKi),      INTENT(   OUT ) :: x(3) ! Section local x unit vector parallel to side A
      REAL(ReKi),      INTENT(   OUT ) :: y(3) ! Section local y unit vector parallel to side B
      x = DirCos(1:3,1)
      y = DirCos(1:3,2)
   END SUBROUTINE GetSectionUnitVectors_Rec

   SUBROUTINE GetSectionFreeSurfaceIntersects_Cyl( pos0, FSPt, k_hat, y_hat, z_hat, n_hat, R, theta1, theta2, secStat)
      REAL(DbKi),      INTENT( In    ) :: pos0(3)
      REAL(DbKi),      INTENT( In    ) :: FSPt(3)
      REAL(ReKi),      INTENT( In    ) :: k_hat(3)
      REAL(ReKi),      INTENT( In    ) :: y_hat(3)
      REAL(ReKi),      INTENT( In    ) :: z_hat(3)
      REAL(ReKi),      INTENT( In    ) :: n_hat(3)
      REAL(DbKi),      INTENT( In    ) :: R
      REAL(DbKi),      INTENT(   OUT ) :: theta1
      REAL(DbKi),      INTENT(   OUT ) :: theta2
      INTEGER(IntKi),  INTENT(   OUT ) :: secStat
      REAL(DbKi)                       :: a, b, c, d, d2
      REAL(DbKi)                       :: alpha, beta
      REAL(DbKi)                       :: tmp
      CHARACTER(*),    PARAMETER       :: RoutineName = 'GetSectionFreeSurfaceIntersects_Cyl'

      a  = R * dot_product(y_hat,n_hat)
      b  = R * dot_product(z_hat,n_hat)
      c  = dot_product(FSPt-pos0,n_hat)
      d2 = a*a+b*b
      IF ( d2 >= c*c ) THEN ! Has intersection
         d = SQRT(d2)
         IF (b>=0.0) THEN
            alpha =  ACOS(a/d)
         ELSE
            alpha = -ACOS(a/d)
         END IF
         beta   = ACOS(c/d)
         theta1 = alpha - beta
         theta2 = alpha + beta
         IF ( dot_product( (cos(theta2)-cos(theta1))*z_hat-(sin(theta2)-sin(theta1))*y_hat, n_hat) < 0.0 ) THEN
            tmp    = theta1
            theta1 = theta2
            theta2 = tmp + 2.0*PI_D
         END IF
         secStat = 1;
      ELSE IF (c > 0.0) THEN ! Section is fully submerged
         theta1 = -1.5*PI_D
         theta2 =  0.5*PI_D
         secStat = 2;
      ELSE ! Section is completely dry
         theta1 = -0.5*PI_D
         theta2 = -0.5*PI_D
         secStat = 0;
      END IF

   END SUBROUTINE GetSectionFreeSurfaceIntersects_Cyl

   SUBROUTINE GetSectionHstLds_Cyl( origin, pos0, k_hat, y_hat, z_hat, R, dRdl, theta1, theta2, dFdl)

      REAL(DbKi),      INTENT( IN    ) :: origin(3)
      REAL(DbKi),      INTENT( IN    ) :: pos0(3)
      REAL(DbKi),      INTENT( IN    ) :: k_hat(3)
      REAL(DbKi),      INTENT( IN    ) :: y_hat(3)
      REAL(DbKi),      INTENT( IN    ) :: z_hat(3)
      REAL(DbKi),      INTENT( IN    ) :: R
      REAL(DbKi),      INTENT( IN    ) :: dRdl
      REAL(DbKi),      INTENT( IN    ) :: theta1
      REAL(DbKi),      INTENT( IN    ) :: theta2
      REAL(DbKi),      INTENT(   OUT ) :: dFdl(6)
      REAL(DbKi)                       :: C0, C1, C2
      REAL(DbKi)                       :: Z0, dTheta, sinTheta1, sinTheta2, cosTheta1, cosTheta2, cosPhi
      
      Z0        = pos0(3)
      dTheta    = theta2 - theta1
      sinTheta1 = SIN(theta1)
      sinTheta2 = SIN(theta2)
      cosTheta1 = COS(theta1)
      cosTheta2 = COS(theta2)
      cosPhi    = SQRT(k_hat(1)**2+k_hat(2)**2)

      C0 = Z0*dTheta                +     R*cosPhi*(cosTheta1   -cosTheta2)
      C1 = Z0*(sinTheta2-sinTheta1) + 0.5*R*cosPhi*(cosTheta2**2-cosTheta1**2)
      C2 = Z0*(cosTheta1-cosTheta2) + 0.5*R*cosPhi*(dTheta-sinTheta2*cosTheta2+sinTheta1*cosTheta1)

      dFdl(1:3) = -R   *dRdl*C0*k_hat + R*C1*y_hat + R*C2*z_hat
      dFdl(4:6) = -R**2*dRdl*C2*y_hat + R**2*dRdl*C1*z_hat      + CROSS_PRODUCT((pos0-origin),dFdl(1:3))
      dFdl = dFdl * p%WaveField%WtrDens * g

   END SUBROUTINE GetSectionHstLds_Cyl

   SUBROUTINE GetSectionHstLds_Rec( origin, pos0, k_hat, x_hat, y_hat, Sa, Sb, dSadl, dSbdl, rFS, nFS, dFdl, secStat)

      REAL(DbKi),      INTENT( IN    ) :: origin(3)
      REAL(DbKi),      INTENT( IN    ) :: pos0(3)
      REAL(DbKi),      INTENT( IN    ) :: k_hat(3)
      REAL(DbKi),      INTENT( IN    ) :: x_hat(3)
      REAL(DbKi),      INTENT( IN    ) :: y_hat(3)
      REAL(DbKi),      INTENT( IN    ) :: Sa, Sb
      REAL(DbKi),      INTENT( IN    ) :: dSadl, dSbdl
      REAL(DbKi),      INTENT( IN    ) :: rFS(3)
      REAL(DbKi),      INTENT( IN    ) :: nFS(3)
      REAL(DbKi),      INTENT(   OUT ) :: dFdl(6)
      INTEGER(IntKi),  INTENT(   OUT ) :: secStat

      INTEGER(IntKi)                   :: i, v1, v2, numVInWtr
      REAL(DbKi)                       :: s, s1, s2, sInt, x0, y0, z0, x1, y1, z1, x2, y2, z2
      REAL(DbKi)                       :: n(3,4), rv(3,4)
      REAL(DbKi)                       :: C(3)
      LOGICAL                          :: vInWtr(4)
      
      ! Origin point about which moment is computed
      x0 = origin(1)
      y0 = origin(2)
      z0 = origin(3)

      ! Global coordinates of the four vertices of the section
      rv(:,1) = pos0 + x_hat * (-0.5*Sa) + y_hat * (-0.5*Sb)
      rv(:,2) = pos0 + x_hat * ( 0.5*Sa) + y_hat * (-0.5*Sb)
      rv(:,3) = pos0 + x_hat * ( 0.5*Sa) + y_hat * ( 0.5*Sb)
      rv(:,4) = pos0 + x_hat * (-0.5*Sa) + y_hat * ( 0.5*Sb)

      ! Unit normal vector of side walls
      n(:,1)  =  y_hat + k_hat * 0.5 * dSbdl
      n(:,2)  = -x_hat + k_hat * 0.5 * dSadl
      n(:,3)  = -y_hat + k_hat * 0.5 * dSbdl
      n(:,4)  =  x_hat + k_hat * 0.5 * dSadl

      ! Check for and count vertices in water
      numVInWtr = 0
      do i = 1,4
         vInWtr(i) = ( Dot_Product(rv(:,i)-rFS,nFS) <= 0.0 )
         if ( vInWtr(i) ) then
            numVInWtr = numVInWtr + 1
         end if
      end do

      ! Return section status
      if (numVInWtr == 0) then            ! Completely dry section
         secStat = 0;
      else if (numVInWtr == 4) then       ! Completely submerged section
         secStat = 2;
      else                                ! Partially wetted section that intersects the free surface
         secStat = 1;
      end if

      ! Compute the hydrostatic force and moment on each of the 4 sides
      dFdl = 0.0
      do i = 1,4
         v1 = i
         if (i == 4) then
            v2 = 1
         else
            v2 = i+1
         end if
         x1 = rv(1,v1)
         y1 = rv(2,v1)
         z1 = rv(3,v1)
         x2 = rv(1,v2)
         y2 = rv(2,v2)
         z2 = rv(3,v2)
         if ( (i==1) .or. (i==3) ) then
            s = Sa
         else
            s = Sb
         end if
         if (EqualRealNos(s,0.0_DbKi)) cycle
         if ( vInWtr(v1) .and. vInWtr(v2) ) then
            ! Side fully submerged
            s1 = 0.0
            s2 = s;
         else if ( vInWtr(v1) .or. vInWtr(v2) ) then
            ! Side partially wetted
            sInt = s * DOT_PRODUCT(rFS-rv(:,v1),nFS) / DOT_PRODUCT(rv(:,v2)-rv(:,v1),nFS)
            if ( vInWtr(v1) ) then
               s1 = 0.0
               s2 = sInt
            else
               s1 = sInt
               s2 = s
            end if
         else
            ! Side fully out of water
            s1 = 0.0;
            s2 = 0.0;
         end if

         dFdl(1:3) = dFdl(1:3) + &
            -n(:,i) * ( z1*(s2-s1) + 0.5*(z2-z1)/s*(s2*s2-s1*s1) )
         C(1) = (z2-z1)*(x2-x1)/3.0/(s*s)*(s2**3-s1**3) + 0.5*((z2-z1)*(x1-x0)/s+(x2-x1)*z1/s)*(s2*s2-s1*s1) + z1*(x1-x0)*(s2-s1)
         C(2) = (z2-z1)*(y2-y1)/3.0/(s*s)*(s2**3-s1**3) + 0.5*((z2-z1)*(y1-y0)/s+(y2-y1)*z1/s)*(s2*s2-s1*s1) + z1*(y1-y0)*(s2-s1)
         C(3) = (z2-z1)*(z2-z1)/3.0/(s*s)*(s2**3-s1**3) + 0.5*((z2-z1)*(z1-z0)/s+(z2-z1)*z1/s)*(s2*s2-s1*s1) + z1*(z1-z0)*(s2-s1)
         dFdl(4:6) = dFdl(4:6) - CROSS_PRODUCT(C,n(:,i))

      end do

      dFdl = dFdl * p%WaveField%WtrDens * g

   END SUBROUTINE GetSectionHstLds_Rec

   SUBROUTINE getElementHstLds_Mod2_Cyl( pos1In, pos2In, FSPtIn, k_hatIn, y_hatIn, z_hatIn, n_hatIn, r1In, r2In, dlIn, F_B1, F_B2, ErrStat, ErrMsg )
      
      REAL(ReKi),      INTENT( IN    ) :: pos1In(3)
      REAL(ReKi),      INTENT( IN    ) :: pos2In(3)
      REAL(ReKi),      INTENT( IN    ) :: FSPtIn(3)
      REAL(ReKi),      INTENT( IN    ) :: k_hatIn(3)
      REAL(ReKi),      INTENT( IN    ) :: y_hatIn(3)
      REAL(ReKi),      INTENT( IN    ) :: z_hatIn(3)
      REAL(ReKi),      INTENT( IN    ) :: n_hatIn(3)
      REAL(ReKi),      INTENT( IN    ) :: r1In
      REAL(ReKi),      INTENT( IN    ) :: r2In
      REAL(ReKi),      INTENT( IN    ) :: dlIn
      REAL(ReKi),      INTENT(   OUT ) :: F_B1(6)
      REAL(ReKi),      INTENT(   OUT ) :: F_B2(6)
      INTEGER(IntKi),  INTENT(   OUT ) :: ErrStat ! Error status of the operation
      CHARACTER(*),    INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None
      REAL(DbKi)                       :: theta1, theta2
      REAL(DbKi)                       :: dFdl1(6), dFdlMid(6), dFdl2(6), F_B(6)
      REAL(DbKi)                       :: i, dl, r1, r2, rMid, dRdl, posMid(3), pos1(3), pos2(3), FSPt(3), k_hat(3), y_hat(3), z_hat(3), n_hat(3)
      INTEGER(IntKi)                   :: secStat1, secStatMid, secStat2
      CHARACTER(*),    PARAMETER       :: routineName = "getElementHstLds_Mod2_Cyl"
      INTEGER(IntKi)                   :: errStat2
      CHARACTER(ErrMsgLen)             :: errMsg2
      ErrStat   = ErrID_None
      ErrMsg    = ""  
  
      pos1   = REAL(pos1In,DbKi)
      pos2   = REAL(pos2In,DbKi)
      r1     = REAL(r1In,DbKi)
      r2     = REAL(r2In,DbKi)
      dl     = REAL(dlIn,DbKi)
      dRdl   = (r2-r1)/dl
      rMid   = 0.5*(  r1+  r2)
      posMid = 0.5*(pos1In+pos2In)
      FSPt   = REAL(FSPtIn,DbKi)
      k_hat  = REAL(k_hatIn,DbKi)
      y_hat  = REAL(y_hatIn,DbKi)
      z_hat  = REAL(z_hatIn,DbKi)
      n_hat  = REAL(n_hatIn,DbKi)
      
      ! Avoid sections coincident with the SWL
      IF ( ABS(k_hat(3)) > 0.999999_ReKi ) THEN ! Vertical member
         IF ( EqualRealNos( pos1(3), 0.0_DbKi ) ) THEN
            pos1(3) = pos1(3) - 1.0E-6 * dl
         END IF
         IF ( EqualRealNos( pos2(3), 0.0_DbKi ) ) THEN
            pos2(3) = pos2(3) - 1.0E-6 * dl
         END IF
         IF ( EqualRealNos( posMid(3), 0.0_DbKi ) ) THEN
            posMid(3) = posMid(3) - 1.0E-6 * dl
         END IF
      END IF

      ! Section load at node 1
      CALL GetSectionFreeSurfaceIntersects_Cyl( pos1,   FSPt, REAL(k_hat,ReKi), REAL(y_hat,ReKi), REAL(z_hat,ReKi), REAL(n_hat,ReKi), r1, theta1, theta2, secStat1)
      CALL GetSectionHstLds_Cyl( pos1, pos1,   k_hat, y_hat, z_hat, r1,   dRdl, theta1, theta2, dFdl1)

      ! Section load at midpoint
      CALL GetSectionFreeSurfaceIntersects_Cyl( posMid, FSPt, REAL(k_hat,ReKi), REAL(y_hat,ReKi), REAL(z_hat,ReKi), REAL(n_hat,ReKi), rMid, theta1, theta2, secStatMid)
      CALL GetSectionHstLds_Cyl( pos1, posMid, k_hat, y_hat, z_hat, rMid, dRdl, theta1, theta2, dFdlMid)

      ! Section load at node 2
      CALL GetSectionFreeSurfaceIntersects_Cyl( pos2,   FSPt, REAL(k_hat,ReKi), REAL(y_hat,ReKi), REAL(z_hat,ReKi), REAL(n_hat,ReKi), r2, theta1, theta2, secStat2)
      CALL GetSectionHstLds_Cyl( pos1, pos2,   k_hat, y_hat, z_hat, r2,   dRdl, theta1, theta2, dFdl2)

      ! Adaptively refine the load integration over the element
      CALL RefineElementHstLds_Cyl(pos1,pos1,posMid,pos2,FSPt,r1,rMid,r2,dl,dRdl,secStat1,secStatMid,secStat2,k_hat,y_hat,z_hat,n_hat,dFdl1,dFdlMid,dFdl2,1,F_B,ErrStat2,ErrMsg2)
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! Distribute the hydrostatic load to the two end nodes
      F_B1(1:3) = 0.5 * F_B(1:3)
      F_B2(1:3) = 0.5 * F_B(1:3)
      F_B(4:6)  = F_B(4:6) - CROSS_PRODUCT(k_hat*dl,F_B2(1:3))
      F_B1(4:6) = 0.5 * F_B(4:6)
      F_B2(4:6) = 0.5 * F_B(4:6)

   END SUBROUTINE getElementHstLds_Mod2_Cyl

   SUBROUTINE getElementHstLds_Mod2_Rec( pos1In, pos2In, FSPtIn, k_hatIn, x_hatIn, y_hatIn, n_hatIn, Sa1In, Sa2In, Sb1In, Sb2In, dlIn, F_B1, F_B2, ErrStat, ErrMsg )
      
      REAL(ReKi),      INTENT( IN    ) :: pos1In(3)
      REAL(ReKi),      INTENT( IN    ) :: pos2In(3)
      REAL(ReKi),      INTENT( IN    ) :: FSPtIn(3)
      REAL(ReKi),      INTENT( IN    ) :: k_hatIn(3)
      REAL(ReKi),      INTENT( IN    ) :: x_hatIn(3)
      REAL(ReKi),      INTENT( IN    ) :: y_hatIn(3)
      REAL(ReKi),      INTENT( IN    ) :: n_hatIn(3)
      REAL(ReKi),      INTENT( IN    ) :: Sa1In, Sb1In
      REAL(ReKi),      INTENT( IN    ) :: Sa2In, Sb2In
      REAL(ReKi),      INTENT( IN    ) :: dlIn
      REAL(ReKi),      INTENT(   OUT ) :: F_B1(6)
      REAL(ReKi),      INTENT(   OUT ) :: F_B2(6)
      INTEGER(IntKi),  INTENT(   OUT ) :: ErrStat ! Error status of the operation
      CHARACTER(*),    INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None

      REAL(DbKi)                       :: dFdl1(6), dFdlMid(6), dFdl2(6), F_B(6)
      REAL(DbKi)                       :: i, dl, Sa1, Sa2, SaMid, Sb1, Sb2, SbMid, dSadl, dSbdl, posMid(3), pos1(3), pos2(3), FSPt(3), k_hat(3), x_hat(3), y_hat(3), n_hat(3)
      INTEGER(IntKi)                   :: secStat1, secStatMid, secStat2
      CHARACTER(*),    PARAMETER       :: routineName = "getElementHstLds_Mod2_Rec"
      INTEGER(IntKi)                   :: errStat2
      CHARACTER(ErrMsgLen)             :: errMsg2
      ErrStat   = ErrID_None
      ErrMsg    = ""  
  
      pos1   = REAL(pos1In,DbKi)
      pos2   = REAL(pos2In,DbKi)
      Sa1    = REAL(Sa1In,DbKi)
      Sa2    = REAL(Sa2In,DbKi)
      Sb1    = REAL(Sb1In,DbKi)
      Sb2    = REAL(Sb2In,DbKi)
      dl     = REAL(dlIn,DbKi)
      dSadl  = (Sa2-Sa1)/dl
      dSbdl  = (Sb2-Sb1)/dl
      SaMid  = 0.5*(  Sa1+  Sa2)
      SbMid  = 0.5*(  Sb1+  Sb2)
      posMid = 0.5*(pos1In+pos2In)
      FSPt   = REAL(FSPtIn,DbKi)
      k_hat  = REAL(k_hatIn,DbKi)
      x_hat  = REAL(x_hatIn,DbKi)
      y_hat  = REAL(y_hatIn,DbKi)
      n_hat  = REAL(n_hatIn,DbKi)
      
      ! Avoid sections coincident with the SWL
      IF ( ABS(k_hat(3)) > 0.999999_ReKi ) THEN ! Vertical member
         IF ( EqualRealNos( pos1(3), 0.0_DbKi ) ) THEN
            pos1(3) = pos1(3) - 1.0E-6 * dl
         END IF
         IF ( EqualRealNos( pos2(3), 0.0_DbKi ) ) THEN
            pos2(3) = pos2(3) - 1.0E-6 * dl
         END IF
         IF ( EqualRealNos( posMid(3), 0.0_DbKi ) ) THEN
            posMid(3) = posMid(3) - 1.0E-6 * dl
         END IF
      END IF

      ! Section load at node 1
      CALL GetSectionHstLds_Rec( pos1, pos1,   k_hat, x_hat, y_hat, Sa1,   Sb1,   dSadl, dSbdl, FSPt, n_hat, dFdl1,   secStat1)

      ! Section load at midpoint
      CALL GetSectionHstLds_Rec( pos1, posMid, k_hat, x_hat, y_hat, SaMid, SbMid, dSadl, dSbdl, FSPt, n_hat, dFdlMid, secStatMid)

      ! Section load at node 2
      CALL GetSectionHstLds_Rec( pos1, pos2,   k_hat, x_hat, y_hat, Sa2,   Sb2,   dSadl, dSbdl, FSPt, n_hat, dFdl2,   secStat2)

      ! Adaptively refine the load integration over the element
      CALL RefineElementHstLds_Rec(pos1,pos1,posMid,pos2,FSPt,Sa1,SaMid,Sa2,Sb1,SbMid,Sb2,dl,dSadl,dSbdl, &
              secStat1,secStatMid,secStat2,k_hat,x_hat,y_hat,n_hat,dFdl1,dFdlMid,dFdl2,1,F_B,ErrStat2,ErrMsg2)
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! Distribute the hydrostatic load to the two end nodes
      F_B1(1:3) = 0.5 * F_B(1:3)
      F_B2(1:3) = 0.5 * F_B(1:3)
      F_B(4:6)  = F_B(4:6) - CROSS_PRODUCT(k_hat*dl,F_B2(1:3))
      F_B1(4:6) = 0.5 * F_B(4:6)
      F_B2(4:6) = 0.5 * F_B(4:6)

   END SUBROUTINE getElementHstLds_Mod2_Rec

   RECURSIVE SUBROUTINE RefineElementHstLds_Cyl( origin, pos1, posMid, pos2, FSPt, r1, rMid, r2, dl, dRdl,secStat1,secStatMid,secStat2, k_hat, y_hat, z_hat, n_hat, dFdl1, dFdlMid, dFdl2, recurLvl, F_B_5pt, ErrStat, ErrMsg)

      REAL(DbKi),      INTENT( IN    ) :: origin(3)
      REAL(DbKi),      INTENT( IN    ) :: pos1(3)
      REAL(DbKi),      INTENT( IN    ) :: posMid(3)
      REAL(DbKi),      INTENT( IN    ) :: pos2(3)
      REAL(DbKi),      INTENT( IN    ) :: FSPt(3)
      REAL(DbKi),      INTENT( IN    ) :: r1
      REAL(DbKi),      INTENT( IN    ) :: rMid
      REAL(DbKi),      INTENT( IN    ) :: r2
      REAL(DbKi),      INTENT( IN    ) :: dl
      REAL(DbKi),      INTENT( IN    ) :: dRdl
      INTEGER(IntKi),  INTENT( IN    ) :: secStat1
      INTEGER(IntKi),  INTENT( IN    ) :: secStatMid
      INTEGER(IntKi),  INTENT( IN    ) :: secStat2
      REAL(DbKi),      INTENT( IN    ) :: k_hat(3)
      REAL(DbKi),      INTENT( IN    ) :: y_hat(3)
      REAL(DbKi),      INTENT( IN    ) :: z_hat(3)
      REAL(DbKi),      INTENT( IN    ) :: n_hat(3)
      REAL(DbKi),      INTENT( IN    ) :: dFdl1(6)
      REAL(DbKi),      INTENT( IN    ) :: dFdlMid(6)
      REAL(DbKi),      INTENT( IN    ) :: dFdl2(6)
      INTEGER(IntKi),  INTENT( IN    ) :: recurLvl
      REAL(DbKi),      INTENT(   OUT ) :: F_B_5pt(6)
      INTEGER(IntKi),  INTENT(   OUT ) :: ErrStat ! Error status of the operation
      CHARACTER(*),    INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None

      REAL(DbKi)                       :: theta1,theta2
      REAL(DbKi)                       :: posMidL(3), posMidR(3)
      REAL(DbKi)                       :: rMidL, rMidR
      REAL(DbKi)                       :: dFdlMidL(6), dFdlMidR(6), F_B_3pt(6)
      REAL(DbKi)                       :: error(6), tmp(6)
      LOGICAL                          :: refine, tolMet
      INTEGER(IntKi)                   :: i
      INTEGER(IntKi)                   :: secStatMidL, secStatMidR
      REAL(DbKi),      PARAMETER       :: RelTol      = 1.0E-6
      REAL(DbKi),      PARAMETER       :: AbsTol      = 1.0E-8
      INTEGER(IntKi),  PARAMETER       :: maxRecurLvl = 50
      CHARACTER(*),    PARAMETER       :: RoutineName = "RefineElementHstLds_Cyl"
      
      ErrStat = ErrID_None
      ErrMsg  = ""

      posMidL = 0.5*(pos1+posMid)
      posMidR = 0.5*(posMid+pos2)
      rMidL   = 0.5*(r1+rMid)
      rMidR   = 0.5*(rMid+r2)

      ! Avoid sections coincident with the SWL
      IF ( ABS(k_hat(3)) > 0.999999_ReKi ) THEN ! Vertical member
         IF ( EqualRealNos( posMidL(3), 0.0_DbKi ) ) THEN
            posMidL(3) = posMidL(3) - 1.0E-6 * dl
         END IF
         IF ( EqualRealNos( posMidR(3), 0.0_DbKi ) ) THEN
            posMidR(3) = posMidR(3) - 1.0E-6 * dl
         END IF
      END IF

      ! Total hydrostatic load on the element (Simpsons Rule)
      F_B_3pt = (dFdl1 + 4.0*dFdlMid + dFdl2) * dl/6.0

      ! Mid point of left section
      CALL GetSectionFreeSurfaceIntersects_Cyl( posMidL, FSPt, REAL(k_hat,ReKi), REAL(y_hat,ReKi), REAL(z_hat,ReKi), REAL(n_hat,ReKi), rMidL, theta1, theta2, secStatMidL)
      CALL GetSectionHstLds_Cyl( origin, posMidL, k_hat, y_hat, z_hat, rMidL, dRdl, theta1, theta2, dFdlMidL)

      ! Mid point of right section
      CALL GetSectionFreeSurfaceIntersects_Cyl( posMidR, FSPt, REAL(k_hat,ReKi), REAL(y_hat,ReKi), REAL(z_hat,ReKi), REAL(n_hat,ReKi), rMidR, theta1, theta2, secStatMidR)
      CALL GetSectionHstLds_Cyl( origin, posMidR, k_hat, y_hat, z_hat, rMidR, dRdl, theta1, theta2, dFdlMidR)
      
      F_B_5pt = (dFdl1 + 4.0*dFdlMidL + 2.0*dFdlMid + 4.0*dFdlMidR + dFdl2) * dl/12.0

      error = ABS(F_B_3pt - F_B_5pt)
      tolMet = .TRUE.
      DO i = 1,6
         IF ( error(i) > MAX(RelTol*ABS(F_B_5pt(i)),AbsTol) ) THEN
            tolMet = .FALSE.
         END IF
      END DO
      refine = .NOT. tolMet
      IF (ABS(secStat1-secStat2)>1) THEN ! (Sub)element bounds the waterplane
         refine = .TRUE. ! Keep refining irrespective of tolMet to avoid premature termination
      END IF
      IF ( recurLvl > maxRecurLvl ) THEN
         refine = .FALSE.
         IF (.NOT. tolMet) THEN
            CALL SetErrStat(ErrID_Warn, 'Tolerance for element hydrostatic load not met after the maximum allowed level of recursion is reached. Consider reducing MDivSize.', ErrStat, ErrMsg, RoutineName )
         ! ELSE
            ! Free surface is likely normal to the element.
         END IF
      END IF
      
      IF (refine) THEN ! Recursively refine the load integration if tolerance not met
         CALL RefineElementHstLds_Cyl(origin,pos1,posMidL,posMid,FSPt,r1,rMidL,rMid,0.5*dl,dRdl,secStat1,secStatMidL,secStatMid,k_hat,y_hat,z_hat,n_hat,dFdl1,dFdlMidL,dFdlMid, recurLvl+1, tmp, ErrStat, ErrMsg)
         CALL RefineElementHstLds_Cyl(origin,posMid,posMidR,pos2,FSPt,rMid,rMidR,r2,0.5*dl,dRdl,secStatMid,secStatMidR,secStat2,k_hat,y_hat,z_hat,n_hat,dFdlMid,dFdlMidR,dFdl2, recurLvl+1, F_B_5pt, ErrStat, ErrMsg)
         F_B_5pt = F_B_5pt + tmp
      END IF

   END SUBROUTINE RefineElementHstLds_Cyl

   RECURSIVE SUBROUTINE RefineElementHstLds_Rec( origin, pos1, posMid, pos2, FSPt, Sa1, SaMid, Sa2, Sb1, SbMid, Sb2, dl, dSadl, dSbdl, &
      secStat1, secStatMid, secStat2, k_hat, x_hat, y_hat, n_hat, dFdl1, dFdlMid, dFdl2, recurLvl, F_B_5pt, ErrStat, ErrMsg)

      REAL(DbKi),      INTENT( IN    ) :: origin(3)
      REAL(DbKi),      INTENT( IN    ) :: pos1(3)
      REAL(DbKi),      INTENT( IN    ) :: posMid(3)
      REAL(DbKi),      INTENT( IN    ) :: pos2(3)
      REAL(DbKi),      INTENT( IN    ) :: FSPt(3)
      REAL(DbKi),      INTENT( IN    ) :: Sa1, Sb1
      REAL(DbKi),      INTENT( IN    ) :: SaMid, SbMid
      REAL(DbKi),      INTENT( IN    ) :: Sa2, Sb2
      REAL(DbKi),      INTENT( IN    ) :: dl
      REAL(DbKi),      INTENT( IN    ) :: dSadl, dSbdl
      INTEGER(IntKi),  INTENT( IN    ) :: secStat1
      INTEGER(IntKi),  INTENT( IN    ) :: secStatMid
      INTEGER(IntKi),  INTENT( IN    ) :: secStat2
      REAL(DbKi),      INTENT( IN    ) :: k_hat(3)
      REAL(DbKi),      INTENT( IN    ) :: x_hat(3)
      REAL(DbKi),      INTENT( IN    ) :: y_hat(3)
      REAL(DbKi),      INTENT( IN    ) :: n_hat(3)
      REAL(DbKi),      INTENT( IN    ) :: dFdl1(6)
      REAL(DbKi),      INTENT( IN    ) :: dFdlMid(6)
      REAL(DbKi),      INTENT( IN    ) :: dFdl2(6)
      INTEGER(IntKi),  INTENT( IN    ) :: recurLvl
      REAL(DbKi),      INTENT(   OUT ) :: F_B_5pt(6)
      INTEGER(IntKi),  INTENT(   OUT ) :: ErrStat ! Error status of the operation
      CHARACTER(*),    INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None

      REAL(DbKi)                       :: posMidL(3), posMidR(3)
      REAL(DbKi)                       :: SaMidL, SbMidL, SaMidR, SbMidR
      REAL(DbKi)                       :: dFdlMidL(6), dFdlMidR(6), F_B_3pt(6)
      REAL(DbKi)                       :: error(6), tmp(6)
      LOGICAL                          :: refine, tolMet
      INTEGER(IntKi)                   :: i
      INTEGER(IntKi)                   :: secStatMidL, secStatMidR
      REAL(DbKi),      PARAMETER       :: RelTol      = 1.0E-6
      REAL(DbKi),      PARAMETER       :: AbsTol      = 1.0E-8
      INTEGER(IntKi),  PARAMETER       :: maxRecurLvl = 50
      CHARACTER(*),    PARAMETER       :: RoutineName = "RefineElementHstLds_Rec"
      
      ErrStat = ErrID_None
      ErrMsg  = ""

      posMidL = 0.5*(pos1+posMid)
      posMidR = 0.5*(posMid+pos2)
      SaMidL  = 0.5*(Sa1+SaMid)
      SbMidL  = 0.5*(Sb1+SbMid)
      SaMidR  = 0.5*(SaMid+Sa2)
      SbMidR  = 0.5*(SbMid+Sb2)

      ! Avoid sections coincident with the SWL
      IF ( ABS(k_hat(3)) > 0.999999_ReKi ) THEN ! Vertical member
         IF ( EqualRealNos( posMidL(3), 0.0_DbKi ) ) THEN
            posMidL(3) = posMidL(3) - 1.0E-6 * dl
         END IF
         IF ( EqualRealNos( posMidR(3), 0.0_DbKi ) ) THEN
            posMidR(3) = posMidR(3) - 1.0E-6 * dl
         END IF
      END IF

      ! Total hydrostatic load on the element (Simpsons Rule)
      F_B_3pt = (dFdl1 + 4.0*dFdlMid + dFdl2) * dl/6.0

      ! Mid point of left section
      CALL GetSectionHstLds_Rec( origin, posMidL, k_hat, x_hat, y_hat, SaMidL, SbMidL, dSadl, dSbdl, FSPt, n_hat, dFdlMidL, secStatMidL)

      ! Mid point of right section
      CALL GetSectionHstLds_Rec( origin, posMidR, k_hat, x_hat, y_hat, SaMidR, SbMidR, dSadl, dSbdl, FSPt, n_hat, dFdlMidR, secStatMidR)
      
      F_B_5pt = (dFdl1 + 4.0*dFdlMidL + 2.0*dFdlMid + 4.0*dFdlMidR + dFdl2) * dl/12.0

      error = ABS(F_B_3pt - F_B_5pt)
      tolMet = .TRUE.
      DO i = 1,6
         IF ( error(i) > MAX(RelTol*ABS(F_B_5pt(i)),AbsTol) ) THEN
            tolMet = .FALSE.
         END IF
      END DO
      refine = .NOT. tolMet
      IF (ABS(secStat1-secStat2)>1) THEN ! (Sub)element bounds the waterplane
         refine = .TRUE. ! Keep refining irrespective of tolMet to avoid premature termination
      END IF
      IF ( recurLvl > maxRecurLvl ) THEN
         refine = .FALSE.
         IF (.NOT. tolMet) THEN
            CALL SetErrStat(ErrID_Warn, 'Tolerance for element hydrostatic load not met after the maximum allowed level of recursion is reached. Consider reducing MDivSize.', ErrStat, ErrMsg, RoutineName )
         ! ELSE
            ! Free surface is likely normal to the element.
         END IF
      END IF
      
      IF (refine) THEN ! Recursively refine the load integration if tolerance not met
         CALL RefineElementHstLds_Rec(origin,pos1,posMidL,posMid,FSPt,Sa1,SaMidL,SaMid,Sb1,SbMidL,SbMid,0.5*dl,dSadl,dSbdl, &
                secStat1,secStatMidL,secStatMid,k_hat,x_hat,y_hat,n_hat,dFdl1,dFdlMidL,dFdlMid,recurLvl+1,tmp,ErrStat,ErrMsg)
         CALL RefineElementHstLds_Rec(origin,posMid,posMidR,pos2,FSPt,SaMid,SaMidR,Sa2,SbMid,SbMidR,Sb2,0.5*dl,dSadl,dSbdl, &
                secStatMid,secStatMidR,secStat2,k_hat,x_hat,y_hat,n_hat,dFdlMid,dFdlMidR,dFdl2,recurLvl+1,F_B_5pt,ErrStat,ErrMsg)
         F_B_5pt = F_B_5pt + tmp
      END IF

   END SUBROUTINE RefineElementHstLds_Rec

   SUBROUTINE GetEndPlateHstLds_Cyl(pos0, k_hat, y_hat, z_hat, R, theta1, theta2, F)

      REAL(ReKi),      INTENT( IN    ) :: pos0(3)
      REAL(ReKi),      INTENT( IN    ) :: k_hat(3)
      REAL(ReKi),      INTENT( IN    ) :: y_hat(3)
      REAL(ReKi),      INTENT( IN    ) :: z_hat(3)
      REAL(ReKi),      INTENT( IN    ) :: R
      REAL(DbKi),      INTENT( IN    ) :: theta1
      REAL(DbKi),      INTENT( IN    ) :: theta2
      REAL(ReKi),      INTENT(   OUT ) :: F(6)
      REAL(DbKi)                       :: C0, C1, C2, a, b, tmp1, tmp2, tmp3
      REAL(DbKi)                       :: Z0, cosPhi, dTheta
      REAL(DbKi)                       :: y1, y2
      REAL(DbKi)                       :: z1, z2, z1_2, z2_2, z1_3, z2_3, z1_4, z2_4
      REAL(DbKi)                       :: dy, dy_3, dz, dz_2, dz_3, dz_4, sz
      REAL(DbKi)                       :: R_2, R_4
      REAL(DbKi)                       :: Fk, My, Mz

      Z0     = pos0(3)
      cosPhi = SQRT(k_hat(1)**2+k_hat(2)**2)
      dTheta = theta2-theta1
      y1     = R*COS(theta1)
      z1     = R*SIN(theta1)
      y2     = R*COS(theta2)
      z2     = R*SIN(theta2)
      z1_2   = z1*z1
      z1_3   = z1*z1_2
      z1_4   = z1*z1_3
      z2_2   = z2*z2
      z2_3   = z2*z2_2
      z2_4   = z2*z2_3
      R_2    = R*R
      R_4    = R_2*R_2
      dy     = y2-y1
      sz     = z1+z2
      dy_3   = y2*y2*y2-y1*y1*y1
      dz_2   = z2_2-z1_2
      dz_3   = z2_3-z1_3
      dz_4   = z2_4-z1_4
      tmp1   = y1*z2-y2*z1
      tmp2   = z1_2+z1*z2+z2_2

      ! End plate force in the k_hat direction
      Fk     = -0.5*Z0*(R_2*dTheta-tmp1) + cosPhi/6.0*( 2.0*dy_3 - z1*z2*dy - z1_2*(y2+2.0*y1) + z2_2*(y1+2.0*y2) )
      F(1:3) = p%WaveField%WtrDens * g * Fk * k_hat

      ! End plate moment in the y_hat and z_hat direction
      My     = Z0/6.0*( 2.0*dy_3 + 2.0*dy*tmp2 + 3.0*tmp1*sz ) &    ! y_hat component
               + cosPhi/24.0*( -3.0*R_4*dTheta + 3.0*y1*z1*(2.0*z1_2-R_2) - 3.0*y2*z2*(2.0*z2_2-R_2) &
                                                            + 6.0*dy*sz*(z1_2+z2_2) + 8.0*tmp1*tmp2  )
      IF (EqualRealNos(z1, z2)) THEN ! z_hat component (Nonzero only when z1 /= z2)
         Mz = 0.0
      ELSE
         dz = z2-z1
         a  = dy/dz
         b  = tmp1/dz
         tmp1 = a*a+1.0
         tmp2 = a*b
         tmp3 = b*b-R_2
         Mz =     -Z0/ 6.0*(    tmp1*dz_3 + 3.0*tmp2*dz_2 + 3.0*tmp3*dz  ) &
              -cosPhi/24.0*(3.0*tmp1*dz_4 + 8.0*tmp2*dz_3 + 6.0*tmp3*dz_2)
      END IF
      F(4:6) = p%WaveField%WtrDens * g * (My*y_hat + Mz*z_hat)

   END SUBROUTINE GetEndPlateHstLds_Cyl

   SUBROUTINE GetEndPlateHstLds_Rec(pos0, k_hat, x_hat, y_hat, Sa, Sb, rFS, nFS, F)

      REAL(ReKi),      INTENT( IN    ) :: pos0(3)
      REAL(ReKi),      INTENT( IN    ) :: k_hat(3)
      REAL(ReKi),      INTENT( IN    ) :: x_hat(3)
      REAL(ReKi),      INTENT( IN    ) :: y_hat(3)
      REAL(ReKi),      INTENT( IN    ) :: Sa, Sb
      REAL(ReKi),      INTENT( IN    ) :: rFS(3)
      REAL(ReKi),      INTENT( IN    ) :: nFS(3)
      REAL(ReKi),      INTENT(   OUT ) :: F(6)

      INTEGER(IntKi)                   :: i, numVInWtr
      REAL(DbKi)                       :: z0, s1, s2, h1s1, h1s2, h2s1, h2s2
      REAL(DbKi)                       :: rv(3,4), Ftmp1(6), Ftmp2(6)
      LOGICAL                          :: vInWtr(4)

      z0     = pos0(3)
      ! Global coordinates of the four vertices of the section
      rv(:,1) = pos0 + x_hat * (-0.5*Sa) + y_hat * (-0.5*Sb)
      rv(:,2) = pos0 + x_hat * ( 0.5*Sa) + y_hat * (-0.5*Sb)
      rv(:,3) = pos0 + x_hat * ( 0.5*Sa) + y_hat * ( 0.5*Sb)
      rv(:,4) = pos0 + x_hat * (-0.5*Sa) + y_hat * ( 0.5*Sb)

      ! Check for and count vertices in water
      numVInWtr = 0
      do i = 1,4
         vInWtr(i) = ( Dot_Product(rv(:,i)-rFS,nFS) <= 0.0 )
         if ( vInWtr(i) ) then
            numVInWtr = numVInWtr + 1
         end if
      end do

      if (numVInWtr == 0) then ! Dry endplate
         F = 0.0
      else if (numVInWtr == 1) then ! Only one vertex in water
         if (vInWtr(1)) then
            ! Sides 4 & 1 intersects the free surface
            h2s1 =  Sb * ( 0.5 - dot_product(rFS-rv(:,4),nFS)/dot_product(rv(:,1)-rv(:,4),nFS) )
            h2s2 = -0.5*Sb
            h1s1 = -0.5*Sb
            h1s2 = -0.5*Sb
            s1   = -0.5*Sa
            s2   =  Sa * (-0.5 + dot_product(rFS-rv(:,1),nFS)/dot_product(rv(:,2)-rv(:,1),nFS) )
         else if (vInWtr(2)) then
            ! Sides 1 & 2 intersects the free surface
            h2s1 = -0.5*Sb
            h2s2 =  Sb * (-0.5 + dot_product(rFS-rv(:,2),nFS)/dot_product(rv(:,3)-rv(:,2),nFS) )
            h1s1 = -0.5*Sb
            h1s2 = -0.5*Sb
            s1   =  Sa * (-0.5 + dot_product(rFS-rv(:,1),nFS)/dot_product(rv(:,2)-rv(:,1),nFS) )
            s2   =  0.5*Sa
         else if (vInWtr(3)) then
            ! Sides 2 & 3 intersects the free surface
            h2s1 =  0.5*Sb
            h2s2 =  0.5*Sb
            h1s1 =  0.5*Sb
            h1s2 =  Sb * (-0.5 + dot_product(rFS-rv(:,2),nFS)/dot_product(rv(:,3)-rv(:,2),nFS) )
            s1   =  Sa * ( 0.5 - dot_product(rFS-rv(:,3),nFS)/dot_product(rv(:,4)-rv(:,3),nFS) )
            s2   =  0.5*Sa
         else if (vInWtr(4)) then
            ! Sides 3 & 4 intersects the free surface
            h2s1 =  0.5*Sb
            h2s2 =  0.5*Sb
            h1s1 =  Sb * ( 0.5 - dot_product(rFS-rv(:,4),nFS)/dot_product(rv(:,1)-rv(:,4),nFS) )
            h1s2 =  0.5*Sb
            s1   = -0.5*Sa
            s2   =  Sa * (0.5 - dot_product(rFS-rv(:,3),nFS)/dot_product(rv(:,4)-rv(:,3),nFS) )
         end if
         call GetHstLdsOnTrapezoid(REAL(pos0,DbKi),s1,s2,h1s1,h1s2,h2s1,h2s2,REAL(k_hat,DbKi),REAL(x_hat,DbKi),REAL(y_hat,DbKi),Ftmp1)
         F = Ftmp1
      else if (numVInWtr == 2) then ! Two neighboring vertices in water
         if (vInWtr(1) .and. vInWtr(2)) then
            ! Sides 2 & 4 intersects the free surface
            ! Side 1 submerged and side 3 dry
            h2s1 =  Sb * ( 0.5 - dot_product(rFS-rv(:,4),nFS)/dot_product(rv(:,1)-rv(:,4),nFS) )
            h2s2 =  Sb * (-0.5 + dot_product(rFS-rv(:,2),nFS)/dot_product(rv(:,3)-rv(:,2),nFS) )
            h1s1 = -0.5*Sb
            h1s2 = -0.5*Sb
            s1   = -0.5*Sa
            s2   =  0.5*Sa
            Ftmp2 = 0.0_DbKi
         else if (vInWtr(2) .and. vInWtr(3)) then
            ! Sides 1 & 3 intersects the free surface
            ! Side 2 submerged and side 4 dry
            ! Integrate in two pieces
            h2s1 = -0.5*Sb
            h2s2 =  0.5*Sb
            h1s1 = -0.5*Sb
            h1s2 = -0.5*Sb
            s1   =  Sa * (-0.5 + dot_product(rFS-rv(:,1),nFS)/dot_product(rv(:,2)-rv(:,1),nFS) )
            s2   =  Sa * (0.5 - dot_product(rFS-rv(:,3),nFS)/dot_product(rv(:,4)-rv(:,3),nFS) )
            call GetHstLdsOnTrapezoid(REAL(pos0,DbKi),s2,0.5_DbKi*Sa,-0.5_DbKi*Sb,-0.5_DbKi*Sb,0.5_DbKi*Sb,0.5_DbKi*Sb,REAL(k_hat,DbKi),REAL(x_hat,DbKi),REAL(y_hat,DbKi),Ftmp2)
         else if (vInWtr(3) .and. vInWtr(4)) then
            ! Sides 2 & 4 intersects the free surface
            ! Side 3 submerged and side 1 dry
            h2s1 =  0.5*Sb
            h2s2 =  0.5*Sb
            h1s1 =  Sb * ( 0.5 - dot_product(rFS-rv(:,4),nFS)/dot_product(rv(:,1)-rv(:,4),nFS) )
            h1s2 =  Sb * (-0.5 + dot_product(rFS-rv(:,2),nFS)/dot_product(rv(:,3)-rv(:,2),nFS) )
            s1   = -0.5*Sa
            s2   =  0.5*Sa
            Ftmp2 = 0.0_DbKi
         else if (vInWtr(4) .and. vInWtr(1)) then
            ! Sides 1 & 3 intersects the free surface
            ! Side 4 submerged and side 2 dry
            ! Integrate in two pieces
            h2s1 =  0.5*Sb
            h2s2 = -0.5*Sb
            h1s1 = -0.5*Sb
            h1s2 = -0.5*Sb
            s1   =  Sa * ( 0.5 - dot_product(rFS-rv(:,3),nFS)/dot_product(rv(:,4)-rv(:,3),nFS) )
            s2   =  Sa * (-0.5 + dot_product(rFS-rv(:,1),nFS)/dot_product(rv(:,2)-rv(:,1),nFS) )
            call GetHstLdsOnTrapezoid(REAL(pos0,DbKi),-0.5_DbKi*Sa,s1,-0.5_DbKi*Sb,-0.5_DbKi*Sb,0.5_DbKi*Sb,0.5_DbKi*Sb,REAL(k_hat,DbKi),REAL(x_hat,DbKi),REAL(y_hat,DbKi),Ftmp2)
         end if
         call GetHstLdsOnTrapezoid(REAL(pos0,DbKi),s1,s2,h1s1,h1s2,h2s1,h2s2,REAL(k_hat,DbKi),REAL(x_hat,DbKi),REAL(y_hat,DbKi),Ftmp1)
         F = Ftmp1 + Ftmp2
      else if (numVInWtr == 3) then ! Only one vertex out of water
         if (.not. vInWtr(1)) then
            ! Sides 4 & 1 intersects the free surface
            h2s1 =  Sb * ( 0.5 - dot_product(rFS-rv(:,4),nFS)/dot_product(rv(:,1)-rv(:,4),nFS) )
            h2s2 = -0.5*Sb
            h1s1 = -0.5*Sb
            h1s2 = -0.5*Sb
            s1   = -0.5*Sa;
            s2   =  Sa * (-0.5 + dot_product(rFS-rv(:,1),nFS)/dot_product(rv(:,2)-rv(:,1),nFS) )
         else if (.not. vInWtr(2)) then
            ! Sides 1 & 2 intersects the free surface
            h2s1 = -0.5*Sb
            h2s2 =  Sb * (-0.5 + dot_product(rFS-rv(:,2),nFS)/dot_product(rv(:,3)-rv(:,2),nFS) )
            h1s1 = -0.5*Sb
            h1s2 = -0.5*Sb
            s1   =  Sa * (-0.5 + dot_product(rFS-rv(:,1),nFS)/dot_product(rv(:,2)-rv(:,1),nFS) )
            s2   =  0.5*Sa
         else if (.not. vInWtr(3)) then
            ! Sides 2 & 3 intersects the free surface
            h2s1 =  0.5*Sb
            h2s2 =  0.5*Sb
            h1s1 =  0.5*Sb
            h1s2 =  Sb * (-0.5 + dot_product(rFS-rv(:,2),nFS)/dot_product(rv(:,3)-rv(:,2),nFS) )
            s1   =  Sa * ( 0.5 - dot_product(rFS-rv(:,3),nFS)/dot_product(rv(:,4)-rv(:,3),nFS) )
            s2   =  0.5*Sa
         else if (.not. vInWtr(4)) then
            ! Sides 3 & 4 intersects the free surface
            h2s1 =  0.5*Sb
            h2s2 =  0.5*Sb
            h1s1 =  Sb * ( 0.5 - dot_product(rFS-rv(:,4),nFS)/dot_product(rv(:,1)-rv(:,4),nFS) )
            h1s2 =  0.5*Sb
            s1   = -0.5*Sa
            s2   =  Sa * ( 0.5 - dot_product(rFS-rv(:,3),nFS)/dot_product(rv(:,4)-rv(:,3),nFS) )
         end if
         call GetHstLdsOnTrapezoid(REAL(pos0,DbKi),s1,s2,h1s1,h1s2,h2s1,h2s2,REAL(k_hat,DbKi),REAL(x_hat,DbKi),REAL(y_hat,DbKi),Ftmp1)
         F(1:3) = -p%WaveField%WtrDens*g*z0*Sa*Sb*k_hat - Ftmp1(1:3)
         F(4:6) =  p%WaveField%WtrDens*g*(Sa**3*Sb*x_hat(3)*y_hat-Sa*Sb**3*y_hat(3)*x_hat)/12.0 - Ftmp1(4:6)
      else if (numVInWtr == 4) then ! Submerged endplate
         F(1:3) = -p%WaveField%WtrDens*g*z0*Sa*Sb*k_hat
         F(4:6) =  p%WaveField%WtrDens*g*(Sa**3*Sb*x_hat(3)*y_hat-Sa*Sb**3*y_hat(3)*x_hat)/12.0
      end if

   END SUBROUTINE GetEndPlateHstLds_Rec

   SUBROUTINE GetHstLdsOnTrapezoid(pos0, s1, s2, h1s1, h1s2, h2s1, h2s2, k_hat, x_hat, y_hat, F)
      REAL(DbKi),      INTENT( IN    ) :: pos0(3)
      REAL(DbKi),      INTENT( IN    ) :: s1, s2, h1s1, h1s2, h2s1, h2s2
      REAL(DbKi),      INTENT( IN    ) :: k_hat(3), x_hat(3), y_hat(3)
      REAL(DbKi),      INTENT(   OUT ) :: F(6)

      REAL(DbKi)                       :: z0
      REAL(DbKi)                       :: ds, ds2, ds3, ds4
      REAL(DbKi)                       :: p1, q1, p2, q2
      REAL(DbKi)                       :: dp, dp2, dp3
      REAL(DbKi)                       :: dq, dq2, dq3
      REAL(DbKi)                       :: dpq, dp2q, dpq2, tmp

      if (EqualRealNos(s1,s2)) then
         F = 0.0
         return;
      end if

      z0   = pos0(3)
      ds   = s2   -s1
      ds2  = s2**2-s1**2
      ds3  = s2**3-s1**3
      ds4  = s2**4-s1**4
      p1   = (h1s2-h1s1)/ds
      q1   = (h1s1*s2-h1s2*s1)/ds
      p2   = (h2s2-h2s1)/ds
      q2   = (h2s1*s2-h2s2*s1)/ds
      dp   = p2-p1
      dq   = q2-q1
      dp2  = p2**2-p1**2
      dq2  = q2**2-q1**2
      dp3  = p2**3-p1**3
      dq3  = q2**3-q1**3
      dpq  = p2*q2-p1*q1
      dp2q = p2**2*q2-p1**2*q1
      dpq2 = p2*q2**2-p1*q1**2
      tmp  = 3.0*dp2*ds4+8.0*dpq*ds3+6.0*dq2*ds2

      F(1:3) = -( 0.5*z0*(dp*ds2+2.0*dq*ds) &
                 +x_hat(3)/6.0*(2.0*dp*ds3+3.0*dq*ds2) &
                 +y_hat(3)/6.0*(dp2*ds3+3.0*dpq*ds2+3.0*(q2**2-q1**2)*ds) &
                ) * k_hat

      F(4:6) = -( z0/6.0*(dp2*ds3+3.0*dpq*ds2+3.0*dq2*ds) &
                 +x_hat(3)/24.0*tmp &
                 +y_hat(3)/12.0*(dp3*ds4+4.0*dp2q*ds3+6.0*dpq2*ds2+4.0*dq3*ds) )*x_hat &
               +( z0/6.0*(2.0*dp*ds3+3.0*dq*ds2) &
                 +x_hat(3)/12.0*(3.0*dp*ds4+4.0*dq*ds3) &
                 +y_hat(3)/24.0*tmp )*y_hat

      F = p%WaveField%WtrDens * g * F

   END SUBROUTINE GetHstLdsOnTrapezoid

   SUBROUTINE getElementHstLds_Mod1( mem, Time, pos1, pos2, Zeta1, Zeta2, k_hat, r1, r2, dl, alphaIn, Is1stElement, F_B0, F_B1, F_B2, ErrStat, ErrMsg )
      
      TYPE(Morison_MemberType), intent(in) :: mem
      REAL(DbKi),      INTENT( IN    ) :: Time
      REAL(ReKi),      INTENT( IN    ) :: pos1(3)
      REAL(ReKi),      INTENT( IN    ) :: pos2(3)
      REAL(ReKi),      INTENT( IN    ) :: Zeta1
      REAL(ReKi),      INTENT( IN    ) :: Zeta2
      REAL(ReKi),      INTENT( IN    ) :: k_hat(3)
      REAL(ReKi),      INTENT( IN    ) :: r1
      REAL(ReKi),      INTENT( IN    ) :: r2
      REAL(ReKi),      INTENT( IN    ) :: dl
      REAL(ReKi),      INTENT( IN    ) :: alphaIn
      LOGICAL,         INTENT( IN    ) :: Is1stElement
      REAL(ReKi),      INTENT(   OUT ) :: F_B0(6) ! Lumped load at the first  node of the last    element
      REAL(ReKi),      INTENT(   OUT ) :: F_B1(6) ! Lumped load at the first  node of the current element
      REAL(ReKi),      INTENT(   OUT ) :: F_B2(6) ! Lumped load at the second node of the current element
      INTEGER(IntKi),  INTENT(   OUT ) :: ErrStat ! Error status of the operation
      CHARACTER(*),    INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None
      REAL(ReKi)                       :: alpha, dRdl, SubRatio
      REAL(ReKi)                       :: Vs, Vrc, Vhc
      REAL(ReKi)                       :: h0, rh, a0, b0, a0b0, s0, C_1, C_2, Z0
      REAL(ReKi)                       :: sinGamma, cosGamma, tanGamma
      REAL(ReKi)                       :: FbVec(3), MbVec(3), FSInt(3), n_hat(3), t_hat(3), s_hat(3), r_hat(3)
      INTEGER(IntKi)                   :: errStat2
      CHARACTER(ErrMsgLen)             :: errMsg2
      INTEGER(IntKi),  PARAMETER       :: pwr = 3 ! Exponent for buoyancy node distribution smoothing
      CHARACTER(*),    PARAMETER       :: RoutineName = "getElementHstLds_Mod1"
      
      F_B0 = 0.0_ReKi
      F_B1 = 0.0_ReKi
      F_B2 = 0.0_ReKi

      dRdl = (r2 - r1)/dl

      IF ( (z1 < Zeta1) .AND. (z2 < Zeta2) ) THEN  ! If element is fully submerged
         ! Compute the waterplane shape, the submerged volume, and it's geometric center
         ! No need to consider tapered and non-tapered elements separately
         Vs   =  Pi*dl   *(r1**2 +     r1*r2 +     r2**2 ) /  3.0_ReKi   ! volume of total submerged portion
         Vhc  =  Pi*dl**2*(r1**2 + 2.0*r1*r2 + 3.0*r2**2 ) / 12.0_ReKi   ! Submerged Volume * h_c

         ! Hydrostatic force on element
         FbVec = (/0.0_ReKi,0.0_ReKi,Vs/) - Pi*( r2*r2*z2 - r1*r1*z1) *k_hat
         FbVec = p%WaveField%WtrDens * g * FbVec

         ! Hydrostatic moment on element about the lower node
         MbVec = (Vhc+0.25*Pi*(r2**4-r1**4)) * Cross_Product(k_hat,(/0.0_ReKi,0.0_ReKi,1.0_ReKi/))
         MbVec = p%WaveField%WtrDens * g * MbVec

         ! Distribute element load to nodes
         alpha = alphaIn*(z2-Zeta2)**pwr/(alphaIn*(z2-Zeta2)**pwr+(1.0_ReKi-alphaIn)*(z1-Zeta1)**pwr)

         ! Hydrostatic force
         F_B1(1:3) = (1-alpha) * FbVec
         F_B2(1:3) =    alpha  * FbVec
         
         ! Hydrostatic moment correction followed by redistribution
         MbVec = MbVec - Cross_Product( k_hat*dl, F_B2(1:3))
         F_B1(4:6) = (1-alpha) * MbVec
         F_B2(4:6) =    alpha  * MbVec
         
      ELSE IF ( (z1 < Zeta1) .AND. (z2 >= Zeta2) ) THEN ! Element is partially submerged
         ! Submergence ratio
         SubRatio = ( Zeta1-pos1(3) ) / ( (Zeta1-pos1(3)) - (Zeta2-pos2(3)) )
         ! The position of the intersection between the free surface and the element
         FSInt    = SubRatio * (pos2-pos1) + pos1
         ! Distances along element centerline from point 1 to the waterplane
         h0       = SubRatio * dl
         ! Scaled radius of element at point where its centerline crosses the waterplane
         rh       = r1 + h0*dRdl
         ! Estimate the free-surface normal at the free-surface intersection, n_hat
         IF ( p%WaveField%WaveStMod > 0_IntKi ) THEN ! If wave stretching is enabled, compute free surface normal
            CALL GetFreeSurfaceNormal( Time, FSInt, rh, n_hat, ErrStat2, ErrMsg2 )
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         ELSE ! Without wave stretching, use the normal of the SWL
            n_hat = (/0.0_ReKi,0.0_ReKi,1.0_ReKi/)
         END IF
         ! Get other relevant unit vectors, t_hat, r_hat, and s_hat
         t_hat    = Cross_Product(k_hat,n_hat)
         sinGamma = SQRT(Dot_Product(t_hat,t_hat))
         cosGamma = Dot_Product(k_hat,n_hat)
         tanGamma = sinGamma/cosGamma
         IF (sinGamma < 0.0001) THEN  ! Free surface normal is aligned with the element
            ! Arbitrary choice for t_hat as long as it is perpendicular to k_hat
            IF ( k_hat(3) < 0.999999_ReKi ) THEN
               t_hat = (/-k_hat(2),k_hat(1),0.0_ReKi/)   
               t_hat = t_hat / SQRT(Dot_Product(t_hat,t_hat))
            ELSE ! k_hat is close to vertical (0,0,1)
               t_hat = (/1.0_ReKi,0.0_ReKi,0.0_ReKi/);
            END IF
         ELSE
            t_hat = t_hat / sinGamma
         END IF
         s_hat = Cross_Product(t_hat,n_hat)
         r_hat = Cross_Product(t_hat,k_hat)

         ! Compute the waterplane shape, the submerged volume, and it's geometric center
         IF (abs(dRdl) < 0.0001) THEN  ! non-tapered member

            IF (cosGamma < 0.0001) THEN
               CALL SetErrStat(ErrID_Fatal, 'Element cannot be parallel to the free surface.  This has happened for Member ID '//trim(num2lstr(mem%MemberID)), errStat, errMsg, RoutineName )
            END IF

            a0   = r1/cosGamma                                       ! Semi major axis of waterplane
            b0   = r1                                                ! Semi minor axis of waterplane
            a0b0 = a0*b0
            s0   = 0.0_ReKi                                          ! Distance from the center of the waterplane to the element centerline
            Vs   =       Pi*r1**2*h0                                 ! volume of total submerged portion
            Vrc  = -0.25*Pi*r1**4*tanGamma                           ! Submerged Volume * r_c
            Vhc  = 0.125*Pi*r1**2* (4.0*h0**2 + r1**2 * tanGamma**2) ! Submerged Volume * h_c

         ELSE  ! tapered member
            C_1  = 1.0_ReKi - dRdl**2 * tanGamma**2
            IF (C_1 < 0.0001) THEN ! The free surface is nearly tangent to the element wall
               CALL SetErrStat(ErrID_Fatal, 'Element cannot be parallel to the free surface.  This has happened for Member ID '//trim(num2lstr(mem%MemberID)), errStat, errMsg, RoutineName )
            END IF

            a0   = rh/(C_1*cosGamma)                                 ! Semi major axis of waterplane
            b0   = rh/sqrt(C_1)                                      ! Semi minor axis of waterplane
            a0b0 = a0*b0
            C_2  = a0b0*rh*cosGamma - r1**3
            s0   = -rh*dRdl*tanGamma/C_1/cosGamma                    ! Distance from the center of the waterplane to the element centerline
            Vs   =  Pi*C_2/(3.0*dRdl)                                ! volume of total submerged portion
            Vrc  = -0.25*Pi *  a0b0*rh**2*sinGamma/C_1                                                 ! Submerged Volume * r_c
            Vhc  =  0.25*Pi * (a0b0*rh**2*cosGamma/C_1 - r1**4 - 4.0_ReKi/3.0_ReKi*r1*C_2 ) /dRdl**2   ! Submerged Volume * h_c

         END IF

         ! z-coordinate of the center of the waterplane in the global earth-fixed system
         Z0 = z1+h0*k_hat(3)+s0*s_hat(3)  

         ! Hydrostatic force on element
         FbVec = (/0.0_ReKi,0.0_ReKi,Vs/) - Pi*a0b0*Z0*n_hat + Pi*r1**2*z1*k_hat
         FbVec = p%WaveField%WtrDens * g * FbVec

         ! Hydrostatic moment on element about the lower node
         MbVec = Cross_Product( Vrc*r_hat+Vhc*k_hat, (/0.0_ReKi,0.0_ReKi,1.0_ReKi/) ) &
                 + 0.25*Pi*a0b0* ( ( s_hat(3)*a0*a0 + 4.0*(s0-h0*sinGamma)*Z0 )*t_hat - t_hat(3)*b0*b0*s_hat ) &
                 - 0.25*Pi*r1**4*(   r_hat(3)                                  *t_hat - t_hat(3)   *   r_hat )
         MbVec = p%WaveField%WtrDens * g * MbVec

         IF ( Is1stElement ) THEN ! This is the 1st element of the member
            ! Assign the element load to the lower (1st) node of the member
            F_B1(1:3) = FbVec
            F_B1(4:6) = MbVec
         ELSE ! This is not the 1st element of the member
            ! Distribute element load to nodes
            alpha = (1.0-alphaIn)*(z1-Zeta1)**pwr / ( -alphaIn*(z2-Zeta2)**pwr + (1.0-alphaIn)*(z1-Zeta1)**pwr )
            ! Hydrostatic force
            F_B0(1:3) = (1-alpha) * FbVec
            F_B1(1:3) =    alpha  * FbVec
            ! Hydrostatic moment correction followed by redistribution
            MbVec = MbVec - Cross_Product( -k_hat*dl, F_B0(1:3))
            F_B0(4:6) = (1-alpha) * MbVec
            F_B1(4:6) =    alpha  * MbVec
         END IF
      END IF
   END SUBROUTINE getElementHstLds_Mod1

   SUBROUTINE YawMember(member, PtfmRefY, ErrStat, ErrMsg)
      Type(Morison_MemberType), intent(inout) :: member
      Real(ReKi),               intent(in   ) :: PtfmRefY
      Integer(IntKi),           intent(  out) :: ErrStat
      Character(*),             intent(  out) :: ErrMsg

      Real(ReKi)                              :: k(3), x_hat(3), y_hat(3)
      Real(ReKi)                              :: kkt(3,3)
      Real(ReKi)                              :: Ak(3,3)
      Integer(IntKi)                          :: ErrStat2
      Character(ErrMsgLen)                    :: ErrMsg2

      Character(*), parameter                 :: RoutineName = 'YawMember'

      ErrStat = ErrID_None
      ErrMsg  = ''

      call hiFrameTransform(h2i,PtfmRefY,member%k,k,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      member%k   = k

      call hiFrameTransform(h2i,PtfmRefY,member%kkt,kkt,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      member%kkt = kkt

      call hiFrameTransform(h2i,PtfmRefY,member%Ak,Ak,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      member%Ak  = Ak

      IF (member%MSecGeom == MSecGeom_Rec) THEN

         call hiFrameTransform(h2i,PtfmRefY,member%x_hat,x_hat,ErrStat2,ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         member%x_hat   = x_hat

         call hiFrameTransform(h2i,PtfmRefY,member%y_hat,y_hat,ErrStat2,ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         member%y_hat   = y_hat

      END IF

   END SUBROUTINE YawMember

   SUBROUTINE YawJoint(JointNo,PtfmRefY,AM_End,An_End,DP_Const_End,I_MG_End,ErrStat,ErrMsg)
      Integer(IntKi),           intent(in   ) :: JointNo
      Real(ReKi),               intent(in   ) :: PtfmRefY
      Real(ReKi),               intent(  out) :: AM_End(3,3)
      Real(ReKi),               intent(  out) :: An_End(3)
      Real(ReKi),               intent(  out) :: DP_Const_End(3)
      Real(ReKi),               intent(  out) :: I_MG_End(3,3)
      Integer(IntKi),           intent(  out) :: ErrStat
      Character(*),             intent(  out) :: ErrMsg
      
      Integer(IntKi)                          :: ErrStat2
      Character(ErrMsgLen)                    :: ErrMsg2

      Character(*), parameter                 :: RoutineName = 'YawJoint'
      
      ErrStat = ErrID_None
      ErrMsg  = ''

      call hiFrameTransform(h2i,PtfmRefY,p%AM_End(:,:,jointNo),AM_End,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call hiFrameTransform(h2i,PtfmRefY,p%An_End(:,jointNo),An_End,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call hiFrameTransform(h2i,PtfmRefY,p%DP_Const_End(:,jointNo),DP_Const_End,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call hiFrameTransform(h2i,PtfmRefY,p%I_MG_End(:,:,jointNo),I_MG_End,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   END SUBROUTINE YawJoint

   SUBROUTINE getMemBallastHiPt(member,z_hi, ErrStat, ErrMsg)
      ! This subroutine returns the highest point of a member's internal ballast
      Type(Morison_MemberType), intent(in   ) :: member
      Real(ReKi),               intent(  out) :: z_hi
      Integer(IntKi),           intent(  out) :: ErrStat
      Character(*),             intent(  out) :: ErrMsg

      Integer(IntKi)                          :: N, elemNo
      Real(ReKi)                              :: CMatrix0(3,3), CMatrix(3,3)
      Real(ReKi)                              :: k_hat(3), x_hat(3), y_hat(3), z_hat(3)
      Real(ReKi)                              :: l, rIn, SaIn, SbIn, z0, pos1(3), pos2(3)

      Character(*), Parameter                 :: RoutineName = 'getMemBallastHiPt'

      ErrStat = ErrId_None
      ErrMsg  = ""

      IF (member%memfloodstatus == 0) THEN  ! Unflooded member
         CALL SetErrStat(ErrId_Fatal," Coding error: getMemBallastHiPt should never be called with an unflooded member with memfloodstatus = 0. ",ErrStat,ErrMsg,RoutineName )
         RETURN
      END IF

      N     = member%NElements

      IF (member%MSecGeom == MSecGeom_Cyl) THEN
         pos1  = m%DispNodePosHst(:,member%NodeIndx(1))
         pos2  = m%DispNodePosHst(:,member%NodeIndx(2))
         k_hat = pos2-pos1
         k_hat = k_hat / SQRT(Dot_Product(k_hat,k_hat))
         CALL GetSectionUnitVectors_Cyl( k_hat, y_hat, z_hat )
         ! Check the starting section
         rIn   = member%Rin(1)
         z_hi  = pos1(3) + rIn * z_hat(3)
         IF (member%memfloodstatus == 1) THEN            ! Fully flooded, check other end
            pos1  = m%DispNodePosHst(:,member%NodeIndx(N  ))
            pos2  = m%DispNodePosHst(:,member%NodeIndx(N+1))
            k_hat = pos2-pos1
            k_hat = k_hat / SQRT(Dot_Product(k_hat,k_hat))
            CALL GetSectionUnitVectors_Cyl( k_hat, y_hat, z_hat )
            ! Check the ending section
            rIn    = member%Rin(N+1)
            z_hi   = MAX( pos2(3) + rIn * z_hat(3), z_hi)
         ELSE IF (member%memfloodstatus == 2) THEN       ! Partially flooded, check the end of the flooded section
            elemNo = member%elem_fill
            pos1  = m%DispNodePosHst(:,member%NodeIndx(elemNo  ))
            pos2  = m%DispNodePosHst(:,member%NodeIndx(elemNo+1))
            k_hat = pos2-pos1
            k_hat = k_hat / SQRT(Dot_Product(k_hat,k_hat))
            CALL GetSectionUnitVectors_Cyl( k_hat, y_hat, z_hat )
            pos2   = pos1 + member%h_fill * k_hat        ! End of flooded section
            l      = member%h_fill/member%dl
            rIn    = member%Rin(elemNo) * (1.0-l) + member%Rin(elemNo+1) * l
            z_hi   = MAX( pos2(3) + rIn * z_hat(3), z_hi)
         END IF
      ELSE IF (member%MSecGeom == MSecGeom_Rec) THEN
         ! DirCos matrix of undisplaced member
         CALL Morison_DirCosMtrx( u%Mesh%Position(:,member%NodeIndx(1  )), u%Mesh%Position(:,member%NodeIndx(N+1)), member%MSpinOrient, CMatrix0 )
         ! Check the vertices of the starting section
         pos1  = m%DispNodePosHst(:,member%NodeIndx(1))
         pos2  = m%DispNodePosHst(:,member%NodeIndx(2))
         z0 = pos1(3)
         CMatrix = matmul(transpose(u%Mesh%Orientation(:,:,member%NodeIndx(1  ))),CMatrix0)
         CALL GetSectionUnitVectors_Rec( CMatrix, x_hat, y_hat )
         SaIn = member%SaIn(1)
         SbIn = member%SbIn(1)
         z_hi =     z0 - 0.5*SaIn*x_hat(3) - 0.5*SbIn*y_hat(3)
         z_hi = MAX(z0 + 0.5*SaIn*x_hat(3) - 0.5*SbIn*y_hat(3), z_hi)
         z_hi = MAX(z0 + 0.5*SaIn*x_hat(3) + 0.5*SbIn*y_hat(3), z_hi)
         z_hi = MAX(z0 - 0.5*SaIn*x_hat(3) + 0.5*SbIn*y_hat(3), z_hi)
         IF (member%memfloodstatus == 1) THEN            ! Fully flooded, check other end
            ! Check the vertices of the ending section
            pos1  = m%DispNodePosHst(:,member%NodeIndx(N  ))
            pos2  = m%DispNodePosHst(:,member%NodeIndx(N+1))
            z0 = pos2(3)
            CMatrix = matmul(transpose(u%Mesh%Orientation(:,:,member%NodeIndx(N+1))),CMatrix0)
            CALL GetSectionUnitVectors_Rec( CMatrix, x_hat, y_hat )
            SaIn = member%SaIn(N+1)
            SbIn = member%SbIn(N+1)
            z_hi = MAX(z0 - 0.5*SaIn*x_hat(3) - 0.5*SbIn*y_hat(3), z_hi)
            z_hi = MAX(z0 + 0.5*SaIn*x_hat(3) - 0.5*SbIn*y_hat(3), z_hi)
            z_hi = MAX(z0 + 0.5*SaIn*x_hat(3) + 0.5*SbIn*y_hat(3), z_hi)
            z_hi = MAX(z0 - 0.5*SaIn*x_hat(3) + 0.5*SbIn*y_hat(3), z_hi)
         ELSE IF (member%memfloodstatus == 2) THEN       ! Partially flooded, check the end of the flooded section
            elemNo = member%elem_fill
            pos1  = m%DispNodePosHst(:,member%NodeIndx(elemNo  ))
            pos2  = m%DispNodePosHst(:,member%NodeIndx(elemNo+1))
            if ( member%h_fill>0.5*member%dl ) then
               CMatrix = matmul(transpose(u%Mesh%Orientation(:,:,member%NodeIndx(elemNo+1))),CMatrix0)
            else
               CMatrix = matmul(transpose(u%Mesh%Orientation(:,:,member%NodeIndx(elemNo  ))),CMatrix0)
            end if
            CALL GetSectionUnitVectors_Rec( CMatrix, x_hat, y_hat )
            k_hat = Cross_Product(x_hat,y_hat)
            pos2 = pos1 + member%h_fill * k_hat          ! End of filled section
            z0   = pos2(3)
            l = member%h_fill/member%dl
            SaIn = member%SaIn(elemNo) * (1.0-l) + member%SaIn(elemNo+1) * l
            SbIn = member%SbIn(elemNo) * (1.0-l) + member%SbIn(elemNo+1) * l
            z_hi = MAX(z0 - 0.5*SaIn*x_hat(3) - 0.5*SbIn*y_hat(3), z_hi)
            z_hi = MAX(z0 + 0.5*SaIn*x_hat(3) - 0.5*SbIn*y_hat(3), z_hi)
            z_hi = MAX(z0 + 0.5*SaIn*x_hat(3) + 0.5*SbIn*y_hat(3), z_hi)
            z_hi = MAX(z0 - 0.5*SaIn*x_hat(3) + 0.5*SbIn*y_hat(3), z_hi)
         END IF
      END IF
   END SUBROUTINE getMemBallastHiPt

   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
      !if (Failed) then
      !   call FailCleanup()
      !endif
   end function Failed

END SUBROUTINE Morison_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
subroutine LumpDistrHydroLoads( f_hydro, k_hat, dl, h_c, lumpedLoad )
   real(ReKi), intent(in   ) :: f_hydro(3)
   real(ReKi), intent(in   ) :: k_hat(3)
   real(ReKi), intent(in   ) :: dl
   real(ReKi), intent(in   ) :: h_c
   real(ReKi), intent(inout) :: lumpedLoad(6)
   lumpedLoad(1:3) = f_hydro*dl
   lumpedLoad(4:6) = cross_product(k_hat*h_c, f_hydro)*dl
end subroutine LumpDistrHydroLoads
!----------------------------------------------------------------------------------------------------------------------------------
! Takes loads on node i in element tilted frame and converts to 6DOF loads at node i and adjacent node
PURE SUBROUTINE DistributeElementLoads(Fl, Fr, M, sinPhi, cosPhi, SinBeta, cosBeta, alpha, F1, F2)
   
   REAL(ReKi),                     INTENT    ( IN    )  :: Fl        ! (N)   axial load about node i
   REAL(ReKi),                     INTENT    ( IN    )  :: Fr        ! (N)   radial load about node i in direction of tilt
   REAL(ReKi),                     INTENT    ( IN    )  :: M         ! (N-m) radial moment about node i, positive in direction of tilt angle
   REAL(ReKi),                     INTENT    ( IN    )  :: sinPhi    ! trig functions of  tilt angle 
   REAL(ReKi),                     INTENT    ( IN    )  :: cosPhi   
   REAL(ReKi),                     INTENT    ( IN    )  :: sinBeta   ! trig functions of heading of tilt
   REAL(ReKi),                     INTENT    ( IN    )  :: cosBeta  
   REAL(ReKi),                     INTENT    ( IN    )  :: alpha     ! fraction of load staying with node i (1-alpha goes to other node)  
   
   REAL(ReKi),                     INTENT    ( OUT   )  :: F1(6)   ! (N, Nm) force/moment vector for node i
   REAL(ReKi),                     INTENT    ( OUT   )  :: F2(6)   ! (N, Nm) force/moment vector for the other node (whether i+1, or i-1)
   REAL(ReKi)                                           :: F(6)
   
   F(1) =  cosBeta*(Fl*sinPhi + Fr*cosPhi)
   F(2) =  sinBeta*(Fl*sinPhi + Fr*cosPhi)
   F(3) =          (Fl*cosPhi - Fr*sinPhi)
   F(4) = -sinBeta * M                    
   F(5) =  cosBeta * M                    
   F(6) =  0.0

   F1 = F*alpha
   F2 = F*(1.0_ReKi-alpha)
   
END SUBROUTINE DistributeElementLoads
!----------------------------------------------------------------------------------------------------------------------------------
! Takes loads on end node i and converts to 6DOF loads, adding to the nodes existing loads
PURE SUBROUTINE AddEndLoad(Fl, M, sinPhi, cosPhi, SinBeta, cosBeta, Fi)
   
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
   Fi(4) = Fi(4) - M*sinBeta
   Fi(5) = Fi(5) + M*cosBeta
   
END SUBROUTINE AddEndLoad


!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for updating discrete states
SUBROUTINE Morison_UpdateDiscState( Time, u, p, x, xd, z, OtherState, m, errStat, errMsg )   
!..................................................................................................................................
   
   REAL(DbKi),                        INTENT(IN   )  :: Time        !< Current simulation time in seconds   
   TYPE(Morison_InputType),           INTENT(IN   )  :: u           !< Inputs at Time                       
   TYPE(Morison_ParameterType),       INTENT(IN   )  :: p           !< Parameters                                 
   TYPE(Morison_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(Morison_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Input:  Discrete states at Time; 
                                                                    !< Output: Discrete states at Time + Interval
   TYPE(Morison_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
   TYPE(Morison_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time          
   TYPE(Morison_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables            
   INTEGER(IntKi),                    INTENT(  OUT)  :: errStat     !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: errMsg      !< Error message if errStat /= ErrID_None
   INTEGER(IntKi)                                    :: J
   INTEGER(IntKi)                                    :: nodeInWater
   REAL(ReKi)                                        :: pos(3), vrel(3), FV(3), vmag, vmagf, An_End(3)
   REAL(SiKi)                                        :: FVTmp(3)
   INTEGER(IntKi)                                    :: errStat2
   CHARACTER(ErrMsgLen)                              :: errMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'Morison_UpdateDiscState'
   
   ! Initialize errStat  
   errStat = ErrID_None         
   errMsg  = ""               

   ! Update state of the relative normal velocity high-pass filter at each joint
   DO J = 1, p%NJoints

      ! Get joint position
      IF (p%WaveDisp == 0 ) THEN
         ! use the initial X,Y location
         pos(1) = u%Mesh%Position(1,J)
         pos(2) = u%Mesh%Position(2,J)
      ELSE
         ! Use current X,Y location
         pos(1) = u%Mesh%TranslationDisp(1,J) + u%Mesh%Position(1,J)
         pos(2) = u%Mesh%TranslationDisp(2,J) + u%Mesh%Position(2,J)
      END IF
      IF (p%WaveField%WaveStMod > 0 .AND. p%WaveDisp /= 0) THEN ! Wave stretching enabled
         pos(3) = u%Mesh%Position(3,J) + u%Mesh%TranslationDisp(3,J) - p%WaveField%MSL2SWL  ! Use the current Z location.
      ELSE ! Wave stretching disabled
         pos(3) = u%Mesh%Position(3,J) - p%WaveField%MSL2SWL  ! We are intentionally using the undisplaced Z position of the node.
      END IF

      ! Get fluid velocity at the joint
      CALL WaveField_GetNodeWaveVel( p%WaveField, m%WaveField_m, Time, pos, .FALSE., .TRUE., nodeInWater, FVTmp, ErrStat2, ErrMsg2 )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      FV   = REAL(FVTmp, ReKi)
      vrel = ( FV - u%Mesh%TranslationVel(:,J) ) * nodeInWater

      ! Transform An_End based on reference yaw offset
      call hiFrameTransform(h2i,u%PtfmRefY,p%An_End(:,j),An_End,ErrStat2,ErrMsg2)
      call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

      ! Compute the dot product of the relative velocity vector with the directional Area of the Joint
      vmag  = vrel(1)*An_End(1) + vrel(2)*An_End(2) + vrel(3)*An_End(3)
      ! High-pass filtering
      vmagf = p%VRelNFiltConst(J) * (vmag + xd%V_rel_n_FiltStat(J))
      ! Update relative normal velocity filter state for joint J 
      xd%V_rel_n_FiltStat(J) = vmagf-vmag

   END DO ! J = 1, p%NJoints

END SUBROUTINE Morison_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE Morison
!**********************************************************************************************************************************
