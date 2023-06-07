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
   use SeaState_Interp
  ! USE HydroDyn_Output_Types
   USE NWTC_Library

   
   IMPLICIT NONE
   
   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: Morison_ProgDesc = ProgDesc( 'Morison', '', '' )

   
      ! ..... Public Subroutines ...................................................................................................
   PUBLIC:: Morison_GenerateSimulationNodes
   
   PUBLIC :: Morison_Init                           ! Initialization routine
   PUBLIC :: Morison_CalcOutput                     ! Routine for computing outputs
   PUBLIC :: Morison_UpdateDiscState                ! Routine for updating discrete states
   
CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
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
      
      DirCos(1, 1) = (z1-z0)/xz
      DirCos(1, 2) = -(x1-x0)*(y1-y0)/(xz*xyz)
      DirCos(1, 3) = (x1-x0)/xyz
      
      DirCos(2, 1) = 0.0
      DirCos(2, 2) = xz/xyz
      DirCos(2, 3) = (y1-y0)/xyz
      
      DirCos(3, 1) = -(x1-x0)/xz
      DirCos(3, 2) = -(y1-y0)*(z1-z0)/(xz*xyz)
      DirCos(3, 3) = (z1-z0)/xyz

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
SUBROUTINE TaperCalc(R1, R2, H, taperV, h_c)
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
   
END SUBROUTINE TaperCalc
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
SUBROUTINE MarineGrowthPartSegment(R1, R2, Rmg1, Rmg2, L, rho,  Vinner, Vouter, m_mg, h_c, Ilmg, Irmg)
   
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
   call TaperCalc(R1, R2, L, Vinner, cVinner) 

   ! get V and CV for marine growth displacement
   call TaperCalc(Rmg1, Rmg2, L, Vouter, cVouter) 
   
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

END SUBROUTINE MarineGrowthPartSegment
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FloodedBallastPartSegment(R1, R2, L, rho, V, m, h_c, Il, Ir)
   
   REAL(ReKi),                     INTENT    ( IN    )  :: R1  ! interior radius of element at node point
   REAL(ReKi),                     INTENT    ( IN    )  :: R2  ! interior radius of other end of part-element
   REAL(ReKi),                     INTENT    ( IN    )  :: L   ! distance positive along axis to end of part-element
   REAL(ReKi),                     INTENT    ( OUT   )  :: V   ! volume from inner radius
   REAL(ReKi),                     INTENT    ( IN    )  :: rho ! density of ballast
   REAL(ReKi),                     INTENT    ( OUT   )  :: m   ! mass of material
   REAL(ReKi),                     INTENT    ( OUT   )  :: h_c    ! center of mass offset from first node
   REAL(ReKi),                     INTENT    ( OUT   )  :: Il   ! moment of inertia about axis
   REAL(ReKi),                     INTENT    ( OUT   )  :: Ir   ! moment of inertia about radial axis from first node
   

   
   ! get V and CV for flooded part of part-element
   call TaperCalc(R1, R2, L, V, h_c) 
   m = rho*V
   
   call CylInertia(R1, R2, L, rho, Il, Ir)  ! inertias for filled section

END SUBROUTINE FloodedBallastPartSegment
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE WriteSummaryFile( UnSum, MSL2SWL, numJoints, numNodes, nodes, numMembers, members, &
                             NOutputs, OutParam, MOutLst, JOutLst, uMesh, yMesh, p, m, errStat, errMsg ) 
                             
   INTEGER,                               INTENT ( IN    )  :: UnSum
   REAL(ReKi),                            INTENT ( IN    )  :: MSL2SWL
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
   TYPE(Morison_NodeType)                      ::  node1, node2
   real(ReKi)                                  :: ptLoad(6)
   logical                                     :: fillFlag
   type(Morison_MemberType)                    :: mem
   REAL(ReKi)                                  :: Cd1, Cd2, Ca1, Ca2, Cp1, Cp2, AxCd1, AxCd2, AxCa1, AxCa2, AxCp1, AxCp2, Cb1, Cb2, JAxCd1, JAxCd2, JAxCa1, JAxCa2, JAxCp1, JAxCp2 ! tmp coefs
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
         
      !  IF ( node2%Position(3) <= MSL2SWL .AND. node1%Position(3) >= -WtrDpth) totalDisplVol = totalDisplVol + elementVol

 
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
         
         if ( yMesh%Position(3,J) <= MSL2SWL ) then  ! need to check relative to MSL2SWL offset because the Mesh Positons are relative to MSL
            
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
      WRITE( UnSum, '(1X,A5,21(2X,A10))' ) '  i  ', '  MbrIndx ', '   Nxi    ', '   Nyi    ', '   Nzi    ', '     R    ', '    t     ', '   tMG    ', '  MGDens  ', ' PropPot  ', 'FilledFlag', 'FilledMass', '    Cd    ', '    Ca    ', '    Cp    ', '    Cb    ', '   AxCd   ',  '   AxCa   ', '   AxCp   ', '   JAxCd  ', '   JAxCa  ', '   JAxCp  '
      WRITE( UnSum, '(1X,A5,21(2X,A10))' ) ' (-) ', '    (-)   ', '   (m)    ', '   (m)    ', '   (m)    ', '    (m)   ', '   (m)    ', '   (m)    ', ' (kg/m^3) ', '   (-)    ', '   (-)    ', '  (kg)    ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ',  '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   '
      
         ! Write the node data
      do I = 1,numJoints   
         ! need to add MSL2SWL offset from this because the Positons are relative to SWL, but we should report them relative to MSL here
         pos = nodes(i)%Position
         pos(3) = pos(3) + MSL2SWL
         write( UnSum, '(1X,I5,(2X,A10),3(2X,F10.4),2(2X,A10),2(2X,ES10.3),10(2X,A10),3(2X,ES10.3))' ) i,'    -     ', pos, '    -     ',  '    -     ',  nodes(i)%tMG,  nodes(i)%MGdensity,  '    -     ',  '    -     ',  '    -     ', '    -     ',  '    -     ',  '    -     ',  '    -     ',  '    -     ',  '    -     ',  '    -     ',  nodes(i)%JAxCd,  nodes(i)%JAxCa, nodes(i)%JAxCp
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
            pos(3) = pos(3) + MSL2SWL
            if (members(j)%flipped) then
               II=members(j)%NElements+2-I
            else
               II=I
            endif
            write( UnSum, '(1X,I5,(2X,I10),3(2X,F10.4),4(2X,ES10.3),2(6X,L6),8(2X,ES10.3),3(7x,A5))' ) c, members(j)%MemberID, pos, members(j)%R(ii),  members(j)%R(ii)-members(j)%Rin(ii),  members(j)%tMG(ii),  members(j)%MGdensity(ii),  members(j)%PropPot,  fillFlag,  members(j)%m_fb_u(ii)+members(j)%m_fb_l(ii),  members(j)%Cd(ii),  members(j)%Ca(ii),  members(j)%Cp(ii),  members(j)%Cb(ii),  members(j)%AxCd(ii),  members(j)%AxCa(ii),  members(j)%AxCp(ii), '  -  ',  '  -  ',  '  -  '
         end do
      end do
      
      
      write( UnSum,  '(//)' ) 
      write( UnSum,  '(A8)' ) 'Members'
      write( UnSum,  '(/)' ) 
      write( UnSum, '(1X,A8,2X,A6,2X,A6,33(2X,A12))' ) 'MemberID', 'joint1','joint2','  Length  ', '   NElem    ', '   Volume   ', '  MGVolume  ', '      R1    ', '     t1     ', '      R2    ', '     t2     ', ' PropPot  ', 'FilledFlag', 'FillDensity', '  FillFSLoc ', '  FillMass  ', '     Cd1    ', '    Ca1   ', '     Cp1    ', '     Cb1    ', '    AxCd1   ', '    AxCa1   ', '    AxCp1   ', '   JAxCd1   ', '   JAxCa1   ', '  JAxCp1    ', '     Cd2    ', '     Ca2    ', '     Cp2    ', '     Cb2    ', '    AxCd2   ', '    AxCa2   ', '    AxCp2   ', '   JAxCd2   ', '   JAxCa2   ', '   JAxCp2   '
      write( UnSum, '(1X,A8,2X,A6,2X,A6,33(2X,A12))' ) '  (-)   ', ' (-)  ',' (-)  ','   (m)    ', '    (-)     ', '   (m^3)    ', '   (m^3)    ', '      (m)   ', '     (m)    ', '      (m)   ', '     (m)    ', '   (-)    ', '   (-)    ', ' (kg/m^3)  ', '     (-)    ', '    (kg)    ', '     (-)    ', '    (-)   ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '    (-)     ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    '
      
      
      
      
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
     
         Cd1   = members(i)%Cd(1)
         Cd2   = members(i)%Cd(N+1)
         Ca1   = members(i)%Ca(1)
         Ca2   = members(i)%Ca(N+1)
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
       
         
         write( UnSum, '(1X,I8,2X,I6,2X,I6,2X,ES12.5,2X,I12, 6(2X,ES12.5),2(2X,L12),23(2X,ES12.5))' )  members(i)%MemberID, &
                       members(i)%NodeIndx(1), members(i)%NodeIndx(N+1), members(i)%RefLength, N, &
                       memberVol, MGvolume, members(i)%Rmg(1), members(i)%Rmg(1)-members(i)%Rin(1), &
                       members(i)%Rmg(N+1), members(i)%Rmg(N+1)-members(i)%Rin(N+1),  &
                       members(i)%PropPot, filledFlag, members(i)%FillDens, members(i)%FillFSLoc, &
                       mass_fill, Cd1, Ca1, Cp1, Cb1, AxCd1, AxCa1, AxCp1, JAxCd1, JAxCa1, JAxCp1, &
                       Cd2, Ca2, Cp2, Cb2, AxCd2, AxCa2, AxCp2, JAxCd2, JAxCa2, JAxCp2
      
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
            pos(3) = pos(3) + MSL2SWL
            pos2 = node2%Position
            pos2(3) = pos2(3) + MSL2SWL
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
            pos(3) = pos(3) + MSL2SWL
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
SUBROUTINE SetDepthBasedCoefs( z, tMG, NCoefDpth, CoefDpths, Cd, Ca, Cp, AxCd, AxCa, AxCp, Cb )
   
   REAL(ReKi), INTENT (IN   )             :: z ! Z location relative to MSL inertial system
   REAL(ReKi), INTENT (IN   )             :: tMG
   INTEGER,    INTENT (IN   )             :: NCoefDpth
   TYPE(Morison_CoefDpths), INTENT (IN   ):: CoefDpths(:)
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
      AxCd   = CoefDpths(indx1)%DpthCd*(1-s)     + CoefDpths(indx2)%DpthAxCd*s
      AxCa   = CoefDpths(indx1)%DpthCa*(1-s)     + CoefDpths(indx2)%DpthAxCa*s
      AxCp   = CoefDpths(indx1)%DpthCp*(1-s)     + CoefDpths(indx2)%DpthAxCp*s
      Cb     = CoefDpths(indx1)%DpthCb*(1-s)     + CoefDpths(indx2)%DpthCb*s
   end if
   

END SUBROUTINE SetDepthBasedCoefs



!====================================================================================================
SUBROUTINE SetExternalHydroCoefs(  MSL2SWL, MCoefMod, MmbrCoefIDIndx, SimplCd, SimplCdMG, SimplCa, SimplCaMG, SimplCp, &
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
   type(Morison_CoefMembers), allocatable, intent(in   )  :: CoefMembers(:)
   integer(IntKi),                         intent(in   )  :: NCoefDpth
   type(Morison_CoefDpths),   allocatable, intent(in   )  :: CoefDpths(:)
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
         CALL SetDepthBasedCoefs( nodes(member%NodeIndx(i))%Position(3)+MSL2SWL,  member%tMG(i), NCoefDpth, CoefDpths, member%Cd(i), member%Ca(i), &
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
  
end subroutine SetExternalHydroCoefs

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
         ! Linearly interpolate the coef values based on depth
         !CALL FindInterpFactor( z, CoefDpths(indx1)%Dpth, CoefDpths(indx2)%Dpth, s )
      
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
   call AllocAry(member%dRdl_mg      , member%NElements,   'member%dRdl_mg'      , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%dRdl_mg_b    , member%NElements,   'member%dRdl_mg_b'    , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%dRdl_in      , member%NElements,   'member%dRdl_in'      , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
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
   call AllocAry(member%I_rfb_l      , member%NElements,   'member%I_rfb_l      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%I_rfb_u      , member%NElements,   'member%I_rfb_u      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%m_mg_l       , member%NElements,   'member%m_mg_l       ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%m_mg_u       , member%NElements,   'member%m_mg_u       ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%h_cmg_l      , member%NElements,   'member%h_cmg_l      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%h_cmg_u      , member%NElements,   'member%h_cmg_u      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%I_lmg_l      , member%NElements,   'member%I_lmg_l      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%I_lmg_u      , member%NElements,   'member%I_lmg_u      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%I_rmg_l      , member%NElements,   'member%I_rmg_l      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%I_rmg_u      , member%NElements,   'member%I_rmg_u      ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%Cfl_fb       , member%NElements,   'member%Cfl_fb       ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%Cfr_fb       , member%NElements,   'member%Cfr_fb       ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%CM0_fb       , member%NElements,   'member%CM0_fb       ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName) 
   call AllocAry(member%R            , member%NElements+1, 'member%R            ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%RMG          , member%NElements+1, 'member%RMG          ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%RMGB         , member%NElements+1, 'member%RMGB         ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%Rin          , member%NElements+1, 'member%Rin          ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%tMG          , member%NElements+1, 'member%tMG          ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%MGdensity    , member%NElements+1, 'member%MGdensity    ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%Cd           , member%NElements+1, 'member%Cd           ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%Ca           , member%NElements+1, 'member%Ca           ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%Cp           , member%NElements+1, 'member%Cp           ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%AxCd         , member%NElements+1, 'member%AxCd         ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%AxCa         , member%NElements+1, 'member%AxCa         ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%AxCp         , member%NElements+1, 'member%AxCp         ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%Cb           , member%NElements+1, 'member%Cb           ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_D    , 6, member%NElements+1, 'memberLoads%F_D'   , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_A    , 6, member%NElements+1, 'memberLoads%F_A'   , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_B    , 6, member%NElements+1, 'memberLoads%F_B'   , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_BF   , 6, member%NElements+1, 'memberLoads%F_BF'  , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_I    , 6, member%NElements+1, 'memberLoads%F_I'   , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_If   , 6, member%NElements+1, 'memberLoads%F_If'  , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_WMG  , 6, member%NElements+1, 'memberLoads%F_WMG' , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_IMG  , 6, member%NElements+1, 'memberLoads%F_IMG' , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   
   if (ErrStat >= AbortErrLev) return

   ! Initialize everything to zero
   member%NodeIndx      = 0.0_ReKi
   member%dRdl_mg       = 0.0_ReKi
   member%dRdl_mg_b     = 0.0_ReKi
   member%dRdl_in       = 0.0_ReKi
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
   member%I_rfb_l       = 0.0_ReKi
   member%I_rfb_u       = 0.0_ReKi
   member%m_mg_l        = 0.0_ReKi
   member%m_mg_u        = 0.0_ReKi
   member%h_cmg_l       = 0.0_ReKi
   member%h_cmg_u       = 0.0_ReKi
   member%I_lmg_l       = 0.0_ReKi
   member%I_lmg_u       = 0.0_ReKi
   member%I_rmg_l       = 0.0_ReKi
   member%I_rmg_u       = 0.0_ReKi
   member%Cfl_fb        = 0.0_ReKi
   member%Cfr_fb        = 0.0_ReKi
   member%CM0_fb        = 0.0_ReKi
   member%R             = 0.0_ReKi
   member%RMG           = 0.0_ReKi
   member%RMGB          = 0.0_ReKi
   member%Rin           = 0.0_ReKi
   member%tMG           = 0.0_ReKi
   member%MGdensity     = 0.0_ReKi
   member%Cd            = 0.0_ReKi
   member%Ca            = 0.0_ReKi
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
      
   end if    
   
end subroutine FlipMemberNodeData
!----------------------------------------------------------------------------------------------------------------------------------
subroutine SetMemberProperties( MSL2SWL, gravity, member, MCoefMod, MmbrCoefIDIndx, MmbrFilledIDIndx, propSet1, propSet2, InitInp, errStat, errMsg )
   real(ReKi),                   intent (in   )  :: MSL2SWL
   real(ReKi),                   intent (in   )  :: gravity
   type(Morison_MemberType),     intent (inout)  :: member
   integer(IntKi),               intent (in   )  :: MCoefMod
   integer(IntKi),               intent (in   )  :: MmbrCoefIDIndx
   integer(IntKi),               intent (in   )  :: MmbrFilledIDIndx
   type(Morison_MemberPropType), intent (in   )  :: propSet1             ! property set of node 1
   type(Morison_MemberPropType), intent (in   )  :: propSet2             ! property set of node N+1
   type(Morison_InitInputType),  intent (in   )  :: InitInp
   integer(IntKi),               intent (  out)  :: errStat              ! returns a non-zero value when an error occurs            
   character(*),                 intent (  out)  :: errMsg               ! Error message if errStat /= ErrID_None

   integer(IntKi) :: N, i
   real(ReKi)     :: WtrDepth,s, dl
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
   
   WtrDepth = InitInp%WtrDpth
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
      call SetNodeMG( InitInp%NMGDepths, InitInp%MGDepths, InitInp%Nodes(member%NodeIndx(i)), InitInp%MSL2SWL, member%tMG(i), member%MGDensity(i) )
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

   call SetExternalHydroCoefs(  MSL2SWL, MCoefMod, MmbrCoefIDIndx, InitInp%SimplCd, InitInp%SimplCdMG, InitInp%SimplCa, InitInp%SimplCaMG, InitInp%SimplCp, &
                                   InitInp%SimplCpMG, InitInp%SimplAxCd, InitInp%SimplAxCdMG, InitInp%SimplAxCa, InitInp%SimplAxCaMG, InitInp%SimplAxCp, InitInp%SimplAxCpMG, &
                                   InitInp%SimplCb, InitInp%SimplCbMG, InitInp%SimplMCF, & 
                                   InitInp%CoefMembers, InitInp%NCoefDpth, InitInp%CoefDpths, InitInp%Nodes, member )
   
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
         CALL SetErrStat(ErrID_Fatal, 'MacCamy-Fuchs members must be surface piercing.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, 'SetMemberProperties' )   
         RETURN
      END IF
      ! Check inclination
      If ( ABS(phi) .GE. 0.174533 ) THEN ! If inclination from vertical is greater than 10 deg
         CALL SetErrStat(ErrID_Fatal, 'MacCamy-Fuchs members must be within 10 degrees from vertical.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, 'SetMemberProperties' )   
         RETURN
      END IF
      ! Check radius
      DO i = 1, member%NElements+1
         IF ( (member%RMG(i) .GT. 1.1_ReKi*REAL(0.5_SiKi*InitInp%MCFD)) .OR. (member%RMG(i) .LT. 0.9_ReKi*REAL(0.5_SiKi*InitInp%MCFD)) ) THEN
            ! Error because MacCamy-Fuchs members must have a diameter within +/-10% of MCFD specified in seastate.
            CALL SetErrStat(ErrID_Fatal, 'MacCamy-Fuchs members must have a diameter within +/-10% of MCFD specified in the SeaState input file.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, 'SetMemberProperties' )   
            RETURN
         END IF
      END DO
      ! Check draft-to-radius ratio
      IF ( (-InitInp%Nodes(member%NodeIndx(1))%Position(3)) < 0.5_SiKi*InitInp%MCFD ) THEN
         CALL SetErrStat(ErrID_Fatal, 'Initial draft of MacCamy-Fuchs members should be at least as large as their radius.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, 'SetMemberProperties' )   
         RETURN
      END IF
   END IF

   ! find fill location of member (previously in SetElementFillProps)
   member%MmbrFilledIDIndx = MmbrFilledIDIndx ! Set this to the parameter version of this member data
   if ( MmbrFilledIDIndx > 0 ) then    
      member%FillDens     =  InitInp%FilledGroups(MmbrFilledIDIndx)%FillDens
      member%FillFSLoc    =  InitInp%FilledGroups(MmbrFilledIDIndx)%FillFSLoc - InitInp%MSL2SWL
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
         if ( Zb <= -InitInp%WtrDpth ) then
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
      if (member%MHstLMod == 1) then
         if ( abs(Zb) < abs(member%Rmg(N+1)*sinPhi) ) then
            call SetErrStat(ErrID_Fatal, 'The upper end-plate of a member must not cross the water plane.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, 'SetMemberProperties' )   
         end if
         if ( abs(Za) < abs(member%Rmg(1)*sinPhi) ) then
            call SetErrStat(ErrID_Fatal, 'The lower end-plate of a member must not cross the water plane.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, 'SetMemberProperties' )   
         end if
      end if
      if ( ( Za < -WtrDepth .and. Zb >= -WtrDepth ) .and. ( phi > 10.0*d2r .or. abs((member%RMG(N+1) - member%RMG(1))/member%RefLength)>0.1 ) ) then
         call SetErrStat(ErrID_Fatal, 'A member which crosses the seabed must not be inclined more than 10 degrees from vertical or have a taper larger than 0.1.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, 'SetMemberProperties' )   
      end if
      
   end if
   

   ! calculate h_floor if seabed-piercing
   member%h_floor = 0.0_ReKi
   member%i_floor = member%NElements+1  ! Default to entire member is below the seabed
   member%doEndBuoyancy = .false.
   if (Za < -WtrDepth) then
      do i= 2, member%NElements+1
         Za = InitInp%Nodes(member%NodeIndx(i))%Position(3)
         if (Za > -WtrDepth) then            ! find the lowest node above the seabed
            
            if (cosPhi < 0.173648178 ) then ! phi > 80 degrees and member is seabed crossing
               call SetErrStat(ErrID_Fatal, 'A seabed crossing member must have an inclination angle of <= 80 degrees from vertical.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, 'SetMemberProperties' )
            end if
            
            member%h_floor = (-WtrDepth-Za)/cosPhi  ! get the distance from the node to the seabed along the member axis (negative value)
            member%i_floor = i-1                    ! record the number of the element that pierces the seabed
            member%doEndBuoyancy = .true.
            exit
         else if ( EqualRealNos(Za, -WtrDepth ) ) then
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
      
      member%alpha(   i) = GetAlpha(member%RMGB(i), member%RMGB(i+1))   ! Only used to distribute external buoyancy load to nodes
      member%alpha_fb(i) = GetAlpha(member%Rin( i), member%Rin( i+1))
      
   end do

   member%Vinner = 0.0_ReKi  ! Total  volume of member without marine growth
   member%Vouter = 0.0_ReKi  ! Total outer volume of member including marine growth
   member%Vballast = 0.0_ReKi ! Total ballasted volume of member
   
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

         CALL MarineGrowthPartSegment(member%R(i  ), Rmid, member%RMG(i  ),RmidMG, Lmid, member%MGDensity(i),  Vinner_l, Vouter_l, member%m_mg_l(i), member%h_cmg_l(i), member%I_lmg_l(i), member%I_rmg_l(i))   ! get precomputed quantities for lower half-segment
         CALL MarineGrowthPartSegment(member%R(i+1), Rmid, member%RMG(i+1),RmidMG,-Lmid, member%MGDensity(i),  Vinner_u, Vouter_u, member%m_mg_u(i), member%h_cmg_u(i), member%I_lmg_u(i), member%I_rmg_u(i))   ! get precomputed quantities for upper half-segment
         
      else if (i == member%i_floor) then         
         ! crossing seabed: get the properties for part-element above the seabed and lump to the upper node      

         Rmid   = (-member%h_floor*member%R(  i) +(dl+member%h_floor)*member%R(  i+1))/dl
         RmidMG = (-member%h_floor*member%RMG(i) +(dl+member%h_floor)*member%RMG(i+1))/dl
         Lmid   = -member%h_floor

         CALL MarineGrowthPartSegment(member%R(i+1), Rmid, member%RMG(i+1),RmidMG, -Lmid, member%MGDensity(i),  Vinner_u, Vouter_u, member%m_mg_u(i), member%h_cmg_u(i), member%I_lmg_u(i), member%I_rmg_u(i))   ! get precomputed quantities for upper half-segment
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
            CALL FloodedBallastPartSegment(member%Rin(i  ), Rmidin,  Lmid, member%FillDens, Vballast_l, member%m_fb_l(i), member%h_cfb_l(i), member%I_lfb_l(i), member%I_rfb_l(i))   ! get precomputed quantities for lower half-segment
            CALL FloodedBallastPartSegment(member%Rin(i+1), Rmidin, -Lmid, member%FillDens, Vballast_u, member%m_fb_u(i), member%h_cfb_u(i), member%I_lfb_u(i), member%I_rfb_u(i))   ! get precomputed quantities for upper half-segment
 
         ! partially filled element, so split at FillFSLoc
         else if ((i > member%i_floor)  .AND. (member%FillFSLoc < Zb)) then

            ! get the properties for each partial-element lumped to the appropriate node
            Lmid   = member%FillFSLoc - Za 
            Rmidin = member%Rin(i)+(Lmid/(Zb-Za))*(member%Rin(i+1)-member%Rin(i))  ! radius of member interior at middle of segment, where division occurs
            CALL FloodedBallastPartSegment(member%Rin(i  ), Rmidin,  Lmid, member%FillDens, Vballast_l, member%m_fb_l(i), member%h_cfb_l(i), member%I_lfb_l(i), member%I_rfb_l(i))   ! get precomputed quantities for lower half-segment
            CALL FloodedBallastPartSegment(member%Rin(i+1), Rmidin, -Lmid, 0.0_ReKi, Vballast_u, member%m_fb_u(i), member%h_cfb_u(i), member%I_lfb_u(i), member%I_rfb_u(i))   ! get precomputed quantities for upper half-segment
 
         else if (i == member%i_floor) then     ! Hopefully we don't have a partially filled element crossing the seabed.
 
            ! crossing seabed: get the properties for part-element above the seabed and lump to the upper node
            RmidMG = (-member%h_floor*member%RMG(i) +(dl+member%h_floor)*member%RMG(i+1))/dl
            Rmidin = (-member%h_floor*member%Rin(i) +(dl+member%h_floor)*member%Rin(i+1))/dl
            Lmid   = -member%h_floor
            CALL FloodedBallastPartSegment(member%Rin(i+1), Rmidin, -Lmid, member%FillDens,  Vballast_u, member%m_fb_u(i), member%h_cfb_u(i), member%I_lfb_u(i), member%I_rfb_u(i))   ! get precomputed quantities for upper half-segment
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
         
         if (Zb < -WtrDepth) then
            ! fully buried element, do not add these volume contributions to totals
         else if (0.0 >= Zb) then   ! Bug fix per OpenFAST issue #844   GJH 2/3/2022
            ! fully submerged elements.  
            ! NOTE: For an element which is fractionaly in the seabed, the entire element volume is added to totals
            member%Vinner = member%Vinner + Vinner_l + Vinner_u
            member%Vouter = member%Vouter + Vouter_l + Vouter_u
            member%Vsubmerged = member%Vsubmerged + Vouter_l + Vouter_u
         else if ((0.0 > Za) .AND. (0.0 < Zb)) then ! Bug fix per OpenFAST issue #844   GJH 2/3/2022
            ! if (i == 1) then
            !    call SetErrStat(ErrID_Fatal, 'The lowest element of a member must not cross the free surface.  This is true for MemberID '//trim(num2lstr(member%MemberID)), errStat, errMsg, 'SetMemberProperties')
            ! end if
            
            ! partially submerged element
            member%Vinner = member%Vinner + Vinner_l + Vinner_u
            member%Vouter = member%Vouter + Vouter_l + Vouter_u
            ! compute volume portion which is submerged
            Lmid = -Za/cosPhi 
            call TaperCalc( member%Rmg(i), member%Rmg(i)+Lmid*member%dRdl_mg(i), Lmid, Vouter_l, h_c)
            
            member%Vsubmerged = member%Vsubmerged + Vouter_l 
            
         else ! fully above the water
            member%Vinner = member%Vinner + Vinner_l + Vinner_u
            member%Vouter = member%Vouter + Vouter_l + Vouter_u
         end if 
      end if
      
      ! ------------------ flooded ballast weight (done) --------------------
      ! NOTE: this section of code is somewhat redundant with "flooded ballast inertia" section above

      li = dl*(i-1)
      ! fully buried element
      if (Zb < -WtrDepth) then
         member%floodstatus(i) = 0
      
      ! fully filled elements 
      else if (member%memfloodstatus > 0 .and. member%FillFSLoc > Zb) then  
         member%floodstatus(i) = 1
         member%Vballast = member%Vballast + Vballast_l + Vballast_u
         ! depth-adjusted force distribution constant
         member%alpha_fb_star(i) = member%alpha_fb(i)*( Zb - member%FillFSLoc )**3 / ( ( (1-member%alpha_fb(i))*(Za - member%FillFSLoc))**3 + member%alpha_fb(i)*(Zb - member%FillFSLoc)**3 )
         
         ! force and moment magnitude constants

         
         member%Cfl_fb(i) = TwoPi * member%dRdl_in(i) * member%FillDens * gravity * dl *( (li - member%l_fill)*member%Rin(i) + 0.5*((li - member%l_fill)* member%dRdl_in(i) + member%Rin(i))*dl + 1.0/3.0* member%dRdl_in(i)*dl**2 )
         member%Cfr_fb(i) =    Pi *                     member%FillDens * gravity * dl *( member%Rin(i)**2 +  member%dRdl_in(i)*member%Rin(i)*dl +1.0/3.0 * member%dRdl_in(i)**2 *dl**2 )
         member%CM0_fb(i) = TwoPi *                     member%FillDens * gravity * dl *( 0.25*dl**3* member%dRdl_in(i)**4 + 0.25*dl**3* member%dRdl_in(i)**2 + dl**2* member%dRdl_in(i)**3*member%Rin(i) + 2.0/3.0*dl**2* member%dRdl_in(i)*member%Rin(i) + 1.5*dl* member%dRdl_in(i)**2*member%Rin(i)**2 + 0.5*dl*member%Rin(i)**2 +  member%dRdl_in(i)*member%Rin(i)**3 )
         
         
      ! partially filled element
      else if ((member%memfloodstatus > 0) .and. (member%FillFSLoc > Za) .AND. (member%FillFSLoc < Zb)) then
         
         ! Need to enforce the modeling requirement that the first/bottom-most element of a member be fully flooded
         if (i == 1) then
            call SetErrStat(ErrID_Fatal,'The modeling of partially flooded/ballested members requires that the first/bottom-most element of a member must be fully flooded. This is not true for MemberID '//trim(num2lstr(member%MemberID)),ErrStat,ErrMsg,'SetMemberProperties')
            return
         end if
         ! Need to enforce the modeling requirement that a partially flooded member must not be close to horizontal
         if ( (InitInp%Nodes(member%NodeIndx(N+1))%Position(3) - member%Rin(N+1)*sinPhi) < member%FillFSLoc ) then
            call SetErrStat(ErrID_Fatal,'The modeling of partially flooded/ballested members requires the the member not be near horizontal.  This is not true for MemberID '//trim(num2lstr(member%MemberID)),ErrStat,ErrMsg,'SetMemberProperties') 
            return
         end if
         
         member%floodstatus(i) = 2
         
         ! length along axis from node i to fill level
         member%h_fill = member%l_fill - (i-1)*dl
         !Since this element is only partially flooded/ballasted, compute the Volume fraction which is filled
         call TaperCalc( member%Rin(i), member%Rin(i)+member%h_fill*member%dRdl_in(i), member%h_fill, Vballast_l, h_c)
         Vballast_u = 0.0
         member%Vballast = member%Vballast + Vballast_l + Vballast_u ! Note: Vballast_l will match calculations above
       
         
         ! depth-adjusted force distribution constant
         member%alpha_fb_star(i) = (1 - member%alpha_fb(i))*( Za - member%FillFSLoc )**3 / ( ( (1-member%alpha_fb(i))*(Za - member%FillFSLoc))**3 - member%alpha_fb(i)*(Zb - member%FillFSLoc)**3 )
         
         ! force and moment magnitude constants
         member%Cfl_fb(i) = TwoPi * member%dRdl_in(i) * member%FillDens * gravity * member%h_fill *( (li - member%l_fill)*member%Rin(i) + 0.5*((li - member%l_fill)*member%dRdl_in(i) + member%Rin(i))*member%h_fill + 1.0/3.0*member%dRdl_in(i)*member%h_fill**2 )
         member%Cfr_fb(i) =    Pi * member%FillDens * gravity * member%h_fill *( member%Rin(i)**2 + member%dRdl_in(i)*member%Rin(i)*member%h_fill +1.0/3.0 *member%dRdl_in(i)**2 *member%h_fill**2 )
         member%CM0_fb(i) = TwoPi * member%FillDens * gravity * member%h_fill *( 0.25*member%h_fill**3*member%dRdl_in(i)**4 + 0.25*member%h_fill**3*member%dRdl_in(i)**2 + member%h_fill**2*member%dRdl_in(i)**3*member%Rin(i) + 2.0/3.0*member%h_fill**2*member%dRdl_in(i)*member%Rin(i)  &
                                                                                 + 1.5*member%h_fill*member%dRdl_in(i)**2*member%Rin(i)**2 + 0.5*member%h_fill*member%Rin(i)**2 + member%dRdl_in(i)*member%Rin(i)**3 ) &
                                    -0.25 * member%FillDens * gravity * Pi * (  member%Rin(i) + member%h_fill*member%dRdl_in(i))**4
      
      ! unflooded element
      else
         member%floodstatus(i) = 0
      
      end if
      

   end do ! end looping through elements   
  
 
end subroutine SetMemberProperties

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
      p%Members(i)%MemberID  = InitInp%InpMembers(i)%MemberID
      p%Members(i)%RefLength = InitInp%InpMembers(i)%RefLength
      p%Members(i)%dl        = InitInp%InpMembers(i)%dl
      p%Members(i)%NElements = InitInp%InpMembers(i)%NElements
      p%Members(i)%PropPot   = InitInp%InpMembers(i)%PropPot
      p%Members(i)%MHstLMod  = InitInp%InpMembers(i)%MHstLMod
      ! p%Members(i)%MCF       = InitInp%InpMembers(i)%MCF
      
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
      call SetMemberProperties( InitInp%MSL2SWL, InitInp%Gravity, p%Members(i), InitInp%InpMembers(i)%MCoefMod, InitInp%InpMembers(i)%MmbrCoefIDIndx, InitInp%InpMembers(i)%MmbrFilledIDIndx, InitInp%MPropSets(prop1Indx), InitInp%MPropSets(prop2Indx), InitInp, errStat2, errMsg2 ) 
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

   TYPE(Morison_MemberType) :: member    ! the current member
   INTEGER                  :: i, j, k
   REAL(ReKi)               :: v2D(3,1), pos(3)
   real(ReKi)               :: An(3), An_drag(3), Vn(3), I_n(3), sgn, Amag, Amag_drag, Vmag, Imag, Ir_MG_end, Il_MG_end, R_I(3,3), IRl_mat(3,3), tMG, MGdens
   integer(IntKi)           :: MemberEndIndx
   INTEGER, ALLOCATABLE       :: commonNodeLst(:)
   LOGICAL, ALLOCATABLE       :: usedJointList(:)
   integer(IntKi)           :: errStat2              ! returns a non-zero value when an error occurs            
   CHARACTER(errMsgLen)     :: errMsg2     ! Error message if errStat2 /= ErrID_None

      
      ! Initialize errStat        
   errStat = ErrID_None         
   errMsg  = ""               
      
  

      ! Define parameters here:  
   p%DT         = Interval
   p%WtrDens    = InitInp%WtrDens
   p%WtrDpth    = InitInp%WtrDpth
   p%Gravity    = InitInp%Gravity
   p%NNodes     = InitInp%NNodes
   p%NJoints    = InitInp%NJoints
   p%NStepWave  = InitInp%NStepWave
   p%NumOuts    = InitInp%NumOuts
   p%NMOutputs  = InitInp%NMOutputs                       ! Number of members to output [ >=0 and <10]
   p%MSL2SWL    = InitInp%MSL2SWL
   p%WaveDisp   = InitInp%WaveDisp
   p%AMMod      = InitInp%AMMod
   p%WaveStMod  = InitInp%WaveStMod

   ! Only compute added-mass force up to the free surface if wave stretching is enabled
   IF ( p%WaveStMod .EQ. 0_IntKi ) THEN
       ! Setting AMMod to zero just in case. Probably redundant.
       p%AMMod = 0_IntKi
   END IF

   p%WaveElev1  => InitInp%WaveElev1
   IF (associated(InitInp%WaveElev2)) THEN
      p%WaveElev2 => InitInp%WaveElev2
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
      call SetNodeMG( InitInp%NMGDepths, InitInp%MGDepths, InitInp%Nodes(i), InitInp%MSL2SWL, InitInp%Nodes(i)%tMG, InitInp%Nodes(i)%MGDensity )
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
      pos(3) = pos(3) + InitInp%MSL2SWL
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
   m%LastIndWave              = 1

   
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
      
      IF ( InitInp%InpJoints(i)%Position(3) >= -p%WtrDpth ) THEN
   
         ! loop through each member attached to the joint, getting the radius of its appropriate end
         DO J = 1, InitInp%InpJoints(I)%NConnections
      
            ! identify attached member and which end to us
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
               sgn = 1.0                                 ! Local coord sys points out of member at ending node, so leave sign of local z vector
            END IF

            ! Account for reordering of what the original node for the end was -- This affects the sign of the An term which can pose a problem for members crossing the waterline
            if (member%Flipped)   sgn = -1.0 * sgn
               
            ! Compute the signed quantities for this member end (for drag regardless of PropPot value), and add them to the joint values
            An_drag = An_drag + sgn* member%k*Pi*(member%RMG(MemberEndIndx))**2     ! area-weighted normal vector
               
            ! For the following quantities, the attached member cannot be modeled using WAMIT if we're to count it
            IF  (.NOT. member%PropPot) THEN

               ! Compute the signed quantities for this member end, and add them to the joint values
               An = An + sgn* member%k*Pi*(member%RMG(MemberEndIndx))**2     ! area-weighted normal vector
               Vn = Vn + sgn* member%k*   (member%RMG(MemberEndIndx))**3     ! r^3-weighted normal vector used for mass
               I_n=I_n + sgn* member%k*Pi*(member%RMG(MemberEndIndx))**4     ! r^4-weighted normal vector used for moments of inertia
               if (tMG == -999.0) then
                  ! All member nodes at this joint will have the same MG thickness and density, so only do this once
                  tMG = member%tMG(MemberEndIndx)
                  MGdens = member%MGdensity(MemberEndIndx) 
               end if
            END IF
         
         END DO   !J = 1, InitInp%InpJoints(I)%NConnections

         Vn = Vn*TwoPi/3.0_ReKi ! Semisphere volume is Vn = 2/3 pi \sum (r_MG^3 k)
         
         p%An_End(:,i) = An_drag 
         Amag_drag = Dot_Product(An_drag ,An_drag)
         Amag = Dot_Product(An ,An)
         IF (EqualRealNos(Amag_drag, 0.0_ReKi)) THEN
            p%DragConst_End(i) =  0.0
         ELSE
            p%DragConst_End(i) = InitInp%Nodes(i)%JAxCd*p%WtrDens / ( 4.0_ReKi * Amag_drag )
         END IF
         ! magnitudes of normal-weighted values
         Amag = sqrt(Amag)
         Vmag = sqrt(Dot_Product(Vn ,Vn))
         Imag = sqrt(Dot_Product(I_n,I_n))
      
         ! Constant part of the external hydrodynamic added mass term
         if ( Vmag > 0.0 ) then
            v2D(:,1) = Vn        
            p%AM_End(:,:,i) = (InitInp%Nodes(I)%JAxCa*InitInp%WtrDens/ Vmag)*matmul(v2D, transpose(v2D))
         end if
         
         ! Constant part of the external hydrodynamic dynamic pressure force
         if ( Amag > 0.0 ) then
            p%DP_Const_End(:,i) = -InitInp%Nodes(i)%JAxCp*An 
         endif
         
         ! marine growth mass/inertia magnitudes
         p%Mass_MG_End(i) = MGdens * tMG * Amag
         p%F_WMG_End(3,i) =        -MGdens * tMG * Amag * InitInp%Gravity  ! Z component of the directional force due to marine growth mass at joint
         Ir_MG_end   =  0.25 * MGdens * tMG * Imag  ! radial moment of inertia magnitude
         Il_MG_end   =  0.5  * MGdens * tMG * Imag  ! axial moment of inertia magnitude
      
         ! get rotation matrix for moment of inertia orientations
         call RodrigMat(I_n, R_I, errStat, errMsg)
         IF ( errStat >= AbortErrLev ) RETURN

         ! globally-oreinted moment of inertia matrix for joint
         Irl_mat = 0.0
         Irl_mat(1,1) = Ir_MG_end
         Irl_mat(2,2) = Ir_MG_end
         Irl_mat(3,3) = Il_MG_end
      
         p%I_MG_End(:,:,i) = MatMul( MatMul(R_I, Irl_mat), Transpose(R_I) ) ! final moment of inertia matrix for node
         

      END IF  ! InitInp%InpJoints(i)%Position(3) >= -p%WtrDpth
   
      p%DragMod_End   (i) = InitInp%Nodes(i)%JAxFDMod
      IF ( InitInp%Nodes(i)%JAxVnCOff .LE. 0.0_ReKi) THEN
         p%VRelNFiltConst(i) = 1.0_ReKi
         p%DragLoFSc_End (i) = 1.0_ReKi
      ELSE
         p%VRelNFiltConst(i) = exp(-2.0*Pi*InitInp%Nodes(i)%JAxVnCOff * p%DT)
         p%DragLoFSc_End (i) = InitInp%Nodes(i)%JAxFDLoFSc
      END IF

   END DO ! looping through nodes that are joints, i
          
         ! Define initial guess for the system inputs here:
         !    u%DummyInput = 0
         ! Define system output initializations (set up mesh) here:  
         ! Define initialization-routine output here:
         
   
   ! Setup the 4D grid information for the Interpolatin Module
   p%seast_interp_p = InitInp%seast_interp_p

   ! Setup 3D SWL grids needed for wave stretching
   
   IF (p%WaveStMod > 0_IntKi) THEN ! Wave stretching enabled
      
      ! Allocate variables for the wave dynamics at the SWL - Needed for wave stretching
      ALLOCATE ( p%WaveDynP0 (0:p%NStepWave,p%seast_interp_p%n(2),p%seast_interp_p%n(3)), STAT=errStat2 )
         IF ( errStat2 /= 0 ) call SetErrStat(ErrID_Fatal,'Error allocating space for p%WaveDynP0.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE ( p%WaveVel0 (0:p%NStepWave,p%seast_interp_p%n(2),p%seast_interp_p%n(3),3), STAT=errStat2 )
         IF ( errStat2 /= 0 ) call SetErrStat(ErrID_Fatal,'Error allocating space for p%WaveVel0.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE ( p%WaveAcc0 (0:p%NStepWave,p%seast_interp_p%n(2),p%seast_interp_p%n(3),3), STAT=errStat2 )
         IF ( errStat2 /= 0 ) call SetErrStat(ErrID_Fatal,'Error allocating space for p%WaveAcc0.', ErrStat, ErrMsg, RoutineName)
      ALLOCATE ( p%WaveAccMCF0 (0:p%NStepWave,p%seast_interp_p%n(2),p%seast_interp_p%n(3),3), STAT=errstat2 )
         IF ( errStat2 /= 0 ) call SetErrStat(ErrID_Fatal,'Error allocating space for p%WaveAccMCF0.', ErrStat, ErrMsg, RoutineName)
         
      if (ErrStat >= AbortErrLev) RETURN
   
      ! Copy the wave data at the SWL
      DO i = 1,p%seast_interp_p%n(2)
        DO j = 1,p%seast_interp_p%n(3)
          p%WaveDynP0(:,i,j) = p%WaveDynP(:,i,j,p%seast_interp_p%n(4))
          DO k = 1,3
            p%WaveVel0(:,i,j,k) = p%WaveVel(:,i,j,p%seast_interp_p%n(4),k)
            p%WaveAcc0(:,i,j,k) = p%WaveAcc(:,i,j,p%seast_interp_p%n(4),k)
          END DO
        END DO
      END DO
   
      ! Also copy the MacCamy-Fuchs scaled wave acceleration at the SWL if available
      IF (ASSOCIATED(p%WaveAccMCF)) THEN
        DO i = 1,p%seast_interp_p%n(2)
          DO j = 1,p%seast_interp_p%n(3)
            DO k = 1,3
              p%WaveAccMCF0(:,i,j,k) = p%WaveAccMCF(:,i,j,p%seast_interp_p%n(4),k)
            END DO
          END DO
        END DO
      END IF
   
   END IF
   
   ! Initialize the outputs      
   CALL MrsnOUT_Init( InitInp, y, p, InitOut, errStat2, errMsg2 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if ( errStat >= AbortErrLev ) return
   
   ! We will call CalcOutput to compute the loads for the initial reference position
   ! Then we can use the computed load components in the Summary File
   ! NOTE: Morison module has no states, otherwise we could no do this. GJH
   
   call Morison_CalcOutput(0.0_DbKi, u, p, x, xd, z, OtherState, y, m, errStat2, errMsg2 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if ( errStat >= AbortErrLev ) return
   
      ! Write Summary information to *HydroDyn* summary file now that everything has been initialized. 
   CALL WriteSummaryFile( InitInp%UnSum, InitInp%MSL2SWL, InitInp%NJoints, InitInp%NNodes, InitInp%Nodes, p%NMembers, p%Members, &
                          p%NumOuts, p%OutParam, p%MOutLst, p%JOutLst, u%Mesh, y%Mesh, p, m, errStat2, errMsg2 )
   call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   if ( errStat >= AbortErrLev ) return
                                                       
   !Contains:
   !   SUBROUTINE CleanUpInitOnErr
   !   IF (ALLOCATED(sw(1)%array))  DEALLOCATE(sw(1)%array, STAT=aviFail)
   !   END SUBROUTINE

END SUBROUTINE Morison_Init
!----------------------------------------------------------------------------------------------------------------------------------
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
FUNCTION GetAlpha(R1,R2)
   ! calculates relative center of volume location for a (tapered) cylindrical element
   real(ReKi)    :: GetAlpha
   REAL(ReKi),                     INTENT    ( IN    )  :: R1  ! interior radius of element at node point
   REAL(ReKi),                     INTENT    ( IN    )  :: R2  ! interior radius of other end of part-element
   
   IF ( EqualRealNos(R1, R2) ) THEN ! To cover the case where R1=R2=0
      GetAlpha = 0.5
   ELSE
      GetAlpha = (R1*R1 + 2.0*R1*R2 + 3.0*R2*R2)/4.0/(R1*R1 + R1*R2 + R2*R2)
   END IF
   
END FUNCTION GetAlpha

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AllocateNodeLoadVariables(InitInp, p, m, NNodes, errStat, errMsg )
   TYPE(Morison_InitInputType),       INTENT(IN   )  :: InitInp     ! Initialization inputs
   TYPE(Morison_ParameterType),       INTENT(INOUT)  :: p           ! parameter variables
   TYPE(Morison_MiscVarType),         INTENT(INOUT)  :: m           ! Misc/optimization variables            
   INTEGER(IntKi),                    INTENT(IN   )  :: NNodes      ! number of nodes in node list
   INTEGER(IntKi),                    INTENT(  OUT)  :: errStat     ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: errMsg      ! Error message if errStat /= ErrID_None
   integer(IntKi)           :: errStat2              ! returns a non-zero value when an error occurs            
   CHARACTER(errMsgLen)     :: errMsg2     ! Error message if errStat2 /= ErrID_None
   character(*), parameter :: routineName = 'AllocateNodeLoadVariables'
   
      ! Initialize errStat
         
   errStat = ErrID_None         
   errMsg  = ""               
      
   call AllocAry( m%nodeInWater        , NNodes   , 'm%nodeInWater'  , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%vrel         ,    3, NNodes   , 'm%vrel'         , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   !call AllocAry( m%F_D          ,    6, NNodes   , 'm%F_D'          , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   !call AllocAry( m%F_A          ,    6, NNodes   , 'm%F_A'          , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   !call AllocAry( m%F_B          ,    6, NNodes   , 'm%F_B'          , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   !call AllocAry( m%F_BF         ,    6, NNodes   , 'm%F_BF'         , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   !call AllocAry( m%F_I          ,    6, NNodes   , 'm%F_I'          , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   !call AllocAry( m%F_If         ,    6, NNodes   , 'm%F_If'         , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   !call AllocAry( m%F_WMG        ,    6, NNodes   , 'm%F_WMG'        , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   !call AllocAry( m%F_IMG        ,    6, NNodes   , 'm%F_IMG'        , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%FV           ,    3, NNodes   , 'm%FV'           , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%FA           ,    3, NNodes   , 'm%FA'           , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%FAMCF        ,    3, NNodes   , 'm%FAMCF'        , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%FDynP        ,       NNodes   , 'm%FDynP'        , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%WaveElev     ,       NNodes   , 'm%WaveElev'        , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%WaveElev1    ,       NNodes   , 'm%WaveElev1'        , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%WaveElev2    ,       NNodes   , 'm%WaveElev2'        , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%An_End       ,    3, p%NJoints, 'p%An_End'       , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%DragConst_End,       p%NJoints, 'p%DragConst_End', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%F_I_End      ,    3, p%NJoints, 'm%F_I_End'      , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%F_BF_End     ,    6, p%NJoints, 'm%F_BF_End'     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%F_A_End      ,    3, p%NJoints, 'm%F_A_End'      , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%F_D_End      ,    3, p%NJoints, 'm%F_D_End'      , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%F_B_End      ,    6, p%NJoints, 'm%F_B_End'      , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%F_IMG_End    ,    6, p%NJoints, 'm%F_IMG_End'    , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%I_MG_End     , 3, 3, p%NJoints, 'p%I_MG_End'     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%F_WMG_End    ,    3, p%NJoints, 'p%F_WMG_End'    , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%Mass_MG_End  ,       p%NJoints, 'p%Mass_MG_End'  , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%AM_End       , 3, 3, p%NJoints, 'p%AM_End'       , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%DP_Const_End ,    3, p%NJoints, 'p%DP_Const_End' , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)

   call AllocAry( m%V_rel_n        ,     p%NJoints, 'm%V_rel_n'        , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( m%V_rel_n_HiPass ,     p%NJoints, 'm%V_rel_n_HiPass' , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)

   call AllocAry( p%DragMod_End ,     p%NJoints, 'p%DragMod_End' , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%DragLoFSc_End ,     p%NJoints, 'p%DragLoFSc_End' , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( p%VRelNFiltConst ,     p%NJoints, 'p%VRelNFiltConst' , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)

   if (errStat >= AbortErrLev) return
   
   m%nodeInWater   = 0
   m%vrel          = 0.0_ReKi
   !m%F_D           = 0.0_ReKi
   !m%F_A           = 0.0_ReKi 
   !m%F_B           = 0.0
   !m%F_BF          = 0.0
   !m%F_I           = 0.0
   !m%F_If          = 0.0
   !m%F_WMG         = 0.0
   !m%F_IMG         = 0.0
   m%FV            = 0.0_ReKi
   m%FA            = 0.0_ReKi
   m%FDynP         = 0.0_ReKi
   p%An_End        = 0.0
   p%DragConst_End = 0.0
   m%F_I_End       = 0.0
   m%F_BF_End      = 0.0
   m%F_A_End       = 0.0
   m%F_D_End       = 0.0
   m%F_B_End       = 0.0
   m%F_IMG_End     = 0.0
   p%DP_Const_End  = 0.0
   p%I_MG_End      = 0.0
   p%Mass_MG_End   = 0.0
   p%F_WMG_End     = 0.0
   p%AM_End        = 0.0
   
   m%V_rel_n        = 0.0_ReKi
   m%V_rel_n_HiPass = 0.0_ReKi

   p%WaveVel    => InitInp%WaveVel      
   p%WaveAcc    => InitInp%WaveAcc
   p%WaveDynP   => InitInp%WaveDynP    
   p%WaveTime   => InitInp%WaveTime   
   p%PWaveVel0  => InitInp%PWaveVel0      
   p%PWaveAcc0  => InitInp%PWaveAcc0
   p%PWaveDynP0 => InitInp%PWaveDynP0
   
   p%WaveAccMCF => InitInp%WaveAccMCF
   p%PWaveAccMCF0 => InitInp%PWaveAccMCF0
   

   


END SUBROUTINE AllocateNodeLoadVariables

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is similar to InterpWrappedStpReal, except it returns only the slope for the interpolation.
!! By returning the slope based on Time, we don't have to calculate this for every variable (Yary) we want to interpolate.
!! NOTE: p%WaveTime (and most arrays here) start with index of 0 instead of 1, so we will subtract 1 from "normal" interpolation
!! schemes.
FUNCTION GetInterpolationSlope(Time, p, m, IntWrapIndx) RESULT( InterpSlope )
      REAL(DbKi),                        INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(Morison_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(Morison_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables            
      INTEGER, OPTIONAL,                 INTENT(  OUT)  :: IntWrapIndx

      REAL(SiKi)                                        :: Time_SiKi
      REAL(SiKi)                                        :: TimeMod
      REAL(ReKi)                                        :: InterpSlope

      Time_SiKi = REAL(Time, SiKi)
      TimeMod = MOD(Time_SiKi, p%WaveTime(p%NStepWave)) !p%WaveTime starts at index 0, so it has p%NStepWave+1 elements
      IF ( TimeMod <= p%WaveTime(1) )  THEN !second element
         m%LastIndWave = 0
      END IF
      
      IF ( TimeMod <= p%WaveTime(0) )  THEN
         m%LastIndWave = 0
         InterpSlope = 0.0_ReKi  ! returns values at m%LastIndWave
         IF(PRESENT(IntWrapIndx)) IntWrapIndx = 0
      ELSE IF ( TimeMod >= p%WaveTime(p%NStepWave) )  THEN
         m%LastIndWave = p%NStepWave-1
         InterpSlope = 1.0_ReKi  ! returns values at p%NStepWave
         IF(PRESENT(IntWrapIndx)) IntWrapIndx = p%NStepWave
      ELSE
         m%LastIndWave = MAX( MIN( m%LastIndWave, p%NStepWave-1 ), 0 )

         DO

            IF ( TimeMod < p%WaveTime(m%LastIndWave) )  THEN

               m%LastIndWave = m%LastIndWave - 1

            ELSE IF ( TimeMod >= p%WaveTime(m%LastIndWave+1) )  THEN

               m%LastIndWave = m%LastIndWave + 1

            ELSE
               IF(PRESENT(IntWrapIndx)) IntWrapIndx = m%LastIndWave
               
               InterpSlope = ( TimeMod - p%WaveTime(m%LastIndWave) )/( p%WaveTime(m%LastIndWave+1) - p%WaveTime(m%LastIndWave) )
               RETURN ! stop checking DO loop
            END IF

         END DO
   
      END IF
      
END FUNCTION GetInterpolationSlope
!----------------------------------------------------------------------------------------------------------------------------------
!> Use in conjunction with GetInterpolationSlope, to replace InterpWrappedStpReal here.
FUNCTION InterpolateWithSlope(InterpSlope, Ind, YAry)
      REAL(ReKi), INTENT(IN)                            :: InterpSlope
      INTEGER(IntKi), INTENT(IN )                       :: Ind           !< Misc/optimization variables
      REAL(SiKi), INTENT(IN)                            :: YAry(0:)
      REAL(ReKi)                                        :: InterpolateWithSlope

      InterpolateWithSlope = ( YAry(Ind+1) - YAry(Ind) )*InterpSlope + YAry(Ind)

END FUNCTION InterpolateWithSlope
!----------------------------------------------------------------------------------------------------------------------------------
!> Use in conjunction with GetInterpolationSlope, to replace InterpWrappedStpReal here.
FUNCTION InterpolateWithSlopeR(InterpSlope, Ind, YAry)
      REAL(ReKi), INTENT(IN)                            :: InterpSlope
      INTEGER(IntKi), INTENT(IN )                       :: Ind           !< Misc/optimization variables
      REAL(ReKi), INTENT(IN)                            :: YAry(0:)
      REAL(ReKi)                                        :: InterpolateWithSlopeR

      InterpolateWithSlopeR = ( YAry(Ind+1) - YAry(Ind) )*InterpSlope + YAry(Ind)

END FUNCTION InterpolateWithSlopeR
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
   INTEGER                                           :: I, J, K
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
   REAL(ReKi)               :: CMatrix(3,3), CTrans(3,3) ! Direction cosine matrix for element, and its transpose
   REAL(ReKi)               :: z1
   REAL(ReKi)               :: z2
   REAL(ReKi)               :: r1
   REAL(ReKi)               :: r2
   REAL(ReKi)               :: r1b
   REAL(ReKi)               :: r2b
   REAL(ReKi)               :: rMidb
   REAL(ReKi)               :: dRdl_mg     ! shorthand for taper including marine growth of element i
   REAL(ReKi)               :: dRdl_mg_b   ! shorthand for taper including marine growth of element i with radius scaling by sqrt(Cb)
   REAL(ReKi)               :: RMGFSInt    ! Member radius with marine growth at the intersection with the instantaneous free surface
   REAL(ReKi)               :: g     ! gravity constant
   REAL(ReKi)               :: k_hat(3), k_hat1(3), k_hat2(3) ! Elemental unit vector pointing from 1st node to 2nd node of the element
   REAL(ReKi)               :: n_hat(3)
   REAL(ReKi)               :: alpha ! final load distribution factor for element
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
   REAL(ReKi)               :: pos1(3), pos2(3), positionXY(2)   
   REAL(ReKi)               :: Imat(3,3)
   REAL(ReKi)               :: iArm(3), iTerm(3), Ioffset, h_c, dRdl_p, dRdl_pp, f_hydro(3), Am(3,3), lstar, deltal, deltalLeft, deltalRight
   REAL(ReKi)               :: h, h_c_AM, deltal_AM
   REAL(ReKi)               :: F_WMG(6), F_IMG(6), F_If(6), F_B0(6), F_B1(6), F_B2(6), F_B_End(6)

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
   REAL(ReKi)               :: FDynPFSInt
   REAL(ReKi)               :: vrelFSInt(3)
   REAL(ReKi)               :: pos1Prime(3)
   REAL(ReKi)               :: WtrDpth
   REAL(ReKi)               :: FAMCFFSInt(3)
   INTEGER(IntKi)           :: MemSubStat, NumFSX
   REAL(DbKi)               :: theta1, theta2
   REAL(ReKi)               :: y_hat(3), z_hat(3), posMid(3), zetaMid, FSPt(3)
   INTEGER(IntKi)           :: secStat

   LOGICAL                  :: Is1stElement

   ! Initialize errStat
   errStat = ErrID_None         
   errMsg  = ""               
   Imat    = 0.0_ReKi   
   g       = p%Gravity
   WtrDpth = p%WtrDpth + p%MSL2SWL ! Water depth measured from the free surface
   
   !InterpolationSlope = GetInterpolationSlope(Time, p, m, IntWrapIndx)

   !===============================================================================================
   ! Calculate the fluid kinematics at all mesh nodes and store for use in the equations below

   DO j = 1, p%NNodes
      IF (p%WaveDisp == 0 ) THEN
         ! use the initial X,Y location
         pos1(1) = u%Mesh%Position(1,j)
         pos1(2) = u%Mesh%Position(2,j)
      ELSE
         ! Use current X,Y location
         pos1(1) = u%Mesh%TranslationDisp(1,j) + u%Mesh%Position(1,j)
         pos1(2) = u%Mesh%TranslationDisp(2,j) + u%Mesh%Position(2,j)
      END IF
      
      IF (p%WaveStMod > 0 .AND. p%WaveDisp /= 0) THEN ! Wave stretching enabled
        pos1(3) = u%Mesh%Position(3,j) + u%Mesh%TranslationDisp(3,j) - p%MSL2SWL  ! Use the current Z location.
      ELSE ! Wave stretching disabled
        pos1(3) = u%Mesh%Position(3,j) - p%MSL2SWL  ! We are intentionally using the undisplaced Z position of the node.
      END IF
            
      ! Compute the free surface elevation at the x/y position of all nodes
      positionXY = (/pos1(1),pos1(2)/)      
      m%WaveElev1(j) = SeaSt_Interp_3D( Time, positionXY, p%WaveElev1, p%seast_interp_p, m%SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (associated(p%WaveElev2)) THEN
        m%WaveElev2(j) = SeaSt_Interp_3D( Time, positionXY, p%WaveElev2, p%seast_interp_p, m%SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        m%WaveElev(j) =  m%WaveElev1(j) + m%WaveElev2(j)
      ELSE
        m%WaveElev(j) =  m%WaveElev1(j) 
      END IF      
      
      IF (p%WaveStMod == 0) THEN ! No wave stretching
    
          IF ( pos1(3) <= 0.0_ReKi) THEN ! Node is at or below the SWL
              ! Use location to obtain interpolated values of kinematics         
              call SeaSt_Interp_Setup( Time, pos1, p%seast_interp_p, m%seast_interp_m, ErrStat2, ErrMsg2 ) 
                call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              m%FV(:,j)  = SeaSt_Interp_4D_Vec( p%WaveVel, m%seast_interp_m, ErrStat2, ErrMsg2 )
                call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              m%FA(:,j) = SeaSt_Interp_4D_Vec( p%WaveAcc, m%seast_interp_m, ErrStat2, ErrMsg2 )
                call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              m%FDynP(j)  = SeaSt_Interp_4D    ( p%WaveDynP, m%seast_interp_m, ErrStat2, ErrMsg2 )
                call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              m%vrel(:,j) = m%FV(:,j) - u%Mesh%TranslationVel(:,j)
              m%nodeInWater(j) = 1_IntKi
          ELSE ! Node is above the SWL
              m%FV(:,j)  = 0.0
              m%FA(:,j)  = 0.0
              m%FDynP(j) = 0.0
              m%vrel(:,j) = 0.0  
              m%nodeInWater(j) = 0_IntKi
          END IF
      
      ELSE ! Wave stretching enabled
      
          IF ( pos1(3) <= m%WaveElev(j)) THEN ! Node is submerged
          
              m%nodeInWater(j) = 1_IntKi
 
              IF (p%WaveStMod <3) THEN ! Vertical or extrapolated wave stretching
          
                  IF ( pos1(3) <= 0.0_SiKi) THEN ! Node is below the SWL - evaluate wave dynamics as usual
          
                      ! Use location to obtain interpolated values of kinematics         
                      call SeaSt_Interp_Setup( Time, pos1, p%seast_interp_p, m%seast_interp_m, ErrStat2, ErrMsg2 ) 
                        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                      m%FV(:,j)  = SeaSt_Interp_4D_Vec( p%WaveVel, m%seast_interp_m, ErrStat2, ErrMsg2 )
                        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                      m%FA(:,j) = SeaSt_Interp_4D_Vec( p%WaveAcc, m%seast_interp_m, ErrStat2, ErrMsg2 )
                        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                      m%FDynP(j)  = SeaSt_Interp_4D    ( p%WaveDynP, m%seast_interp_m, ErrStat2, ErrMsg2 )
                        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
          
                  ELSE ! Node is above SWL - need wave stretching
          
                      ! Vertical wave stretching
                      m%FV(:,j)  = SeaSt_Interp_3D_vec( Time, positionXY, p%WaveVel0, p%seast_interp_p,  m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
                        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                      m%FA(:,j)  = SeaSt_Interp_3D_vec( Time, positionXY, p%WaveAcc0, p%seast_interp_p,  m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
                        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                      m%FDynP(j) = SeaSt_Interp_3D    ( Time, positionXY, p%WaveDynP0, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
                        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                      
                      ! Extrapoled wave stretching
                      IF (p%WaveStMod == 2) THEN 
                        m%FV(:,j)  = m%FV(:,j)  + SeaSt_Interp_3D_vec( Time, positionXY, p%PWaveVel0,  p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos1(3)
                          call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        m%FA(:,j)  = m%FA(:,j)  + SeaSt_Interp_3D_vec( Time, positionXY, p%PWaveAcc0,  p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos1(3)
                          call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        m%FDynP(j) = m%FDynP(j) + SeaSt_Interp_3D    ( Time, positionXY, p%PWaveDynP0, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos1(3)
                          call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                      END IF
          
                  END IF ! Node is submerged
 
              ELSE ! Wheeler stretching - no need to check whether the node is above or below SWL
                  
                  ! Map the node z-position linearly from [-WtrDpth,m%WaveElev(j)] to [-WtrDpth,0] 
                  pos1Prime = pos1
                  pos1Prime(3) = WtrDpth*(WtrDpth+pos1(3))/(WtrDpth+m%WaveElev(j))-WtrDpth
                  
                  ! Obtain the wave-field variables by interpolation with the mapped position.
                  call SeaSt_Interp_Setup( Time, pos1Prime, p%seast_interp_p, m%seast_interp_m, ErrStat2, ErrMsg2 ) 
                    call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  m%FV(:,j)  = SeaSt_Interp_4D_Vec( p%WaveVel, m%seast_interp_m, ErrStat2, ErrMsg2 )
                    call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  m%FA(:,j) = SeaSt_Interp_4D_Vec( p%WaveAcc, m%seast_interp_m, ErrStat2, ErrMsg2 )
                    call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  m%FDynP(j)  = SeaSt_Interp_4D    ( p%WaveDynP, m%seast_interp_m, ErrStat2, ErrMsg2 )
                    call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              
              END IF
          
              m%vrel(:,j) = m%FV(:,j) - u%Mesh%TranslationVel(:,j)
        
          ELSE ! Node is out of water - zero-out all wave dynamics
          
              m%nodeInWater(j) = 0_IntKi  
              m%FV(:,j)  = 0.0
              m%FA(:,j)  = 0.0
              m%FDynP(j) = 0.0
              m%vrel(:,j) = 0.0 
          
          END IF ! If node is in or out of water
      
      END IF ! If wave stretching is on or off

   END DO ! j = 1, p%NNodes
   
   ! Scaled fluid acceleration for the MacCamy-Fuchs model
   IF ( ASSOCIATED(p%WaveAccMCF) ) THEN
      DO im = 1,p%NMembers
         IF ( p%Members(im)%PropMCF .AND. ( .NOT. p%Members(im)%PropPot ) ) THEN
            DO i = 1,p%Members(im)%NElements+1
               j = p%Members(im)%NodeIndx(i)
               
               IF (p%WaveDisp == 0 ) THEN
                  ! use the initial X,Y location
                  pos1(1) = u%Mesh%Position(1,j)
                  pos1(2) = u%Mesh%Position(2,j)
               ELSE
                  ! Use current X,Y location
                  pos1(1) = u%Mesh%TranslationDisp(1,j) + u%Mesh%Position(1,j)
                  pos1(2) = u%Mesh%TranslationDisp(2,j) + u%Mesh%Position(2,j)
               END IF
      
               IF (p%WaveStMod > 0 .AND. p%WaveDisp /= 0) THEN ! Wave stretching enabled
                  pos1(3) = u%Mesh%Position(3,j) + u%Mesh%TranslationDisp(3,j) - p%MSL2SWL  ! Use the current Z location.
               ELSE ! Wave stretching disabled
                  pos1(3) = u%Mesh%Position(3,j) - p%MSL2SWL  ! We are intentionally using the undisplaced Z position of the node.
               END IF
            
               ! Compute the free surface elevation at the x/y position of all nodes
               positionXY = (/pos1(1),pos1(2)/)
               
               IF (p%WaveStMod == 0) THEN ! No wave stretching
    
                  IF ( pos1(3) <= 0.0_ReKi) THEN ! Node is at or below the SWL
                     ! Use location to obtain interpolated values of kinematics         
                     call SeaSt_Interp_Setup( Time, pos1, p%seast_interp_p, m%seast_interp_m, ErrStat2, ErrMsg2 ) 
                        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     m%FAMCF(:,j) = SeaSt_Interp_4D_Vec( p%WaveAccMCF, m%seast_interp_m, ErrStat2, ErrMsg2 )
                        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  ELSE ! Node is above the SWL
                     m%FAMCF(:,j)  = 0.0
                  END IF
      
               ELSE ! Wave stretching enabled
      
                  IF ( pos1(3) <= m%WaveElev(j)) THEN ! Node is submerged
      
                     IF (p%WaveStMod <3) THEN ! Vertical or extrapolated wave stretching

                        IF ( pos1(3) <= 0.0_SiKi) THEN ! Node is below the SWL - evaluate wave dynamics as usual
                           ! Use location to obtain interpolated values of kinematics         
                           call SeaSt_Interp_Setup( Time, pos1, p%seast_interp_p, m%seast_interp_m, ErrStat2, ErrMsg2 ) 
                              call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                           m%FAMCF(:,j) = SeaSt_Interp_4D_Vec( p%WaveAccMCF, m%seast_interp_m, ErrStat2, ErrMsg2 )
                              call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        ELSE ! Node is above SWL - need wave stretching
                           
                           
                           ! Vertical wave stretching
                           m%FAMCF(:,j)  = SeaSt_Interp_3D_vec( Time, positionXY, p%WaveAccMCF0, p%seast_interp_p, m%SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
                              call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                           
                           ! Extrapoled wave stretching
                           IF (p%WaveStMod == 2) THEN 
                              m%FAMCF(:,j)  = m%FAMCF(:,j)  + SeaSt_Interp_3D_vec( Time, positionXY, p%PWaveAccMCF0,  p%seast_interp_p, m%SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos1(3)
                                 call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                           END IF
          
                        END IF ! Node is submerged
 
                     ELSE ! Wheeler stretching - no need to check whether the node is above or below SWL
                  
                        ! Map the node z-position linearly from [-WtrDpth,m%WaveElev(j)] to [-WtrDpth,0] 
                        pos1Prime = pos1
                        pos1Prime(3) = WtrDpth*(WtrDpth+pos1(3))/(WtrDpth+m%WaveElev(j))-WtrDpth
                  
                        ! Obtain the wave-field variables by interpolation with the mapped position.
                        call SeaSt_Interp_Setup( Time, pos1Prime, p%seast_interp_p, m%seast_interp_m, ErrStat2, ErrMsg2 ) 
                           call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                        m%FAMCF(:,j) = SeaSt_Interp_4D_Vec( p%WaveAccMCF, m%seast_interp_m, ErrStat2, ErrMsg2 )
                           call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              
                    END IF
        
                  ELSE ! Node is out of water - zero-out all wave dynamics
          
                     m%FAMCF(:,j)  = 0.0
          
                  END IF ! If node is in or out of water
      
               END IF ! If wave stretching is on or off

            END DO
         END IF
      END DO
   END IF

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
      N   = p%Members(im)%NElements      
      mem = p%Members(im)   !@mhall: does this have much overhead?

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
      IF ( p%WaveStMod .EQ. 0_IntKi ) THEN ! No wave stretching - Only need to check the two ends
         IF ( m%nodeInWater(mem%NodeIndx(1)) .NE. m%nodeInWater(mem%NodeIndx(N+1)) ) THEN
            MemSubStat = 1_IntKi  ! Member centerline crosses the SWL once
         ELSE IF ( m%nodeInWater(mem%NodeIndx(1)) .EQ. 0_IntKi ) THEN
            MemSubStat = 3_IntKi  ! Member centerline completely above water
         ELSE
            MemSubStat = 0_IntKi  ! Member centerline fully submerged
         END IF 
      ELSE IF ( p%WaveStMod > 0_IntKi ) THEN ! Has wave stretching - Need to check every node
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
            ! the first and last NodeIndx values point to the corresponding Joint nodes idices which are at the start of the Mesh

            pos1    = u%Mesh%TranslationDisp(:, mem%NodeIndx(i))   + u%Mesh%Position(:, mem%NodeIndx(i)) 
            pos1(3) = pos1(3) - p%MSL2SWL
            pos2    = u%Mesh%TranslationDisp(:, mem%NodeIndx(i+1)) + u%Mesh%Position(:, mem%NodeIndx(i+1)) 
            pos2(3) = pos2(3) - p%MSL2SWL

            call GetOrientationAngles( pos1, pos2, phi, sinPhi, cosPhi, tanPhi, sinBeta, cosBeta, k_hat, errStat2, errMsg2 )
              call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            call Morison_DirCosMtrx( pos1, pos2, CMatrix )
            CTrans  = transpose(CMatrix)
            ! save some commonly used variables   
            dl        = mem%dl
            z1        = pos1(3)          ! get node z locations from input mesh
            z2        = pos2(3)
            r1        = mem%RMG(i  )     ! outer radius at element nodes including marine growth
            r2        = mem%RMG(i+1)
            r1b       = mem%RMGB(i  )    ! outer radius at element nodes including marine growth scaled by sqrt(Cb)
            r2b       = mem%RMGB(i+1)
            dRdl_mg   = mem%dRdl_mg(i)   ! Taper of element including marine growth
            dRdl_mg_b = mem%dRdl_mg_b(i) ! Taper of element including marine growth with radius scaling by sqrt(Cb)
            a_s1      = u%Mesh%TranslationAcc(:, mem%NodeIndx(i  ))
            alpha_s1  = u%Mesh%RotationAcc   (:, mem%NodeIndx(i  ))
            omega_s1  = u%Mesh%RotationVel   (:, mem%NodeIndx(i  ))
            a_s2      = u%Mesh%TranslationAcc(:, mem%NodeIndx(i+1))
            alpha_s2  = u%Mesh%RotationAcc   (:, mem%NodeIndx(i+1))
            omega_s2  = u%Mesh%RotationVel   (:, mem%NodeIndx(i+1))
              
            ! ------------------ marine growth: Sides: Section 4.1.2 --------------------  
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
            
            ! lower node
            Ioffset   = mem%h_cmg_l(i)*mem%h_cmg_l(i)*mem%m_mg_l(i)
            Imat(1,1) = mem%I_rmg_l(i) - Ioffset
            Imat(2,2) = mem%I_rmg_l(i) - Ioffset
            Imat(3,3) = mem%I_lmg_l(i) - Ioffset
            Imat      =  matmul(matmul(CMatrix, Imat), CTrans)
            iArm = mem%h_cmg_l(i) * k_hat
            iTerm     = ( -a_s1 - cross_product(omega_s1, cross_product(omega_s1,iArm )) - cross_product(alpha_s1,iArm) ) * mem%m_mg_l(i)
            F_IMG(1:3) = iTerm
            F_IMG(4:6) = - cross_product(a_s1 * mem%m_mg_l(i), mem%h_cmg_l(i) * k_hat) + matmul(Imat, alpha_s1)  &
                         - cross_product(omega_s1,matmul(Imat,omega_s1))
            m%memberLoads(im)%F_IMG(:,i) = m%memberLoads(im)%F_IMG(:,i) + F_IMG
            y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + F_IMG(1:3)
            y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + F_IMG(4:6)
            
            ! upper node
            Ioffset   = mem%h_cmg_u(i)*mem%h_cmg_u(i)*mem%m_mg_u(i)
            Imat(1,1) = mem%I_rmg_u(i) - Ioffset
            Imat(2,2) = mem%I_rmg_u(i) - Ioffset
            Imat(3,3) = mem%I_lmg_u(i) - Ioffset
            Imat      =  matmul(matmul(CMatrix, Imat), CTrans)
            iArm = mem%h_cmg_u(i) * k_hat
            iTerm     = ( -a_s2 - cross_product(omega_s2, cross_product(omega_s2,iArm )) - cross_product(alpha_s2,iArm) ) * mem%m_mg_u(i)
            F_IMG(1:3) = iTerm
            F_IMG(4:6) = - cross_product(a_s2 * mem%m_mg_u(i), mem%h_cmg_u(i) * k_hat) + matmul(Imat, alpha_s2) &
                         - cross_product(omega_s2,matmul(Imat,omega_s2))
            m%memberLoads(im)%F_IMG(:,i+1) = m%memberLoads(im)%F_IMG(:,i+1) + F_IMG
            y%Mesh%Force (:,mem%NodeIndx(i+1)) = y%Mesh%Force (:,mem%NodeIndx(i+1)) + F_IMG(1:3)
            y%Mesh%Moment(:,mem%NodeIndx(i+1)) = y%Mesh%Moment(:,mem%NodeIndx(i+1)) + F_IMG(4:6)

            ! ------------------- buoyancy loads: sides: Sections 3.1 and 3.2 ------------------------
            IF (mem%MHstLMod == 1) THEN
               IF ( p%WaveStMod > 0_IntKi ) THEN ! If wave stretching is enabled, compute buoyancy up to free surface
                  CALL GetTotalWaveElev( Time, pos1, Zeta1, ErrStat2, ErrMsg2 )
                  CALL GetTotalWaveElev( Time, pos2, Zeta2, ErrStat2, ErrMsg2 )
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               ELSE ! Without wave stretching, compute buoyancy based on SWL
                  Zeta1 = 0.0_ReKi
                  Zeta2 = 0.0_ReKi
               END IF
               Is1stElement = ( i .EQ. 1)
               CALL getElementHstLds_Mod1( Time, pos1, pos2, Zeta1, Zeta2, k_hat, r1b, r2b, dl, mem%alpha(i), Is1stElement, F_B0, F_B1, F_B2, ErrStat2, ErrMsg2 )
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
               rMidb  = 0.5 * (r1b +r2b )
               IF (p%WaveStMod > 0) THEN
                  CALL GetTotalWaveElev( Time, posMid, ZetaMid, ErrStat2, ErrMsg2 )
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  CALL GetFreeSurfaceNormal( Time, posMid, rMidb, n_hat, ErrStat2, ErrMsg2 )
                    CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  FSPt = (/posMid(1),posMid(2),ZetaMid/) ! Reference point on the free surface
               ELSE
                  FSPt = (/posMid(1),posMid(2),0.0/)
                  n_hat = (/0.0,0.0,1.0/)
               END IF    
               CALL GetSectionUnitVectors( k_hat, y_hat, z_hat )
               CALL getElementHstLds_Mod2( pos1, pos2, FSPt, k_hat, y_hat, z_hat, n_hat, r1b, r2b, dl, F_B1, F_B2, ErrStat2, ErrMsg2)
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
      DO i = max(mem%i_floor,1), N    ! loop through member elements that are not completely buried in the seabed
         
         ! calculate instantaneous incline angle and heading, and related trig values
         ! the first and last NodeIndx values point to the corresponding Joint nodes idices which are at the start of the Mesh

         pos1    = u%Mesh%TranslationDisp(:, mem%NodeIndx(i))   + u%Mesh%Position(:, mem%NodeIndx(i)) 
         pos1(3) = pos1(3) - p%MSL2SWL
         pos2    = u%Mesh%TranslationDisp(:, mem%NodeIndx(i+1)) + u%Mesh%Position(:, mem%NodeIndx(i+1)) 
         pos2(3) = pos2(3) - p%MSL2SWL

         call GetOrientationAngles( pos1, pos2, phi, sinPhi, cosPhi, tanPhi, sinBeta, cosBeta, k_hat, errStat2, errMsg2 )
           call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         call Morison_DirCosMtrx( pos1, pos2, CMatrix )
         CTrans  = transpose(CMatrix)
         ! save some commonly used variables   
         dl        = mem%dl
         z1        = pos1(3)          ! get node z locations from input mesh
         z2        = pos2(3)
         r1        = mem%RMG(i  )     ! outer radius at element nodes including marine growth
         r2        = mem%RMG(i+1)
         r1b       = mem%RMGB(i  )    ! outer radius at element nodes including marine growth scaled by sqrt(Cb)
         r2b       = mem%RMGB(i+1)
         dRdl_mg   = mem%dRdl_mg(i)   ! Taper of element including marine growth
         dRdl_mg_b = mem%dRdl_mg_b(i) ! Taper of element including marine growth with radius scaling by sqrt(Cb)
         a_s1      = u%Mesh%TranslationAcc(:, mem%NodeIndx(i  ))
         alpha_s1  = u%Mesh%RotationAcc   (:, mem%NodeIndx(i  ))
         omega_s1  = u%Mesh%RotationVel   (:, mem%NodeIndx(i  ))
         a_s2      = u%Mesh%TranslationAcc(:, mem%NodeIndx(i+1))
         alpha_s2  = u%Mesh%RotationAcc   (:, mem%NodeIndx(i+1))
         omega_s2  = u%Mesh%RotationVel   (:, mem%NodeIndx(i+1))

         ! ------------------ flooded ballast inertia: sides: Section 6.1.1 : Always compute regardless of PropPot setting ---------------------
         ! lower node
         Ioffset   = mem%h_cfb_l(i)*mem%h_cfb_l(i)*mem%m_fb_l(i)
         Imat(1,1) = mem%I_rfb_l(i) - Ioffset
         Imat(2,2) = mem%I_rfb_l(i) - Ioffset
         Imat(3,3) = mem%I_lfb_l(i) - Ioffset
         iArm = mem%h_cfb_l(i) * k_hat
         iTerm     = ( -a_s1  - cross_product(omega_s1, cross_product(omega_s1,iArm ))  -  cross_product(alpha_s1,iArm) ) * mem%m_fb_l(i)
         F_If(1:3) =  iTerm
         F_If(4:6) =  - cross_product(a_s1 * mem%m_fb_l(i), mem%h_cfb_l(i) * k_hat) + matmul(Imat, alpha_s1) &
                      - cross_product(omega_s1,matmul(Imat,omega_s1)) 
         m%memberLoads(im)%F_If(:,i) = m%memberLoads(im)%F_If(:,i) + F_If
         y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + F_If(1:3)
         y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + F_If(4:6)
         
         ! upper node
         Ioffset   = mem%h_cfb_u(i)*mem%h_cfb_u(i)*mem%m_fb_u(i)
         Imat(1,1) = mem%I_rfb_u(i) - Ioffset
         Imat(2,2) = mem%I_rfb_u(i) - Ioffset
         Imat(3,3) = mem%I_lfb_u(i) - Ioffset
         iArm = mem%h_cfb_u(i) * k_hat
         iTerm     = ( -a_s2  - cross_product(omega_s2, cross_product(omega_s2,iArm ))  -  cross_product(alpha_s2,iArm) ) * mem%m_fb_u(i)
         F_If(1:3) = iTerm
         F_If(4:6) = - cross_product(a_s2 * mem%m_fb_u(i), mem%h_cfb_u(i) * k_hat) + matmul(Imat, alpha_s2) &
                     - cross_product(omega_s2,matmul(Imat,omega_s2)) 
         m%memberLoads(im)%F_If(:,i+1) = m%memberLoads(im)%F_If(:,i+1) + F_If
         y%Mesh%Force (:,mem%NodeIndx(i+1)) = y%Mesh%Force (:,mem%NodeIndx(i+1)) + F_If(1:3)
         y%Mesh%Moment(:,mem%NodeIndx(i+1)) = y%Mesh%Moment(:,mem%NodeIndx(i+1)) + F_If(4:6)  
         
         ! ------------------ flooded ballast weight : sides : Section 5.1.2 & 5.2.2  : Always compute regardless of PropPot setting ---------------------
         
         ! NOTE: For memfloodstatus and floodstatus: 0 = fully buried or not ballasted, 1 = fully flooded, 2 = partially flooded
         
         ! fully filled elements
         if (mem%floodstatus(i) == 1) then
            
            ! Compute lstar
            if ( mem%memfloodstatus == 2) then  
               ! partially flooded MEMBER
               lstar = dl*(i-1) - mem%l_fill
            elseif (cosPhi >= 0.0 ) then
               lstar = dl*(i-N-1) 
            else
               lstar = dl*(i-1)
            end if
            Fl =TwoPi * mem%dRdl_in(i) * mem%FillDens * p%gravity * dl *( -( mem%Rin(i) + 0.5* mem%dRdl_in(i)*dl )*mem%z_overfill +  &
                        ( lstar*mem%Rin(i) + 0.5*(lstar*mem%dRdl_in(i) + mem%Rin(i) )*dl + mem%dRdl_in(i)*dl**2/3.0 )*cosphi )

            ! forces and moment in tilted coordinates about node i
            Fr = mem%Cfr_fb(i)*sinPhi     
            Moment  = mem%CM0_fb(i)*sinPhi - Fr*mem%alpha_fb_star(i)*dl
           
            ! calculate full vector and distribute to nodes
            call DistributeElementLoads(Fl, Fr, Moment, sinPhi, cosPhi, sinBeta, cosBeta, (1-mem%alpha_fb_star(i)), F_B1, F_B2)
            m%memberLoads(im)%F_BF(:, i)   = m%memberLoads(im)%F_BF(:, i) + F_B2  ! 1-alpha
            m%memberLoads(im)%F_BF(:, i+1) = m%memberLoads(im)%F_BF(:, i+1) + F_B1 ! alpha
            y%Mesh%Force (:,mem%NodeIndx(i  )) = y%Mesh%Force (:,mem%NodeIndx(i  )) + F_B2(1:3)
            y%Mesh%Moment(:,mem%NodeIndx(i  )) = y%Mesh%Moment(:,mem%NodeIndx(i  )) + F_B2(4:6)
            y%Mesh%Force (:,mem%NodeIndx(i+1)) = y%Mesh%Force (:,mem%NodeIndx(i+1)) + F_B1(1:3)
            y%Mesh%Moment(:,mem%NodeIndx(i+1)) = y%Mesh%Moment(:,mem%NodeIndx(i+1)) + F_B1(4:6)
           
         ! partially filled element
         else if (mem%floodstatus(i) == 2) then
           
            ! forces and moment in tilted coordinates about node i
            Fl = mem%Cfl_fb(i)*cosPhi     
            Fr = mem%Cfr_fb(i)*sinPhi     
            Moment  = mem%CM0_fb(i)*sinPhi + Fr*(1 - mem%alpha_fb_star(i))*dl
        
            ! calculate full vector and distribute to nodes
            call DistributeElementLoads(Fl, Fr, Moment, sinPhi, cosPhi, sinBeta, cosBeta, mem%alpha_fb_star(i), F_B1, F_B2)
            m%memberLoads(im)%F_BF(:, i) = m%memberLoads(im)%F_BF(:, i) + F_B1     ! alpha
            m%memberLoads(im)%F_BF(:, i-1) = m%memberLoads(im)%F_BF(:, i-1) + F_B2 ! 1- alpha
            y%Mesh%Force (:,mem%NodeIndx(i  )) = y%Mesh%Force (:,mem%NodeIndx(i  )) + F_B1(1:3)
            y%Mesh%Moment(:,mem%NodeIndx(i  )) = y%Mesh%Moment(:,mem%NodeIndx(i  )) + F_B1(4:6)
            y%Mesh%Force (:,mem%NodeIndx(i-1)) = y%Mesh%Force (:,mem%NodeIndx(i-1)) + F_B2(1:3)
            y%Mesh%Moment(:,mem%NodeIndx(i-1)) = y%Mesh%Moment(:,mem%NodeIndx(i-1)) + F_B2(4:6)    
        
         ! no load for unflooded element or element fully below seabed
        
         end if
   
      END DO ! i = max(mem%i_floor,1), N    ! loop through member elements that are not fully buried in the seabed   

      !-----------------------------------------------------------------------------------------------------!
      !                               External Hydrodynamic Side Loads - Start                              !
      !-----------------------------------------------------------------------------------------------------!
      IF ( p%WaveStMod > 0 .AND. MemSubStat == 1 .AND. (m%NodeInWater(mem%NodeIndx(N+1)).EQ.0_IntKi) ) THEN 
      !----------------------------Apply load smoothing----------------------------!
      ! only when:
      ! 1. wave stretching is enabled
      ! 2. member centerline crosses the free surface exactly once
      ! 3. the last node is out of water, which implies the first node is in water
                
        FSElem = -1 ! Initialize the No. of the partially wetted element as -1
      
        DO i = mem%i_floor+1,N ! loop through member nodes starting from the first node above seabed, but skip the last node which should not be submerged anyways
           
           ! Get positions of node i and i+1
           IF (p%WaveDisp /= 0) THEN ! Use current position
              pos1    = u%Mesh%TranslationDisp(:, mem%NodeIndx(i))   + u%Mesh%Position(:, mem%NodeIndx(i))  
              pos2    = u%Mesh%TranslationDisp(:, mem%NodeIndx(i+1)) + u%Mesh%Position(:, mem%NodeIndx(i+1))
           ELSE ! Use initial position
              pos1    = u%Mesh%Position(:, mem%NodeIndx(i))  
              pos2    = u%Mesh%Position(:, mem%NodeIndx(i+1))
           END if
           ! We need to subtract the MSL2SWL offset to place this in the SWL reference system
           pos1(3) = pos1(3) - p%MSL2SWL
           pos2(3) = pos2(3) - p%MSL2SWL 
           
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
         
           ! Compute the slope of member radius
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
       
           !-------------------- hydrodynamic drag loads: sides: Section 7.1.2 ------------------------!
           vec = matmul( mem%Ak,m%vrel(:,mem%NodeIndx(i)) )
           f_hydro = mem%Cd(i)*p%WtrDens*mem%RMG(i)*TwoNorm(vec)*vec  +  &
                     0.5*mem%AxCd(i)*p%WtrDens*pi*mem%RMG(i)*dRdl_p * abs(dot_product( mem%k, m%vrel(:,mem%NodeIndx(i)) )) * matmul( mem%kkt, m%vrel(:,mem%NodeIndx(i)) )
           CALL LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%memberLoads(im)%F_D(:, i) )
           y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + m%memberLoads(im)%F_D(1:3, i)
           y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + m%memberLoads(im)%F_D(4:6, i)
           IF (i == FSElem) THEN ! Save the distributed load at the first node below the free surface
             F_D0 = f_hydro
           END IF
           
           IF ( .NOT. mem%PropPot ) THEN
              !-------------------- hydrodynamic added mass loads: sides: Section 7.1.3 ------------------------!
              Am = mem%Ca(i)*p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*mem%Ak + 2.0*mem%AxCa(i)*p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*dRdl_p*mem%kkt
              f_hydro = -matmul( Am, u%Mesh%TranslationAcc(:,mem%NodeIndx(i)) )

              IF ( p%AMMod .EQ. 0_IntKi ) THEN ! Compute added-mass force up to the SWL
                 z1 = u%Mesh%Position(3, mem%NodeIndx(i)) - p%MSL2SWL ! Undisplaced z-position of the current node
                 IF ( z1 > 0.0_ReKi ) THEN ! Node is above SWL undisplaced; zero added-mass force
                    f_hydro = 0.0_ReKi
                    CALL LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%memberLoads(im)%F_A(:, i) )
                 ELSE
                    ! Need to compute deltal_AM and h_c_AM based on the formulation without wave stretching.
                    z2 = u%Mesh%Position(3, mem%NodeIndx(i+1)) - p%MSL2SWL ! Undisplaced z-position of the next node
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
              IF (mem%PropMCF) THEN
                 f_hydro=                     p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)       * matmul( mem%Ak,  m%FAMCF(:,mem%NodeIndx(i)) ) + &
                              2.0*mem%AxCa(i)*p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*dRdl_p * matmul( mem%kkt, m%FA(:,mem%NodeIndx(i)) ) + &
                              2.0*m%FDynP(mem%NodeIndx(i))*mem%AxCp(i)*pi*mem%RMG(i)*dRdl_pp*mem%k 
              ELSE
                 f_hydro=(mem%Ca(i)+mem%Cp(i))*p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)       * matmul( mem%Ak,  m%FA(:,mem%NodeIndx(i)) ) + &
                              2.0*mem%AxCa(i)*p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*dRdl_p * matmul( mem%kkt, m%FA(:,mem%NodeIndx(i)) ) + &
                              2.0*m%FDynP(mem%NodeIndx(i))*mem%AxCp(i)*pi*mem%RMG(i)*dRdl_pp*mem%k 
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
        ! Get wave dynamics at the free surface intersection
        IF (p%WaveStMod <3) THEN ! Vertical or extrapolated stretching
           
           IF ( FSInt(3) <= 0.0_ReKi) THEN ! Intersection is below SWL - evaluate wave dynamics as usual
              
              ! Use location to obtain interpolated values of kinematics
              CALL SeaSt_Interp_Setup( Time, FSInt, p%seast_interp_p, m%seast_interp_m, ErrStat2, ErrMsg2 ) 
                CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              FVFSInt = SeaSt_Interp_4D_Vec( p%WaveVel, m%seast_interp_m, ErrStat2, ErrMsg2 )
                CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              FAFSInt = SeaSt_Interp_4D_Vec( p%WaveAcc, m%seast_interp_m, ErrStat2, ErrMsg2 )
                CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              FDynPFSInt = SeaSt_Interp_4D    ( p%WaveDynP, m%seast_interp_m, ErrStat2, ErrMsg2 )
                CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
           
           ELSE ! Intersection is above SWL - need wave stretching
              
              ! Vertical wave stretching
              FVFSInt    = SeaSt_Interp_3D_vec( Time, FSInt(1:2), p%WaveVel0, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
                CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              FAFSInt    = SeaSt_Interp_3D_vec( Time, FSInt(1:2), p%WaveAcc0, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
                CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              FDynPFSInt = SeaSt_Interp_3D    ( Time, FSInt(1:2), p%WaveDynP0, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
                CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              
              ! Extrapolated wave stretching
              IF (p%WaveStMod == 2) THEN 
                FVFSInt    = FVFSInt    + SeaSt_Interp_3D_vec( Time, FSInt(1:2), p%PWaveVel0,  p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * FSInt(3)
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                FAFSInt    = FAFSInt    + SeaSt_Interp_3D_vec( Time, FSInt(1:2), p%PWaveAcc0,  p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * FSInt(3)
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                FDynPFSInt = FDynPFSInt + SeaSt_Interp_3D    ( Time, FSInt(1:2), p%PWaveDynP0, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * FSInt(3)
                  CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              END IF
              
           END IF
           
        ELSE ! Wheeler stretching
           
           ! Points on the free surface is always mapped back to z=0 of the unstretched wave field
           ! Can evaluate the wave-field variables in the same way as vertical stretching
           FVFSInt = SeaSt_Interp_3D_vec( Time, FSInt(1:2), p%WaveVel0, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
             CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
           FAFSInt = SeaSt_Interp_3D_vec( Time, FSInt(1:2), p%WaveAcc0, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
             CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
           FDynPFSInt = SeaSt_Interp_3D( Time, FSInt(1:2), p%WaveDynP0, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
             CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
           
        END IF


        IF ( mem%PropMCF .AND. ( .NOT. mem%PropPot ) ) THEN
           IF (p%WaveStMod <3) THEN ! Vertical or extrapolated stretching
           
              IF ( FSInt(3) <= 0.0_ReKi) THEN ! Intersection is below SWL - evaluate wave dynamics as usual
              
              ! Use location to obtain interpolated values of kinematics
              CALL SeaSt_Interp_Setup( Time, FSInt, p%seast_interp_p, m%seast_interp_m, ErrStat2, ErrMsg2 ) 
                CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              FAMCFFSInt = SeaSt_Interp_4D_Vec( p%WaveAccMCF, m%seast_interp_m, ErrStat2, ErrMsg2 )
                CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

              ELSE ! Intersection is above SWL - need wave stretching
              
                 ! Vertical wave stretching
                 FAMCFFSInt    = SeaSt_Interp_3D_vec( Time, FSInt(1:2), p%WaveAccMCF0, p%seast_interp_p, m%SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
                   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

                 ! Extrapolated wave stretching
                 IF (p%WaveStMod == 2) THEN 
                    FAMCFFSInt    = FAMCFFSInt    + SeaSt_Interp_3D_vec( Time, FSInt(1:2), p%PWaveAccMCF0,  p%seast_interp_p, m%SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * FSInt(3)
                      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                 END IF
              
              END IF
           
           ELSE ! Wheeler stretching
           
              ! Points on the free surface is always mapped back to z=0 of the unstretched wave field
              ! Can evaluate the wave-field variables in the same way as vertical stretching
              FAMCFFSInt = SeaSt_Interp_3D_vec( Time, FSInt(1:2), p%WaveAccMCF0, p%seast_interp_p, m%SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
                CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
           END IF
        END IF

        ! Viscous drag:
        ! Compute relative velocity at the free surface intersection. 
        ! Linear interpolation between the two nodes of the element is used to estimate velocity of the structure
        vrelFSInt = FVFSInt - ( & 
               SubRatio  * u%Mesh%TranslationVel(:,mem%NodeIndx(FSElem+1)) + &
          (1.0-SubRatio) * u%Mesh%TranslationVel(:,mem%NodeIndx(FSElem  ))   &
        )
        dRdl_p  = abs(mem%dRdl_mg(FSElem))
        RMGFSInt = SubRatio * mem%RMG(FSElem+1) + (1.0-SubRatio) * mem%RMG(FSElem)

        vec = matmul( mem%Ak,vrelFSInt )
        F_DS = mem%Cd(FSElem)*p%WtrDens*RMGFSInt*TwoNorm(vec)*vec  +  &
                  0.5*mem%AxCd(FSElem)*p%WtrDens*pi*RMGFSInt*dRdl_p * & 
                  abs(dot_product( mem%k, vrelFSInt )) * matmul( mem%kkt, vrelFSInt )

        ! Hydrodynamic added mass and inertia loads
        IF ( .NOT. mem%PropPot ) THEN
           
           ! ------------------- hydrodynamic added mass loads: sides: Section 7.1.3 ------------------------
           IF (p%AMMod > 0_IntKi) THEN
              Am =      mem%Ca(FSElem)*p%WtrDens*pi*RMGFSInt*RMGFSInt*mem%Ak + &
                  2.0*mem%AxCa(FSElem)*p%WtrDens*pi*RMGFSInt*RMGFSInt*dRdl_p*mem%kkt
              F_AS = -matmul( Am, &
                         SubRatio  * u%Mesh%TranslationAcc(:,mem%NodeIndx(FSElem+1)) + &
                    (1.0-SubRatio) * u%Mesh%TranslationAcc(:,mem%NodeIndx(FSElem  )) )
           END IF
         
           ! ------------------- hydrodynamic inertia loads: sides: Section 7.1.4 ------------------------
           IF ( mem%PropMCF) THEN
              F_IS=                             p%WtrDens*pi*RMGFSInt*RMGFSInt   * matmul( mem%Ak,  FAMCFFSInt ) + &
                           2.0*mem%AxCa(FSElem)*p%WtrDens*pi*RMGFSInt*RMGFSInt*dRdl_p  * matmul( mem%kkt, FAFSInt ) + &
                           2.0*mem%AxCp(FSElem)          *pi*RMGFSInt                *dRdl_pp * FDynPFSInt*mem%k
           ELSE
              F_IS=(mem%Ca(FSElem)+mem%Cp(FSElem))*p%WtrDens*pi*RMGFSInt*RMGFSInt   * matmul( mem%Ak,  FAFSInt ) + &
                           2.0*mem%AxCa(FSElem)*p%WtrDens*pi*RMGFSInt*RMGFSInt*dRdl_p  * matmul( mem%kkt, FAFSInt ) + &
                           2.0*mem%AxCp(FSElem)          *pi*RMGFSInt                *dRdl_pp * FDynPFSInt*mem%k
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
           ! We need to subtract the MSL2SWL offset to place this in the SWL reference system
           ! Using the initial z-position to be consistent with the evaluation of wave kinematics
           IF (p%WaveStMod > 0 .AND. p%WaveDisp /= 0) THEN
              ! Use current z-position
              z1 = u%Mesh%Position(3, mem%NodeIndx(i)) + u%Mesh%TranslationDisp(3, mem%NodeIndx(i)) - p%MSL2SWL
           ELSE
              ! Use initial z-position
              z1 = u%Mesh%Position(3, mem%NodeIndx(i)) - p%MSL2SWL 
           END IF         
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
                    IF (p%WaveStMod > 0 .AND. p%WaveDisp /= 0) THEN ! Use current z-position
                       z2 = u%Mesh%Position(3, mem%NodeIndx(i-1)) + u%Mesh%TranslationDisp(3, mem%NodeIndx(i-1)) - p%MSL2SWL
                    ELSE ! Use initial z-position
                       z2 = u%Mesh%Position(3, mem%NodeIndx(i-1)) - p%MSL2SWL 
                    END IF
                    IF ( p%WaveStMod > 0_IntKi ) THEN ! Wave stretching enabled
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
                    IF (p%WaveStMod > 0 .AND. p%WaveDisp /= 0) THEN ! Use current z-position
                       z2 = u%Mesh%Position(3, mem%NodeIndx(i+1)) + u%Mesh%TranslationDisp(3, mem%NodeIndx(i+1)) - p%MSL2SWL
                    ELSE ! Use initial z-position
                       z2 = u%Mesh%Position(3, mem%NodeIndx(i+1)) - p%MSL2SWL 
                    END IF
                    IF ( p%WaveStMod > 0_IntKi ) THEN ! Wave stretching enabled
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

           ! Compute the slope of the member radius
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
         
           !--------------------- hydrodynamic drag loads: sides: Section 7.1.2 --------------------------------! 
           vec = matmul( mem%Ak,m%vrel(:,mem%NodeIndx(i)) )
           f_hydro = mem%Cd(i)*p%WtrDens*mem%RMG(i)*TwoNorm(vec)*vec  +  &
                     0.5*mem%AxCd(i)*p%WtrDens*pi*mem%RMG(i)*dRdl_p * abs(dot_product( mem%k, m%vrel(:,mem%NodeIndx(i)) )) * matmul( mem%kkt, m%vrel(:,mem%NodeIndx(i)) )
           CALL LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%memberLoads(im)%F_D(:, i) )
           y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + m%memberLoads(im)%F_D(1:3, i)
           y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + m%memberLoads(im)%F_D(4:6, i)
            
           IF ( .NOT. mem%PropPot ) THEN
              !-------------------- hydrodynamic added mass loads: sides: Section 7.1.3 ------------------------!
              Am = mem%Ca(i)*p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*mem%Ak + 2.0*mem%AxCa(i)*p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*dRdl_p*mem%kkt
              f_hydro = -matmul( Am, u%Mesh%TranslationAcc(:,mem%NodeIndx(i)) )
              IF ( p%AMMod .EQ. 0_IntKi ) THEN ! Always compute added-mass force on nodes below SWL when undisplaced
                 z1 = u%Mesh%Position(3, mem%NodeIndx(i)) - p%MSL2SWL ! Undisplaced z-position of the current node
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
                       z2 = u%Mesh%Position(3, mem%NodeIndx(i+1)) - p%MSL2SWL
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
              IF ( mem%PropMCF ) THEN
                 f_hydro=                     p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)        * matmul( mem%Ak,  m%FAMCF(:,mem%NodeIndx(i)) ) + &
                              2.0*mem%AxCa(i)*p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*dRdl_p * matmul( mem%kkt, m%FA(:,mem%NodeIndx(i)) ) + &
                              2.0*m%FDynP(mem%NodeIndx(i))*mem%AxCp(i)*pi*mem%RMG(i)*dRdl_pp*mem%k 
              ELSE
                 f_hydro=(mem%Ca(i)+mem%Cp(i))*p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)       * matmul( mem%Ak,  m%FA(:,mem%NodeIndx(i)) ) + &
                              2.0*mem%AxCa(i)*p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*dRdl_p * matmul( mem%kkt, m%FA(:,mem%NodeIndx(i)) ) + &
                              2.0*m%FDynP(mem%NodeIndx(i))*mem%AxCp(i)*pi*mem%RMG(i)*dRdl_pp*mem%k 
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
      pos1    = u%Mesh%TranslationDisp(:, mem%NodeIndx(1)) + u%Mesh%Position(:, mem%NodeIndx(1)) 
      pos1(3) = pos1(3) - p%MSL2SWL
      pos2    = u%Mesh%TranslationDisp(:, mem%NodeIndx(2)) + u%Mesh%Position(:, mem%NodeIndx(2)) 
      pos2(3) = pos2(3) - p%MSL2SWL
      z1 = pos1(3)
      
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
         pos1    = u%Mesh%TranslationDisp(:, mem%NodeIndx(N))   + u%Mesh%Position(:, mem%NodeIndx(N))
         pos1(3) = pos1(3) - p%MSL2SWL
         pos2    = u%Mesh%TranslationDisp(:, mem%NodeIndx(N+1)) + u%Mesh%Position(:, mem%NodeIndx(N+1))
         pos2(3) = pos2(3) - p%MSL2SWL
         call GetOrientationAngles( pos1, pos2, phi2, sinPhi2, cosPhi2, tanPhi, sinBeta2, cosBeta2, k_hat2, errStat2, errMsg2 )
           call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      end if
      
      ! We need to subtract the MSL2SWL offset to place this  in the SWL reference system
      pos2    = u%Mesh%TranslationDisp(:, mem%NodeIndx(N+1)) + u%Mesh%Position(:, mem%NodeIndx(N+1))
      pos2(3) = pos2(3) - p%MSL2SWL
      z2 = pos2(3)
      
      !----------------------------------- filled buoyancy loads: starts -----------------------------------!
      !TODO: Do the equations below still work if z1 > z2 ?
      !TODO: Should not have to test seabed crossing in time-marching loop
      if ( mem%i_floor == 0 ) then   ! both ends are above seabed
         !--- Water ballast buoyancy ---
         ! if member is fully flooded
         if (mem%memfloodstatus == 1) then
         !if (mem%z_overfill >= 0) then 
            Fl      = -mem%FillDens * g * pi *mem%Rin(  1)**2* (mem%z_overfill + max(z2-z1, 0.0_ReKi))
            Moment  =  mem%FillDens * g * pi *0.25*mem%Rin(  1)**4*sinPhi
            call AddEndLoad(Fl, Moment, sinPhi1, cosPhi1, sinBeta1, cosBeta1, m%F_BF_End(:, mem%NodeIndx(1)))
            
            Fl      =   mem%FillDens * g * pi *mem%Rin(N+1)**2* (mem%z_overfill + max(z1-z2, 0.0_ReKi))
            Moment  =  -mem%FillDens * g * pi *0.25*mem%Rin(N+1)**4*sinPhi            
            call AddEndLoad(Fl, Moment, sinPhi2, cosPhi2, sinBeta2, cosBeta2, m%F_BF_End(:, mem%NodeIndx(N+1)))
            
         ! if member is partially flooded
         else if (mem%l_fill > 0) then 
            Fl      = -mem%FillDens * g * pi *mem%Rin(1)**2*mem%l_fill*cosPhi
            Moment  =  mem%FillDens * g * pi *0.25*mem%Rin(1)**4*sinPhi
            call AddEndLoad(Fl, Moment, sinPhi1, cosPhi1, sinBeta1, cosBeta1, m%F_BF_End(:, mem%NodeIndx(1)))
         else
            ! no load if member is not flooded at all
         end if
         
      elseif ( mem%i_floor < mem%NElements+1 ) then ! upper node is still above the seabed, but lower node is below seabed
         !if (mem%z_overfill >= 0) then 
         if (mem%memfloodstatus == 1) then
            Fl      =   mem%FillDens * g * pi *mem%Rin(N+1)**2* (mem%z_overfill + max(z1-z2, 0.0_ReKi))
            Moment  =  -mem%FillDens * g * pi *0.25*mem%Rin(N+1)**4*sinPhi            
            call AddEndLoad(Fl, Moment, sinPhi2, cosPhi2, sinBeta2, cosBeta2, m%F_BF_End(:, mem%NodeIndx(N+1)))
         end if
         
      else    
         ! no loads because both end nodes are below seabed
      end if

      !------------------------------------ filled buoyancy loads: ends ------------------------------------!

      ! --- no inertia loads from water ballast modeled on ends

      !---------------------------------- external buoyancy loads: starts ----------------------------------!
      if ( (.not. mem%PropPot) .AND. (mem%MHstLMod /= 0) ) then

         ! Get positions of member end nodes
         pos1    = u%Mesh%TranslationDisp(:, mem%NodeIndx(1  )) + u%Mesh%Position(:, mem%NodeIndx(1  )) 
         pos1(3) = pos1(3) - p%MSL2SWL
         r1      = mem%RMGB(1  )
         pos2    = u%Mesh%TranslationDisp(:, mem%NodeIndx(N+1)) + u%Mesh%Position(:, mem%NodeIndx(N+1))
         pos2(3) = pos2(3) - p%MSL2SWL
         r2      = mem%RMGB(N+1)
         if (mem%i_floor == 0) then  ! both ends above or at seabed
            ! Compute loads on the end plate of node 1
            IF (p%WaveStMod > 0) THEN
               CALL GetTotalWaveElev( Time, pos1, Zeta1, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               CALL GetFreeSurfaceNormal( Time, pos1, r1, n_hat, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FSPt = (/pos1(1),pos1(2),Zeta1/) ! Reference point on the free surface
            ELSE
               FSPt = (/pos1(1),pos1(2),0.0/)
               n_hat = (/0.0,0.0,1.0/)
            END IF
            CALL GetSectionUnitVectors( k_hat1, y_hat, z_hat )
            CALL GetSectionFreeSurfaceIntersects( REAL(pos1,DbKi), REAL(FSPt,DbKi), k_hat1, y_hat, z_hat, n_hat, REAL(r1,DbKi), theta1, theta2, secStat)
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL GetEndPlateHstLds(pos1, k_hat1, y_hat, z_hat, r1, theta1, theta2, F_B_End)
            m%F_B_End(:, mem%NodeIndx(  1)) = m%F_B_End(:, mem%NodeIndx(  1)) + F_B_End
            IF (mem%MHstLMod == 1) THEN ! Check for partially wetted end plates
               IF ( .NOT.( EqualRealNos((theta2-theta1),0.0_DbKi) .OR. EqualRealNos((theta2-theta1),2.0_DbKi*PI_D) ) ) THEN
                   CALL SetErrStat(ErrID_Warn, 'End plate is partially wetted with MHstLMod = 1. The buoyancy load and distribution potentially have large error. This has happened to the first node of Member ID ' //trim(num2lstr(mem%MemberID)), errStat, errMsg, RoutineName )
               END IF
            END IF
            ! Compute loads on the end plate of node N+1
            IF (p%WaveStMod > 0) THEN
               CALL GetTotalWaveElev( Time, pos2, Zeta2, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               CALL GetFreeSurfaceNormal( Time, pos2, r2, n_hat, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FSPt = (/pos2(1),pos2(2),Zeta2/) ! Reference point on the free surface
            ELSE
               FSPt = (/pos2(1),pos2(2),0.0/)
               n_hat = (/0.0,0.0,1.0/)
            END IF
            CALL GetSectionUnitVectors( k_hat2, y_hat, z_hat )
            CALL GetSectionFreeSurfaceIntersects( REAL(pos2,DbKi), REAL(FSPt,DbKi), k_hat2, y_hat, z_hat, n_hat, REAL(r2,DbKi), theta1, theta2, secStat)
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL GetEndPlateHstLds(pos2, k_hat2, y_hat, z_hat, r2, theta1, theta2, F_B_End)
            m%F_B_End(:, mem%NodeIndx(N+1)) = m%F_B_End(:, mem%NodeIndx(N+1)) - F_B_End
            IF (mem%MHstLMod == 1) THEN ! Check for partially wetted end plates
               IF ( .NOT.( EqualRealNos((theta2-theta1),0.0_DbKi) .OR. EqualRealNos((theta2-theta1),2.0_DbKi*PI_D) ) ) THEN
                   CALL SetErrStat(ErrID_Warn, 'End plate is partially wetted with MHstLMod = 1. The buoyancy load and distribution potentially have large error. This has happened to the last node of Member ID ' //trim(num2lstr(mem%MemberID)), errStat, errMsg, RoutineName )
               END IF
            END IF
         elseif ( mem%doEndBuoyancy ) then ! The member crosses the seabed line so only the upper end potentially have hydrostatic load
            ! Only compute the loads on the end plate of node N+1
            IF (p%WaveStMod > 0) THEN
               CALL GetTotalWaveElev( Time, pos2, Zeta2, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               CALL GetFreeSurfaceNormal( Time, pos2, r2, n_hat, ErrStat2, ErrMsg2 )
                 CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               FSPt = (/pos2(1),pos2(2),Zeta2/) ! Reference point on the free surface
            ELSE
               FSPt = (/pos2(1),pos2(2),0.0/)
               n_hat = (/0.0,0.0,1.0/)
            END IF
            CALL GetSectionUnitVectors( k_hat2, y_hat, z_hat )
            CALL GetSectionFreeSurfaceIntersects( REAL(pos2,DbKi), REAL(FSPt,DbKi), k_hat2, y_hat, z_hat, n_hat, REAL(r2,DbKi), theta1, theta2, secStat)
              CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            CALL GetEndPlateHstLds(pos2, k_hat2, y_hat, z_hat, r2, theta1, theta2, F_B_End)
            m%F_B_End(:, mem%NodeIndx(N+1)) = m%F_B_End(:, mem%NodeIndx(N+1)) - F_B_End
            IF (mem%MHstLMod == 1) THEN ! Check for partially wetted end plates
               IF ( .NOT.( EqualRealNos((theta2-theta1),0.0_DbKi) .OR. EqualRealNos((theta2-theta1),2.0_DbKi*PI_D) ) ) THEN
                   CALL SetErrStat(ErrID_Warn, 'End plate is partially wetted with MHstLMod = 1. The buoyancy load and distribution potentially have large error. This has happened to the last node of Member ID ' //trim(num2lstr(mem%MemberID)), errStat, errMsg, RoutineName )
               END IF
            END IF
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
   !TODO: Where's F_WMF_End computed?
      
   DO J = 1, p%NJoints
     
      ! Obtain the node index because WaveVel, WaveAcc, and WaveDynP are defined in the node indexing scheme, not the markers (No longer relevant?)
      ! The first NJoints nodes are all the joints with the rest being the internal nodes. See Morison_GenerateSimulationNodes.
            
      ! NOTE: 
      ! The PropPot values are only for members, and when the p%AM_End, p%DP_Const_End, p%Mass_MG_End, and p%I_MG_End are computed at init,
      ! contributions to these values are added only if the member connecting to the joint is NOT modeled with potential flow theory
      ! However, the p%An_End term used data from ALL members attached to a node, regardless of the PropPot setting, because the drag force is alway on.
      ! Therefore, no need to check PropPot here.
      
      ! Effect of wave stretching already baked into m%FDynP, m%FA, and m%vrel. No additional modification needed.
         
      ! Lumped added mass loads
      qdotdot                 = reshape((/u%Mesh%TranslationAcc(:,J),u%Mesh%RotationAcc(:,J)/),(/6/)) 
      m%F_A_End(:,J)          = m%nodeInWater(j) * matmul( p%AM_End(:,:,J) , ( - qdotdot(1:3)) )
         
      ! TODO: The original code did not multiply by nodeInWater, but should we? GJH
      ! Should be ok because m%FDynP and m%FA are both zeroed above the SWL (when WaveStMod=0) or the instantaneous free surface (when WaveStMod>0)
      m%F_I_End(:,J) =   (p%DP_Const_End(:,j) * m%FDynP(j) + matmul(p%AM_End(:,:,j),m%FA(:,j)))
         
      ! Marine growth inertia: ends: Section 4.2.2
      ! With wave stretching, m%nodeInWater is based on the instantaneous free surface and the current body position if (WaveDisp/=0).
      ! This should still be ok because with wave stretching, we do not allow joints to come out of water if initially submerged or 
      ! enter water if initially out of water. This is enforced when computing the side loads above.
      m%F_IMG_End(1:3,j) = -m%nodeInWater(j) * p%Mass_MG_End(j)*qdotdot(1:3)
      m%F_IMG_End(4:6,j) = -m%nodeInWater(j) * (matmul(p%I_MG_End(:,:,j),qdotdot(4:6)) - cross_product(u%Mesh%RotationVel(:,J),matmul(p%I_MG_End(:,:,j),u%Mesh%RotationVel(:,J))))

      ! Compute the dot product of the relative velocity vector with the directional Area of the Joint
      ! m%nodeInWater(j) is probably not necessary because m%vrel is zeroed when the node is out of water
      vmag  = m%nodeInWater(j) * ( m%vrel(1,j)*p%An_End(1,J) + m%vrel(2,j)*p%An_End(2,J) + m%vrel(3,j)*p%An_End(3,J) )
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
               m%F_D_End(i,j) = (1.0_ReKi - p%DragLoFSc_End(j)) * p%An_End(i,j) * p%DragConst_End(j) * abs(vmagf)*vmagf &   
                                          + p%DragLoFSc_End(j)  * p%An_End(i,j) * p%DragConst_End(j) * abs(vmag )*vmag  
            ELSE IF (p%DragMod_End(J) .EQ. 1_IntKi) THEN
               ! Note: vmag is zero if node is not in the water
               m%F_D_End(i,j) = (1.0_ReKi - p%DragLoFSc_End(j)) * p%An_End(i,j) * p%DragConst_End(j) * abs(vmagf)*max(vmagf,0.0_ReKi) &
                                          + p%DragLoFSc_End(j)  * p%An_End(i,j) * p%DragConst_End(j) * abs(vmag) *max(vmag, 0.0_ReKi)
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
        


   CONTAINS
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
      Zeta = SeaSt_Interp_3D( Time, pos(1:2), p%WaveElev1, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (associated(p%WaveElev2)) THEN
         Zeta = Zeta + SeaSt_Interp_3D( Time, pos(1:2), p%WaveElev2, p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
           CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF
      
   END SUBROUTINE GetTotalWaveElev

   SUBROUTINE GetFreeSurfaceNormal( Time, pos, r, n, ErrStat, ErrMsg)
      REAL(DbKi),      INTENT( In    ) :: Time
      REAL(ReKi),      INTENT( In    ) :: pos(*)  ! Position at which free-surface normal is to be calculated. Third entry ignored if present.
      REAL(ReKi),      INTENT( In    ) :: r       ! Distance for central differencing
      REAL(ReKi),      INTENT(   OUT ) :: n(3)    ! Free-surface normal vector
      INTEGER(IntKi),  INTENT(   OUT ) :: ErrStat ! Error status of the operation
      CHARACTER(*),    INTENT(   OUT ) :: ErrMsg  ! Error message if errStat /= ErrID_None
      REAL(ReKi)                       :: r1,ZetaP,ZetaM,dZetadx,dZetady
      CHARACTER(*),    PARAMETER       :: RoutineName = 'GetFreeSurfaceNormal'
      INTEGER(IntKi)                   :: errStat2
      CHARACTER(ErrMsgLen)             :: errMsg2
      ErrStat   = ErrID_None
      ErrMsg    = ""

      r1 = MAX(r,1.0e-6) ! In case r is zero

      CALL GetTotalWaveElev( Time, (/pos(1)+r1,pos(2)/), ZetaP, ErrStat2, ErrMsg2 )
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL GetTotalWaveElev( Time, (/pos(1)-r1,pos(2)/), ZetaM, ErrStat2, ErrMsg2 )
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      dZetadx = (ZetaP-ZetaM)/(2.0_ReKi*r1)
      
      CALL GetTotalWaveElev( Time, (/pos(1),pos(2)+r1/), ZetaP, ErrStat2, ErrMsg2 )
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      CALL GetTotalWaveElev( Time, (/pos(1),pos(2)-r1/), ZetaM, ErrStat2, ErrMsg2 )
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      dZetady = (ZetaP-ZetaM)/(2.0_ReKi*r1)
      
      n = (/-dZetadx,-dZetady,1.0_ReKi/)
      n = n / SQRT(Dot_Product(n,n))

   END SUBROUTINE GetFreeSurfaceNormal

   SUBROUTINE GetSectionUnitVectors( k, y, z )
      REAL(ReKi),      INTENT( In    ) :: k(3) ! Member axial unit vector
      REAL(ReKi),      INTENT(   OUT ) :: y(3) ! Horizontal unit vector perpendicular to k
      REAL(ReKi),      INTENT(   OUT ) :: z(3) ! Unit vector perpendicular to k and y with positive vertical component
      IF ( ABS(k(3)) > 0.999999_ReKi ) THEN ! k is effectively vertical
         y = (/0.0,1.0,0.0/)
      ELSE
         y = (/-k(2),k(1),0.0/)
         y = y / SQRT(Dot_Product(y,y))      
      ENDIF
      z = cross_product(k,y)
      IF ( z(3) < 0.0 ) THEN ! Flip y and z so z points upward
         y = -y;
         z = -z;
      END IF
   END SUBROUTINE GetSectionUnitVectors

   SUBROUTINE GetSectionFreeSurfaceIntersects( pos0, FSPt, k_hat, y_hat, z_hat, n_hat, R, theta1, theta2, secStat)
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
      CHARACTER(*),    PARAMETER       :: RoutineName = 'GetSectionFreeSurfaceIntersects'

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

   END SUBROUTINE GetSectionFreeSurfaceIntersects

   SUBROUTINE GetSectionHstLds( origin, pos0, k_hat, y_hat, z_hat, R, dRdl, theta1, theta2, dFdl)

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
      dFdl = dFdl * p%WtrDens * g

   END SUBROUTINE GetSectionHstLds

   SUBROUTINE getElementHstLds_Mod2( pos1In, pos2In, FSPtIn, k_hatIn, y_hatIn, z_hatIn, n_hatIn, r1In, r2In, dlIn, F_B1, F_B2, ErrStat, ErrMsg )
      
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
      CHARACTER(*),    PARAMETER       :: routineName = "getElementHstLds_Mod2"
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
      CALL GetSectionFreeSurfaceIntersects( pos1,   FSPt, REAL(k_hat,ReKi), REAL(y_hat,ReKi), REAL(z_hat,ReKi), REAL(n_hat,ReKi), r1, theta1, theta2, secStat1)
      CALL GetSectionHstLds( pos1, pos1,   k_hat, y_hat, z_hat, r1,   dRdl, theta1, theta2, dFdl1)

      ! Section load at midpoint
      CALL GetSectionFreeSurfaceIntersects( posMid, FSPt, REAL(k_hat,ReKi), REAL(y_hat,ReKi), REAL(z_hat,ReKi), REAL(n_hat,ReKi), rMid, theta1, theta2, secStatMid)
      CALL GetSectionHstLds( pos1, posMid, k_hat, y_hat, z_hat, rMid, dRdl, theta1, theta2, dFdlMid)

      ! Section load at node 2
      CALL GetSectionFreeSurfaceIntersects( pos2,   FSPt, REAL(k_hat,ReKi), REAL(y_hat,ReKi), REAL(z_hat,ReKi), REAL(n_hat,ReKi), r2, theta1, theta2, secStat2)
      CALL GetSectionHstLds( pos1, pos2,   k_hat, y_hat, z_hat, r2,   dRdl, theta1, theta2, dFdl2)

      ! Adaptively refine the load integration over the element
      CALL RefineElementHstLds(pos1,pos1,posMid,pos2,FSPt,r1,rMid,r2,dl,dRdl,secStat1,secStatMid,secStat2,k_hat,y_hat,z_hat,n_hat,dFdl1,dFdlMid,dFdl2,1,F_B,ErrStat2,ErrMsg2)
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! Distribute the hydrostatic load to the two end nodes
      F_B1(1:3) = 0.5 * F_B(1:3)
      F_B2(1:3) = 0.5 * F_B(1:3)
      F_B(4:6)  = F_B(4:6) - CROSS_PRODUCT(k_hat*dl,F_B2(1:3))
      F_B1(4:6) = 0.5 * F_B(4:6)
      F_B2(4:6) = 0.5 * F_B(4:6)

   END SUBROUTINE getElementHstLds_Mod2

   RECURSIVE SUBROUTINE RefineElementHstLds( origin, pos1, posMid, pos2, FSPt, r1, rMid, r2, dl, dRdl,secStat1,secStatMid,secStat2, k_hat, y_hat, z_hat, n_hat, dFdl1, dFdlMid, dFdl2, recurLvl, F_B_5pt, ErrStat, ErrMsg)

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
      CHARACTER(*),    PARAMETER       :: RoutineName = "RefineElementHstLds"
      
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
      CALL GetSectionFreeSurfaceIntersects( posMidL, FSPt, REAL(k_hat,ReKi), REAL(y_hat,ReKi), REAL(z_hat,ReKi), REAL(n_hat,ReKi), rMidL, theta1, theta2, secStatMidL)
      CALL GetSectionHstLds( origin, posMidL, k_hat, y_hat, z_hat, rMidL, dRdl, theta1, theta2, dFdlMidL)

      ! Mid point of right section
      CALL GetSectionFreeSurfaceIntersects( posMidR, FSPt, REAL(k_hat,ReKi), REAL(y_hat,ReKi), REAL(z_hat,ReKi), REAL(n_hat,ReKi), rMidR, theta1, theta2, secStatMidR)
      CALL GetSectionHstLds( origin, posMidR, k_hat, y_hat, z_hat, rMidR, dRdl, theta1, theta2, dFdlMidR)
      
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
         CALL RefineElementHstLds(origin,pos1,posMidL,posMid,FSPt,r1,rMidL,rMid,0.5*dl,dRdl,secStat1,secStatMidL,secStatMid,k_hat,y_hat,z_hat,n_hat,dFdl1,dFdlMidL,dFdlMid, recurLvl+1, tmp, ErrStat, ErrMsg)
         CALL RefineElementHstLds(origin,posMid,posMidR,pos2,FSPt,rMid,rMidR,r2,0.5*dl,dRdl,secStatMid,secStatMidR,secStat2,k_hat,y_hat,z_hat,n_hat,dFdlMid,dFdlMidR,dFdl2, recurLvl+1, F_B_5pt, ErrStat, ErrMsg)
         F_B_5pt = F_B_5pt + tmp
      END IF

   END SUBROUTINE RefineElementHstLds

   SUBROUTINE GetEndPlateHstLds(pos0, k_hat, y_hat, z_hat, R, theta1, theta2, F)

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
      F(1:3) = p%WtrDens * g * Fk * k_hat

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
      F(4:6) = p%WtrDens * g * (My*y_hat + Mz*z_hat)

   END SUBROUTINE GetEndPlateHstLds

   SUBROUTINE getElementHstLds_Mod1( Time, pos1, pos2, Zeta1, Zeta2, k_hat, r1, r2, dl, alphaIn, Is1stElement, F_B0, F_B1, F_B2, ErrStat, ErrMsg )
      
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
         FbVec = p%WtrDens * g * FbVec

         ! Hydrostatic moment on element about the lower node
         MbVec = (Vhc+0.25*Pi*(r2**4-r1**4)) * Cross_Product(k_hat,(/0.0_ReKi,0.0_ReKi,1.0_ReKi/))
         MbVec = p%WtrDens * g * MbVec

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
         IF ( p%WaveStMod > 0_IntKi ) THEN ! If wave stretching is enabled, compute free surface normal
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
         FbVec = p%WtrDens * g * FbVec

         ! Hydrostatic moment on element about the lower node
         MbVec = Cross_Product( Vrc*r_hat+Vhc*k_hat, (/0.0_ReKi,0.0_ReKi,1.0_ReKi/) ) &
                 + 0.25*Pi*a0b0* ( ( s_hat(3)*a0*a0 + 4.0*(s0-h0*sinGamma)*Z0 )*t_hat - t_hat(3)*b0*b0*s_hat ) &
                 - 0.25*Pi*r1**4*(   r_hat(3)                                  *t_hat - t_hat(3)   *   r_hat )
         MbVec = p%WtrDens * g * MbVec

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
SUBROUTINE DistributeElementLoads(Fl, Fr, M, sinPhi, cosPhi, SinBeta, cosBeta, alpha, F1, F2)
   
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
   
   F1(1) =  cosBeta*(Fl*sinPhi + Fr*cosPhi)*alpha
   F1(2) =  sinBeta*(Fl*sinPhi + Fr*cosPhi)*alpha
   F1(3) =          (Fl*cosPhi - Fr*sinPhi)*alpha
   F1(4) = -sinBeta * M                    *alpha
   F1(5) =  cosBeta * M                    *alpha
   F1(6) =  0.0
      
   F2(1) =  cosBeta*(Fl*sinPhi + Fr*cosPhi)*(1-alpha)
   F2(2) =  sinBeta*(Fl*sinPhi + Fr*cosPhi)*(1-alpha)
   F2(3) =          (Fl*cosPhi - Fr*sinPhi)*(1-alpha)
   F2(4) = -sinBeta * M                    *(1-alpha)
   F2(5) =  cosBeta * M                    *(1-alpha)
   F2(6) =  0.0
   
END SUBROUTINE DistributeElementLoads
!----------------------------------------------------------------------------------------------------------------------------------
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
   REAL(ReKi)                                        :: WtrDpth
   REAL(ReKi)                                        :: pos(3), posPrime(3), positionXY(2)
   REAL(SiKi)                                        :: WaveElev, WaveElev1, WaveElev2
   REAL(ReKi)                                        :: vrel(3), FV(3), vmag, vmagf
   INTEGER(IntKi)                                    :: errStat2
   CHARACTER(ErrMsgLen)                              :: errMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'Morison_UpdateDiscState'
   
   ! Initialize errStat  
   errStat = ErrID_None         
   errMsg  = ""               
   
   ! Water depth measured from the free surface
   WtrDpth = p%WtrDpth + p%MSL2SWL

   ! Update state of the relative normal velocity high-pass filter at each joint
   DO j = 1, p%NJoints 
      ! Get joint position
      IF (p%WaveDisp == 0 ) THEN
         ! use the initial X,Y location
         pos(1) = u%Mesh%Position(1,j)
         pos(2) = u%Mesh%Position(2,j)
      ELSE
         ! Use current X,Y location
         pos(1) = u%Mesh%TranslationDisp(1,j) + u%Mesh%Position(1,j)
         pos(2) = u%Mesh%TranslationDisp(2,j) + u%Mesh%Position(2,j)
      END IF
      IF (p%WaveStMod > 0 .AND. p%WaveDisp /= 0) THEN ! Wave stretching enabled
        pos(3) = u%Mesh%Position(3,j) + u%Mesh%TranslationDisp(3,j) - p%MSL2SWL  ! Use the current Z location.
      ELSE ! Wave stretching disabled
        pos(3) = u%Mesh%Position(3,j) - p%MSL2SWL  ! We are intentionally using the undisplaced Z position of the node.
      END IF      
      ! Compute the free surface elevation at the x/y position of the joint
      positionXY = (/pos(1),pos(2)/)      
      WaveElev1 = SeaSt_Interp_3D( Time, positionXY, p%WaveElev1, p%seast_interp_p, m%SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
        CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (associated(p%WaveElev2)) THEN
        WaveElev2 = SeaSt_Interp_3D( Time, positionXY, p%WaveElev2, p%seast_interp_p, m%SeaSt_Interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        WaveElev  = WaveElev1 + WaveElev2
      ELSE
        WaveElev  = WaveElev1 
      END IF
      ! Compute fluid and relative velocity at the joint
      IF (p%WaveStMod == 0) THEN ! No wave stretching
          IF ( pos(3) <= 0.0_ReKi) THEN ! Node is at or below the SWL
              ! Use location to obtain interpolated values of kinematics         
              call SeaSt_Interp_Setup( Time, pos, p%seast_interp_p, m%seast_interp_m, ErrStat2, ErrMsg2 ) 
                call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              FV   = SeaSt_Interp_4D_Vec( p%WaveVel, m%seast_interp_m, ErrStat2, ErrMsg2 )
                call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              vrel = FV - u%Mesh%TranslationVel(:,j)
          ELSE ! Node is above the SWL
              vrel = 0.0_ReKi
          END IF
      ELSE ! Wave stretching enabled
          IF ( pos(3) <= WaveElev ) THEN ! Node is submerged
              IF ( p%WaveStMod < 3 ) THEN ! Vertical or extrapolated wave stretching
                  IF ( pos(3) <= 0.0_ReKi) THEN ! Node is below the SWL - evaluate wave dynamics as usual
                      ! Use location to obtain interpolated values of kinematics         
                      call SeaSt_Interp_Setup( Time, pos, p%seast_interp_p, m%seast_interp_m, ErrStat2, ErrMsg2 ) 
                        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                      FV = SeaSt_Interp_4D_Vec( p%WaveVel, m%seast_interp_m, ErrStat2, ErrMsg2 )
                        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  ELSE ! Node is above SWL - need wave stretching
                      ! Vertical wave stretching
                      FV = SeaSt_Interp_3D_vec( Time, positionXY, p%WaveVel0, p%seast_interp_p,  m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 )
                        call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                      ! Extrapoled wave stretching
                      IF (p%WaveStMod == 2) THEN 
                        FV = FV + SeaSt_Interp_3D_vec( Time, positionXY, p%PWaveVel0,  p%seast_interp_p, m%seast_interp_m%FirstWarn_Clamp, ErrStat2, ErrMsg2 ) * pos(3)
                          call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                      END IF
                  END IF ! Node is submerged
              ELSE ! Wheeler stretching - no need to check whether the node is above or below SWL
                  ! Map the node z-position linearly from [-WtrDpth,WaveElev] to [-WtrDpth,0] 
                  posPrime = pos
                  posPrime(3) = WtrDpth*(WtrDpth+pos(3))/(WtrDpth+WaveElev)-WtrDpth
                  ! Obtain the wave-field variables by interpolation with the mapped position.
                  call SeaSt_Interp_Setup( Time, posPrime, p%seast_interp_p, m%seast_interp_m, ErrStat2, ErrMsg2 ) 
                    call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                  FV = SeaSt_Interp_4D_Vec( p%WaveVel, m%seast_interp_m, ErrStat2, ErrMsg2 )
                    call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
              END IF
              vrel = FV - u%Mesh%TranslationVel(:,j)
          ELSE ! Node is out of water - zero-out all wave dynamics
              vrel = 0.0_ReKi
          END IF ! If node is in or out of water
      END IF ! If wave stretching is on or off
      ! Compute the dot product of the relative velocity vector with the directional Area of the Joint
      vmag  = vrel(1)*p%An_End(1,J) + vrel(2)*p%An_End(2,J) + vrel(3)*p%An_End(3,J)
      ! High-pass filtering
      vmagf = p%VRelNFiltConst(J) * (vmag + xd%V_rel_n_FiltStat(J))
      ! Update relative normal velocity filter state for joint J 
      xd%V_rel_n_FiltStat(J) = vmagf-vmag
   END DO ! j = 1, p%NJoints
END SUBROUTINE Morison_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE Morison
!**********************************************************************************************************************************
