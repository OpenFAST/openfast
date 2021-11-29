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
   PUBLIC:: Morison_GenerateSimulationNodes
   
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

subroutine GetOrientationAngles(p1, p2, phi, sinPhi, cosPhi, tanPhi, sinBeta, cosBeta, k_hat, errStat, errMsg)
   real(ReKi),   intent(in   ) :: p1(3),p2(3)
   real(ReKi),   intent(  out) :: phi, sinPhi, cosPhi, tanPhi, sinBeta, cosBeta, k_hat(3)
   integer,      intent(  out) :: errStat              ! returns a non-zero value when an error occurs  
   character(*), intent(  out) :: errMsg               ! Error message if errStat /= ErrID_None
   
   
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
            call SeterrStat(ErrID_Fatal, 'An element of the Morison structure has co-located endpoints!  This should never occur.  Please review your model.', errStat, errMsg, 'Morison_CalcOutput' )
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


SUBROUTINE WriteSummaryFile( UnSum, g, MSL2SWL, WtrDpth, numJoints, numNodes, nodes, numMembers, members, &
                             NOutputs, OutParam, NMOutputs, MOutLst,  NJOutputs, JOutLst, uMesh, yMesh, &
                             p, m, errStat, errMsg ) 
                             
   INTEGER,                               INTENT ( IN    )  :: UnSum
   REAL(ReKi),                            INTENT ( IN    )  :: g                    ! gravity
   REAL(ReKi),                            INTENT ( IN    )  :: MSL2SWL
   REAL(ReKi),                            INTENT ( IN    )  :: WtrDpth
   INTEGER,                               INTENT ( IN    )  :: numJoints
   INTEGER,                               INTENT ( IN    )  :: numNodes
   TYPE(Morison_NodeType),   ALLOCATABLE, INTENT ( IN    )  :: nodes(:)  
   INTEGER,                               INTENT ( IN    )  :: numMembers
   TYPE(Morison_MemberType), ALLOCATABLE, INTENT ( IN    )  :: members(:)
   INTEGER,                               INTENT ( IN    )  :: NOutputs
   TYPE(OutParmType),        ALLOCATABLE, INTENT ( IN    )  :: OutParam(:)
   INTEGER,                               INTENT ( IN    )  :: NMOutputs
   TYPE(Morison_MOutput),    ALLOCATABLE, INTENT ( IN    )  :: MOutLst(:)
   INTEGER,                               INTENT ( IN    )  :: NJOutputs
   TYPE(Morison_JOutput),    ALLOCATABLE, INTENT ( IN    )  :: JOutLst(:)
   TYPE(MeshType),                        INTENT ( INOUT )  :: uMesh
   TYPE(MeshType),                        INTENT ( INOUT )  :: yMesh
   TYPE(Morison_ParameterType),           INTENT ( IN    ) :: p
   TYPE(Morison_MiscVarType),             INTENT ( IN    ) :: m 
   INTEGER,                               INTENT (   OUT )  :: errStat             ! returns a non-zero value when an error occurs  
   CHARACTER(*),                          INTENT (   OUT )  :: errMsg              ! Error message if errStat /= ErrID_None

   INTEGER                                     :: I, J
   REAL(ReKi)                                  :: l                   ! length of an element
   LOGICAL                                     :: filledFlag          ! flag indicating if element is filled/flooded
   CHARACTER(2)                                :: strFmt
   CHARACTER(ChanLen)                          :: strNodeType         ! string indicating type of node: End, Interior, Super
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
   real(ReKi)                                  :: pos(3), pos2(3)     ! Position of a node or joint in the MSL inertial system
   INTEGER                                     :: mbrIndx, nodeIndx, c, N
   CHARACTER(ChanLen)                          :: tmpName
   REAL(ReKi)                                  :: totalFillMass, mass_fill, fillVol, memberVol
   REAL(ReKi)                                  :: totalMGMass, mass_MG
   TYPE(Morison_NodeType)                      ::  node1, node2
   real(ReKi)                                  :: ptLoad(6)
   logical                                     :: fillFlag
   type(Morison_MemberType)                    :: mem
   REAL(ReKi)                                  :: Cd1, Cd2, Ca1, Ca2, Cp1, Cp2, AxCd1, AxCd2, AxCa1, AxCa2, AxCp1, AxCp2, JAxCd1, JAxCd2, JAxCa1, JAxCa2, JAxCp1, JAxCp2 ! tmp coefs
   real(ReKi)                                  :: F_B(6, numNodes), F_BF(6, numNodes), F_WMG(6, numNodes)
      ! Initialize data
   errStat       = ErrID_None
   errMsg        = ""
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
   
   IF ( UnSum > 0 ) THEN
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
                     ,errStat           = errStat           &
                     ,ErrMess           = errMsg            &
                     ,Force             = .TRUE.            &
                     ,Moment            = .TRUE.            &
                     )
         ! Create the node on the mesh
 
      CALL MeshPositionNode (WRP_Mesh                              &
                              , 1                                  &
                              , (/0.0_ReKi, 0.0_ReKi, 0.0_ReKi/)   &  
                              , errStat                            &
                              , errMsg                             &
                              )
      
      IF ( errStat /= 0 ) RETURN
       
      
         ! Create the mesh element
      CALL MeshConstructElement (  WRP_Mesh            &
                                  , ELEMENT_POINT      &                         
                                  , errStat            &
                                  , errMsg             &
                                  , 1                  &
                                )
      CALL MeshCommit ( WRP_Mesh           &
                      , errStat            &
                      , errMsg             )
   
      IF ( errStat /= ErrID_None ) RETURN
            
         ! we need the translation displacement mesh for loads transfer:
      CALL MeshCopy ( SrcMesh  = WRP_Mesh            &
                    , DestMesh = WRP_Mesh_position   &
                    , CtrlCode = MESH_SIBLING        &
                    , IOS      = COMPONENT_INPUT     &
                    , TranslationDisp = .TRUE.       &
                    , errStat  = errStat             &
                    , ErrMess  = errMsg              )  ! automatically sets    DestMesh%RemapFlag = .TRUE.
                    
      IF ( errStat /= ErrID_None ) RETURN
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
         
      CALL MeshMapCreate           ( yMesh, WRP_Mesh, M_P_2_P, errStat, errMsg                )
        !CALL CheckError( errStat, 'Message from MeshMapCreate HD_M_L_2_ED_P: '//NewLine//errMsg )
      CALL Transfer_Point_to_Point( yMesh, WRP_Mesh, M_P_2_P, errStat, errMsg, uMesh, WRP_Mesh_position )
      
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
      CALL Transfer_Point_to_Point( yMesh, WRP_Mesh, M_P_2_P, errStat, errMsg, uMesh, WRP_Mesh_position )
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
      CALL Transfer_Point_to_Point( yMesh, WRP_Mesh, M_P_2_P, errStat, errMsg, uMesh, WRP_Mesh_position )
      MG_Wt(1:3) = WRP_Mesh%Force(:,1)
      MG_Wt(4:6) = WRP_Mesh%Moment(:,1)
      !
       CALL MeshMapDestroy( M_P_2_P, errStat, errMsg ); IF ( errStat /= ErrID_None ) CALL WrScr(TRIM(errMsg))
     
           
      WRITE( UnSum,  '(//)' ) 
      WRITE( UnSum, '(A36)' ) 'Weight loads about ( 0.0, 0.0, 0.0 )'
      WRITE( UnSum, '(A36)' ) '------------------------------------'
      WRITE( UnSum, '(18x,6(2X,A20))' ) '  MGFxi  ', '  MGFyi  ', '  MGFzi  ', '  MGMxi  ', '  MGMyi  ', '  MGMzi  '
      WRITE( UnSum, '(18x,6(2X,A20))' ) '   (N)   ', '   (N)   ', '   (N)   ', '  (N-m)  ', '  (N-m)  ', '  (N-m)  '
      WRITE( UnSum, '(A18,6(2X,ES20.6))') 'Marine Growth:   ', MG_Wt(1), MG_Wt(2), MG_Wt(3), MG_Wt(4), MG_Wt(5), MG_Wt(6)

      
      CALL MeshDestroy(WRP_Mesh, errStat, errMsg ); IF ( errStat /= ErrID_None ) CALL WrScr(TRIM(errMsg))
      CALL MeshDestroy(WRP_Mesh_position, errStat, errMsg ); IF ( errStat /= ErrID_None ) CALL WrScr(TRIM(errMsg))
      !
      !   ! Write the header for this section
      WRITE( UnSum,  '(//)' ) 
      WRITE( UnSum,  '(A14,I4,A44)' ) 'Nodes (first [',numJoints,'] are joints, remainder are internal nodes)'
      WRITE( UnSum,  '(/)' ) 
      WRITE( UnSum, '(1X,A5,20(2X,A10))' ) '  i  ', '  MbrIndx ', '   Nxi    ', '   Nyi    ', '   Nzi    ', '     R    ', '    t     ', '   tMG    ', '  MGDens  ', ' PropPot  ', 'FilledFlag', 'FilledMass', '    Cd    ', '    Ca    ', '    Cp    ', '   AxCd   ',  '   AxCa   ', '   AxCp   ', '   JAxCd  ', '   JAxCa  ', '   JAxCp  '
      WRITE( UnSum, '(1X,A5,20(2X,A10))' ) ' (-) ', '    (-)   ', '   (m)    ', '   (m)    ', '   (m)    ', '    (m)   ', '   (m)    ', '   (m)    ', ' (kg/m^3) ', '   (-)    ', '   (-)    ', '  (kg)    ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ',  '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   ', '    (-)   '
      
         ! Write the node data
      do I = 1,numJoints   
         ! need to add MSL2SWL offset from this because the Positons are relative to SWL, but we should report them relative to MSL here
         pos = nodes(i)%Position
         pos(3) = pos(3) + MSL2SWL
         write( UnSum, '(1X,I5,(2X,A10),3(2X,F10.4),2(2X,A10),2(2X,ES10.3),9(2X,A10),3(2X,ES10.3))' ) i,'    -     ', pos, '    -     ',  '    -     ',  nodes(i)%tMG,  nodes(i)%MGdensity,  '    -     ',  '    -     ',  '    -     ', '    -     ',  '    -     ',  '    -     ',  '    -     ',  '    -     ',  '    -     ',  nodes(i)%JAxCd,  nodes(i)%JAxCa, nodes(i)%JAxCp
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
            write( UnSum, '(1X,I5,(2X,I10),3(2X,F10.4),4(2X,ES10.3),2(6X,L6),7(2X,ES10.3),3(7x,A5))' ) c, members(j)%MemberID, pos, members(j)%R(i),  members(j)%R(i)-members(j)%Rin(i),  members(j)%tMG(i),  members(j)%MGdensity(i),  members(j)%PropPot,  fillFlag,  members(j)%m_fb_u(i)+members(j)%m_fb_l(i),  members(j)%Cd(i),  members(j)%Ca(i),  members(j)%Cp(i),  members(j)%AxCd(i),  members(j)%AxCa(i),  members(j)%AxCp(i), '  -  ',  '  -  ',  '  -  '
         end do
      end do
      
      
      write( UnSum,  '(//)' ) 
      write( UnSum,  '(A8)' ) 'Members'
      write( UnSum,  '(/)' ) 
      write( UnSum, '(1X,A8,2X,A6,2X,A6,31(2X,A12))' ) 'MemberID', 'joint1','joint2','  Length  ', '   NElem    ', '   Volume   ', '  MGVolume  ', '      R1    ', '     t1     ', '      R2    ', '     t2     ', ' PropPot  ', 'FilledFlag', 'FillDensity', '  FillFSLoc ', '  FillMass  ', '     Cd1    ', '    Ca1   ', '     Cp1    ', '    AxCd1   ', '    AxCa1   ', '    AxCp1   ', '   JAxCd1   ', '   JAxCa1   ', '  JAxCp1    ', '     Cd2    ', '     Ca2    ', '     Cp2    ', '    AxCd2   ', '    AxCa2   ', '    AxCp2   ', '   JAxCd2   ', '   JAxCa2   ', '   JAxCp2   '
      write( UnSum, '(1X,A8,2X,A6,2X,A6,31(2X,A12))' ) '  (-)   ', ' (-)  ',' (-)  ','   (m)    ', '    (-)     ', '   (m^3)    ', '   (m^3)    ', '      (m)   ', '     (m)    ', '      (m)   ', '     (m)    ', '   (-)    ', '   (-)    ', ' (kg/m^3)  ', '     (-)    ', '    (kg)    ', '     (-)    ', '    (-)   ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '    (-)     ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    ', '     (-)    '
      
      
      
      
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
     
         JAxCd1 = nodes(members(i)%NodeIndx(1  ))%JAxCd
         JAxCd2 = nodes(members(i)%NodeIndx(1+N))%JAxCd
         JAxCa1 = nodes(members(i)%NodeIndx(1  ))%JAxCa
         JAxCa2 = nodes(members(i)%NodeIndx(1+N))%JAxCa
         JAxCp1 = nodes(members(i)%NodeIndx(1  ))%JAxCp
         JAxCp2 = nodes(members(i)%NodeIndx(1+N))%JAxCp
       
         
         write( UnSum, '(1X,I8,2X,I6,2X,I6,2X,ES12.5,2X,I12, 6(2X,ES12.5),2(2X,L12),21(2X,ES12.5))' )  members(i)%MemberID, &
                       members(i)%NodeIndx(1), members(i)%NodeIndx(N+1), members(i)%RefLength, N, &
                       memberVol, MGvolume, members(i)%Rmg(1), members(i)%Rmg(1)-members(i)%Rin(1), &
                       members(i)%Rmg(N+1), members(i)%Rmg(N+1)-members(i)%Rin(N+1),  &
                       members(i)%PropPot, filledFlag, members(i)%FillDens, members(i)%FillFSLoc, &
                       mass_fill, Cd1, Ca1, Cp1, AxCd1, AxCa1, AxCp1, JAxCd1, JAxCa1, JAxCp1, &
                       Cd2, Ca2, Cp2, AxCd2, AxCa2, AxCp2, JAxCd2, JAxCa2, JAxCp2
      
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
      
   
   END IF

END SUBROUTINE WriteSummaryFile


        
      
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
      if ( EqualRealNos(memLength, 0.0_ReKi) )then
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
   allocate ( nodes(maxNodes), STAT = errStat )
      if ( errStat /= 0 ) then
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
SUBROUTINE SetDepthBasedCoefs( z, tMG, NCoefDpth, CoefDpths, Cd, Ca, Cp, AxCd, AxCa, AxCp )
   
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
   else
      Cd     = CoefDpths(indx1)%DpthCd*(1-s)     + CoefDpths(indx2)%DpthCd*s
      Ca     = CoefDpths(indx1)%DpthCa*(1-s)     + CoefDpths(indx2)%DpthCa*s
      Cp     = CoefDpths(indx1)%DpthCp*(1-s)     + CoefDpths(indx2)%DpthCp*s
      AxCd   = CoefDpths(indx1)%DpthCd*(1-s)     + CoefDpths(indx2)%DpthAxCd*s
      AxCa   = CoefDpths(indx1)%DpthCa*(1-s)     + CoefDpths(indx2)%DpthAxCa*s
      AxCp   = CoefDpths(indx1)%DpthCp*(1-s)     + CoefDpths(indx2)%DpthAxCp*s
   end if
   

END SUBROUTINE SetDepthBasedCoefs



!====================================================================================================
!SUBROUTINE SetExternalHydroCoefs
SUBROUTINE SetExternalHydroCoefs(  MSL2SWL, MCoefMod, MmbrCoefIDIndx, SimplCd, SimplCdMG, SimplCa, SimplCaMG, SimplCp, &
                                   SimplCpMG, SimplAxCd, SimplAxCdMG, SimplAxCa, SimplAxCaMG, SimplAxCp, SimplAxCpMG, CoefMembers,    &
                                   NCoefDpth, CoefDpths, numNodes, nodes, member )   
!     This private subroutine generates the Cd, Ca, Cp, CdMG, CaMG and CpMG coefs for the member based on
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
   type(Morison_CoefMembers), allocatable, intent(in   )  :: CoefMembers(:)
   integer(IntKi),                         intent(in   )  :: NCoefDpth
   type(Morison_CoefDpths),   allocatable, intent(in   )  :: CoefDpths(:)
   integer(IntKi),                         intent(in   )  :: numNodes
   type(Morison_NodeType),    allocatable, intent(in   )  :: nodes(:)
   type(Morison_MemberType),               intent(inout)  :: member
   
   type(Morison_NodeType)                      :: node, node1, node2
   integer(IntKi)                              :: i, j
   real(ReKi)                                  :: s, Cd, CdMG, Ca, CaMG, Cp, CpMG, AxCa, AxCp, AxCaMG, AxCpMG
  
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
         else
            member%Cd    (i) = SimplCd
            member%Ca    (i) = SimplCa
            member%Cp    (i) = SimplCp
            member%AxCd  (i) = SimplAxCd
            member%AxCa  (i) = SimplAxCa
            member%AxCp  (i) = SimplAxCp
         end if
      end do
         
   CASE (2) ! Depth-based model: coefficients are set using depth-based table data
      do i = 1, member%NElements + 1
         CALL SetDepthBasedCoefs( nodes(member%NodeIndx(i))%Position(3)+MSL2SWL,  member%tMG(i), NCoefDpth, CoefDpths, member%Cd(i), member%Ca(i), &
                                    member%Cp(i), member%AxCd(i), member%AxCa(i), member%AxCp(i) )
      end do
         
   CASE (3) ! Member-based model: coefficients set using member-specific coefficient tables
       do i = 1, member%NElements + 1
         ! Pull member  end-node data from the tables and then linearly interpolate it onto the interior member nodes    
         s = (real(i,ReKi)-1.0) / real(member%NElements,ReKi)
         if ( member%tMG(i) > 0.0_ReKi ) then
            member%Cd    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCdMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCdMG2*s
            member%Ca    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCaMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCaMG2*s
            member%Cp    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCpMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCpMG2*s 
            member%AxCd  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCaMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCdMG2*s
            member%AxCa  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCaMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCaMG2*s
            member%AxCp  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCpMG1*(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCpMG2*s
         else
            member%Cd    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCd1 *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCd2 *s
            member%Ca    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCa1 *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCa2 *s
            member%Cp    (i) = CoefMembers(MmbrCoefIDIndx)%MemberCp1 *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberCp2 *s
            member%AxCd  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCd1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCd2  *s
            member%AxCa  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCa1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCa2  *s
            member%AxCp  (i) = CoefMembers(MmbrCoefIDIndx)%MemberAxCp1  *(1-s) + CoefMembers(MmbrCoefIDIndx)%MemberAxCp2  *s
         end if
      end do
   end select
  
end subroutine SetExternalHydroCoefs


SUBROUTINE SetNodeMG( numMGDepths, MGDepths, node, MSL2SWL, tMG, MGdensity )
   ! sets the margine growth thickness of a single node (previously all nodes)
   INTEGER,                                  INTENT( IN    )  :: numMGDepths
   TYPE(Morison_MGDepthsType), ALLOCATABLE,  INTENT( IN    )  :: MGDepths(:)
   TYPE(Morison_NodeType),                   INTENT( IN    )  :: node
   real(ReKi),                               intent( in    )  :: MSL2SWL
   real(ReKi),                               intent( inout )  :: tMG
   real(ReKi),                               intent( inout )  :: MGdensity
   
   INTEGER                 :: I, J
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
   call AllocAry(member%NodeIndx     , member%NElements,   'member%NodeIndx'     , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%dRdl_mg      , member%NElements,   'member%dRdl_mg'      , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
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
   call AllocAry(member%Rin          , member%NElements+1, 'member%Rin          ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%tMG          , member%NElements+1, 'member%tMG          ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%MGdensity    , member%NElements+1, 'member%MGdensity    ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%Cd           , member%NElements+1, 'member%Cd           ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%Ca           , member%NElements+1, 'member%Ca           ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%Cp           , member%NElements+1, 'member%Cp           ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%AxCd         , member%NElements+1, 'member%AxCd         ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%AxCa         , member%NElements+1, 'member%AxCa         ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry(member%AxCp         , member%NElements+1, 'member%AxCp         ', errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_D    , 6, member%NElements+1, 'memberLoads%F_D'   , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_A    , 6, member%NElements+1, 'memberLoads%F_A'   , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_B    , 6, member%NElements+1, 'memberLoads%F_B'   , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_BF   , 6, member%NElements+1, 'memberLoads%F_BF'  , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_I    , 6, member%NElements+1, 'memberLoads%F_I'   , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_If   , 6, member%NElements+1, 'memberLoads%F_If'  , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_WMG  , 6, member%NElements+1, 'memberLoads%F_WMG' , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
   call AllocAry( memberLoads%F_IMG  , 6, member%NElements+1, 'memberLoads%F_IMG' , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)

   ! Initialize everything to zero
   member%NodeIndx      = 0.0_ReKi
   member%dRdl_mg       = 0.0_ReKi
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
   member%Rin           = 0.0_ReKi
   member%tMG           = 0.0_ReKi
   member%MGdensity     = 0.0_ReKi
   member%Cd            = 0.0_ReKi
   member%Ca            = 0.0_ReKi
   member%Cp            = 0.0_ReKi
   member%AxCd          = 0.0_ReKi
   member%AxCa          = 0.0_ReKi
   member%AxCp          = 0.0_ReKi
   memberLoads%F_D      = 0.0_ReKi
   memberLoads%F_A      = 0.0_ReKi
   memberLoads%F_B      = 0.0_ReKi
   memberLoads%F_BF     = 0.0_ReKi
   memberLoads%F_I      = 0.0_ReKi
   memberLoads%F_If     = 0.0_ReKi
   memberLoads%F_WMG    = 0.0_ReKi
   memberLoads%F_IMG    = 0.0_ReKi

end subroutine AllocateMemberDataArrays

subroutine FlipMemberNodeData( member, nodes, doSwap, errStat, errMsg )
   type(Morison_MemberType),     intent (inout)  :: member
   type(Morison_NodeType),       intent (in   )  :: nodes(:)
   logical,                      intent (  out)  :: doSwap
   integer(IntKi),               intent (  out)  :: errStat              ! returns a non-zero value when an error occurs            
   character(*),                 intent (  out)  :: errMsg               ! Error message if errStat /= ErrID_None
   
   integer(IntKi) :: i, j1, j2, numMemNodes, indx
   
   errStat = ErrID_None
   errMSg  = ''
   
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
   type(Morison_NodeType) :: node1, node2
   real(ReKi)     :: vec(3), vecLen
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
   phi = acos(vec(3)/memLength)  ! incline angle   
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
                                   InitInp%SimplCpMG, InitInp%SimplAxCd, InitInp%SimplAxCdMG, InitInp%SimplAxCa, InitInp%SimplAxCaMG, InitInp%SimplAxCp, InitInp%SimplAxCpMG, InitInp%CoefMembers,    &
                                   InitInp%NCoefDpth, InitInp%CoefDpths, InitInp%NNodes, InitInp%Nodes, member )       
   
   ! calculate reference incline angle and heading, and related trig values.  Note: members are straight to start
   Za = InitInp%Nodes(member%NodeIndx(1  ))%Position(3) 
   Zb = InitInp%Nodes(member%NodeIndx(N+1))%Position(3)

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
      if ( abs(Zb) < abs(member%Rmg(N+1)*sinPhi) ) then
         call SetErrStat(ErrID_Fatal, 'The upper end-plate of a member must not cross the water plane.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, 'SetMemberProperties' )   
      end if
      if ( abs(Za) < abs(member%Rmg(1)*sinPhi) ) then
         call SetErrStat(ErrID_Fatal, 'The lower end-plate of a member must not cross the water plane.  This is not true for Member ID '//trim(num2lstr(member%MemberID)), errStat, errMsg, 'SetMemberProperties' )   
      end if
      
      if ( ( Za < -WtrDepth .and. Zb >= -WtrDepth ) .and. ( phi > 10.0*d2r .or. abs((member%RMG(N+1) - member%RMG(i))/member%RefLength)>0.1 ) ) then
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
      member%dRdl_mg(i) = (member%RMG(i+1) - member%RMG(i))/dl
      member%dRdl_in(i) = (member%Rin(i+1) - member%Rin(i))/dl
      
      member%alpha(   i) = GetAlpha(member%RMG(i), member%RMG(i+1))
      member%alpha_fb(i) = GetAlpha(member%Rin(i), member%Rin(i+1))
      
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
            CALL FloodedBallastPartSegment(member%Rin(i+1), Rmidin, -Lmid, 0.0, Vballast_u, member%m_fb_u(i), member%h_cfb_u(i), member%I_lfb_u(i), member%I_rfb_u(i))   ! get precomputed quantities for upper half-segment
 
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
         else if (0.0 > Zb) then 
            ! fully submerged elements.  
            ! NOTE: For an element which is fractionaly in the seabed, the entire element volume is added to totals
            member%Vinner = member%Vinner + Vinner_l + Vinner_u
            member%Vouter = member%Vouter + Vouter_l + Vouter_u
            member%Vsubmerged = member%Vsubmerged + Vouter_l + Vouter_u
         else if ((0.0 > Za) .AND. (0.0 <= Zb)) then
            if (i == 1) then
               call SetErrStat(ErrID_Fatal, 'The lowest element of a member must not cross the free surface.  This is true for MemberID '//trim(num2lstr(member%MemberID)), errStat, errMsg, 'SetMemberProperties')
            end if
            
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
   ALLOCATE ( p%Members(p%NMembers), STAT = errStat )
   IF ( errStat /= ErrID_None ) THEN
      errMsg  = ' Error allocating space for the members array.'
      errStat = ErrID_Fatal
      RETURN
   END IF   
   
   ALLOCATE ( m%MemberLoads(p%NMembers), STAT = errStat )
   IF ( errStat /= ErrID_None ) THEN
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
      
      call AllocateMemberDataArrays(p%Members(i), m%MemberLoads(i), errStat2, errMsg2) ; call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'SetupMembers')
        
      p%Members(i)%NodeIndx  = InitInp%InpMembers(i)%NodeIndx ! now that the parameter version is allocated, copy the data from the InitInp version
      
      ! only reorder the nodes if the end nodes do not follow the necessary coordinate ordering rules
      call FlipMemberNodeData(p%Members(i), InitInp%nodes, doSwap, errStat2, errMsg2) ; call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'SetupMembers')
      if (doSwap) then
            prop2Indx = InitInp%InpMembers(I)%MPropSetID1Indx
            prop1Indx = InitInp%InpMembers(I)%MPropSetID2Indx
      else
            prop1Indx = InitInp%InpMembers(I)%MPropSetID1Indx
            prop2Indx = InitInp%InpMembers(I)%MPropSetID2Indx
      end if
      ! Now populate the various member data arrays using the HydroDyn input file data
      call SetMemberProperties( InitInp%MSL2SWL, InitInp%Gravity, p%Members(i), InitInp%InpMembers(i)%MCoefMod, InitInp%InpMembers(i)%MmbrCoefIDIndx, InitInp%InpMembers(i)%MmbrFilledIDIndx, InitInp%MPropSets(prop1Indx), InitInp%MPropSets(prop2Indx), InitInp, errStat2, errMsg2 ) ; call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'SetupMembers')
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
   REAL(DbKi),                        INTENT(INOUT)  :: Interval    !< Coupling interval in seconds: the rate that 
                                                                     !!   (1) Morison_UpdateStates() is called in loose coupling &
                                                                     !!   (2) Morison_UpdateDiscState() is called in tight coupling.
                                                                     !!   Input is the suggested time from the glue code; 
                                                                     !!   Output is the actual coupling interval that will be used 
                                                                     !!   by the glue code.
   TYPE(Morison_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
   INTEGER(IntKi),                    INTENT(  OUT)  :: errStat     !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: errMsg      !< Error message if errStat /= ErrID_None

   TYPE(Morison_MemberType) :: member    ! the current member
   type(Morison_MemberInputType) :: inpMember ! current input file-based member
   INTEGER                  :: N, i, j, count
   REAL(ReKi)               :: dl
   REAL(ReKi)               :: vec(3),v2D(3,1), pos(3)
   REAL(ReKi)               :: phi    ! member tilt angle
   REAL(ReKi)               :: beta   ! member tilt heading
   REAL(ReKi)               :: cosPhi
   REAL(ReKi)               :: sinPhi
   REAL(ReKi)               :: tanPhi
   REAL(ReKi)               :: sinBeta
   REAL(ReKi)               :: cosBeta
   REAL(ReKi)               :: Za
   REAL(ReKi)               :: Zb
   real(ReKi)               :: memLength ! reference member length
   real(ReKi)               :: An(3), An_drag(3), Vn(3), I_n(3), Z0, sgn, Amag, Amag_drag, Vmag, Imag, Ir_MG_end, Il_MG_end, R_I(3,3), IRl_mat(3,3), tMG, MGdens, F_I(3), F_DP(3), af(3), VnDotAf
   integer(IntKi)           :: MemberEndIndx, ncommon
   INTEGER, ALLOCATABLE       :: commonNodeLst(:)
   LOGICAL, ALLOCATABLE       :: usedJointList(:)
   integer(IntKi)           :: errStat2              ! returns a non-zero value when an error occurs            
   CHARACTER(errMsgLen)     :: errMsg2     ! Error message if errStat2 /= ErrID_None

      
      ! Initialize errStat        
   errStat = ErrID_None         
   errMsg  = ""               
      
  
      
      ! Initialize the NWTC Subroutine Library         
   CALL NWTC_Init(  )

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
   p%OutSwtch   = InitInp%OutSwtch
   p%MSL2SWL    = InitInp%MSL2SWL
   
   ALLOCATE ( p%MOutLst(p%NMOutputs), STAT = errStat )
   IF ( errStat /= ErrID_None ) THEN
      errMsg  = ' Error allocating space for MOutLst array.'
      errStat = ErrID_Fatal
      RETURN
   END IF
   IF (ALLOCATED(InitInp%MOutLst) ) &
      p%MOutLst =    InitInp%MOutLst           ! Member output data
      
   p%NJOutputs = InitInp%NJOutputs                        ! Number of joints to output [ >=0 and <10]
      
   ALLOCATE ( p%JOutLst(p%NJOutputs), STAT = errStat )
   IF ( errStat /= ErrID_None ) THEN
      errMsg  = ' Error allocating space for JOutLst array.'
      errStat = ErrID_Fatal
      RETURN
   END IF
   IF (ALLOCATED(InitInp%JOutLst) ) &
      p%JOutLst =    InitInp%JOutLst            ! Joint output data
 
   ! ----------------------- set up the members -----------------------
   call SetupMembers( InitInp, p, m, errStat2, errMsg2 ) ; call SetErrStat( errStat2, errMsg2, errStat, errMsg, 'Morison_Init' )
   if ( errStat >= AbortErrLev ) return
   
   !------------------------ set up joint (or joint-node) properties --
   do i = 1, InitInp%NJoints
      InitInp%Nodes(i)%JAxCd = InitInp%AxialCoefs(InitInp%InpJoints(i)%JointAxIDIndx)%AxCd
      InitInp%Nodes(i)%JAxCa = InitInp%AxialCoefs(InitInp%InpJoints(i)%JointAxIDIndx)%AxCa
      InitInp%Nodes(i)%JAxCp = InitInp%AxialCoefs(InitInp%InpJoints(i)%JointAxIDIndx)%AxCp
      InitInp%Nodes(i)%JAxCd = InitInp%AxialCoefs(InitInp%InpJoints(i)%JointAxIDIndx)%AxCd
      InitInp%Nodes(i)%JAxCa = InitInp%AxialCoefs(InitInp%InpJoints(i)%JointAxIDIndx)%AxCa
      InitInp%Nodes(i)%JAxCp = InitInp%AxialCoefs(InitInp%InpJoints(i)%JointAxIDIndx)%AxCp  
      ! Redundant work (these are already assigned to the member data arrays, 
      ! but is needed on the joint data because we report the tMG, and MGDensity at each Joint node in the Summary File
      call SetNodeMG( InitInp%NMGDepths, InitInp%MGDepths, InitInp%Nodes(i), InitInp%MSL2SWL, InitInp%Nodes(i)%tMG, InitInp%Nodes(i)%MGDensity )
   end do

      ! allocate and copy in node-based load and hydrodynamic arrays
   call AllocateNodeLoadVariables(InitInp, p, m, p%NNodes, errStat, errMsg )
   call MOVE_ALLOC( InitInp%nodeInWater, p%nodeInWater )   

 
   
   ! Create the input and output meshes associated with loads at the nodes
      
   CALL MeshCreate( BlankMesh      = u%Mesh          &
                     ,IOS          = COMPONENT_INPUT        &
                     ,Nnodes       = p%NNodes      &
                     ,errStat      = errStat                &
                     ,ErrMess      = errMsg                 &
                     ,TranslationDisp = .TRUE.              &
                     ,Orientation     = .TRUE.              &
                     ,TranslationVel  = .TRUE.              &
                     ,RotationVel     = .TRUE.              &
                     ,TranslationAcc  = .TRUE.              &
                     ,RotationAcc     = .TRUE.               )

   IF ( errStat >= AbortErrLev ) RETURN

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
                        , errStat                  &
                        , errMsg                   &
                        ) !, transpose(p%Nodes(I)%R_LToG)          )
      IF ( errStat /= 0 ) RETURN

!TODO: Do we still need this for visualization?  How is it used? GJH 3/26/2020  Actually need a line mesh to properly visualize the members
     ! Morison_Rad(count) = p%Nodes(I)%R   ! set this for FAST visualization
      
     
   
         ! Create the mesh element
   
      CALL MeshConstructElement (u%Mesh   &
                            , ELEMENT_POINT      &                                  
                            , errStat            &
                            , errMsg  &
                            , i                  &
                                        )

   END DO

   CALL MeshCommit ( u%Mesh   &
                      , errStat            &
                      , errMsg             )
   
   IF ( errStat /= 0 ) THEN
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
                     ,errStat      = errStat                &
                     ,ErrMess      = errMsg                 &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 )
   u%Mesh%RemapFlag = .TRUE.
   y%Mesh%RemapFlag = .TRUE.

         ! Define initial system states here:

   x%DummyContState           = 0
   xd%DummyDiscState          = 0
   z%DummyConstrState         = 0
   OtherState%DummyOtherState = 0
   m%LastIndWave              = 1

   ! IF ( p%OutSwtch > 0 ) THEN  @mhall: I think the below need to be allocated in all cases


   
   ! allocate and initialize joint-specific arrays   
      
   ALLOCATE ( commonNodeLst(10), STAT = errStat )
   IF ( errStat /= ErrID_None ) THEN
      errMsg  = ' Error allocating space for the commonNodeLst array.'
      errStat = ErrID_Fatal
      RETURN
   END IF 
   commonNodeLst = -1
   
   ALLOCATE ( usedJointList(p%NJoints), STAT = errStat )
   IF ( errStat /= ErrID_None ) THEN
      errMsg  = ' Error allocating space for the UsedJointList array.'
      errStat = ErrID_Fatal
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
            p%AM_End(:,:,i) = (InitInp%Nodes(I)%JAxCa*InitInp%WtrDens/ Vmag)*matmul(transpose(v2D), v2D) 
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
   
   END DO ! looping through nodes that are joints, i
          
         ! Define initial guess for the system inputs here:
         !    u%DummyInput = 0
         ! Define system output initializations (set up mesh) here:  
         ! Define initialization-routine output here:
         
         ! Initialize the outputs      
   IF ( p%OutSwtch > 0) then  !@mhall: moved this "if" to after allocations
   
      CALL MrsnOUT_Init( InitInp, y, p, InitOut, errStat, errMsg )
      IF ( errStat > AbortErrLev ) RETURN
      
         ! Determine if we need to perform output file handling
      
      IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN  
         CALL MrsnOUT_OpenOutput( Morison_ProgDesc%Name, TRIM(InitInp%OutRootName)//'.HD', p, InitOut, errStat, errMsg )
         IF ( errStat > AbortErrLev ) RETURN
      END IF
      
   END IF  
   
   ! We will call CalcOutput to compute the loads for the initial reference position
   ! Then we can use the computed load components in the Summary File
   ! NOTE: Morison module has no states, otherwise we could no do this. GJH
   
   call Morison_CalcOutput(0.0_DbKi, u, p, x, xd, z, OtherState, y, m, errStat, errMsg )
   
      ! Write Summary information now that everything has been initialized. 
   CALL WriteSummaryFile( InitInp%UnSum, InitInp%Gravity, InitInp%MSL2SWL, InitInp%WtrDpth, InitInp%NJoints, InitInp%NNodes, InitInp%Nodes, p%NMembers, p%Members, &
                          p%NumOuts, p%OutParam, p%NMOutputs, p%MOutLst,  p%NJOutputs, p%JOutLst, u%Mesh, y%Mesh, &
                          p, m, errStat, errMsg )
   IF ( errStat > AbortErrLev ) RETURN  
                                                       
   !Contains:
   !   SUBROUTINE CleanUpInitOnErr
   !   IF (ALLOCATED(sw(1)%array))  DEALLOCATE(sw(1)%array, STAT=aviFail)
   !   END SUBROUTINE

END SUBROUTINE Morison_Init


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


FUNCTION GetAlpha(R1,R2)
   ! calculates relative center of volume location for a (tapered) cylindrical element
   real(ReKi)    :: GetAlpha
   REAL(ReKi),                     INTENT    ( IN    )  :: R1  ! interior radius of element at node point
   REAL(ReKi),                     INTENT    ( IN    )  :: R2  ! interior radius of other end of part-element
   
      
   GetAlpha = (R1*R1 + 2.0*R1*R2 + 3.0*R2*R2)/4.0/(R1*R1 + R1*R2 + R2*R2)

   
END FUNCTION GetAlpha


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
   call AllocAry( m%FDynP        ,       NNodes   , 'm%FDynP'        , errStat2, errMsg2); call SetErrStat(errStat2, errMsg2, errStat, errMsg, routineName)
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
   if (errStat == ErrID_Fatal) return
   
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
   
   allocate( p%WaveVel(0:p%NStepWave, p%NNodes, 3), STAT = errStat )
   IF ( errStat /= ErrID_None ) THEN
      errMsg  = ' Error allocating space for wave velocities array.'
      errStat = ErrID_Fatal
      RETURN
   END IF
   p%WaveVel = InitInp%WaveVel      
      
   allocate( p%WaveAcc(0:p%NStepWave, p%NNodes, 3), STAT = errStat )
   IF ( errStat /= ErrID_None ) THEN
      errMsg  = ' Error allocating space for wave accelerations array.'
      errStat = ErrID_Fatal
      RETURN
   END IF
   p%WaveAcc = InitInp%WaveAcc
      
   allocate( p%WaveDynP(0:p%NStepWave, p%NNodes), STAT = errStat )
   IF ( errStat /= ErrID_None ) THEN
      errMsg  = ' Error allocating space for wave dynamic pressure array.'
      errStat = ErrID_Fatal
      RETURN
   END IF
   p%WaveDynP = InitInp%WaveDynP
      
   allocate( p%WaveTime(0:p%NStepWave), STAT = errStat )
   IF ( errStat /= ErrID_None ) THEN
      errMsg  = ' Error allocating space for wave time array.'
      errStat = ErrID_Fatal
      RETURN
   END IF     
   p%WaveTime     = InitInp%WaveTime   
   
   
END SUBROUTINE AllocateNodeLoadVariables



!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE Morison_End( u, p, x, xd, z, OtherState, y, m, errStat, errMsg )
!..................................................................................................................................

      TYPE(Morison_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(Morison_ParameterType),       INTENT(INOUT)  :: p           !< Parameters     
      TYPE(Morison_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(Morison_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(Morison_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(Morison_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states            
      TYPE(Morison_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(Morison_MiscVarType),         INTENT(INOUT)  :: m           !< Misc/optimization variables            
      INTEGER(IntKi),                    INTENT(  OUT)  :: errStat     !< Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: errMsg      !< Error message if errStat /= ErrID_None



         ! Initialize errStat
         
      errStat = ErrID_None         
      errMsg  = ""               
      
      
         ! Place any last minute operations or calculations here:


         ! Close files here:     
                  
 

         ! Destroy the input data:
         
      CALL Morison_DestroyInput( u, errStat, errMsg )


         ! Determine if we need to close the output file
         
      IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN   
         CALL MrsnOut_CloseOutput( p, errStat, errMsg )         
      END IF 
         
         ! Destroy the parameter data:
         
      
      CALL Morison_DestroyParam( p, errStat, errMsg )


         ! Destroy the state data:
         
      CALL Morison_DestroyContState(   x,           errStat, errMsg )
      CALL Morison_DestroyDiscState(   xd,          errStat, errMsg )
      CALL Morison_DestroyConstrState( z,           errStat, errMsg )
      CALL Morison_DestroyOtherState(  OtherState,  errStat, errMsg )
         
      CALL Morison_DestroyMisc( m, errStat, errMsg )

         ! Destroy the output data:
         
      CALL Morison_DestroyOutput( y, errStat, errMsg )


      

END SUBROUTINE Morison_End
!----------------------------------------------------------------------------------------------------------------------------------
!> This is a loose coupling routine for solving constraint states, integrating continuous states, and updating discrete and other 
!! states. Continuous, constraint, discrete, and other states are updated to values at t + Interval.
SUBROUTINE Morison_UpdateStates( Time, u, p, x, xd, z, OtherState, m, errStat, errMsg )
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
      INTEGER(IntKi),                     INTENT(  OUT) :: errStat     !< Error status of the operation     
      CHARACTER(*),                       INTENT(  OUT) :: errMsg      !< Error message if errStat /= ErrID_None

         ! Local variables
                  
      INTEGER(IntKi)                                    :: errStat2    ! Error status of the operation (occurs after initial error)
      CHARACTER(errMsgLen)                              :: errMsg2     ! Error message if errStat2 /= ErrID_None
                        
         ! Initialize errStat
         
      errStat = ErrID_None         
      errMsg  = ""               
      
           
      

      
END SUBROUTINE Morison_UpdateStates
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
!> Use in conjunction with GetInterpolationSlope, to replace InterpWrappedStpReal here.
FUNCTION InterpolateWithSlope(InterpSlope, Ind, YAry)
      REAL(ReKi), INTENT(IN)                            :: InterpSlope
      INTEGER(IntKi), INTENT(IN )                       :: Ind           !< Misc/optimization variables
      REAL(SiKi), INTENT(IN)                            :: YAry(0:)
      REAL(ReKi)                                        :: InterpolateWithSlope

      InterpolateWithSlope = ( YAry(Ind+1) - YAry(Ind) )*InterpSlope + YAry(Ind)

END FUNCTION InterpolateWithSlope
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
      
   REAL(ReKi)                                        :: F_DP(6), kvec(3), v(3),  vf(3), vrel(3), vmag
   INTEGER                                           :: I, J, K, nodeIndx, IntWrapIndx
   REAL(ReKi)                                        :: AllOuts(MaxMrsnOutputs)
   REAL(ReKi)                                        :: qdotdot(6) ,qdotdot2(3)     ! The structural acceleration of a mesh node
   !REAL(ReKi)                                        :: accel_fluid(6) ! Acceleration of fluid at the mesh node
   REAL(ReKi)                                        :: dragFactor     ! The lumped drag factor
   REAL(ReKi)                                        :: AnProd         ! Dot product of the directional area of the joint
   REAL(ReKi)                                        :: C(3,3)
   REAL(ReKi)                                        :: sgn
   REAL(ReKi)                                        :: D_AM_M(6,6)
   REAL(ReKi)                                        :: nodeInWater
   REAL(ReKi)                                        :: D_dragConst     ! The distributed drag factor
   REAL(ReKi)                                        :: InterpolationSlope 


      
   TYPE(Morison_MemberType) :: mem     ! the current member
   INTEGER                  :: N       ! Number of elements within a given member
   REAL(ReKi)               :: dl      ! Element length within a given member, m
   REAL(ReKi)               :: vec(3)  ! Vector pointing from a member's 1st node to its last node
   REAL(ReKi)               :: phi, phi1, phi2     ! member tilt angle
   REAL(ReKi)               :: beta    ! member tilt heading
   real(ReKi)               :: vecLen  ! distance between member end nodes (joints) [this should never be zero but we test for it just in case]
   REAL(ReKi)               :: cosPhi, cosPhi1, cosPhi2
   REAL(ReKi)               :: sinPhi, sinPhi1, sinPhi2
   REAL(ReKi)               :: tanPhi
   REAL(ReKi)               :: sinBeta, sinBeta1, sinBeta2
   REAL(ReKi)               :: cosBeta, cosBeta1, cosBeta2
   real(ReKi)               :: CMatrix(3,3), CTrans(3,3) ! Direction cosine matrix for element, and its transpose
   REAL(ReKi)               :: z1
   REAL(ReKi)               :: z2
   REAL(ReKi)               :: r1
   REAL(ReKi)               :: r2
   real(ReKi)               :: p1(3), p2(3)
   REAL(ReKi)               :: dRdl_mg     ! shorthand for taper including marine growth of element i
   REAL(ReKi)               :: Rmid  
   REAL(ReKi)               :: RmidMG
   REAL(ReKi)               :: Rmidin
   REAL(ReKi)               :: Lmid  
   real(ReKi)               :: g     ! gravity constant
   REAL(ReKi)               :: h0    ! distances along cylinder centerline from point 1 to the waterplane
   real(ReKi)               :: k_hat(3), k_hat1(3), k_hat2(3) ! Elemental unit vector pointing from 1st node to 2nd node of the element
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
   REAL(ReKi)               :: Moment     !moment induced about the center of the cylinder's bottom face
   REAL(ReKi)               :: BuoyF(3) ! buoyancy force vector aligned with an element
   REAL(ReKi)               :: BuoyM(3) ! buoyancy moment vector aligned with an element
   integer(IntKi)           :: im    ! counter   
   real(ReKi)               :: a_s1(3)       
   real(ReKi)               :: alpha_s1(3)
   real(ReKi)               :: omega_s1(3)
   real(ReKi)               :: a_s2(3)       
   real(ReKi)               :: alpha_s2(3)
   real(ReKi)               :: omega_s2(3)
   real(ReKi)               :: pos1(3), pos2(3)   
   real(ReKi)               :: Imat(3,3)
   real(ReKi)               :: iArm(3), iTerm(3), Ioffset, h_c, dRdl_p, dRdl_pp, f_hydro(3), Am(3,3), lstar, deltal
   real(ReKi)               :: C_1, C_2, a0b0, z1d, z2d, h
   real(ReKi)               :: F_WMG(6), F_IMG(6), F_If(6), F_A(6), F_I(6), F_D(6), F_B1(6), F_B2(6)

      ! Initialize errStat
         
   errStat = ErrID_None         
   errMsg  = ""               
   Imat    = 0.0_ReKi   
   g       = p%Gravity
   
   InterpolationSlope = GetInterpolationSlope(Time, p, m, IntWrapIndx)

   !===============================================================================================
   ! Calculate the fluid kinematics at all mesh nodes and store for use in the equations below
   
   do j = 1, p%NNodes
      m%nodeInWater(j) = REAL( p%nodeInWater(IntWrapIndx,j), ReKi )
      
         ! Determine the dynamic pressure at the node
      m%FDynP(j) = InterpolateWithSlope(InterpolationSlope, m%LastIndWave, p%WaveDynP(:,j))
      do i=1,3
            ! Determine the fluid acceleration and velocity and relative structural velocity at the node
         m%FA(i,j) = InterpolateWithSlope(InterpolationSlope, m%LastIndWave, p%WaveAcc(:,j,i)) 
               
         m%FV(i,j) = InterpolateWithSlope(InterpolationSlope, m%LastIndWave, p%WaveVel(:,j,i)) 
         m%vrel(i,j) = m%FV(i,j) - u%Mesh%TranslationVel(i,j)
      end do
   end do
   
   ! ==============================================================================================
   ! Calculate instantaneous loads on each member except for the hydrodynamic loads on member ends.
   ! This covers aspects of the load calculations previously in CreateDistributedMesh.  

   ! Zero out previous time-steps loads (these are loads which are computed at the member-level and summed onto a node, 
   !    so they need to be zeroed out before the summations happen)
   !m%F_WMG   = 0.0_ReKi
   !m%F_IMG   = 0.0_ReKi
   m%F_BF_End= 0.0_ReKi
   !m%F_If    = 0.0_ReKi
   !m%F_D     = 0.0_ReKi
   !m%F_A     = 0.0_ReKi
   !m%F_I     = 0.0_ReKi
   !m%F_B     = 0.0_ReKi
   !m%F_BF    = 0.0_ReKi
   m%F_B_End = 0.0_ReKi
   y%Mesh%Force  = 0.0_ReKi
   y%Mesh%Moment = 0.0_ReKi
   
   ! Loop through each member
   DO im = 1, p%NMembers    
      N   = p%Members(im)%NElements      
      mem = p%Members(im)   !@mhall: does this have much overhead?

      !zero member loads
      m%memberLoads(im)%F_B = 0.0_ReKi
      m%memberLoads(im)%F_BF = 0.0_ReKi
      m%memberLoads(im)%F_D = 0.0_ReKi
      m%memberLoads(im)%F_A = 0.0_ReKi
      m%memberLoads(im)%F_I = 0.0_ReKi
      m%memberLoads(im)%F_WMG = 0.0_ReKi
      m%memberLoads(im)%F_IMG = 0.0_ReKi
      m%memberLoads(im)%F_If = 0.0_ReKi

      DO i =1,N    ! loop through member elements

            ! calculate isntantaneous incline angle and heading, and related trig values
         ! the first and last NodeIndx values point to the corresponding Joint nodes idices which are at the start of the Mesh

         pos1    = u%Mesh%TranslationDisp(:, mem%NodeIndx(i))   + u%Mesh%Position(:, mem%NodeIndx(i)) 
         pos1(3) = pos1(3) - p%MSL2SWL
         pos2    = u%Mesh%TranslationDisp(:, mem%NodeIndx(i+1)) + u%Mesh%Position(:, mem%NodeIndx(i+1)) 
         pos2(3) = pos2(3) - p%MSL2SWL

         call GetOrientationAngles( pos1, pos2, phi, sinPhi, cosPhi, tanPhi, sinBeta, cosBeta, k_hat, errStat2, errMsg2 )
         call Morison_DirCosMtrx( pos1, pos2, CMatrix )
         CTrans  = transpose(CMatrix)
         ! save some commonly used variables   
         dl      = mem%dl
         z1      = pos1(3)   ! get node z locations from input mesh
         z2      = pos2(3)
         r1      = mem%RMG(i  )                         ! outer radius element nodes including marine growth
         r2      = mem%RMG(i+1)
         dRdl_mg = mem%dRdl_mg(i)                                    ! Taper of element including marine growth
         a_s1    = u%Mesh%TranslationAcc(:, mem%NodeIndx(i  ))
         alpha_s1= u%Mesh%RotationAcc   (:, mem%NodeIndx(i  ))
         omega_s1= u%Mesh%RotationVel   (:, mem%NodeIndx(i  ))
         a_s2    = u%Mesh%TranslationAcc(:, mem%NodeIndx(i+1))
         alpha_s2= u%Mesh%RotationAcc   (:, mem%NodeIndx(i+1))
         omega_s2= u%Mesh%RotationVel   (:, mem%NodeIndx(i+1))
        
         if ( .not. mem%PropPot )  then ! Member is NOT modeled with Potential Flow Theory
         
            ! should i_floor theshold be applied to below calculations to avoid wasting time on computing zero-valued things? <<<<<
            ! should lumped half-element coefficients get combined at initialization? <<<
              
            ! ------------------ marine growth: Sides: Section 4.1.2 --------------------  
            F_WMG = 0.0_ReKi

            ! lower node
            !m%F_WMG(3, mem%NodeIndx(i  )) = m%F_WMG(3, mem%NodeIndx(i  )) - mem%m_mg_l(i)*g ! weight force  : Note: this is a constant
            !m%F_WMG(4, mem%NodeIndx(i  )) = m%F_WMG(4, mem%NodeIndx(i  )) - mem%m_mg_l(i)*g * mem%h_cmg_l(i)* sinPhi * sinBeta! weight force
            !m%F_WMG(5, mem%NodeIndx(i  )) = m%F_WMG(5, mem%NodeIndx(i  )) + mem%m_mg_l(i)*g * mem%h_cmg_l(i)* sinPhi * cosBeta! weight force
            
            F_WMG(3) = - mem%m_mg_l(i)*g ! weight force  : Note: this is a constant
            F_WMG(4) = - mem%m_mg_l(i)*g * mem%h_cmg_l(i)* sinPhi * sinBeta! weight force
            F_WMG(5) =   mem%m_mg_l(i)*g * mem%h_cmg_l(i)* sinPhi * cosBeta! weight force
            m%memberLoads(im)%F_WMG(:,i) = m%memberLoads(im)%F_WMG(:,i) + F_WMG
            y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + F_WMG(1:3)
            y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + F_WMG(4:6)
            
            ! upper node
            !m%F_WMG(3, mem%NodeIndx(i+1)) = m%F_WMG(3, mem%NodeIndx(i+1)) - mem%m_mg_u(i)*g ! weight force  : Note: this is a constant 
            !m%F_WMG(4, mem%NodeIndx(i+1)) = m%F_WMG(4, mem%NodeIndx(i+1)) - mem%m_mg_u(i)*g * mem%h_cmg_u(i)* sinPhi * sinBeta! weight force
            !m%F_WMG(5, mem%NodeIndx(i+1)) = m%F_WMG(5, mem%NodeIndx(i+1)) + mem%m_mg_u(i)*g * mem%h_cmg_u(i)* sinPhi * cosBeta! weight force
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
            !m%F_IMG(1:3, mem%NodeIndx(i  )) = m%F_IMG(1:3, mem%NodeIndx(i  )) + iTerm
            !m%F_IMG(4:6, mem%NodeIndx(i  )) = m%F_IMG(4:6, mem%NodeIndx(i  )) &
            !                                  - cross_product(a_s1 * mem%m_mg_l(i), mem%h_cmg_l(i) * k_hat) &
            !                                  + matmul(Imat, alpha_s1)  &
            !                                  - cross_product(omega_s1,matmul(Imat,omega_s1))
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
            !m%F_IMG(1:3, mem%NodeIndx(i+1)) = m%F_IMG(1:3, mem%NodeIndx(i+1)) + iTerm
            !m%F_IMG(4:6, mem%NodeIndx(i+1)) = m%F_IMG(4:6, mem%NodeIndx(i+1)) &
            !                                  - cross_product(a_s2 * mem%m_mg_u(i), mem%h_cmg_u(i) * k_hat) &
            !                                  + matmul(Imat, alpha_s2) &
            !                                  - cross_product(omega_s2,matmul(Imat,omega_s2))
            F_IMG(1:3) = iTerm
            F_IMG(4:6) = - cross_product(a_s2 * mem%m_mg_u(i), mem%h_cmg_u(i) * k_hat) + matmul(Imat, alpha_s2) &
                         - cross_product(omega_s2,matmul(Imat,omega_s2))
            m%memberLoads(im)%F_IMG(:,i+1) = m%memberLoads(im)%F_IMG(:,i+1) + F_IMG
            y%Mesh%Force (:,mem%NodeIndx(i+1)) = y%Mesh%Force (:,mem%NodeIndx(i+1)) + F_IMG(1:3)
            y%Mesh%Moment(:,mem%NodeIndx(i+1)) = y%Mesh%Moment(:,mem%NodeIndx(i+1)) + F_IMG(4:6)

            ! ------------------- buoyancy loads: sides: Sections 3.1 and 3.2 ------------------------

!TODO: What about elements which are buried in the seabed?  This doesn't seem to be tested for
            if (z1 < 0) then   ! if segment is at least partially submerged ...
              
              
               if (z1*z2 <= 0) then ! special calculation if the slice is partially submerged
                  
                  ! Check that this is not the 1st element of the member
                  if ( i == 1 ) then
                     call SeterrStat(ErrID_Fatal, 'The lowest element of a Morison member has become partially submerged!  This is not allowed.  Please review your model and create a discretization such that even with displacements, the lowest element of a member does not become partially submerged.', errStat, errMsg, 'Morison_CalcOutput' )                  
                     return
                  end if
                  
                  h0 = -z1/cosPhi             ! distances along element centerline from point 1 to the waterplane
              
              
                  if (abs(dRdl_mg) < 0.0001) then      ! untapered cylinder case

                     Vs =    Pi*r1*r1*h0   ! volume of total submerged portion
                     if ( EqualRealNos(Vs, 0.0_ReKi) ) then
                        cx = 0.0_ReKi  ! Avoid singularity, but continue to provide the correct solution
                     else
                        cr = 0.25*r1*r1*tanPhi/h0
                        cl = 0.5*h0 + 0.125*r1*r1*tanPhi*tanPhi/h0
                        cx = cr*cosPhi + cl*sinPhi
                     end if
                    
                     !alpha0 = 0.5*h0/dl            ! force distribution between end nodes
                 
                  else       ! inclined tapered cylinder case (note I've renamed r0 to rh here!!)
                     !===================
                     !Per plan equations
                     ! NOTE:  Variable changes of Plan     vs       Code
                     !---------------------------------------------------
                     !                             V                 Vs
                     !                             a_h               a0
                     !                             b_h               b0
                     !                             x_c               cx
                     !                             h                 h0
                     !                             r1                r_MG,i
                     !                             r_c               cr
                     !                             h_c               cl
                     ! NOTE: a0 and b0 always appear as a0b0, never separately.
                     rh   = r1 + h0*dRdl_mg    ! radius of element at point where its centerline crosses the waterplane
                     C_1  = 1.0_ReKi - dRdl_mg**2 * tanPhi**2
                     ! waterplane ellipse shape
                     b0   = rh/sqrt(C_1)
                     a0   = rh/((C_1)*cosPhi)             ! simplified from what's in ConicalCalcs.ipynb
                     a0b0 = a0*b0
                     C_2  = a0b0*rh*cosPhi - r1**3
                     cl   = -(-0.75*a0b0*rh**2*cosPhi + 0.75*r1**4*C_1 + r1*C_1*C_2) / (dRdl_mg*C_1*C_2)
                     cr   = (0.75*a0b0*dRdl_mg*rh**2*sinPhi)/(C_1*C_2)
                     cx   = cr*cosPhi + cl*sinPhi 
                     Vs   = pi*(a0b0*rh*cosPhi - r1**3)/(3.0*dRdl_mg)       
                  
                     ! End per plan equations
                     !===================
                  
                     !rh = r1 + h0*dRdl_mg    ! radius of element at point where its centerline crosses the waterplane
                     !l1 = r1/dRdl_mg  ! distance from cone end to bottom node
                     !              
                     !! waterplane ellipse shape
                     !b0 = rh/sqrt(1 - dRdl_mg**2 * tanPhi**2)
                     !a0 = rh/((1 - dRdl_mg**2*tanPhi**2)*cosPhi)             ! simplified from what's in ConicalCalcs.ipynb
                     !
                     !! segment submerged volume
                     !!Vs = pi*(a0*b0*rh*cosPhi - l1**3*dRdl_mg**3)/(3*dRdl_mg) !Original code
                     !Vs = pi*(a0*b0*rh*cosPhi - r1**3)/(3*dRdl_mg)        !Plan doc
                     !
                     !! centroid of segment submerged volume (relative to bottom node)
                     !cx = -0.25*(3*a0*b0*rh*rh*(dRdl_mg**2 + 1)*cosPhi + 3.0*l1**4*dRdl_mg**4*(dRdl_mg**2*tanPhi**2 - 1) + 4*l1*dRdl_mg*(dRdl_mg**2*tanPhi**2 - 1)*(a0*b0*rh*cosPhi - 1.0*l1**3*dRdl_mg**3))*sin(phi)/(dRdl_mg*(dRdl_mg**2*tanPhi**2 - 1)*(a0*b0*rh*cosPhi - l1**3*dRdl_mg**3))
                              
                     !alpha0 = (r1*r1 + 2*r1*r2 + 3*r2**2)/4/(r1*r1 + r1*r2 + r2**2)  ! this can be precomputed
              
                  end if

                  pwr = 3
                  alpha    = (1.0-mem%alpha(i))*z1**pwr/(-mem%alpha(i)*z2**pwr + (1.0-mem%alpha(i))*z1**pwr)

                  Fb  = Vs*p%WtrDens*g       !buoyant force
                  Fr  = -Fb*sinPhi     !radial component of buoyant force
                  Fl  = Fb*cosPhi      !axial component of buoyant force
                  Moment = -Fb*cx      !This was matt's code        !moment induced about the center of the cylinder's bottom face

                  ! calculate (imaginary) bottom plate forces/moment to subtract from displacement-based values
                  Fl  = Fl  + p%WtrDens*g*z1* Pi *r1*r1        
                  Moment  = Moment  + p%WtrDens*g* sinPhi * Pi/4.0*r1**4       


                  ! reduce taper-based moment to remove (not double count) radial force distribution to each node 
                  Moment  = Moment + Fr*(1.0_ReKi-alpha)*dl
                  !call DistributeElementLoads(Fl, Fr, Moment, sinPhi, cosPhi, sinBeta, cosBeta, alpha, m%F_B(:, mem%NodeIndx(i)), m%F_B(:, mem%NodeIndx(i-1)))
                  call DistributeElementLoads(Fl, Fr, Moment, sinPhi, cosPhi, sinBeta, cosBeta, alpha, F_B1, F_B2)
                  m%memberLoads(im)%F_B(:, i) = m%memberLoads(im)%F_B(:, i) + F_B1      ! alpha
                  m%memberLoads(im)%F_B(:, i-1) = m%memberLoads(im)%F_B(:, i-1) + F_B2  ! 1-alpha
                  y%Mesh%Force (:,mem%NodeIndx(i  )) = y%Mesh%Force (:,mem%NodeIndx(i  )) + F_B1(1:3)
                  y%Mesh%Moment(:,mem%NodeIndx(i  )) = y%Mesh%Moment(:,mem%NodeIndx(i  )) + F_B1(4:6)
                  y%Mesh%Force (:,mem%NodeIndx(i-1)) = y%Mesh%Force (:,mem%NodeIndx(i-1)) + F_B2(1:3)
                  y%Mesh%Moment(:,mem%NodeIndx(i-1)) = y%Mesh%Moment(:,mem%NodeIndx(i-1)) + F_B2(4:6)
               else ! normal, fully submerged case
              
                  Fl = -2.0*Pi*dRdl_mg*p%WtrDens*g*dl*( z1*r1 + 0.5*(z1*dRdl_mg + r1*cosPhi)*dl + 1.0/3.0*(dRdl_mg*cosPhi*dl*dl) )   ! from CylinderCalculationsR1.ipynb
              
                  Fr = -Pi*p%WtrDens*g*dl*(r1*r1 + dRdl_mg*r1*dl + (dRdl_mg**2*dl**2)/3.0)*sinPhi                          ! from CylinderCalculationsR1.ipynb
                  Moment = -Pi*dl*g*p%WtrDens*(3.0*dl**3*dRdl_mg**4 + 3.0*dl**3*dRdl_mg**2 + 12.0*dl**2*dRdl_mg**3*r1 + 8.0*dl**2*dRdl_mg*r1 + 18.0*dl*dRdl_mg**2*r1*r1 + 6.0*dl*r1*r1 + 12.0*dRdl_mg*r1**3)*sinPhi/12.0   ! latest from CylinderCalculationsR1.ipynb

                  ! precomputed as mem%alpha(i) ... alpha0 = (r1*r1 + 2*r1*r2 + 3*r2**2)/4/(r1*r1 + r1*r2 + r2**2)
      !TODO: Review the below alpha eqn, GJH           
                  z1d = -min(0.0_ReKi,z1)
                  z2d = -min(0.0_ReKi,z2)
                   
                  pwr = 3
                  alpha = mem%alpha(i)*z2d**pwr/(mem%alpha(i)*z2d**pwr+(1-mem%alpha(i))*z1d**pwr)
                             
              
                  ! reduce moment to remove (not double count) radial force distribution to each node
                  Moment = Moment - Fr*alpha*dl
                  ! TODO: Should the order be, i, i+1 GJH
                  !call DistributeElementLoads(Fl, Fr, Moment, sinPhi, cosPhi, sinBeta, cosBeta, alpha, m%F_B(:, mem%NodeIndx(i+1)), m%F_B(:, mem%NodeIndx(i)))
                  call DistributeElementLoads(Fl, Fr, Moment, sinPhi, cosPhi, sinBeta, cosBeta, alpha, F_B1, F_B2)
                  m%memberLoads(im)%F_B(:,i+1) = m%memberLoads(im)%F_B(:,i+1) + F_B1  ! alpha
                  m%memberLoads(im)%F_B(:, i)  = m%memberLoads(im)%F_B(:, i)  + F_B2  ! 1-alpha
                  y%Mesh%Force (:,mem%NodeIndx(i  )) = y%Mesh%Force (:,mem%NodeIndx(i  )) + F_B2(1:3)
                  y%Mesh%Moment(:,mem%NodeIndx(i  )) = y%Mesh%Moment(:,mem%NodeIndx(i  )) + F_B2(4:6)
                  y%Mesh%Force (:,mem%NodeIndx(i+1)) = y%Mesh%Force (:,mem%NodeIndx(i+1)) + F_B1(1:3)
                  y%Mesh%Moment(:,mem%NodeIndx(i+1)) = y%Mesh%Moment(:,mem%NodeIndx(i+1)) + F_B1(4:6)
               end if  ! submergence cases
             
            end if ! element at least partially submerged
         
         end if ! NOT Modeled with Potential flow theory
      
         ! ------------------ flooded ballast inertia: sides: Section 6.1.1 : Always compute regardless of PropPot setting ---------------------

         ! lower node
         Ioffset   = mem%h_cfb_l(i)*mem%h_cfb_l(i)*mem%m_fb_l(i)
         Imat(1,1) = mem%I_rfb_l(i) - Ioffset
         Imat(2,2) = mem%I_rfb_l(i) - Ioffset
         Imat(3,3) = mem%I_lfb_l(i) - Ioffset
         iArm = mem%h_cfb_l(i) * k_hat
         iTerm     = ( -a_s1  - cross_product(omega_s1, cross_product(omega_s1,iArm ))  -  cross_product(alpha_s1,iArm) ) * mem%m_fb_l(i)
         !m%F_If(1:3, mem%NodeIndx(i  )) = m%F_If(1:3, mem%NodeIndx(i  )) + iTerm
         !m%F_If(4:6, mem%NodeIndx(i  )) = m%F_If(4:6, mem%NodeIndx(i  )) &
         !                                 - cross_product(a_s1 * mem%m_fb_l(i), mem%h_cfb_l(i) * k_hat) &
         !                                 + matmul(Imat, alpha_s1) &
         !                                 - cross_product(omega_s1,matmul(Imat,omega_s1)) 
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
         !m%F_If(1:3, mem%NodeIndx(i+1)) = m%F_If(1:3, mem%NodeIndx(i+1)) + iTerm
         !m%F_If(4:6, mem%NodeIndx(i+1)) = m%F_If(4:6, mem%NodeIndx(i+1)) &
         !                                 - cross_product(a_s2 * mem%m_fb_u(i), mem%h_cfb_u(i) * k_hat) &
         !                                 + matmul(Imat, alpha_s2) &
         !                                 - cross_product(omega_s2,matmul(Imat,omega_s2)) 
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
            !Fl = mem%Cfl_fb(i)*cosPhi     
            Fr = mem%Cfr_fb(i)*sinPhi     
            Moment  = mem%CM0_fb(i)*sinPhi - Fr*mem%alpha_fb_star(i)*dl
           
            ! calculate full vector and distribute to nodes
            !call DistributeElementLoads(Fl, Fr, Moment, sinPhi, cosPhi, sinBeta, cosBeta, (1-mem%alpha_fb_star(i)), m%F_BF(:, mem%NodeIndx(i)), m%F_BF(:, mem%NodeIndx(i+1)))
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
            !call DistributeElementLoads(Fl, Fr, Moment, sinPhi, cosPhi, sinBeta, cosBeta, mem%alpha_fb_star(i), m%F_BF(:, mem%NodeIndx(i)), m%F_BF(:, mem%NodeIndx(i-1)))
            call DistributeElementLoads(Fl, Fr, Moment, sinPhi, cosPhi, sinBeta, cosBeta, mem%alpha_fb_star(i), F_B1, F_B2)
            m%memberLoads(im)%F_BF(:, i) = m%memberLoads(im)%F_BF(:, i) + F_B1     ! alpha
            m%memberLoads(im)%F_BF(:, i-1) = m%memberLoads(im)%F_BF(:, i-1) + F_B2 ! 1- alpha
            y%Mesh%Force (:,mem%NodeIndx(i  )) = y%Mesh%Force (:,mem%NodeIndx(i  )) + F_B1(1:3)
            y%Mesh%Moment(:,mem%NodeIndx(i  )) = y%Mesh%Moment(:,mem%NodeIndx(i  )) + F_B1(4:6)
            y%Mesh%Force (:,mem%NodeIndx(i-1)) = y%Mesh%Force (:,mem%NodeIndx(i-1)) + F_B2(1:3)
            y%Mesh%Moment(:,mem%NodeIndx(i-1)) = y%Mesh%Moment(:,mem%NodeIndx(i-1)) + F_B2(4:6)    
        
         ! no load for unflooded element or element fully below seabed
        
         end if
        
           
         
   
      END DO ! i =1,N    ! loop through member elements       

       
      ! External Hydrodynamic Side Loads
         ! NOTE: All geometry-related calculations are based on the undisplaced configuration of the structure
      
      DO i =1,N+1    ! loop through member nodes
         ! We need to subtract the MSL2SWL offset to place this in the SWL reference system
         z1 = u%Mesh%Position(3, mem%NodeIndx(i)) - p%MSL2SWL
         if ( i > mem%i_floor .and. z1 <= 0.0 ) then  ! node is above (or at? TODO: check) seabed and below or at free-surface)
            ! TODO: Note that for computational efficiency, we could precompute h_c and deltal for each element when we are NOT using wave stretching
            ! We would still need to test at time marching for nodes just below the free surface because that uses the current locations not the reference locations
            ! see table in Section 7.1.1
            if ( i == 1 ) then
               deltal = mem%dl/2.0_ReKi
               h_c    = mem%dl/4.0_ReKi
            elseif (i == N+1) then
               deltal =  mem%dl/2.0_ReKi
               h_c    = -mem%dl/4.0_ReKi
            elseif ( mem%i_floor == i+1 ) then ! This node is the upper node of an element which crosses the seabed
               deltal = mem%dl/2.0_ReKi - mem%h_floor  ! TODO: h_floor is negative valued, should we be subrtracting it from dl/2? GJH
               h_c    = 0.5_ReKi*(mem%dl/2.0_ReKi + mem%h_floor)
            else
               ! We need to subtract the MSL2SWL offset to place this  in the SWL reference system
               pos1 =   u%Mesh%Position(:, mem%NodeIndx(i))
               pos1(3) = pos1(3) - p%MSL2SWL
               pos2 =   u%Mesh%Position(:, mem%NodeIndx(i+1))
               pos2(3) = pos2(3) - p%MSL2SWL
               if (pos1(3) <= 0.0 .and. 0.0 < pos2(3) ) then ! This node is just below the free surface !TODO: Needs to be augmented for wave stretching
                  ! We need to subtract the MSL2SWL offset to place this  in the SWL reference system
                  !TODO: Fix this one
                  pos1 =  u%Mesh%Position(:, mem%NodeIndx(i)) ! use reference position for following equation
                  pos1(3) = pos1(3) - p%MSL2SWL
                  h = (  pos1(3) ) / mem%cosPhi_ref !TODO: Needs to be augmented for wave stretching
                  deltal = mem%dl/2.0 + h
                  h_c    = 0.5*(h-mem%dl/2.0)
               else
                  ! This node is a fully submerged interior node
                  deltal = mem%dl
                  h_c    = 0.0_ReKi
               end if
            
            end if

            if (i == 1) then
               dRdl_p  = abs(mem%dRdl_mg(i))
               dRdl_pp = mem%dRdl_mg(i)   
            elseif ( i > 1 .and. i < (N+1)) then
               dRdl_p  = 0.5*( abs(mem%dRdl_mg(i-1)) + abs(mem%dRdl_mg(i)) )
               dRdl_pp = 0.5*( mem%dRdl_mg(i-1) + mem%dRdl_mg(i) )
            else 
               dRdl_p  = abs(mem%dRdl_mg(N))
               dRdl_pp = mem%dRdl_mg(N)
            end if
         
            ! ------------------- hydrodynamic drag loads: sides: Section 7.1.2 ------------------------ 
            vec = matmul( mem%Ak,m%vrel(:,mem%NodeIndx(i)) )
            f_hydro = mem%Cd(i)*p%WtrDens*mem%RMG(i)*TwoNorm(vec)*vec  +  &
                      0.5*mem%AxCd(i)*p%WtrDens*pi*mem%RMG(i)*dRdl_p * matmul( dot_product( mem%k, m%vrel(:,mem%NodeIndx(i)) )*mem%kkt, m%vrel(:,mem%NodeIndx(i)) )
!            call LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%F_D(:, mem%NodeIndx(i)) )
            call LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%memberLoads(im)%F_D(:, i) )
            y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + m%memberLoads(im)%F_D(1:3, i)
            y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + m%memberLoads(im)%F_D(4:6, i)
            
            if ( .not. mem%PropPot ) then
               ! ------------------- hydrodynamic added mass loads: sides: Section 7.1.3 ------------------------
               Am = mem%Ca(i)*p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*mem%Ak + 2.0*mem%AxCa(i)*p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*dRdl_p*mem%kkt
               f_hydro = -matmul( Am, u%Mesh%TranslationAcc(:,mem%NodeIndx(i)) )
               !call LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%F_A(:, mem%NodeIndx(i)) )
               call LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%memberLoads(im)%F_A(:, i) )
               y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + m%memberLoads(im)%F_A(1:3, i)
               y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + m%memberLoads(im)%F_A(4:6, i)
         
               ! ------------------- hydrodynamic inertia loads: sides: Section 7.1.4 ------------------------
               f_hydro=(mem%Ca(i)+mem%Cp(i))*p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)       * matmul( mem%Ak,  m%FA(:,mem%NodeIndx(i)) ) + &
                            2.0*mem%AxCa(i)*p%WtrDens*pi*mem%RMG(i)*mem%RMG(i)*dRdl_p * matmul( mem%kkt, m%FA(:,mem%NodeIndx(i)) ) + &
                            2.0*m%FDynP(mem%NodeIndx(i))*mem%AxCp(i)*pi*mem%RMG(i)*dRdl_pp*mem%k 
               !call LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%F_I(:, mem%NodeIndx(i)) )
               call LumpDistrHydroLoads( f_hydro, mem%k, deltal, h_c, m%memberLoads(im)%F_I(:, i) )
               y%Mesh%Force (:,mem%NodeIndx(i)) = y%Mesh%Force (:,mem%NodeIndx(i)) + m%memberLoads(im)%F_I(1:3, i)
               y%Mesh%Moment(:,mem%NodeIndx(i)) = y%Mesh%Moment(:,mem%NodeIndx(i)) + m%memberLoads(im)%F_I(4:6, i)
            end if
         end if ! ( i > mem%i_floor .and. Zi <= 0.0 )
         
      END DO ! i =1,N+1    ! loop through member nodes       
      
      
      ! Any end plate loads that are modeled on a per-member basis
      
      ! reassign convenience variables to correspond to member ends
      ! We need to subtract the MSL2SWL offset to place this  in the SWL reference system
      pos1    = u%Mesh%TranslationDisp(:, mem%NodeIndx(1)) + u%Mesh%Position(:, mem%NodeIndx(1)) 
      pos1(3) = pos1(3) - p%MSL2SWL
      pos2    = u%Mesh%TranslationDisp(:, mem%NodeIndx(2)) + u%Mesh%Position(:, mem%NodeIndx(2)) 
      pos2(3) = pos2(3) - p%MSL2SWL
      z1 = pos1(3)
      
      call GetOrientationAngles( pos1, pos2, phi1, sinPhi1, cosPhi1, tanPhi, sinBeta1, cosBeta1, k_hat1, errStat2, errMsg2 )
      if ( N == 1 ) then       ! Only one element in member
         sinPhi2 = sinPhi1
         cosPhi2 = cosPhi1
         sinBeta2  = sinBeta1
         cosBeta2  = cosBeta1
      else
         !  We need to subtract the MSL2SWL offset to place this  in the SWL reference system
         pos1    = u%Mesh%TranslationDisp(:, mem%NodeIndx(N))   + u%Mesh%Position(:, mem%NodeIndx(N))
         pos1(3) = pos1(3) - p%MSL2SWL
         pos2    = u%Mesh%TranslationDisp(:, mem%NodeIndx(N+1)) + u%Mesh%Position(:, mem%NodeIndx(N+1))
         pos2(3) = pos2(3) - p%MSL2SWL
         call GetOrientationAngles( pos1, pos2, phi2, sinPhi2, cosPhi2, tanPhi, sinBeta2, cosBeta2, k_hat2, errStat2, errMsg2 )
      end if
      ! We need to subtract the MSL2SWL offset to place this  in the SWL reference system
      pos2    = u%Mesh%TranslationDisp(:, mem%NodeIndx(N+1)) + u%Mesh%Position(:, mem%NodeIndx(N+1))
      pos2(3) = pos2(3) - p%MSL2SWL
      z2 = pos2(3)
      
      ! Check the member does not exhibit any of the following conditions
      if (.not. mem%PropPot) then 
         if ( abs(z2) < abs(mem%Rmg(N+1)*sinPhi2) ) then
            call SetErrStat(ErrID_Fatal, 'The upper end-plate of a member must not cross the water plane.  This is not true for Member ID '//trim(num2lstr(mem%MemberID)), errStat, errMsg, 'Morison_CalcOutput' )   
         end if
         if ( abs(z1) < abs(mem%Rmg(1)*sinPhi1) ) then
            call SetErrStat(ErrID_Fatal, 'The lower end-plate of a member must not cross the water plane.  This is not true for Member ID '//trim(num2lstr(mem%MemberID)), errStat, errMsg, 'Morison_CalcOutput' )   
         end if
      end if

! TODO: Do the equations below still work if z1 > z2 ?
 !TODO, should not have to test seabed crossing in time-marching loop

      
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

      ! --- no inertia loads from water ballast modeled on ends

      ! --- external buoyancy loads: ends ---

      if ( .not. mem%PropPot ) then
         ! We need to subtract the MSL2SWL offset to place this  in the SWL reference system
         pos1    = u%Mesh%TranslationDisp(:, mem%NodeIndx(1))   + u%Mesh%Position(:, mem%NodeIndx(1))
         pos1(3) = pos1(3) - p%MSL2SWL
         pos2    = u%Mesh%TranslationDisp(:, mem%NodeIndx(N+1)) + u%Mesh%Position(:, mem%NodeIndx(N+1))
         pos2(3) = pos2(3) - p%MSL2SWL
         z1 = pos1(3)
         z2 = pos2(3)
         if (mem%i_floor == 0) then  ! both ends above or at seabed
            if (z2<= 0.0_ReKi) then
               ! Compute loads on both ends
               Fl      = -p%WtrDens * g * pi *mem%RMG(1)**2*z1
               Moment  = -p%WtrDens * g * pi *0.25*mem%RMG(1)**4*sinPhi
               call AddEndLoad(Fl, Moment, sinPhi1, cosPhi1, sinBeta1, cosBeta1, m%F_B_End(:, mem%NodeIndx(1))) 
               Fl      = p%WtrDens * g * pi *mem%RMG(N+1)**2*z2
               Moment  = p%WtrDens * g * pi *0.25*mem%RMG(N+1)**4*sinPhi
               call AddEndLoad(Fl, Moment, sinPhi2, cosPhi2, sinBeta2, cosBeta2, m%F_B_End(:, mem%NodeIndx(N+1)))
            elseif ( z1< 0.0_ReKi ) then
               ! Compute loads only on lower end
               Fl      = -p%WtrDens * g * pi *mem%RMG(1)**2*z1
               Moment  = -p%WtrDens * g * pi *0.25*mem%RMG(1)**4*sinPhi
               call AddEndLoad(Fl, Moment, sinPhi1, cosPhi1, sinBeta1, cosBeta1, m%F_B_End(:, mem%NodeIndx(1)))
            else
               ! Entire member is above the still water line
            end if

    !     elseif ( (mem%i_floor < mem%NElements) .and. (z2<= 0.0_ReKi) ) then ! The member crosses the seabed line so only the upper end could have bouyancy effects, if at or below free surface
         elseif ( (mem%doEndBuoyancy) .and. (z2<= 0.0_ReKi) ) then ! The member crosses the seabed line so only the upper end could have bouyancy effects, if at or below free surface
            ! Only compute the buoyancy contribution from the upper end
            Fl      = p%WtrDens * g * pi *mem%RMG(N+1)**2*z2
            Moment  = p%WtrDens * g * pi *0.25*mem%RMG(N+1)**4*sinPhi
            call AddEndLoad(Fl, Moment, sinPhi2, cosPhi2, sinBeta2, cosBeta2, m%F_B_End(:, mem%NodeIndx(N+1)))
         else
            ! entire member is buried below the seabed         
         end if
         
      end if   ! PropPot
      
   end do ! im - looping through members
      
   !do j = 1, p%NNodes    
   !   ! Sum side load components onto output mesh
   !   DO i=1,6
   !      IF (i < 4 ) THEN 
   !         y%Mesh%Force(I,J)    =  m%F_D(I,J) + m%F_A(I,J) + m%F_I(I,J) + m%F_B(I,J) + m%F_BF(I,J) + m%F_If(i,j) + m%F_WMG(i,j) + m%F_IMG(i,j)
   !      ELSE 
   !         y%Mesh%Moment(I-3,J) =  m%F_D(I,J) + m%F_A(I,J) + m%F_I(I,J) + m%F_B(I,J) + m%F_BF(I,J) + m%F_If(i,j) + m%F_WMG(i,j) + m%F_IMG(i,j)   
   !      END IF
   !   END DO  ! 
   !end do
 

   ! --- Hydrodynamic drag loads: joints
      
      ! NOTE:  All wave kinematics have already been zeroed out above the SWL or instantaneous wave height (for WaveStMod > 0), so loads derived from the kinematics will be correct
      !        without the use of a nodeInWater value, but other loads need to be multiplied by nodeInWater to zero them out above the SWL or instantaneous wave height.
      
      DO J = 1, p%NJoints
         
            ! Obtain the node index because WaveVel, WaveAcc, and WaveDynP are defined in the node indexing scheme, not the markers

         
            ! Compute the dot product of the relative velocity vector with the directional Area of the Joint
         vmag =  m%nodeInWater(j) * ( m%vrel(1,j)*p%An_End(1,J) + m%vrel(2,j)*p%An_End(2,J) + m%vrel(3,j)*p%An_End(3,J) )
         
  !NOTE: The PropPot values are only for members, and when the p%AM_End, p%DP_Const_End, p%Mass_MG_End, and p%I_MG_End are computed at init,
  !      contributions to these values are added only if the member connecting to the joint is NOT modeled with potential flow theory
  !      However, the p%An_End term used data from ALL members attached to a node, regardless of the PropPot setting.
         
            ! Lumped added mass loads
         qdotdot                 = reshape((/u%Mesh%TranslationAcc(:,J),u%Mesh%RotationAcc(:,J)/),(/6/)) 
         m%F_A_End(:,J)          = m%nodeInWater(j) * matmul( p%AM_End(:,:,J) , ( - qdotdot(1:3)) )
         
         ! TODO: The original code did not multiply by nodeInWater, but should we? GJH
         m%F_I_End(:,J) =   (p%DP_Const_End(:,j) * m%FDynP(j) + matmul(p%AM_End(:,:,j),m%FA(:,j)))
         
         ! Marine growth inertia: ends: Section 4.2.2  
         m%F_IMG_End(1:3,j) = -m%nodeInWater(j) * p%Mass_MG_End(j)*qdotdot(1:3)
         m%F_IMG_End(4:6,j) = -m%nodeInWater(j) * (matmul(p%I_MG_End(:,:,j),qdotdot(4:6)) - cross_product(u%Mesh%RotationVel(:,J),matmul(p%I_MG_End(:,:,j),u%Mesh%RotationVel(:,J))))
         
         DO I=1,6
                        
            ! We are now combining the dynamic pressure term into the inertia term
            
            
            IF (I < 4 ) THEN
              
              
               m%F_D_End(i,j) =  p%An_End(i,j)*p%DragConst_End(j)*abs(vmag)*vmag  ! Note: vmag is zero if node is not in the water
               y%Mesh%Force(i,j)    = y%Mesh%Force(i,j)    + m%F_D_End(i,j) + m%F_I_End(i,j) + p%F_WMG_End(i,j) + m%F_B_End(i,j) + m%F_BF_End(i,j) + m%F_A_End(i,j) + m%F_IMG_End(i,j)
            ELSE
               y%Mesh%Moment(i-3,j) = y%Mesh%Moment(i-3,j) + m%F_B_End(i,j) + m%F_BF_End(i,j)  + m%F_IMG_End(i,j)
            END IF
         END DO      ! I=1,6
      ENDDO          ! J = 1, p%NJoints
     
         ! OutSwtch determines whether or not to actually output results via the WriteOutput array
         ! 1 = Morison will generate an output file of its own.  2 = the caller will handle the outputs, but
         ! Morison needs to provide them.  3 = Both 1 and 2, 0 = No one needs the Morison outputs provided
         ! via the WriteOutput array.
         
      IF ( p%OutSwtch > 0 ) THEN
     
            ! Map calculated results into the AllOuts Array
         CALL MrsnOut_MapOutputs(Time, y, p, u, m, AllOuts, errStat, errMsg)
               
      
            ! Put the output data in the WriteOutput array
   
         DO I = 1,p%NumOuts

            y%WriteOutput(I) = p%OutParam(I)%SignM * AllOuts( p%OutParam(I)%Indx )
      
         END DO
         
         
            ! Generate output into the output file
            
         IF ( p%OutSwtch == 1 .OR. p%OutSwtch == 3 ) THEN
            CALL MrsnOut_WriteOutputs( p%UnOutFile, Time, y, p, errStat, errMsg )         
         END IF
      END IF
      
   
END SUBROUTINE Morison_CalcOutput

subroutine LumpDistrHydroLoads( f_hydro, k_hat, dl, h_c, lumpedLoad )
   real(ReKi), intent(in   ) :: f_hydro(3)
   real(ReKi), intent(in   ) :: k_hat(3)
   real(ReKi), intent(in   ) :: dl
   real(ReKi), intent(in   ) :: h_c
   real(ReKi), intent(inout) :: lumpedLoad(6)
   !lumpedLoad(1:3) = lumpedLoad(1:3) + f_hydro*dl
   !lumpedLoad(4:6) = lumpedLoad(4:6) + cross_product(k_hat*h_c, f_hydro)*dl
   lumpedLoad(1:3) = f_hydro*dl
   lumpedLoad(4:6) = cross_product(k_hat*h_c, f_hydro)*dl
end subroutine LumpDistrHydroLoads

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
         

   !F1(1) = F1(1) +  cosBeta*(Fl*sinPhi + Fr*cosPhi)*alpha
   !F1(2) = F1(2) -  sinBeta*(Fl*sinPhi + Fr*cosPhi)*alpha
   !F1(3) = F1(3) +          (Fl*cosPhi - Fr*sinPhi)*alpha
   !F1(4) = F1(4) +  sinBeta * M                    *alpha
   !F1(5) = F1(5) +  cosBeta * M                    *alpha
   !!F1(6) = F1(6) + 0.0
   !   
   !F2(1) = F2(1) +  cosBeta*(Fl*sinPhi + Fr*cosPhi)*(1-alpha)
   !F2(2) = F2(2) -  sinBeta*(Fl*sinPhi + Fr*cosPhi)*(1-alpha)
   !F2(3) = F2(3) +          (Fl*cosPhi - Fr*sinPhi)*(1-alpha)
   !F2(4) = F2(4) +  sinBeta * M                    *(1-alpha)
   !F2(5) = F2(5) +  cosBeta * M                    *(1-alpha)
   !!F2(6) = F2(6) + 0.0
   
   F1(1) =  cosBeta*(Fl*sinPhi + Fr*cosPhi)*alpha
   F1(2) =  sinBeta*(Fl*sinPhi + Fr*cosPhi)*alpha
   F1(3) =         (Fl*cosPhi - Fr*sinPhi)*alpha
   F1(4) =  -sinBeta * M                    *alpha
   F1(5) =  cosBeta * M                    *alpha
   F1(6) =  0.0
      
   F2(1) =  cosBeta*(Fl*sinPhi + Fr*cosPhi)*(1-alpha)
   F2(2) =  sinBeta*(Fl*sinPhi + Fr*cosPhi)*(1-alpha)
   F2(3) =          (Fl*cosPhi - Fr*sinPhi)*(1-alpha)
   F2(4) =  -sinBeta * M                    *(1-alpha)
   F2(5) =  cosBeta * M                    *(1-alpha)
   F2(6) = 0.0
   
   !F1(1) =  cosBeta*(-Fl*sinPhi + Fr*cosPhi)*alpha
   !F1(2) =  sinBeta*(-Fl*sinPhi + Fr*cosPhi)*alpha
   !F1(3) =         (Fl*cosPhi + Fr*sinPhi)*alpha
   !F1(4) =  -sinBeta * M                    *alpha
   !F1(5) =  cosBeta * M                    *alpha
   !F1(6) =  0.0
   !   
   !F2(1) =  cosBeta*(-Fl*sinPhi + Fr*cosPhi)*(1-alpha)
   !F2(2) =  sinBeta*(-Fl*sinPhi + Fr*cosPhi)*(1-alpha)
   !F2(3) =          (Fl*cosPhi + Fr*sinPhi)*(1-alpha)
   !F2(4) =  -sinBeta * M                    *(1-alpha)
   !F2(5) =  cosBeta * M                    *(1-alpha)
   !F2(6) = 0.0

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
   Fi(4) = Fi(4) - M*sinBeta
   Fi(5) = Fi(5) + M*cosBeta
   
END SUBROUTINE AddEndLoad


!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
SUBROUTINE Morison_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, errStat, errMsg )  
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
      INTEGER(IntKi),                    INTENT(  OUT)  :: errStat     !< Error status of the operation     
      CHARACTER(*),                      INTENT(  OUT)  :: errMsg      !< Error message if errStat /= ErrID_None

               
         ! Initialize errStat
         
      errStat = ErrID_None         
      errMsg  = ""               
      
      
         ! Compute the first time derivatives of the continuous states here:
      
      dxdt%DummyContState = 0.0
         

END SUBROUTINE Morison_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for updating discrete states
SUBROUTINE Morison_UpdateDiscState( Time, u, p, x, xd, z, OtherState, m, errStat, errMsg )   
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
      INTEGER(IntKi),                    INTENT(  OUT)  :: errStat     !< Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: errMsg      !< Error message if errStat /= ErrID_None

               
         ! Initialize errStat
         
      errStat = ErrID_None         
      errMsg  = ""               
      
      
         ! Update discrete states here:
      
      ! StateData%DiscState = 

END SUBROUTINE Morison_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
SUBROUTINE Morison_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, errStat, errMsg )   
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
      INTEGER(IntKi),                    INTENT(  OUT)  :: errStat     !< Error status of the operation
      CHARACTER(*),                      INTENT(  OUT)  :: errMsg      !< Error message if errStat /= ErrID_None

               
         ! Initialize errStat
         
      errStat = ErrID_None         
      errMsg  = ""               
      
      
         ! Solve for the constraint states here:
      
      z_residual%DummyConstrState = 0.0_ReKi

END SUBROUTINE Morison_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
   
END MODULE Morison
!**********************************************************************************************************************************
