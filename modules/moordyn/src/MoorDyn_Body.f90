!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2020-2021 Alliance for Sustainable Energy, LLC
! Copyright (C) 2015-2019 Matthew Hall
!
!    This file is part of MoorDyn.
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
MODULE MoorDyn_Body

   USE MoorDyn_Types
   USE MoorDyn_IO
   USE NWTC_Library
   USE MoorDyn_Misc
   !USE MoorDyn_Line,           only : Line_SetEndKinematics, Line_GetEndStuff
   USE MoorDyn_Point,           only : Point_SetKinematics, Point_GetNetForceAndMass
   USE MoorDyn_Rod,             only : Rod_Initialize, Rod_SetKinematics, Rod_GetNetForceAndMass
   
   IMPLICIT NONE

   PRIVATE

   INTEGER(IntKi), PARAMETER            :: wordy = 0   ! verbosity level. >1 = more console output

   PUBLIC :: Body_Setup
   PUBLIC :: Body_Initialize
   PUBLIC :: Body_InitializeUnfree
   PUBLIC :: Body_SetKinematics
   PUBLIC :: Body_SetState
   PUBLIC :: Body_SetDependentKin
   PUBLIC :: Body_GetStateDeriv
   PUBLIC :: Body_DoRHS
   PUBLIC :: Body_GetCoupledForce
   PUBLIC :: Body_AddPoint
   PUBLIC :: Body_AddRod
   
   

CONTAINS


   SUBROUTINE Body_Setup( Body, tempArray, p, ErrStat, ErrMsg)

      TYPE(MD_Body),     INTENT(INOUT)    :: Body          ! the single body object of interest
      REAL(DbKi),        INTENT(IN)       :: tempArray(6)  ! initial pose of body
      TYPE(MD_ParameterType), INTENT(IN   ) :: p       ! Parameters
      INTEGER,           INTENT(INOUT )   :: ErrStat       ! returns a non-zero value when an error occurs
      CHARACTER(*),      INTENT(INOUT )   :: ErrMsg        ! Error message if ErrStat /= ErrID_None

      INTEGER(4)                          :: J             ! Generic index
!      INTEGER(4)                          :: K             ! Generic index
!      INTEGER(IntKi)                      :: N

      REAL(DbKi)                          :: Mtemp(6,6) = 0.0_DbKi  ! temporary mass matrix
      
      ErrStat = ErrID_None
      ErrMsg = ""

      ! set initial velocity to zero
      Body%v6 = 0.0_DbKi

      ! set external load to zero
      Body%FextG = 0.0_DbKi
      Body%FextL = 0.0_DbKi

      ! set external damping to zero
      Body%BlinG  = 0.0_DbKi
      Body%BquadG = 0.0_DbKi
      Body%BlinL  = 0.0_DbKi
      Body%BquadL = 0.0_DbKi

      !also set number of attached rods and points to zero initially
      Body%nAttachedC = 0
      Body%nAttachedR = 0

      ! set up body initial mass matrix (excluding any rods or attachements)
      DO J=1,3
         Mtemp(J,J) = Body%BodyM          ! fill in mass
         Mtemp(3+J,3+J) = Body%BodyI(J)   ! fill in inertia   
      END DO
      
      CALL TranslateMass6to6DOF(Body%rCG, Mtemp, Body%M0)  ! account for potential CG offset <<< is the direction right? <<<
        
      DO J=1,3
         Body%M0(J,J) = Body%M0(J,J) + Body%BodyV*Body%BodyCa(J)* p%rhow ! add added mass in each direction about ref point (so only diagonals) <<< eventually expand to multi D
      END DO
   
      ! --------------- if this is an independent body (not coupled) ----------
      ! set initial position and orientation of body from input file 
      Body%r6(1:3) = tempArray(1:3)
      Body%r6(4:6) = tempArray(4:6) * (pi/180)

      ! calculate orientation matrix based on latest angles
      !RotMat(r6[3], r6[4], r6[5], OrMat);
      Body%OrMat = TRANSPOSE( EulerConstruct( Body%r6(4:6) ) )  ! full Euler angle approach <<<< need to check order 

      IF (wordy > 0) print *, "Set up Body ",Body%IdNum, ", type ", Body%typeNum

      ! need to add cleanup sub <<<

   END SUBROUTINE Body_Setup

!   ! used to initialize bodies that aren't free i.e. don't have states
!   !--------------------------------------------------------------
!   SUBROUTINE Body_InitializeUnfree(Body, r6_in, mesh, mesh_index, m)
!
!      Type(MD_Body),         INTENT(INOUT)  :: Body        ! the Body object
!      Real(DbKi),            INTENT(IN   )  :: r6_in(6)    ! state vector section for this line
!      TYPE(MeshType),        INTENT(INOUT)  :: mesh        !
!      Integer(IntKi),        INTENT(IN   )  :: mesh_index  ! index of the node in the mesh for the current object being initialized
!      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m           ! passing along all mooring objects
!
!      INTEGER(IntKi)                        :: l           ! index of segments or nodes along line
!      REAL(DbKi)                            :: rRef(3)     ! reference position of mesh node
!      REAL(DbKi)                            :: OrMat(3,3)  ! DCM for body orientation based on r6_in
!      REAL(DbKi)                            :: dummyStates(12) 
!   
!   
!      rRef = 0.0_DbKi   ! <<< maybe this should be the offsets of the local platform origins from the global origins in future? And that's what's specificed by the Body input coordinates?
!      
!      CALL MeshPositionNode(mesh, mesh+index, rRef,ErrStat2,ErrMsg2)! "assign the coordinates (u%PtFairleadDisplacement%Position) of each node in the global coordinate space"
!
!      CALL CheckError( ErrStat2, ErrMsg2 )
!      IF (ErrStat >= AbortErrLev) RETURN
!
!      ! Apply offsets due to initial platform rotations and translations (fixed Jun 19, 2015)
!      CALL SmllRotTrans('body rotation matrix', r6_in(4),r6_in(5),r6_in(6), OrMat, '', ErrStat2, ErrMsg2)
!      mesh%TranslationDisp(1, mesh_index) = r6_in(1) + OrMat(1,1)*rRef(1) + OrMat(2,1)*rRef(2) + OrMat(3,1)*rRef(3) - rRef(1)
!      mesh%TranslationDisp(2, mesh_index) = r6_in(2) + OrMat(1,2)*rRef(1) + OrMat(2,2)*rRef(2) + OrMat(3,2)*rRef(3) - rRef(2)
!      mesh%TranslationDisp(3, mesh_index) = r6_in(3) + OrMat(1,3)*rRef(1) + OrMat(2,3)*rRef(2) + OrMat(3,3)*rRef(3) - rRef(3)
!
!      ! what about node point orientation ???
!
!      ! If any Rod is fixed to the body (not pinned), initialize it now because otherwise it won't be initialized
!      DO l=1, Body%nAttachedR
!         if (m%RodList(Body%attachedR(l))%typeNum == 2)  CALL Rod_Initialize(m%RodList(Body%attachedR(l)), dummyStates, m%LineList)
!      END DO
!      
!      ! Note: Points don't need any initialization
!      
!   END SUBROUTINE Body_InitializeUnfree
!   !--------------------------------------------------------------


   ! used to initialize bodies that are free
   !--------------------------------------------------------------
   SUBROUTINE Body_Initialize(Body, states, m)

      Type(MD_Body),         INTENT(INOUT)  :: Body            ! the Body object
      Real(DbKi),            INTENT(INOUT)  :: states(:)       ! state vector section for this Body
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m               ! passing along all mooring objects

      INTEGER(IntKi)                        :: l               ! index of segments or nodes along line
      REAL(DbKi)                            :: dummyStates(12) ! dummy vector to mimic states when initializing a rigidly attached rod
   
      IF (wordy > 0) print *, "initializing Body ", Body%idNum

      ! the r6 and v6 vectors should have already been set
      ! r and rd of ends have already been set by setup function or by parent object   <<<<< right? <<<<<


      if (Body%typeNum == 0) then               ! free body type
      
         ! assign initial body kinematics to state vector
         states(1:6 ) = Body%v6 ! zero velocities for initialization (set to 0 in Body_Setup)
         states(7:12) = Body%r6
      
      else if (Body%typeNum ==2 ) then           ! pinned rod type (coupled or attached to something previously via setPinKin)
      
         states(1:3)   = Body%v6(4:6) ! zero velocities for initialization (set to 0 in Body_Setup)
         states(4:6)   = Body%r6(4:6) ! body orentations
         
      end if
      
      ! set positions of any dependent points and rods now (before they are initialized)
      CALL Body_SetDependentKin(Body, 0.0_DbKi, m)
            
      ! If any Rod is fixed to the body (not pinned), initialize it now because otherwise it won't be initialized
      DO l=1, Body%nAttachedR
         if (m%RodList(Body%attachedR(l))%typeNum == 2)  CALL Rod_Initialize(m%RodList(Body%attachedR(l)), dummyStates,  m)
      END DO
      
      ! Note: Points don't need any initialization
      
   END SUBROUTINE Body_Initialize
   !--------------------------------------------------------------
   
   ! used to initialize bodies that are coupled or fixed
   !--------------------------------------------------------------
   SUBROUTINE Body_InitializeUnfree(Body, m)

      Type(MD_Body),         INTENT(INOUT)  :: Body            ! the Body object
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m               ! passing along all mooring objects

      INTEGER(IntKi)                        :: l               ! index of segments or nodes along line
      REAL(DbKi)                            :: dummyStates(12) ! dummy vector to mimic states when initializing a rigidly attached rod
   
   
      ! set positions of any dependent points and rods now (before they are initialized)
      CALL Body_SetDependentKin(Body, 0.0_DbKi, m)
            
      ! If any Rod is fixed to the body (not pinned), initialize it now because otherwise it won't be initialized
      DO l=1, Body%nAttachedR
         if (m%RodList(Body%attachedR(l))%typeNum == 2)  CALL Rod_Initialize(m%RodList(Body%attachedR(l)), dummyStates,  m)
      END DO
      
      ! Note: Points don't need any initialization
      
   END SUBROUTINE Body_InitializeUnfree
   !--------------------------------------------------------------




   ! set kinematics for Bodies if they are coupled (or ground)
   !--------------------------------------------------------------
   SUBROUTINE Body_SetKinematics(Body, r6_in, v6_in, a6_in, t, m)

      Type(MD_Body),         INTENT(INOUT)  :: Body       ! the Body object
      Real(DbKi),            INTENT(IN   )  :: r6_in(6)   ! 6-DOF position
      Real(DbKi),            INTENT(IN   )  :: v6_in(6)   ! 6-DOF velocity
      Real(DbKi),             INTENT(IN   ) :: a6_in(6)       ! 6-DOF acceleration 
      Real(DbKi),            INTENT(IN   )  :: t         ! instantaneous time
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m         ! passing along all mooring objects (for simplicity, since Bodies deal with Rods and Points)


!      INTEGER(IntKi)                   :: l

      ! store current time
      Body%time = t

      if (Body%typeNum == 2) then ! body pinned to coupling point
      
         ! set Body translational kinematics based on BCs (linear model for now) 
         Body%r6(1:3) = r6_in(1:3)
         Body%v6(1:3) = v6_in(1:3)
         Body%a6(1:3) = a6_in(1:3)

         ! Body rotations are left alone and will be handled, along with passing kinematics to dependent objects, by separate call to setState

      else ! body rigidly coupled to coupling point
         Body%r6 = r6_in
         Body%v6 = v6_in
         Body%a6 = a6_in
                  
         ! since this body has no states and all DOFs have been set, pass its kinematics to dependent attachments
         CALL Body_SetDependentKin(Body, t, m)

      end if

   END SUBROUTINE Body_SetKinematics
   !--------------------------------------------------------------


   !--------------------------------------------------------------
   SUBROUTINE Body_SetState(Body, X, t, m)

      Type(MD_Body),         INTENT(INOUT)  :: Body           ! the Body object
      Real(DbKi),            INTENT(IN   )  :: X(:)           ! state vector section for this line
      Real(DbKi),            INTENT(IN   )  :: t              ! instantaneous time
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m              ! passing along all mooring objects

!      INTEGER(IntKi)                        :: l              ! index of segments or nodes along line
!      INTEGER(IntKi)                        :: J              ! index
   
      ! store current time
      Body%time = t
      
      if (Body%typeNum == 0) then ! free Body type
      
         Body%r6 = X(7:12)   ! get positions      
         Body%v6 = X(1:6)    ! get velocities

         ! set positions of any dependent points and rods
         CALL Body_SetDependentKin(Body, t, m) 

      else if (Body%typeNum == 2) then

         Body%r6(4:6) = X(4:6) ! get positions
         Body%v6(4:6) = X(1:3) ! get velocities


         ! set positions of any dependent points and rods
         CALL Body_SetDependentKin(Body, t, m) 

      else 
         Call WrScr("Error: Body::setState called for a non-free Body type in MoorDyn")   ! <<<
      end if
      
   END SUBROUTINE Body_SetState
   !--------------------------------------------------------------


   ! set the states (positions and velocities) of any points or rods that are part of this body
   ! also computes the orientation matrix (never skip this sub!)
   !--------------------------------------------------------------
   SUBROUTINE Body_SetDependentKin(Body, t, m)

      Type(MD_Body),         INTENT(INOUT)  :: Body        ! the Bodyion object
      REAL(DbKi),            INTENT(IN   )  :: t
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m           ! passing along all mooring objects (for simplicity, since Bodies deal with Rods and Points)

      INTEGER(IntKi)                        :: l              ! index of attached objects
   
      Real(DbKi)                            :: rPoint(3)
      Real(DbKi)                            :: rdPoint(3)
      Real(DbKi)                            :: rRod(6)
      Real(DbKi)                            :: vRod(6)
      Real(DbKi)                            :: aRod(6)

      

      ! calculate orientation matrix based on latest angles
      !CALL SmllRotTrans('', Body%r6(4), Body%r6(5), Body%r6(6), Body%TransMat, '', ErrStat2, ErrMsg2)
      Body%OrMat = TRANSPOSE( EulerConstruct( Body%r6(4:6) ) ) ! full Euler angle approach <<<< need to check order 
  
      ! set kinematics of any dependent points
      do l = 1,Body%nAttachedC
      
         CALL transformKinematics(Body%rPointRel(:,l), Body%r6, Body%OrMat, Body%v6, rPoint, rdPoint) !<<< should double check this function
                  
         ! >>> need to add acceleration terms here too? <<<
                  
         ! pass above to the point and get it to calculate the forces
         CALL Point_SetKinematics( m%PointList(Body%attachedC(l)), rPoint, rdPoint, m%zeros6(1:3), t, m)
      end do
      
      ! set kinematics of any dependent Rods
      do l=1,Body%nAttachedR
      
         ! calculate displaced coordinates/orientation and velocities of each rod <<<<<<<<<<<<<
         ! do 3d details of Rod ref point
         CALL TransformKinematicsA( Body%r6RodRel(1:3,l), Body%r6(1:3), Body%OrMat, Body%v6, Body%a6, rRod(1:3), vRod(1:3), aRod(1:3))  ! set first three entires (end A translation) of rRod and rdRod
         ! does the above function need to take in all 6 elements of r6RodRel??
         
         ! do rotational stuff
         rRod(4:6) = MATMUL(Body%OrMat, Body%r6RodRel(4:6,l))    !<<<<<< correct? <<<<< rotateVector3(r6RodRel[i]+3, OrMat, rRod+3);   ! rotate rod relative unit vector by OrMat to get unit vec in reference coords
         vRod(4:6) = Body%v6(4:6)  ! transformed rotational velocity.  <<< is this okay as is? <<<<
         aRod(4:6) = Body%a6(4:6) 
         
         ! pass above to the rod and get it to calculate the forces
         CALL Rod_SetKinematics(m%RodList(Body%attachedR(l)), rRod, vRod, aRod, t, m)
      end do

   END SUBROUTINE Body_SetDependentKin
   !--------------------------------------------------------------
   

   !--------------------------------------------------------------
   SUBROUTINE Body_GetStateDeriv(Body, Xd, m, p)

      Type(MD_Body),         INTENT(INOUT)  :: Body          ! the Bodyion object
      Real(DbKi),            INTENT(INOUT)  :: Xd(:)            ! state derivative vector section for this line
      
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m           ! passing along all mooring objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p                ! Parameters
      
      INTEGER(IntKi)                        :: J                ! index
      
      Real(DbKi)                            :: Fnet     (6)     ! net force and moment about reference point

      Real(DbKi)                            :: acc(6)           ! 6DOF acceleration vector
      
      Real(DbKi)                            :: y_temp (6)       ! temporary vector for LU decomposition
      Real(DbKi)                            :: LU_temp(6,6)     ! temporary matrix for LU decomposition
      

      ! Initialize temp variables
      y_temp   = 0.0_DbKi
! FIXME: should LU_temp be set to M_out before calling LUsolve?????
      LU_temp  = 0.0_DbKi

      CALL Body_DoRHS(Body, m, p)

      IF (Body%typeNum == 0) THEN ! Free body

         ! solve for accelerations in [M]{a}={f} using LU decomposition
         CALL LUsolve(6, Body%M, LU_temp, Body%F6net, y_temp, acc)

         ! fill in state derivatives
         Xd(7:12) = Body%v6       ! dxdt = V   (velocities)
         Xd(1:6)  = acc           ! dVdt = a   (accelerations) 

         ! store accelerations in case they're useful as output
         Body%a6 = acc

      ELSE ! Pinned Body, 3 states (rotational only)

         ! Account for moment response due to inertial coupling
         Fnet = Body%F6net
         Fnet(4:6) = Fnet(4:6) - MATMUL(Body%M(4:6,1:3), Body%a6(1:3))  

         ! solve for accelerations in [M]{a}={f} using LU decomposition
         CALL LUsolve(3, Body%M(4:6,4:6), LU_temp(4:6,4:6), Fnet(4:6), y_temp(4:6), acc(4:6))

         ! fill in state derivatives
         Xd(4:6) = Body%v6(4:6)       ! dxdt = V   (velocities)
         Xd(1:3)  = acc(4:6)           ! dVdt = a   (accelerations) 

         ! store accelerations in case they're useful as output
         Body%a6(4:6) = acc(4:6)

      ENDIF
   
      ! check for NaNs (should check all state derivatives, not just first 6)
      DO J = 1, 6
         IF (Is_NaN(Xd(J))) THEN
            CALL WrScr("NaN detected at time "//trim(Num2LStr(Body%time))//" in Body "//trim(Int2LStr(Body%IdNum))//" in MoorDyn,")
            IF (wordy > 0) print *, "state derivatives:"
            IF (wordy > 0) print *, Xd
            EXIT
         END IF
      END DO


   END SUBROUTINE Body_GetStateDeriv
   !--------------------------------------------------------------

   !--------------------------------------------------------------
   SUBROUTINE Body_DoRHS(Body, m, p)

      Type(MD_Body),         INTENT(INOUT)  :: Body        ! the Bodyion object
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m           ! passing along all mooring objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p           ! Parameters
      
      !TYPE(MD_MiscVarType), INTENT(INOUT)  :: m       ! misc/optimization variables

      INTEGER(IntKi)             :: l         ! index of attached lines
      INTEGER(IntKi)             :: i         ! Generic loop counter

      Real(DbKi)                 :: Fgrav(3)           ! body weight force
      Real(DbKi)                 :: body_rCGrotated(3) ! instantaneous vector from body ref point to CG
      Real(DbKi)                 :: U(3)               ! water velocity - zero for now
!      Real(DbKi)                 :: Ud(3)              ! water acceleration- zero for now
      Real(DbKi)                 :: vi(6)              ! relative water velocity (last 3 terms are rotatonal and will be set to zero
      Real(DbKi)                 :: F6_i(6)            ! net force and moments from an attached object
      Real(DbKi)                 :: M6_i(6,6)          ! mass and inertia from an attached object
      Real(DbKi)                 :: cda(6)             ! body drag coefficients
      Real(DbKi)                 :: cda_t(3,3) = 0.0         ! matrix with translational drag coefficients as diagonals
      Real(DbKi)                 :: cda_r(3,3) = 0.0         ! matrix with rotational drag coefficients as diagonals
      Real(DbKi)                 :: w(3)                     ! body angular velocity vector
      Real(DbKi)                 :: Fcentripetal(3)        ! centripetal force
      Real(DbKi)                 :: Mcentripetal(3)        ! centripetal moment     
      Real(DbKi)                 :: v3L(3)             ! Body translational velocity in the local body-fixed coordinate system
      Real(DbKi)                 :: FDL(3)             ! Part of user-defined damping force defined in the local body-fixed coordinate system


      ! Initialize variables
      U = 0.0_DbKi      ! Set to zero for now
      Body%F6net = 0.0_DbKi

      ! First, the body's own mass matrix must be adjusted based on its orientation so that 
      ! we have a mass matrix in the global orientation frame
      Body%M = RotateM6(Body%M0, Body%OrMat)

      !gravity on core body
      Fgrav(1) = 0.0_DbKi
      Fgrav(2) = 0.0_DbKi
      Fgrav(3) = Body%bodyV * p%rhow * p%g - Body%bodyM * p%g ! weight+buoyancy vector

      body_rCGrotated = MATMUL(Body%OrMat, Body%rCG) ! rotateVector3(body_rCG, OrMat, body_rCGrotated); ! relative vector to body CG in inertial orientation
      CALL translateForce3to6DOF(body_rCGrotated, Fgrav, Body%F6net)  ! gravity forces and moments about body ref point given CG location

      ! Add user-defined external force and damping on body defined in the global earth-fixed coordinate system (assumed to be applied at the body ref point)
      Body%F6net(1:3) = Body%F6net(1:3) + Body%FextG
      do i = 1,3
         Body%F6net(i) = Body%F6net(i) - Body%BlinG(i) * Body%v6(i) - Body%BquadG(i) * ABS(Body%v6(i)) * Body%v6(i)
      end do

      ! Add user-defined external force and damping on body defined in the local body-fixed coordinate system (assumed to be applied at the body ref point)
      Body%F6net(1:3) = Body%F6net(1:3) + MATMUL( Body%OrMat, Body%FextL)
      v3L = MATMUL( TRANSPOSE(Body%OrMat), Body%v6(1:3) )
      do i = 1,3
         FDL(i) = - Body%BlinL(i) * v3L(i) - Body%BquadL(i) * ABS(v3L(i)) * v3L(i)
      end do
      Body%F6net(1:3) = Body%F6net(1:3) + MATMUL( Body%OrMat, FDL )

      ! Centripetal force and moment due to COM not being at body origin plus gyroscopic moment
      w = Body%v6(4:6)
      Fcentripetal = - MATMUL(Body%M(1:3,1:3), CROSS_PRODUCT(w, CROSS_PRODUCT(w, body_rCGrotated)))
      Mcentripetal = - CROSS_PRODUCT(w, MATMUL(Body%M(4:6,4:6), w)) 

      Body%F6net(1:3) = Body%F6net(1:3) + Fcentripetal
      Body%F6net(4:6) = Body%F6net(4:6) + Mcentripetal

      ! --------------------------------- apply wave kinematics ------------------------------------
      !env->waves->getU(r6, t, U); ! call generic function to get water velocities <<<<<<<<< all needs updating

      !   for (int J=0; J<3; J++)
      !      Ud[J] = 0.0;                 ! set water accelerations as zero for now
      ! ------------------------------------------------------------------------------------------

      ! viscous drag calculation (on core body)
      vi(1:3) = U - Body%v6(1:3)  ! relative flow velocity over body ref point
      vi(4:6) =   - Body%v6(4:6)  ! for rotation, this is just the negative of the body's rotation for now (not allowing flow rotation)

      cda_t(1,1) = Body%bodyCdA(1)
      cda_t(2,2) = Body%bodyCdA(2)
      cda_t(3,3) = Body%bodyCdA(3)
      cda_r(1,1) = Body%bodyCdA(4)
      cda_r(2,2) = Body%bodyCdA(5)
      cda_r(3,3) = Body%bodyCdA(6)

      cda(1:3) = MATMUL( MATMUL( MATMUL(Body%OrMat,cda_t) , transpose(Body%OrMat) ) , vi(1:3) * norm2(vi(1:3)) );
      cda(4:6) = MATMUL( MATMUL( MATMUL(Body%OrMat,cda_r) , transpose(Body%OrMat) ) , vi(4:6) * norm2(vi(4:6)) );
      Body%F6net = Body%F6net + 0.5*p%rhoW*cda

      

   
   
      ! Get contributions from any dependent points
      do l = 1,Body%nAttachedC
      
         ! get net force and mass from Point on body ref point (global orientation)
         CALL Point_GetNetForceAndMass( m%PointList(Body%attachedC(l)), Body%r6(1:3), Body%v6(4:6), F6_i, M6_i, m, p)
         
         if (ABS(F6_i(5)) > 1.0E12) then
            Call WrScr( "Warning: extreme pitch moment from body-attached Point "//trim(num2lstr(l)))
         end if
         
         ! sum quantitites
         Body%F6net = Body%F6net + F6_i
         Body%M     = Body%M     + M6_i
       end do
      
      ! Get contributions from any dependent Rods
      do l=1,Body%nAttachedR
      
         ! get net force and mass from Rod on body ref point (global orientation)
         CALL Rod_GetNetForceAndMass(m%RodList(Body%attachedR(l)), Body%r6(1:3), Body%v6(4:6), F6_i, M6_i, m, p)
         
         if (ABS(F6_i(5)) > 1.0E12) then
            Call WrScr("Warning: extreme pitch moment from body-attached Rod "//trim(num2lstr(l)))
         end if
         
         ! sum quantitites
         Body%F6net = Body%F6net + F6_i
         Body%M     = Body%M     + M6_i
      end do


   END SUBROUTINE Body_DoRHS
   !=====================================================================


      ! calculate the aggregate 3/6DOF rigid-body loads of a coupled rod including inertial loads
   !--------------------------------------------------------------
   SUBROUTINE Body_GetCoupledForce(t, Body, Fnet_out, m, p)
      real(R8Ki),            intent(in   )  :: t           ! time - for ramping inertial loading
      Type(MD_Body),         INTENT(INOUT)  :: Body        ! the Body object
      Real(DbKi),            INTENT(  OUT)  :: Fnet_out(6) ! force and moment vector
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m           ! passing along all mooring objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p           ! Parameters
      
      Real(DbKi)                            :: F6_iner(6)  ! inertial reaction force
      
      ! do calculations of forces and masses on the body
      CALL Body_DoRHS(Body, m, p)

      ! add inertial loads as appropriate 
      if (Body%typeNum == -1) then                          
      
         if (p%inertialF == 1) then      ! include inertial components 
            F6_iner = -MATMUL(Body%M, Body%a6)     ! unstable in OpenFAST v4 and below becasue of loose coupling with ED and SD. Transients in acceleration can cause issues
         elseif (p%inertialF == 2) then  ! include inertial components, but ramp up load
            F6_iner = -MATMUL(Body%M, Body%a6)
            if (t < p%inertialF_rampT) F6_iner = F6_iner * t / p%inertialF_rampT
         else
            ! When OpenFAST v5 is released w/ tight coupling, remove this hack and just use the inertial term above 
            F6_iner = 0.0
         endif

         Body%F6net = Body%F6net + F6_iner        ! add inertial loads
         Fnet_out = Body%F6net
                 
      else if (Body%typeNum == 2) then  ! pinned coupled body  

         if (p%inertialF == 1) then      ! include inertial components  
            ! inertial loads ... from input translational ... and solved rotational ... acceleration
            F6_iner(1:3)  = -MATMUL(Body%M(1:3,1:3), Body%a6(1:3)) - MATMUL(Body%M(1:3,4:6), Body%a6(4:6))
         elseif (p%inertialF == 2) then
            F6_iner(1:3)  = -MATMUL(Body%M(1:3,1:3), Body%a6(1:3)) - MATMUL(Body%M(1:3,4:6), Body%a6(4:6))
            if (t < p%inertialF_rampT) F6_iner = F6_iner * t / p%inertialF_rampT
         else 
            F6_iner(1:3) = 0.0
         endif
         
         Body%F6net(1:3) = Body%F6net(1:3) + F6_iner(1:3)     ! add translational inertial loads
         Body%F6net(4:6) = 0.0_DbKi
         Fnet_out = Body%F6net

      else
         Call WrScr( "ERROR, Body_GetCoupledForce called for wrong (non-coupled) body type in MoorDyn!")
      end if
   
   END SUBROUTINE Body_GetCoupledForce
   !--------------------------------------------------------------
   


   ! this function handles assigning a point to a body
   !--------------------------------------------------------------
   SUBROUTINE Body_AddPoint(Body, pointID, coords)

      Type(MD_Body),      INTENT(INOUT)  :: Body        ! the Point object
      Integer(IntKi),     INTENT(IN   )  :: pointID
      REAL(DbKi),         INTENT(IN   )  :: coords(3)


      IF (wordy > 0) Print*, "P", pointID, "->B", Body%IdNum
      
      IF(Body%nAttachedC < 30) THEN                ! this is currently just a maximum imposed by a fixed array size.  could be improved.
         Body%nAttachedC = Body%nAttachedC + 1     ! increment the number pointed
         Body%AttachedC(Body%nAttachedC) = pointID
         Body%rPointRel(:,Body%nAttachedC) = coords  ! store relative position of point on body
      ELSE
         call WrScr("too many Points attached to Body "//trim(num2lstr(Body%IdNum))//" in MoorDyn!")
      END IF

   END SUBROUTINE Body_AddPoint


   ! this function handles assigning a rod to a body
   !--------------------------------------------------------------
   SUBROUTINE Body_AddRod(Body, rodID, coords)

      Type(MD_Body),      INTENT(INOUT)  :: Body        ! the Point object
      Integer(IntKi),     INTENT(IN   )  :: rodID
      REAL(DbKi),         INTENT(IN   )  :: coords(6)  ! positions of rod ends A and B relative to body
      
      REAL(DbKi)                         :: tempUnitVec(3)
      REAL(DbKi)                         :: dummyLength

      IF (wordy > 0) Print*, "R", rodID, "->B", Body%IdNum
      
      IF(Body%nAttachedR < 30) THEN                ! this is currently just a maximum imposed by a fixed array size.  could be improved.
         Body%nAttachedR = Body%nAttachedR + 1     ! increment the number connected
         
         ! store rod ID
         Body%AttachedR(Body%nAttachedR) = rodID   
         
         ! store Rod end A relative position and unit vector from end A to B
         CALL UnitVector(coords(1:3), coords(4:6), tempUnitVec, dummyLength)
         Body%r6RodRel(1:3, Body%nAttachedR) = coords(1:3)
         Body%r6RodRel(4:6, Body%nAttachedR) = tempUnitVec
         
      ELSE
         call WrScr("too many rods attached to Body "//trim(num2lstr(Body%IdNum))//" in MoorDyn")
      END IF

   END SUBROUTINE Body_AddRod



END MODULE MoorDyn_Body
