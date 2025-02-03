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
MODULE MoorDyn_Point

   USE MoorDyn_Types
   USE MoorDyn_IO
   USE NWTC_Library
   USE MoorDyn_Misc
   USE MoorDyn_Line,           only : Line_SetEndKinematics, Line_GetEndStuff
   
   IMPLICIT NONE

   PRIVATE

   INTEGER(IntKi), PARAMETER            :: wordy = 0   ! verbosity level. >1 = more console output

   PUBLIC :: Point_Initialize
   PUBLIC :: Point_SetKinematics
   PUBLIC :: Point_SetState
   PUBLIC :: Point_GetStateDeriv
   PUBLIC :: Point_DoRHS
   PUBLIC :: Point_GetCoupledForce
   PUBLIC :: Point_GetNetForceAndMass
   PUBLIC :: Point_AddLine
   PUBLIC :: Point_RemoveLine
   

CONTAINS


   !--------------------------------------------------------------
   SUBROUTINE Point_Initialize(Point, states, m)

      Type(MD_Point), INTENT(INOUT)  :: Point        ! the Point object
      Real(DbKi),       INTENT(INOUT)  :: states(6)      ! state vector section for this Point
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m          ! passing along all mooring objects

      INTEGER(IntKi)                   :: l


      if (Point%typeNum == 0) then  ! error check
      
         ! pass kinematics to any attached lines so they have initial positions at this initialization stage
         DO l=1,Point%nAttached
            IF (wordy > 1) print *, "Point ",  Point%IdNum, " setting end kinematics of line ", Point%attached(l), " to ", Point%r
            CALL Line_SetEndKinematics(m%LineList(Point%attached(l)), Point%r, Point%rd, 0.0_DbKi, Point%Top(l))
         END DO


         ! assign initial node kinematics to state vector
         states(4:6) = Point%r
         states(1:3) = Point%rd
         
         
         IF (wordy > 0) print *, "Initialized Point ", Point%IdNum
      
      else 
         CALL WrScr("   Error: wrong Point type given to Point_Initialize for number "//trim(Int2Lstr(Point%idNum)))
      end if
      
   END SUBROUTINE Point_Initialize
   !--------------------------------------------------------------


   !--------------------------------------------------------------
   SUBROUTINE Point_SetKinematics(Point, r_in, rd_in, a_in, t, m)

      Type(MD_Point), INTENT(INOUT)  :: Point        ! the Point object
      Real(DbKi),       INTENT(IN   )  :: r_in( 3)       ! position
      Real(DbKi),       INTENT(IN   )  :: rd_in(3)       ! velocity
      Real(DbKi),       INTENT(IN   )  :: a_in(3)        ! acceleration (only used for coupled points)
      Real(DbKi),       INTENT(IN   )  :: t              ! instantaneous time
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m          ! passing along all mooring objects
      

      INTEGER(IntKi)                   :: l

      ! store current time
      Point%time = t

      
                     
      ! set position and velocity
      Point%r  = r_in
      Point%rd = rd_in
      Point%a = a_in
               
      ! pass latest kinematics to any attached lines
      DO l=1,Point%nAttached
         CALL Line_SetEndKinematics(m%LineList(Point%attached(l)), Point%r, Point%rd, t, Point%Top(l))
      END DO
   
   
         
   END SUBROUTINE Point_SetKinematics
   !--------------------------------------------------------------

   !--------------------------------------------------------------
   SUBROUTINE Point_SetState(Point, X, t, m)

      Type(MD_Point),      INTENT(INOUT)  :: Point        ! the Point object
      Real(DbKi),            INTENT(IN   )  :: X(:)           ! state vector section for this line
      Real(DbKi),            INTENT(IN   )  :: t              ! instantaneous time
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m              ! passing along all mooring objects

      INTEGER(IntKi)                        :: l              ! index of segments or nodes along line
      INTEGER(IntKi)                        :: J              ! index
   

      ! store current time
      Point%time = t
      
      ! from state values, get r and rdot values
      DO J=1,3
         Point%r( J) = X(3 + J)  ! get positions
         Point%rd(J) = X(    J)  ! get velocities
      END DO
           
     ! pass latest kinematics to any attached lines
      DO l=1,Point%nAttached
         CALL Line_SetEndKinematics(m%LineList(Point%attached(l)), Point%r, Point%rd, t, Point%Top(l))
      END DO
      
   END SUBROUTINE Point_SetState
   !--------------------------------------------------------------

   !--------------------------------------------------------------
   SUBROUTINE Point_GetStateDeriv(Point, Xd, m, p)

      Type(MD_Point),      INTENT(INOUT)  :: Point          ! the Point object
      Real(DbKi),            INTENT(INOUT)  :: Xd(:)            ! state derivative vector section for this line
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m              ! passing along all mooring objects
      
      TYPE(MD_ParameterType),INTENT(IN   )  :: p                ! Parameters
      
      !TYPE(MD_MiscVarType), INTENT(INOUT)  :: m       ! misc/optimization variables

      !INTEGER(IntKi)             :: l         ! index of attached lines
      INTEGER(IntKi)                        :: J                ! index
!      INTEGER(IntKi)                        :: K                ! index      
!      Real(DbKi)                            :: Sum1             ! for adding things
      
      Real(DbKi)                            :: S(3,3)           ! inverse mass matrix


      CALL Point_DoRHS(Point, m, p)
                  
!      // solve for accelerations in [M]{a}={f} using LU decomposition
!      double M_tot[9];                     // serialize total mass matrix for easy processing
!      for (int I=0; I<3; I++) for (int J=0; J<3; J++) M_tot[3*I+J]=M[I][J];
!      double LU[9];                        // serialized matrix that will hold LU matrices combined
!      Crout(3, M_tot, LU);                  // perform LU decomposition on mass matrix
!      double acc[3];                        // acceleration vector to solve for
!      solveCrout(3, LU, Fnet, acc);     // solve for acceleration vector

      ! solve for accelerations in [M]{a}={f} using LU decomposition
!      CALL LUsolve(6, M_out, LU_temp, Fnet_out, y_temp, acc)
   
                  
      ! invert node mass matrix
      CALL Inverse3by3(S, Point%M)

      ! accelerations 
      Point%a = MATMUL(S, Point%Fnet)

      ! fill in state derivatives
      Xd(4:6) = Point%rd           ! dxdt = V    (velocities)
      Xd(1:3) = Point%a            ! dVdt = RHS * A  (accelerations)
      

      ! check for NaNs
      DO J = 1, 6
         IF (Is_NaN(Xd(J))) THEN
            CALL WrScr("NaN detected at time "//trim(Num2LStr(Point%time))//" in Point "//trim(Int2LStr(Point%IdNum))//" in MoorDyn.")
            IF (wordy > 1) print *, "state derivatives:"
            IF (wordy > 1) print *, Xd
            EXIT
         END IF
      END DO

   END SUBROUTINE Point_GetStateDeriv
   !--------------------------------------------------------------

   !--------------------------------------------------------------
   SUBROUTINE Point_DoRHS(Point, m, p)

      Type(MD_Point),      INTENT(INOUT)  :: Point ! the Point object      
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m       ! passing along all mooring objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p       ! Parameters
      
      !TYPE(MD_MiscVarType), INTENT(INOUT)  :: m       ! misc/optimization variables

      INTEGER(IntKi)             :: l         ! index of attached lines
!      INTEGER(IntKi)             :: I         ! index
      INTEGER(IntKi)             :: J         ! index
!      INTEGER(IntKi)             :: K         ! index

      Real(DbKi)                 :: Fnet_i(3) ! force from an attached line
      Real(DbKi)                 :: Moment_dummy(3) ! dummy vector to hold unused line end moments
      Real(DbKi)                 :: M_i(3,3)  ! mass from an attached line


      ! start with the Point's own forces including buoyancy and weight, and its own mass
      Point%Fnet(1) = Point%pointFX
      Point%Fnet(2) = Point%pointFY
      Point%Fnet(3) = Point%pointFZ + Point%pointV*p%rhoW*p%g - Point%pointM*p%g

      Point%M    = 0.0_DbKi  ! clear (zero) the point mass matrix
      
      DO J = 1,3
        Point%M   (J,J) = Point%pointM  ! set the diagonals to the self-mass (to start with)
      END DO


   !   print *, "point number", Point%IdNum
   !   print *, "attached lines: ", Point%attached
   !   print *, "size of line list" , size(m%LineList)

      ! loop through attached lines, adding force and mass contributions
      DO l=1,Point%nAttached
         
      !   print *, "  l", l
      !   print *, Point%attached(l)
      !   print *, m%LineList(Point%attached(l))%Fnet
      !   
      !   
      !   print *, "  attached line ID", m%LineList(Point%attached(l))%IdNum
         
         CALL Line_GetEndStuff(m%LineList(Point%attached(l)), Fnet_i, Moment_dummy, M_i, Point%Top(l))
         
         ! sum quantitites
         Point%Fnet = Point%Fnet + Fnet_i
         Point%M    = Point%M    + M_i
         
      END DO


      ! XXXWhen this sub is called, any self weight, buoyancy, or external forcing should have already been
      ! added by the calling subroutine.  The only thing left is any added mass or drag forces from the point (e.g. float)
      ! itself, which will be added below.XXX


   !   IF (EqualRealNos(t, 0.0_DbKi)) THEN  ! this is old: with current IC gen approach, we skip the first call to the line objects, because they're set AFTER the call to the points
   !
   !      DO J = 1,3
   !         Xd(3+J) = X(J)        ! velocities - these are unused in integration
   !         Xd(J) = 0.0_DbKi           ! accelerations - these are unused in integration
   !      END DO
   !   ELSE
   !      ! from state values, get r and rdot values
   !      DO J = 1,3
   !         Point%r(J)  = X(3 + J)   ! get positions
   !         Point%rd(J) = X(J)       ! get velocities
   !      END DO
   !   END IF
      

      ! add any added mass and drag forces from the Point body itself
      DO J = 1,3
         Point%Fnet(J)   = Point%Fnet(J) - 0.5 * p%rhoW * Point%rd(J) * abs(Point%rd(J)) * Point%pointCdA;  ! add drag forces - corrected Nov 24
         Point%M   (J,J) = Point%M   (J,J) + Point%pointV*p%rhoW*Point%pointCa;                               ! add added mass

      END DO

      ! Added user-defined external force and damping on point
      Point%Fnet = Point%Fnet + Point%Fext
      DO J = 1, 3
         Point%Fnet(J) = Point%Fnet(J) - Point%Blin(J) * Point%rd(J) - Point%Bquad(J) * ABS(Point%rd(J)) * Point%rd(J)
      END DO

      ! would this sub ever need to include the m*a inertial term?  Is it ever called for coupled points? <<<

   END SUBROUTINE Point_DoRHS
   !=====================================================================


   ! calculate the force including inertial loads on point that is coupled
   !--------------------------------------------------------------
   SUBROUTINE Point_GetCoupledForce(Point,  Fnet_out, m, p)
   
      Type(MD_Point),      INTENT(INOUT)  :: Point     ! the Point object
      Real(DbKi),            INTENT(  OUT)  :: Fnet_out(3) ! force and moment vector about rRef
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m              ! passing along all mooring objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p           ! Parameters

      Real(DbKi)                            :: F_iner(3)   ! inertial force

      IF (Point%typeNum == -1) then
         ! calculate forces and masses of point
         CALL Point_DoRHS(Point, m, p)

         ! add inertial loads as appropriate
         F_iner = -MATMUL(Point%M, Point%a)    ! inertial loads
         Fnet_out = Point%Fnet + F_iner          ! add inertial loads

      ELSE
         CALL WrScr("Point_GetCoupledForce called for wrong (uncoupled) Point type in MoorDyn!")
      END IF
      
   END SUBROUTINE Point_GetCoupledForce


   ! calculate the force and mass contributions of the point on the parent body (only for type 3 points?)
   !--------------------------------------------------------------
   SUBROUTINE Point_GetNetForceAndMass(Point, rRef, wRef, Fnet_out, M_out, m, p)
   
      Type(MD_Point),      INTENT(INOUT)  :: Point     ! the Point object
      Real(DbKi),            INTENT(IN   )  :: rRef(3)     ! global coordinates of reference point (i.e. the parent body)
      Real(DbKi),            INTENT(IN   )  :: wRef(3)     ! global angular velocities of reference point (i.e. the parent body)
      Real(DbKi),            INTENT(  OUT)  :: Fnet_out(6) ! force and moment vector about rRef
      Real(DbKi),            INTENT(  OUT)  :: M_out(6,6)  ! mass and inertia matrix about rRef
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m           ! passing along all mooring objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p           ! Parameters

      Real(DbKi)                            :: rRel(  3)   ! position of point relative to the body reference point (global orientation frame)
      Real(DbKi)                            :: Fcentripetal(3)        ! centripetal force
      Real(DbKi)                            :: Mcentripetal(3)        ! centripetal moment     


      CALL Point_DoRHS(Point, m, p)

      rRel = Point%r - rRef    ! vector from body reference point to node

      ! convert net force into 6dof force about body ref point
      CALL translateForce3to6DOF(rRel, Point%Fnet, Fnet_out)

      ! convert mass matrix to 6by6 mass matrix about body ref point
      CALL translateMass3to6DOF(rRel, Point%M, M_out)

      ! add in the centripetal force and moment on the body. If rRel is zero there will be no translational centripetal component
      Fcentripetal = - MATMUL(M_out(1:3,1:3), CROSS_PRODUCT(wRef, CROSS_PRODUCT(wRef,rRel)))
      Mcentripetal = - CROSS_PRODUCT(wRef, MATMUL(M_out(4:6,4:6), wRef))

      Fnet_out(1:3) = Fnet_out(1:3) + Fcentripetal
      Fnet_out(4:6) = Fnet_out(4:6) + Mcentripetal

   END SUBROUTINE Point_GetNetForceAndMass
   
   
 
 
   ! this function handles assigning a line to a connection node
   !--------------------------------------------------------------
   SUBROUTINE Point_AddLine(Point, lineID, TopOfLine)

      Type(MD_Point), INTENT (INOUT)   :: Point        ! the Point object
      Integer(IntKi),   INTENT( IN )     :: lineID
      Integer(IntKi),   INTENT( IN )     :: TopOfLine

      IF (wordy > 0) Print*, "L", lineID, "->C", Point%IdNum
      
      IF (Point%nAttached <10) THEN ! this is currently just a maximum imposed by a fixed array size.  could be improved.
         Point%nAttached = Point%nAttached + 1  ! add the line to the number connected
         Point%Attached(Point%nAttached) = lineID
         Point%Top(Point%nAttached) = TopOfLine  ! attached to line ... 1 = top/fairlead(end B), 0 = bottom/anchor(end A)
      ELSE
         call WrScr("Too many lines connected to Point "//trim(num2lstr(Point%IdNum))//" in MoorDyn!")
      END IF

   END SUBROUTINE Point_AddLine


   ! this function handles removing a line from a connection node
   !--------------------------------------------------------------
   SUBROUTINE Point_RemoveLine(Point, lineID, TopOfLine, rEnd, rdEnd)

      Type(MD_Point), INTENT (INOUT)   :: Point        ! the Point object
      Integer(IntKi),   INTENT( IN )     :: lineID
      Integer(IntKi),   INTENT(  OUT)    :: TopOfLine
      REAL(DbKi),       INTENT(INOUT)    :: rEnd(3)
      REAL(DbKi),       INTENT(INOUT)    :: rdEnd(3)
      
      Integer(IntKi)    :: l,m,J
      logical           :: found

      found = .false.
      
      DO l = 1,Point%nAttached    ! look through attached lines
      
         IF (Point%Attached(l) == lineID) THEN   ! if this is the line's entry in the attachment list
         
            TopOfLine = Point%Top(l);                ! record which end of the line was attached
            
            DO m = l,Point%nAttached 
            
               Point%Attached(m) = Point%Attached(m+1)  ! move subsequent line links forward one spot in the list to eliminate this line link
               Point%Top(     m) =      Point%Top(m+1) 
            
               Point%nAttached = Point%nAttached - 1                      ! reduce attached line counter by 1
            
               ! also pass back the kinematics at the end
               DO J = 1,3
                  rEnd( J) = Point%r( J)
                  rdEnd(J) = Point%rd(J)
               END DO
               
               call WrScr( "Detached line "//trim(num2lstr(lineID))//" from Point "//trim(num2lstr(Point%IdNum)))
               
               EXIT
            END DO
            
            found = .true.
         
         END IF
         
      END DO

      IF (.not. found) THEN   ! detect if line not found TODO: fix this, its wrong. If pointNnattached is oprginally 2, then it will be 1 after one run of the loop and l will also be 1
         CALL WrScr("Error: failed to find line to remove during RemoveLine call to Point "//trim(num2lstr(Point%IdNum))//". Line "//trim(num2lstr(lineID)))
      END IF
      
   END SUBROUTINE Point_RemoveLine



END MODULE MoorDyn_Point
