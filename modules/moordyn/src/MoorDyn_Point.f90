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

   PUBLIC :: Connect_Initialize
   PUBLIC :: Connect_SetKinematics
   PUBLIC :: Connect_SetState
   PUBLIC :: Connect_GetStateDeriv
   PUBLIC :: Connect_DoRHS
   PUBLIC :: Connect_GetCoupledForce
   PUBLIC :: Connect_GetNetForceAndMass
   PUBLIC :: Connect_AddLine
   PUBLIC :: Connect_RemoveLine
   

CONTAINS


   !--------------------------------------------------------------
   SUBROUTINE Connect_Initialize(Connect, states, m)

      Type(MD_Connect), INTENT(INOUT)  :: Connect        ! the Connection object
      Real(DbKi),       INTENT(INOUT)  :: states(6)      ! state vector section for this Connection
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m          ! passing along all mooring objects

      INTEGER(IntKi)                   :: l


      if (Connect%typeNum == 0) then  ! error check
      
         ! pass kinematics to any attached lines so they have initial positions at this initialization stage
         DO l=1,Connect%nAttached
            IF (wordy > 1) print *, "Connect ",  Connect%IdNum, " setting end kinematics of line ", Connect%attached(l), " to ", Connect%r
            CALL Line_SetEndKinematics(m%LineList(Connect%attached(l)), Connect%r, Connect%rd, 0.0_DbKi, Connect%Top(l))
         END DO


         ! assign initial node kinematics to state vector
         states(4:6) = Connect%r
         states(1:3) = Connect%rd
         
         
         IF (wordy > 0) print *, "Initialized Connection ", Connect%IdNum
      
      else 
         CALL WrScr("   Error: wrong Point type given to Connect_Initialize for number "//trim(Int2Lstr(Connect%idNum)))
      end if
      
   END SUBROUTINE Connect_Initialize
   !--------------------------------------------------------------


   !--------------------------------------------------------------
   SUBROUTINE Connect_SetKinematics(Connect, r_in, rd_in, a_in, t, m)

      Type(MD_Connect), INTENT(INOUT)  :: Connect        ! the Connection object
      Real(DbKi),       INTENT(IN   )  :: r_in( 3)       ! position
      Real(DbKi),       INTENT(IN   )  :: rd_in(3)       ! velocity
      Real(DbKi),       INTENT(IN   )  :: a_in(3)        ! acceleration (only used for coupled connects)
      Real(DbKi),       INTENT(IN   )  :: t              ! instantaneous time
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m          ! passing along all mooring objects
      

      INTEGER(IntKi)                   :: l

      ! store current time
      Connect%time = t

      
    !  if (Connect%typeNum==0) THEN ! anchor ( <<< to be changed/expanded) ... in MoorDyn F also used for coupled connections
                        
         ! set position and velocity
         Connect%r  = r_in
         Connect%rd = rd_in
         Connect%a = a_in
                 
         ! pass latest kinematics to any attached lines
         DO l=1,Connect%nAttached
            CALL Line_SetEndKinematics(m%LineList(Connect%attached(l)), Connect%r, Connect%rd, t, Connect%Top(l))
         END DO
      
     ! else
     !    
     !    PRINT*,"Error: setKinematics called for wrong Connection type. Connection ", Connect%IdNum, " type ", Connect%typeNum
         
   !  END IF
      
         
   END SUBROUTINE Connect_SetKinematics
   !--------------------------------------------------------------

   !--------------------------------------------------------------
   SUBROUTINE Connect_SetState(Connect, X, t, m)

      Type(MD_Connect),      INTENT(INOUT)  :: Connect        ! the Connection object
      Real(DbKi),            INTENT(IN   )  :: X(:)           ! state vector section for this line
      Real(DbKi),            INTENT(IN   )  :: t              ! instantaneous time
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m              ! passing along all mooring objects

      INTEGER(IntKi)                        :: l              ! index of segments or nodes along line
      INTEGER(IntKi)                        :: J              ! index
   

      ! store current time
      Connect%time = t
      
      ! from state values, get r and rdot values
      DO J=1,3
         Connect%r( J) = X(3 + J)  ! get positions
         Connect%rd(J) = X(    J)  ! get velocities
      END DO
           
     ! pass latest kinematics to any attached lines
      DO l=1,Connect%nAttached
         CALL Line_SetEndKinematics(m%LineList(Connect%attached(l)), Connect%r, Connect%rd, t, Connect%Top(l))
      END DO
      
   END SUBROUTINE Connect_SetState
   !--------------------------------------------------------------

   !--------------------------------------------------------------
   SUBROUTINE Connect_GetStateDeriv(Connect, Xd, m, p)

      Type(MD_Connect),      INTENT(INOUT)  :: Connect          ! the Connection object
      Real(DbKi),            INTENT(INOUT)  :: Xd(:)            ! state derivative vector section for this line
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m              ! passing along all mooring objects
      
      TYPE(MD_ParameterType),INTENT(IN   )  :: p                ! Parameters
      
      !TYPE(MD_MiscVarType), INTENT(INOUT)  :: m       ! misc/optimization variables

      !INTEGER(IntKi)             :: l         ! index of attached lines
      INTEGER(IntKi)                        :: J                ! index
!      INTEGER(IntKi)                        :: K                ! index      
!      Real(DbKi)                            :: Sum1             ! for adding things
      
      Real(DbKi)                            :: S(3,3)           ! inverse mass matrix


      CALL Connect_DoRHS(Connect, m, p)
                  
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
      CALL Inverse3by3(S, Connect%M)

      ! accelerations 
      Connect%a = MATMUL(S, Connect%Fnet)

      ! fill in state derivatives
      Xd(4:6) = Connect%rd           ! dxdt = V    (velocities)
      Xd(1:3) = Connect%a            ! dVdt = RHS * A  (accelerations)
      

      ! check for NaNs
      DO J = 1, 6
         IF (Is_NaN(Xd(J))) THEN
            CALL WrScr("NaN detected at time "//trim(Num2LStr(Connect%time))//" in Point "//trim(Int2LStr(Connect%IdNum))//" in MoorDyn.")
            IF (wordy > 1) print *, "state derivatives:"
            IF (wordy > 1) print *, Xd
            EXIT
         END IF
      END DO

   END SUBROUTINE Connect_GetStateDeriv
   !--------------------------------------------------------------

   !--------------------------------------------------------------
   SUBROUTINE Connect_DoRHS(Connect, m, p)

      Type(MD_Connect),      INTENT(INOUT)  :: Connect ! the Connection object      
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


      ! start with the Connection's own forces including buoyancy and weight, and its own mass
      Connect%Fnet(1) = Connect%conFX
      Connect%Fnet(2) = Connect%conFY
      Connect%Fnet(3) = Connect%conFZ + Connect%conV*p%rhoW*p%g - Connect%conM*p%g

      Connect%M    = 0.0_DbKi  ! clear (zero) the connect mass matrix
      
      DO J = 1,3
        Connect%M   (J,J) = Connect%conM  ! set the diagonals to the self-mass (to start with)
      END DO


   !   print *, "connection number", Connect%IdNum
   !   print *, "attached lines: ", Connect%attached
   !   print *, "size of line list" , size(m%LineList)

      ! loop through attached lines, adding force and mass contributions
      DO l=1,Connect%nAttached
         
      !   print *, "  l", l
      !   print *, Connect%attached(l)
      !   print *, m%LineList(Connect%attached(l))%Fnet
      !   
      !   
      !   print *, "  attached line ID", m%LineList(Connect%attached(l))%IdNum
         
         CALL Line_GetEndStuff(m%LineList(Connect%attached(l)), Fnet_i, Moment_dummy, M_i, Connect%Top(l))
         
         ! sum quantitites
         Connect%Fnet = Connect%Fnet + Fnet_i
         Connect%M    = Connect%M    + M_i
         
      END DO


      ! XXXWhen this sub is called, any self weight, buoyancy, or external forcing should have already been
      ! added by the calling subroutine.  The only thing left is any added mass or drag forces from the connection (e.g. float)
      ! itself, which will be added below.XXX


   !   IF (EqualRealNos(t, 0.0_DbKi)) THEN  ! this is old: with current IC gen approach, we skip the first call to the line objects, because they're set AFTER the call to the connects
   !
   !      DO J = 1,3
   !         Xd(3+J) = X(J)        ! velocities - these are unused in integration
   !         Xd(J) = 0.0_DbKi           ! accelerations - these are unused in integration
   !      END DO
   !   ELSE
   !      ! from state values, get r and rdot values
   !      DO J = 1,3
   !         Connect%r(J)  = X(3 + J)   ! get positions
   !         Connect%rd(J) = X(J)       ! get velocities
   !      END DO
   !   END IF
      

      ! add any added mass and drag forces from the Connect body itself
      DO J = 1,3
         Connect%Fnet(J)   = Connect%Fnet(J) - 0.5 * p%rhoW * Connect%rd(J) * abs(Connect%rd(J)) * Connect%conCdA;  ! add drag forces - corrected Nov 24
         Connect%M   (J,J) = Connect%M   (J,J) + Connect%conV*p%rhoW*Connect%conCa;                               ! add added mass

      END DO
      
      ! would this sub ever need to include the m*a inertial term?  Is it ever called for coupled connects? <<<

   END SUBROUTINE Connect_DoRHS
   !=====================================================================


   ! calculate the force including inertial loads on connect that is coupled
   !--------------------------------------------------------------
   SUBROUTINE Connect_GetCoupledForce(Connect,  Fnet_out, m, p)
   
      Type(MD_Connect),      INTENT(INOUT)  :: Connect     ! the Connect object
      Real(DbKi),            INTENT(  OUT)  :: Fnet_out(3) ! force and moment vector about rRef
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m              ! passing along all mooring objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p           ! Parameters

      Real(DbKi)                            :: F_iner(3)   ! inertial force

      IF (Connect%typeNum == -1) then
         ! calculate forces and masses of connect
         CALL Connect_DoRHS(Connect, m, p)

         ! add inertial loads as appropriate
         F_iner = -MATMUL(Connect%M, Connect%a)    ! inertial loads
         Fnet_out = Connect%Fnet + F_iner          ! add inertial loads

      ELSE
         CALL WrScr("Connect_GetCoupledForce called for wrong (uncoupled) Point type in MoorDyn!")
      END IF
      
   END SUBROUTINE Connect_GetCoupledForce


   ! calculate the force and mass contributions of the connect on the parent body (only for type 3 connects?)
   !--------------------------------------------------------------
   SUBROUTINE Connect_GetNetForceAndMass(Connect, rRef, Fnet_out, M_out, m, p)
   
      Type(MD_Connect),      INTENT(INOUT)  :: Connect     ! the Connect object
      Real(DbKi),            INTENT(IN   )  :: rRef(3)     ! global coordinates of reference point (i.e. the parent body)
      Real(DbKi),            INTENT(  OUT)  :: Fnet_out(6) ! force and moment vector about rRef
      Real(DbKi),            INTENT(  OUT)  :: M_out(6,6)  ! mass and inertia matrix about rRef
      TYPE(MD_MiscVarType),  INTENT(INOUT)  :: m           ! passing along all mooring objects
      TYPE(MD_ParameterType),INTENT(IN   )  :: p           ! Parameters

      Real(DbKi)                            :: rRel(  3)   ! position of connection relative to the body reference point (global orientation frame)


      CALL Connect_DoRHS(Connect, m, p)

      rRel = Connect%r - rRef    ! vector from body reference point to node

      ! convert net force into 6dof force about body ref point
      CALL translateForce3to6DOF(rRel, Connect%Fnet, Fnet_out)
      
      ! convert mass matrix to 6by6 mass matrix about body ref point
      CALL translateMass3to6DOF(rRel, Connect%M, M_out)

   END SUBROUTINE Connect_GetNetForceAndMass
   
   
 
 
   ! this function handles assigning a line to a connection node
   !--------------------------------------------------------------
   SUBROUTINE Connect_AddLine(Connect, lineID, TopOfLine)

      Type(MD_Connect), INTENT (INOUT)   :: Connect        ! the Connection object
      Integer(IntKi),   INTENT( IN )     :: lineID
      Integer(IntKi),   INTENT( IN )     :: TopOfLine

      IF (wordy > 0) Print*, "L", lineID, "->C", Connect%IdNum
      
      IF (Connect%nAttached <10) THEN ! this is currently just a maximum imposed by a fixed array size.  could be improved.
         Connect%nAttached = Connect%nAttached + 1  ! add the line to the number connected
         Connect%Attached(Connect%nAttached) = lineID
         Connect%Top(Connect%nAttached) = TopOfLine  ! attached to line ... 1 = top/fairlead(end B), 0 = bottom/anchor(end A)
      ELSE
         Print*, "Too many lines connected to Point ", Connect%IdNum, " in MoorDyn!"
      END IF

   END SUBROUTINE Connect_AddLine


   ! this function handles removing a line from a connection node
   !--------------------------------------------------------------
   SUBROUTINE Connect_RemoveLine(Connect, lineID, TopOfLine, rEnd, rdEnd)

      Type(MD_Connect), INTENT (INOUT)   :: Connect        ! the Connection object
      Integer(IntKi),   INTENT( IN )     :: lineID
      Integer(IntKi),   INTENT(  OUT)    :: TopOfLine
      REAL(DbKi),       INTENT(INOUT)    :: rEnd(3)
      REAL(DbKi),       INTENT(INOUT)    :: rdEnd(3)
      
      Integer(IntKi)    :: l,m,J
      
      DO l = 1,Connect%nAttached    ! look through attached lines
      
         IF (Connect%Attached(l) == lineID) THEN   ! if this is the line's entry in the attachment list
         
            TopOfLine = Connect%Top(l);                ! record which end of the line was attached
            
            DO m = l,Connect%nAttached-1 
            
               Connect%Attached(m) = Connect%Attached(m+1)  ! move subsequent line links forward one spot in the list to eliminate this line link
               Connect%Top(     m) =      Connect%Top(m+1) 
            
               Connect%nAttached = Connect%nAttached - 1                      ! reduce attached line counter by 1
            
               ! also pass back the kinematics at the end
               DO J = 1,3
                  rEnd( J) = Connect%r( J)
                  rdEnd(J) = Connect%rd(J)
               END DO
               
               print*, "Detached line ", lineID, " from Connection ", Connect%IdNum
               
               EXIT
            END DO
            
            IF (l == Connect%nAttached) THEN   ! detect if line not found
               print *, "Error: failed to find line to remove during removeLineFromConnect call to connection ", Connect%IdNum, ". Line ", lineID
            END IF
         
         END IF
         
      END DO
      
   END SUBROUTINE Connect_RemoveLine



END MODULE MoorDyn_Point
