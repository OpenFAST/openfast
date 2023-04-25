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
MODULE MoorDyn_Line

   USE MoorDyn_Types
   USE MoorDyn_IO
   USE NWTC_Library
   USE MoorDyn_Misc
   
   IMPLICIT NONE

   PRIVATE

   INTEGER(IntKi), PARAMETER            :: wordy = 0   ! verbosity level. >1 = more console output

   PUBLIC :: SetupLine
   PUBLIC :: Line_Initialize
   PUBLIC :: Line_SetState
   PUBLIC :: Line_GetStateDeriv
   PUBLIC :: Line_SetEndKinematics
   PUBLIC :: Line_GetEndStuff
   PUBLIC :: Line_GetEndSegmentInfo
   PUBLIC :: Line_SetEndOrientation
   
   

CONTAINS


   !-----------------------------------------------------------------------
   ! >>>>>>>>>>>>>> rename/reorganize this subroutine >>>>>>>>>>>>>
   SUBROUTINE SetupLine (Line, LineProp, p, ErrStat, ErrMsg)
      ! allocate arrays in line object

      TYPE(MD_Line), INTENT(INOUT)       :: Line          ! the single line object of interest
      TYPE(MD_LineProp), INTENT(INOUT)   :: LineProp      ! the single line property set for the line of interest
      TYPE(MD_ParameterType), INTENT(IN   ) :: p       ! Parameters
      INTEGER,       INTENT(   INOUT )   :: ErrStat       ! returns a non-zero value when an error occurs
      CHARACTER(*),  INTENT(   INOUT )   :: ErrMsg        ! Error message if ErrStat /= ErrID_None

      INTEGER(4)                         :: I, J          ! Generic index
      INTEGER(IntKi)                     :: N
      REAL(DbKi)                         :: temp


      N = Line%N  ! number of segments in this line (for code readability)

      ! -------------- save some section properties to the line object itself -----------------

      Line%d   = LineProp%d
      Line%rho = LineProp%w/(Pi/4.0 * Line%d * Line%d)
      
      Line%EA   = LineProp%EA
      ! note: Line%BA is set later
      Line%EA_D = LineProp%EA_D
      Line%BA_D = LineProp%BA_D
      Line%EI   = LineProp%EI  !<<< for bending stiffness
      
      Line%Can   = LineProp%Can
      Line%Cat   = LineProp%Cat
      Line%Cdn   = LineProp%Cdn
      Line%Cdt   = LineProp%Cdt      
      
      ! copy over elasticity data
      Line%ElasticMod = LineProp%ElasticMod
      
      Line%nEApoints = LineProp%nEApoints
      DO I = 1,Line%nEApoints
         Line%stiffXs(I) = LineProp%stiffXs(I)
         Line%stiffYs(I) = LineProp%stiffYs(I)  ! note: this does not convert to E (not EA) like done in C version
      END DO
      
      Line%nBApoints = LineProp%nBApoints
      DO I = 1,Line%nBApoints
         Line%dampXs(I) = LineProp%dampXs(I)
         Line%dampYs(I) = LineProp%dampYs(I)
      END DO
      
      Line%nEIpoints = LineProp%nEIpoints
      DO I = 1,Line%nEIpoints
         Line%bstiffXs(I) = LineProp%bstiffXs(I)
         Line%bstiffYs(I) = LineProp%bstiffYs(I)                     ! copy over
      END DO
      
      
      
      ! Specify specific internal damping coefficient (BA) for this line.
      ! Will be equal to inputted BA of LineType if input value is positive.
      ! If input value is negative, it is considered to be desired damping ratio (zeta)
      ! from which the line's BA can be calculated based on the segment natural frequency.
      IF (LineProp%BA < 0) THEN
         ! - we assume desired damping coefficient is zeta = -LineProp%BA
         ! - highest axial vibration mode of a segment is wn = sqrt(k/m) = 2N/UnstrLen*sqrt(EA/w)
         Line%BA = -LineProp%BA * Line%UnstrLen / Line%N * SQRT(LineProp%EA * LineProp%w)
         IF (wordy > 1) print *, 'Based on zeta, BA set to ', Line%BA
         
         IF (wordy > 1) print *, 'Negative BA input detected, treating as -zeta.  For zeta = ', -LineProp%BA, ', setting BA to ', Line%BA
         
      ELSE
         Line%BA = LineProp%BA
         IF (wordy > 1) temp = Line%N * Line%BA / Line%UnstrLen * SQRT(1.0/(LineProp%EA * LineProp%w))
         IF (wordy > 1) print *, 'BA set as input to ', Line%BA, '. Corresponding zeta is ', temp
      END IF
      
      !temp = 2*Line%N / Line%UnstrLen * sqrt( LineProp%EA / LineProp%w) / TwoPi
      !print *, 'Segment natural frequency is ', temp, ' Hz'
      
      
      !print *, "Line ElasticMod is ", Line%ElasticMod
      !print *, "EA (static value) is", Line%EA
      !print *, "EA_D              is", Line%EA_D 
      !print *, "BA                is", Line%BA 
      !print *, "BA_D              is", Line%BA_D 
      

      ! allocate node positions and velocities (NOTE: these arrays start at ZERO)
      ALLOCATE ( Line%r(3, 0:N), Line%rd(3, 0:N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating r and rd arrays.'
         !CALL CleanUp()
         RETURN
      END IF
      
      ! if using viscoelastic model, allocate additional state quantities
      if (Line%ElasticMod == 2) then
         ALLOCATE ( Line%dl_1(N), STAT = ErrStat )
         IF ( ErrStat /= ErrID_None ) THEN
            ErrMsg  = ' Error allocating dl_1 array.'
            !CALL CleanUp()
            RETURN
         END IF
         ! initialize to zero
         Line%dl_1 = 0.0_DbKi
      end if
      
      ! allocate node and segment tangent vectors
      ALLOCATE ( Line%q(3, 0:N), Line%qs(3, N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating q or qs array.'
         !CALL CleanUp()
         RETURN
      END IF

      ! allocate segment scalar quantities
      ALLOCATE ( Line%l(N), Line%ld(N), Line%lstr(N), Line%lstrd(N), Line%Kurv(0:N), Line%V(N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating segment scalar quantity arrays.'
         !CALL CleanUp()
         RETURN
      END IF
      Line%Kurv = 0.0_DbKi

      ! assign values for l, ld, and V
      DO J=1,N
         Line%l(J) = Line%UnstrLen/REAL(N, DbKi)
         Line%ld(J)= 0.0_DbKi
         Line%V(J) = Line%l(J)*0.25*Pi*LineProp%d*LineProp%d
      END DO
      
      ! allocate water related vectors
      ALLOCATE ( Line%U(3, 0:N), Line%Ud(3, 0:N), Line%zeta(0:N), Line%PDyn(0:N), STAT = ErrStat )
      ! set to zero initially (important of wave kinematics are not being used)
      Line%U    = 0.0_DbKi
      Line%Ud   = 0.0_DbKi
      Line%zeta = 0.0_DbKi
      Line%PDyn = 0.0_DbKi

      ! allocate segment tension and internal damping force vectors
      ALLOCATE ( Line%T(3, N), Line%Td(3, N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating T and Td arrays.'
         !CALL CleanUp()
         RETURN
      END IF

      ! allocate node force vectors
      ALLOCATE ( Line%W(3, 0:N), Line%Dp(3, 0:N), Line%Dq(3, 0:N), Line%Ap(3, 0:N), &
         Line%Aq(3, 0:N), Line%B(3, 0:N), Line%Bs(3, 0:N), Line%Fnet(3, 0:N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating node force arrays.'
         !CALL CleanUp()
         RETURN
      END IF

      ! set gravity and bottom contact forces to zero initially (because the horizontal components should remain at zero)
      Line%W = 0.0_DbKi
      Line%B = 0.0_DbKi
      
      ! allocate mass and inverse mass matrices for each node (including ends)
      ALLOCATE ( Line%S(3, 3, 0:N), Line%M(3, 3, 0:N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating T and Td arrays.'
         !CALL CleanUp()
         RETURN
      END IF
      
      
      if (p%writeLog > 1) then
         write(p%UnLog, '(A)') "  - Line"//trim(num2lstr(Line%IdNum))
         write(p%UnLog, '(A)') "    ID: "//trim(num2lstr(Line%IdNum))
         write(p%UnLog, '(A)') "    UnstrLen: "//trim(num2lstr(Line%UnstrLen))
         write(p%UnLog, '(A)') "    N   : "//trim(num2lstr(Line%N   ))
         write(p%UnLog, '(A)') "    d   : "//trim(num2lstr(Line%d   ))
         write(p%UnLog, '(A)') "    rho : "//trim(num2lstr(Line%rho ))
         write(p%UnLog, '(A)') "    E   : "//trim(num2lstr(Line%EA  ))
         write(p%UnLog, '(A)') "    EI  : "//trim(num2lstr(Line%EI  ))
         !write(p%UnLog, '(A)') "    BAin: "//trim(num2lstr(Line%BAin))
         write(p%UnLog, '(A)') "    Can : "//trim(num2lstr(Line%Can ))
         write(p%UnLog, '(A)') "    Cat : "//trim(num2lstr(Line%Cat ))
         write(p%UnLog, '(A)') "    Cdn : "//trim(num2lstr(Line%Cdn ))
         write(p%UnLog, '(A)') "    Cdt : "//trim(num2lstr(Line%Cdt ))
         !write(p%UnLog, '(A)') "    ww_l: " << ( (rho - env->rho_w)*(pi/4.*d*d) )*9.81 << endl;	
      end if
    
      
      ! need to add cleanup sub <<<


   END SUBROUTINE SetupLine
   !--------------------------------------------------------------





   !----------------------------------------------------------------------------------------=======
   SUBROUTINE Line_Initialize (Line, LineProp, rhoW, ErrStat, ErrMsg)
      ! calculate initial profile of the line using quasi-static model

      TYPE(MD_Line),     INTENT(INOUT)       :: Line        ! the single line object of interest
      TYPE(MD_LineProp), INTENT(INOUT)       :: LineProp    ! the single line property set for the line of interest
      REAL(DbKi),        INTENT(IN)          :: rhoW
      INTEGER,           INTENT(   INOUT )   :: ErrStat     ! returns a non-zero value when an error occurs
      CHARACTER(*),      INTENT(   INOUT )   :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      REAL(DbKi)                             :: COSPhi      ! Cosine of the angle between the xi-axis of the inertia frame and the X-axis of the local coordinate system of the current mooring line (-)
      REAL(DbKi)                             :: SINPhi      ! Sine   of the angle between the xi-axis of the inertia frame and the X-axis of the local coordinate system of the current mooring line (-)
      REAL(DbKi)                             :: XF          ! Horizontal distance between anchor and fairlead of the current mooring line (meters)
      REAL(DbKi)                             :: ZF          ! Vertical   distance between anchor and fairlead of the current mooring line (meters)
      INTEGER(4)                             :: I           ! Generic index
      INTEGER(4)                             :: J           ! Generic index


      INTEGER(IntKi)                         :: ErrStat2      ! Error status of the operation
      CHARACTER(ErrMsgLen)                   :: ErrMsg2       ! Error message if ErrStat2 /= ErrID_None
      REAL(DbKi)                             :: WetWeight
      REAL(DbKi)                             :: SeabedCD = 0.0_DbKi
      REAL(DbKi)                             :: TenTol = 0.0001_DbKi
      REAL(DbKi), ALLOCATABLE                :: LSNodes(:)
      REAL(DbKi), ALLOCATABLE                :: LNodesX(:)
      REAL(DbKi), ALLOCATABLE                :: LNodesZ(:)
      INTEGER(IntKi)                         :: N


      N = Line%N ! for convenience

       ! try to calculate initial line profile using catenary routine (from FAST v.7)
       ! note: much of this function is adapted from the FAST source code

          ! Transform the fairlead location from the inertial frame coordinate system
          !   to the local coordinate system of the current line (this coordinate
          !   system lies at the current anchor, Z being vertical, and X directed from
          !   current anchor to the current fairlead).  Also, compute the orientation
          !   of this local coordinate system:

             XF         = SQRT( ( Line%r(1,N) - Line%r(1,0) )**2.0 + ( Line%r(2,N) - Line%r(2,0) )**2.0 )
             ZF         =         Line%r(3,N) - Line%r(3,0)

             IF ( XF == 0.0 )  THEN  ! .TRUE. if the current mooring line is exactly vertical; thus, the solution below is ill-conditioned because the orientation is undefined; so set it such that the tensions and nodal positions are only vertical
                COSPhi  = 0.0_DbKi
                SINPhi  = 0.0_DbKi
             ELSE                    ! The current mooring line must not be vertical; use simple trigonometry
                COSPhi  =       ( Line%r(1,N) - Line%r(1,0) )/XF
                SINPhi  =       ( Line%r(2,N) - Line%r(2,0) )/XF
             ENDIF

        WetWeight = LineProp%w - 0.25*Pi*LineProp%d*LineProp%d*rhoW

        !LineNodes = Line%N + 1  ! number of nodes in line for catenary model to worry about

        ! allocate temporary arrays for catenary routine
        ALLOCATE ( LSNodes(N+1), STAT = ErrStat )
        IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Error allocating LSNodes array.'
          CALL CleanUp()
          RETURN
        END IF

        ALLOCATE ( LNodesX(N+1), STAT = ErrStat )
        IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Error allocating LNodesX array.'
          CALL CleanUp()
          RETURN
        END IF

        ALLOCATE ( LNodesZ(N+1), STAT = ErrStat )
        IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Error allocating LNodesZ array.'
          CALL CleanUp()
          RETURN
        END IF

        ! Assign node arc length locations
        LSNodes(1) = 0.0_DbKi
        DO I=2,N
          LSNodes(I) = LSNodes(I-1) + Line%l(I-1)  ! note: l index is because line segment indices start at 1
        END DO
        LSNodes(N+1) = Line%UnstrLen  ! ensure the last node length isn't longer than the line due to numerical error

          ! Solve the analytical, static equilibrium equations for a catenary (or
          !   taut) mooring line with seabed interaction in order to find the
          !   horizontal and vertical tensions at the fairlead in the local coordinate
          !   system of the current line:
          ! NOTE: The values for the horizontal and vertical tensions at the fairlead
          !       from the previous time step are used as the initial guess values at
          !       at this time step (because the LAnchHTe(:) and LAnchVTe(:) arrays
          !       are stored in a module and thus their values are saved from CALL to
          !       CALL).


             CALL Catenary ( XF           , ZF          , Line%UnstrLen, LineProp%EA  , &
                             WetWeight    , SeabedCD,    TenTol,     (N+1)     , &
                             LSNodes, LNodesX, LNodesZ , ErrStat2, ErrMsg2)

      IF (ErrStat2 == ErrID_None) THEN ! if it worked, use it
          ! Transform the positions of each node on the current line from the local
          !   coordinate system of the current line to the inertial frame coordinate
          !   system:

          DO J = 0,N ! Loop through all nodes per line where the line position and tension can be output
             Line%r(1,J) = Line%r(1,0) + LNodesX(J+1)*COSPhi
             Line%r(2,J) = Line%r(2,0) + LNodesX(J+1)*SINPhi
             Line%r(3,J) = Line%r(3,0) + LNodesZ(J+1)
          ENDDO              ! J - All nodes per line where the line position and tension can be output


      ELSE ! if there is a problem with the catenary approach, just stretch the nodes linearly between fairlead and anchor

          !CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Line_Initialize')
          call WrScr(" Catenary solve of Line "//trim(Num2LStr(Line%IdNum))//" unsuccessful. Initializing as linear.")

!          print *, "Node positions: "

          DO J = 0,N ! Loop through all nodes per line where the line position and tension can be output
             Line%r(1,J) = Line%r(1,0) + (Line%r(1,N) - Line%r(1,0))*REAL(J, DbKi)/REAL(N, DbKi)
             Line%r(2,J) = Line%r(2,0) + (Line%r(2,N) - Line%r(2,0))*REAL(J, DbKi)/REAL(N, DbKi)
             Line%r(3,J) = Line%r(3,0) + (Line%r(3,N) - Line%r(3,0))*REAL(J, DbKi)/REAL(N, DbKi)
             
!             print*, Line%r(:,J)
          ENDDO
          
!          print*,"FYI line end A and B node coords are"
!          print*, Line%r(:,0)
!          print*, Line%r(:,N)
      ENDIF



      CALL CleanUp()  ! deallocate temporary arrays



   CONTAINS


      !-----------------------------------------------------------------------
      SUBROUTINE CleanUp()
           ! deallocate temporary arrays

           IF (ALLOCATED(LSNodes))  DEALLOCATE(LSNodes)
           IF (ALLOCATED(LNodesX))  DEALLOCATE(LNodesX)
           IF (ALLOCATED(LNodesZ))  DEALLOCATE(LNodesZ)

        END SUBROUTINE CleanUp
      !-----------------------------------------------------------------------


      !-----------------------------------------------------------------------
      SUBROUTINE Catenary ( XF_In, ZF_In, L_In  , EA_In, &
                            W_In , CB_In, Tol_In, N    , &
                            s_In , X_In , Z_In , ErrStat, ErrMsg    )

         ! This subroutine is copied from FAST v7 with minor modifications

         ! This routine solves the analytical, static equilibrium equations
         ! for a catenary (or taut) mooring line with seabed interaction.
         ! Stretching of the line is accounted for, but bending stiffness
         ! is not.  Given the mooring line properties and the fairlead
         ! position relative to the anchor, this routine finds the line
         ! configuration and tensions.  Since the analytical solution
         ! involves two nonlinear equations (XF and  ZF) in two unknowns
         ! (HF and VF), a Newton-Raphson iteration scheme is implemented in
         ! order to solve for the solution.  The values of HF and VF that
         ! are passed into this routine are used as the initial guess in
         ! the iteration.  The Newton-Raphson iteration is only accurate in
         ! double precision, so all of the input/output arguments are
         ! converteds to/from double precision from/to default precision.

         ! >>>> TO DO: streamline this function, if it's still to be used at all <<<<

         !     USE                             Precision


         IMPLICIT                        NONE


            ! Passed Variables:

         INTEGER(4), INTENT(IN   )    :: N                                               ! Number of nodes where the line position and tension can be output (-)

         REAL(DbKi), INTENT(IN   )    :: CB_In                                           ! Coefficient of seabed static friction drag (a negative value indicates no seabed) (-)
         REAL(DbKi), INTENT(IN   )    :: EA_In                                           ! Extensional stiffness of line (N)
     !    REAL(DbKi), INTENT(  OUT)    :: HA_In                                           ! Effective horizontal tension in line at the anchor   (N)
     !    REAL(DbKi), INTENT(INOUT)    :: HF_In                                           ! Effective horizontal tension in line at the fairlead (N)
         REAL(DbKi), INTENT(IN   )    :: L_In                                            ! Unstretched length of line (meters)
         REAL(DbKi), INTENT(IN   )    :: s_In     (N)                                    ! Unstretched arc distance along line from anchor to each node where the line position and tension can be output (meters)
     !    REAL(DbKi), INTENT(  OUT)    :: Te_In    (N)                                    ! Effective line tensions at each node (N)
         REAL(DbKi), INTENT(IN   )    :: Tol_In                                          ! Convergence tolerance within Newton-Raphson iteration specified as a fraction of tension (-)
     !    REAL(DbKi), INTENT(  OUT)    :: VA_In                                           ! Effective vertical   tension in line at the anchor   (N)
     !    REAL(DbKi), INTENT(INOUT)    :: VF_In                                           ! Effective vertical   tension in line at the fairlead (N)
         REAL(DbKi), INTENT(IN   )    :: W_In                                            ! Weight of line in fluid per unit length (N/m)
         REAL(DbKi), INTENT(  OUT)    :: X_In     (N)                                    ! Horizontal locations of each line node relative to the anchor (meters)
         REAL(DbKi), INTENT(IN   )    :: XF_In                                           ! Horizontal distance between anchor and fairlead (meters)
         REAL(DbKi), INTENT(  OUT)    :: Z_In     (N)                                    ! Vertical   locations of each line node relative to the anchor (meters)
         REAL(DbKi), INTENT(IN   )    :: ZF_In                                           ! Vertical   distance between anchor and fairlead (meters)
             INTEGER,                      INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs
             CHARACTER(*),                 INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None


            ! Local Variables:

         REAL(DbKi)                   :: CB                                              ! Coefficient of seabed static friction (a negative value indicates no seabed) (-)
         REAL(DbKi)                   :: CBOvrEA                                         ! = CB/EA
         REAL(DbKi)                   :: DET                                             ! Determinant of the Jacobian matrix (m^2/N^2)
         REAL(DbKi)                   :: dHF                                             ! Increment in HF predicted by Newton-Raphson (N)
         REAL(DbKi)                   :: dVF                                             ! Increment in VF predicted by Newton-Raphson (N)
         REAL(DbKi)                   :: dXFdHF                                          ! Partial derivative of the calculated horizontal distance with respect to the horizontal fairlead tension (m/N): dXF(HF,VF)/dHF
         REAL(DbKi)                   :: dXFdVF                                          ! Partial derivative of the calculated horizontal distance with respect to the vertical   fairlead tension (m/N): dXF(HF,VF)/dVF
         REAL(DbKi)                   :: dZFdHF                                          ! Partial derivative of the calculated vertical   distance with respect to the horizontal fairlead tension (m/N): dZF(HF,VF)/dHF
         REAL(DbKi)                   :: dZFdVF                                          ! Partial derivative of the calculated vertical   distance with respect to the vertical   fairlead tension (m/N): dZF(HF,VF)/dVF
         REAL(DbKi)                   :: EA                                              ! Extensional stiffness of line (N)
         REAL(DbKi)                   :: EXF                                             ! Error function between calculated and known horizontal distance (meters): XF(HF,VF) - XF
         REAL(DbKi)                   :: EZF                                             ! Error function between calculated and known vertical   distance (meters): ZF(HF,VF) - ZF
         REAL(DbKi)                   :: HA                                              ! Effective horizontal tension in line at the anchor   (N)
         REAL(DbKi)                   :: HF                                              ! Effective horizontal tension in line at the fairlead (N)
         REAL(DbKi)                   :: HFOvrW                                          ! = HF/W
         REAL(DbKi)                   :: HFOvrWEA                                        ! = HF/WEA
         REAL(DbKi)                   :: L                                               ! Unstretched length of line (meters)
         REAL(DbKi)                   :: Lamda0                                          ! Catenary parameter used to generate the initial guesses of the horizontal and vertical tensions at the fairlead for the Newton-Raphson iteration (-)
         REAL(DbKi)                   :: LMax                                            ! Maximum stretched length of the line with seabed interaction beyond which the line would have to double-back on itself; here the line forms an "L" between the anchor and fairlead (i.e. it is horizontal along the seabed from the anchor, then vertical to the fairlead) (meters)
         REAL(DbKi)                   :: LMinVFOvrW                                      ! = L - VF/W
         REAL(DbKi)                   :: LOvrEA                                          ! = L/EA
         REAL(DbKi)                   :: s        (N)                                    ! Unstretched arc distance along line from anchor to each node where the line position and tension can be output (meters)
         REAL(DbKi)                   :: sOvrEA                                          ! = s(I)/EA
         REAL(DbKi)                   :: SQRT1VFOvrHF2                                   ! = SQRT( 1.0_DbKi + VFOvrHF2      )
         REAL(DbKi)                   :: SQRT1VFMinWLOvrHF2                              ! = SQRT( 1.0_DbKi + VFMinWLOvrHF2 )
         REAL(DbKi)                   :: SQRT1VFMinWLsOvrHF2                             ! = SQRT( 1.0_DbKi + VFMinWLsOvrHF*VFMinWLsOvrHF )
         REAL(DbKi)                   :: Te       (N)                                    ! Effective line tensions at each node (N)
         REAL(DbKi)                   :: Tol                                             ! Convergence tolerance within Newton-Raphson iteration specified as a fraction of tension (-)
         REAL(DbKi)                   :: VA                                              ! Effective vertical   tension in line at the anchor   (N)
         REAL(DbKi)                   :: VF                                              ! Effective vertical   tension in line at the fairlead (N)
         REAL(DbKi)                   :: VFMinWL                                         ! = VF - WL
         REAL(DbKi)                   :: VFMinWLOvrHF                                    ! = VFMinWL/HF
         REAL(DbKi)                   :: VFMinWLOvrHF2                                   ! = VFMinWLOvrHF*VFMinWLOvrHF
         REAL(DbKi)                   :: VFMinWLs                                        ! = VFMinWL + Ws
         REAL(DbKi)                   :: VFMinWLsOvrHF                                   ! = VFMinWLs/HF
         REAL(DbKi)                   :: VFOvrHF                                         ! = VF/HF
         REAL(DbKi)                   :: VFOvrHF2                                        ! = VFOvrHF*VFOvrHF
         REAL(DbKi)                   :: VFOvrWEA                                        ! = VF/WEA
         REAL(DbKi)                   :: W                                               ! Weight of line in fluid per unit length (N/m)
         REAL(DbKi)                   :: WEA                                             ! = W*EA
         REAL(DbKi)                   :: WL                                              ! Total weight of line in fluid (N): W*L
         REAL(DbKi)                   :: Ws                                              ! = W*s(I)
         REAL(DbKi)                   :: X        (N)                                    ! Horizontal locations of each line node relative to the anchor (meters)
         REAL(DbKi)                   :: XF                                              ! Horizontal distance between anchor and fairlead (meters)
         REAL(DbKi)                   :: XF2                                             ! = XF*XF
         REAL(DbKi)                   :: Z        (N)                                    ! Vertical   locations of each line node relative to the anchor (meters)
         REAL(DbKi)                   :: ZF                                              ! Vertical   distance between anchor and fairlead (meters)
         REAL(DbKi)                   :: ZF2                                             ! = ZF*ZF

         INTEGER(4)                   :: I                                               ! Index for counting iterations or looping through line nodes (-)
         INTEGER(4)                   :: MaxIter                                         ! Maximum number of Newton-Raphson iterations possible before giving up (-)

         LOGICAL                      :: FirstIter                                       ! Flag to determine whether or not this is the first time through the Newton-Raphson interation (flag)


         ErrStat = ERrId_None


            ! The Newton-Raphson iteration is only accurate in double precision, so
            !   convert the input arguments into double precision:

         CB     = REAL( CB_In    , DbKi )
         EA     = REAL( EA_In    , DbKi )
         HF = 0.0_DbKi  !    = REAL( HF_In    , DbKi )
         L      = REAL( L_In     , DbKi )
         s  (:) = REAL( s_In  (:), DbKi )
         Tol    = REAL( Tol_In   , DbKi )
         VF = 0.0_DbKi   ! keeping this for some error catching functionality? (at first glance)  ! VF     = REAL( VF_In    , DbKi )
         W      = REAL( W_In     , DbKi )
         XF     = REAL( XF_In    , DbKi )
         ZF     = REAL( ZF_In    , DbKi )


         
      !  HF and VF cannot be initialized to zero when a  portion of the line rests on the seabed and the anchor tension is nonzero
         
      ! Generate the initial guess values for the horizontal and vertical tensions
      !   at the fairlead in the Newton-Raphson iteration for the catenary mooring
      !   line solution.  Use starting values documented in: Peyrot, Alain H. and
      !   Goulois, A. M., "Analysis Of Cable Structures," Computers & Structures,
      !   Vol. 10, 1979, pp. 805-813:
         XF2     = XF*XF
         ZF2     = ZF*ZF

         IF     ( XF           == 0.0_DbKi    )  THEN ! .TRUE. if the current mooring line is exactly vertical
            Lamda0 = 1.0D+06
         ELSEIF ( L <= SQRT( XF2 + ZF2 ) )  THEN ! .TRUE. if the current mooring line is taut
            Lamda0 = 0.2_DbKi
         ELSE                                    ! The current mooring line must be slack and not vertical
            Lamda0 = SQRT( 3.0_DbKi*( ( L**2 - ZF2 )/XF2 - 1.0_DbKi ) )
         ENDIF

         HF = ABS( 0.5_DbKi*W*  XF/     Lamda0      )
         VF =      0.5_DbKi*W*( ZF/TANH(Lamda0) + L )         
                                    

            ! Abort when there is no solution or when the only possible solution is
            !   illogical:

         IF (    Tol <= EPSILON(TOL) )  THEN   ! .TRUE. when the convergence tolerance is specified incorrectly
           ErrStat = ErrID_Warn
           ErrMsg = ' Convergence tolerance must be greater than zero in routine Catenary().'
           return
         ELSEIF ( XF <  0.0_DbKi )  THEN   ! .TRUE. only when the local coordinate system is not computed correctly
           ErrStat = ErrID_Warn
           ErrMsg =  ' The horizontal distance between an anchor and its'// &
                         ' fairlead must not be less than zero in routine Catenary().'
           return

         ELSEIF ( ZF <  0.0_DbKi )  THEN   ! .TRUE. if the fairlead has passed below its anchor
           ErrStat = ErrID_Warn
           ErrMsg =  " A line's fairlead is defined as below its anchor. You may need to swap a line's fairlead and anchor end nodes."
           return

         ELSEIF ( L  <= 0.0_DbKi )  THEN   ! .TRUE. when the unstretched line length is specified incorrectly
           ErrStat = ErrID_Warn
           ErrMsg =  ' Unstretched length of line must be greater than zero in routine Catenary().'
           return

         ELSEIF ( EA <= 0.0_DbKi )  THEN   ! .TRUE. when the unstretched line length is specified incorrectly
           ErrStat = ErrID_Warn
           ErrMsg =  ' Extensional stiffness of line must be greater than zero in routine Catenary().'
           return

         ELSEIF ( W  == 0.0_DbKi )  THEN   ! .TRUE. when the weight of the line in fluid is zero so that catenary solution is ill-conditioned
           ErrStat = ErrID_Warn
           ErrMsg = ' The weight of the line in fluid must not be zero. '// &
                         ' Routine Catenary() cannot solve quasi-static mooring line solution.'
           return


         ELSEIF ( W  >  0.0_DbKi )  THEN   ! .TRUE. when the line will sink in fluid

            LMax      = XF - EA/W + SQRT( (EA/W)*(EA/W) + 2.0_DbKi*ZF*EA/W )  ! Compute the maximum stretched length of the line with seabed interaction beyond which the line would have to double-back on itself; here the line forms an "L" between the anchor and fairlead (i.e. it is horizontal along the seabed from the anchor, then vertical to the fairlead)

            IF ( ( L  >=  LMax   ) .AND. ( CB >= 0.0_DbKi ) )  then  ! .TRUE. if the line is as long or longer than its maximum possible value with seabed interaction
               ErrStat = ErrID_Warn
               !ErrMsg =  ' Unstretched mooring line length too large. '// &
               !             ' Routine Catenary() cannot solve quasi-static mooring line solution.'
               return
            END IF

         ENDIF


            ! Initialize some commonly used terms that don't depend on the iteration:

         WL      =          W  *L
         WEA     =          W  *EA
         LOvrEA  =          L  /EA
         CBOvrEA =          CB /EA
         MaxIter = INT(1.0_DbKi/Tol)   ! Smaller tolerances may take more iterations, so choose a maximum inversely proportional to the tolerance



            ! To avoid an ill-conditioned situation, ensure that the initial guess for
            !   HF is not less than or equal to zero.  Similarly, avoid the problems
            !   associated with having exactly vertical (so that HF is zero) or exactly
            !   horizontal (so that VF is zero) lines by setting the minimum values
            !   equal to the tolerance.  This prevents us from needing to implement
            !   the known limiting solutions for vertical or horizontal lines (and thus
            !   complicating this routine):

         HF = MAX( HF, Tol )
         XF = MAX( XF, Tol )
         ZF = MAX( ZF, TOl )



            ! Solve the analytical, static equilibrium equations for a catenary (or
            !   taut) mooring line with seabed interaction:

            ! Begin Newton-Raphson iteration:

         I         = 1        ! Initialize iteration counter
         FirstIter = .TRUE.   ! Initialize iteration flag

         DO


            ! Initialize some commonly used terms that depend on HF and VF:

            VFMinWL            = VF - WL
            LMinVFOvrW         = L  - VF/W
            HFOvrW             =      HF/W
            HFOvrWEA           =      HF/WEA
            VFOvrWEA           =      VF/WEA
            VFOvrHF            =      VF/HF
            VFMinWLOvrHF       = VFMinWL/HF
            VFOvrHF2           = VFOvrHF     *VFOvrHF
            VFMinWLOvrHF2      = VFMinWLOvrHF*VFMinWLOvrHF
            SQRT1VFOvrHF2      = SQRT( 1.0_DbKi + VFOvrHF2      )
            SQRT1VFMinWLOvrHF2 = SQRT( 1.0_DbKi + VFMinWLOvrHF2 )


            ! Compute the error functions (to be zeroed) and the Jacobian matrix
            !   (these depend on the anticipated configuration of the mooring line):

            IF ( ( CB <  0.0_DbKi ) .OR. ( W  <  0.0_DbKi ) .OR. ( VFMinWL >  0.0_DbKi ) )  THEN   ! .TRUE. when no portion of the line      rests on the seabed

               EXF    = (   LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                       &
                          - LOG( VFMinWLOvrHF +               SQRT1VFMinWLOvrHF2 )                                         )*HFOvrW &
                      + LOvrEA*  HF                         - XF
               EZF    = (                                     SQRT1VFOvrHF2                                              &
                          -                                   SQRT1VFMinWLOvrHF2                                           )*HFOvrW &
                      + LOvrEA*( VF - 0.5_DbKi*WL )         - ZF

               dXFdHF = (   LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                       &
                          - LOG( VFMinWLOvrHF +               SQRT1VFMinWLOvrHF2 )                                         )/     W &
                      - (      ( VFOvrHF      + VFOvrHF2     /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      ) &
                          -    ( VFMinWLOvrHF + VFMinWLOvrHF2/SQRT1VFMinWLOvrHF2 )/( VFMinWLOvrHF + SQRT1VFMinWLOvrHF2 )   )/     W &
                      + LOvrEA
               dXFdVF = (      ( 1.0_DbKi     + VFOvrHF      /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      ) &
                          -    ( 1.0_DbKi     + VFMinWLOvrHF /SQRT1VFMinWLOvrHF2 )/( VFMinWLOvrHF + SQRT1VFMinWLOvrHF2 )   )/     W
               dZFdHF = (                                     SQRT1VFOvrHF2                                              &
                          -                                   SQRT1VFMinWLOvrHF2                                           )/     W &
                      - (                       VFOvrHF2     /SQRT1VFOvrHF2                                              &
                          -                     VFMinWLOvrHF2/SQRT1VFMinWLOvrHF2                                           )/     W
               dZFdVF = (                       VFOvrHF      /SQRT1VFOvrHF2                                              &
                          -                     VFMinWLOvrHF /SQRT1VFMinWLOvrHF2                                           )/     W &
                      + LOvrEA


            ELSEIF (                                           -CB*VFMinWL <  HF         )  THEN   ! .TRUE. when a  portion of the line      rests on the seabed and the anchor tension is nonzero

               EXF    =     LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          *HFOvrW &
                      - 0.5_DbKi*CBOvrEA*W*  LMinVFOvrW*LMinVFOvrW                                                                  &
                      + LOvrEA*  HF           + LMinVFOvrW  - XF
               EZF    = (                                     SQRT1VFOvrHF2                                   - 1.0_DbKi   )*HFOvrW &
                      + 0.5_DbKi*VF*VFOvrWEA                - ZF

               dXFdHF =     LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          /     W &
                      - (      ( VFOvrHF      + VFOvrHF2     /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2        ) )/     W &
                      + LOvrEA
               dXFdVF = (      ( 1.0_DbKi     + VFOvrHF      /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2        ) )/     W &
                      + CBOvrEA*LMinVFOvrW - 1.0_DbKi/W
               dZFdHF = (                                     SQRT1VFOvrHF2                                   - 1.0_DbKi &
                          -                     VFOvrHF2     /SQRT1VFOvrHF2                                                )/     W
               dZFdVF = (                       VFOvrHF      /SQRT1VFOvrHF2                                                )/     W &
                      + VFOvrWEA


            ELSE                                                ! 0.0_DbKi <  HF  <= -CB*VFMinWL   !             A  portion of the line must rest  on the seabed and the anchor tension is    zero

               EXF    =     LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          *HFOvrW &
                      - 0.5_DbKi*CBOvrEA*W*( LMinVFOvrW*LMinVFOvrW - ( LMinVFOvrW - HFOvrW/CB )*( LMinVFOvrW - HFOvrW/CB ) )        &
                      + LOvrEA*  HF           + LMinVFOvrW  - XF
               EZF    = (                                     SQRT1VFOvrHF2                                   - 1.0_DbKi   )*HFOvrW &
                      + 0.5_DbKi*VF*VFOvrWEA                - ZF

               dXFdHF =     LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          /     W &
                      - (      ( VFOvrHF      + VFOvrHF2     /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      )   )/     W &
                      + LOvrEA - ( LMinVFOvrW - HFOvrW/CB )/EA
               dXFdVF = (      ( 1.0_DbKi     + VFOvrHF      /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      )   )/     W &
                      + HFOvrWEA           - 1.0_DbKi/W
               dZFdHF = (                                     SQRT1VFOvrHF2                                   - 1.0_DbKi &
                          -                     VFOvrHF2     /SQRT1VFOvrHF2                                                )/     W
               dZFdVF = (                       VFOvrHF      /SQRT1VFOvrHF2                                                )/     W &
                      + VFOvrWEA


            ENDIF


            ! Compute the determinant of the Jacobian matrix and the incremental
            !   tensions predicted by Newton-Raphson:
            
            
            DET = dXFdHF*dZFdVF - dXFdVF*dZFdHF
            
            if ( EqualRealNos( DET, 0.0_DbKi ) ) then               
!bjj: there is a serious problem with the debugger here when DET = 0
                ErrStat = ErrID_Warn
                ErrMsg =  ' Iteration not convergent (DET is 0). '// &
                          ' Routine Catenary() cannot solve quasi-static mooring line solution.'
                return
            endif

               
            dHF = ( -dZFdVF*EXF + dXFdVF*EZF )/DET    ! This is the incremental change in horizontal tension at the fairlead as predicted by Newton-Raphson
            dVF = (  dZFdHF*EXF - dXFdHF*EZF )/DET    ! This is the incremental change in vertical   tension at the fairlead as predicted by Newton-Raphson

            dHF = dHF*( 1.0_DbKi - Tol*I )            ! Reduce dHF by factor (between 1 at I = 1 and 0 at I = MaxIter) that reduces linearly with iteration count to ensure that we converge on a solution even in the case were we obtain a nonconvergent cycle about the correct solution (this happens, for example, if we jump to quickly between a taut and slack catenary)
            dVF = dVF*( 1.0_DbKi - Tol*I )            ! Reduce dHF by factor (between 1 at I = 1 and 0 at I = MaxIter) that reduces linearly with iteration count to ensure that we converge on a solution even in the case were we obtain a nonconvergent cycle about the correct solution (this happens, for example, if we jump to quickly between a taut and slack catenary)

            dHF = MAX( dHF, ( Tol - 1.0_DbKi )*HF )   ! To avoid an ill-conditioned situation, make sure HF does not go less than or equal to zero by having a lower limit of Tol*HF [NOTE: the value of dHF = ( Tol - 1.0_DbKi )*HF comes from: HF = HF + dHF = Tol*HF when dHF = ( Tol - 1.0_DbKi )*HF]

            ! Check if we have converged on a solution, or restart the iteration, or
            !   Abort if we cannot find a solution:

            IF ( ( ABS(dHF) <= ABS(Tol*HF) ) .AND. ( ABS(dVF) <= ABS(Tol*VF) ) )  THEN ! .TRUE. if we have converged; stop iterating! [The converge tolerance, Tol, is a fraction of tension]

               EXIT


            ELSEIF ( ( I == MaxIter )        .AND. (       FirstIter         ) )  THEN ! .TRUE. if we've iterated MaxIter-times for the first time;

            ! Perhaps we failed to converge because our initial guess was too far off.
            !   (This could happen, for example, while linearizing a model via large
            !   pertubations in the DOFs.)  Instead, use starting values documented in:
            !   Peyrot, Alain H. and Goulois, A. M., "Analysis Of Cable Structures,"
            !   Computers & Structures, Vol. 10, 1979, pp. 805-813:
            ! NOTE: We don't need to check if the current mooring line is exactly
            !       vertical (i.e., we don't need to check if XF == 0.0), because XF is
            !       limited by the tolerance above.

               XF2 = XF*XF
               ZF2 = ZF*ZF

               IF ( L <= SQRT( XF2 + ZF2 ) )  THEN ! .TRUE. if the current mooring line is taut
                  Lamda0 = 0.2_DbKi
               ELSE                                ! The current mooring line must be slack and not vertical
                  Lamda0 = SQRT( 3.0_DbKi*( ( L*L - ZF2 )/XF2 - 1.0_DbKi ) )
               ENDIF

               HF  = MAX( ABS( 0.5_DbKi*W*  XF/     Lamda0      ), Tol )   ! As above, set the lower limit of the guess value of HF to the tolerance
               VF  =           0.5_DbKi*W*( ZF/TANH(Lamda0) + L )


            ! Restart Newton-Raphson iteration:

               I         = 0
               FirstIter = .FALSE.
               dHF       = 0.0_DbKi
               dVF       = 0.0_DbKi


           ELSEIF ( ( I == MaxIter )        .AND. ( .NOT. FirstIter         ) )  THEN ! .TRUE. if we've iterated as much as we can take without finding a solution; Abort
             ErrStat = ErrID_Warn
             ErrMsg =  ' Iteration not convergent. '// &
                       ' Routine Catenary() cannot solve quasi-static mooring line solution.'
             RETURN


           ENDIF


            ! Increment fairlead tensions and iteration counter so we can try again:

            HF = HF + dHF
            VF = VF + dVF

            I  = I  + 1


         ENDDO



            ! We have found a solution for the tensions at the fairlead!

            ! Now compute the tensions at the anchor and the line position and tension
            !   at each node (again, these depend on the configuration of the mooring
            !   line):

         IF ( ( CB <  0.0_DbKi ) .OR. ( W  <  0.0_DbKi ) .OR. ( VFMinWL >  0.0_DbKi ) )  THEN   ! .TRUE. when no portion of the line      rests on the seabed

            ! Anchor tensions:

            HA = HF
            VA = VFMinWL


            ! Line position and tension at each node:

            DO I = 1,N  ! Loop through all nodes where the line position and tension are to be computed

               IF ( ( s(I) <  0.0_DbKi ) .OR. ( s(I) >  L ) )  THEN
                 ErrStat = ErrID_Warn
                 ErrMsg = ' All line nodes must be located between the anchor ' &
                                 //'and fairlead (inclusive) in routine Catenary().'
                 RETURN
               END IF

               Ws                  = W       *s(I)                                  ! Initialize
               VFMinWLs            = VFMinWL + Ws                                   ! some commonly
               VFMinWLsOvrHF       = VFMinWLs/HF                                    ! used terms
               sOvrEA              = s(I)    /EA                                    ! that depend
               SQRT1VFMinWLsOvrHF2 = SQRT( 1.0_DbKi + VFMinWLsOvrHF*VFMinWLsOvrHF ) ! on s(I)

               X (I)    = (   LOG( VFMinWLsOvrHF + SQRT1VFMinWLsOvrHF2 ) &
                            - LOG( VFMinWLOvrHF  + SQRT1VFMinWLOvrHF2  )   )*HFOvrW                     &
                        + sOvrEA*  HF
               Z (I)    = (                        SQRT1VFMinWLsOvrHF2   &
                            -                      SQRT1VFMinWLOvrHF2      )*HFOvrW                     &
                        + sOvrEA*(         VFMinWL + 0.5_DbKi*Ws    )
               Te(I)    = SQRT( HF*HF +    VFMinWLs*VFMinWLs )

            ENDDO       ! I - All nodes where the line position and tension are to be computed


         ELSEIF (                                           -CB*VFMinWL <  HF         )  THEN   ! .TRUE. when a  portion of the line      rests on the seabed and the anchor tension is nonzero

            ! Anchor tensions:

            HA = HF + CB*VFMinWL
            VA = 0.0_DbKi


            ! Line position and tension at each node:

            DO I = 1,N  ! Loop through all nodes where the line position and tension are to be computed

               IF ( ( s(I) <  0.0_DbKi ) .OR. ( s(I) >  L ) )  THEN
                 ErrStat = ErrID_Warn
                 ErrMsg =  ' All line nodes must be located between the anchor ' &
                                 //'and fairlead (inclusive) in routine Catenary().'
                 RETURN
               END IF

               Ws                  = W       *s(I)                                  ! Initialize
               VFMinWLs            = VFMinWL + Ws                                   ! some commonly
               VFMinWLsOvrHF       = VFMinWLs/HF                                    ! used terms
               sOvrEA              = s(I)    /EA                                    ! that depend
               SQRT1VFMinWLsOvrHF2 = SQRT( 1.0_DbKi + VFMinWLsOvrHF*VFMinWLsOvrHF ) ! on s(I)

               IF (     s(I) <= LMinVFOvrW             )  THEN ! .TRUE. if this node rests on the seabed and the tension is nonzero

                  X (I) = s(I)                                                                          &
                        + sOvrEA*( HF + CB*VFMinWL + 0.5_DbKi*Ws*CB )
                  Z (I) = 0.0_DbKi
                  Te(I) =       HF    + CB*VFMinWLs

               ELSE           ! LMinVFOvrW < s <= L            !           This node must be above the seabed

                  X (I) =     LOG( VFMinWLsOvrHF + SQRT1VFMinWLsOvrHF2 )    *HFOvrW                     &
                        + sOvrEA*  HF + LMinVFOvrW                    - 0.5_DbKi*CB*VFMinWL*VFMinWL/WEA
                  Z (I) = ( - 1.0_DbKi           + SQRT1VFMinWLsOvrHF2     )*HFOvrW                     &
                        + sOvrEA*(         VFMinWL + 0.5_DbKi*Ws    ) + 0.5_DbKi*   VFMinWL*VFMinWL/WEA
                  Te(I) = SQRT( HF*HF +    VFMinWLs*VFMinWLs )

               ENDIF

            ENDDO       ! I - All nodes where the line position and tension are to be computed


         ELSE                                                ! 0.0_DbKi <  HF  <= -CB*VFMinWL   !             A  portion of the line must rest  on the seabed and the anchor tension is    zero

            ! Anchor tensions:

            HA = 0.0_DbKi
            VA = 0.0_DbKi


            ! Line position and tension at each node:

            DO I = 1,N  ! Loop through all nodes where the line position and tension are to be computed

               IF ( ( s(I) <  0.0_DbKi ) .OR. ( s(I) >  L ) )  THEN
                  ErrStat = ErrID_Warn
                   ErrMsg =  ' All line nodes must be located between the anchor ' &
                                 //'and fairlead (inclusive) in routine Catenary().'
                   RETURN
               END IF

               Ws                  = W       *s(I)                                  ! Initialize
               VFMinWLs            = VFMinWL + Ws                                   ! some commonly
               VFMinWLsOvrHF       = VFMinWLs/HF                                    ! used terms
               sOvrEA              = s(I)    /EA                                    ! that depend
               SQRT1VFMinWLsOvrHF2 = SQRT( 1.0_DbKi + VFMinWLsOvrHF*VFMinWLsOvrHF ) ! on s(I)

               IF (     s(I) <= LMinVFOvrW - HFOvrW/CB )  THEN ! .TRUE. if this node rests on the seabed and the tension is    zero

                  X (I) = s(I)
                  Z (I) = 0.0_DbKi
                  Te(I) = 0.0_DbKi

               ELSEIF ( s(I) <= LMinVFOvrW             )  THEN ! .TRUE. if this node rests on the seabed and the tension is nonzero

                  X (I) = s(I)                     - ( LMinVFOvrW - 0.5_DbKi*HFOvrW/CB )*HF/EA          &
                        + sOvrEA*( HF + CB*VFMinWL + 0.5_DbKi*Ws*CB ) + 0.5_DbKi*CB*VFMinWL*VFMinWL/WEA
                  Z (I) = 0.0_DbKi
                  Te(I) =       HF    + CB*VFMinWLs

               ELSE           ! LMinVFOvrW < s <= L            !           This node must be above the seabed

                  X (I) =     LOG( VFMinWLsOvrHF + SQRT1VFMinWLsOvrHF2 )    *HFOvrW                     &
                        + sOvrEA*  HF + LMinVFOvrW - ( LMinVFOvrW - 0.5_DbKi*HFOvrW/CB )*HF/EA
                  Z (I) = ( - 1.0_DbKi           + SQRT1VFMinWLsOvrHF2     )*HFOvrW                     &
                        + sOvrEA*(         VFMinWL + 0.5_DbKi*Ws    ) + 0.5_DbKi*   VFMinWL*VFMinWL/WEA
                  Te(I) = SQRT( HF*HF +    VFMinWLs*VFMinWLs )

               ENDIF

            ENDDO       ! I - All nodes where the line position and tension are to be computed


         ENDIF



            ! The Newton-Raphson iteration is only accurate in double precision, so
            !   convert the output arguments back into the default precision for real
            !   numbers:

         !HA_In    = REAL( HA   , DbKi )  !mth: for this I only care about returning node positions
         !HF_In    = REAL( HF   , DbKi )
         !Te_In(:) = REAL( Te(:), DbKi )
         !VA_In    = REAL( VA   , DbKi )
         !VF_In    = REAL( VF   , DbKi )
         X_In (:) = REAL( X (:), DbKi )
         Z_In (:) = REAL( Z (:), DbKi )

      END SUBROUTINE Catenary
      !-----------------------------------------------------------------------


   END SUBROUTINE Line_Initialize
   !--------------------------------------------------------------

   
   !--------------------------------------------------------------
   SUBROUTINE Line_SetState(Line, X, t)

      TYPE(MD_Line),    INTENT(INOUT)  :: Line           ! the current Line object
      Real(DbKi),       INTENT(IN   )  :: X(:)           ! state vector section for this line
      Real(DbKi),       INTENT(IN   )  :: t              ! instantaneous time

      INTEGER(IntKi)                   :: i              ! index of segments or nodes along line
      INTEGER(IntKi)                   :: J              ! index
   

      ! store current time
      Line%time = t
      
      ! set interior node positions and velocities based on state vector
      DO I=1,Line%N-1
         DO J=1,3
         
            Line%r( J,I) = X( 3*Line%N-3 + 3*I-3 + J)  ! get positions
            Line%rd(J,I) = X(              3*I-3 + J)  ! get velocities
            
         END DO
      END DO
      
      ! if using viscoelastic model, also set the static stiffness stretch
      if (Line%ElasticMod == 2) then
         do I=1,Line%N
            Line%dl_1(I) = X( 6*Line%N-6 + I)   ! these will be the last N entries in the state vector
         end do
      end if
         
   END SUBROUTINE Line_SetState
   !--------------------------------------------------------------

   !--------------------------------------------------------------
   SUBROUTINE Line_GetStateDeriv(Line, Xd, m, p)  !, FairFtot, FairMtot, AnchFtot, AnchMtot)

      TYPE(MD_Line),          INTENT(INOUT) :: Line    ! the current Line object
      Real(DbKi),             INTENT(INOUT) :: Xd(:)   ! state derivative vector section for this line
      TYPE(MD_MiscVarType),   INTENT(INOUT) :: m       ! passing along all mooring objects
      TYPE(MD_ParameterType), INTENT(IN   ) :: p       ! Parameters
      
      !   Real(DbKi), INTENT( IN )      :: X(:)           ! state vector, provided
      !   Real(DbKi), INTENT( INOUT )   :: Xd(:)          ! derivative of state vector, returned ! cahnged to INOUT
      !   Real(DbKi), INTENT (IN)       :: t              ! instantaneous time
      !   TYPE(MD_Line), INTENT (INOUT) :: Line           ! label for the current line, for convenience
      !   TYPE(MD_LineProp), INTENT(IN) :: LineProp       ! the single line property set for the line of interest
      !    Real(DbKi), INTENT(INOUT)     :: FairFtot(:)    ! total force on Connect top of line is attached to
      !    Real(DbKi), INTENT(INOUT)     :: FairMtot(:,:)  ! total mass of Connect top of line is attached to
      !    Real(DbKi), INTENT(INOUT)     :: AnchFtot(:)    ! total force on Connect bottom of line is attached to
      !    Real(DbKi), INTENT(INOUT)     :: AnchMtot(:,:)  ! total mass of Connect bottom of line is attached to


      INTEGER(IntKi)                   :: i              ! index of segments or nodes along line
      INTEGER(IntKi)                   :: J              ! index
      INTEGER(IntKi)                   :: K              ! index
      INTEGER(IntKi)                   :: N              ! number of segments in line
      Real(DbKi)                       :: d              ! line diameter
      Real(DbKi)                       :: rho            ! line material density [kg/m^3]
      Real(DbKi)                       :: Sum1           ! for summing squares
      Real(DbKi)                       :: dummyLength    ! 
      Real(DbKi)                       :: m_i            ! node mass
      Real(DbKi)                       :: v_i            ! node submerged volume
      Real(DbKi)                       :: Vi(3)          ! relative water velocity at a given node
      Real(DbKi)                       :: Vp(3)          ! transverse relative water velocity component at a given node
      Real(DbKi)                       :: Vq(3)          ! tangential relative water velocity component at a given node
      Real(DbKi)                       :: SumSqVp        !
      Real(DbKi)                       :: SumSqVq        !
      Real(DbKi)                       :: MagVp          !
      Real(DbKi)                       :: MagVq          !
      Real(DbKi)                       :: MagT           ! tension stiffness force magnitude
      Real(DbKi)                       :: MagTd          ! tension damping force magnitude
      Real(DbKi)                       :: Xi             ! used in interpolating from lookup table
      Real(DbKi)                       :: Yi             ! used in interpolating from lookup table
      Real(DbKi)                       :: dl             ! stretch of a segment [m]
      Real(DbKi)                       :: ld_1           ! rate of change of static stiffness portion of segment [m/s]
      Real(DbKi)                       :: EA_1           ! stiffness of 'static stiffness' portion of segment, combines with dynamic stiffness to give static stiffnes [m/s]

      Real(DbKi)                       :: Kurvi          ! temporary curvature value [1/m]
      Real(DbKi)                       :: pvec(3)        ! the p vector used in bending stiffness calcs
      Real(DbKi)                       :: Mforce_im1(3)  ! force vector for a contributor to the effect of a bending moment [N]
      Real(DbKi)                       :: Mforce_ip1(3)  ! force vector for a contributor to the effect of a bending moment [N]
      Real(DbKi)                       :: Mforce_i(  3)  ! force vector for a contributor to the effect of a bending moment [N]

      Real(DbKi)                       :: depth          ! local water depth interpolated from bathymetry grid [m]
      Real(DbKi)                       :: nvec(3)        ! local seabed surface normal vector (positive out)
      Real(DbKi)                       :: Fn(3)          ! seabed contact normal force vector
      Real(DbKi)                       :: Vn(3)          ! normal velocity of a line node relative to the seabed slope [m/s]
      Real(DbKi)                       :: Vsb(3)         ! tangent velocity of a line node relative to the seabed slope [m/s]
      Real(DbKi)                       :: Va(3)          ! velocity of a line node in the axial or "in-line" direction [m/s]
      Real(DbKi)                       :: Vt(3)          ! velocity of a line node in the transverse direction [m/s]
      Real(DbKi)                       :: VtMag          ! magnitude of the transverse velocity of a line node [m/s]
      Real(DbKi)                       :: VaMag          ! magnitude of the axial velocity of a line node [m/s]
      Real(DbKi)                       :: FkTmax         ! maximum kinetic friction force in the transverse direction (scalar)
      Real(DbKi)                       :: FkAmax         ! maximum kinetic friction force in the axial direction (scalar)
      Real(DbKi)                       :: FkT(3)         ! kinetic friction force in the transverse direction (vector)
      Real(DbKi)                       :: FkA(3)         ! kinetic friction force in the axial direction (vector)
      !Real(DbKi)                       :: mc_T           ! ratio of the transverse static friction coefficient to the transverse kinetic friction coefficient
      !Real(DbKi)                       :: mc_A           ! ratio of the axial static friction coefficient to the axial kinetic friction coefficient
      Real(DbKi)                       :: FfT(3)         ! total friction force in the transverse direction
      Real(DbKi)                       :: FfA(3)         ! total friction force in the axial direction
      Real(DbKi)                       :: Ff(3)          ! total friction force on the line node


      N = Line%N                      ! for convenience
      d = Line%d    
      rho = Line%rho

      ! note that end node kinematics should have already been set by attached objects

      !   ! set end node positions and velocities from connect objects' states
      !   DO J = 1, 3
      !      Line%r( J,N) = m%ConnectList(Line%FairConnect)%r(J)
      !      Line%r( J,0) = m%ConnectList(Line%AnchConnect)%r(J)
      !      Line%rd(J,N) = m%ConnectList(Line%FairConnect)%rd(J)
      !      Line%rd(J,0) = m%ConnectList(Line%AnchConnect)%rd(J)
      !   END DO



      ! -------------------- calculate various kinematic quantities ---------------------------
      DO I = 1, N
         
         
         ! calculate current (Stretched) segment lengths and unit tangent vectors (qs) for each segment (this is used for bending calculations)
         CALL UnitVector(Line%r(:,I-1), Line%r(:,I), Line%qs(:,I), Line%lstr(I))
         
         ! should add catch here for if lstr is ever zero
   
         Sum1 = 0.0_DbKi
         DO J = 1, 3
            Sum1 = Sum1 + (Line%r(J,I) - Line%r(J,I-1))*(Line%rd(J,I) - Line%rd(J,I-1))
         END DO
         Line%lstrd(I) = Sum1/Line%lstr(I)                          ! segment stretched length rate of change

      !       Line%V(I) = Pi/4.0 * d*d*Line%l(I)                        !volume attributed to segment
      END DO

      !calculate unit tangent vectors (q) for each internal node based on adjacent node positions
      DO I = 1, N-1
        CALL UnitVector(Line%r(:,I-1), Line%r(:,I+1), Line%q(:,I), dummyLength)
      END DO
      
      ! calculate unit tangent vectors for either end node if the line has no bending stiffness of if either end is pinned (otherwise it's already been set via setEndStateFromRod)
      if ((Line%endTypeA==0) .or. (Line%EI==0.0)) then 
         CALL UnitVector(Line%r(:,0), Line%r(:,1), Line%q(:,0), dummyLength)
      end if
      if ((Line%endTypeB==0) .or. (Line%EI==0.0)) then 
         CALL UnitVector(Line%r(:,N-1), Line%r(:,N), Line%q(:,N), dummyLength)
      end if
      
      ! apply wave kinematics (if there are any) 
      DO i=0,N
         CALL getWaterKin(p, Line%r(1,i), Line%r(2,i), Line%r(3,i), Line%time, m%WaveTi, Line%U(:,i), Line%Ud(:,i), Line%zeta(i), Line%PDyn(i))
      END DO
      

      ! --------------- calculate mass (including added mass) matrix for each node -----------------
      DO I = 0, N
         IF (I==0) THEN
            m_i = Pi/8.0 *d*d*Line%l(1)*rho
            v_i = 0.5 *Line%V(1)
         ELSE IF (I==N) THEN
            m_i = pi/8.0 *d*d*Line%l(N)*rho;
            v_i = 0.5*Line%V(N)
         ELSE
            m_i = pi/8.0 * d*d*rho*(Line%l(I) + Line%l(I+1))
            v_i = 0.5 *(Line%V(I) + Line%V(I+1))
         END IF

         DO J=1,3
            DO K=1,3
               IF (J==K) THEN
                  Line%M(K,J,I) = m_i + p%rhoW*v_i*( Line%Can*(1 - Line%q(J,I)*Line%q(K,I)) + Line%Cat*Line%q(J,I)*Line%q(K,I) )
               ELSE
                  Line%M(K,J,I) = p%rhoW*v_i*( Line%Can*(-Line%q(J,I)*Line%q(K,I)) + Line%Cat*Line%q(J,I)*Line%q(K,I) )
               END IF
            END DO
         END DO

         CALL Inverse3by3(Line%S(:,:,I), Line%M(:,:,I))             ! invert mass matrix
      END DO


      ! ------------------  CALCULATE FORCES ON EACH NODE ----------------------------

      ! loop through the segments
      DO I = 1, N
      
         ! handle nonlinear stiffness if needed
         if (Line%nEApoints > 0) then
            
            Xi = Line%lstr(I)/Line%l(I) - 1.0                   ! strain rate based on inputs
            Yi = 0.0_DbKi

            ! find stress based on strain
            if (Xi < 0.0) then                                  ! if negative strain (compression), zero stress
               Yi = 0.0_DbKi
            else if (Xi < Line%stiffXs(1)) then                 ! if strain below first data point, interpolate from zero
               Yi = Xi * Line%stiffYs(1)/Line%stiffXs(1)
            else if (Xi >= Line%stiffXs(Line%nEApoints)) then   ! if strain exceeds last data point, use last data point
               Yi = Line%stiffYs(Line%nEApoints)
            else                                                ! otherwise we're in range of the table so interpolate!
               do J=1, Line%nEApoints-1                         ! go through lookup table until next entry exceeds inputted strain rate
                  if (Line%stiffXs(J+1) > Xi) then
                     Yi = Line%stiffYs(J) + (Xi-Line%stiffXs(J)) * (Line%stiffYs(J+1)-Line%stiffYs(J))/(Line%stiffXs(J+1)-Line%stiffXs(J))
                     exit
                  end if
               end do
            end if

            ! calculate a young's modulus equivalent value based on stress/strain
            Line%EA = Yi/Xi
         end if
         
         
         ! >>>> could do similar as above for nonlinear damping or bending stiffness <<<<         
         if (Line%nBApoints > 0) print *, 'Nonlinear elastic damping not yet implemented'
         if (Line%nEIpoints > 0) print *, 'Nonlinear bending stiffness not yet implemented'
            
            
         ! basic elasticity model
         if (Line%ElasticMod == 1) then
            ! line tension, inherently including possibility of dynamic length changes in l term
            if (Line%lstr(I)/Line%l(I) > 1.0) then
               MagT  =  Line%EA *( Line%lstr(I)/Line%l(I) - 1.0 )
            else
               MagT = 0.0_DbKi                              ! cable can't "push"
            end if
            ! line internal damping force based on line-specific BA value, including possibility of dynamic length changes in l and ld terms
            MagTd = Line%BA* ( Line%lstrd(I) -  Line%lstr(I)*Line%ld(I)/Line%l(I) )/Line%l(I)
         
         ! viscoelastic model
         else if (Line%ElasticMod == 2) then
         
            EA_1 = Line%EA_D*Line%EA/(Line%EA_D - Line%EA)! calculated EA_1 which is the stiffness in series with EA_D that will result in the desired static stiffness of EA_S
         
            dl = Line%lstr(I) - Line%l(I) ! delta l of this segment
         
            ld_1 = (Line%EA_D*dl - (Line%EA_D + EA_1)*Line%dl_1(I) + Line%BA_D*Line%lstrd(I)) /( Line%BA_D + Line%BA) ! rate of change of static stiffness portion [m/s]
            
            !MagT = (Line%EA*Line%dl_S(I) + Line%BA*ld_S)/ Line%l(I)   ! compute tension based on static portion (dynamic portion would give same)
            MagT  = EA_1*Line%dl_1(I)/ Line%l(I)  
            MagTd = Line%BA*ld_1        / Line%l(I)  
            
            ! update state derivative for static stiffness stretch (last N entries in the state vector)
            Xd( 6*N-6 + I) = ld_1
         
         end if

         
         do J = 1, 3
            !Line%T(J,I) = Line%EA *( 1.0/Line%l(I) - 1.0/Line%lstr(I) ) * (Line%r(J,I)-Line%r(J,I-1))
            Line%T(J,I)  = MagT * Line%qs(J,I)
            !Line%Td(J,I) = Line%BA* ( Line%lstrd(I) / Line%l(I) ) * (Line%r(J,I)-Line%r(J,I-1)) / Line%lstr(I)  ! note new form of damping coefficient, BA rather than Cint
            Line%Td(J,I) = MagTd * Line%qs(J,I)
         end do
      end do


      ! Bending loads   
      Line%Bs = 0.0_DbKi         ! zero bending forces
         
      if (Line%EI > 0) then
         ! loop through all nodes to calculate bending forces due to bending stiffness
         do i=0,N
         
            ! end node A case (only if attached to a Rod, i.e. a cantilever rather than pinned connection)
            if (i==0) then
            
               if (Line%endTypeA > 0) then ! if attached to Rod i.e. cantilever connection
               
                  Kurvi = GetCurvature(Line%lstr(1), Line%q(:,0), Line%qs(:,1))  ! curvature (assuming rod angle is node angle which is middle of if there was a segment -1/2)
         
                  pvec = cross_product(Line%q(:,0), Line%qs(:,1))                ! get direction of bending radius axis
                  
                  Mforce_ip1 = cross_product(Line%qs(:,1), pvec)                 ! get direction of resulting force from bending to apply on node i+1

                  call scalevector(pvec, Kurvi*Line%EI, Line%endMomentA)         ! record bending moment at end for potential application to attached object

                  call scalevector(Mforce_ip1, Kurvi*Line%EI/Line%lstr(1), Mforce_ip1)  ! scale force direction vectors by desired moment force magnitudes to get resulting forces on adjacent nodes
                     
                  Mforce_i = -Mforce_ip1                                         ! set force on node i to cancel out forces on adjacent nodes
                  
                  ! apply these forces to the node forces
                  Line%Bs(:,i  ) = Line%Bs(:,i  ) +  Mforce_i
                  Line%Bs(:,i+1) = Line%Bs(:,i+1) +  Mforce_ip1
                  
               end if
            
            ! end node A case (only if attached to a Rod, i.e. a cantilever rather than pinned connection)
            else if (i==N) then
            
               if (Line%endTypeB > 0) then ! if attached to Rod i.e. cantilever connection
               
                  Kurvi = GetCurvature(Line%lstr(N), Line%qs(:,N), Line%q(:,N))  ! curvature (assuming rod angle is node angle which is middle of if there was a segment -1/2
                  
                  pvec = cross_product(Line%qs(:,N), Line%q(:,N))                ! get direction of bending radius axis
                  
                  Mforce_im1 = cross_product(Line%qs(:,N), pvec)                 ! get direction of resulting force from bending to apply on node i-1
                  
                  call scalevector(pvec, -Kurvi*Line%EI, Line%endMomentB )       ! record bending moment at end (note end B is oposite sign as end A)
                  
                  call scalevector(Mforce_im1, Kurvi*Line%EI/Line%lstr(N), Mforce_im1)  ! scale force direction vectors by desired moment force magnitudes to get resulting forces on adjacent nodes
                  
                  Mforce_i =  -Mforce_im1                                         ! set force on node i to cancel out forces on adjacent nodes
                  
                  ! apply these forces to the node forces
                  Line%Bs(:,i-1) = Line%Bs(:,i-1) +  Mforce_im1
                  Line%Bs(:,i  ) = Line%Bs(:,i  ) +  Mforce_i
                  
               end if
            
            else   ! internal node
            
               Kurvi = GetCurvature(Line%lstr(i)+Line%lstr(i+1), Line%qs(:,i), Line%qs(:,i+1))  ! curvature
               
               pvec = cross_product(Line%qs(:,i), Line%qs(:,i+1))  ! get direction of bending radius axis
               
               Mforce_im1 = cross_product(Line%qs(:,i  ), pvec)    ! get direction of resulting force from bending to apply on node i-1
               Mforce_ip1 = cross_product(Line%qs(:,i+1), pvec)    ! get direction of resulting force from bending to apply on node i+1
               
               ! scale force direction vectors by desired moment force magnitudes to get resulting forces on adjacent nodes
               call scalevector(Mforce_im1, Kurvi*Line%EI/Line%lstr(i  ), Mforce_im1)
               call scalevector(Mforce_ip1, Kurvi*Line%EI/Line%lstr(i+1), Mforce_ip1)
        
               Mforce_i = -Mforce_im1 - Mforce_ip1                 ! set force on node i to cancel out forces on adjacent nodes
               
               ! apply these forces to the node forces
               Line%Bs(:,i-1) = Line%Bs(:,i-1) + Mforce_im1 
               Line%Bs(:,i  ) = Line%Bs(:,i  ) + Mforce_i
               Line%Bs(:,i+1) = Line%Bs(:,i+1) + Mforce_ip1
               
            end if
                       
            ! record curvature at node
            Line%Kurv(i) = Kurvi
            
         end do   ! for i=0,N (looping through nodes)
      end if  ! if EI > 0
         
         


      ! loop through the nodes
      DO I = 0, N

         !submerged weight (including buoyancy)
         IF (I==0) THEN
            Line%W(3,I) = Pi/8.0*d*d* Line%l(1)*(rho - p%rhoW) *(-p%g)   ! assuming g is positive
         ELSE IF (i==N)  THEN
            Line%W(3,I) = pi/8.0*d*d* Line%l(N)*(rho - p%rhoW) *(-p%g)
         ELSE
            Line%W(3,I) = pi/8.0*d*d* (Line%l(I)*(rho - p%rhoW) + Line%l(I+1)*(rho - p%rhoW) )*(-p%g)  ! left in this form for future free surface handling
         END IF

         ! relative flow velocities
         DO J = 1, 3
            Vi(J) = Line%U(J,I) - Line%rd(J,I)                               ! relative flow velocity over node -- this is where wave velicites would be added
         END DO

         ! decomponse relative flow into components
         SumSqVp = 0.0_DbKi                                         ! start sums of squares at zero
         SumSqVq = 0.0_DbKi
         DO J = 1, 3
            Vq(J) = DOT_PRODUCT( Vi , Line%q(:,I) ) * Line%q(J,I);   ! tangential relative flow component
            Vp(J) = Vi(J) - Vq(J)                                    ! transverse relative flow component
            SumSqVq = SumSqVq + Vq(J)*Vq(J)
            SumSqVp = SumSqVp + Vp(J)*Vp(J)
         END DO
         MagVp = sqrt(SumSqVp)                                      ! get magnitudes of flow components
         MagVq = sqrt(SumSqVq)

         ! transverse and tangenential drag
         IF (I==0) THEN
            Line%Dp(:,I) = 0.25*p%rhoW*Line%Cdn*    d*Line%l(1) * MagVp * Vp
            Line%Dq(:,I) = 0.25*p%rhoW*Line%Cdt* Pi*d*Line%l(1) * MagVq * Vq
         ELSE IF (I==N)  THEN
            Line%Dp(:,I) = 0.25*p%rhoW*Line%Cdn*    d*Line%l(N) * MagVp * Vp
            Line%Dq(:,I) = 0.25*p%rhoW*Line%Cdt* Pi*d*Line%l(N) * MagVq * Vq
         ELSE
            Line%Dp(:,I) = 0.25*p%rhoW*Line%Cdn*    d*(Line%l(I) + Line%l(I+1)) * MagVp * vp
            Line%Dq(:,I) = 0.25*p%rhoW*Line%Cdt* Pi*d*(Line%l(I) + Line%l(I+1)) * MagVq * vq
         END IF

         ! F-K force from fluid acceleration not implemented yet

         ! bottom contact (stiffness and damping, vertical-only for now)  - updated Nov 24 for general case where anchor and fairlead ends may deal with bottom contact forces
         ! bottom contact - updated throughout October 2021 for seabed bathymetry and friction models
         
         ! interpolate the local depth from the bathymetry grid and return the vector normal to the seabed slope
         CALL getDepthFromBathymetry(m%BathymetryGrid, m%BathGrid_Xs, m%BathGrid_Ys, Line%r(1,I), Line%r(2,I), depth, nvec)

         IF (Line%r(3,I) < -depth) THEN   ! for every line node at or below the seabed
         
            ! calculate the velocity components of the node relative to the seabed
            Vn = DOT_PRODUCT( Line%rd(:,I), nvec) * nvec  ! velocity component normal to the local seabed slope
            Vsb = Line%rd(:,I) - Vn                       ! velocity component along (tangent to) the seabed
            
            ! calculate the normal contact force on the line node
            IF (I==0) THEN
               Fn = ( (-depth - Line%r(3,I))*nvec(3)*nvec*p%kBot - Vn*p%cBot) * 0.5*d*(            Line%l(I+1) )
            ELSE IF (I==N) THEN
               Fn = ( (-depth - Line%r(3,I))*nvec(3)*nvec*p%kBot - Vn*p%cBot) * 0.5*d*(Line%l(I)               )
            ELSE
               Fn = ( (-depth - Line%r(3,I))*nvec(3)*nvec*p%kBot - Vn*p%cBot) * 0.5*d*(Line%l(I) + Line%l(I+1) )
            END IF
         
            ! calculate the axial and transverse components of the node velocity vector along the seabed
            Va = DOT_PRODUCT( Vsb , Line%q(:,I) ) * Line%q(:,I)
            Vt = Vsb - Va
            
            ! calculate the magnitudes of each velocity
            VaMag = SQRT(Va(1)**2+Va(2)**2+Va(3)**2)
            VtMag = SQRT(Vt(1)**2+Vt(2)**2+Vt(3)**2)

            ! find the maximum possible kinetic friction force using transverse and axial kinetic friction coefficients
            FkTmax = p%mu_kT*SQRT(Fn(1)**2+Fn(2)**2+Fn(3)**2)
            FkAmax = p%mu_kA*SQRT(Fn(1)**2+Fn(2)**2+Fn(3)**2)
            ! turn the maximum kinetic friction forces into vectors in the direction of their velocities
            DO J = 1, 3
               IF (VtMag==0) THEN
                  FkT(J) = 0.0_DbKi
               ELSE
                  FkT(J) = FkTmax*Vt(J)/VtMag
               END IF
               IF (VaMag==0) THEN
                  FkA(J) = 0.0_DbKi
               ELSE
                  FkA(J) = FkAmax*Va(J)/VaMag
               END IF
            END DO
            ! calculate the ratio between the static and kinetic coefficients of friction
            !mc_T = p%mu_sT/p%mu_kT
            !mc_A = p%mu_sA/p%mu_kA
            
            ! calculate the transverse friction force
            IF (p%mu_kT*p%cv*VtMag > p%mc*FkTmax) THEN   ! if the friction force of the linear curve is greater than the maximum friction force allowed adjusted for static friction,
               FfT = -FkT                                ! then the friction force is the maximum kinetic friction force vector (constant part of the curve)
            ELSE                                         ! if the friction force of the linear curve is less than the maximum friction force allowed adjusted for static friction,
               FfT = -p%mu_kT*p%cv*Vt                    ! then the friction force is the calculated value of the linear line
            END IF
            ! calculate the axial friction force
            IF (p%mu_kA*p%cv*VaMag > p%mc*FkAmax) THEN   ! if the friction force of the linear curve is greater than the maximum friction force allowed adjusted for static friction,
               FfA = -FkA                                ! then the friction force is the maximum kinetic friction force vector (constant part of the curve)
            ELSE                                         ! if the friction force of the linear curve is less than the maximum friction force allowed adjusted for static friction,
               FfA = -p%mu_kA*p%cv*Va                    ! then the friction force is the calculated value of the linear line
            END IF
            ! NOTE: these friction forces have a negative sign here to indicate a force in the opposite direction of motion

            ! the total friction force is along the plane of the seabed slope, which is just the vector sum of the transverse and axial components
            Ff = FfT + FfA
            
         ELSE
            Fn = 0.0_DbKi
            Ff = 0.0_DbKi
         END IF
         
         
         ! the total force from bottom contact on the line node is the sum of the normal contact force and the friction force
         Line%B(:,I) = Fn + Ff

         ! total forces
         IF (I==0)  THEN
            Line%Fnet(:,I) = Line%T(:,1)                 + Line%Td(:,1)                  + Line%W(:,I) + Line%Dp(:,I) + Line%Dq(:,I) + Line%B(:,I) + Line%Bs(:,I)
         ELSE IF (I==N)  THEN
            Line%Fnet(:,I) =                -Line%T(:,N)                  - Line%Td(:,N) + Line%W(:,I) + Line%Dp(:,I) + Line%Dq(:,I) + Line%B(:,I) + Line%Bs(:,I)
         ELSE
            Line%Fnet(:,I) = Line%T(:,I+1) - Line%T(:,I) + Line%Td(:,I+1) - Line%Td(:,I) + Line%W(:,I) + Line%Dp(:,I) + Line%Dq(:,I) + Line%B(:,I) + Line%Bs(:,I)
         END IF

      END DO  ! I  - done looping through nodes

      ! loop through internal nodes and update their states  <<< should/could convert to matrix operations instead of all these loops
      DO I=1, N-1
         DO J=1,3

            ! calculate RHS constant (premultiplying force vector by inverse of mass matrix  ... i.e. rhs = S*Forces)
            Sum1 = 0.0_DbKi                               ! reset temporary accumulator <<< could turn this into a Line%a array to save and output node accelerations
            DO K = 1, 3
              Sum1 = Sum1 + Line%S(K,J,I) * Line%Fnet(K,I)   ! matrix-vector multiplication [S i]{Forces i}  << double check indices
            END DO ! K
            
            ! update states
            Xd(3*N-3 + 3*I-3 + J) = Line%rd(J,I);       ! dxdt = V  (velocities)
            Xd(        3*I-3 + J) = Sum1                ! dVdt = RHS * A  (accelerations)
            
         END DO ! J
      END DO  ! I


      ! check for NaNs
      DO J = 1, 6*(N-1)
         IF (Is_NaN(Xd(J))) THEN
            print *, "NaN detected at time ", Line%time, " in Line ", Line%IdNum, " in MoorDyn."
            IF (wordy > 1) THEN
               print *, "state derivatives:"
               print *, Xd
               
               
               
               print *, "m_i  p%rhoW   v_i Line%Can  Line%Cat"
               print *, m_i 
               print *, p%rhoW
               print *, v_i
               print *, Line%Can
               print *, Line%Cat
               
               print *, "Line%q"
               print *, Line%q
               
               print *, "Line%r"
               print *, Line%r
               
               
               print *, "Here is the mass matrix set"
               print *, Line%M
               
               print *, "Here is the inverted mass matrix set"
               print *, Line%S
               
               print *, "Here is the net force set"
               print *, Line%Fnet
            END IF            
            
            EXIT
         END IF
      END DO


      !   ! add force and mass of end nodes to the Connects they correspond to <<<<<<<<<<<< do this from Connection instead now!
      !   DO J = 1,3
      !      FairFtot(J) = FairFtot(J) + Line%F(J,N)
      !      AnchFtot(J) = AnchFtot(J) + Line%F(J,0)
      !      DO K = 1,3
      !         FairMtot(K,J) = FairMtot(K,J) + Line%M(K,J,N)
      !         AnchMtot(K,J) = AnchMtot(K,J) + Line%M(K,J,0)
      !      END DO
      !   END DO

   END SUBROUTINE Line_GetStateDeriv
   !=====================================================================


   !--------------------------------------------------------------
   SUBROUTINE Line_SetEndKinematics(Line, r_in, rd_in, t, topOfLine)

      TYPE(MD_Line),    INTENT(INOUT)  :: Line           ! the current Line object
      Real(DbKi),       INTENT(IN   )  :: r_in( 3)       ! state vector section for this line
      Real(DbKi),       INTENT(IN   )  :: rd_in(3)       ! state vector section for this line
      Real(DbKi),       INTENT(IN   )  :: t              ! instantaneous time
      INTEGER(IntKi),   INTENT(IN   )  :: topOfLine      ! 0 for end A (Node 0), 1 for end B (node N)

      Integer(IntKi)                   :: J      
      INTEGER(IntKi)                   :: inode
      
      IF (topOfLine==1) THEN
         inode = Line%N      
         Line%endTypeB = 0   ! set as ball rather than rigid connection (unless changed later by SetEndOrientation)
      ELSE
         inode = 0
         Line%endTypeA = 0   ! set as ball rather than rigid connection (unless changed later by SetEndOrientation)
      END IF 
      
      !Line%r( :,inode) = r_in
      !Line%rd(:,inode) = rd_in
      
      DO J = 1,3
         Line%r( :,inode) = r_in
         Line%rd(:,inode) = rd_in
      END DO
      
   !   print *, "SetEndKinematics of line ", Line%idNum, " top?:", topOfLine
   !   print *, r_in
   !   print *, Line%r( :,inode), "  - confirming, node ", inode
   !   print *, rd_in
  
      Line%time = t
   
   END SUBROUTINE Line_SetEndKinematics
   !--------------------------------------------------------------
   

   ! get force, moment, and mass of line at line end node
   !--------------------------------------------------------------
   SUBROUTINE Line_GetEndStuff(Line, Fnet_out, Moment_out, M_out, topOfLine)
   
      TYPE(MD_Line),    INTENT(INOUT)  :: Line           ! label for the current line, for convenience
      REAL(DbKi),       INTENT(  OUT)  :: Fnet_out(3)    ! net force on end node
      REAL(DbKi),       INTENT(  OUT)  :: Moment_out(3)  ! moment on end node (future capability)
      REAL(DbKi),       INTENT(  OUT)  :: M_out(3,3)     ! mass matrix of end node
      INTEGER(IntKi),   INTENT(IN   )  :: topOfLine      ! 0 for end A (Node 0), 1 for end B (node N)
      
      Integer(IntKi)                   :: J
!      INTEGER(IntKi)                   :: inode
      
      IF (topOfLine==1) THEN           ! end B of line
         Fnet_out   = Line%Fnet(:, Line%N)
         Moment_out = Line%endMomentB
         M_out      = Line%M(:,:, Line%N)
      ELSE                             ! end A of line
         Fnet_out   = Line%Fnet(:, 0)
         Moment_out = Line%endMomentA
         M_out      = Line%M(:,:, 0)
      END IF

   END SUBROUTINE Line_GetEndStuff
   !--------------------------------------------------------------
  
   ! Get bending stiffness vector from line end for use in computing orientation of zero-length rods
   !--------------------------------------------------------------
   SUBROUTINE Line_GetEndSegmentInfo(Line, q_EI_dl, topOfLine, rodEndB)

      TYPE(MD_Line),    INTENT(INOUT)  :: Line           ! label for the current line, for convenience
      REAL(DbKi),       INTENT(  OUT)  :: q_EI_dl(3)     ! EI/dl of the line end segment multiplied by the axis unit vector with the correct sign
      INTEGER(IntKi),   INTENT(IN   )  :: topOfLine      ! 0 for end A (Node 0), 1 for end B (node N)
      INTEGER(IntKi),   INTENT(IN   )  :: rodEndB        ! rodEndB=0 means the line is attached to Rod end A, =1 means attached to Rod end B (implication for unit vector sign)
      
      REAL(DbKi)   :: qEnd(3)
      REAL(DbKi)   :: dlEnd
            
      if (topOfLine==1) then
         CALL UnitVector(Line%r(:,Line%N-1), Line%r(:,Line%N), qEnd, dlEnd)  ! unit vector of last line segment
         if (rodEndB == 0) then
            q_EI_dl =  qEnd*Line%EI/dlEnd  ! -----line----->[A==ROD==>B]
         else
            q_EI_dl = -qEnd*Line%EI/dlEnd  ! -----line----->[B==ROD==>A]
         end if
      else
         CALL UnitVector(Line%r(:,0       ), Line%r(:,1     ), qEnd, dlEnd)  ! unit vector of first line segment
         if (rodEndB == 0) then
            q_EI_dl = -qEnd*Line%EI/dlEnd  ! <----line-----[A==ROD==>B]
         else
            q_EI_dl =  qEnd*Line%EI/dlEnd  ! <----line-----[B==ROD==>A]
         end if
      end if
      
   END SUBROUTINE Line_GetEndSegmentInfo
   !--------------------------------------------------------------


   ! set end node unit vector of a line (this is called when attached to a Rod, only applicable for bending stiffness)
   !--------------------------------------------------------------
   SUBROUTINE Line_SetEndOrientation(Line, qin, topOfLine, rodEndB)

      TYPE(MD_Line),    INTENT(INOUT)  :: Line           ! label for the current line, for convenience
      REAL(DbKi),       INTENT(IN   )  :: qin(3)         ! the rod's axis unit vector
      INTEGER(IntKi),   INTENT(IN   )  :: topOfLine      ! 0 for end A (Node 0), 1 for end B (node N)
      INTEGER(IntKi),   INTENT(IN   )  :: rodEndB        ! =0 means the line is attached to Rod end A, =1 means attached to Rod end B (implication for unit vector sign)

      if (topOfLine==1) then
      
         Line%endTypeB = 1                  ! indicate attached to Rod (at every time step, just in case line gets detached)
         
         if (rodEndB==1) then
            Line%q(:,Line%N) = -qin   ! -----line----->[B<==ROD==A]
         else
            Line%q(:,Line%N) =  qin   ! -----line----->[A==ROD==>B]
         end if
      else
         
         Line%endTypeA = 1                  ! indicate attached to Rod (at every time step, just in case line gets detached)                 ! indicate attached to Rod
         
         if (rodEndB==1) then
            Line%q(:,0     ) =  qin   ! [A==ROD==>B]-----line----->
         else
            Line%q(:,0     ) = -qin   ! [B<==ROD==A]-----line----->
         end if
      end if

   END SUBROUTINE Line_SetEndOrientation
   !--------------------------------------------------------------


END MODULE MoorDyn_Line
