!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015  Matthew Hall
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
MODULE MoorDyn

   USE MoorDyn_Types
   USE MoorDyn_IO
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: MD_ProgDesc = ProgDesc( 'MoorDyn', 'v1.01.02F', '8-Apr-2016' )


   PUBLIC :: MD_Init
   PUBLIC :: MD_UpdateStates
   PUBLIC :: MD_CalcOutput
   PUBLIC :: MD_CalcContStateDeriv
   PUBLIC :: MD_End

CONTAINS

   !=========================================   MD_Init   ===================================
   SUBROUTINE MD_Init(InitInp, u, p, x, xd, z, other, y, m, DTcoupling, InitOut, ErrStat, ErrMsg)

      IMPLICIT NONE

      TYPE(MD_InitInputType),       INTENT(INOUT)  :: InitInp     ! INTENT(INOUT) : Input data for initialization routine
      TYPE(MD_InputType),           INTENT(  OUT)  :: u           ! INTENT( OUT) : An initial guess for the input; input mesh must be defined
      TYPE(MD_ParameterType),       INTENT(  OUT)  :: p           ! INTENT( OUT) : Parameters
      TYPE(MD_ContinuousStateType), INTENT(  OUT)  :: x           ! INTENT( OUT) : Initial continuous states
      TYPE(MD_DiscreteStateType),   INTENT(  OUT)  :: xd          ! INTENT( OUT) : Initial discrete states
      TYPE(MD_ConstraintStateType), INTENT(  OUT)  :: z           ! INTENT( OUT) : Initial guess of the constraint states
      TYPE(MD_OtherStateType),      INTENT(  OUT)  :: other       ! INTENT( OUT) : Initial other states
      TYPE(MD_OutputType),          INTENT(  OUT)  :: y           ! INTENT( OUT) : Initial system outputs (outputs are not calculated; only the output mesh is initialized)
      TYPE(MD_MiscVarType),         INTENT(  OUT)  :: m           ! INTENT( OUT) : Initial misc/optimization variables
      REAL(DbKi),                   INTENT(INOUT)  :: DTcoupling  ! Coupling interval in seconds: the rate that Output is the actual coupling interval
      TYPE(MD_InitOutputType),      INTENT(INOUT)  :: InitOut     ! Output for initialization routine
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
      REAL(DbKi)                                   :: t              ! instantaneous time, to be used during IC generation
      INTEGER(IntKi)                               :: I              ! index
      INTEGER(IntKi)                               :: J              ! index
      INTEGER(IntKi)                               :: K              ! index
      INTEGER(IntKi)                               :: Converged      ! flag indicating whether the dynamic relaxation has converged
      INTEGER(IntKi)                               :: N              ! convenience integer for readability: number of segments in the line
      REAL(ReKi)                                   :: Pos(3)         ! array for setting absolute fairlead positions in mesh
      REAL(DbKi)                                   :: TransMat(3,3)  ! rotation matrix for setting fairlead positions correctly if there is initial platform rotation
      REAL(DbKi), ALLOCATABLE                      :: FairTensIC(:,:)! array of size Nfairs, 3 to store three latest fairlead tensions of each line
      CHARACTER(20)                                :: TempString     ! temporary string for incidental use
      INTEGER(IntKi)                               :: ErrStat2       ! Error status of the operation
      CHARACTER(ErrMsgLen)                         :: ErrMsg2        ! Error message if ErrStat2 /= ErrID_None
      
      TYPE(MD_InputType)    :: uArray(1)    ! a size-one array for u to make call to TimeStep happy
      REAL(DbKi)            :: utimes(1)    ! a size-one array saying time is 0 to make call to TimeStep happy  

      CHARACTER(MaxWrScrLen)                       :: Message


      ErrStat = ErrID_None
      ErrMsg  = ""

      ! Initialize the NWTC Subroutine Library
      CALL NWTC_Init( )

      ! Display the module information
      CALL DispNVD( MD_ProgDesc )
      InitOut%Ver = MD_ProgDesc


      !---------------------------------------------------------------------------------------------
      !                   Get all the inputs taken care of
      !---------------------------------------------------------------------------------------------


      ! set environmental parameters from input data and error check
      ! (should remove these values as options from MoorDyn input file for consistency?)

      p%g        = InitInp%g
      p%WtrDpth  = InitInp%WtrDepth
      p%rhoW     = InitInp%rhoW

      p%RootName = TRIM(InitInp%RootName)//'.MD'  ! all files written from this module will have this root name


      ! call function that reads input file and creates cross-referenced Connect and Line objects
      CALL MDIO_ReadInput(InitInp, p, m, ErrStat2, ErrMsg2)
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN

      ! process the OutList array and set up the index arrays for the requested output quantities
      CALL MDIO_ProcessOutList(InitInp%OutList, p, m, y, InitOut, ErrStat2, ErrMsg2 )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN


      !-------------------------------------------------------------------------------------------------
      !          Connect mooring system together and make necessary allocations
      !-------------------------------------------------------------------------------------------------

      CALL WrNR( '   Creating mooring system.  ' )

      p%NFairs = 0   ! this is the number of "vessel" type Connections.  being consistent with MAP terminology
      p%NConns = 0   ! this is the number of "connect" type Connections.  not to be confused with NConnects, the number of Connections
      p%NAnchs = 0   ! this is the number of "fixed" type Connections.

      ! cycle through Connects and identify Connect types
      DO I = 1, p%NConnects
               
         TempString = m%ConnectList(I)%type
         CALL Conv2UC(TempString)
         if (TempString == 'FIXED') then
            m%ConnectList(I)%TypeNum = 0
            p%NAnchs = p%NAnchs + 1
         else if (TempString == 'VESSEL') then
            m%ConnectList(I)%TypeNum = 1
            p%NFairs = p%NFairs + 1             ! if a vessel connection, increment fairlead counter
         else if (TempString == 'CONNECT') then
            m%ConnectList(I)%TypeNum = 2
            p%NConns = p%NConns + 1
         else
            CALL CheckError( ErrID_Fatal, 'Error in provided Connect type.  Must be fixed, vessel, or connect.' )
            RETURN
         END IF
      END DO

      CALL WrScr(trim(Num2LStr(p%NFairs))//' fairleads, '//trim(Num2LStr(p%NAnchs))//' anchors, '//trim(Num2LStr(p%NConns))//' connects.')


      ! allocate fairleads list
      ALLOCATE ( m%FairIdList(p%NFairs), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         CALL CheckError( ErrID_Fatal, 'Error allocating space for FairIdList array.')
         RETURN
      END IF

      ! allocate connect list
      ALLOCATE ( m%ConnIdList(p%NConns), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         CALL CheckError( ErrID_Fatal, 'Error allocating space for ConnIdList array.')
         RETURN
      END IF


      ! now go back through and record the fairlead Id numbers (this is all the "connecting" that's required)
      J = 1  ! counter for fairlead number
      K = 1  ! counter for connect number
      DO I = 1,p%NConnects
         IF (m%ConnectList(I)%TypeNum == 1) THEN
           m%FairIdList(J) = I             ! if a vessel connection, add ID to list
           J = J + 1
         ELSE IF (m%ConnectList(I)%TypeNum == 2) THEN
           m%ConnIdList(K) = I             ! if a connect connection, add ID to list
           K = K + 1
         END IF
      END DO


      ! go through lines and allocate variables
      DO I = 1, p%NLines
         CALL SetupLine( m%LineList(I), m%LineTypeList(m%LineList(I)%PropsIdNum), p%rhoW ,  ErrStat2, ErrMsg2)
            CALL CheckError( ErrStat2, ErrMsg2 )
            IF (ErrStat >= AbortErrLev) RETURN
      END DO


      !------------------------------------------------------------------------------------
      !                               prepare state vector
      !------------------------------------------------------------------------------------

      ! allocate list of starting state vector indices for each line  - does this belong elsewhere?
      ALLOCATE ( m%LineStateIndList(p%NLines), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
        CALL CheckError(ErrID_Fatal, ' Error allocating LineStateIndList array.')
        RETURN
      END IF


      ! figure out required size of state vector and how it will be apportioned to Connect and Lines (J is keeping track of the growing size of the state vector)
      J = p%NConns*6   ! start index of first line's states (added six state variables for each "connect"-type connection)

      DO I = 1, p%NLines
         m%LineStateIndList(I) = J+1            ! assign start index of each line
         J = J + 6*(m%LineList(I)%N - 1)  !add 6 state variables for each internal node
      END DO


      ! allocate state vector for RK2 based on size just calculated
      ALLOCATE ( x%states(J), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
        ErrMsg  = ' Error allocating state vector.'
        !CALL CleanUp()
        RETURN
      END IF


      ! get header information for the FAST output file   <<< what does this mean?


      !--------------------------------------------------------------------------
      !             create i/o meshes for fairlead positions and forces
      !--------------------------------------------------------------------------

      ! create input mesh for fairlead kinematics
      CALL MeshCreate(BlankMesh=u%PtFairleadDisplacement , &
                    IOS= COMPONENT_INPUT           , &
                    Nnodes=p%NFairs                , &
                    TranslationDisp=.TRUE.         , &
                    TranslationVel=.TRUE.          , &
                    ErrStat=ErrStat2                , &
                    ErrMess=ErrMsg2)

      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN


      !  --------------------------- set up initial condition of each fairlead -------------------------------
      DO i = 1,p%NFairs

         Pos(1) = m%ConnectList(m%FairIdList(i))%conX ! set relative position of each fairlead i (I'm pretty sure this is just relative to ptfm origin)
         Pos(2) = m%ConnectList(m%FairIdList(i))%conY
         Pos(3) = m%ConnectList(m%FairIdList(i))%conZ

         CALL MeshPositionNode(u%PtFairleadDisplacement,i,Pos,ErrStat2,ErrMsg2)! "assign the coordinates of each node in the global coordinate space"

         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN


         ! set offset position of each node to according to initial platform position
         CALL SmllRotTrans('initial fairlead positions due to platform rotation', InitInp%PtfmInit(4),InitInp%PtfmInit(5),InitInp%PtfmInit(6), TransMat, '', ErrStat2, ErrMsg2)  ! account for possible platform rotation

         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN

         ! Apply initial platform rotations and translations (fixed Jun 19, 2015)
         u%PtFairleadDisplacement%TranslationDisp(1,i) = InitInp%PtfmInit(1) + Transmat(1,1)*Pos(1) + Transmat(2,1)*Pos(2) + TransMat(3,1)*Pos(3) - Pos(1)
         u%PtFairleadDisplacement%TranslationDisp(2,i) = InitInp%PtfmInit(2) + Transmat(1,2)*Pos(1) + Transmat(2,2)*Pos(2) + TransMat(3,2)*Pos(3) - Pos(2)
         u%PtFairleadDisplacement%TranslationDisp(3,i) = InitInp%PtfmInit(3) + Transmat(1,3)*Pos(1) + Transmat(2,3)*Pos(2) + TransMat(3,3)*Pos(3) - Pos(3)

         ! set velocity of each node to zero
         u%PtFairleadDisplacement%TranslationVel(1,i) = 0.0_DbKi
         u%PtFairleadDisplacement%TranslationVel(2,i) = 0.0_DbKi
         u%PtFairleadDisplacement%TranslationVel(3,i) = 0.0_DbKi
         
         !print *, 'Fairlead ', i, ' z TranslationDisp at start is ', u%PtFairleadDisplacement%TranslationDisp(3,i)
         !print *, 'Fairlead ', i, ' z Position at start is ', u%PtFairleadDisplacement%Position(3,i)


         ! set each node as a point element
         CALL MeshConstructElement(u%PtFairleadDisplacement, ELEMENT_POINT, ErrStat2, ErrMsg2, i)

         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN

      END DO    ! I


      CALL MeshCommit ( u%PtFairleadDisplacement, ErrStat, ErrMsg )

      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN


      ! copy the input fairlead kinematics mesh to make the output mesh for fairlead loads, PtFairleadLoad
      CALL MeshCopy ( SrcMesh  = u%PtFairleadDisplacement,   DestMesh = y%PtFairleadLoad, &
                      CtrlCode = MESH_SIBLING,               IOS      = COMPONENT_OUTPUT, &
                      Force    = .TRUE.,                     ErrStat  = ErrStat2, ErrMess=ErrMsg2 )

      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN


      ! --------------------------------------------------------------------
      !   go through all Connects and set position based on input file
      ! --------------------------------------------------------------------

      ! first do it for all connections (connect and anchor types will be saved)
      DO I = 1, p%NConnects
         m%ConnectList(I)%r(1) = m%ConnectList(I)%conX
         m%ConnectList(I)%r(2) = m%ConnectList(I)%conY
         m%ConnectList(I)%r(3) = m%ConnectList(I)%conZ
         m%ConnectList(I)%rd(1) = 0.0_DbKi
         m%ConnectList(I)%rd(2) = 0.0_DbKi
         m%ConnectList(I)%rd(3) = 0.0_DbKi
      END DO

      ! then do it for fairlead types
      DO I = 1,p%NFairs
         DO J = 1, 3
            m%ConnectList(m%FairIdList(I))%r(J)  = u%PtFairleadDisplacement%Position(J,I) + u%PtFairleadDisplacement%TranslationDisp(J,I)
            m%ConnectList(m%FairIdList(I))%rd(J) = 0.0_DbKi
         END DO
      END DO

      ! for connect types, write the coordinates to the state vector
      DO I = 1,p%NConns
         x%states(6*I-2:6*I)   = m%ConnectList(m%ConnIdList(I))%r  ! double check order of r vs rd
         x%states(6*I-5:6*I-3) = m%ConnectList(m%ConnIdList(I))%rd
      END DO

      ! --------------------------------------------------------------------
      !          open output file(s) and write header lines
      CALL MDIO_OpenOutput( InitInp%FileName, p, m, InitOut, ErrStat2, ErrMsg2 )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN
      ! --------------------------------------------------------------------



      ! --------------------------------------------------------------------
      ! size active tensioning inputs arrays based on highest channel number read from input file for now <<<<<<<
      ! --------------------------------------------------------------------
      
      ! find the highest channel number
      N = 0
      DO I = 1, p%NLines
         IF ( m%LineList(I)%CtrlChan > N ) then
            N = m%LineList(I)%CtrlChan       
         END IF
      END DO   
      
      ! allocate the input arrays (if any requested)
      if (N > 0) then
         call AllocAry( u%DeltaL, N, 'u%DeltaL', ErrStat2, ErrMsg2 )
            call CheckError( ErrStat2, ErrMsg2 )
            if (ErrStat >= AbortErrLev) return
            u%DeltaL =  0.0_ReKi
         call AllocAry( u%DeltaLdot, N, 'u%DeltaLdot', ErrStat2, ErrMsg2 )
            call CheckError( ErrStat2, ErrMsg2 )
            if (ErrStat >= AbortErrLev) return
            u%DeltaLdot =  0.0_ReKi
         call AllocAry( InitOut%CableCChanRqst, N, 'CableCChanRqst', ErrStat2, ErrMsg2 )
            call CheckError( ErrStat2, ErrMsg2 )
            if (ErrStat >= AbortErrLev) return
         InitOut%CableCChanRqst = .FALSE.    ! Initialize to false
         do J=1,p%NLines
            if (m%LineList(J)%CtrlChan > 0)  InitOut%CableCChanRqst(m%LineList(J)%CtrlChan) = .TRUE.
         enddo
      endif


      ! --------------------------------------------------------------------
      ! go through lines and initialize internal node positions using Catenary()
      ! --------------------------------------------------------------------
      DO I = 1, p%NLines

         N = m%LineList(I)%N ! for convenience
         
         !TODO: apply any initial adjustment of line length from active tensioning <<<<<<<<<<<<
         ! >>> maybe this should be skipped <<<<

         ! set end node positions and velocities from connect objects
         m%LineList(I)%r(:,N) = m%ConnectList(m%LineList(I)%FairConnect)%r
         m%LineList(I)%r(:,0) = m%ConnectList(m%LineList(I)%AnchConnect)%r
         m%LineList(I)%rd(:,N) = (/ 0.0, 0.0, 0.0 /)  ! set anchor end velocities to zero
         m%LineList(I)%rd(:,0) = (/ 0.0, 0.0, 0.0 /)  ! set fairlead end velocities to zero

         ! set initial line internal node positions using quasi-static model or straight-line interpolation from anchor to fairlead
         CALL InitializeLine( m%LineList(I), m%LineTypeList(m%LineList(I)%PropsIdNum), p%rhoW ,  ErrStat2, ErrMsg2)
            CALL CheckError( ErrStat2, ErrMsg2 )
            IF (ErrStat >= AbortErrLev) RETURN
            IF (ErrStat >= ErrId_Warn) CALL WrScr("   Catenary solver failed for one or more lines.  Using linear node spacing.")  ! make this statement more accurate

         ! assign the resulting internal node positions to the integrator initial state vector! (velocities leave at 0)
         DO J = 1, N-1
           DO K = 1, 3
             x%states(m%LineStateIndList(I) + 3*N-3 + 3*J-3 + K-1 ) = m%LineList(I)%r(K,J) ! assign position
             x%states(m%LineStateIndList(I)         + 3*J-3 + K-1 ) = 0.0_DbKi ! assign velocities (of zero)
           END DO
         END DO

      END DO    !I = 1, p%NLines


!      ! try writing output for troubleshooting purposes (TEMPORARY)
!      CALL MDIO_WriteOutputs(-1.0_DbKi, p, m, y, ErrStat, ErrMsg)
!      IF ( ErrStat >= AbortErrLev ) THEN
!         ErrMsg = ' Error in MDIO_WriteOutputs: '//TRIM(ErrMsg)
!         RETURN
!      END IF


      ! --------------------------------------------------------------------
      !           do dynamic relaxation to get ICs
      ! --------------------------------------------------------------------

      CALL WrScr("   Finalizing ICs using dynamic relaxation."//NewLine)  ! newline because next line writes over itself

      ! boost drag coefficient of each line type
      DO I = 1, p%NTypes
         m%LineTypeList(I)%Cdn = m%LineTypeList(I)%Cdn * InitInp%CdScaleIC
         m%LineTypeList(I)%Cdt = m%LineTypeList(I)%Cdt * InitInp%CdScaleIC
      END DO

      ! allocate array holding three latest fairlead tensions
      ALLOCATE ( FairTensIC(p%NFairs,3), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         CALL CheckError( ErrID_Fatal, ErrMsg2 )
         RETURN
      END IF

      ! initialize fairlead tension memory at zero
      DO J = 1,p%NFairs
         DO I = 1, 3
            FairTensIC(J,I) = 0.0_DbKi
         END DO
      END DO

      t = 0.0_DbKi     ! start time at zero

      ! because TimeStep wants an array...
      call MD_CopyInput( u, uArray(1), MESH_NEWCOPY, ErrStat2, ErrMsg2 )


      DO I = 1, ceiling(InitInp%TMaxIC/InitInp%DTIC)   ! loop through IC gen time steps, up to maximum

         ! integrate the EOMs one DTIC s time step
         CALL TimeStep ( t, InitInp%DTIC, uArray, utimes, p, x, xd, z, other, m, ErrStat, ErrMsg )
            CALL CheckError( ErrStat2, ErrMsg2 )
            IF (ErrStat >= AbortErrLev) RETURN

         ! store new fairlead tension (and previous fairlead tensions for comparison)
         DO J = 1, p%NFairs
            FairTensIC(J,3) = FairTensIC(J,2)
            FairTensIC(J,2) = FairTensIC(J,1)
            FairTensIC(J,1) = TwoNorm(m%ConnectList(m%FairIdList(J))%Ftot(:))
         END DO

         ! provide status message
         ! bjj: putting this in a string so we get blanks to cover up previous values (if current string is shorter than previous one)
         Message = '   t='//trim(Num2LStr(t))//'  FairTen 1: '//trim(Num2LStr(FairTensIC(1,1)))// &
                        ', '//trim(Num2LStr(FairTensIC(1,2)))//', '//trim(Num2LStr(FairTensIC(1,3))) 
         CALL WrOver( Message )

         ! check for convergence (compare current tension at each fairlead with previous two values)
         IF (I > 2) THEN
            Converged = 1
            DO J = 1, p%NFairs   ! check for non-convergence
               IF (( abs( FairTensIC(J,1)/FairTensIC(J,2) - 1.0 ) > InitInp%threshIC ) .OR. ( abs( FairTensIC(J,1)/FairTensIC(J,3) - 1.0 ) > InitInp%threshIC ) ) THEN
                  Converged = 0
                  EXIT
               END IF
            END DO

            IF (Converged == 1)  THEN  ! (J == p%NFairs) THEN   ! if we made it with all cases satisfying the threshold
               CALL WrScr('   Fairlead tensions converged to '//trim(Num2LStr(100.0*InitInp%threshIC))//'% after '//trim(Num2LStr(t))//' seconds.')
               EXIT  ! break out of the time stepping loop
            END IF
         END IF

         IF (I == ceiling(InitInp%TMaxIC/InitInp%DTIC) ) THEN
            CALL WrScr('   Fairlead tensions did not converge within TMaxIC='//trim(Num2LStr(InitInp%TMaxIC))//' seconds.')
            !ErrStat = ErrID_Warn
            !ErrMsg = '  MD_Init: ran dynamic convergence to TMaxIC without convergence'
         END IF

      END DO ! I ... looping through time steps

      CALL MD_DestroyInput( uArray(1), ErrStat2, ErrMsg2 )

      ! UNboost drag coefficient of each line type
      DO I = 1, p%NTypes
         m%LineTypeList(I)%Cdn = m%LineTypeList(I)%Cdn / InitInp%CdScaleIC
         m%LineTypeList(I)%Cdt = m%LineTypeList(I)%Cdt / InitInp%CdScaleIC
      END DO


      p%dtCoupling = DTcoupling  ! store coupling time step for use in updatestates

      other%dummy = 0
      xd%dummy    = 0
      z%dummy     = 0

   CONTAINS

      SUBROUTINE CheckError(ErrID,Msg)
         ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev

            ! Passed arguments
         INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
         CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
         CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)

         ! Set error status/message;
         IF ( ErrID /= ErrID_None ) THEN

            IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine   ! if there's a pre-existing warning/error, retain the message and start a new line

            ErrMsg = TRIM(ErrMsg)//' MD_Init:'//TRIM(Msg)
            ErrStat = MAX(ErrStat, ErrID)

            ! Clean up if we're going to return on error: close files, deallocate local arrays


            IF ( ErrStat >= AbortErrLev ) THEN                
               IF (ALLOCATED(m%FairIdList       ))  DEALLOCATE(m%FairIdList       )
               IF (ALLOCATED(m%ConnIdList       ))  DEALLOCATE(m%ConnIdList       )
               IF (ALLOCATED(m%LineStateIndList ))  DEALLOCATE(m%LineStateIndList )
               IF (ALLOCATED(x%states           ))  DEALLOCATE(x%states           )
               IF (ALLOCATED(FairTensIC         ))  DEALLOCATE(FairTensIC         ) 
            END IF
         END IF

      END SUBROUTINE CheckError

   END SUBROUTINE MD_Init
   !==============================================================================================




   !==============================================================================================
   SUBROUTINE MD_UpdateStates( t, n, u, utimes, p, x, xd, z, other, m, ErrStat, ErrMsg)

      REAL(DbKi)                      , INTENT(IN   ) :: t
      INTEGER(IntKi)                  , INTENT(IN   ) :: n
      TYPE(MD_InputType)              , INTENT(INOUT) :: u(:)       ! INTENT(INOUT) ! had to change this to INOUT
      REAL(DbKi)                      , INTENT(IN   ) :: utimes(:)
      TYPE(MD_ParameterType)          , INTENT(IN   ) :: p          ! INTENT(IN   )
      TYPE(MD_ContinuousStateType)    , INTENT(INOUT) :: x          ! INTENT(INOUT)
      TYPE(MD_DiscreteStateType)      , INTENT(INOUT) :: xd         ! INTENT(INOUT)
      TYPE(MD_ConstraintStateType)    , INTENT(INOUT) :: z          ! INTENT(INOUT)
      TYPE(MD_OtherStateType)         , INTENT(INOUT) :: other      ! INTENT(INOUT)
      TYPE(MD_MiscVarType)            , INTENT(INOUT) :: m          ! INTENT(INOUT)
      INTEGER(IntKi)                  , INTENT(  OUT) :: ErrStat    ! Error status of the operation
      CHARACTER(*)                    , INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

      INTEGER(IntKi)                                  :: ErrStat2   ! Error status of the operation
      CHARACTER(ErrMsgLen)                            :: ErrMsg2    ! Error message if ErrStat2 /= ErrID_None

! moved to TimeStep      TYPE(MD_InputType)                              :: u_interp   !
      INTEGER(IntKi)                                  :: nTime

      REAL(DbKi)                                      :: t2         ! copy of time passed to TimeStep


      nTime = size(u) ! the number of times of input data provided?

      t2 = t

! >>> removing this section and putting it inside loop of TimeStep (to be done at every time step) <<<
!      ! create space for arrays/meshes in u_interp
!      CALL MD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
!         CALL CheckError( ErrStat2, ErrMsg2 )
!         IF (ErrStat >= AbortErrLev) RETURN
!
!      ! interpolate input mesh to correct time
!      CALL MD_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat2, ErrMsg2)
!         CALL CheckError( ErrStat2, ErrMsg2 )
!         IF (ErrStat >= AbortErrLev) RETURN
!
!
!      ! go through fairleads and apply motions from driver
!      DO I = 1, p%NFairs
!         DO J = 1,3
!            m%ConnectList(m%FairIdList(I))%r(J)  = u_interp%PtFairleadDisplacement%Position(J,I) + u_interp%PtFairleadDisplacement%TranslationDisp(J,I)
!            m%ConnectList(m%FairIdList(I))%rd(J) = u_interp%PtFairleadDisplacement%TranslationVel(J,I)  ! is this right? <<<
!         END DO
!      END DO
!

      ! call function that loops through mooring model time steps
      CALL TimeStep ( t2, p%dtCoupling, u, utimes, p, x, xd, z, other, m, ErrStat2, ErrMsg2 )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN


      ! clean up input interpolation stuff
  ! moved to TimeStep    CALL MD_DestroyInput(u_interp, ErrStat, ErrMsg)


   CONTAINS

      SUBROUTINE CheckError(ErrId, Msg)
        ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev

         INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
         CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


         IF ( ErrID /= ErrID_None ) THEN

            IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine  ! keep existing error message if there is one
            ErrMsg = TRIM(ErrMsg)//' MD_UpdateStates:'//TRIM(Msg)      ! add current error message
            ErrStat = MAX(ErrStat, ErrID)

            CALL WrScr( ErrMsg )  ! do this always or only if warning level?

            IF( ErrStat > ErrID_Warn ) THEN
       !         CALL MD_DestroyInput( u_interp, ErrStat, ErrMsg )
                RETURN
            END IF
         END IF

      END SUBROUTINE CheckError

   END SUBROUTINE MD_UpdateStates
   !========================================================================================



   !========================================================================================
   SUBROUTINE MD_CalcOutput( t, u, p, x, xd, z, other, y, m, ErrStat, ErrMsg )

      REAL(DbKi)                     , INTENT(IN   ) :: t
      TYPE( MD_InputType )           , INTENT(IN   ) :: u       ! INTENT(IN   )
      TYPE( MD_ParameterType )       , INTENT(IN   ) :: p       ! INTENT(IN   )
      TYPE( MD_ContinuousStateType ) , INTENT(IN   ) :: x       ! INTENT(IN   )
      TYPE( MD_DiscreteStateType )   , INTENT(IN   ) :: xd      ! INTENT(IN   )
      TYPE( MD_ConstraintStateType ) , INTENT(IN   ) :: z       ! INTENT(IN   )
      TYPE( MD_OtherStateType )      , INTENT(IN   ) :: other   ! INTENT(IN   )
      TYPE( MD_OutputType )          , INTENT(INOUT) :: y       ! INTENT(INOUT)
      TYPE(MD_MiscVarType)           , INTENT(INOUT) :: m       ! INTENT(INOUT)
      INTEGER(IntKi)                 , INTENT(INOUT) :: ErrStat
      CHARACTER(*)                   , INTENT(INOUT) :: ErrMsg

      TYPE(MD_ContinuousStateType)                   :: dxdt    ! time derivatives of continuous states (initialized in CalcContStateDeriv)
      INTEGER(IntKi)                                 :: I       ! counter
      INTEGER(IntKi)                                 :: J       ! counter

      INTEGER(IntKi)                                 :: ErrStat2   ! Error status of the operation
      CHARACTER(ErrMsgLen)                           :: ErrMsg2    ! Error message if ErrStat2 /= ErrID_None


      ! below updated to make sure outputs are current (based on provided x and u)  - similar to what's in UpdateStates

      ! go through fairleads and apply motions from driver
      DO I = 1, p%NFairs
         DO J = 1,3
            m%ConnectList(m%FairIdList(I))%r(J)  = u%PtFairleadDisplacement%Position(J,I) + u%PtFairleadDisplacement%TranslationDisp(J,I)
            m%ConnectList(m%FairIdList(I))%rd(J) = u%PtFairleadDisplacement%TranslationVel(J,I)  ! is this right? <<<
         END DO
      END DO

      ! call CalcContStateDeriv in order to run model and calculate dynamics with provided x and u
      CALL MD_CalcContStateDeriv( t, u, p, x, xd, z, other, m, dxdt, ErrStat, ErrMsg )


      ! assign net force on fairlead Connects to the output mesh
      DO i = 1, p%NFairs
         DO J=1,3
            y%PtFairleadLoad%Force(J,I) = m%ConnectList(m%FairIdList(I))%Ftot(J)
         END DO
      END DO


      ! calculate outputs (y%WriteOutput) for glue code and write any m outputs to MoorDyn output files
      CALL MDIO_WriteOutputs(t, p, m, y, ErrStat2, ErrMsg2)
      CALL CheckError(ErrStat2, 'In MDIO_WriteOutputs: '//trim(ErrMsg2))
      IF ( ErrStat >= AbortErrLev ) RETURN


      ! destroy dxdt
      CALL MD_DestroyContState( dxdt, ErrStat2, ErrMsg2)
      CALL CheckError(ErrStat2, 'When destroying dxdt: '//trim(ErrMsg2))
      IF ( ErrStat >= AbortErrLev ) RETURN



   CONTAINS

      SUBROUTINE CheckError(ErrId, Msg)
        ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev

         INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
         CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


         IF ( ErrID /= ErrID_None ) THEN

            IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine  ! keep existing error message if there is one
            ErrMsg = TRIM(ErrMsg)//' MD_CalcOutput:'//TRIM(Msg)      ! add current error message
            ErrStat = MAX(ErrStat, ErrID)

            CALL WrScr( ErrMsg )  ! do this always or only if warning level? <<<<<<<<<<<<<<<<<<<<<< probably should remove all instances

            IF( ErrStat > ErrID_Warn ) THEN
                CALL MD_DestroyContState( dxdt, ErrStat2, ErrMsg2)
            END IF
         END IF

      END SUBROUTINE CheckError

   END SUBROUTINE MD_CalcOutput
   !=============================================================================================


   !=============================================================================================
   SUBROUTINE MD_CalcContStateDeriv( t, u, p, x, xd, z, other, m, dxdt, ErrStat, ErrMsg )
   ! Tight coupling routine for computing derivatives of continuous states
   ! this is modelled off what used to be subroutine DoRHSmaster

      REAL(DbKi),                         INTENT(IN )    :: t       ! Current simulation time in seconds
      TYPE(MD_InputType),                 INTENT(IN )    :: u       ! Inputs at t
      TYPE(MD_ParameterType),             INTENT(IN )    :: p       ! Parameters
      TYPE(MD_ContinuousStateType),       INTENT(IN )    :: x       ! Continuous states at t
      TYPE(MD_DiscreteStateType),         INTENT(IN )    :: xd      ! Discrete states at t
      TYPE(MD_ConstraintStateType),       INTENT(IN )    :: z       ! Constraint states at t
      TYPE(MD_OtherStateType),            INTENT(IN )    :: other   ! Other states at t
      TYPE(MD_MiscVarType),               INTENT(INOUT)  :: m       ! misc/optimization variables
      TYPE(MD_ContinuousStateType),       INTENT(  OUT)  :: dxdt    ! Continuous state derivatives at t
      INTEGER(IntKi),                     INTENT( OUT)   :: ErrStat ! Error status of the operation
      CHARACTER(*),                       INTENT( OUT)   :: ErrMsg  ! Error message if ErrStat /= ErrID_None


      INTEGER(IntKi)                                     :: L       ! index
      INTEGER(IntKi)                                     :: J       ! index
      INTEGER(IntKi)                                     :: K       ! index
      INTEGER(IntKi)                                     :: Istart  ! start index of line/connect in state vector
      INTEGER(IntKi)                                     :: Iend    ! end index of line/connect in state vector


      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg = ""

      ! allocations of dxdt (as in SubDyn.  "INTENT(OUT) automatically deallocates the arrays on entry, we have to allocate them here"  is this right/efficient?)
      ALLOCATE ( dxdt%states(size(x%states)), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
          ErrMsg  = ' Error allocating dxdt%states array.'
          RETURN
      END IF


      ! clear connection force and mass values
      DO L = 1, p%NConnects
        DO J = 1,3
          m%ConnectList(L)%Ftot(J) = 0.0_DbKi
          m%ConnectList(L)%Ftot(J) = 0.0_DbKi
          DO K = 1,3
            m%ConnectList(L)%Mtot(K,J) = 0.0_DbKi
            m%ConnectList(L)%Mtot(K,J) = 0.0_DbKi
          END DO
        END DO
      END DO
      
      ! update fairlead positions for instantaneous values (fixed 2015-06-22)    
      DO K = 1, p%NFairs
         DO J = 1,3
            m%ConnectList(m%FairIdList(K))%r(J)  = u%PtFairleadDisplacement%Position(J,K) + u%PtFairleadDisplacement%TranslationDisp(J,K)
            m%ConnectList(m%FairIdList(K))%rd(J) = u%PtFairleadDisplacement%TranslationVel(J,K)  ! is this right? <<<
         END DO
      END DO

      ! apply line length changes from active tensioning if applicable
      DO L = 1, p%NLines
         IF (m%LineList(L)%CtrlChan > 0) then
            
            ! do a bounds check to prohibit excessive segment length changes (until a method to add/remove segments is created)
            IF ( u%DeltaL(m%LineList(L)%CtrlChan) > m%LineList(L)%UnstrLen / m%LineList(L)%N ) then
                ErrStat = ErrID_Fatal
                ErrMsg  = ' Active tension command will make a segment longer than the limit of twice its original length.'
                print *, u%DeltaL(m%LineList(L)%CtrlChan), " is an increase of more than ", (m%LineList(L)%UnstrLen / m%LineList(L)%N)
                print *, u%DeltaL
                print*, m%LineList(L)%CtrlChan
                RETURN
            END IF
            IF ( u%DeltaL(m%LineList(L)%CtrlChan) < -0.5 * m%LineList(L)%UnstrLen / m%LineList(L)%N ) then
             ErrStat = ErrID_Fatal
                ErrMsg  = ' Active tension command will make a segment shorter than the limit of half its original length.'
                print *, u%DeltaL(m%LineList(L)%CtrlChan), " is a reduction of more than half of ", (m%LineList(L)%UnstrLen / m%LineList(L)%N)
                print *, u%DeltaL
                print*, m%LineList(L)%CtrlChan
                RETURN
            END IF                
            
            ! for now this approach only acts on the fairlead end segment, and assumes all segment lengths are otherwise equal size
            m%LineList(L)%l( m%LineList(L)%N) = m%LineList(L)%UnstrLen/m%LineList(L)%N + u%DeltaL(m%LineList(L)%CtrlChan)       
            m%LineList(L)%ld(m%LineList(L)%N) =                                       u%DeltaLdot(m%LineList(L)%CtrlChan)       
         END IF
      END DO      

      ! do Line force and acceleration calculations, also add end masses/forces to respective Connects
      DO L = 1, p%NLines
        Istart = m%LineStateIndList(L)
        Iend = Istart + 6*(m%LineList(L)%N - 1) - 1
        CALL DoLineRHS(x%states(Istart:Iend), dxdt%states(Istart:Iend), t, m%LineList(L), &
          m%LineTypeList(m%LineList(L)%PropsIdNum), &
          m%ConnectList(m%LineList(L)%FairConnect)%Ftot, m%ConnectList(m%LineList(L)%FairConnect)%Mtot, &
          m%ConnectList(m%LineList(L)%AnchConnect)%Ftot, m%ConnectList(m%LineList(L)%AnchConnect)%Mtot )
      END DO


      ! perform connection force and mass calculations (done to all connects for sake of calculating fairlead/anchor loads)
      DO L = 1, p%NConnects
         ! add Connect's own forces including buoyancy and weight
         m%ConnectList(L)%Ftot(1) =m%ConnectList(L)%Ftot(1) + m%ConnectList(L)%conFX
         m%ConnectList(L)%Ftot(2) =m%ConnectList(L)%Ftot(2) + m%ConnectList(L)%conFY
         m%ConnectList(L)%Ftot(3) =m%ConnectList(L)%Ftot(3) + m%ConnectList(L)%conFZ + m%ConnectList(L)%conV*p%rhoW*p%g - m%ConnectList(L)%conM*p%g

         ! add Connect's own mass
         DO J = 1,3
            m%ConnectList(L)%Mtot(J,J) = m%ConnectList(L)%Mtot(J,J) + m%ConnectList(L)%conM
         END DO
      END DO  ! L


      ! do Connect acceleration calculations - changed to do only connect types
      DO L = 1, p%NConns
        Istart = L*6-5
        Iend = L*6
        CALL DoConnectRHS(x%states(Istart:Iend), dxdt%states(Istart:Iend), t, m%ConnectList(m%ConnIDList(L)))
      END DO


   CONTAINS


      !======================================================================
      SUBROUTINE DoLineRHS (X, Xd, t, Line, LineProp, FairFtot, FairMtot, AnchFtot, AnchMtot)

         Real(DbKi), INTENT( IN )      :: X(:)           ! state vector, provided
         Real(DbKi), INTENT( INOUT )   :: Xd(:)          ! derivative of state vector, returned ! cahnged to INOUT
         Real(DbKi), INTENT (IN)       :: t              ! instantaneous time
         TYPE(MD_Line), INTENT (INOUT) :: Line           ! label for the current line, for convenience
         TYPE(MD_LineProp), INTENT(IN) :: LineProp       ! the single line property set for the line of interest
         Real(DbKi), INTENT(INOUT)     :: FairFtot(:)    ! total force on Connect top of line is attached to
         Real(DbKi), INTENT(INOUT)     :: FairMtot(:,:)  ! total mass of Connect top of line is attached to
         Real(DbKi), INTENT(INOUT)     :: AnchFtot(:)    ! total force on Connect bottom of line is attached to
         Real(DbKi), INTENT(INOUT)     :: AnchMtot(:,:)  ! total mass of Connect bottom of line is attached to


         INTEGER(IntKi)                :: I              ! index of segments or nodes along line
         INTEGER(IntKi)                :: J              ! index
         INTEGER(IntKi)                :: K              ! index
         INTEGER(IntKi)                :: N              ! number of segments in line
         Real(DbKi)                    :: d              ! line diameter
         Real(DbKi)                    :: rho            ! line material density [kg/m^3]
         Real(DbKi)                    :: Sum1           ! for summing squares
         Real(DbKi)                    :: m_i            ! node mass
         Real(DbKi)                    :: v_i            ! node submerged volume
         Real(DbKi)                    :: Vi(3)          ! relative water velocity at a given node
         Real(DbKi)                    :: Vp(3)          ! transverse relative water velocity component at a given node
         Real(DbKi)                    :: Vq(3)          ! tangential relative water velocity component at a given node
         Real(DbKi)                    :: SumSqVp        !
         Real(DbKi)                    :: SumSqVq        !
         Real(DbKi)                    :: MagVp          !
         Real(DbKi)                    :: MagVq          !


         N = Line%N                      ! for convenience
         d = LineProp%d                  ! for convenience
         rho = LineProp%w/(Pi/4.0*d*d)



         ! set end node positions and velocities from connect objects' states
         DO J = 1, 3
            Line%r( J,N) = m%ConnectList(Line%FairConnect)%r(J)
            Line%r( J,0) = m%ConnectList(Line%AnchConnect)%r(J)
            Line%rd(J,N) = m%ConnectList(Line%FairConnect)%rd(J)
            Line%rd(J,0) = m%ConnectList(Line%AnchConnect)%rd(J)
         END DO

         ! set interior node positions and velocities
         DO I = 1, N-1
            DO J = 1, 3
               Line%r( J,I) = X( 3*N-3 + 3*I-3 + J)      ! r(J,I)  = X[3*N-3 + 3*i-3 + J]; // get positions  .. used to start from m%LineStateIndList(Line%IdNum) in whole state vector
               Line%rd(J,I) = X(         3*I-3 + J)      ! rd(J,I) = X[        3*i-3 + J]; // get velocities
            END DO
         END DO

         ! calculate instantaneous (stretched) segment lengths and rates << should add catch here for if lstr is ever zero
         DO I = 1, N
            Sum1 = 0.0_DbKi
            DO J = 1, 3
               Sum1 = Sum1 + (Line%r(J,I) - Line%r(J,I-1)) * (Line%r(J,I) - Line%r(J,I-1))
            END DO
            Line%lstr(I) = sqrt(Sum1)                                  ! stretched segment length

            Sum1 = 0.0_DbKi
            DO J = 1, 3
               Sum1 = Sum1 + (Line%r(J,I) - Line%r(J,I-1))*(Line%rd(J,I) - Line%rd(J,I-1))
            END DO
            Line%lstrd(I) = Sum1/Line%lstr(I)                          ! segment stretched length rate of change

    !       Line%V(I) = Pi/4.0 * d*d*Line%l(I)                        !volume attributed to segment
         END DO

         !calculate unit tangent vectors (q) for each node (including ends) note: I think these are pointing toward 0 rather than N!
         CALL UnitVector(Line%q(:,0), Line%r(:,1), Line%r(:,0)) ! compute unit vector q
         DO I = 1, N-1
           CALL UnitVector(Line%q(:,I), Line%r(:,I+1), Line%r(:,I-1)) ! compute unit vector q ... using adjacent two nodes!
         END DO
         CALL UnitVector(Line%q(:,N), Line%r(:,N), Line%r(:,N-1))     ! compute unit vector q


         ! wave kinematics not implemented yet


         !calculate mass (including added mass) matrix for each node
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
                     Line%M(K,J,I) = m_i + p%rhoW*v_i*( LineProp%Can*(1 - Line%q(J,I)*Line%q(K,I)) + LineProp%Cat*Line%q(J,I)*Line%q(K,I) )
                  ELSE
                     Line%M(K,J,I) = p%rhoW*v_i*( LineProp%Can*(-Line%q(J,I)*Line%q(K,I)) + LineProp%Cat*Line%q(J,I)*Line%q(K,I) )
                  END IF
               END DO
            END DO

            CALL Inverse3by3(Line%S(:,:,I), Line%M(:,:,I))             ! invert mass matrix
         END DO


         ! ------------------  CALCULATE FORCES ON EACH NODE ----------------------------

         ! loop through the segments
         DO I = 1, N

            ! line tension, inherently including possibility of dynamic length changes in l term
            IF (Line%lstr(I)/Line%l(I) > 1.0) THEN
               DO J = 1, 3
                  Line%T(J,I) = LineProp%EA *( 1.0/Line%l(I) - 1.0/Line%lstr(I) ) * (Line%r(J,I)-Line%r(J,I-1))
               END DO
            ELSE
               DO J = 1, 3
                  Line%T(J,I) = 0.0_DbKi                              ! cable can't "push"
               END DO
            END if

            ! line internal damping force based on line-specific BA value, including possibility of dynamic length changes in l and ld terms
            DO J = 1, 3
               Line%Td(J,I) = Line%BA* ( Line%lstrd(I) -  Line%lstr(I)*Line%ld(I)/Line%l(I) )/Line%l(I)  * (Line%r(J,I)-Line%r(J,I-1)) / Line%lstr(I)
            END DO
         END DO



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

            !relative flow velocities
            DO J = 1, 3
               Vi(J) = 0.0 - Line%rd(J,I)                               ! relative flow velocity over node -- this is where wave velicites would be added
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
               DO J = 1, 3
                  Line%Dp(J,I) = 0.25*p%rhoW*LineProp%Cdn*    d*Line%l(1) * MagVp * Vp(J)
                  Line%Dq(J,I) = 0.25*p%rhoW*LineProp%Cdt* Pi*d*Line%l(1) * MagVq * Vq(J)
               END DO
            ELSE IF (I==N)  THEN
               DO J = 1, 3
                  Line%Dp(J,I) = 0.25*p%rhoW*LineProp%Cdn*    d*Line%l(N) * MagVp * Vp(J);
                  Line%Dq(J,I) = 0.25*p%rhoW*LineProp%Cdt* Pi*d*Line%l(N) * MagVq * Vq(J)
               END DO
            ELSE
               DO J = 1, 3
                  Line%Dp(J,I) = 0.25*p%rhoW*LineProp%Cdn*    d*(Line%l(I) + Line%l(I+1)) * MagVp * vp(J);
                  Line%Dq(J,I) = 0.25*p%rhoW*LineProp%Cdt* Pi*d*(Line%l(I) + Line%l(I+1)) * MagVq * vq(J);
               END DO
            END IF

            ! F-K force from fluid acceleration not implemented yet

            ! bottom contact (stiffness and damping, vertical-only for now)  - updated Nov 24 for general case where anchor and fairlead ends may deal with bottom contact forces

            IF (Line%r(3,I) < -p%WtrDpth) THEN
               IF (I==0) THEN
                  Line%B(3,I) = ( (-p%WtrDpth - Line%r(3,I))*p%kBot - Line%rd(3,I)*p%cBot) * 0.5*d*(            Line%l(I+1) ) 
               ELSE IF (I==N) THEN
                  Line%B(3,I) = ( (-p%WtrDpth - Line%r(3,I))*p%kBot - Line%rd(3,I)*p%cBot) * 0.5*d*(Line%l(I)               ) 
               ELSE
                  Line%B(3,I) = ( (-p%WtrDpth - Line%r(3,I))*p%kBot - Line%rd(3,I)*p%cBot) * 0.5*d*(Line%l(I) + Line%l(I+1) ) 



               END IF
            ELSE
               Line%B(3,I) = 0.0_DbKi
            END IF

            ! total forces
            IF (I==0)  THEN
               DO J = 1, 3
                  Line%F(J,I) = Line%T(J,1)                 + Line%Td(J,1)                  + Line%W(J,I) + Line%Dp(J,I) + Line%Dq(J,I) + Line%B(J,I)
               END DO
            ELSE IF (I==N)  THEN
               DO J = 1, 3
                  Line%F(J,I) =                -Line%T(J,N)                  - Line%Td(J,N) + Line%W(J,I) + Line%Dp(J,I) + Line%Dq(J,I) + Line%B(J,I)
               END DO
            ELSE
               DO J = 1, 3
                  Line%F(J,I) = Line%T(J,I+1) - Line%T(J,I) + Line%Td(J,I+1) - Line%Td(J,I) + Line%W(J,I) + Line%Dp(J,I) + Line%Dq(J,I) + Line%B(J,I)
               END DO
            END IF

         END DO  ! I  - done looping through nodes


         ! loop through internal nodes and update their states
         DO I=1, N-1
            DO J=1,3

               ! calculate RHS constant (premultiplying force vector by inverse of mass matrix  ... i.e. rhs = S*Forces)
               Sum1 = 0.0_DbKi                               ! reset temporary accumulator
               DO K = 1, 3
                 Sum1 = Sum1 + Line%S(K,J,I) * Line%F(K,I)   ! matrix-vector multiplication [S i]{Forces i}  << double check indices
               END DO ! K

               ! update states
               Xd(3*N-3 + 3*I-3 + J) = X(3*I-3 + J);         ! dxdt = V  (velocities)
               Xd(        3*I-3 + J) = Sum1                  ! dVdt = RHS * A  (accelerations)

            END DO ! J
         END DO  ! I


         ! add force and mass of end nodes to the Connects they correspond to
         DO J = 1,3
            FairFtot(J) = FairFtot(J) + Line%F(J,N)
            AnchFtot(J) = AnchFtot(J) + Line%F(J,0)
            DO K = 1,3
               FairMtot(K,J) = FairMtot(K,J) + Line%M(K,J,N)
               AnchMtot(K,J) = AnchMtot(K,J) + Line%M(K,J,0)
            END DO
         END DO

      END SUBROUTINE DoLineRHS
      !=====================================================================


      !======================================================================
      SUBROUTINE DoConnectRHS (X, Xd, t, Connect)

         ! This subroutine is for the "Connect" type of Connections only.  Other types don't have their own state variables.
      
         Real(DbKi),       INTENT( IN )    :: X(:)           ! state vector for this connect, provided
         Real(DbKi),       INTENT( OUT )   :: Xd(:)          ! derivative of state vector for this connect, returned
         Real(DbKi),       INTENT (IN)     :: t              ! instantaneous time
         Type(MD_Connect), INTENT (INOUT)  :: Connect        ! Connect number


         !INTEGER(IntKi)             :: I         ! index of segments or nodes along line
         INTEGER(IntKi)             :: J         ! index
         INTEGER(IntKi)             :: K         ! index
         Real(DbKi)                 :: Sum1      ! for adding things

         ! When this sub is called, the force and mass contributions from the attached Lines should already have been added to
         ! Fto and Mtot by the Line RHS function.  Also, any self weight, buoyancy, or external forcing should have already been
         ! added by the calling subroutine.  The only thing left is any added mass or drag forces from the connection (e.g. float)
         ! itself, which will be added below.


         IF (EqualRealNos(t, 0.0_DbKi)) THEN  ! this is old: with current IC gen approach, we skip the first call to the line objects, because they're set AFTER the call to the connects

            DO J = 1,3
               Xd(3+J) = X(J)        ! velocities - these are unused in integration
               Xd(J) = 0.0_DbKi           ! accelerations - these are unused in integration
            END DO
         ELSE
            ! from state values, get r and rdot values
            DO J = 1,3
               Connect%r(J)  = X(3 + J)   ! get positions
               Connect%rd(J) = X(J)       ! get velocities
            END DO
         END IF
         

         ! add any added mass and drag forces from the Connect body itself
         DO J = 1,3
            Connect%Ftot(J) = Connect%Ftot(J) - 0.5 * p%rhoW * Connect%rd(J) * abs(Connect%rd(J)) * Connect%conCdA;  ! add drag forces - corrected Nov 24
            Connect%Mtot(J,J) = Connect%Mtot(J,J) + Connect%conV*p%rhoW*Connect%conCa;                               ! add added mass
         END DO
                     
         ! invert node mass matrix
         CALL Inverse3by3(Connect%S, Connect%Mtot)

         DO J = 1,3
            ! RHS constant - (premultiplying force vector by inverse of mass matrix  ... i.e. rhs = S*Forces
            Sum1 = 0.0_DbKi   ! reset accumulator
            DO K = 1, 3
               Sum1 = Sum1 + Connect%S(K,J) * Connect%Ftot(K)   !  matrix multiplication [S i]{Forces i}
            END DO

            ! update states
            Xd(3 + J) = X(J)          ! dxdt = V    (velocities)
            Xd(J) = Sum1              ! dVdt = RHS * A  (accelerations)
         END DO

      END SUBROUTINE DoConnectRHS
      !=====================================================================

   END SUBROUTINE MD_CalcContStateDeriv
   !=============================================================================================


   !===============================================================================================
   SUBROUTINE MD_End(u, p, x, xd, z, other, y, m, ErrStat , ErrMsg)
      TYPE(MD_InputType) ,            INTENT(INOUT) :: u
      TYPE(MD_ParameterType) ,        INTENT(INOUT) :: p
      TYPE(MD_ContinuousStateType) ,  INTENT(INOUT) :: x
      TYPE(MD_DiscreteStateType) ,    INTENT(INOUT) :: xd
      TYPE(MD_ConstraintStateType) ,  INTENT(INOUT) :: z
      TYPE(MD_OtherStateType) ,       INTENT(INOUT) :: other
      TYPE(MD_OutputType) ,           INTENT(INOUT) :: y
      TYPE(MD_MiscVarType),           INTENT(INOUT) :: m      
      INTEGER(IntKi),                 INTENT(  OUT) :: ErrStat
      CHARACTER(*),                   INTENT(  OUT) :: ErrMsg

!      INTEGER(IntKi)                                :: i=0

      INTEGER(IntKi)                               :: ErrStat2      ! Error status of the operation
      CHARACTER(ErrMsgLen)                         :: ErrMsg2       ! Error message if ErrStat2 /= ErrID_None

      ErrStat = ErrID_None
      ErrMsg  = ""



      ! deallocate data associated with file output
      CALL MDIO_CloseOutput ( p, m, ErrStat2, ErrMsg2 )
         CALL CheckError( ErrStat2, ErrMsg2 )
         !IF (ErrStat >= AbortErrLev) RETURN


      ! deallocate FAST data structures
      CALL MD_DestroyInput(u, ErrStat2, ErrMsg2)
         CALL CheckError( ErrStat2, ErrMsg2 )
      CALL MD_DestroyParam(p, ErrStat2, ErrMsg2)
         CALL CheckError( ErrStat2, ErrMsg2 )
      CALL MD_DestroyContState(x, ErrStat2, ErrMsg2)  ! <--- getting access violation
         CALL CheckError( ErrStat2, ErrMsg2 )
      CALL MD_DestroyDiscState(xd, ErrStat2, ErrMsg2)
         CALL CheckError( ErrStat2, ErrMsg2 )
      CALL MD_DestroyConstrState(z, ErrStat2, ErrMsg2)
         CALL CheckError( ErrStat2, ErrMsg2 )
      CALL MD_DestroyOtherState(other,ErrStat2,ErrMsg2) ! <--- getting access violation
         CALL CheckError( ErrStat2, ErrMsg2 )
      CALL MD_DestroyOutput(y, ErrStat2, ErrMsg2)
         CALL CheckError( ErrStat2, ErrMsg2 )
      CALL MD_DestroyMisc(m, ErrStat2, ErrMsg2)
         CALL CheckError( ErrStat2, ErrMsg2 )

 !     IF ( ErrStat==ErrID_None) THEN
 !        CALL WrScr('MoorDyn closed without errors')
 !     ELSE
 !        CALL WrScr('MoorDyn closed with errors')
 !     END IF


   CONTAINS

      SUBROUTINE CheckError(ErrId, Msg)
        ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev

         INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
         CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)


         IF ( ErrID /= ErrID_None ) THEN

            IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine  ! keep existing error message if there is one
            ErrMsg = TRIM(ErrMsg)//' MD_End:'//TRIM(Msg)      ! add current error message
            ErrStat = MAX(ErrStat, ErrID)

            CALL WrScr( ErrMsg )  ! do this always or only if warning level?

         END IF

      END SUBROUTINE CheckError


   END SUBROUTINE MD_End                                                                         !   -------+
   !==========================================================================================================


!!==========   MD_CheckError   =======     <---------------------------------------------------------------+
! SUBROUTINE MD_CheckError(InMsg,OutMsg)
!   ! Passed arguments
!!   CHARACTER(*), INTENT(IN   ) :: InMsg       ! The input string
!   CHARACTER(*), INTENT(INOUT) :: OutMsg      ! The error message (ErrMsg)!
!
 !  OutMsg = InMsg
 !  RETURN
 !END SUBROUTINE MD_CheckError                                                                  !   -------+
 !==========================================================================================================




   !========================================================================================================
   SUBROUTINE TimeStep ( t, dtStep, u, utimes, p, x, xd, z, other, m, ErrStat, ErrMsg )
      REAL(DbKi)                     , INTENT(INOUT)      :: t
      REAL(ReKi)                     , INTENT(IN   )      :: dtStep     ! how long to advance the time for
      TYPE( MD_InputType )           , INTENT(INOUT)      :: u(:)       ! INTENT(IN   )
      REAL(DbKi)                     , INTENT(IN   )      :: utimes(:)  ! times corresponding to elements of u(:)?
      TYPE( MD_ParameterType )       , INTENT(IN   )      :: p          ! INTENT(IN   )
      TYPE( MD_ContinuousStateType ) , INTENT(INOUT)      :: x
      TYPE( MD_DiscreteStateType )   , INTENT(IN   )      :: xd         ! INTENT(IN   )
      TYPE( MD_ConstraintStateType ) , INTENT(IN   )      :: z          ! INTENT(IN   )
      TYPE( MD_OtherStateType )      , INTENT(IN   )      :: other      ! INTENT(INOUT)
      TYPE(MD_MiscVarType)           , INTENT(INOUT)      :: m          ! INTENT(INOUT)
      INTEGER(IntKi)                 , INTENT(  OUT)      :: ErrStat
      CHARACTER(*)                   , INTENT(  OUT)      :: ErrMsg


      TYPE(MD_ContinuousStateType)                        :: dxdt       ! time derivatives of continuous states (initialized in CalcContStateDeriv)
      TYPE(MD_ContinuousStateType)                        :: x2         ! temporary copy of continuous states used in RK2 calculations
      INTEGER(IntKi)                                      :: NdtM       ! the number of time steps to make with the mooring model
      Real(DbKi)                                          :: dtM        ! the actual time step size to use
      INTEGER(IntKi)                                      :: Nx         ! size of states vector
      INTEGER(IntKi)                                      :: I          ! counter
      INTEGER(IntKi)                                      :: J          ! counter
      TYPE(MD_InputType)                                  :: u_interp   ! interpolated instantaneous input values to be calculated for each mooring time step

      Real(DbKi)                                          :: tDbKi   ! double version because that's what MD_Input_ExtrapInterp needs.
      
      
      ! allocate space for x2
      CALL MD_CopyContState( x, x2, 0, ErrStat, ErrMsg)
         
      ! create space for arrays/meshes in u_interp   ... is it efficient to do this every time step???
      CALL MD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg)
         

      Nx = size(x%states)


      ! round dt to integer number of time steps
      NdtM = ceiling(dtStep/p%dtM0)                  ! get number of mooring time steps to do based on desired time step size
      dtM = dtStep/float(NdtM)                       ! adjust desired time step to satisfy dt with an integer number of time steps


      !loop through line integration time steps
      DO I = 1, NdtM                                 ! for (double ts=t; ts<=t+ICdt-dts; ts+=dts)
      
      
         !tDbKi = t        ! get DbKi version of current time (why does ExtrapInterp except different time type than UpdateStates?)
         
      
         ! -------------------------------------------------------------------------------
         !       RK2 integrator written here, now calling CalcContStateDeriv
         !--------------------------------------------------------------------------------

         ! step 1

         CALL MD_Input_ExtrapInterp(u, utimes, u_interp, t          , ErrStat, ErrMsg)   ! interpolate input mesh to correct time (t)
      
         CALL MD_CalcContStateDeriv( t, u_interp, p, x, xd, z, other, m, dxdt, ErrStat, ErrMsg )
         DO J = 1, Nx
            x2%states(J) = x%states(J) + 0.5*dtM*dxdt%states(J)                                           !x1 = x0 + dt*f0/2.0;
         END DO

         ! step 2
   
         CALL MD_Input_ExtrapInterp(u, utimes, u_interp, t + 0.5_DbKi*dtM, ErrStat, ErrMsg)   ! interpolate input mesh to correct time (t+0.5*dtM)
            
         CALL MD_CalcContStateDeriv( (t + 0.5_DbKi*dtM), u_interp, p, x2, xd, z, other, m, dxdt, ErrStat, ErrMsg )       !called with updated states x2 and time = t + dt/2.0
         DO J = 1, Nx
            x%states(J) = x%states(J) + dtM*dxdt%states(J)
         END DO

         t = t + dtM  ! update time
         
         !print *, " In TimeStep t=", t, ",  L1N8Pz=", M%LineList(1)%r(3,8), ", dL1=", u_interp%DeltaL(1)

         !----------------------------------------------------------------------------------

   ! >>> below should no longer be necessary thanks to using ExtrapInterp of u(:) within the mooring time stepping loop.. <<<
   !      ! update Fairlead positions by integrating velocity and last position (do this AFTER the processing of the time step rather than before)
   !      DO J = 1, p%NFairs
   !         DO K = 1, 3
   !          m%ConnectList(m%FairIdList(J))%r(K) = m%ConnectList(m%FairIdList(J))%r(K) + m%ConnectList(m%FairIdList(J))%rd(K)*dtM
   !         END DO
   !      END DO
      
   
      END DO  ! I  time steps


      ! destroy dxdt and x2, and u_interp
      CALL MD_DestroyContState( dxdt, ErrStat, ErrMsg)
      CALL MD_DestroyContState( x2, ErrStat, ErrMsg)
      CALL MD_DestroyInput(u_interp, ErrStat, ErrMsg)
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error destroying dxdt or x2.'
      END IF
      !   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MD_UpdateStates')

      
      ! check for NaNs - is this a good place/way to do it?
      DO J = 1, Nx
         IF (Is_NaN(REAL(x%states(J),DbKi))) THEN
            ErrStat = ErrID_Fatal
            ErrMsg = ' NaN state detected.'
         END IF
      END DO
 

   END SUBROUTINE TimeStep
   !======================================================================



   !=======================================================================
   SUBROUTINE SetupLine (Line, LineProp, rhoW, ErrStat, ErrMsg)
      ! allocate arrays in line object

      TYPE(MD_Line), INTENT(INOUT)       :: Line          ! the single line object of interest
      TYPE(MD_LineProp), INTENT(INOUT)   :: LineProp      ! the single line property set for the line of interest
      REAL(ReKi),    INTENT(IN)          :: rhoW
      INTEGER,       INTENT(   INOUT )   :: ErrStat       ! returns a non-zero value when an error occurs
      CHARACTER(*),  INTENT(   INOUT )   :: ErrMsg        ! Error message if ErrStat /= ErrID_None

      INTEGER(4)                         :: J             ! Generic index
      INTEGER(4)                         :: K             ! Generic index
      INTEGER(IntKi)                     :: N

      N = Line%N  ! number of segments in this line (for code readability)

      ! allocate node positions and velocities (NOTE: these arrays start at ZERO)
      ALLOCATE ( Line%r(3, 0:N), Line%rd(3, 0:N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating r and rd arrays.'
         !CALL CleanUp()
         RETURN
      END IF

      ! allocate node tangent vectors
      ALLOCATE ( Line%q(3, 0:N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating q array.'
         !CALL CleanUp()
         RETURN
      END IF

      ! allocate segment scalar quantities
      ALLOCATE ( Line%l(N), Line%ld(N), Line%lstr(N), Line%lstrd(N), Line%V(N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating segment scalar quantity arrays.'
         !CALL CleanUp()
         RETURN
      END IF

      ! assign values for l and V
      DO J=1,N
         Line%l(J) = Line%UnstrLen/REAL(N, DbKi)
         Line%ld(J)= 0.0_DbKi
         Line%V(J) = Line%l(J)*0.25*Pi*LineProp%d*LineProp%d
      END DO

      ! allocate segment tension and internal damping force vectors
      ALLOCATE ( Line%T(3, N), Line%Td(3, N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating T and Td arrays.'
         !CALL CleanUp()
         RETURN
      END IF

      ! allocate node force vectors
      ALLOCATE ( Line%W(3, 0:N), Line%Dp(3, 0:N), Line%Dq(3, 0:N), Line%Ap(3, 0:N), &
         Line%Aq(3, 0:N), Line%B(3, 0:N), Line%F(3, 0:N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating node force arrays.'
         !CALL CleanUp()
         RETURN
      END IF

      ! set gravity and bottom contact forces to zero initially (because the horizontal components should remain at zero)
      DO J = 0,N
         DO K = 1,3
            Line%W(K,J) = 0.0_DbKi
            Line%B(K,J) = 0.0_DbKi
         END DO
      END DO

      ! allocate mass and inverse mass matrices for each node (including ends)
      ALLOCATE ( Line%S(3, 3, 0:N), Line%M(3, 3, 0:N), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating T and Td arrays.'
         !CALL CleanUp()
         RETURN
      END IF
      
      ! Specify specific internal damping coefficient (BA) for this line.
      ! Will be equal to inputted BA of LineType if input value is positive.
      ! If input value is negative, it is considered to be desired damping ratio (zeta)
      ! from which the line's BA can be calculated based on the segment natural frequency.
      IF (LineProp%BA < 0) THEN
         ! - we assume desired damping coefficient is zeta = -LineProp%BA
         ! - highest axial vibration mode of a segment is wn = sqrt(k/m) = 2N/UnstrLen*sqrt(EA/w)
         Line%BA = -LineProp%BA * Line%UnstrLen / Line%N * SQRT(LineProp%EA * LineProp%w)
      !  print *, 'Based on zeta, BA set to ', Line%BA
         
      !  print *, 'Negative BA input detected, treating as -zeta.  For zeta = ', -LineProp%BA, ', setting BA to ', Line%BA
         
      ELSE
         Line%BA = LineProp%BA
      !  temp = Line%N * Line%BA / Line%UnstrLen * SQRT(1.0/(LineProp%EA * LineProp%w))
      !  print *, 'BA set as input to ', Line%BA, '. Corresponding zeta is ', temp
      END IF
      
      !temp = 2*Line%N / Line%UnstrLen * sqrt( LineProp%EA / LineProp%w) / TwoPi
      !print *, 'Segment natural frequency is ', temp, ' Hz'
      
      
      ! need to add cleanup sub <<<


   END SUBROUTINE SetupLine
   !======================================================================




   !===============================================================================================
   SUBROUTINE InitializeLine (Line, LineProp, rhoW, ErrStat, ErrMsg)
      ! calculate initial profile of the line using quasi-static model

      TYPE(MD_Line),     INTENT(INOUT)       :: Line        ! the single line object of interest
      TYPE(MD_LineProp), INTENT(INOUT)       :: LineProp    ! the single line property set for the line of interest
      REAL(ReKi),        INTENT(IN)          :: rhoW
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

          DO J = 0,Line%N ! Loop through all nodes per line where the line position and tension can be output
             Line%r(1,J) = Line%r(1,0) + LNodesX(J+1)*COSPhi
             Line%r(2,J) = Line%r(2,0) + LNodesX(J+1)*SINPhi
             Line%r(3,J) = Line%r(3,0) + LNodesZ(J+1)
          ENDDO              ! J - All nodes per line where the line position and tension can be output


      ELSE ! if there is a problem with the catenary approach, just stretch the nodes linearly between fairlead and anchor

          CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'InitializeLine')

          DO J = 0,Line%N ! Loop through all nodes per line where the line position and tension can be output
             Line%r(1,J) = Line%r(1,0) + (Line%r(1,N) - Line%r(1,0))*REAL(J, DbKi)/REAL(N, DbKi)
             Line%r(2,J) = Line%r(2,0) + (Line%r(2,N) - Line%r(2,0))*REAL(J, DbKi)/REAL(N, DbKi)
             Line%r(3,J) = Line%r(3,0) + (Line%r(3,N) - Line%r(3,0))*REAL(J, DbKi)/REAL(N, DbKi)
          ENDDO
      ENDIF



      CALL CleanUp()  ! deallocate temporary arrays



   CONTAINS


      !=======================================================================
      SUBROUTINE CleanUp()
           ! deallocate temporary arrays

           IF (ALLOCATED(LSNodes))  DEALLOCATE(LSNodes)
           IF (ALLOCATED(LNodesX))  DEALLOCATE(LNodesX)
           IF (ALLOCATED(LNodesZ))  DEALLOCATE(LNodesZ)

        END SUBROUTINE CleanUp
      !=======================================================================


      !=======================================================================
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
           ErrMsg =  ' A fairlead has passed below its anchor.'
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
               ErrMsg =  ' Unstretched mooring line length too large. '// &
                            ' Routine Catenary() cannot solve quasi-static mooring line solution.'
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
      !=======================================================================


   END SUBROUTINE InitializeLine
   !======================================================================



   ! ============ below are some math convenience functions ===============
   ! should add error checking if I keep these, but hopefully there are existing NWTCLib functions to replace them


   ! return unit vector (u) in direction from r1 to r2
   !=======================================================================
   SUBROUTINE UnitVector( u, r1, r2 )
      REAL(DbKi), INTENT(OUT)   :: u(:)
      REAL(DbKi), INTENT(IN)    :: r1(:)
      REAL(DbKi), INTENT(IN)    :: r2(:)

      REAL(DbKi)                :: Length

      u = r2 - r1
      Length = TwoNorm(u)

      if ( .NOT. EqualRealNos(length, 0.0_DbKi ) ) THEN
        u = u / Length
      END IF

   END SUBROUTINE UnitVector
   !=======================================================================


   !compute the inverse of a 3-by-3 matrix m
   !=======================================================================
   SUBROUTINE Inverse3by3( Minv, M )
      Real(DbKi), INTENT(OUT)   :: Minv(:,:)  ! returned inverse matrix
      Real(DbKi), INTENT(IN)    :: M(:,:)     ! inputted matrix

      Real(DbKi)                :: det        ! the determinant
      Real(DbKi)                :: invdet     ! inverse of the determinant

      det = M(1, 1) * (M(2, 2) * M(3, 3) - M(3, 2) * M(2, 3)) - &
            M(1, 2) * (M(2, 1) * M(3, 3) - M(2, 3) * M(3, 1)) + &
            M(1, 3) * (M(2, 1) * M(3, 2) - M(2, 2) * M(3, 1));

      invdet = 1.0 / det   ! because multiplying is faster than dividing

      Minv(1, 1) = (M(2, 2) * M(3, 3) - M(3, 2) * M(2, 3)) * invdet
      Minv(1, 2) = (M(1, 3) * M(3, 2) - M(1, 2) * M(3, 3)) * invdet
      Minv(1, 3) = (M(1, 2) * M(2, 3) - M(1, 3) * M(2, 2)) * invdet
      Minv(2, 1) = (M(2, 3) * M(3, 1) - M(2, 1) * M(3, 3)) * invdet
      Minv(2, 2) = (M(1, 1) * M(3, 3) - M(1, 3) * M(3, 1)) * invdet
      Minv(2, 3) = (M(2, 1) * M(1, 3) - M(1, 1) * M(2, 3)) * invdet
      Minv(3, 1) = (M(2, 1) * M(3, 2) - M(3, 1) * M(2, 2)) * invdet
      Minv(3, 2) = (M(3, 1) * M(1, 2) - M(1, 1) * M(3, 2)) * invdet
      Minv(3, 3) = (M(1, 1) * M(2, 2) - M(2, 1) * M(1, 2)) * invdet

   END SUBROUTINE Inverse3by3
   !=======================================================================


END MODULE MoorDyn
