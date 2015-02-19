MODULE MoorDyn

   USE MoorDyn_Types
   USE MoorDyn_IO
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: MD_ProgDesc = ProgDesc( 'MoorDyn', 'v0.9.01-mth', '10-Dec-2014' )


   PUBLIC :: MD_Init
   PUBLIC :: MD_UpdateStates
   PUBLIC :: MD_CalcOutput
   PUBLIC :: MD_End

CONTAINS

   !==========   MD_Init   ======     <----------------------------------------------------------------------+
   SUBROUTINE MD_Init(InitInp, u, p, x, xd, z, other, y, DTcoupling, InitOut, ErrStat, ErrMsg)
      IMPLICIT NONE
      TYPE(MD_InitInputType),       INTENT(INOUT)  :: InitInp     ! INTENT(INOUT) : Input data for initialization routine
      TYPE(MD_InputType),           INTENT(  OUT)  :: u           ! INTENT( OUT) : An initial guess for the input; input mesh must be defined
      TYPE(MD_ParameterType),       INTENT(  OUT)  :: p           ! INTENT( OUT) : Parameters
      TYPE(MD_ContinuousStateType), INTENT(  OUT)  :: x           ! INTENT( OUT) : Initial continuous states
      TYPE(MD_DiscreteStateType),   INTENT(  OUT)  :: xd          ! INTENT( OUT) : Initial discrete states
      TYPE(MD_ConstraintStateType), INTENT(  OUT)  :: z           ! INTENT( OUT) : Initial guess of the constraint states
      TYPE(MD_OtherStateType),      INTENT(  OUT)  :: other       ! INTENT( OUT) : Initial other/optimization states
      TYPE(MD_OutputType),          INTENT(  OUT)  :: y           ! INTENT( OUT) : Initial system outputs (outputs are not calculated; only the output mesh is initialized)
      REAL(DbKi),                   INTENT(INOUT)  :: DTcoupling  ! Coupling interval in seconds: the rate that Output is the actual coupling interval
      TYPE(MD_InitOutputType),      INTENT(INOUT)  :: InitOut     ! Output for initialization routine
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
      REAL(ReKi)                                   :: t      ! instantaneous time, to be used during IC generation
      INTEGER(IntKi)                               :: I      ! index
      INTEGER(IntKi)                               :: J      ! index
      INTEGER(IntKi)                               :: K      ! index
      INTEGER(IntKi)                               :: N      ! convenience integer for readability: number of segments in the line
      REAL(ReKi)                                   :: Pos(3)      ! temporary array for setting absolute fairlead positions in mesh
      REAL(ReKi), ALLOCATABLE                      :: FairTensIC(:,:)      ! array of size Nfairs, 3 to store three latest fairlead tensions of each line
      CHARACTER(20)                     :: TempString       ! temporary string for incidental use
      !    INTEGER(IntKi)                               :: NumNodes = 0
      INTEGER(IntKi)                               :: ErrStat2      ! Error status of the operation
      CHARACTER(LEN(ErrMsg))                       :: ErrMsg2       ! Error message if ErrStat /= ErrID_None

      ErrStat = ErrID_None
      ErrMsg  = ""

      ! Initialize the NWTC Subroutine Library
      CALL NWTC_Init( )

      ! Display the module information
      CALL DispNVD( MD_ProgDesc )
      InitOut%Ver = MD_ProgDesc



      ! -------------------------------------------------------------------------------------------

      ! call function that reads input file and creates cross-referenced Connect and Line objects

      ! -------------------------------------------------------------------------------------------
      CALL MDIO_ReadInput(InitInp, p, other, ErrStat2, ErrMsg2)

      CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MD_Init' )         ! should I use this with all calls???
      IF ( ErrStat >= AbortErrLev ) THEN
         !CALL CleanUp()
         RETURN
      END IF


      ! process the OutList array and set up the index arrays for the requested output quantities
      CALL MDIO_SetOutParam(InitInp%OutList, p, other, y, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         RETURN
      END IF

      print *, 'SetOutParam is done.  Size of y%WriteOutput is ', size(y%WriteOutput)


      CALL WrScr( 'moving on to the next step in MD_Init' )

      !-------------------------------------------------------------------------------------------------
      !  Connect mooring system together and make necessary allocations
      !-------------------------------------------------------------------------------------------------

      ! cycle through Connects and identify Vessel types
      DO I = 1, p%NConnects
      TempString = other%ConnectList(I)%type
      CALL Conv2UC(TempString)
      if (TempString == 'FIXED') then
         other%ConnectList(I)%TypeNum = 0
      else if (TempString == 'VESSEL') then
         other%ConnectList(I)%TypeNum = 1
         p%NFairs = p%NFairs + 1             ! if a vessel connection, increment fairlead counter
      else if (TempString == 'CONNECT') then
         other%ConnectList(I)%TypeNum = 2
      else
         ErrStat = ErrID_Fatal
         ErrMsg = ' Error in provided Connect type.  Must be fixed, vessel, or connect.'
         !CALL CleanUp()
         RETURN
      END IF
      END DO

      ! allocate fairleads list
      ALLOCATE ( other%FairIdList(p%NFairs), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating space for FairIdList array.'
         !CALL CleanUp()
         RETURN
      END IF

      ! now go back through and record the fairlead Id numbers (this is all the "connecting" that's required)
      J = 1  ! counter for fairlead number
      DO I = 1,p%NConnects
         IF (other%ConnectList(I)%TypeNum == 1) THEN
           other%FairIdList(J) = I             ! if a vessel connection, add ID to fairlead list
           J = J + 1
         END IF
      END DO


      ! go through lines and allocate variables
      DO I = 1, p%NLines
        CALL SetupLine( other%LineList(I), other%LineTypeList(other%LineList(I)%PropsIdNum), p%rhoW ,  ErrStat, ErrMsg)
      END DO


      ! now prepare state vectors

      ! allocate list of starting state vector indices for each line  - does this belong elsewhere?
      ALLOCATE ( other%LineStateIndList(p%NLines), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
        ErrMsg  = ' Error allocating LineStateIndList array.'
        !CALL CleanUp()
        RETURN
      END IF

      ! J is keeping track of the growing size of the state vector
      J = p%NConnects*6   ! start index of first line's states (added six state variables for each connection) - could work out a way to avoid states for fairleads and anchors

      DO I = 1, p%NLines
         other%LineStateIndList(I) = J+1            ! assign start index of each line
         J = J + 6*(other%LineList(I)%N - 1)  !add 6 state variables for each internal node
      END DO

      ! allocate state vector and F vector for RK2 based on size just calculated
      ALLOCATE ( x%states(J), other%F(J), STAT = ErrStat )

      IF ( ErrStat /= ErrID_None ) THEN
        ErrMsg  = ' Error allocating state vector.'
        !CALL CleanUp()
        RETURN
      END IF

      CALL WrScr( 'Done allocating variables' )


    ! create meshes for fairlead positions and forces

   ! below is for fairlead i/o
    !==========   MD Mesh initialization   ======      <--------------------------+
    ! get header information for the FAST output file <-- matt: what does this mean?
    ! Create the input mesh                                            !          |
                                                                       !          |
    CALL MeshCreate(BlankMesh=u%PtFairleadDisplacement , &             !          |
                    IOS= COMPONENT_INPUT           , &                 !          |
                    NNodes=p%NFairs                , &                 !          |
                    TranslationDisp=.TRUE.         , &                 !          |
                    TranslationVel=.TRUE.          , &                 !          |
                    ErrStat=ErrStat                , &                 !          |
                    ErrMess=ErrMsg)                                    !          |
    IF(ErrStat/=ErrID_None) CALL WrScr(TRIM(ErrMsg))                   !          |
                                                                       !          |
    DO i = 1,p%NFairs

       ! set position of each node (this way rather than manually... for error checking?  seems more tedious)
       Pos(1) = other%ConnectList(other%FairIdList(i))%conX ! set initial position of each fairlead i (IFP)
       Pos(2) = other%ConnectList(other%FairIdList(i))%conY                                                 !          |
       Pos(3) = other%ConnectList(other%FairIdList(i))%conZ                                                !          |
                                                                       !          |
       CALL MeshPositionNode(u%PtFairleadDisplacement,i,Pos,ErrStat,ErrMsg)! "assign the coordinates of each node in the global coordinate space"
       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))             !           |

       ! also set velocity of each node to zero (manually)
       u%PtFairleadDisplacement%TranslationVel(1,i) = 0.0_ReKi
       u%PtFairleadDisplacement%TranslationVel(2,i) = 0.0_ReKi
       u%PtFairleadDisplacement%TranslationVel(3,i) = 0.0_ReKi

       ! set each node as a point element
       CALL MeshConstructElement(u%PtFairleadDisplacement , &              !          |
                                 ELEMENT_POINT        , &              !          |
                                 ErrStat              , &              !          |
                                 ErrMsg               , &              !          |
                                 i)                                    !          |
       IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))             !          |
    END DO    ! I                                                      !          |
                                                                       !          |
    CALL MeshCommit ( u%PtFairleadDisplacement, ErrStat, ErrMsg )      !          |
    IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                !          |
                                                                       !          |
    ! now, copy the input PtFairleadDisplacement to output    PtFairleadLoad
    CALL MeshCopy ( SrcMesh  = u%PtFairleadDisplacement , &                !          |
                    DestMesh = y%PtFairleadLoad     , &                !          |
                    CtrlCode = MESH_SIBLING         , &                !          |
                    Force    = .TRUE.               , &                !          |
                    ErrStat  = ErrStat              , &                !          |
                    ErrMess  = ErrMsg                 )                !          |
    IF (ErrStat /= ErrID_None) CALL WrScr(TRIM(ErrMsg))                !          |
                                                                       !          |
    y%PtFairleadLoad%IOS = COMPONENT_OUTPUT                            !          |
    ! End mesh initialization                                          !   -------+
    !==============================================================================



    ! go through all Connects and set position based on input file

    DO I = 1, p%NConnects

      IF ( other%ConnectList(I)%TypeNum == 0 ) THEN !fixed
        other%ConnectList(I)%r(1) = other%ConnectList(I)%conX
        other%ConnectList(I)%r(2) = other%ConnectList(I)%conY
        other%ConnectList(I)%r(3) = other%ConnectList(I)%conZ
        other%ConnectList(I)%rd(1) = 0.0_ReKi
        other%ConnectList(I)%rd(2) = 0.0_ReKi
        other%ConnectList(I)%rd(3) = 0.0_ReKi
      ELSE IF ( other%ConnectList(I)%TypeNum == 1 ) THEN !vessel
        other%ConnectList(I)%r(1) = other%ConnectList(I)%conX  ! this is temporary - doesn't account for platform initial position
        other%ConnectList(I)%r(2) = other%ConnectList(I)%conY
        other%ConnectList(I)%r(3) = other%ConnectList(I)%conZ
        other%ConnectList(I)%rd(1) = 0.0_ReKi
        other%ConnectList(I)%rd(2) = 0.0_ReKi
        other%ConnectList(I)%rd(3) = 0.0_ReKi

!        r[0] = TransMat[0]*conX + TransMat[1]*conY + TransMat[2]*conZ + pX[0];  // x
!        r[1] = TransMat[3]*conX + TransMat[4]*conY + TransMat[5]*conZ + pX[1];  // y
!        r[2] = TransMat[6]*conX + TransMat[7]*conY + TransMat[8]*conZ + pX[2];  // z
!
!        for (int I=0; I<3; I++) rd(I) = 0.0;
      ELSE IF ( other%ConnectList(I)%TypeNum == 2 ) THEN !connect
        other%ConnectList(I)%r(1) = other%ConnectList(I)%conX        ! put node at initial guess!!!
        other%ConnectList(I)%r(2) = other%ConnectList(I)%conY
        other%ConnectList(I)%r(3) = other%ConnectList(I)%conZ
        other%ConnectList(I)%rd(1) = 0.0_ReKi
        other%ConnectList(I)%rd(2) = 0.0_ReKi
        other%ConnectList(I)%rd(3) = 0.0_ReKi
      END IF

      ! assign to state vector
      !DO J = 1, 3
      x%states(6*I-5:6*I-3) = other%ConnectList(I)%r !X[3 + I] = r(I);
      x%states(6*I-2:6*I) = other%ConnectList(I)%rd  !X[    I] = rd(I);

    END DO  !I = 1, p%NConnects


    ! ------------------------------- open output file(s) and write header lines ---------------------------
    CALL MDIO_OpenOutput( InitInp%FileName,  p, other, InitOut, ErrStat, ErrMsg )



    ! go through lines and initialize internal node positions using Catenary()
    DO I = 1, p%NLines

      !get ICs for line using quasi-static approach

      N = other%LineList(I)%N ! for convenience

      ! set end node positions and velocities from connect objects
      other%LineList(I)%r(:,N) = other%ConnectList(other%LineList(I)%FairConnect)%r
      other%LineList(I)%r(:,0) = other%ConnectList(other%LineList(I)%AnchConnect)%r
      other%LineList(I)%rd(:,N) = (/ 0.0, 0.0, 0.0 /)  ! set anchor end velocities to zero
      other%LineList(I)%rd(:,0) = (/ 0.0, 0.0, 0.0 /)  ! set fairlead end velocities to zero

      print *, 'Line', other%LineList(I)%IdNum, ' PropsIdNum is ', other%LineList(I)%PropsIdNum, ' ... ', other%LineTypeList(other%LineList(I)%PropsIdNum)%IdNum

      ! set initial line internal node positions using quasi-static model or straight-line interpolation from anchor to fairlead
      CALL InitializeLine( other%LineList(I), other%LineTypeList(other%LineList(I)%PropsIdNum), p%rhoW ,  ErrStat, ErrMsg)


      ! assign the resulting internal node positions to the integrator initial state vector! (velocities leave at 0)
      DO J = 1, N-1
        DO K = 1, 3
          x%states(other%LineStateIndList(I) + 3*N-3 + 3*J-3 + K-1 ) = other%LineList(I)%r(K,J) ! assign position
          x%states(other%LineStateIndList(I) +         3*J-3 + K-1 ) = 0.0_ReKi ! assign velocities (of zero)
        END DO
      END DO

    END DO    !I = 1, p%NLines


    ! try writing output for troubleshooting purposes (TEMPORARY)
    !CALL MDIO_WriteOutputs(-1.0_DbKi, p, other, y, ErrStat, ErrMsg)


    ! ----------------- do dynamic relaxation to get ICs ------------------------

    CALL WrScr("   Finalizing ICs using dynamic relaxation.")

    ! boost drag coefficient of each line type
    DO I = 1, p%NTypes
      other%LineTypeList(I)%Cdn = other%LineTypeList(I)%Cdn * InitInp%CdScaleIC
      other%LineTypeList(I)%Cdt = other%LineTypeList(I)%Cdt * InitInp%CdScaleIC
    END DO

    ! allocate array holding three latest fairlead tensions
    ALLOCATE ( FairTensIC(p%NFairs,3), STAT = ErrStat )
    IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Error allocating FairTensIC array.'
      ! handle how?
    END IF

    ! initialize fairlead tension memory at zero
    DO J = 1,p%NFairs
      DO I = 1, 3
        FairTensIC(J,I) = 0.0_ReKi
      END DO
    END DO

    t = 0.0_ReKi     ! start time at zero
    p%Ndt = ceiling(InitInp%DTIC/InitInp%DTmooring)  ! get number of mooring time steps per IC convergence analysis time step
    p%dt = InitInp%DTIC/p%Ndt                       ! adjusted mooring time step size for IC generation

    print *, InitInp%DTIC
    print *, InitInp%DTmooring
    print *, p%Ndt

    print *, 'about to start IC time stepping'

    print *, InitInp%TMaxIC
    print *, p%dt
    print *, ceiling(InitInp%TMaxIC/InitInp%DTIC)

    DO I = 1, ceiling(InitInp%TMaxIC/InitInp%DTIC)   ! loop through IC gen time steps, up to maximum

    !  print *, 'Time step ', I, ', t = ', t

      !double t = iic*ICdt;         // IC gen time (s).

      !// loop through line integration time steps
      !for (double ts=t; ts<=t+ICdt-dts; ts+=dts)
      !  rk2 (states, ts, dts )

      CALL TimeStep ( t, InitInp%DTIC, u, p, x, xd, z, other, ErrStat, ErrMsg )

    ! try writing output for troubleshooting purposes (TEMPORARY)
!    CALL MDIO_WriteOutputs(REAL(t,DbKi) , p, other, ErrStat, ErrMsg)


      ! store previous fairlead tensions for comparison
      DO J = 1, p%NFairs
        FairTensIC(J,3) = FairTensIC(J,2)
        FairTensIC(J,2) = FairTensIC(J,1)
        FairTensIC(J,1) = TwoNorm(other%ConnectList(other%FairIdList(J))%Ftot(:))

    !    print *, "   t = ", t , " s, tension at Node " , other%FairIdList(J), " is ", FairTensIC(J,1)

      END DO

      ! check for convergence (compare current tension at each fairlead with previous two values)
      IF (I > 2) THEN
        DO J = 1, p%NFairs
          IF (( abs( FairTensIC(J,1)/FairTensIC(J,2) - 1.0 ) > InitInp%threshIC ) .OR. ( abs( FairTensIC(J,1)/FairTensIC(J,3) - 1.0 ) > InitInp%threshIC ) ) THEN
            EXIT
          END IF
        END DO

        IF (J == p%NFairs) THEN   ! if we made it with all cases satisfying the threshold

          print *,"   Fairlead tensions converged to ", 100.0*InitInp%threshIC, " after ", t , " seconds."
          EXIT  ! break out of the time stepping loop
        END IF
      END IF

      ! should have a 'timed out without convergence' message also

    END DO ! I ... looping through time steps


     ! UNboost drag coefficient of each line type
    DO I = 1, p%NTypes
      other%LineTypeList(I)%Cdn = other%LineTypeList(I)%Cdn / InitInp%CdScaleIC
      other%LineTypeList(I)%Cdt = other%LineTypeList(I)%Cdt / InitInp%CdScaleIC
    END DO



    ! trying writing an actual output line - temporary
    !CALL MDIO_WriteOutputs( 0.0_DbKi, p, other, y, ErrStat, ErrMsg )


    ! set time step stuff for actual simulation
    p%Ndt = ceiling(DTcoupling/InitInp%DTmooring)  ! get number of mooring time steps per IC convergence analysis time step
    p%dt = DTcoupling/p%Ndt                       ! adjusted mooring time step size for IC generation

    p%dtCoupling = DTcoupling  ! not sure if this is "proper"  store coupling time step for use in updatestates

    ! Give the program description (name, version number, date)
    !InitOut%Ver = ProgDesc('MoorDyn',TRIM(InitOut%version),TRIM(InitOut%compilingData)) ! rhis is a duplicate of above


    CALL WrScr( 'Done MD_Init' )


  END SUBROUTINE MD_Init                                                                        !   -------+
  !==========================================================================================================




  !==========   MD_UpdateStates   ======     <-------------------------------------------------------------+
  SUBROUTINE MD_UpdateStates( t, n, u, utimes, p, x, xd, z, other, ErrStat, ErrMsg)
    REAL(DbKi)                      , INTENT(IN   ) :: t
    INTEGER(IntKi)                  , INTENT(IN   ) :: n
    REAL(DbKi)                      , INTENT(IN   ) :: utimes(:)
    TYPE(MD_InputType)              , INTENT(INOUT) :: u(:)       ! INTENT(INOUT) ! had to change this to INOUT
    TYPE(MD_ParameterType)          , INTENT(IN   ) :: p          ! INTENT(IN   )
    TYPE(MD_ContinuousStateType)    , INTENT(INOUT) :: x          ! INTENT(INOUT)
    TYPE(MD_DiscreteStateType)      , INTENT(INOUT) :: xd         ! INTENT(INOUT)
    TYPE(MD_ConstraintStateType)    , INTENT(INOUT) :: z          ! INTENT(INOUT)
    TYPE(MD_OtherStateType)         , INTENT(INOUT) :: other      ! INTENT(INOUT)
    INTEGER(IntKi)                  , INTENT(  OUT) :: ErrStat    ! Error status of the operation
    CHARACTER(*)                    , INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

    INTEGER(IntKi)  :: I ! counter
    INTEGER(IntKi)  :: J ! counter
    INTEGER(IntKi)  :: K ! counter

    TYPE(MD_InputType)         :: u_interp  !

    REAL(ReKi) :: t2  ! why did I do this?


    t2 = real(t, ReKi)

    ! create space for arrays/meshes in u_interp
    CALL MD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg)
    CALL MD_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat, ErrMsg)
    !mth: I don't know what this is but it seems cool!


    !-------------------------------------------------------------------------------
    ! go through fairleads and apply motions from driver
    DO I = 1, p%NFairs
      DO J = 1,3

        other%ConnectList(other%FairIdList(I))%r(J) =  u_interp%PtFairleadDisplacement%Position(J,I) ! assign fairlead position
       !other%ConnectList(other%FairIdList(I))%r(J) =  u_interp%PtFairleadDisplacement%TranslationDisp(J,I) ! is this the same as above?

        other%ConnectList(other%FairIdList(I))%rd(J) =  u_interp%PtFairleadDisplacement%TranslationVel(J,I) ! assign fairlead velocity

      END DO
    END DO  ! I



    !-------------------------------------------------------------------------------
    ! loop through mooring model time steps
    ! DO I = 1, p%Ndt  <<< the looping through mooring time steps is done within TimeStep()!
      CALL TimeStep ( t2, p%dtCoupling, u_interp, p, x, xd, z, other, ErrStat, ErrMsg )
    ! END DO ! I



    !-------------------------------------------------------------------------------
    ! check for errors (return) or warnings (print)  - this seems redundant
    IF( ErrStat /= ErrID_None ) THEN
       IF( ErrStat > ErrID_Warn ) THEN
          CALL MD_DestroyInput( u_interp, ErrStat, ErrMsg )
          RETURN
       ELSE
          CALL WrScr( ErrMsg )
       END IF
    END IF

    CALL MD_DestroyInput(u_interp, ErrStat, ErrMsg)


  END SUBROUTINE MD_UpdateStates                                                                !   -------+
  !==========================================================================================================


  !==========   MD_CalcOutput   ======     <---------------------------------------------------------------+
  SUBROUTINE MD_CalcOutput( t, u, p, x, xd, z, other, y, ErrStat, ErrMsg )
    REAL(DbKi)                     , INTENT(IN   ) :: t
    TYPE( MD_InputType )           , INTENT(IN   ) :: u       ! INTENT(IN   )
    TYPE( MD_ParameterType )       , INTENT(IN   ) :: p       ! INTENT(IN   )
    TYPE( MD_ContinuousStateType ) , INTENT(IN   ) :: x       ! INTENT(IN   )
    TYPE( MD_DiscreteStateType )   , INTENT(IN   ) :: xd      ! INTENT(IN   )
    TYPE( MD_ConstraintStateType ) , INTENT(IN   ) :: z       ! INTENT(IN   )
    TYPE( MD_OtherStateType )      , INTENT(INOUT) :: other   ! INTENT(INOUT)
    TYPE( MD_OutputType )          , INTENT(INOUT) :: y       ! INTENT(INOUT)
    INTEGER(IntKi)                 , INTENT(  OUT) :: ErrStat
    CHARACTER(*)                   , INTENT(  OUT) :: ErrMsg

    INTEGER(IntKi)             :: I  ! counter
    INTEGER(IntKi)             :: J  ! counter


    ! assign net force on fairlead Connects to the outputs
    DO i = 1, p%NFairs
      DO J=1,3
        y%PtFairleadLoad%Force(J,I) = other%ConnectList(other%FairIdList(I))%Ftot(J)
      END DO
    END DO


    ! calculate outputs (y%WriteOutput) for glue code and write any other outputs to MoorDyn output files
     CALL MDIO_WriteOutputs(REAL(t,DbKi) , p, other, y, ErrStat, ErrMsg)


    ! check for errors (return) or warnings (print)  <<< ?
    IF( ErrStat /= ErrID_None ) THEN
       IF( ErrStat > ErrID_Warn ) THEN
          RETURN
       ELSE
          CALL WrScr( ErrMsg )
       END IF
    END IF

  END SUBROUTINE MD_CalcOutput                                                                  !   -------+
  !==========================================================================================================


  !==========   MD_End   ======     <----------------------------------------------------------------------+
  SUBROUTINE MD_End(u, p, x, xd, z, other, y, ErrStat , ErrMsg)
    TYPE(MD_InputType) ,            INTENT(INOUT) :: u
    TYPE(MD_ParameterType) ,        INTENT(INOUT) :: p
    TYPE(MD_ContinuousStateType) ,  INTENT(INOUT) :: x
    TYPE(MD_DiscreteStateType) ,    INTENT(INOUT) :: xd
    TYPE(MD_ConstraintStateType) ,  INTENT(INOUT) :: z
    TYPE(MD_OtherStateType) ,       INTENT(INOUT) :: other
    TYPE(MD_OutputType) ,           INTENT(INOUT) :: y
    INTEGER(IntKi),                 INTENT(  OUT) :: ErrStat
    CHARACTER(*),                   INTENT(  OUT) :: ErrMsg
    INTEGER(IntKi)                                :: i=0

    ErrStat = ErrID_None
    ErrMsg  = ""



    ! FOR YOU MATT
    !
    ! do deallocations here
    ! ....
    ! Done


    ! close output, for now
    CALL MDIO_CloseOutput ( p, ErrStat, ErrMsg )


    ! check for errors (return) or warnings (print)  - this seems redundant
    IF( ErrStat /= ErrID_None ) THEN
       IF( ErrStat > ErrID_Warn ) THEN
          RETURN
       ELSE
          CALL WrScr( ErrMsg )
       END IF
    END IF





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


! DO I NEED SUBROUTINES FOR COPYING AND DESTROYIG THE DERIVED DATA TYPES, AS STATED IN P13 OF THE PROGRAMMING HANDBOOK?
! I.E. COPY AND DELETE FUNCTIONS FOR CONNECTS LINES AND PROPTYPES?


  !========================================================================================================
  SUBROUTINE TimeStep ( t, dt, u, p, x, xd, z, other, ErrStat, ErrMsg )
    REAL(ReKi)                     , INTENT(INOUT) :: t
    REAL(ReKi)                     , INTENT(IN   ) :: dt      ! - DOES THIS EVEN NEED TO BE PASSED? - OUTER time step size (mooring time step size is stored in p)
    TYPE( MD_InputType )           , INTENT(IN   ) :: u       ! INTENT(IN   )
    TYPE( MD_ParameterType )       , INTENT(IN   ) :: p       ! INTENT(IN   )
    TYPE( MD_ContinuousStateType ) , INTENT(INOUT) :: x
    TYPE( MD_DiscreteStateType )   , INTENT(IN   ) :: xd      ! INTENT(IN   )
    TYPE( MD_ConstraintStateType ) , INTENT(IN   ) :: z       ! INTENT(IN   )
    TYPE( MD_OtherStateType )      , INTENT(INOUT) :: other   ! INTENT(INOUT)
 !   TYPE( MD_OutputType )          , INTENT(INOUT) :: y       ! INTENT(INOUT)
    INTEGER(IntKi)                 , INTENT(  OUT) :: ErrStat
    CHARACTER(*)                   , INTENT(  OUT) :: ErrMsg


    INTEGER(IntKi)  :: I ! counter
    INTEGER(IntKi)  :: J ! counter
    INTEGER(IntKi)  :: K ! counter


    !loop through line integration time steps
    DO I = 1, p%Ndt  !for (double ts=t; ts<=t+ICdt-dts; ts+=dts)

      ! try writing output for troubleshooting purposes (TEMPORARY)
      ! CALL MDIO_WriteOutputs(REAL(t,DbKi) , p, other, ErrStat, ErrMsg)


      CALL RK2 (x%states, other%F, t, p%dt, ErrStat, ErrMsg)  ! pass state vector, temporary vector for calculations, instantaneous time, and time step size to integration
        !rk2 (states, ts, dts )

      ! update Fairlead positions by integrating velocity and last position (do this AFTER the processing of the time step rather than before)
      DO J = 1, p%NFairs
        DO K = 1, 3
          other%ConnectList(other%FairIdList(J))%r(K) = other%ConnectList(other%FairIdList(J))%r(K) + other%ConnectList(other%FairIdList(J))%rd(K)*p%dt
        END DO
      END DO

    END DO

    !is anything else required here at the end???

  CONTAINS


    !=====================================================================
    SUBROUTINE RK2 ( X, F, t, dt, ErrStat, ErrMsg )
    ! Do RK2 integration

      Real(ReKi),     INTENT( INOUT )   :: X(:) ! state vector
      Real(ReKi),     INTENT( INOUT )   :: F(:) ! vector used for integration calculations (passed to avoid allocating here)
      Real(ReKi),     INTENT (INOUT)    :: t    ! time
      Real(ReKi),     INTENT (IN)       :: dt    ! time step size

      INTEGER,        INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs
      CHARACTER(*),   INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None

      !Real(ReKi)                        :: time

      INTEGER(IntKi)                     :: I ! index

      ErrStat = ErrID_None
      ErrMsg  = ""

      CALL RHSmaster(X, F, t)                         !f0 = f ( t, x0 );

      DO I = 1, size(X)
        X(I) = X(I) + 0.5*dt*F(I)                  !x1 = x0 + dt*f0/2.0;
      END DO

      !time = t + dt/2.0

      CALL RHSmaster(X, F, (t + 0.5*dt) )             !f1 = f ( t1, x1 );

      DO I = 1, size(X)
        X(I) = X(I) + dt*F(I)
      END DO

      t = t + dt  ! update time

    END SUBROUTINE RK2
    !=================================================================


    !======================================================================
    SUBROUTINE RHSMaster( X,  Xd, t)
      ! currently, this subroutine does all the calculations for the connections and lines

      Real(ReKi), INTENT( IN )  :: X(:)  ! state vector, provided
      Real(ReKi), INTENT( OUT ) :: Xd(:) ! derivative of state vector, returned
      Real(ReKi), INTENT (IN)   :: t    ! instantaneous time

      INTEGER(IntKi)                     :: L ! index
      INTEGER(IntKi)                     :: Istart ! start index of line/connect in state vector
      INTEGER(IntKi)                     :: Iend ! end index of line/connect in state vector


      ! clear connection force and mass values
      DO L = 1, p%NConnects
        DO J = 1,3
          other%ConnectList(L)%Ftot(J) = 0.0_ReKi
          other%ConnectList(L)%Ftot(J) = 0.0_ReKi
          DO K = 1,3
            other%ConnectList(L)%Mtot(K,J) = 0.0_ReKi
            other%ConnectList(L)%Mtot(K,J) = 0.0_ReKi
          END DO
        END DO
      END DO

      ! do Line force and acceleration calculations, also add end masses/forces to respective Connects
      DO L = 1, p%NLines
        Istart = other%LineStateIndList(L)
        Iend = Istart + 6*(other%LineList(L)%N-1) - 1
        CALL DoLineRHS(X(Istart:Iend), Xd(Istart:Iend), t, other%LineList(L), &
          other%LineTypeList(other%LineList(L)%PropsIdNum), &
          other%ConnectList(other%LineList(L)%FairConnect)%Ftot, other%ConnectList(other%LineList(L)%FairConnect)%Mtot, &
          other%ConnectList(other%LineList(L)%AnchConnect)%Ftot, other%ConnectList(other%LineList(L)%AnchConnect)%Mtot )
      END DO

      ! do Connect acceleration calculations
      DO L = 1, p%NConnects
        Istart = L*6-5
        Iend = L*6
        CALL DoConnectRHS(X(Istart:Iend), Xd(Istart:Iend), t, other%ConnectList(L))
      END DO


    END SUBROUTINE RHSMaster
    !=====================================================================



    !======================================================================
    SUBROUTINE DoLineRHS (X, Xd, t, Line, LineProp, FairFtot, FairMtot, AnchFtot, AnchMtot)

      Real(ReKi), INTENT( IN )  :: X(:)               ! state vector, provided
      Real(ReKi), INTENT( OUT ) :: Xd(:)              ! derivative of state vector, returned
      Real(ReKi), INTENT (IN)   :: t                  ! instantaneous time
      TYPE(MD_Line), INTENT (INOUT)  :: Line          ! label for the current line, for convenience
      TYPE(MD_LineProp), INTENT(IN) :: LineProp       ! the single line property set for the line of interest
      Real(ReKi), INTENT(INOUT)     :: FairFtot(:)    ! total force on Connect top of line is attached to
      Real(ReKi), INTENT(INOUT)     :: FairMtot(:,:)  ! total mass of Connect top of line is attached to
      Real(ReKi), INTENT(INOUT)     :: AnchFtot(:)    ! total force on Connect bottom of line is attached to
      Real(ReKi), INTENT(INOUT)     :: AnchMtot(:,:)  ! total mass of Connect bottom of line is attached to


      INTEGER(IntKi)             :: I         ! index of segments or nodes along line
      INTEGER(IntKi)             :: J         ! index
      INTEGER(IntKi)             :: K         ! index
      INTEGER(IntKi)             :: N         ! number of segments in line
      Real(ReKi)                 :: d         ! line diameter
      Real(ReKi)                 :: rho       ! line material density
      Real(ReKi)                 :: Sum1      ! for summing squares
      Real(ReKi)                 :: m_i       ! node mass
      Real(ReKi)                 :: v_i       ! node submerged volume
      Real(ReKi)                 :: Vi(3)     ! relative water velocity at a given node
      Real(ReKi)                 :: Vp(3)     ! transverse relative water velocity component at a given node
      Real(ReKi)                 :: Vq(3)     ! tangential relative water velocity component at a given node
      Real(ReKi)                 :: SumSqVp   !
      Real(ReKi)                 :: SumSqVq   !
      Real(ReKi)                 :: MagVp   !
      Real(ReKi)                 :: MagVq   !

      N = Line%N                      ! for convenience
      d = LineProp%d                  ! for convenience
      rho = LineProp%w/(Pi/4.0*d*d)


      ! set end node positions and velocities from connect objects' states
      DO J = 1, 3
        Line%r(J,N) = other%ConnectList(Line%FairConnect)%r(J)
        Line%r(J,0) = other%ConnectList(Line%AnchConnect)%r(J)
      END DO

      ! set interior node positions and velocities
      DO I = 1, N-1
        DO J = 1, 3
          Line%r(J,I) = X( 3*N-3 + 3*I-3 + J)      ! r(J,I)  = X[3*N-3 + 3*i-3 + J]; // get positions  .. used to start from other%LineStateIndList(Line%IdNum) in whole state vector
          Line%rd(J,I) = X(         3*I-3 + J)      ! rd(J,I) = X[        3*i-3 + J]; // get velocities
        END DO
      END DO

 !     print *, ' r(1,...) are ', Line%r(1,:)

      ! calculate instantaneous (stretched) segment lengths and rates << should add catch here for if lstr is ever zero
      DO I = 0, N-1
        Sum1 = 0.0_ReKi
        DO J = 1, 3
          Sum1 = Sum1 + (Line%r(J,I+1) - Line%r(J,I))*(Line%r(J,I+1) - Line%r(J,I))
        END DO
        Line%lstr(I) = sqrt(Sum1)                                ! stretched segment length

 !   print *, '  lstr() is ', Line%lstr(I)

        Sum1 = 0.0_ReKi
        DO J = 1, 3
          Sum1 = Sum1 + (Line%r(J,I+1) - Line%r(J,I))*(Line%rd(J,I+1) - Line%rd(J,I))
        END DO
        Line%lstrd(I) = Sum1/Line%lstr(I)                        ! strain rate of segment

    !    Line%V(I) = Pi/4.0 * d*d*Line%l(I)     !volume attributed to segment
      END DO

      !calculate unit tangent vectors (q) for each node (including ends)
      CALL UnitVector(Line%q(:,0), Line%r(:,1), Line%r(:,0)) ! compute unit vector q
      DO I = 1, N-1
        CALL UnitVector(Line%q(:,I), Line%r(:,I+1), Line%r(:,I-1)) ! compute unit vector q ... using adjacent two nodes!
      END DO
      CALL UnitVector(Line%q(:,N), Line%r(:,N), Line%r(:,N-1))    ! compute unit vector q


      ! wave kinematics not implemented yet


      !calculate mass (including added mass) matrix for each node
      DO I = 0, N
        IF (I==0) THEN
          m_i = Pi/8.0 *d*d*Line%l(0)*rho
          v_i = 0.5 *Line%V(I)
        ELSE IF (I==N) THEN
          m_i = pi/8.0 *d*d*Line%l(N-1)*rho;
          v_i = 0.5*Line%V(I-1)
        ELSE
          m_i = pi/8.0 * d*d*rho*(Line%l(I) + Line%l(I+1))
          v_i = 0.5 *(Line%V(I-1) + Line%V(I))
        END IF

  !      print *, '  m() is ', m_i

        DO J=1,3
          DO K=1,3
            IF (J==K) THEN
              Line%M(K,J,I) = m_i + p%rhoW*v_i*( LineProp%Can*(1 - Line%q(J,I)*Line%q(K,I)) + LineProp%Cat*Line%q(J,I)*Line%q(K,I) )
            ELSE
              Line%M(K,J,I) = p%rhoW*v_i*( LineProp%Can*(-Line%q(J,I)*Line%q(K,I)) + LineProp%Cat*Line%q(J,I)*Line%q(K,I) )
            END IF
          END DO
        END DO

        CALL Inverse3by3(Line%S(:,:,I), Line%M(:,:,I))   ! invert mass matrix
      END DO


      ! ------------------  CALCULATE FORCES ON EACH NODE ----------------------------

      ! loop through the segments
      DO I = 0, N-1

        ! line tension
        IF (Line%lstr(I)/Line%l(I) > 1.0) THEN
          DO J = 1, 3
            Line%T(J,I) = LineProp%EA *( 1.0/Line%l(I) - 1.0/Line%lstr(I) ) * (Line%r(J,I+1)-Line%r(J,I))
          END DO
        ELSE
          DO J = 1, 3
            Line%T(J,I) = 0.0_ReKi  ! cable can't "push"
          END DO
        END if

  !    print *, '  T is ', Line%T(:,I)

        ! line internal damping force
        DO J = 1, 3
          Line%Td(J,I) = LineProp%BA* ( Line%lstrd(I) / Line%l(I) ) * (Line%r(J,I+1)-Line%r(J,I)) / Line%lstr(I)  ! note new form of damping coefficient, BA rather than Cint
        END DO
      END DO



      ! loop through the nodes
      DO I = 0, N

        !submerged weight (including buoyancy)
        IF (I==0) THEN
          Line%W(3,I) = Pi/8.0*d*d* Line%l(I)*(rho - p%rhoW) *(-p%g)   ! assuming g is positive
        ELSE IF (i==N)  THEN
          Line%W(3,I) = pi/8.0*d*d* Line%l(I)*(rho - p%rhoW) *(-p%g)
        ELSE
          Line%W(3,I) = pi/8.0*d*d* (Line%l(I)*(rho - p%rhoW) + Line%l(I+1)*(rho - p%rhoW) )*(-p%g)  ! left in this form for future free surface handling
        END IF


 !     print *, '  W is ', Line%W(:,I)

        !relative flow velocities
        DO J = 1, 3
          Vi(J) = 0.0 - Line%rd(J,I)                          ! relative flow velocity over node -- this is where wave velicites would be added
        END DO

        ! decomponse relative flow into components
        SumSqVp = 0.0_ReKi                                         ! start sums of squares at zero
        SumSqVq = 0.0_ReKi
        DO J = 1, 3
          Vq(J) = DOT_PRODUCT( Vi , Line%q(:,I) ) * Line%q(J,I);   ! tangential relative flow component
          Vp(J) = Vi(J) - Vq(J)                                    ! transverse relative flow component
          SumSqVq = SumSqVq + Vq(J)*Vq(J)
          SumSqVp = SumSqVp + Vp(J)*Vp(J)
        END DO
        MagVp = sqrt(SumSqVp)                                   ! get magnitudes of flow components
        MagVq = sqrt(SumSqVq)

        ! transverse and tangenential drag
        IF (I==0) THEN
          DO J = 1, 3
            Line%Dp(J,I) = 0.25*p%rhoW*LineProp%Cdn*    d*Line%l(I) * MagVp * Vp(J)
            Line%Dq(J,I) = 0.25*p%rhoW*LineProp%Cdt* Pi*d*Line%l(I) * MagVq * Vq(J)
          END DO
        ELSE IF (I==N)  THEN
          DO J = 1, 3
            Line%Dp(J,I) = 0.25*p%rhoW*LineProp%Cdn*    d*Line%l(I-1) * MagVp * Vp(J);
            Line%Dq(J,I) = 0.25*p%rhoW*LineProp%Cdt* Pi*d*Line%l(I-1) * MagVq * Vq(J)
          END DO
        ELSE
          DO J = 1, 3
            Line%Dp(J,I) = 0.25*p%rhoW*LineProp%Cdn*    d*(Line%l(I) + Line%l(I-1)) * MagVp * vp(J);
            Line%Dq(J,I) = 0.25*p%rhoW*LineProp%Cdt* Pi*d*(Line%l(I) + Line%l(I-1)) * MagVq * vq(J);
          END DO
        END IF

        ! F-K force from fluid acceleration not implemented yet

        ! bottom contact (stiffness and damping)
        IF (Line%r(3,I) < -p%WtrDpth) THEN
          Line%B(3,I) = ( (-p%WtrDpth - Line%r(3,I))*p%kBot - Line%rd(3,I)*p%cBot) * 0.5*d*(Line%l(I) + Line%l(I-1) ) ! vertical only for now
        ELSE
          Line%B(3,I) = 0.0_ReKi
        END IF

        ! total forces
        IF (I==0)  THEN
          DO J = 1, 3
            Line%F(J,I) = Line%T(J,I)                 + Line%Td(J,I)                  + Line%W(J,I) + Line%Dp(J,I) + Line%Dq(J,I) + Line%B(J,I)
          END DO
        ELSE IF (I==N)  THEN
          DO J = 1, 3
            Line%F(J,I) =              -Line%T(J,I-1)                - Line%Td(J,I-1) + Line%W(J,I) + Line%Dp(J,I) + Line%Dq(J,I) + Line%B(J,I)
          END DO
        ELSE
          DO J = 1, 3
            Line%F(J,I) = Line%T(J,I) - Line%T(J,I-1) + Line%Td(J,I) - Line%Td(J,I-1) + Line%W(J,I) + Line%Dp(J,I) + Line%Dq(J,I) + Line%B(J,I)
          END DO
        END IF

      END DO  ! I  - done looping through nodes


    ! some checks

!    I = 4
!    print *, 'N6 (1) l is ', Line%l(I)
!    print *, 'N6 (1) M11 is ', Line%M(1,1,I)
!    print *, 'N6 (1) T is ', Line%T(1,I)
!    print *, 'N6 (1) Td is ', Line%Td(1,I)
!    print *, 'N6 (1) W is ', Line%W(1,I)
!    print *, 'N6 (1) Dp is ', Line%Dp(1,I)
!    print *, 'N6 (1) Dq is ', Line%Dq(1,I)
!    print *, 'N6 (1) B is ', Line%B(1,I)


      ! loop through internal nodes and update their states
      DO I=1, N-1
        DO J=1,3
          ! calculate RHS constant (premultiplying force vector by inverse of mass matrix  ... i.e. rhs = S*Forces)
          Sum1 = 0.0_ReKi                              ! reset temporary accumulator
          DO K = 1, 3
            Sum1 = Sum1 + Line%S(K,J,I) * Line%F(K,I)   ! matrix-vector multiplication [S i]{Forces i}  << double check indices
          END DO ! K

          ! update states
          Xd(3*N-3 + 3*I-3 + J) = X(3*I-3 + J);         ! dxdt = V  (velocities)
          Xd(        3*I-3 + J) = Sum1                ! dVdt = RHS * A  (accelerations)
        END DO ! J
      END DO  ! I

      ! add force and mass of end nodes to the Connects they correspond to!
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
    SUBROUTINE DoConnectRHS (X, Xd, t, Connect)  ! by passing with these only am I somehow increasing the compartmentalization/safety of the code???  can I ??
      !-------------------- RHS calculations for Connects -------------------

      Real(ReKi),       INTENT( IN )    :: X(:)           ! state vector, provided
      Real(ReKi),       INTENT( OUT )   :: Xd(:)          ! derivative of state vector, returned
      Real(ReKi),       INTENT (IN)     :: t              ! instantaneous time
      Type(MD_Connect), INTENT (INOUT)  :: Connect        ! Connect number


      INTEGER(IntKi)             :: I         ! index of segments or nodes along line
      INTEGER(IntKi)             :: J         ! index
      INTEGER(IntKi)             :: K         ! index
      Real(ReKi)                 :: Sum1      ! for adding things


      ! ----------------------------------------------------------------------------------------------------------
      ! the force and mass contributions from the attached Lines should already have been added to
      ! Fto and Mtot by the Line RHS function

      ! add Connect's own forces including buoyancy and weight  (should this only be done for Connect type Connects?)
      Connect%Ftot(1) = Connect%Ftot(1) + Connect%conFX
      Connect%Ftot(2) = Connect%Ftot(2) + Connect%conFY
      Connect%Ftot(3) = Connect%Ftot(3) + Connect%conFZ + Connect%conV*p%rhoW*p%g - Connect%conM*p%g

      ! add Connect's own mass
      DO J = 1,3
        Connect%Mtot(J,J) = Connect%Mtot(J,J) + Connect%conM
      END DO


      ! ------ behavior dependant on connect type -------

      IF (Connect%TypeNum==0)  THEN ! fixed type
        Connect%r(1) = Connect%conX
        Connect%r(2) = Connect%conY
        Connect%r(3) = Connect%conZ
        DO J = 1,3
          Connect%rd(J) = 0.0_ReKi
        END DO
      ELSE IF (Connect%TypeNum==1)  THEN ! vessel type (moves with platform)

        ! fairlead positions are updated from the previous value by integrating the velocity (assumed constant during the driver time step)
        ! this is done in subroutine TimeStep

      ELSE IF (Connect%TypeNum==2)  THEN ! "connect" type

        IF (EqualRealNos(t, 0.0)) THEN  ! this is old: with current IC gen approach, we skip the first call to the line objects, because they're set AFTER the call to the connects

          DO J = 1,3
            Xd(3+I) = X(I)        ! velocities - these are unused in integration
            Xd(I) = 0.0_ReKi           ! accelerations - these are unused in integration
          END DO
        ELSE
          ! from state values, get r and rdot values
          DO J = 1,3
            Connect%r(J)  = X(3 + J)   ! get positions
            Connect%rd(J) = X(J)       ! get velocities
          END DO
        END IF

        ! invert node mass matrix
        CALL Inverse3by3(Connect%S, Connect%Mtot)

        DO J = 1,3
          ! RHS constant - (premultiplying force vector by inverse of mass matrix  ... i.e. rhs = S*Forces
          Sum1 = 0.0_ReKi   ! reset accumulator
          DO K = 1, 3
            Sum1 = Sum1 + Connect%S(K,J) * Connect%Ftot(K)   !  matrix multiplication [S i]{Forces i}
          END DO

          ! update states
          Xd(3 + J) = X(J)          ! dxdt = V    (velocities)
          Xd(I) = Sum1              ! dVdt = RHS * A  (accelerations)
        END DO
      END IF

    END SUBROUTINE DoConnectRHS
    !=====================================================================

  END SUBROUTINE TimeStep
  !======================================================================



  !=======================================================================
  SUBROUTINE SetupLine (Line, LineProp, rhoW, ErrStat, ErrMsg)
    ! calculate initial profile of the line using quasi-static model

    TYPE(MD_Line), INTENT(INOUT)       :: Line            ! the single line object of interest
    TYPE(MD_LineProp), INTENT(INOUT)   :: LineProp    ! the single line property set for the line of interest
    REAL(ReKi),    INTENT(IN)          :: rhoW
    INTEGER,       INTENT(   INOUT )   :: ErrStat              ! returns a non-zero value when an error occurs
    CHARACTER(*),  INTENT(   INOUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None

     INTEGER(4)                   :: I                                               ! Generic index
     INTEGER(4)                   :: J                                               ! Generic index
     INTEGER(4)                   :: K                                               ! Generic index

    INTEGER(IntKi)         :: N

      N = Line%N  ! number of segments in this line (for code readability)

        ! allocate node positions and velocities (NOTE: these arrays start at ZERO)
        ALLOCATE ( Line%r(3, 0:N), Line%rd(3, 0:N), STAT = ErrStat )
        IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating r and rd arrays.'
         !CALL CleanUp()
         RETURN
        END IF

        ! allocate node tangent vectors
        ALLOCATE ( Line%q(3, N+1), STAT = ErrStat )
        IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating q array.'
         !CALL CleanUp()
         RETURN
        END IF

        ! allocate segment scalar quantities
        ALLOCATE ( Line%l(N), Line%lstr(N), Line%lstrd(N), Line%V(N), STAT = ErrStat )
        IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating segment scalar quantity arrays.'
         !CALL CleanUp()
         RETURN
        END IF

        ! assign values for l and V
        DO J=0,N
          Line%l(J) = Line%UnstrLen/REAL(N, DbKi)
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
        ALLOCATE ( Line%W(3, N+1), Line%Dp(3, N+1), Line%Dq(3, N+1), Line%Ap(3, N+1), &
         Line%Aq(3, N+1), Line%B(3, N+1), Line%F(3, N+1), STAT = ErrStat )
        IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating node force arrays.'
         !CALL CleanUp()
         RETURN
        END IF

        ! set gravity and bottom contact forces to zero initially (because the horizontal components should remain at zero)
        DO J = 1,N+1
          DO K = 1,3
            Line%W(K,J) = 0.0_ReKi
            Line%B(K,J) = 0.0_ReKi
          END DO
        END DO

        ! allocate mass and inverse mass matrices for each node (including ends)
        ALLOCATE ( Line%S(3, 3, N+1), Line%M(3, 3, N+1), STAT = ErrStat )
        IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error allocating T and Td arrays.'
         !CALL CleanUp()
         RETURN
        END IF

  END SUBROUTINE SetupLine




  !===============================================================================================
  SUBROUTINE InitializeLine (Line, LineProp, rhoW, ErrStat, ErrMsg)
    ! calculate initial profile of the line using quasi-static model

    TYPE(MD_Line), INTENT(INOUT)       :: Line            ! the single line object of interest
    TYPE(MD_LineProp), INTENT(INOUT)   :: LineProp    ! the single line property set for the line of interest
    REAL(ReKi),    INTENT(IN)          :: rhoW
    INTEGER,       INTENT(   INOUT )   :: ErrStat              ! returns a non-zero value when an error occurs
    CHARACTER(*),  INTENT(   INOUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None

     REAL(ReKi)                   :: COSPhi                                          ! Cosine of the angle between the xi-axis of the inertia frame and the X-axis of the local coordinate system of the current mooring line (-)
     REAL(ReKi)                   :: SINPhi                                          ! Sine   of the angle between the xi-axis of the inertia frame and the X-axis of the local coordinate system of the current mooring line (-)
     REAL(ReKi)                   :: XF                                              ! Horizontal distance between anchor and fairlead of the current mooring line (meters)
     REAL(ReKi)                   :: ZF                                              ! Vertical   distance between anchor and fairlead of the current mooring line (meters)

     INTEGER(4)                   :: I                                               ! Generic index
  !   INTEGER(4)                   :: IndRdtn                                         ! Generic index for the radiation problem
     INTEGER(4)                   :: J                                               ! Generic index
  !   INTEGER(4)                   :: JNode                                           ! The index of the current platform node / element (-) [1 to PtfmNodes]
  !   INTEGER(4)                   :: K                                               ! Generic index
  !   INTEGER(4), SAVE             :: LastIndRdtn                                     ! Index into the radiation     arrays saved from the last call as a starting point for this call.
  !   INTEGER(4), SAVE             :: LastIndRdtn2                                    ! Index into the radiation     arrays saved from the last call as a starting point for this call.
  !   INTEGER(4), SAVE             :: LastIndWave = 1                                 ! Index into the incident wave arrays saved from the last call as a starting point for this call.


    REAL(ReKi)                   :: WetWeight
    REAL(ReKi)                   :: SeabedCD = 0.0_ReKi
    REAL(ReKi)                   :: TenTol = 0.0001_ReKi
    REAL(ReKi), ALLOCATABLE      ::LSNodes(:)
    REAL(ReKi), ALLOCATABLE      ::LNodesX(:)
    REAL(ReKi), ALLOCATABLE      ::LNodesZ(:)
    INTEGER(IntKi)         :: N


  !  print *, 'initializing Line', Line%IdNum, ' of LineType ', LineProp%IdNum

  !  print *, 'LineProp EA is ', LineProp%EA


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
              COSPhi  = 0.0_ReKi
              SINPhi  = 0.0_ReKi
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
     !   CALL CleanUp()     << should handle deallocating arrays in all cases
        RETURN
      END IF

      ALLOCATE ( LNodesX(N+1), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
        ErrMsg  = ' Error allocating LNodesX array.'
     !   CALL CleanUp()     << should handle deallocating arrays in all cases
        RETURN
      END IF

      ALLOCATE ( LNodesZ(N+1), STAT = ErrStat )
      IF ( ErrStat /= ErrID_None ) THEN
        ErrMsg  = ' Error allocating LNodesZ array.'
     !   CALL CleanUp()     << should handle deallocating arrays in all cases
        RETURN
      END IF

      ! Assign node arc length locations
      LSNodes(1) = 0.0_ReKi
      DO I=2,N
        LSNodes(I) = LSNodes(I-1) + Line%l(I-2)  ! note: l index is because line and segment indices start at 0
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

   !       print *, 'about to call Catenary. Args are ', XF, ZF, Line%UnstrLen, LineProp%EA  , &
    !                       WetWeight    , SeabedCD,    TenTol,     (N+1)     ,  LSNodes

           CALL Catenary ( XF           , ZF          , Line%UnstrLen, LineProp%EA  , &
                           WetWeight    , SeabedCD,    TenTol,     (N+1)     , &
                           LSNodes, LNodesX, LNodesZ , ErrStat, ErrMsg)


    IF (ErrStat == ErrID_None) THEN ! if it worked, use it
        ! Transform the positions of each node on the current line from the local
        !   coordinate system of the current line to the inertial frame coordinate
        !   system:

           DO J = 0,Line%N ! Loop through all nodes per line where the line position and tension can be output
              Line%r(1,J) = Line%r(1,0) + LNodesX(J+1)*COSPhi
              Line%r(2,J) = Line%r(2,0) + LNodesX(J+1)*SINPhi
              Line%r(3,J) = Line%r(3,0) + LNodesZ(J+1)
           ENDDO              ! J - All nodes per line where the line position and tension can be output


    ELSE ! if there is a problem with the catenary approach, just stretch the nodes linearly between fairlead and anchor

  !      }
  !      else
  !      {  ! otherwise just stretch the nodes between the endpoints linearly and hope for the best
          print *, 'Catenary IC generation failed so using straight-line approach instead.'

           DO J = 0,Line%N ! Loop through all nodes per line where the line position and tension can be output
              Line%r(1,J) = Line%r(1,0) + (Line%r(1,N) - Line%r(1,0))*REAL(J, ReKi)/REAL(N, ReKi)
              Line%r(2,J) = Line%r(2,0) + (Line%r(2,N) - Line%r(2,0))*REAL(J, ReKi)/REAL(N, ReKi)
              Line%r(3,J) = Line%r(3,0) + (Line%r(3,N) - Line%r(3,0))*REAL(J, ReKi)/REAL(N, ReKi)
              print *, Line%r(1,J)
           ENDDO              ! J - All nodes per line where the line position and tension can be output
    ENDIF



    CALL CleanUp()  ! deallocate temporary arrays

  ! ------------------------------------------------------------------------------

    CONTAINS

      SUBROUTINE CleanUp()

      ! add things to deallocate temporary arrays here

      END SUBROUTINE CleanUp




    !=======================================================================
    !JASON: SHOULD THIS ROUTINE (Catenary) BE PLACED IN NWTC_Subs OR IN ITS OWN DLL?
    SUBROUTINE Catenary ( XF_In, ZF_In, L_In  , EA_In, &
                          W_In , CB_In, Tol_In, N    , &
                          s_In , X_In , Z_In , ErrStat, ErrMsg    )


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

    REAL(ReKi), INTENT(IN   )    :: CB_In                                           ! Coefficient of seabed static friction drag (a negative value indicates no seabed) (-)
    REAL(ReKi), INTENT(IN   )    :: EA_In                                           ! Extensional stiffness of line (N)
  !  REAL(ReKi), INTENT(  OUT)    :: HA_In                                           ! Effective horizontal tension in line at the anchor   (N)
  !  REAL(ReKi), INTENT(INOUT)    :: HF_In                                           ! Effective horizontal tension in line at the fairlead (N)
    REAL(ReKi), INTENT(IN   )    :: L_In                                            ! Unstretched length of line (meters)
    REAL(ReKi), INTENT(IN   )    :: s_In     (N)                                    ! Unstretched arc distance along line from anchor to each node where the line position and tension can be output (meters)
  !  REAL(ReKi), INTENT(  OUT)    :: Te_In    (N)                                    ! Effective line tensions at each node (N)
    REAL(ReKi), INTENT(IN   )    :: Tol_In                                          ! Convergence tolerance within Newton-Raphson iteration specified as a fraction of tension (-)
  !  REAL(ReKi), INTENT(  OUT)    :: VA_In                                           ! Effective vertical   tension in line at the anchor   (N)
  !  REAL(ReKi), INTENT(INOUT)    :: VF_In                                           ! Effective vertical   tension in line at the fairlead (N)
    REAL(ReKi), INTENT(IN   )    :: W_In                                            ! Weight of line in fluid per unit length (N/m)
    REAL(ReKi), INTENT(  OUT)    :: X_In     (N)                                    ! Horizontal locations of each line node relative to the anchor (meters)
    REAL(ReKi), INTENT(IN   )    :: XF_In                                           ! Horizontal distance between anchor and fairlead (meters)
    REAL(ReKi), INTENT(  OUT)    :: Z_In     (N)                                    ! Vertical   locations of each line node relative to the anchor (meters)
    REAL(ReKi), INTENT(IN   )    :: ZF_In                                           ! Vertical   distance between anchor and fairlead (meters)
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



       ! Abort when there is no solution or when the only possible solution is
       !   illogical:

    IF (    Tol <= 0.0_DbKi )  THEN   ! .TRUE. when the convergence tolerance is specified incorrectly
      ErrStat = ErrID_Warn
      ErrMsg = ' Convergence tolerance must be greater than zero in routine Catenary().'

    ELSEIF ( XF <  0.0_DbKi )  THEN   ! .TRUE. only when the local coordinate system is not computed correctly
      ErrStat = ErrID_Warn
      ErrMsg =  ' The horizontal distance between an anchor and its'// &
                    ' fairlead must not be less than zero in routine Catenary().'

    ELSEIF ( ZF <  0.0_DbKi )  THEN   ! .TRUE. if the fairlead has passed below its anchor
      ErrStat = ErrID_Warn
      ErrMsg =  ' A fairlead has passed below its anchor.'

    ELSEIF ( L  <= 0.0_DbKi )  THEN   ! .TRUE. when the unstretched line length is specified incorrectly
      ErrStat = ErrID_Warn
      ErrMsg =  ' Unstretched length of line must be greater than zero in routine Catenary().'

    ELSEIF ( EA <= 0.0_DbKi )  THEN   ! .TRUE. when the unstretched line length is specified incorrectly
      ErrStat = ErrID_Warn
      ErrMsg =  ' Extensional stiffness of line must be greater than zero in routine Catenary().'

    ELSEIF ( W  == 0.0_DbKi )  THEN   ! .TRUE. when the weight of the line in fluid is zero so that catenary solution is ill-conditioned
      ErrStat = ErrID_Warn
      ErrMsg = ' The weight of the line in fluid must not be zero. '// &
                    ' Routine Catenary() cannot solve quasi-static mooring line solution.'


    ELSEIF ( W  >  0.0_DbKi )  THEN   ! .TRUE. when the line will sink in fluid

       LMax      = XF - EA/W + SQRT( (EA/W)*(EA/W) + 2.0_DbKi*ZF*EA/W )  ! Compute the maximum stretched length of the line with seabed interaction beyond which the line would have to double-back on itself; here the line forms an "L" between the anchor and fairlead (i.e. it is horizontal along the seabed from the anchor, then vertical to the fairlead)

       IF ( ( L  >=  LMax   ) .AND. ( CB >= 0.0_DbKi ) )  &  ! .TRUE. if the line is as long or longer than its maximum possible value with seabed interaction
      ErrStat = ErrID_Warn
      ErrMsg =  ' Unstretched mooring line length too large. '// &
                       ' Routine Catenary() cannot solve quasi-static mooring line solution.'


    ENDIF

    IF (ErrStat > ErrID_None) then
      print *, ErrMsg
      RETURN
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

    !HA_In    = REAL( HA   , ReKi )  !mth: for this I only care about returning node positions
    !HF_In    = REAL( HF   , ReKi )
    !Te_In(:) = REAL( Te(:), ReKi )
    !VA_In    = REAL( VA   , ReKi )
    !VF_In    = REAL( VF   , ReKi )
    X_In (:) = REAL( X (:), ReKi )
    Z_In (:) = REAL( Z (:), ReKi )

    print *, 'done Catenary.  N is ', N, ', first three x positions are ', X_In(1), X_In(2), X_In(3)


    END SUBROUTINE Catenary


  END SUBROUTINE InitializeLine



  ! ============ below are some math convenience functions ===============
  ! should add error checking if I keep these, but hopefully there are existing NWTCLib functions to replace them

  ! return unit vector (u) in direction from r1 to r2
  !-----------------------------------------------------------------------
  SUBROUTINE UnitVector( u, r1, r2 )
    REAL(ReKi), INTENT(OUT)   :: u(:)
    REAL(ReKi), INTENT(IN)    :: r1(:)
    REAL(ReKi), INTENT(IN)    :: r2(:)

    REAL(ReKi)               :: Length

    u = r2 - r1    
    Length = TwoNorm(u) 

    if ( .NOT. EqualRealNos(length, 0.0_ReKi ) ) THEN
      u = u / Length
    END IF
    

   END SUBROUTINE UnitVector



  !computes the inverse of a matrix m
  !-----------------------------------------------------------------------
  SUBROUTINE Inverse3by3( Minv, M )
    Real(ReKi), INTENT(OUT)   :: Minv(:,:)  ! returned inverse matrix
    Real(ReKi), INTENT(IN)    :: M(:,:)  ! inputted matrix

    Real(ReKi)                :: det  ! the determinant
    Real(ReKi)                :: invdet  ! inverse of the determinant

    det = M(1, 1) * (M(2, 2) * M(3, 3) - M(3, 2) * M(2, 3)) - &
             M(1, 2) * (M(2, 1) * M(3, 3) - M(2, 3) * M(3, 1)) + &
             M(1, 3) * (M(2, 1) * M(3, 2) - M(2, 2) * M(3, 1));

    invdet = 1 / det   ! because multiplying is faster than dividing

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



END MODULE MoorDyn
