!**********************************************************************************************************************************
! WLaCava (WGL), Matt Lackner (MAL),  Meghan Glade (MEG), and Semyung Park (SP)
! Tuned Mass Damper Module
!**********************************************************************************************************************************
MODULE StrucCtrl

   USE StrucCtrl_Types
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE


   TYPE(ProgDesc), PARAMETER            :: StC_Ver = ProgDesc( 'StrucCtrl', '', '' )




      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: StC_Init                           ! Initialization routine
   PUBLIC :: StC_End                            ! Ending routine (includes clean up)

   PUBLIC :: StC_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                    !   continuous states, and updating discrete states
   PUBLIC :: StC_CalcOutput                     ! Routine for computing outputs

  ! PUBLIC :: StC_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: StC_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states

   !PUBLIC :: StC_UpdateDiscState                ! Tight coupling routine for updating discrete states

   !PUBLIC :: StC_JacobianPInput                 ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the inputs (u)
   !PUBLIC :: StC_JacobianPContState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the continuous
   !                                                 !   states (x)
   !PUBLIC :: StC_JacobianPDiscState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
   !                                                 !   (Xd), and constraint-state (Z) equations all with respect to the discrete
   !                                                 !   states (xd)
   !PUBLIC :: StC_JacobianPConstrState           ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete-
                                                    !   (Xd), and constraint-state (Z) equations all with respect to the constraint
                                                    !   states (z)


   INTEGER(IntKi), PRIVATE, PARAMETER :: ControlMode_NONE      = 0          !< The (StC-universal) control code for not using a particular type of control

   INTEGER(IntKi), PRIVATE, PARAMETER :: DOFMode_Indept        = 1          !< independent DOFs
   INTEGER(IntKi), PRIVATE, PARAMETER :: DOFMode_Omni          = 2          !< omni-directional
   INTEGER(IntKi), PRIVATE, PARAMETER :: DOFMode_TLCD          = 3          !< tuned liquid column dampers !MEG & SP
   INTEGER(IntKi), PRIVATE, PARAMETER :: DOFMode_Prescribed    = 4          !< prescribed force series
   INTEGER(IntKi), PRIVATE, PARAMETER :: DOFMode_ForceDLL      = 5          !< prescribed force series

   INTEGER(IntKi), PRIVATE, PARAMETER :: CMODE_Semi            = 1          !< semi-active control
   INTEGER(IntKi), PRIVATE, PARAMETER :: CMODE_ActiveEXTERN    = 4          !< active control
   INTEGER(IntKi), PRIVATE, PARAMETER :: CMODE_ActiveDLL       = 5          !< active control

   INTEGER(IntKi), PRIVATE, PARAMETER :: SA_CMODE_GH_vel       = 1          !< 1: velocity-based ground hook control;
   INTEGER(IntKi), PRIVATE, PARAMETER :: SA_CMODE_GH_invVel    = 2          !< 2: Inverse velocity-based ground hook control
   INTEGER(IntKi), PRIVATE, PARAMETER :: SA_CMODE_GH_disp      = 3          !< 3: displacement-based ground hook control
   INTEGER(IntKi), PRIVATE, PARAMETER :: SA_CMODE_Ph_FF        = 4          !< 4: Phase difference Algorithm with Friction Force
   INTEGER(IntKi), PRIVATE, PARAMETER :: SA_CMODE_Ph_DF        = 5          !< 5: Phase difference Algorithm with Damping Force

   integer(IntKi), private, parameter :: PRESCRIBED_FORCE_GLOBAL  = 1_IntKi   !< Prescribed forces are in global coords
   integer(IntKi), private, parameter :: PRESCRIBED_FORCE_LOCAL   = 2_IntKi   !< Prescribed forces are in local  coords

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
SUBROUTINE StC_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(StC_InitInputType),       INTENT(INOUT)  :: InitInp     !< Input data for initialization routine.
      TYPE(StC_InputType),           INTENT(  OUT)  :: u           !< An initial guess for the input; input mesh must be defined
      TYPE(StC_ParameterType),       INTENT(  OUT)  :: p           !< Parameters
      TYPE(StC_ContinuousStateType), INTENT(  OUT)  :: x           !< Initial continuous states
      TYPE(StC_DiscreteStateType),   INTENT(  OUT)  :: xd          !< Initial discrete states
      TYPE(StC_ConstraintStateType), INTENT(  OUT)  :: z           !< Initial guess of the constraint states
      TYPE(StC_OtherStateType),      INTENT(  OUT)  :: OtherState  !< Initial other states
      TYPE(StC_OutputType),          INTENT(INOUT)  :: y           !< Initial system outputs (outputs are not calculated;
                                                                   !!   only the output mesh is initialized)
      TYPE(StC_MiscVarType),         INTENT(  OUT)  :: m           !< Misc (optimization) variables
      REAL(DbKi),                    INTENT(INOUT)  :: Interval    !< Coupling interval in seconds: the rate that
                                                                   !!   (1) StC_UpdateStates() is called in loose coupling &
                                                                   !!   (2) StC_UpdateDiscState() is called in tight coupling.
                                                                   !!   Input is the suggested time from the glue code;
                                                                   !!   Output is the actual coupling interval that will be used
                                                                   !!   by the glue code.
      TYPE(StC_InitOutputType),      INTENT(  OUT)  :: InitOut     !< Output for initialization routine
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


         ! Local variables
      INTEGER(IntKi)                                :: NumOuts
      TYPE(StC_InputFile)                           :: InputFileData ! Data stored in the module's input file
      INTEGER(IntKi)                                :: i_pt          ! Generic counter for mesh point
      INTEGER(IntKi)                                :: i             ! Generic counter for mesh point
      REAL(ReKi), allocatable, dimension(:,:)       :: RefPosGlobal

      type(FileInfoType)                            :: FileInfo_In               !< The derived type for holding the full input file for parsing -- we may pass this in the future
      type(FileInfoType)                            :: FileInfo_In_PrescribeFrc  !< The derived type for holding the prescribed forces input file for parsing -- we may pass this in the future
      character(1024)                               :: PriPath       !< Primary path
      integer(IntKi)                                :: UnEcho
      INTEGER(IntKi)                                :: ErrStat2      ! local error status
      CHARACTER(ErrMsgLen)                          :: ErrMsg2       ! local error message

      CHARACTER(*), PARAMETER                       :: RoutineName = 'StC_Init'

         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ''
      NumOuts = 0
      UnEcho  = -1   ! will be > 0 if echo file is opened

     ! Initialize the NWTC Subroutine Library
   CALL NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information
   CALL DispNVD( StC_Ver )
    !............................................................................................
    ! Read the input file and validate the data
    !............................................................................................

   CALL GetPath( InitInp%InputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.

   if (InitInp%UseInputFile) then
      ! Read the entire input file, minus any comment lines, into the FileInfo_In
      ! data structure in memory for further processing.
      call ProcessComFile( InitInp%InputFile, FileInfo_In, ErrStat2, ErrMsg2 )
   else
         ! put passed string info into the FileInfo_In -- FileInfo structure
      call NWTC_Library_CopyFileInfoType( InitInp%PassedPrimaryInputData, FileInfo_In, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
   endif
   if (Failed())  return;

   ! For diagnostic purposes, the following can be used to display the contents
   ! of the FileInfo_In data structure.
   ! call Print_FileInfo_Struct( CU, FileInfo_In ) ! CU is the screen -- different number on different systems.

      !  Parse the FileInfo_In structure of data from the inputfile into the InitInp%InputFile structure
   CALL StC_ParseInputFileInfo( PriPath, InitInp%InputFile, TRIM(InitInp%RootName), InitInp%NumMeshPts, FileInfo_In, InputFileData, UnEcho, ErrStat2, ErrMsg2 )
   if (Failed())  return;

      ! Using the InputFileData structure, check that it makes sense
   CALL StC_ValidatePrimaryData( InputFileData, InitInp, ErrStat2, ErrMsg2 )
   if (Failed())  return;

      ! read in the prescribed forces file
   if ( InputFileData%StC_DOF_MODE == DOFMode_Prescribed ) then
      if (InitInp%UseInputFile_PrescribeFrc) then
         ! Read the entire input file, minus any comment lines, into the FileInfo_In
         ! data structure in memory for further processing.
         call ProcessComFile( InputFileData%PrescribedForcesFile, FileInfo_In_PrescribeFrc, ErrStat2, ErrMsg2 )
      else
            ! put passed string info into the FileInfo_In -- FileInfo structure
         call NWTC_Library_CopyFileInfoType( InitInp%PassedPrescribeFrcData, FileInfo_In_PrescribeFrc, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
      endif
      if (Failed())  return;
      ! For diagnostic purposes, the following can be used to display the contents
      ! of the FileInfo_In data structure.
      !call Print_FileInfo_Struct( CU, FileInfo_In_PrescribeFrc ) ! CU is the screen -- different number on different systems.
      !  Parse the FileInfo_In_PrescribeFrc structure of data from the inputfile into the InitInp%InputFile structure
      CALL StC_ParseTimeSeriesFileInfo( InputFileData%PrescribedForcesFile, FileInfo_In_PrescribeFrc, InputFileData, UnEcho, ErrStat2, ErrMsg2 )
      if (Failed())  return;
   endif

      !............................................................................................
      ! Define parameters here:
      !............................................................................................
   CALL StC_SetParameters( InputFileData, InitInp, p, Interval, ErrStat2, ErrMsg2 )
   if (Failed())  return;

      !............................................................................................
      ! Define initial system states here:
      !............................................................................................

   xd%DummyDiscState = 0
   z%DummyConstrState = 0

   ! Initialize other states here:
   OtherState%DummyOtherState = 0

   call Init_Misc( p, m, ErrStat2, ErrMsg2 )
   if (Failed())  return;

   ! Allocate continuous states (x)
   call AllocAry(x%StC_x, 6, p%NumMeshPts, 'x%StC_x',  ErrStat2,ErrMsg2)
   if (Failed())  return;

   ! Define initial guess for the system states here:
   do i_pt=1,p%NumMeshPts
      x%StC_x(1,i_pt) = InputFileData%StC_X_DSP
      x%StC_x(2,i_pt) = 0
      x%StC_x(3,i_pt) = InputFileData%StC_Y_DSP
      x%StC_x(4,i_pt) = 0
      if ((p%StC_DOF_MODE == DOFMode_Indept) .and. p%StC_Z_DOF) then    ! Should be zero for omni and TLCD
         x%StC_x(5,i_pt) = InputFileData%StC_Z_DSP
      else
         x%StC_x(5,i_pt) = 0.0_ReKi
      endif
      x%StC_x(6,i_pt) = 0
   enddo


   ! set positions and orientations for tuned mass dampers's
   call AllocAry(InitOut%RelPosition,  3, p%NumMeshPts, 'RelPosition',     ErrStat2,ErrMsg2);  if (Failed())  return;
   call AllocAry(RefPosGlobal,         3, p%NumMeshPts, 'RefPosGlobal',    ErrStat2,ErrMsg2);  if (Failed())  return;

   ! Set the initial positions and orientations for each point (Ref coords)
   do i_pt = 1,p%NumMeshPts
      InitOut%RelPosition(:,i_pt)   = (/ InputFileData%StC_P_X, InputFileData%StC_P_Y, InputFileData%StC_P_Z /)
      RefPosGlobal(:,i_pt)          = InitInp%InitRefPos(:,i_pt) + real( matmul(InitOut%RelPosition(:,i_pt),InitInp%InitRefOrient(:,:,i_pt)), ReKi)
   enddo

    ! Define system output initializations (set up mesh) here:
    ! Create the input and output meshes associated with lumped loads

   ALLOCATE (u%Mesh(p%NumMeshPts), STAT=ErrStat2)
   IF (ErrStat2/=0) THEN
      CALL SetErrStat(ErrID_Fatal,"Error allocating u%Mesh.",ErrStat,ErrMsg,RoutineName)
      CALL Cleanup()
      RETURN
   END IF
   ALLOCATE (y%Mesh(p%NumMeshPts), STAT=ErrStat2)
   IF (ErrStat2/=0) THEN
      CALL SetErrStat(ErrID_Fatal,"Error allocating y%Mesh.",ErrStat,ErrMsg,RoutineName)
      CALL Cleanup()
      RETURN
   END IF

   ! Create Mesh(i_pt)
   DO i_pt = 1,p%NumMeshPts

      CALL MeshCreate( BlankMesh        = u%Mesh(i_pt)      &
                     ,IOS               = COMPONENT_INPUT   &
                     ,Nnodes            = 1                 &
                     ,ErrStat           = ErrStat2          &
                     ,ErrMess           = ErrMsg2           &
                     ,TranslationDisp   = .TRUE.            &
                     ,Orientation       = .TRUE.            &
                     ,TranslationVel    = .TRUE.            &
                     ,RotationVel       = .TRUE.            &
                     ,TranslationAcc    = .TRUE.            &
                     ,RotationAcc       = .TRUE.)
      if (Failed())  return;


         ! Create the node on the mesh
      CALL MeshPositionNode ( u%Mesh(i_pt),1, RefPosGlobal(:,i_pt), ErrStat2, ErrMsg2, InitInp%InitRefOrient(:,:,i_pt) )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         ! Create the mesh element
      CALL MeshConstructElement (  u%Mesh(i_pt)        &
                                  , ELEMENT_POINT      &
                                  , ErrStat2           &
                                  , ErrMsg2            &
                                  , 1                  )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      CALL MeshCommit ( u%Mesh(i_pt)        &
                      , ErrStat2            &
                      , ErrMsg2             )
      if (Failed())  return;

      CALL MeshCopy ( SrcMesh      = u%Mesh(i_pt)           &
                     ,DestMesh     = y%Mesh(i_pt)           &
                     ,CtrlCode     = MESH_SIBLING           &
                     ,IOS          = COMPONENT_OUTPUT       &
                     ,ErrStat      = ErrStat2               &
                     ,ErrMess      = ErrMsg2                &
                     ,Force        = .TRUE.                 &
                     ,Moment       = .TRUE.                 )

      if (Failed())  return;

      u%Mesh(i_pt)%RemapFlag  = .TRUE.
      y%Mesh(i_pt)%RemapFlag  = .TRUE.

      ! Set initial displacements
      u%Mesh(i_pt)%Orientation(1:3,1:3,1) = InitInp%InitOrient(:,:,i_pt)
      u%Mesh(i_pt)%TranslationDisp(1:3,1) = InitInp%InitTransDisp(:,i_pt)
   enddo


   !bjj: removed for now; output handled in ServoDyn
    !IF (NumOuts > 0) THEN
    !   ALLOCATE( y%WriteOutput(NumOuts), STAT = ErrStat )
    !   IF ( ErrStat/= 0 ) THEN
    !      CALL SetErrStat(ErrID_Fatal,'Error allocating output array.',ErrStat,ErrMsg,'StC_Init')
    !      CALL Cleanup()
    !      RETURN
    !   END IF
    !   y%WriteOutput = 0
    !
    !   ! Define initialization-routine output here:
    !   ALLOCATE( InitOut%WriteOutputHdr(NumOuts), InitOut%WriteOutputUnt(NumOuts), STAT = ErrStat )
    !   IF ( ErrStat/= 0 ) THEN
    !      CALL SetErrStat(ErrID_Fatal,'Error allocating output header and units arrays.',ErrStat,ErrMsg,'StC_Init')
    !      CALL Cleanup()
    !      RETURN
    !   END IF
    !
    !   DO i=1,NumOuts
    !        InitOut%WriteOutputHdr(i) = "Heading"//trim(num2lstr(i))
    !        InitOut%WriteOutputUnt(i) = "(-)"
    !   END DO
    !
    !END IF

    !bjj: need to initialize headers/units

      ! Set the interval value to tell ServoDyn we are using (we don't actually change this in StC)
   Interval = p%DT


      ! Initialize the input and output arrays for control channels
      !     NOTE: these will get resized later in ServoDyn!!!!
   if (maxval(p%StC_CChan) > 0) then
      call AllocAry( u%CmdStiff, 3, maxval(p%StC_CChan), 'u%CmdStiff', ErrStat2, ErrMsg2 )
         if (Failed())  return;
      call AllocAry( u%CmdDamp,  3, maxval(p%StC_CChan), 'u%CmdDamp',  ErrStat2, ErrMsg2 )
         if (Failed())  return;
      call AllocAry( u%CmdBrake, 3, maxval(p%StC_CChan), 'u%CmdBrake', ErrStat2, ErrMsg2 )
         if (Failed())  return;
      call AllocAry( u%CmdForce, 3, maxval(p%StC_CChan), 'u%CmdForce', ErrStat2, ErrMsg2 )
         if (Failed())  return;
      call AllocAry( y%MeasDisp, 3, maxval(p%StC_CChan), 'y%MeasDisp', ErrStat2, ErrMsg2 )
         if (Failed())  return;
      call AllocAry( y%MeasVel,  3, maxval(p%StC_CChan), 'y%MeasVel',  ErrStat2, ErrMsg2 )
         if (Failed())  return;
      ! Initialize to zero (if we asked for channel 9, the first 8 channel entries are zero)
      u%CmdStiff  =  0.0_ReKi
      u%CmdDamp   =  0.0_ReKi
      u%CmdBrake  =  0.0_ReKi
      u%CmdForce  =  0.0_ReKi
      y%MeasDisp  =  0.0_ReKi
      y%MeasVel   =  0.0_ReKi
      ! Check that dimensions of x are what we expect
      if (size(p%StC_CChan) /= size(x%StC_x,2)) then
         ErrStat2 = ErrID_Fatal
         ErrMsg2  = "Error in setup of StC_Chan array -- does not match expected dimensionality of x%StC_x"
         if (Failed())  return;
      endif
      ! Set actual values for channels requested
      do i=1,size(p%StC_CChan)
         if (p%StC_CChan(i) > 0) then
            u%CmdStiff(1:3,p%StC_CChan(i))  = (/ p%K_X, p%K_Y, p%K_Z /)
            u%CmdDamp( 1:3,p%StC_CChan(i))  = (/ p%C_X, p%C_Y, p%C_Z /)
            !u%CmdBrake and u%CmdForce--- leave these at zero for now (no input file method to set it)
            !  The states are sized by (6,NumMeshPts).  NumMeshPts is then used to set
            !  size of StC_CChan as well.  For safety, we will check it here.
            y%MeasDisp(1:3,p%StC_CChan(i))  = (/ x%StC_x(1,i), x%StC_x(3,i), x%StC_x(5,i) /)
            y%MeasVel( 1:3,p%StC_CChan(i))  = (/ x%StC_x(2,i), x%StC_x(4,i), x%StC_x(6,i) /)
         endif
      enddo
   endif

   call cleanup()
!................................
CONTAINS
   subroutine Init_Misc( p, m, ErrStat, ErrMsg )
      type(StC_ParameterType),intent(in   )  :: p        !< Parameters
      type(StC_MiscVarType),  intent(inout)  :: m        !< Misc (optimization) variables
      integer(IntKi),         intent(  out) :: ErrStat   ! The error identifier (ErrStat)
      character(ErrMsgLen),   intent(  out) :: ErrMsg    ! The error message (ErrMsg)

      !  Accelerations, velocities, and resultant forces -- used in all tuned mass calcs (so we don't reallocate all the time)
      !  Note: these variables had been allocated multiple places before and sometimes passed between routines. So
      !        they have been moved into MiscVars so that we don so we don't reallocate all the time
      call AllocAry(m%a_G    , 3, p%NumMeshPts,'a_G'     , ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;
      call AllocAry(m%rdisp_P, 3, p%NumMeshPts,'rdisp_P' , ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;
      call AllocAry(m%rdot_P , 3, p%NumMeshPts,'rdot_P'  , ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;
      call AllocAry(m%rddot_P, 3, p%NumMeshPts,'rddot_P' , ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;
      call AllocAry(m%omega_P, 3, p%NumMeshPts,'omega_P' , ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;
      call AllocAry(m%alpha_P, 3, p%NumMeshPts,'alpha_P' , ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;
      call AllocAry(m%Acc    , 3, p%NumMeshPts,'Acc'     , ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;    ! Summed accelerations
      !  Note: the following two were added to misc so that we have the option of outputting the forces and moments
      !        from each tuned mass system at some later point
      call AllocAry(m%F_P    , 3, p%NumMeshPts,'F_P'     , ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;
      call AllocAry(m%M_P    , 3, p%NumMeshPts,'M_P'     , ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;

      !  External and stop forces
      !  Note: these variables had been allocated multiple places before and sometimes passed between routines. So
      !        they have been moved into MiscVars so that we don so we don't reallocate all the time.
      call AllocAry(m%F_stop , 3, p%NumMeshPts, 'F_stop' , ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;  m%F_stop  = 0.0_ReKi
      call AllocAry(m%F_ext  , 3, p%NumMeshPts, 'F_ext'  , ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;  m%F_ext   = 0.0_ReKi
      call AllocAry(m%F_fr   , 3, p%NumMeshPts, 'F_fr'   , ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;  m%F_fr    = 0.0_ReKi
      call AllocAry(m%C_ctrl , 3, p%NumMeshPts, 'C_ctrl' , ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;  m%C_ctrl  = 0.0_ReKi
      call AllocAry(m%C_Brake, 3, p%NumMeshPts, 'C_Brake', ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;  m%C_Brake = 0.0_ReKi
      call AllocAry(m%F_table, 3, p%NumMeshPts, 'F_table', ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;  m%F_table = 0.0_ReKi
      call AllocAry(m%F_k    , 3, p%NumMeshPts, 'F_k'    , ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return;  m%F_k     = 0.0_ReKi

      ! Set spring constants to value from input file
      call AllocAry(m%K      , 3, p%NumMeshPts, 'K'      , ErrStat, ErrMsg);  if (ErrStat >= AbortErrLev) return
      do i_pt=1,p%NumMeshPts
         m%K(1,i_pt) = p%K_X
         m%K(2,i_pt) = p%K_Y
         m%K(3,i_pt) = p%K_Z
      enddo

      ! indexing
      m%PrescribedInterpIdx = 0_IntKi ! index tracker for PrescribedForce option

   end subroutine Init_Misc
   !.........................................
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'StC_Init' )
      Failed = ErrStat >= AbortErrLev
      if (Failed)    call cleanup()
   end function Failed
   !.........................................
   SUBROUTINE cleanup()
      if (UnEcho > 0)                     close(UnEcho)                    ! Close echo file
      if (allocated(RefPosGlobal    ))    deallocate(RefPosGlobal    )
      CALL StC_DestroyInputFile( InputFileData, ErrStat2, ErrMsg2)      ! Ignore warnings here.
   END SUBROUTINE cleanup
!.........................................
END SUBROUTINE StC_Init
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
SUBROUTINE StC_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      TYPE(StC_InputType),           INTENT(INOUT)  :: u           !< System inputs
      TYPE(StC_ParameterType),       INTENT(INOUT)  :: p           !< Parameters
      TYPE(StC_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states
      TYPE(StC_DiscreteStateType),   INTENT(INOUT)  :: xd          !< Discrete states
      TYPE(StC_ConstraintStateType), INTENT(INOUT)  :: z           !< Constraint states
      TYPE(StC_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states
      TYPE(StC_OutputType),          INTENT(INOUT)  :: y           !< System outputs
      TYPE(StC_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Place any last minute operations or calculations here:



         ! Write the StrucCtrl-level output file data if the user requested module-level output
         ! and the current time has advanced since the last stored time step.



         ! Close files here:


         ! Destroy the input data:

      CALL StC_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:

      CALL StC_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:

      CALL StC_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL StC_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL StC_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL StC_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )

      CALL StC_DestroyMisc(  m,  ErrStat, ErrMsg )

         ! Destroy the output data:

      CALL StC_DestroyOutput( y, ErrStat, ErrMsg )

END SUBROUTINE StC_End
!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving constraint states, integrating continuous states, and updating discrete states.
!! Continuous, constraint, and discrete states are updated to values at t + Interval.
SUBROUTINE StC_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   )  :: t               !< Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   )  :: n               !< Current step of the simulation: t = n*Interval
      TYPE(StC_InputType),                INTENT(INOUT)  :: Inputs(:)       !< Inputs at InputTimes
      REAL(DbKi),                         INTENT(IN   )  :: InputTimes(:)   !< Times in seconds associated with Inputs
      TYPE(StC_ParameterType),            INTENT(IN   )  :: p               !< Parameters
      TYPE(StC_ContinuousStateType),      INTENT(INOUT)  :: x               !< Input: Continuous states at t;
                                                                            !!   Output: Continuous states at t + Interval
      TYPE(StC_DiscreteStateType),        INTENT(INOUT)  :: xd              !< Input: Discrete states at t;
                                                                            !!   Output: Discrete states at t + Interval
      TYPE(StC_ConstraintStateType),      INTENT(INOUT)  :: z               !< Input: Constraint states at t;
                                                                            !!   Output: Constraint states at t + Interval
      TYPE(StC_OtherStateType),           INTENT(INOUT)  :: OtherState      !< Input: Other states at t;
                                                                            !!   Output: Other states at t + Interval
      TYPE(StC_MiscVarType),              INTENT(INOUT)  :: m               !< Misc (optimization) variables
      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat         !< Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg          !< Error message if ErrStat /= ErrID_None

         ! Local variables
      !INTEGER                                            :: I               ! Generic loop counter
      !TYPE(StC_ContinuousStateType)                      :: dxdt            ! Continuous state derivatives at t
      !TYPE(StC_DiscreteStateType)                        :: xd_t            ! Discrete states at t (copy)
      !TYPE(StC_ConstraintStateType)                      :: z_Residual      ! Residual of the constraint state functions (Z)
      !TYPE(StC_InputType)                                :: u               ! Instantaneous inputs
      !INTEGER(IntKi)                                     :: ErrStat2        ! Error status of the operation (secondary error)
      !CHARACTER(ErrMsgLen)                               :: ErrMsg2         ! Error message if ErrStat2 /= ErrID_None
      !INTEGER                                            :: nTime           ! number of inputs


      CALL StC_RK4( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

END SUBROUTINE StC_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!!
!!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x).
!!   Define constants k1, k2, k3, and k4 as
!!        k1 = dt * f(t        , x_t        )
!!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!!        k4 = dt * f(t + dt   , x_t + k3   ).
!!   Then the continuous states at t = t + dt are
!!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!!
!! For details, see:
!! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for
!!   Runge-Kutta." Sections 16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England:
!!   Cambridge University Press, pp. 704-716, 1992.
SUBROUTINE StC_RK4( t, n, u, utimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                    INTENT(IN   )  :: t           !< Current simulation time in seconds
      INTEGER(IntKi),                INTENT(IN   )  :: n           !< time step number
      TYPE(StC_InputType),           INTENT(INOUT)  :: u(:)        !< Inputs at t (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                    INTENT(IN   )  :: utimes(:)   !< times of input
      TYPE(StC_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(StC_ContinuousStateType), INTENT(INOUT)  :: x           !< Continuous states at t on input at t + dt on output
      TYPE(StC_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
      TYPE(StC_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t (possibly a guess)
      TYPE(StC_OtherStateType),      INTENT(INOUT)  :: OtherState  !< Other states at t
      TYPE(StC_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(StC_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states
      TYPE(StC_ContinuousStateType)                 :: k1          ! RK4 constant; see above
      TYPE(StC_ContinuousStateType)                 :: k2          ! RK4 constant; see above
      TYPE(StC_ContinuousStateType)                 :: k3          ! RK4 constant; see above
      TYPE(StC_ContinuousStateType)                 :: k4          ! RK4 constant; see above

      TYPE(StC_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
      TYPE(StC_InputType)                           :: u_interp    ! interpolated value of inputs
      integer(IntKi)                                :: i_pt        ! Generic counter for mesh point

      INTEGER(IntKi)                                :: ErrStat2    ! local error status
      CHARACTER(ErrMsgLen)                          :: ErrMsg2     ! local error message (ErrMsg)


         ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

      ! if prescribed forces, there are no states to advance, so return
      if ( p%StC_DOF_MODE == DOFMode_Prescribed ) then
         do i_pt=1,p%NumMeshPts
            x%StC_x(1,i_pt) = 0
            x%StC_x(2,i_pt) = 0
            x%StC_x(3,i_pt) = 0
            x%StC_x(4,i_pt) = 0
            x%StC_x(5,i_pt) = 0
            x%StC_x(6,i_pt) = 0
         enddo
         return
      endif

      CALL StC_CopyContState( x, k1, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL StC_CopyContState( x, k2, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL StC_CopyContState( x, k3, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL StC_CopyContState( x, k4,    MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
      CALL StC_CopyContState( x, x_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      CALL StC_CopyInput( u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ! interpolate u to find u_interp = u(t)
      CALL StC_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ! find xdot at t
      CALL StC_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      do i_pt=1,p%NumMeshPts
         k1%StC_x(:,i_pt)     = p%dt * xdot%StC_x(:,i_pt)
         x_tmp%StC_x(:,i_pt)  = x%StC_x(:,i_pt)  + 0.5 * k1%StC_x(:,i_pt)
      enddo


      ! interpolate u to find u_interp = u(t + dt/2)
      CALL StC_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%dt, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ! find xdot at t + dt/2
      CALL StC_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      do i_pt=1,p%NumMeshPts
         k2%StC_x(:,i_pt)     = p%dt * xdot%StC_x(:,i_pt)
         x_tmp%StC_x(:,i_pt)  = x%StC_x(:,i_pt)  + 0.5 * k2%StC_x(:,i_pt)
      enddo


      ! find xdot at t + dt/2
      CALL StC_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      do i_pt=1,p%NumMeshPts
         k3%StC_x(:,i_pt)     = p%dt * xdot%StC_x(:,i_pt)
         x_tmp%StC_x(:,i_pt)  = x%StC_x(:,i_pt)  + k3%StC_x(:,i_pt)
      enddo


      ! interpolate u to find u_interp = u(t + dt)
      CALL StC_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat2, ErrMsg2)
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      ! find xdot at t + dt
      CALL StC_CalcContStateDeriv( t + p%dt, u_interp, p, x_tmp, xd, z, OtherState, m, xdot, ErrStat2, ErrMsg2 )
         CALL CheckError(ErrStat2,ErrMsg2)
         IF ( ErrStat >= AbortErrLev ) RETURN

      do i_pt=1,p%NumMeshPts
         k4%StC_x(:,i_pt) = p%dt * xdot%StC_x(:,i_pt)
         x%StC_x(:,i_pt)   = x%StC_x(:,i_pt)  +  ( k1%StC_x(:,i_pt)  + 2. * k2%StC_x(:,i_pt)  + 2. * k3%StC_x(:,i_pt)  + k4%StC_x(:,i_pt)  ) / 6.
         ! x%StC_dxdt = x%StC_dxdt +  ( k1%StC_dxdt + 2. * k2%StC_dxdt + 2. * k3%StC_dxdt + k4%StC_dxdt ) / 6.
      enddo

         ! clean up local variables:
      CALL ExitThisRoutine(  )

CONTAINS
   !...............................................................................................................................
   SUBROUTINE ExitThisRoutine()
   ! This subroutine destroys all the local variables
   !...............................................................................................................................

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)


      CALL StC_DestroyContState( xdot,     ErrStat3, ErrMsg3 )
      CALL StC_DestroyContState( k1,       ErrStat3, ErrMsg3 )
      CALL StC_DestroyContState( k2,       ErrStat3, ErrMsg3 )
      CALL StC_DestroyContState( k3,       ErrStat3, ErrMsg3 )
      CALL StC_DestroyContState( k4,       ErrStat3, ErrMsg3 )
      CALL StC_DestroyContState( x_tmp,    ErrStat3, ErrMsg3 )

      CALL StC_DestroyInput(     u_interp, ErrStat3, ErrMsg3 )

   END SUBROUTINE ExitThisRoutine
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         ! local variables
      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(ErrMsgLen)       :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//'StC_RK4:'//TRIM(Msg)
         ErrStat = MAX(ErrStat,ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................

         IF ( ErrStat >= AbortErrLev ) CALL ExitThisRoutine( )


      END IF

   END SUBROUTINE CheckError

END SUBROUTINE StC_RK4
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
SUBROUTINE StC_CalcOutput( Time, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                    INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(StC_InputType),           INTENT(IN   )  :: u           !< Inputs at Time
      TYPE(StC_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(StC_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(StC_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(StC_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(StC_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time
      TYPE(StC_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at Time (Input only so that mesh con-
                                                                   !!  nectivity information does not have to be recalculated)
      TYPE(StC_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      !  local variables for force calcualtions in X-DOF, Y-DOF, and XY-DOF
      real(ReKi), dimension(3)   :: F_X_P
      real(ReKi), dimension(3)   :: F_Y_P
      real(ReKi), dimension(3)   :: F_Z_P
      real(ReKi), dimension(3)   :: F_XY_P

      !  NOTE: the following two sets of variables could likely be combined into arrays
      !        that could be more easily used with array functions like MATMUL, cross_product,
      !        dot_product etc.
      ! Fore-aft TLCD reactionary forces !MEG & SP
      Real(ReKi)                 :: F_x_tlcd_WR_N
      Real(ReKi)                 :: F_y_tlcd_WR_N
      Real(ReKi)                 :: F_x_tlcd_WL_N
      Real(ReKi)                 :: F_y_tlcd_WL_N
      Real(ReKi)                 :: F_y_tlcd_WH_N
      Real(ReKi)                 :: F_z_tlcd_WH_N

      ! Side-side orthogonal TLCD reactionary forces !MEG & SP
      Real(ReKi)                 :: F_x_otlcd_WB_N
      Real(ReKi)                 :: F_y_otlcd_WB_N
      Real(ReKi)                 :: F_x_otlcd_WF_N
      Real(ReKi)                 :: F_y_otlcd_WF_N
      Real(ReKi)                 :: F_x_otlcd_WH_N
      Real(ReKi)                 :: F_z_otlcd_WH_N

      TYPE(StC_ContinuousStateType)              :: dxdt    ! first time derivative of continuous states

      integer(IntKi)       :: i,j         !< generic counter
      integer(IntKi)       :: i_pt        ! Generic counter for mesh point

      ! Local error handling
      integer(IntKi)       :: ErrStat2
      character(ErrMsgLen) :: ErrMsg2


      ErrStat = ErrID_None
      ErrMsg  = ""


      ! Compute accelerations and velocities in local coordinates
      do i_pt=1,p%NumMeshPts
         m%a_G(:,i_pt)     = matmul(u%Mesh(i_pt)%Orientation(:,:,1),p%Gravity)
         m%rdisp_P(:,i_pt) = matmul(u%Mesh(i_pt)%Orientation(:,:,1),u%Mesh(i_pt)%TranslationDisp(:,1))   ! for ground StC_GroundHookDamp
         m%rdot_P(:,i_pt)  = matmul(u%Mesh(i_pt)%Orientation(:,:,1),u%Mesh(i_pt)%TranslationVel(:,1))    ! for ground StC_GroundHookDamp
         m%rddot_P(:,i_pt) = matmul(u%Mesh(i_pt)%Orientation(:,:,1),u%Mesh(i_pt)%TranslationAcc(:,1))
         m%omega_P(:,i_pt) = matmul(u%Mesh(i_pt)%Orientation(:,:,1),u%Mesh(i_pt)%RotationVel(:,1))
         m%alpha_P(:,i_pt) = matmul(u%Mesh(i_pt)%Orientation(:,:,1),u%Mesh(i_pt)%RotationAcc(:,1))
      enddo


         ! calculate the derivative, only to get updated values of m, which are used in the equations below
      CALL StC_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat2, ErrMsg2 ); if (Failed()) return;


      IF (p%StC_DOF_MODE == ControlMode_None) THEN
         do i_pt=1,p%NumMeshPts
            y%Mesh(i_pt)%Force(:,1)  = 0.0_ReKi
            y%Mesh(i_pt)%Moment(:,1) = 0.0_ReKi
            m%F_P(1:3,i_pt)          = 0.0_ReKi
            m%M_P(1:3,i_pt)          = 0.0_ReKi
         enddo
      ELSEIF (p%StC_DOF_MODE == DOFMode_Indept) THEN

         ! StrucCtrl external forces of dependent degrees:
         do i_pt=1,p%NumMeshPts
            F_X_P(2) = - p%M_X * ( m%a_G(2,i_pt) - m%rddot_P(2,i_pt) - (m%alpha_P(3,i_pt) + m%omega_P(1,i_pt)*m%omega_P(2,i_pt))*x%StC_x(1,i_pt) - 2*m%omega_P(3,i_pt)*x%StC_x(2,i_pt) )
            F_X_P(3) = - p%M_X * ( m%a_G(3,i_pt) - m%rddot_P(3,i_pt) + (m%alpha_P(2,i_pt) - m%omega_P(1,i_pt)*m%omega_P(3,i_pt))*x%StC_x(1,i_pt) + 2*m%omega_P(2,i_pt)*x%StC_x(2,i_pt) )

            F_Y_P(1) = - p%M_Y * ( m%a_G(1,i_pt) - m%rddot_P(1,i_pt) + (m%alpha_P(3,i_pt) - m%omega_P(1,i_pt)*m%omega_P(2,i_pt))*x%StC_x(3,i_pt) + 2*m%omega_P(3,i_pt)*x%StC_x(4,i_pt) )
            F_Y_P(3) = - p%M_Y * ( m%a_G(3,i_pt) - m%rddot_P(3,i_pt) - (m%alpha_P(1,i_pt) + m%omega_P(2,i_pt)*m%omega_P(3,i_pt))*x%StC_x(3,i_pt) - 2*m%omega_P(1,i_pt)*x%StC_x(4,i_pt) )

            F_Z_P(1) = - p%M_Z * ( m%a_G(1,i_pt) - m%rddot_P(1,i_pt) - (m%alpha_P(2,i_pt) + m%omega_P(1,i_pt)*m%omega_P(3,i_pt))*x%StC_x(5,i_pt) - 2*m%omega_P(2,i_pt)*x%StC_x(6,i_pt) )
            F_Z_P(2) = - p%M_Z * ( m%a_G(2,i_pt) - m%rddot_P(2,i_pt) + (m%alpha_P(1,i_pt) - m%omega_P(2,i_pt)*m%omega_P(3,i_pt))*x%StC_x(5,i_pt) + 2*m%omega_P(1,i_pt)*x%StC_x(6,i_pt) )

            ! inertial contributions from mass of tuned mass dampers and acceleration of point
            ! forces and moments in local coordinates
            m%F_P(1,i_pt) =  m%K(1,i_pt) * x%StC_x(1,i_pt) + m%C_ctrl(1,i_pt) * x%StC_x(2,i_pt) + m%C_Brake(1,i_pt) * x%StC_x(2,i_pt) - m%F_stop(1,i_pt) - m%F_ext(1,i_pt) - m%F_fr(1,i_pt) - F_Y_P(1) - F_Z_P(1) + m%F_table(1,i_pt)
            m%F_P(2,i_pt) =  m%K(2,i_pt) * x%StC_x(3,i_pt) + m%C_ctrl(2,i_pt) * x%StC_x(4,i_pt) + m%C_Brake(2,i_pt) * x%StC_x(4,i_pt) - m%F_stop(2,i_pt) - m%F_ext(2,i_pt) - m%F_fr(2,i_pt) - F_X_P(2) - F_Z_P(2) + m%F_table(2,i_pt)
            m%F_P(3,i_pt) =  m%K(3,i_pt) * x%StC_x(5,i_pt) + m%C_ctrl(3,i_pt) * x%StC_x(6,i_pt) + m%C_Brake(3,i_pt) * x%StC_x(6,i_pt) - m%F_stop(3,i_pt) - m%F_ext(3,i_pt) - m%F_fr(3,i_pt) - F_X_P(3) - F_Y_P(3) + m%F_table(3,i_pt) - p%StC_Z_PreLd

            m%M_P(1,i_pt) =  - F_Y_P(3)  * x%StC_x(3,i_pt)  +  F_Z_P(2) * x%StC_x(5,i_pt)
            m%M_P(2,i_pt) =    F_X_P(3)  * x%StC_x(1,i_pt)  -  F_Z_P(1) * x%StC_x(5,i_pt)
            m%M_P(3,i_pt) =  - F_X_P(2)  * x%StC_x(1,i_pt)  +  F_Y_P(1) * x%StC_x(3,i_pt)    ! NOTE signs match document, but are changed from prior value

            ! forces and moments in global coordinates
            y%Mesh(i_pt)%Force(:,1) =  real(matmul(transpose(u%Mesh(i_pt)%Orientation(:,:,1)),m%F_P(1:3,i_pt)),ReKi)
            y%Mesh(i_pt)%Moment(:,1) = real(matmul(transpose(u%Mesh(i_pt)%Orientation(:,:,1)),m%M_P(1:3,i_pt)),ReKi)
         enddo

      ELSE IF (p%StC_DOF_MODE == DOFMode_Omni) THEN

         !note: m%F_k is computed earlier in StC_CalcContStateDeriv

         ! StrucCtrl external forces of dependent degrees:
         do i_pt=1,p%NumMeshPts
            F_XY_P(1) = 0
            F_XY_P(2) = 0
            F_XY_P(3) = - p%M_XY * (  m%a_G(3,i_pt) - m%rddot_P(3,i_pt)                                                       &
                                                - (m%alpha_P(1,i_pt) + m%omega_P(2,i_pt)*m%omega_P(3,i_pt))*x%StC_x(3,i_pt)   &
                                                + (m%alpha_P(2,i_pt) - m%omega_P(1,i_pt)*m%omega_P(3,i_pt))*x%StC_x(1,i_pt)   &
                                                - 2*m%omega_P(1,i_pt)*x%StC_x(4,i_pt)                                         &
                                                + 2*m%omega_P(2,i_pt)*x%StC_x(2,i_pt)       )

            ! inertial contributions from mass of tuned mass dampers and acceleration of point
            ! forces and moments in local coordinates
            m%F_P(1,i_pt) =  m%K(1,i_pt) * x%StC_x(1,i_pt) + m%C_ctrl(1,i_pt) * x%StC_x(2,i_pt) + m%C_Brake(1,i_pt) * x%StC_x(2,i_pt) - m%F_stop(1,i_pt) - m%F_ext(1,i_pt) - m%F_fr(1,i_pt) - F_XY_P(1) + m%F_table(1,i_pt)*(m%F_k(1,i_pt))
            m%F_P(2,i_pt) =  m%K(2,i_pt) * x%StC_x(3,i_pt) + m%C_ctrl(2,i_pt) * x%StC_x(4,i_pt) + m%C_Brake(2,i_pt) * x%StC_x(4,i_pt) - m%F_stop(2,i_pt) - m%F_ext(2,i_pt) - m%F_fr(2,i_pt) - F_XY_P(2) + m%F_table(2,i_pt)*(m%F_k(2,i_pt))
            m%F_P(3,i_pt) = - F_XY_P(3)

            m%M_P(1,i_pt) = - F_XY_P(3) * x%StC_x(3,i_pt)
            m%M_P(2,i_pt) =   F_XY_P(3) * x%StC_x(1,i_pt)
            m%M_P(3,i_pt) = - F_XY_P(1) * x%StC_x(3,i_pt) + F_XY_P(2) * x%StC_x(1,i_pt)

            ! forces and moments in global coordinates
            y%Mesh(i_pt)%Force(:,1) =  real(matmul(transpose(u%Mesh(i_pt)%Orientation(:,:,1)),m%F_P(1:3,i_pt)),ReKi)
            y%Mesh(i_pt)%Moment(:,1) = real(matmul(transpose(u%Mesh(i_pt)%Orientation(:,:,1)),m%M_P(1:3,i_pt)),ReKi)
         enddo

      ELSE IF (p%StC_DOF_MODE == DOFMode_TLCD) THEN

         do i_pt=1,p%NumMeshPts
            !fore-aft TLCD external forces of dependent degrees
            F_x_tlcd_WR_N = p%rho_X*p%area_X*((p%L_X-p%B_X)/2+x%StC_x(1,i_pt))*(                               &
                                       m%rddot_P(1,i_pt)                                                       &
                                    +2*m%omega_P(2,i_pt)*x%StC_x(2,i_pt)                                       &
                                      +m%alpha_P(2,i_pt)*((p%L_X-p%B_X)/2+x%StC_x(1,i_pt))                     &
                                      -m%omega_P(2,i_pt)*m%omega_P(2,i_pt)*p%B_X*.5                            &
                                      -m%omega_P(3,i_pt)*m%omega_P(3,i_pt)*p%B_X*.5                            &
                                      +m%omega_P(3,i_pt)*m%omega_P(1,i_pt)*((p%L_X-p%B_X)/2+x%StC_x(1,i_pt))   &
                                      -m%a_G(1,i_pt)  )
            F_y_tlcd_WR_N = p%rho_X*p%area_X*((p%L_X-p%B_X)/2+x%StC_x(1,i_pt))*(                               &
                                       m%rddot_P(2,i_pt)                                                       &
                                    -2*m%omega_P(1,i_pt)*x%StC_x(2,i_pt)                                       &
                                      +m%alpha_P(3,i_pt)*p%B_X*.5                                              &
                                      -m%alpha_P(1,i_pt)*((p%L_X-p%B_X)/2+x%StC_x(1,i_pt))                     &
                                      +m%omega_P(3,i_pt)*m%omega_P(2,i_pt)*((p%L_X-p%B_X)/2+x%StC_x(1,i_pt))   &
                                      +m%omega_P(1,i_pt)*m%omega_P(2,i_pt)*p%B_X*.5                            &
                                      -m%a_G(2,i_pt)  )
            F_x_tlcd_WL_N = p%rho_X*p%area_X*((p%L_X-p%B_X)/2-x%StC_x(1,i_pt))*(                               &
                                       m%rddot_P(1,i_pt)                                                       &
                                    -2*m%omega_P(2,i_pt)*x%StC_x(2,i_pt)                                       &
                                      +m%alpha_P(2,i_pt)*((p%L_X-p%B_X)/2-x%StC_x(1,i_pt))                     &
                                      +m%omega_P(2,i_pt)*m%omega_P(2,i_pt)*p%B_X*.5                            &
                                      +m%omega_P(3,i_pt)*m%omega_P(3,i_pt)*p%B_X*.5                            &
                                      +m%omega_P(3,i_pt)*m%omega_P(1,i_pt)*((p%L_X-p%B_X)/2-x%StC_x(1,i_pt))   &
                                      -m%a_G(1,i_pt)  )
            F_y_tlcd_WL_N = p%rho_X*p%area_X*((p%L_X-p%B_X)/2-x%StC_x(1,i_pt))*(                               &
                                       m%rddot_P(2,i_pt)                                                       &
                                    +2*m%omega_P(1,i_pt)*x%StC_x(2,i_pt)                                       &
                                      -m%alpha_P(3,i_pt)*p%B_X*.5                                              &
                                      -m%alpha_P(1,i_pt)*((p%L_X-p%B_X)/2-x%StC_x(1,i_pt))                     &
                                      +m%omega_P(3,i_pt)*m%omega_P(2,i_pt)*((p%L_X-p%B_X)/2-x%StC_x(1,i_pt))   &
                                      -m%omega_P(1,i_pt)*m%omega_P(2,i_pt)*p%B_X*.5                            &
                                      -m%a_G(2,i_pt)  )
            F_y_tlcd_WH_N = p%rho_X*p%area_X/p%area_ratio_X*p%B_X*(                       &
                                       m%rddot_P(2,i_pt)                                  &
                                    +2*m%omega_P(3,i_pt)*p%area_ratio_X*x%StC_x(2,i_pt)   &
                                      -m%a_G(2,i_pt)  )
            F_z_tlcd_WH_N = p%rho_X*p%area_X/p%area_ratio_X*p%B_X*(                       &
                                       m%rddot_P(3,i_pt)                                  &
                                    -2*m%omega_P(2,i_pt)*p%area_ratio_X*x%StC_x(2,i_pt)   &
                                      -m%a_G(3,i_pt)  )

            !side-to-side TLCD external forces of dependent degrees
            F_x_otlcd_WB_N = p%rho_Y*p%area_Y*((p%L_Y-p%B_Y)/2+x%StC_x(3,i_pt))*(                              &
                                         m%rddot_P(1,i_pt)                                                     &
                                      +2*m%omega_P(2,i_pt)*x%StC_x(4,i_pt)                                     &
                                        +m%alpha_P(2,i_pt)*((p%L_Y-p%B_Y)/2+x%StC_x(3,i_pt))                   &
                                        +m%alpha_P(3,i_pt)*p%B_Y/2-m%omega_P(2,i_pt)*m%omega_P(1,i_pt)*p%B_Y/2 &
                                        +m%omega_P(3,i_pt)*m%omega_P(1,i_pt)*((p%L_Y-p%B_Y)/2+x%StC_x(3,i_pt)) &
                                        -m%a_G(1,i_pt)   )
            F_y_otlcd_WB_N = p%rho_Y*p%area_Y*((p%L_Y-p%B_Y)/2+x%StC_x(3,i_pt))*(                              &
                                         m%rddot_P(2,i_pt)                                                     &
                                      -2*m%omega_P(1,i_pt)*x%StC_x(4,i_pt)                                     &
                                        -m%alpha_P(1,i_pt)*((p%L_Y-p%B_Y)/2+x%StC_x(3,i_pt))                   &
                                        +m%omega_P(3,i_pt)*m%omega_P(2,i_pt)*((p%L_Y-p%B_Y)/2+x%StC_x(3,i_pt)) &
                                        +m%omega_P(3,i_pt)*m%omega_P(3,i_pt)*p%B_Y/2                           &
                                        +m%omega_P(1,i_pt)*m%omega_P(1,i_pt)*p%B_Y/2                           &
                                        -m%a_G(2,i_pt)   )
            F_x_otlcd_WF_N = p%rho_Y*p%area_Y*((p%L_Y-p%B_Y)/2-x%StC_x(3,i_pt))*(                              &
                                         m%rddot_P(1,i_pt)                                                     &
                                      -2*m%omega_P(2,i_pt)*x%StC_x(4,i_pt)                                     &
                                        +m%alpha_P(2,i_pt)*((p%L_Y-p%B_Y)/2-x%StC_x(3,i_pt))                   &
                                        -m%alpha_P(2,i_pt)*p%B_Y/2                                             &
                                        +m%omega_P(2,i_pt)*m%omega_P(1,i_pt)*p%B_Y/2                           &
                                        +m%omega_P(3,i_pt)*m%omega_P(1,i_pt)*((p%L_Y-p%B_Y)/2-x%StC_x(3,i_pt)) &
                                        -m%a_G(1,i_pt)   )
            F_y_otlcd_WF_N = p%rho_Y*p%area_Y*((p%L_Y-p%B_Y)/2-x%StC_x(3,i_pt))*(                              &
                                         m%rddot_P(2,i_pt)                                                     &
                                      +2*m%omega_P(1,i_pt)*x%StC_x(4,i_pt)                                     &
                                        -m%alpha_P(1,i_pt)*((p%L_Y-p%B_Y)/2-x%StC_x(3,i_pt))                   &
                                        +m%omega_P(3,i_pt)*m%omega_P(2,i_pt)*((p%L_Y-p%B_Y)/2-x%StC_x(3,i_pt)) &
                                        -m%omega_P(3,i_pt)*m%omega_P(3,i_pt)*p%B_Y/2                           &
                                        -m%omega_P(1,i_pt)*m%omega_P(1,i_pt)*p%B_Y/2                           &
                                        -m%a_G(2,i_pt)   )
            F_x_otlcd_WH_N = p%rho_Y*p%area_Y/p%area_ratio_Y*p%B_Y*(                         &
                                          m%rddot_P(1,i_pt)                                  &
                                       -2*m%omega_P(3,i_pt)*p%area_ratio_Y*x%StC_x(4,i_pt)   &
                                         -m%a_G(1,i_pt)  )
            F_z_otlcd_WH_N = p%rho_Y*p%area_Y/p%area_ratio_Y*p%B_Y*(                         &
                                          m%rddot_P(3,i_pt)                                  &
                                       +2*m%omega_P(1,i_pt)*p%area_ratio_Y*x%StC_x(4,i_pt)   &
                                         -m%a_G(3,i_pt)  )

            ! forces and moments in local coordinates (from fore-aft and side-to-side TLCDs)
            m%F_P(1,i_pt) = -F_x_tlcd_WR_N - F_x_tlcd_WL_N - p%rho_X*(p%area_X/p%area_ratio_X)*p%B_X*dxdt%StC_x(2,i_pt)*p%area_ratio_X + F_x_otlcd_WB_N + F_x_otlcd_WF_N + F_x_otlcd_WH_N
            m%F_P(2,i_pt) = +F_y_tlcd_WR_N + F_y_tlcd_WL_N - p%rho_Y*(p%area_Y/p%area_ratio_Y)*p%B_Y*dxdt%StC_x(4,i_pt)*p%area_ratio_Y + F_y_tlcd_WH_N  - F_y_otlcd_WB_N - F_y_otlcd_WF_N
            m%F_P(3,i_pt) = -F_z_tlcd_WH_N - F_z_otlcd_WH_N

            m%M_P(1,i_pt) =  F_y_tlcd_WR_N*((p%L_X-p%B_X)/2+x%StC_x(1,i_pt)) + F_y_tlcd_WL_N*((p%L_X-p%B_X)/2-x%StC_x(1,i_pt)) - F_y_otlcd_WB_N*((p%L_Y-p%B_Y)/2+x%StC_x(3,i_pt)) - F_y_otlcd_WF_N*((p%L_Y-p%B_Y)/2-x%StC_x(3,i_pt))
            m%M_P(2,i_pt) = -F_x_tlcd_WR_N*((p%L_X-p%B_X)/2+x%StC_x(1,i_pt)) - F_x_tlcd_WL_N*((p%L_X-p%B_X)/2-x%StC_x(1,i_pt)) + F_x_otlcd_WB_N*((p%L_Y-p%B_Y)/2+x%StC_x(3,i_pt)) + F_x_otlcd_WF_N*((p%L_Y-p%B_Y)/2-x%StC_x(3,i_pt))
            m%M_P(3,i_pt) =  F_y_tlcd_WR_N*p%B_X*.5 - F_y_tlcd_WL_N*p%B_X*.5 + F_x_otlcd_WB_N*p%B_Y*.5 - F_x_otlcd_WF_N*p%B_Y*.5

            ! forces and moments in global coordinates
            y%Mesh(i_pt)%Force(:,1)  = real(matmul(transpose(u%Mesh(i_pt)%Orientation(:,:,1)), m%F_P(1:3,i_pt)),ReKi)
            y%Mesh(i_pt)%Moment(:,1) = real(matmul(transpose(u%Mesh(i_pt)%Orientation(:,:,1)), m%M_P(1:3,i_pt)),ReKi)
         enddo
      ELSEIF ( p%StC_DOF_MODE == DOFMode_Prescribed ) THEN
         !  Note that the prescribed force is applied the same to all Mesh pts
         !  that are passed into this instance of the StC
         do i=1,3
            ! Get interpolated force   -- this is not in any particular coordinate system yet
            m%F_P(i,:)    = InterpStp( real(Time,ReKi), p%StC_PrescribedForce(1,:),p%StC_PrescribedForce(i+1,:),m%PrescribedInterpIdx, size(p%StC_PrescribedForce,2))
            ! Get interpolated moment  -- this is not in any particular coordinate system yet
            m%M_P(i,:)    = InterpStp( real(Time,ReKi), p%StC_PrescribedForce(1,:),p%StC_PrescribedForce(i+4,:),m%PrescribedInterpIdx, size(p%StC_PrescribedForce,2))
         enddo
         if (p%PrescribedForcesCoordSys == PRESCRIBED_FORCE_GLOBAL) then
            ! Global coords
            do i_pt=1,p%NumMeshPts
               y%Mesh(i_pt)%Force(1:3,1)  =  m%F_P(1:3,i_pt)
               y%Mesh(i_pt)%Moment(1:3,1) =  m%M_P(1:3,i_pt)
            enddo
         elseif (p%PrescribedForcesCoordSys == PRESCRIBED_FORCE_LOCAL) then
            ! local coords
            do i_pt=1,p%NumMeshPts
               y%Mesh(i_pt)%Force(1:3,1)  =  matmul(transpose(u%Mesh(i_pt)%Orientation(:,:,1)), m%F_P(1:3,i_pt))
               y%Mesh(i_pt)%Moment(1:3,1) =  matmul(transpose(u%Mesh(i_pt)%Orientation(:,:,1)), m%M_P(1:3,i_pt))
            enddo
         endif
      ELSEIF ( p%StC_DOF_MODE == DOFMode_ForceDLL ) THEN
         !  Note that the prescribed force is applied the same to all Mesh pts
         !  that are passed into this instance of the StC
         if (p%PrescribedForcesCoordSys == PRESCRIBED_FORCE_GLOBAL) then
            ! Global coords
            do i_pt=1,p%NumMeshPts
               y%Mesh(i_pt)%Force(1:3,1)  =  m%F_ext(1:3,i_pt)
               y%Mesh(i_pt)%Moment(1:3,1) =  0
            enddo
         ! Leave in for now just in case we decide there is a use case for a follower force from the DLL
         ! elseif (p%PrescribedForcesCoordSys == PRESCRIBED_FORCE_LOCAL) then
         !    ! local coords
         !    do i_pt=1,p%NumMeshPts
         !       y%Mesh(i_pt)%Force(1:3,1)  =  matmul(transpose(u%Mesh(i_pt)%Orientation(:,:,1)), m%F_P(1:3,i_pt))
         !       y%Mesh(i_pt)%Moment(1:3,1) =  matmul(transpose(u%Mesh(i_pt)%Orientation(:,:,1)), m%M_P(1:3,i_pt))
         !    enddo
         endif
      END IF

      ! Set output values for the measured displacements for  
      do i=1,size(p%StC_CChan)
         if (p%StC_CChan(i) > 0) then
            if (p%StC_DOF_MODE == DOFMode_Indept .or. p%StC_DOF_MODE == DOFMode_Omni) then
               !  The states are sized by (6,NumMeshPts).  NumMeshPts is then used to set
               !  size of StC_CChan as well.  For safety, we will check it here.
               y%MeasDisp(1:3,p%StC_CChan(i))  = (/ x%StC_x(1,i), x%StC_x(3,i), x%StC_x(5,i) /)
               y%MeasVel( 1:3,p%StC_CChan(i))  = (/ x%StC_x(2,i), x%StC_x(4,i), x%StC_x(6,i) /)
            else
               y%MeasDisp(1:3,p%StC_CChan(i))  = 0.0_ReKi
               y%MeasVel( 1:3,p%StC_CChan(i))  = 0.0_ReKi
            endif
         endif
      enddo

      call CleanUp()

CONTAINS
   subroutine CleanUp()
      call StC_DestroyContState(dxdt,ErrStat2,ErrMsg2)    !Ignore error status
   end subroutine CleanUp
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'StC_CalcOutput')
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed
END SUBROUTINE StC_CalcOutput

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
SUBROUTINE StC_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
!..................................................................................................................................

      REAL(DbKi),                    INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(StC_InputType),           INTENT(IN   )  :: u           !< Inputs at Time
      TYPE(StC_ParameterType),       INTENT(IN   )  :: p           !< Parameters
      TYPE(StC_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(StC_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(StC_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(StC_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time
      TYPE(StC_ContinuousStateType), INTENT(  OUT)  :: dxdt        !< Continuous state derivatives at Time
      TYPE(StC_MiscVarType),         INTENT(INOUT)  :: m           !< Misc (optimization) variables
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

      REAL(ReKi), dimension(3,p%NumMeshPts)           :: K          ! tuned mass damper stiffness
      Real(ReKi)                                      :: denom      ! denominator for omni-direction factors
      integer(IntKi)                                  :: i_pt       ! Generic counter for mesh point

      ! Local error handling
      integer(IntKi)       :: ErrStat2
      character(ErrMsgLen) :: ErrMsg2

         ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""


      call AllocAry(dxdt%StC_x,6, p%NumMeshPts,'dxdt%StC_x',  ErrStat2,ErrMsg2); if (Failed()) return;

         ! compute stop force (m%F_stop)
      IF (p%Use_F_TBL) THEN
         m%F_stop = 0.0_ReKi
      ELSE
         CALL StC_CalcStopForce(x,p,m%F_stop)
      END IF

      ! Compute stiffness -- Note that this value may be overwritten by controller
      IF (p%Use_F_TBL) THEN ! use stiffness table
         CALL SpringForceExtrapInterp(x,p,m%F_table,ErrStat2,ErrMsg2);  if (Failed()) return;
         K = 0.0_ReKi
      ELSE ! use preset values
         K(1,:) = p%K_X
         K(2,:) = p%K_Y
         K(3,:) = p%K_Z
      END IF


      ! Compute accelerations and velocities in local coordinates
      do i_pt=1,p%NumMeshPts
         m%a_G(:,i_pt)     = matmul(u%Mesh(i_pt)%Orientation(:,:,1),p%Gravity)
         m%rdisp_P(:,i_pt) = matmul(u%Mesh(i_pt)%Orientation(:,:,1),u%Mesh(i_pt)%TranslationDisp(:,1))   ! for ground StC_GroundHookDamp
         m%rdot_P(:,i_pt)  = matmul(u%Mesh(i_pt)%Orientation(:,:,1),u%Mesh(i_pt)%TranslationVel(:,1))    ! for ground StC_GroundHookDamp
         m%rddot_P(:,i_pt) = matmul(u%Mesh(i_pt)%Orientation(:,:,1),u%Mesh(i_pt)%TranslationAcc(:,1))
         m%omega_P(:,i_pt) = matmul(u%Mesh(i_pt)%Orientation(:,:,1),u%Mesh(i_pt)%RotationVel(:,1))
         m%alpha_P(:,i_pt) = matmul(u%Mesh(i_pt)%Orientation(:,:,1),u%Mesh(i_pt)%RotationAcc(:,1))
      enddo

      ! NOTE: m%F_stop and m%F_table are calculated earlier
      IF ((p%StC_DOF_MODE == ControlMode_None) .or. (p%StC_DOF_MODE == DOFMode_Prescribed)) THEN
         do i_pt=1,p%NumMeshPts
            ! Aggregate acceleration terms
            m%Acc(1:3,i_pt) = 0.0_ReKi 
         enddo

      ELSEIF (p%StC_DOF_MODE == DOFMode_Indept) THEN

         do i_pt=1,p%NumMeshPts
            ! Aggregate acceleration terms
            m%Acc(1,i_pt) = - m%rddot_P(1,i_pt) + m%a_G(1,i_pt) + 1 / p%M_X * ( m%F_ext(1,i_pt) + m%F_stop(1,i_pt) - m%F_table(1,i_pt) )
            m%Acc(2,i_pt) = - m%rddot_P(2,i_pt) + m%a_G(2,i_pt) + 1 / p%M_Y * ( m%F_ext(2,i_pt) + m%F_stop(2,i_pt) - m%F_table(2,i_pt) )
            m%Acc(3,i_pt) = - m%rddot_P(3,i_pt) + m%a_G(3,i_pt) + 1 / p%M_Z * ( m%F_ext(3,i_pt) + m%F_stop(3,i_pt) - m%F_table(3,i_pt) + p%StC_Z_PreLd )
         enddo

      ELSE IF (p%StC_DOF_MODE == DOFMode_Omni) THEN

         do i_pt=1,p%NumMeshPts
            denom = SQRT(x%StC_x(1,i_pt)**2+x%StC_x(3,i_pt)**2)
            IF ( EqualRealNos( denom, 0.0_ReKi) ) THEN
                m%F_k(1,i_pt) = 0.0
                m%F_k(2,i_pt) = 0.0
            ELSE
                  m%F_k(1,i_pt) = x%StC_x(1,i_pt)/denom
                  m%F_k(2,i_pt) = x%StC_x(3,i_pt)/denom
            END IF
            m%F_k(3,i_pt) = 0.0

            ! Aggregate acceleration terms
            m%Acc(1,i_pt) = - m%rddot_P(1,i_pt) + m%a_G(1,i_pt) + 1 / p%M_XY * ( m%F_ext(1,i_pt) + m%F_stop(1,i_pt) - m%F_table(1,i_pt)*(m%F_k(1,i_pt)) )
            m%Acc(2,i_pt) = - m%rddot_P(2,i_pt) + m%a_G(2,i_pt) + 1 / p%M_XY * ( m%F_ext(2,i_pt) + m%F_stop(2,i_pt) - m%F_table(2,i_pt)*(m%F_k(2,i_pt)) )
            m%Acc(3,i_pt) = 0.0_ReKi
         enddo

      ENDIF


      ! Compute the first time derivatives, dxdt%StC_x(1) and dxdt%StC_x(3), of the continuous states,:
      ! Compute elements 1 and 3 of dxdt%StC_x so that we can compute m%C_ctrl,m%C_Brake, and m%F_fr in StC_GroundHookDamp if necessary
      IF ((p%StC_DOF_MODE == ControlMode_None) .or. (p%StC_DOF_MODE == DOFMode_Prescribed)) THEN

         dxdt%StC_x = 0.0_ReKi ! Whole array

      ELSE

         IF (p%StC_DOF_MODE == DOFMode_Indept .AND. .NOT. p%StC_X_DOF) THEN
            do i_pt=1,p%NumMeshPts
               dxdt%StC_x(1,i_pt) = 0.0_ReKi
            enddo
         ELSE
            do i_pt=1,p%NumMeshPts
               dxdt%StC_x(1,i_pt) = x%StC_x(2,i_pt)
            enddo
         END IF

         IF (p%StC_DOF_MODE == DOFMode_Indept .AND. .NOT. p%StC_Y_DOF) THEN
            do i_pt=1,p%NumMeshPts
               dxdt%StC_x(3,i_pt) = 0.0_ReKi
            enddo
         ELSE
            do i_pt=1,p%NumMeshPts
               dxdt%StC_x(3,i_pt) = x%StC_x(4,i_pt)
            enddo
         END IF

         IF (p%StC_DOF_MODE == DOFMode_Indept .AND. .NOT. p%StC_Z_DOF) THEN
            do i_pt=1,p%NumMeshPts
               dxdt%StC_x(5,i_pt) = 0.0_ReKi
            enddo
         ELSE
            do i_pt=1,p%NumMeshPts
               dxdt%StC_x(5,i_pt) = x%StC_x(6,i_pt)
            enddo
         END IF

         if ( .not. (p%StC_DOF_MODE == DOFMode_Indept .AND. p%StC_Z_DOF)) then      ! z not used in any other configuration
            do i_pt=1,p%NumMeshPts
               dxdt%StC_x(5,i_pt) = 0.0_ReKi
            enddo
         endif

      ENDIF


      ! compute damping for dxdt%StC_x(2), dxdt%StC_x(4), and dxdt%StC_x(6)
      IF (p%StC_CMODE == ControlMode_None) THEN
         m%C_ctrl(1,:) = p%C_X
         m%C_ctrl(2,:) = p%C_Y
         m%C_ctrl(3,:) = p%C_Z
         m%C_Brake = 0.0_ReKi
         m%F_fr    = 0.0_ReKi
      ELSE IF (p%StC_CMODE == CMODE_Semi) THEN ! ground hook control
         CALL StC_GroundHookDamp(dxdt,x,u,p,m%rdisp_P,m%rdot_P,m%C_ctrl,m%C_Brake,m%F_fr)
      ELSE IF (p%StC_CMODE == CMODE_ActiveDLL) THEN   ! Active control from DLL
         call StC_ActiveCtrl_StiffDamp(u,p,m%K,m%C_ctrl,m%C_Brake,m%F_ext)
         m%F_fr    = 0.0_ReKi
         if (.not. p%Use_F_TBL) then
            K(1:3,:) = m%K(1:3,:)
!FIXME: for the derivative, I don't know how this should be handled.
!!!   Defaulting to how the stiffness table operates for now, but leaving
!!!   this next code chunk in case it shuold be handed this way instead
         !else
         !   ! Make commanded stiffness a perturbation about the table value (to avoid double counting the table)
         !   ! NOTE: This has not been verified and may have unintended consequences.
         !   do i_pt=1,p%NumMeshPts
         !      K(1:3,i_pt) = m%F_table(1:3,i_pt) - m%K(1:3,i_pt)
         !   enddo
         endif
      END IF


      ! Compute the first time derivatives, dxdt%StC_x(2), dxdt%StC_x(4), and dxdt%StC_x(6), of the continuous states,:
      IF (p%StC_DOF_MODE == DOFMode_Indept) THEN

         IF (p%StC_X_DOF) THEN
            do i_pt=1,p%NumMeshPts
               dxdt%StC_x(2,i_pt) =  ( m%omega_P(2,i_pt)**2 + m%omega_P(3,i_pt)**2 - K(1,i_pt) / p%M_X) * x%StC_x(1,i_pt) &
                                   - ( m%C_ctrl( 1,i_pt)/p%M_X ) * x%StC_x(2,i_pt)                                   &
                                   - ( m%C_Brake(1,i_pt)/p%M_X ) * x%StC_x(2,i_pt)                                   &
                                   + m%Acc(1,i_pt) + m%F_fr(1,i_pt) / p%M_X
            enddo
         ELSE
            do i_pt=1,p%NumMeshPts
               dxdt%StC_x(2,i_pt) = 0.0_ReKi
            enddo
         END IF
         IF (p%StC_Y_DOF) THEN
            do i_pt=1,p%NumMeshPts
               dxdt%StC_x(4,i_pt) =  ( m%omega_P(1,i_pt)**2 + m%omega_P(3,i_pt)**2 - K(2,i_pt) / p%M_Y) * x%StC_x(3,i_pt) &
                                   - ( m%C_ctrl( 2,i_pt)/p%M_Y ) * x%StC_x(4,i_pt)                                   &
                                   - ( m%C_Brake(2,i_pt)/p%M_Y ) * x%StC_x(4,i_pt)                                   &
                                   + m%Acc(2,i_pt) + m%F_fr(2,i_pt) / p%M_Y
            enddo
         ELSE
            do i_pt=1,p%NumMeshPts
               dxdt%StC_x(4,i_pt) = 0.0_ReKi
            enddo
         END IF
         IF (p%StC_Z_DOF) THEN
            do i_pt=1,p%NumMeshPts
               dxdt%StC_x(6,i_pt) =  ( m%omega_P(1,i_pt)**2 + m%omega_P(2,i_pt)**2 - K(3,i_pt) / p%M_Z) * x%StC_x(5,i_pt) &
                                   - ( m%C_ctrl( 3,i_pt)/p%M_Z ) * x%StC_x(6,i_pt)                                   &
                                   - ( m%C_Brake(3,i_pt)/p%M_Z ) * x%StC_x(6,i_pt)                                   &
                                   + m%Acc(3,i_pt) + m%F_fr(3,i_pt) / p%M_Z
            enddo
         ELSE
            do i_pt=1,p%NumMeshPts
               dxdt%StC_x(6,i_pt) = 0.0_ReKi
            enddo
         END IF

      ELSE IF (p%StC_DOF_MODE == DOFMode_Omni) THEN   ! Only includes X and Y
               ! Compute the first time derivatives of the continuous states of Omnidirectional tuned masse damper mode by sm 2015-0904
         do i_pt=1,p%NumMeshPts
            dxdt%StC_x(2,i_pt) =  ( m%omega_P(2,i_pt)**2 + m%omega_P(3,i_pt)**2 - K(1,i_pt) / p%M_XY) * x%StC_x(1,i_pt)   &
                                - ( m%C_ctrl( 1,i_pt)/p%M_XY ) * x%StC_x(2,i_pt)                                     &
                                - ( m%C_Brake(1,i_pt)/p%M_XY ) * x%StC_x(2,i_pt)                                     &
                                +  m%Acc(1,i_pt) + 1/p%M_XY * ( m%F_fr(1,i_pt) )                                     &
                                - ( m%omega_P(1,i_pt)*m%omega_P(2,i_pt) - m%alpha_P(3,i_pt) ) * x%StC_x(3,i_pt)      &
                               +2 * m%omega_P(3,i_pt) * x%StC_x(4,i_pt)
            dxdt%StC_x(4,i_pt) =  ( m%omega_P(1,i_pt)**2 + m%omega_P(3,i_pt)**2 - K(2,i_pt) / p%M_XY) * x%StC_x(3,i_pt)   &
                                - ( m%C_ctrl( 2,i_pt)/p%M_XY ) * x%StC_x(4,i_pt)                                     &
                                - ( m%C_Brake(2,i_pt)/p%M_XY ) * x%StC_x(4,i_pt)                                     &
                                +  m%Acc(2,i_pt) + 1/p%M_XY * ( m%F_fr(2,i_pt) )                                     &
                                - ( m%omega_P(1,i_pt)*m%omega_P(2,i_pt) + m%alpha_P(3,i_pt) ) * x%StC_x(1,i_pt)      &
                               -2 * m%omega_P(3,i_pt) * x%StC_x(2,i_pt)
            dxdt%StC_x(6,i_pt) = 0.0_ReKi ! Z is off
         enddo

      ELSE IF (p%StC_DOF_MODE == DOFMode_TLCD) THEN !MEG & SP
         ! Compute the first time derivatives of the continuous states of TLCD mode
         do i_pt=1,p%NumMeshPts
            dxdt%StC_x(2,i_pt) = (2*p%rho_X*p%area_X*x%StC_x(1,i_pt)*m%rddot_P(3,i_pt)                                  &
                                   +p%rho_X*p%area_X*p%B_X*m%alpha_P(2,i_pt)*((p%L_X-p%B_X)/2)                          &
                                   -p%rho_X*p%area_X*p%B_X*m%omega_P(1,i_pt)*m%omega_P(3,i_pt)*((p%L_X-p%B_X)/2)        &
                                 +2*p%rho_X*p%area_X*m%omega_P(1,i_pt)*m%omega_P(1,i_pt)*x%StC_x(1,i_pt)*(p%L_X-p%B_X)  &
                                 +2*p%rho_X*p%area_X*m%omega_P(2,i_pt)*m%omega_P(2,i_pt)*x%StC_x(1,i_pt)*(p%L_X-p%B_X)  &
                                 +2*p%rho_X*p%area_X*x%StC_x(1,i_pt)*m%a_G(3,i_pt)                                      &
                                   -p%rho_X*p%area_X*p%B_X*m%rddot_P(1,i_pt)                                            &
                                   +p%rho_X*p%area_X*p%B_X*m%a_G(1,i_pt)                                                &
                                -.5*p%rho_X*p%area_X*p%headLossCoeff_X*p%area_ratio_X*p%area_ratio_X*x%StC_x(2,i_pt)    &
                                       *ABS(x%StC_x(2,i_pt)))/(p%rho_X*p%area_X*(p%L_X-p%B_X+p%area_ratio_X*p%B_X))        
            dxdt%StC_x(4,i_pt) = (2*p%rho_Y*p%area_Y*x%StC_x(3,i_pt)*m%rddot_P(3,i_pt)                                     &
                                   +p%rho_Y*p%area_Y*p%B_Y*m%alpha_P(1,i_pt)*((p%L_Y-p%B_Y)/2)                             &
                                   -p%rho_Y*p%area_Y*p%B_Y*m%omega_P(2,i_pt)*m%omega_P(3,i_pt)*((p%L_Y-p%B_Y)/2)           &
                                 +2*p%rho_Y*p%area_Y*x%StC_x(3,i_pt)*m%omega_P(1,i_pt)*m%omega_P(1,i_pt)*(p%L_Y-p%B_Y)     &
                                 +2*p%rho_Y*p%area_Y*x%StC_x(3,i_pt)*m%omega_P(2,i_pt)*m%omega_P(2,i_pt)*(p%L_Y-p%B_Y)     &
                                 +2*p%rho_Y*p%area_Y*x%StC_x(3,i_pt)*m%a_G(3,i_pt)-p%rho_Y*p%area_Y*p%B_Y*m%rddot_P(2,i_pt)&
                                   +p%rho_Y*p%area_Y*p%B_Y*m%a_G(2,i_pt)                                                   &
                                -.5*p%rho_Y*p%area_Y*p%headLossCoeff_Y*p%area_ratio_Y*p%area_ratio_Y*x%StC_x(4,i_pt)       &
                                       *ABS(x%StC_x(4,i_pt)))/(p%rho_Y*p%area_Y*(p%L_Y-p%B_Y+p%area_ratio_Y*p%B_Y))
            dxdt%StC_x(6,i_pt) = 0.0_ReKi ! Z is off
         enddo

      ELSE IF ( p%StC_DOF_MODE == DOFMode_Prescribed .or. p%StC_DOF_MODE == DOFMode_ForceDLL) THEN
      ! if prescribed forces, there are no states to advance, so return
         do i_pt=1,p%NumMeshPts
            dxdt%StC_x(1,i_pt) = 0
            dxdt%StC_x(2,i_pt) = 0
            dxdt%StC_x(3,i_pt) = 0
            dxdt%StC_x(4,i_pt) = 0
            dxdt%StC_x(5,i_pt) = 0
            dxdt%StC_x(6,i_pt) = 0
         enddo
         return
      END IF

      call CleanUp()
      return

CONTAINS
   subroutine CleanUp()
   end subroutine CleanUp
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'StC_CalcContStateDeriv')
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed
END SUBROUTINE StC_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE StC_CalcStopForce(x,p,F_stop)
   TYPE(StC_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(StC_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   Real(ReKi), dimension(:,:),    INTENT(INOUT)  :: F_stop      !< stop forces
   ! local variables
   Real(ReKi), dimension(3)                      :: F_SK      !stop spring forces
   Real(ReKi), dimension(3)                      :: F_SD      !stop damping forces
   INTEGER(IntKi)                                :: i         ! counter
   INTEGER(IntKi)                                :: i_pt      ! counter for mesh points
   INTEGER(IntKi)                                :: j         ! counter for index into x%StC_x
   do i_pt=1,p%NumMeshPts
      DO i=1,3 ! X, Y, and Z
         j=2*(i-1)+1
         IF ( x%StC_x(j,i_pt) > p%P_SP(i) ) THEN
            F_SK(i) = p%K_S(i) *( p%P_SP(i) - x%StC_x(j,i_pt)  )
         ELSEIF ( x%StC_x(j,i_pt) < p%N_SP(i) ) THEN
            F_SK(i) = p%K_S(i) * ( p%N_SP(i) - x%StC_x(j,i_pt) )
         ELSE
            F_SK(i)  = 0.0_ReKi
         ENDIF
         IF ( (x%StC_x(j,i_pt) > p%P_SP(i)) .AND. (x%StC_x(j+1,i_pt) > 0) ) THEN
            F_SD(i) = -p%C_S(i) *( x%StC_x(j+1,i_pt)  )
         ELSEIF ( (x%StC_x(j,i_pt) < p%N_SP(i)) .AND. (x%StC_x(j+1,i_pt) < 0) ) THEN
            F_SD(i) = -p%C_S(i) *( x%StC_x(j+1,i_pt)  )
         ELSE
            F_SD(i)  = 0.0_ReKi
         ENDIF
         F_stop(i,i_pt) = F_SK(i) + F_SD(i)
      END DO
   enddo
END SUBROUTINE StC_CalcStopForce
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE StC_GroundHookDamp(dxdt,x,u,p,rdisp_P,rdot_P,C_ctrl,C_Brake,F_fr)
   TYPE(StC_ContinuousStateType),         INTENT(IN   )     :: dxdt        !< Derivative of continuous states at Time (needs elements 1 and 3 only)
   TYPE(StC_ContinuousStateType),         INTENT(IN   )     :: x           !< Continuous states at Time
   TYPE(StC_InputType),                   INTENT(IN   )     :: u           !< Inputs at Time
   TYPE(StC_ParameterType),               INTENT(IN)        :: p           !< The module's parameter data
   REAL(ReKi), dimension(:,:),            INTENT(IN   )     :: rdisp_P     !< translational displacement in local coordinates
   REAL(ReKi), dimension(:,:),            INTENT(IN   )     :: rdot_P      !< translational velocity     in local coordinates
   REAL(ReKi), dimension(:,:),            INTENT(INOUT)     :: C_ctrl      !< extrapolated/interpolated stiffness values
   REAL(ReKi), dimension(:,:),            INTENT(INOUT)     :: C_Brake     !< extrapolated/interpolated braking   values
   REAL(ReKi), dimension(:,:),            INTENT(INOUT)     :: F_fr        !< Friction forces
   INTEGER(IntKi)                                           :: i_pt        !< generic counter for mesh points


   do i_pt=1,p%NumMeshPts
      IF (p%StC_CMODE == CMODE_Semi .AND. p%StC_SA_MODE == SA_CMODE_GH_vel) THEN ! velocity-based ground hook control with high damping for braking

         !X
         IF (dxdt%StC_x(1,i_pt) * rdot_P(1,i_pt) <= 0 ) THEN
            C_ctrl(1,i_pt) = p%StC_X_C_HIGH
         ELSE
            C_ctrl(1,i_pt) = p%StC_X_C_LOW
         END IF

         !Brake X
         IF      ( (x%StC_x(1,i_pt) > p%P_SP(1)-0.2) .AND. (x%StC_x(2,i_pt) > 0) ) THEN
            C_Brake(1,i_pt) = p%StC_X_C_BRAKE
         ELSE IF ( (x%StC_x(1,i_pt) < p%N_SP(1)+0.2) .AND. (x%StC_x(2,i_pt) < 0) ) THEN
            C_Brake(1,i_pt) = p%StC_X_C_BRAKE
         ELSE
            C_Brake(1,i_pt) = 0
         END IF


         ! Y
         IF (dxdt%StC_x(3,i_pt) * rdot_P(2,i_pt) <= 0 ) THEN
            C_ctrl(2,i_pt) = p%StC_Y_C_HIGH
         ELSE
            C_ctrl(2,i_pt) = p%StC_Y_C_LOW
         END IF

         !Brake Y
         IF      ( (x%StC_x(3,i_pt) > p%P_SP(2)-0.2) .AND. (x%StC_x(4,i_pt) > 0) ) THEN
            C_Brake(2,i_pt) = p%StC_Y_C_BRAKE
         ELSE IF ( (x%StC_x(3,i_pt) < p%N_SP(2)+0.2) .AND. (x%StC_x(4,i_pt) < 0) ) THEN
            C_Brake(2,i_pt) = p%StC_Y_C_BRAKE
         ELSE
            C_Brake(2,i_pt) = 0
         END IF


         ! Z
         IF (dxdt%StC_x(5,i_pt) * rdot_P(3,i_pt) <= 0 ) THEN
            C_ctrl(3,i_pt) = p%StC_Z_C_HIGH
         ELSE
            C_ctrl(3,i_pt) = p%StC_Z_C_LOW
         END IF

         !Brake Z
         IF      ( (x%StC_x(5,i_pt) > p%P_SP(3)-0.2) .AND. (x%StC_x(6,i_pt) > 0) ) THEN
            C_Brake(3,i_pt) = p%StC_Z_C_BRAKE
         ELSE IF ( (x%StC_x(5,i_pt) < p%N_SP(3)+0.2) .AND. (x%StC_x(6,i_pt) < 0) ) THEN
            C_Brake(3,i_pt) = p%StC_Z_C_BRAKE
         ELSE
            C_Brake(3,i_pt) = 0
         END IF

      ELSE IF (p%StC_CMODE == CMODE_Semi .AND. p%StC_SA_MODE == SA_CMODE_GH_invVel) THEN ! Inverse velocity-based ground hook control with high damping for braking

         ! X
         IF (dxdt%StC_x(1,i_pt) * rdot_P(1,i_pt) >= 0 ) THEN
            C_ctrl(1,i_pt) = p%StC_X_C_HIGH
         ELSE
            C_ctrl(1,i_pt) = p%StC_X_C_LOW
         END IF

         !Brake X
         IF      ( (x%StC_x(1,i_pt) > p%P_SP(1)-0.2) .AND. (x%StC_x(2,i_pt) > 0) ) THEN
            C_Brake(1,i_pt) = p%StC_X_C_BRAKE
         ELSE IF ( (x%StC_x(1,i_pt) < p%N_SP(1)+0.2) .AND. (x%StC_x(2,i_pt) < 0) ) THEN
            C_Brake(1,i_pt) = p%StC_X_C_BRAKE
         ELSE
            C_Brake(1,i_pt) = 0
         END IF

         ! Y
         IF (dxdt%StC_x(3,i_pt) * rdot_P(2,i_pt) >= 0 ) THEN
            C_ctrl(2,i_pt) = p%StC_Y_C_HIGH
         ELSE
            C_ctrl(2,i_pt) = p%StC_Y_C_LOW
         END IF

         !Brake Y
         IF      ( (x%StC_x(3,i_pt) > p%P_SP(2)-0.2) .AND. (x%StC_x(4,i_pt) > 0) ) THEN
            C_Brake(2,i_pt) = p%StC_Y_C_BRAKE
         ELSE IF ( (x%StC_x(3,i_pt) < p%N_SP(2)+0.2) .AND. (x%StC_x(4,i_pt) < 0) ) THEN
            C_Brake(2,i_pt) = p%StC_Y_C_BRAKE
         ELSE
            C_Brake(2,i_pt) = 0
         END IF

         ! Z
         IF (dxdt%StC_x(5,i_pt) * rdot_P(3,i_pt) >= 0 ) THEN
            C_ctrl(3,i_pt) = p%StC_Z_C_HIGH
         ELSE
            C_ctrl(3,i_pt) = p%StC_Z_C_LOW
         END IF

         !Brake Z
         IF      ( (x%StC_x(5,i_pt) > p%P_SP(3)-0.2) .AND. (x%StC_x(6,i_pt) > 0) ) THEN
            C_Brake(3,i_pt) = p%StC_Z_C_BRAKE
         ELSE IF ( (x%StC_x(5,i_pt) < p%N_SP(3)+0.2) .AND. (x%StC_x(6,i_pt) < 0) ) THEN
            C_Brake(3,i_pt) = p%StC_Z_C_BRAKE
         ELSE
            C_Brake(3,i_pt) = 0
         END IF

      ELSE IF (p%StC_CMODE == CMODE_Semi .AND. p%StC_SA_MODE == SA_CMODE_GH_disp) THEN ! displacement-based ground hook control with high damping for braking

         ! X
         IF (dxdt%StC_x(1,i_pt) * rdisp_P(1,i_pt) <= 0 ) THEN
            C_ctrl(1,i_pt) = p%StC_X_C_HIGH
         ELSE
            C_ctrl(1,i_pt) = p%StC_X_C_LOW
         END IF

         !Brake X
         IF      ( (x%StC_x(1,i_pt) > p%P_SP(1)-0.2) .AND. (x%StC_x(2,i_pt) > 0) ) THEN
            C_Brake(1,i_pt) = p%StC_X_C_BRAKE
         ELSE IF ( (x%StC_x(1,i_pt) < p%N_SP(1)+0.2) .AND. (x%StC_x(2,i_pt) < 0) ) THEN
            C_Brake(1,i_pt) = p%StC_X_C_BRAKE
         ELSE
            C_Brake(1,i_pt) = 0
         END IF

         ! Y
         IF (dxdt%StC_x(3,i_pt) * rdisp_P(2,i_pt) <= 0 ) THEN
            C_ctrl(2,i_pt) = p%StC_Y_C_HIGH
         ELSE
            C_ctrl(2,i_pt) = p%StC_Y_C_LOW
         END IF

         !Brake Y
         IF      ( (x%StC_x(3,i_pt) > p%P_SP(2)-0.2) .AND. (x%StC_x(4,i_pt) > 0) ) THEN
            C_Brake(2,i_pt) = p%StC_Y_C_BRAKE
         ELSE IF ( (x%StC_x(3,i_pt) < p%N_SP(2)+0.2) .AND. (x%StC_x(4,i_pt) < 0) ) THEN
            C_Brake(2,i_pt) = p%StC_Y_C_BRAKE
         ELSE
            C_Brake(2,i_pt) = 0
         END IF

         ! Z
         IF (dxdt%StC_x(5,i_pt) * rdisp_P(3,i_pt) <= 0 ) THEN
            C_ctrl(3,i_pt) = p%StC_Z_C_HIGH
         ELSE
            C_ctrl(3,i_pt) = p%StC_Z_C_LOW
         END IF

         !Brake Z
         IF      ( (x%StC_x(5,i_pt) > p%P_SP(3)-0.2) .AND. (x%StC_x(6,i_pt) > 0) ) THEN
            C_Brake(3,i_pt) = p%StC_Z_C_BRAKE
         ELSE IF ( (x%StC_x(3,i_pt) < p%N_SP(3)+0.2) .AND. (x%StC_x(6,i_pt) < 0) ) THEN
            C_Brake(3,i_pt) = p%StC_Z_C_BRAKE
         ELSE
            C_Brake(3,i_pt) = 0
         END IF

      ELSE IF (p%StC_CMODE == CMODE_Semi .AND. p%StC_SA_MODE == SA_CMODE_Ph_FF) THEN ! Phase Difference Algorithm with Friction Force
            ! X
            ! (a)
         IF      (rdisp_P(1,i_pt) > 0 .AND. rdot_P(1,i_pt) < 0 .AND. x%StC_x(1,i_pt) > 0 .AND. dxdt%StC_x(1,i_pt) < 0) THEN
            F_fr(1,i_pt) = p%StC_X_C_HIGH
            ! (b)
         ELSE IF (rdisp_P(1,i_pt) < 0 .AND. rdot_P(1,i_pt) > 0 .AND. x%StC_x(1,i_pt) < 0 .AND. dxdt%StC_x(1,i_pt) > 0) THEN
            F_fr(1,i_pt) = -p%StC_X_C_HIGH
            ! (c)
         ELSE IF (rdisp_P(1,i_pt) < 0 .AND. rdot_P(1,i_pt) < 0 .AND. x%StC_x(1,i_pt) > 0 .AND. dxdt%StC_x(1,i_pt) > 0) THEN
            F_fr(1,i_pt) = -p%StC_X_C_HIGH
         ELSE IF (rdisp_P(1,i_pt) > 0 .AND. rdot_P(1,i_pt) > 0 .AND. x%StC_x(1,i_pt) < 0 .AND. dxdt%StC_x(1,i_pt) < 0) THEN
            F_fr(1,i_pt) = p%StC_X_C_HIGH
         ELSE
            F_fr(1,i_pt) = p%StC_X_C_LOW
         END IF

         !Brake X
         IF ( (x%StC_x(1,i_pt) > p%P_SP(1)-0.2) .AND. (x%StC_x(2,i_pt) > 0) ) THEN
            C_Brake(1,i_pt) = p%StC_X_C_BRAKE
         ELSE IF ( (x%StC_x(1,i_pt) < p%N_SP(1)+0.2) .AND. (x%StC_x(2,i_pt) < 0) ) THEN
            C_Brake(1,i_pt) = p%StC_X_C_BRAKE
         ELSE
            C_Brake(1,i_pt) = 0
         END IF

            ! Y
            ! (a)
         IF      (rdisp_P(2,i_pt) > 0 .AND. rdot_P(2,i_pt) < 0 .AND. x%StC_x(3,i_pt) > 0 .AND. dxdt%StC_x(3,i_pt) < 0) THEN
            F_fr(2,i_pt) = p%StC_Y_C_HIGH
            ! (b)
         ELSE IF (rdisp_P(2,i_pt) < 0 .AND. rdot_P(2,i_pt) > 0 .AND. x%StC_x(3,i_pt) < 0 .AND. dxdt%StC_x(3,i_pt) > 0) THEN
            F_fr(2,i_pt) = -p%StC_Y_C_HIGH
            ! (c)
         ELSE IF (rdisp_P(2,i_pt) < 0 .AND. rdot_P(2,i_pt) < 0 .AND. x%StC_x(3,i_pt) > 0 .AND. dxdt%StC_x(3,i_pt) > 0) THEN
            F_fr(2,i_pt) = -p%StC_Y_C_HIGH
         ELSE IF (rdisp_P(2,i_pt) > 0 .AND. rdot_P(2,i_pt) > 0 .AND. x%StC_x(3,i_pt) < 0 .AND. dxdt%StC_x(3,i_pt) < 0) THEN
            F_fr(2,i_pt) = p%StC_Y_C_HIGH
         ELSE
            F_fr(2,i_pt) = p%StC_Y_C_LOW
         END IF

         !Brake Y
         IF      ( (x%StC_x(3,i_pt) > p%P_SP(2)-0.2) .AND. (x%StC_x(4,i_pt) > 0) ) THEN
            C_Brake(2,i_pt) = p%StC_Y_C_BRAKE
         ELSE IF ( (x%StC_x(3,i_pt) < p%N_SP(2)+0.2) .AND. (x%StC_x(4,i_pt) < 0) ) THEN
            C_Brake(2,i_pt) = p%StC_Y_C_BRAKE
         ELSE
            C_Brake(2,i_pt) = 0
         END IF

            ! Z
            ! (a)
         IF      (rdisp_P(3,i_pt) > 0 .AND. rdot_P(3,i_pt) < 0 .AND. x%StC_x(5,i_pt) > 0 .AND. dxdt%StC_x(5,i_pt) < 0) THEN
            F_fr(3,i_pt) = p%StC_Z_C_HIGH
            ! (b)
         ELSE IF (rdisp_P(3,i_pt) < 0 .AND. rdot_P(3,i_pt) > 0 .AND. x%StC_x(5,i_pt) < 0 .AND. dxdt%StC_x(5,i_pt) > 0) THEN
            F_fr(3,i_pt) = -p%StC_Z_C_HIGH
            ! (c)
         ELSE IF (rdisp_P(3,i_pt) < 0 .AND. rdot_P(3,i_pt) < 0 .AND. x%StC_x(5,i_pt) > 0 .AND. dxdt%StC_x(5,i_pt) > 0) THEN
            F_fr(3,i_pt) = -p%StC_Z_C_HIGH
         ELSE IF (rdisp_P(3,i_pt) > 0 .AND. rdot_P(3,i_pt) > 0 .AND. x%StC_x(5,i_pt) < 0 .AND. dxdt%StC_x(5,i_pt) < 0) THEN
            F_fr(3,i_pt) = p%StC_Z_C_HIGH
         ELSE
            F_fr(3,i_pt) = p%StC_Z_C_LOW
         END IF

         !Brake Z
         IF      ( (x%StC_x(5,i_pt) > p%P_SP(3)-0.2) .AND. (x%StC_x(6,i_pt) > 0) ) THEN
            C_Brake(3,i_pt) = p%StC_Z_C_BRAKE
         ELSE IF ( (x%StC_x(5,i_pt) < p%N_SP(3)+0.2) .AND. (x%StC_x(6,i_pt) < 0) ) THEN
            C_Brake(3,i_pt) = p%StC_Z_C_BRAKE
         ELSE
            C_Brake(3,i_pt) = 0
         END IF

      ELSE IF (p%StC_CMODE == CMODE_Semi .AND. p%StC_SA_MODE == SA_CMODE_Ph_DF) THEN ! Phase Difference Algorithm with Damping On/Off
            ! X
            ! (a)
         IF      (rdisp_P(1,i_pt) > 0 .AND. rdot_P(1,i_pt) < 0 .AND. x%StC_x(1,i_pt) > 0 .AND. dxdt%StC_x(1,i_pt) < 0) THEN
            C_ctrl(1,i_pt) = p%StC_X_C_HIGH
            ! (b)
         ELSE IF (rdisp_P(1,i_pt) < 0 .AND. rdot_P(1,i_pt) > 0 .AND. x%StC_x(1,i_pt) < 0 .AND. dxdt%StC_x(1,i_pt) > 0) THEN
            C_ctrl(1,i_pt) = p%StC_X_C_HIGH
            ! (c)
         ELSE IF (rdisp_P(1,i_pt) < 0 .AND. rdot_P(1,i_pt) < 0 .AND. x%StC_x(1,i_pt) > 0 .AND. dxdt%StC_x(1,i_pt) > 0) THEN
            C_ctrl(1,i_pt) = p%StC_X_C_HIGH
         ELSE IF (rdisp_P(1,i_pt) > 0 .AND. rdot_P(1,i_pt) > 0 .AND. x%StC_x(1,i_pt) < 0 .AND. dxdt%StC_x(1,i_pt) < 0) THEN
            C_ctrl(1,i_pt) = p%StC_X_C_HIGH
         ELSE
            C_ctrl(1,i_pt) = p%StC_X_C_LOW
         END IF

         !Brake X
         IF      ( (x%StC_x(1,i_pt) > p%P_SP(1)-0.2) .AND. (x%StC_x(2,i_pt) > 0) ) THEN
            C_Brake(1,i_pt) = p%StC_X_C_BRAKE
         ELSE IF ( (x%StC_x(1,i_pt) < p%N_SP(1)+0.2) .AND. (x%StC_x(2,i_pt) < 0) ) THEN
            C_Brake(1,i_pt) = p%StC_X_C_BRAKE
         ELSE
            C_Brake(1,i_pt) = 0
         END IF

            ! Y
            ! (a)
         IF      (rdisp_P(2,i_pt) > 0 .AND. rdot_P(2,i_pt) < 0 .AND. x%StC_x(3,i_pt) > 0 .AND. dxdt%StC_x(3,i_pt) < 0) THEN
            C_ctrl(2,i_pt) = p%StC_Y_C_HIGH
            ! (b)
         ELSE IF (rdisp_P(2,i_pt) < 0 .AND. rdot_P(2,i_pt) > 0 .AND. x%StC_x(3,i_pt) < 0 .AND. dxdt%StC_x(3,i_pt) > 0) THEN
            C_ctrl(2,i_pt) = p%StC_Y_C_HIGH
            ! (c)
         ELSE IF (rdisp_P(2,i_pt) < 0 .AND. rdot_P(2,i_pt) < 0 .AND. x%StC_x(3,i_pt) > 0 .AND. dxdt%StC_x(3,i_pt) > 0) THEN
            C_ctrl(2,i_pt) = p%StC_Y_C_HIGH
         ELSE IF (rdisp_P(2,i_pt) > 0 .AND. rdot_P(2,i_pt) > 0 .AND. x%StC_x(3,i_pt) < 0 .AND. dxdt%StC_x(3,i_pt) < 0) THEN
            C_ctrl(2,i_pt) = p%StC_Y_C_HIGH
         ELSE
            C_ctrl(2,i_pt) = p%StC_Y_C_LOW
         END IF

         !Brake Y
         IF      ( (x%StC_x(3,i_pt) > p%P_SP(2)-0.2) .AND. (x%StC_x(4,i_pt) > 0) ) THEN
            C_Brake(2,i_pt) = p%StC_Y_C_BRAKE
         ELSE IF ( (x%StC_x(3,i_pt) < p%N_SP(2)+0.2) .AND. (x%StC_x(4,i_pt) < 0) ) THEN
            C_Brake(2,i_pt) = p%StC_Y_C_BRAKE
         ELSE
            C_Brake(2,i_pt) = 0
         END IF

            ! Z
            ! (a)
         IF      (rdisp_P(3,i_pt) > 0 .AND. rdot_P(3,i_pt) < 0 .AND. x%StC_x(5,i_pt) > 0 .AND. dxdt%StC_x(5,i_pt) < 0) THEN
            C_ctrl(3,i_pt) = p%StC_Z_C_HIGH
            ! (b)
         ELSE IF (rdisp_P(3,i_pt) < 0 .AND. rdot_P(3,i_pt) > 0 .AND. x%StC_x(5,i_pt) < 0 .AND. dxdt%StC_x(5,i_pt) > 0) THEN
            C_ctrl(3,i_pt) = p%StC_Z_C_HIGH
            ! (c)
         ELSE IF (rdisp_P(3,i_pt) < 0 .AND. rdot_P(3,i_pt) < 0 .AND. x%StC_x(5,i_pt) > 0 .AND. dxdt%StC_x(5,i_pt) > 0) THEN
            C_ctrl(3,i_pt) = p%StC_Z_C_HIGH
         ELSE IF (rdisp_P(3,i_pt) > 0 .AND. rdot_P(3,i_pt) > 0 .AND. x%StC_x(5,i_pt) < 0 .AND. dxdt%StC_x(5,i_pt) < 0) THEN
            C_ctrl(3,i_pt) = p%StC_Z_C_HIGH
         ELSE
            C_ctrl(3,i_pt) = p%StC_Z_C_LOW
         END IF

         !Brake Z
         IF      ( (x%StC_x(5,i_pt) > p%P_SP(3)-0.2) .AND. (x%StC_x(6,i_pt) > 0) ) THEN
            C_Brake(3,i_pt) = p%StC_Z_C_BRAKE
         ELSE IF ( (x%StC_x(5,i_pt) < p%N_SP(3)+0.2) .AND. (x%StC_x(6,i_pt) < 0) ) THEN
            C_Brake(3,i_pt) = p%StC_Z_C_BRAKE
         ELSE
            C_Brake(3,i_pt) = 0
         END IF

      END IF
   enddo
END SUBROUTINE StC_GroundHookDamp
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE StC_ActiveCtrl_StiffDamp(u,p,K_ctrl,C_ctrl,C_Brake,F_ctrl)
   TYPE(StC_InputType),                   INTENT(IN   )  :: u              !< Inputs at Time
   TYPE(StC_ParameterType),               INTENT(IN   )  :: p              !< The module's parameter data
   real(ReKi),                            intent(inout)  :: K_ctrl(:,:)    !< stiffness commanded by dll -- leave alone if no ctrl
   real(ReKi),                            intent(inout)  :: C_ctrl(:,:)    !< damping   commanded by dll
   real(ReKi),                            intent(inout)  :: C_Brake(:,:)   !< brake     commanded by dll
   real(ReKi),                            intent(inout)  :: F_ctrl(:,:)    !< brake     commanded by dll
   integer(IntKi)                                        :: i_pt           ! counter for mesh points
   do i_pt=1,p%NumMeshPts
      if (p%StC_CChan(i_pt) > 0) then     ! This index should have been checked at init, so won't check bounds here
         K_ctrl( 1:3,i_pt) = u%CmdStiff(1:3,p%StC_CChan(i_pt))
         C_ctrl( 1:3,i_pt) = u%CmdDamp( 1:3,p%StC_CChan(i_pt))
         C_Brake(1:3,i_pt) = u%CmdBrake(1:3,p%StC_CChan(i_pt))
         F_ctrl(1:3,i_pt)  = u%CmdForce(1:3,p%StC_CChan(i_pt))
      else  ! Use parameters from file (as if no control) -- leave K value as that may be set by table prior
         C_ctrl(1,:) = p%C_X
         C_ctrl(2,:) = p%C_Y
         C_ctrl(3,:) = p%C_Z
         C_Brake     = 0.0_ReKi
         F_ctrl      = 0.0_ReKi
      endif
   enddo
END SUBROUTINE StC_ActiveCtrl_StiffDamp
!----------------------------------------------------------------------------------------------------------------------------------
!> Extrapolate or interpolate stiffness value based on stiffness table.
SUBROUTINE SpringForceExtrapInterp(x, p, F_table,ErrStat,ErrMsg)
   TYPE(StC_ContinuousStateType),         INTENT(IN   )     :: x           !< Continuous states at Time
   TYPE(StC_ParameterType),               INTENT(IN)        :: p           !< The module's parameter data
   REAL(ReKi), dimension(:,:),            INTENT(INOUT)     :: F_table     !< extrapolated/interpolated stiffness values

   INTEGER(IntKi),                        INTENT(OUT)      :: ErrStat        ! The error status code
   CHARACTER(*),                          INTENT(OUT)      :: ErrMsg         ! The error message, if an error occurred

   ! local variables
   INTEGER(IntKi)                                           :: ErrStat2       ! error status
   INTEGER(IntKi)                                           :: I              ! Loop counter
   INTEGER(IntKi), DIMENSION(3)                             :: J = (/1, 3, 5/) ! Index to StC_x for TMD displacement in each dimension 
   INTEGER(IntKi)                                           :: M              ! location of closest table position
   INTEGER(IntKi)                                           :: Nrows          ! Number of rows in F_TBL
   REAL(ReKi)                                               :: Slope          !
   REAL(ReKi)                                               :: DX             !
   REAL(ReKi)                                               :: Disp(3)        ! Current displacement
   REAL(ReKi), ALLOCATABLE                                  :: TmpRAry(:)
   INTEGER(IntKi)                                           :: i_pt           !< generic counter for mesh point

   ErrStat = ErrID_None
   ErrMsg  = ''

   Nrows = SIZE(p%F_TBL,1)
   ALLOCATE(TmpRAry(Nrows),STAT=ErrStat2)
   IF (ErrStat2 /= 0) then
       call SetErrStat(ErrID_Fatal,'Error allocating temp array.',ErrStat,ErrMsg,'SpringForceExtrapInterp')
      RETURN
   END IF

   do i_pt=1,p%NumMeshPts

      IF (p%StC_DOF_MODE == DOFMode_Indept .OR. p%StC_DOF_MODE == DOFMode_Omni) THEN

         IF (p%StC_DOF_MODE == DOFMode_Indept) THEN
            DO I = 1,3
               Disp(I) = x%StC_x(J(I),i_pt)
            END DO
         ELSE !IF (p%StC_DOF_MODE == DOFMode_Omni) THEN  ! Only X and Y
            Disp = SQRT(x%StC_x(1,i_pt)**2+x%StC_x(3,i_pt)**2) ! constant assignment to vector
            Disp(3) = 0.0_ReKi
         END IF

         
         DO I = 1,3
            TmpRAry = p%F_TBL(:,J(I))-Disp(I)
            TmpRAry = ABS(TmpRAry)
            M = MINLOC(TmpRAry,1)

            !interpolate
            IF ( (Disp(I) > p%F_TBL(M,J(I)) .AND. M /= Nrows) .OR. (Disp(I) < p%F_TBL(M,J(I)) .AND. M == 1) ) THEN
            ! for displacements higher than the closest table value or lower than the lower bound
               Slope = ( p%F_TBL(M+1,J(I)+1) - p%F_TBL(M,J(I)+1) ) / ( p%F_TBL(M+1,J(I)) - p%F_TBL(M,J(I)) )

            ELSE IF ( (Disp(I) < p%F_TBL(M,J(I)) .AND. M /= 1 ) .OR. (Disp(I) > p%F_TBL(M,J(I)) .AND. M == Nrows) ) THEN ! lower
            ! for displacements lower than the closest table value or higher than the upper bound
               Slope = ( p%F_TBL(M,J(I)+1) - p%F_TBL(M-1,J(I)+1) ) / ( p%F_TBL(M,J(I)) - p%F_TBL(M-1,J(I)) )

            ELSE ! equal
               Slope = 0
            END IF

            F_table(I,i_pt) = p%F_TBL(M,J(I)+1) + Slope * ( Disp(I) - p%F_TBL(M,J(I)) )

         END DO

      END IF
   enddo ! Loop over p%NumMeshPts

   DEALLOCATE(TmpRAry)

END SUBROUTINE SpringForceExtrapInterp
!----------------------------------------------------------------------------------------------------------------------------------
!> Parse the inputfile info stored in FileInfo_In.
SUBROUTINE StC_ParseInputFileInfo( PriPath, InputFile, RootName, NumMeshPts, FileInfo_In, InputFileData, UnEcho, ErrStat, ErrMsg )

   implicit    none

      ! Passed variables
   character(*),                    intent(in   )  :: PriPath           !< primary path
   CHARACTER(*),                    intent(in   )  :: InputFile         !< Name of the file containing the primary input data
   CHARACTER(*),                    intent(in   )  :: RootName          !< The rootname of the echo file, possibly opened in this routine
   integer(IntKi),                  intent(in   )  :: NumMeshPts        !< The number of mesh points passed in
   type(StC_InputFile),             intent(inout)  :: InputFileData     !< All the data in the StrucCtrl input file
   type(FileInfoType),              intent(in   )  :: FileInfo_In       !< The derived type for holding the file information.
   integer(IntKi),                  intent(  out)  :: UnEcho            !< The local unit number for this module's echo file
   integer(IntKi),                  intent(  out)  :: ErrStat           !< Error status
   CHARACTER(ErrMsgLen),            intent(  out)  :: ErrMsg            !< Error message

      ! Local variables:
   integer(IntKi)                                  :: i                 !< generic counter
   integer(IntKi)                                  :: ErrStat2          !< Temporary Error status
   character(ErrMsgLen)                            :: ErrMsg2           !< Temporary Error message
   integer(IntKi)                                  :: CurLine           !< current entry in FileInfo_In%Lines array
   real(ReKi)                                      :: TmpRe6(6)         !< temporary 6 number array for reading values in


   ! Initialization
   ErrStat  =  0
   ErrMsg   =  ""
   UnEcho   = -1     ! Echo file unit.  >0 when used


   !-------------------------------------------------------------------------------------------------
   ! General settings
   !-------------------------------------------------------------------------------------------------

   CurLine = 4    ! Skip the first three lines as they are known to be header lines and separators
   call ParseVar( FileInfo_In, CurLine, 'Echo', InputFileData%Echo, ErrStat2, ErrMsg2 )
         if (Failed()) return;

   if ( InputFileData%Echo ) then
      CALL OpenEcho ( UnEcho, TRIM(RootName)//'.ech', ErrStat2, ErrMsg2 )
         if (Failed()) return;
      WRITE(UnEcho, '(A)') 'Echo file for StructCtrl input file: '//trim(InputFile)
      ! Write the first three lines into the echo file
      WRITE(UnEcho, '(A)') FileInfo_In%Lines(1)
      WRITE(UnEcho, '(A)') FileInfo_In%Lines(2)
      WRITE(UnEcho, '(A)') FileInfo_In%Lines(3)

      CurLine = 4
      call ParseVar( FileInfo_In, CurLine, 'Echo', InputFileData%Echo, ErrStat2, ErrMsg2, UnEcho )
            if (Failed()) return
   endif

   !-------------------------------------------------------------------------------------------------
   ! StC DEGREES OF FREEDOM
   !-------------------------------------------------------------------------------------------------

   ! Section break
   if ( InputFileData%Echo )   WRITE(UnEcho, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      !  DOF mode (switch) {  0: No StC or TLCD DOF; 
      !                       1: StC_X_DOF, StC_Y_DOF, and/or StC_Z_DOF (three independent StC DOFs);
      !                       2: StC_XY_DOF (Omni-Directional StC);
      !                       3: TLCD;
      !                       4: Prescribed force/moment time series}
   call ParseVar( FileInfo_In, Curline, 'StC_DOF_MODE', InputFileData%StC_DOF_MODE, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  DOF on or off for StC X (flag) [Used only when StC_DOF_MODE=1]
   call ParseVar( FileInfo_In, Curline, 'StC_X_DOF', InputFileData%StC_X_DOF, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  DOF on or off for StC Y (flag) [Used only when StC_DOF_MODE=1]
   call ParseVar( FileInfo_In, Curline, 'StC_Y_DOF', InputFileData%StC_Y_DOF, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  DOF on or off for StC Z (flag) [Used only when StC_DOF_MODE=1]
   call ParseVar( FileInfo_In, Curline, 'StC_Z_DOF', InputFileData%StC_Z_DOF, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;


   !-------------------------------------------------------------------------------------------------
   ! StC LOCATION [relative to the reference origin of component attached to]
   !-------------------------------------------------------------------------------------------------

   ! Section break
   if ( InputFileData%Echo )   WRITE(UnEcho, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      !  At rest X position of StC(s) (m) [relative to reference origin of the component]
   call ParseVar( FileInfo_In, Curline, 'StC_P_X', InputFileData%StC_P_X, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  At rest Y position of StC(s) (m) [relative to reference origin of the component]
   call ParseVar( FileInfo_In, Curline, 'StC_P_Y', InputFileData%StC_P_Y, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  At rest Z position of StC(s) (m) [relative to reference origin of the component]
   call ParseVar( FileInfo_In, Curline, 'StC_P_Z', InputFileData%StC_P_Z, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;

   !-------------------------------------------------------------------------------------------------
   ! StC INITIAL CONDITIONS [used only when StC_DOF_MODE=1 or 2]
   !-------------------------------------------------------------------------------------------------

   ! Section break
   if ( InputFileData%Echo )   WRITE(UnEcho, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      !  StC X initial displacement (m) [relative to at rest position]
   call ParseVar( FileInfo_In, Curline, 'StC_X_DSP', InputFileData%StC_X_DSP, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  StC Y initial displacement (m) [relative to at rest position]
   call ParseVar( FileInfo_In, Curline, 'StC_Y_DSP', InputFileData%StC_Y_DSP, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  StC Z initial displacement (m) [relative to at rest position; used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
   call ParseVar( FileInfo_In, Curline, 'StC_Z_DSP', InputFileData%StC_Z_DSP, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  StC Z pre-load (N) {"gravity" to offset for gravity load; "none" or 0 to turn off} [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
   call ParseVar( FileInfo_In, Curline, 'StC_Z_PreLd', InputFileData%StC_Z_PreLdC, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;

   !-------------------------------------------------------------------------------------------------
   ! StC CONFIGURATION  [used only when StC_DOF_MODE=1 or 2]
   !-------------------------------------------------------------------------------------------------

   ! Section break
   if ( InputFileData%Echo )   WRITE(UnEcho, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      !  Positive stop position (maximum X mass displacement) (m)
   call ParseVar( FileInfo_In, Curline, 'StC_X_PSP', InputFileData%StC_X_PSP, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  Negative stop position (minimum X mass displacement) (m)
   call ParseVar( FileInfo_In, Curline, 'StC_X_NSP', InputFileData%StC_X_NSP, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  Positive stop position (maximum Y mass displacement) (m)
   call ParseVar( FileInfo_In, Curline, 'StC_Y_PSP', InputFileData%StC_Y_PSP, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  Negative stop position (minimum Y mass displacement) (m)
   call ParseVar( FileInfo_In, Curline, 'StC_Y_NSP', InputFileData%StC_Y_NSP, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  Positive stop position (maximum Z mass displacement) (m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
   call ParseVar( FileInfo_In, Curline, 'StC_Z_PSP', InputFileData%StC_Z_PSP, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  Negative stop position (minimum Z mass displacement) (m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
   call ParseVar( FileInfo_In, Curline, 'StC_Z_NSP', InputFileData%StC_Z_NSP, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;

   !-------------------------------------------------------------------------------------------------
   ! StC MASS, STIFFNESS, & DAMPING  [used only when StC_DOF_MODE=1 or 2]
   !-------------------------------------------------------------------------------------------------

   ! Section break
   if ( InputFileData%Echo )   WRITE(UnEcho, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      !  StC X mass (kg) [must equal StC_Y_M for StC_DOF_MODE = 2]
   call ParseVar( FileInfo_In, Curline, 'StC_X_M', InputFileData%StC_X_M, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  StC Y mass (kg) [must equal StC_X_M for StC_DOF_MODE = 2]
   call ParseVar( FileInfo_In, Curline, 'StC_Y_M', InputFileData%StC_Y_M, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  StC Z mass (kg) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
   call ParseVar( FileInfo_In, Curline, 'StC_Z_M', InputFileData%StC_Z_M, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  StC Z mass (kg) [used only when StC_DOF_MODE=2]
   call ParseVar( FileInfo_In, Curline, 'StC_XY_M', InputFileData%StC_XY_M, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  StC X stiffness (N/m)
   call ParseVar( FileInfo_In, Curline, 'StC_X_K', InputFileData%StC_X_K, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  StC Y stiffness (N/m)
   call ParseVar( FileInfo_In, Curline, 'StC_Y_K', InputFileData%StC_Y_K, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  StC Z stiffness (N/m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
   call ParseVar( FileInfo_In, Curline, 'StC_Z_K', InputFileData%StC_Z_K, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  StC X damping (N/(m/s))
   call ParseVar( FileInfo_In, Curline, 'StC_X_C', InputFileData%StC_X_C, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  StC Y damping (N/(m/s))
   call ParseVar( FileInfo_In, Curline, 'StC_Y_C', InputFileData%StC_Y_C, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  StC Z damping (N/(m/s)) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
   call ParseVar( FileInfo_In, Curline, 'StC_Z_C', InputFileData%StC_Z_C, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  Stop spring X stiffness (N/m)
   call ParseVar( FileInfo_In, Curline, 'StC_X_KS', InputFileData%StC_X_KS, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  Stop spring Y stiffness (N/m)
   call ParseVar( FileInfo_In, Curline, 'StC_Y_KS', InputFileData%StC_Y_KS, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  Stop spring Z stiffness (N/m) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
   call ParseVar( FileInfo_In, Curline, 'StC_Z_KS', InputFileData%StC_Z_KS, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  Stop spring X damping (N/(m/s))
   call ParseVar( FileInfo_In, Curline, 'StC_X_CS', InputFileData%StC_X_CS, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  Stop spring Y damping (N/(m/s))
   call ParseVar( FileInfo_In, Curline, 'StC_Y_CS', InputFileData%StC_Y_CS, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;
      !  Stop spring Z damping (N/(m/s)) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
   call ParseVar( FileInfo_In, Curline, 'StC_Z_CS', InputFileData%StC_Z_CS, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;

   !-------------------------------------------------------------------------------------------------
   ! StC USER-DEFINED SPRING FORCES  [used only when StC_DOF_MODE=1 or 2]
   !-------------------------------------------------------------------------------------------------

   ! Section break
   if ( InputFileData%Echo )   WRITE(UnEcho, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      !  Use spring force from user-defined table (flag)
   call ParseVar( FileInfo_In, Curline, 'Use_F_TBL', InputFileData%Use_F_TBL, ErrStat2, ErrMsg2, UnEcho )
      If (Failed()) return;

      ! NKInpSt      - Number of spring force input stations
   call ParseVar( FileInfo_In, CurLine, 'NKInpSt', InputFileData%NKInpSt, ErrStat2, ErrMsg2, UnEcho)
         if (Failed()) return

   ! Section break --  X  K_X   Y  K_Y   Z  K_Z
   if ( InputFileData%Echo )   WRITE(UnEcho, '(A)') '#TABLE: '//FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
   if ( InputFileData%Echo )   WRITE(UnEcho, '(A)') ' Table Header: '//FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1
   if ( InputFileData%Echo )   WRITE(UnEcho, '(A)') ' Table Units: '//FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

   if (InputFileData%NKInpSt > 0) then
      CALL AllocAry( InputFileData%F_TBL, InputFileData%NKInpSt, 6, 'F_TBL', ErrStat2, ErrMsg2 )
            if (Failed()) return;
         ! TABLE read
      do i=1,InputFileData%NKInpSt
         call ParseAry ( FileInfo_In, CurLine, 'Coordinates', TmpRe6, 6, ErrStat2, ErrMsg2, UnEcho )
               if (Failed()) return;
         InputFileData%F_TBL(i,1) = TmpRe6(1) ! X
         InputFileData%F_TBL(i,2) = TmpRe6(2) ! K_X
         InputFileData%F_TBL(i,3) = TmpRe6(3) ! Y
         InputFileData%F_TBL(i,4) = TmpRe6(4) ! K_Y
         InputFileData%F_TBL(i,5) = TmpRe6(5) ! Z
         InputFileData%F_TBL(i,6) = TmpRe6(6) ! K_Z
      enddo
   endif


   !-------------------------------------------------------------------------------------------------
   ! StructCtrl CONTROL  [used only when StC_DOF_MODE=1 or 2]
   !-------------------------------------------------------------------------------------------------

   ! Section break
   if ( InputFileData%Echo )   WRITE(UnEcho, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      !  Control mode (switch) {
      !     0:none;
      !     1: Semi-Active Control Mode;
      !     4: Active Control Mode through Simulink (not available);
      !     5: Active Control Mode through Bladed interface} (-)
   call ParseVar( FileInfo_In, Curline, 'StC_CMODE', InputFileData%StC_CMODE, ErrStat2, ErrMsg2 )
      If (Failed()) return;
   ! Control channels -- there may be multiple if there are multiple StC mesh points (blade TMD case), but we also allow a single
      ! StC_CChan     - Control channel group for stiffness and damping (StC_[XYZ]_K, StC_[XYZ]_C, and StC_[XYZ]_Brake) [used only when StC_CMODE=4 or StC_CMODE=5]
   allocate( InputFileData%StC_CChan(NumMeshPts), STAT=ErrStat2 )    ! Blade TMD will possibly have independent TMD's for each instance
      if (ErrStat2 /= ErrID_None) ErrMsg2="Error allocating InputFileData%StC_CChan(NumMeshPts)"
      If (Failed()) return;
   call ParseAry( FileInfo_In, CurLine, 'StC_CChan', InputFileData%StC_CChan, NumMeshPts, ErrStat2, ErrMsg2 )
   if ( ErrStat2 /= ErrID_None) then      ! If we didn't read a full array, then try reading just one input
      ! If there was an error, CurLine didn't advance, so no resetting of it needed
      call ParseVar( FileInfo_In, CurLine, 'StC_CChan', InputFileData%StC_CChan(1), ErrStat2, ErrMsg2 )
      InputFileData%StC_CChan(:) = InputFileData%StC_CChan(1)     ! Assign all the same, will check in validation
   endif
      If (Failed()) return;

      !  Semi-Active control mode {
      !     1: velocity-based ground hook control; 
      !     2: Inverse velocity-based ground hook control;
      !     3: displacement-based ground hook control;
      !     4: Phase difference Algorithm with Friction Force;
      !     5: Phase difference Algorithm with Damping Force} (-)
   call ParseVar( FileInfo_In, Curline, 'StC_SA_MODE', InputFileData%StC_SA_MODE, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  StC X high damping for ground hook control
   call ParseVar( FileInfo_In, Curline, 'StC_X_C_HIGH', InputFileData%StC_X_C_HIGH, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  StC X low damping for ground hook control
   call ParseVar( FileInfo_In, Curline, 'StC_X_C_LOW', InputFileData%StC_X_C_LOW, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  StC Y high damping for ground hook control
   call ParseVar( FileInfo_In, Curline, 'StC_Y_C_HIGH', InputFileData%StC_Y_C_HIGH, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  StC Y low damping for ground hook control
   call ParseVar( FileInfo_In, Curline, 'StC_Y_C_LOW', InputFileData%StC_Y_C_LOW, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  StC Z high damping for ground hook control [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
   call ParseVar( FileInfo_In, Curline, 'StC_Z_C_HIGH', InputFileData%StC_Z_C_HIGH, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  StC Z low damping for ground hook control  [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
   call ParseVar( FileInfo_In, Curline, 'StC_Z_C_LOW', InputFileData%StC_Z_C_LOW, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  StC X high damping for braking the StC (Don't use it now. should be zero)
   call ParseVar( FileInfo_In, Curline, 'StC_X_C_BRAKE', InputFileData%StC_X_C_BRAKE, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  StC Y high damping for braking the StC (Don't use it now. should be zero)
   call ParseVar( FileInfo_In, Curline, 'StC_Y_C_BRAKE', InputFileData%StC_Y_C_BRAKE, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  StC Z high damping for braking the StC (Don't use it now. should be zero) [used only when StC_DOF_MODE=1 and StC_Z_DOF=TRUE]
   call ParseVar( FileInfo_In, Curline, 'StC_Z_C_BRAKE', InputFileData%StC_Z_C_BRAKE, ErrStat2, ErrMsg2 )
      If (Failed()) return;

   !-------------------------------------------------------------------------------------------------
   ! TLCD  [used only when StC_DOF_MODE=3]
   !-------------------------------------------------------------------------------------------------

   ! Section break
   if ( InputFileData%Echo )   WRITE(UnEcho, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      !  X TLCD total length (m)
   call ParseVar( FileInfo_In, Curline, 'L_X', InputFileData%L_X, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  X TLCD horizontal length (m)
   call ParseVar( FileInfo_In, Curline, 'B_X', InputFileData%B_X, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  X TLCD cross-sectional area of vertical column (m^2)
   call ParseVar( FileInfo_In, Curline, 'area_X', InputFileData%area_X, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  X TLCD cross-sectional area ratio (vertical column area divided by horizontal column area) (-)
   call ParseVar( FileInfo_In, Curline, 'area_ratio_X', InputFileData%area_ratio_X, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  X TLCD head loss coeff (-)
   call ParseVar( FileInfo_In, Curline, 'headLossCoeff_X', InputFileData%headLossCoeff_X, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  X TLCD liquid density (kg/m^3)
   call ParseVar( FileInfo_In, Curline, 'rho_X', InputFileData%rho_X, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  Y TLCD total length (m)
   call ParseVar( FileInfo_In, Curline, 'L_Y', InputFileData%L_Y, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  Y TLCD horizontal length (m)
   call ParseVar( FileInfo_In, Curline, 'B_Y', InputFileData%B_Y, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  Y TLCD cross-sectional area of vertical column (m^2)
   call ParseVar( FileInfo_In, Curline, 'area_Y', InputFileData%area_Y, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  Y TLCD cross-sectional area ratio (vertical column area divided by horizontal column area) (-)
   call ParseVar( FileInfo_In, Curline, 'area_ratio_Y', InputFileData%area_ratio_Y, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  Y TLCD head loss coeff (-)
   call ParseVar( FileInfo_In, Curline, 'headLossCoeff_Y', InputFileData%headLossCoeff_Y, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      !  Y TLCD liquid density (kg/m^3)
   call ParseVar( FileInfo_In, Curline, 'rho_Y', InputFileData%rho_Y, ErrStat2, ErrMsg2 )
      If (Failed()) return;

   !-------------------------------------------------------------------------------------------------
   ! PRESCRIBED TIME SERIES  [used only when StC_DOF_MODE=4]
   !-------------------------------------------------------------------------------------------------

   ! Section break
   if ( InputFileData%Echo )   WRITE(UnEcho, '(A)') FileInfo_In%Lines(CurLine)    ! Write section break to echo
   CurLine = CurLine + 1

      ! Prescribed forces coordinate system
   call ParseVar( FileInfo_In, Curline, 'PrescribedForcesCoordSys', InputFileData%PrescribedForcesCoordSys, ErrStat2, ErrMsg2 )
      If (Failed()) return;
      ! Prescribed input time series
   call ParseVar( FileInfo_In, Curline, 'PrescribedForcesFile', InputFileData%PrescribedForcesFile, ErrStat2, ErrMsg2 )
      if (Failed()) return;
      if ( PathIsRelative( InputFileData%PrescribedForcesFile ) ) InputFileData%PrescribedForcesFile = TRIM(PriPath)//TRIM(InputFileData%PrescribedForcesFile)


CONTAINS
   !-------------------------------------------------------------------------------------------------
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'StC_ParseInputFileInfo' )
      Failed = ErrStat >= AbortErrLev
      if (Failed) call Cleanup()
   end function Failed
   !-------------------------------------------------------------------------------------------------
   SUBROUTINE Cleanup()
      if (UnEcho  > -1_IntKi)     CLOSE( UnEcho  )
   END SUBROUTINE Cleanup
   !-------------------------------------------------------------------------------------------------
END SUBROUTINE StC_ParseInputFileInfo

!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine checks the data handed in.  If all is good, no errors reported. 
subroutine    StC_ValidatePrimaryData( InputFileData, InitInp, ErrStat, ErrMsg )
   TYPE(StC_InputFile),      INTENT(IN)      :: InputFileData  !< Data stored in the module's input file
   TYPE(StC_InitInputType),  INTENT(IN   )   :: InitInp        !< Input data for initialization routine.
   INTEGER(IntKi),           INTENT(  OUT)   :: ErrStat        !< The error status code
   CHARACTER(ErrMsgLen),     INTENT(  OUT)   :: ErrMsg         !< The error message, if an error occurred

   integer(IntKi)                            :: i              !< generic loop counter
   real(ReKi)                                :: TmpRe
   character(10)                             :: TmpCh
   integer(IntKi)                            :: ErrStat2
   CHARACTER(*), PARAMETER                   :: RoutineName = 'StC_ValidatePrimaryData'

      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ''

      ! Check DOF modes
   IF (  InputFileData%StC_DOF_MODE /= ControlMode_None     .and. &
         InputFileData%StC_DOF_MODE /= DOFMode_Indept       .and. &
         InputFileData%StC_DOF_MODE /= DOFMode_Omni         .and. &
         InputFileData%StC_DOF_MODE /= DOFMode_TLCD         .and. &
         InputFileData%StC_DOF_MODE /= DOFMode_Prescribed   .and. &
         InputFileData%StC_DOF_MODE /= DOFMode_ForceDLL) &
      CALL SetErrStat( ErrID_Fatal, 'DOF mode (StC_DOF_MODE) must be 0 (no DOF), 1 (two independent DOFs), '// &
               'or 2 (omni-directional), or 3 (TLCD), 4 (prescribed force time-series), or 5 (force from external DLL).', ErrStat, ErrMsg, RoutineName )

      ! Check control modes
   IF (  InputFileData%StC_CMODE /= ControlMode_None     .and. &
         InputFileData%StC_CMODE /= CMODE_Semi           .and. &
         InputFileData%StC_CMODE /= CMODE_ActiveDLL ) &
         !InputFileData%StC_CMode /= CMODE_ActiveEXTERN   .and. &    ! Not an option at the moment --> 4 (active with Simulink control),
      CALL SetErrStat( ErrID_Fatal, 'Control mode (StC_CMode) must be 0 (none), 1 (semi-active),'// &
            ' or 5 (active with DLL control) in this version of StrucCtrl.', ErrStat, ErrMsg, RoutineName )

      ! Check control channel
   if ( InputFileData%StC_CMode == CMODE_ActiveDLL ) then
      if ( InputFileData%StC_DOF_MODE /= DOFMode_Indept .and. &
           InputFileData%StC_DOF_MODE /= DOFMode_Omni   .and. &
           InputFileData%StC_DOF_MODE /= DOFMode_ForceDLL) then
         call SetErrStat( ErrID_Fatal, 'Control mode 4 (active with Simulink control), or 5 (active with DLL control) '// &
               'can only be used with independent or omni DOF (StC_DOF_Mode=1 or 2) or force from external DLL '// &
               '(StC_DOF_Mode = 5) in this version of StrucCtrl.', ErrStat, ErrMsg, RoutineName )
      endif
      if (InitInp%NumMeshPts > 1) then
         do i=2,InitInp%NumMeshPts  ! Warn if controlling multiple mesh points with single instance (blade TMD)
            if ( InputFileData%StC_CChan(i) == InputFileData%StC_CChan(1) ) then
               call SetErrStat( ErrID_Warn, 'Same control channel is used for multiple StC instances within input file '// &
                        trim(InitInp%InputFile)//'.  This may not be desired in some cases such as blade TMD active controls.', &
                        ErrStat, ErrMsg, RoutineName )
            endif
         enddo
      endif
      do i=1,InitInp%NumMeshPts     ! Check we are in range of number of control channel groups
         if ( InputFileData%StC_CChan(i) < 0 .or. InputFileData%StC_CChan(i) > 10 ) then
            call SetErrStat( ErrID_Fatal, 'Control channel (StC_CChan) must be between 0 (off) and 10 when StC_CMode=5.', ErrStat, ErrMsg, RoutineName )
         endif
      enddo
   endif

   IF ( InputFileData%StC_SA_MODE /= SA_CMODE_GH_vel    .and. &
        InputFileData%StC_SA_MODE /= SA_CMODE_GH_invVel .and. &
        InputFileData%StC_SA_MODE /= SA_CMODE_GH_disp   .and. &
        InputFileData%StC_SA_MODE /= SA_CMODE_Ph_FF     .and. &
        InputFileData%StC_SA_MODE /= SA_CMODE_Ph_DF     ) then
      CALL SetErrStat( ErrID_Fatal, 'Semi-active control mode (StC_SA_MODE) must be 1 (velocity-based ground hook control), '// &
                   '2 (inverse velocity-based ground hook control), 3 (displacement-based ground hook control), '// &
                   '4 (phase difference algorithm with friction force), or 5 (phase difference algorithm with damping force).', ErrStat, ErrMsg, RoutineName )
   END IF

      ! Prescribed forces
   if (InputFileData%StC_DOF_MODE == DOFMode_Prescribed) then
      if (InputFileData%PrescribedForcesCoordSys /= PRESCRIBED_FORCE_GLOBAL .and. InputFileData%PrescribedForcesCoordSys /= PRESCRIBED_FORCE_LOCAL) then
         call SetErrStat( ErrID_Fatal, 'PrescribedForcesCoordSys must be '//trim(Num2LStr(PRESCRIBED_FORCE_GLOBAL))//   &
                                 ' (Global) or '//trim(Num2LStr(PRESCRIBED_FORCE_LOCAL))//' (local)', ErrStat, ErrMsg, RoutineName )
      endif
   endif

      ! DLL Force - not sure if necessary, but nothing happens if these inputs are incorrect
   if (InputFileData%StC_DOF_MODE == DOFMode_ForceDLL) then
      
      ! Need global force coord
      if (InputFileData%PrescribedForcesCoordSys /= PRESCRIBED_FORCE_GLOBAL) THEN
         call SetErrStat( ErrID_Fatal, 'PrescribedForcesCoordSys must be global ('//trim(Num2LStr(PRESCRIBED_FORCE_GLOBAL))//   &
                                 ') when StC_DOF_MODE is '//trim(Num2LStr(DOFMode_ForceDLL)) , ErrStat, ErrMsg, RoutineName )
      endif

      ! Need active DLL control
      if (InputFileData%StC_CMODE /= CMODE_ActiveDLL) THEN
         call SetErrStat( ErrID_Fatal, 'StC_CMODE must be '//trim(Num2LStr(CMODE_ActiveDLL))//   &
                                 ' when StC_DOF_MODE is '//trim(Num2LStr(DOFMode_ForceDLL)) , ErrStat, ErrMsg, RoutineName )
      endif
   endif

      ! Check masses make some kind of sense
   if (InputFileData%StC_DOF_MODE == DOFMode_Indept .and. InputFileData%StC_X_DOF .and. (InputFileData%StC_X_M <= 0.0_ReKi) )    & 
      call SetErrStat(ErrID_Fatal,'StC_X_M must be > 0 when StC_X_DOF is enabled', ErrStat,ErrMsg,RoutineName)
   if (InputFileData%StC_DOF_MODE == DOFMode_Indept .and. InputFileData%StC_X_DOF .and. (InputFileData%StC_X_K <= 0.0_ReKi) )    & 
      call SetErrStat(ErrID_Fatal,'StC_X_K must be > 0 when StC_X_DOF is enabled', ErrStat,ErrMsg,RoutineName)

   if (InputFileData%StC_DOF_MODE == DOFMode_Indept .and. InputFileData%StC_Y_DOF .and. (InputFileData%StC_Y_M <= 0.0_ReKi) )    & 
      call SetErrStat(ErrID_Fatal,'StC_Y_M must be > 0 when StC_Y_DOF is enabled', ErrStat,ErrMsg,RoutineName)
   if (InputFileData%StC_DOF_MODE == DOFMode_Indept .and. InputFileData%StC_Y_DOF .and. (InputFileData%StC_Y_K <= 0.0_ReKi) )    & 
      call SetErrStat(ErrID_Fatal,'StC_Y_K must be > 0 when StC_Y_DOF is enabled', ErrStat,ErrMsg,RoutineName)

   if (InputFileData%StC_DOF_MODE == DOFMode_Indept .and. InputFileData%StC_Z_DOF .and. (InputFileData%StC_Z_M <= 0.0_ReKi) )    & 
      call SetErrStat(ErrID_Fatal,'StC_Z_M must be > 0 when StC_Z_DOF is enabled', ErrStat,ErrMsg,RoutineName)
   if (InputFileData%StC_DOF_MODE == DOFMode_Indept .and. InputFileData%StC_Z_DOF .and. (InputFileData%StC_Z_K <= 0.0_ReKi) )    & 
      call SetErrStat(ErrID_Fatal,'StC_Z_K must be > 0 when StC_Z_DOF is enabled', ErrStat,ErrMsg,RoutineName)

   if (InputFileData%StC_DOF_MODE == DOFMode_Omni .and. (InputFileData%StC_XY_M <= 0.0_ReKi) )    & 
      call SetErrStat(ErrID_Fatal,'StC_XY_M must be > 0 when DOF mode 2 (omni-directional) is used', ErrStat,ErrMsg,RoutineName)
   if (InputFileData%StC_DOF_MODE == DOFMode_Omni .and. (InputFileData%StC_X_K <= 0.0_ReKi) )    & 
      call SetErrStat(ErrID_Fatal,'StC_X_K must be > 0 when DOF mode 2 (omni-directional) is used', ErrStat,ErrMsg,RoutineName)
   if (InputFileData%StC_DOF_MODE == DOFMode_Omni .and. (InputFileData%StC_Y_K <= 0.0_ReKi) )    & 
      call SetErrStat(ErrID_Fatal,'StC_Y_K must be > 0 when DOF mode 2 (omni-directional) is used', ErrStat,ErrMsg,RoutineName)

      ! Check spring preload in ZDof
   TmpCh = trim(InputFileData%StC_Z_PreLdC)
   call Conv2UC(TmpCh)
   if ( INDEX(TmpCh, "GRAVITY") /= 1 ) then  ! if not gravity, check that it reads ok
      if ( INDEX(TmpCh, "NONE") /= 1 ) then  ! if not NONE, check that it reads ok
         READ (InputFileData%StC_Z_PreLdC,*,IOSTAT=ErrStat2)   TmpRe    ! attempt to read real number
         if (ErrStat2 /= 0) then
            call SetErrStat(ErrID_Fatal,'StC_Z_PreLd must be "gravity", "none" or a real number', ErrStat,ErrMsg,RoutineName)
         endif
      endif
   endif

      ! Sanity checks for the TLCD option
!FIXME: add some sanity checks here

end subroutine StC_ValidatePrimaryData
!----------------------------------------------------------------------------------------------------------------------------------
!> This subroutine sets the parameters, based on the data stored in InputFileData.
SUBROUTINE StC_SetParameters( InputFileData, InitInp, p, Interval, ErrStat, ErrMsg )
!..................................................................................................................................

   TYPE(StC_InputFile),      INTENT(IN   )   :: InputFileData  !< Data stored in the module's input file
   TYPE(StC_InitInputType),  INTENT(IN   )   :: InitInp        !< Input data for initialization routine.
   TYPE(StC_ParameterType),  INTENT(INOUT)   :: p              !< The module's parameter data
   REAL(DbKi),               INTENT(IN   )   :: Interval       !< Coupling interval in seconds: the rate that
   INTEGER(IntKi),           INTENT(  OUT)   :: ErrStat        !< The error status code
   CHARACTER(ErrMsgLen),     INTENT(  OUT)   :: ErrMsg         !< The error message, if an error occurred

      ! Local variables
   character(10)                             :: TmpCh          ! Temporary string for figuring out what to do with preload on Z spring
   INTEGER(IntKi)                            :: ErrStat2       ! Temporary error ID
   CHARACTER(ErrMsgLen)                      :: ErrMsg2        ! Temporary message describing error
   CHARACTER(*), PARAMETER                   :: RoutineName = 'StC_SetParameters'


      ! Initialize variables
   ErrStat = ErrID_None
   ErrMsg  = ''

      ! Filenames
   p%RootName     =  TRIM(InitInp%RootName)     ! Already includes NStC, TStC, or BStC

      ! Constants
   p%DT  = Interval
   p%Gravity = InitInp%Gravity      ! Gravity vector pointed in negative global Z-axis (/0,0,-g/)
   p%NumMeshPts   =  InitInp%NumMeshPts

      ! DOF controls
   p%StC_DOF_MODE = InputFileData%StC_DOF_MODE

   !p%DT = InputFileData%DT
   !p%RootName = 'StC'
   ! DOFs

   p%StC_X_DOF = InputFileData%StC_X_DOF
   p%StC_Y_DOF = InputFileData%StC_Y_DOF
   if (p%StC_DOF_MODE == DOFMode_Indept) then
      p%StC_Z_DOF = InputFileData%StC_Z_DOF
   else
      p%StC_Z_DOF = .false.
   endif

   ! StC X parameters
   p%M_X = InputFileData%StC_X_M
   p%K_X = InputFileData%StC_X_K
   p%C_X = InputFileData%StC_X_C

   ! StC Y parameters
   p%M_Y = InputFileData%StC_Y_M
   p%K_Y = InputFileData%StC_Y_K
   p%C_Y = InputFileData%StC_Y_C

   ! StC Z parameters
   p%M_Z = InputFileData%StC_Z_M
   p%K_Z = InputFileData%StC_Z_K
   p%C_Z = InputFileData%StC_Z_C

   ! StC Omni parameters
   p%M_XY = InputFileData%StC_XY_M

   ! Fore-Aft TLCD Parameters ! MEG & SP
   p%L_X = InputFileData%L_X
   p%B_X = InputFileData%B_X
   p%area_X = InputFileData%area_X
   p%area_ratio_X = InputFileData%area_ratio_X
   p%headLossCoeff_X = InputFileData%headLossCoeff_X
   p%rho_X = InputFileData%rho_X

   !Side-Side TLCD Parameters
   p%L_Y = InputFileData%L_Y
   p%B_Y = InputFileData%B_Y
   p%area_Y = InputFileData%area_Y
   p%area_ratio_Y = InputFileData%area_ratio_Y
   p%headLossCoeff_Y = InputFileData%headLossCoeff_Y
   p%rho_Y = InputFileData%rho_Y ! MEG & SP

     ! vector parameters
   ! stop positions
   p%P_SP(1) = InputFileData%StC_X_PSP
   p%P_SP(2) = InputFileData%StC_Y_PSP
   p%P_SP(3) = InputFileData%StC_Z_PSP
   p%N_SP(1) = InputFileData%StC_X_NSP
   p%N_SP(2) = InputFileData%StC_Y_NSP
   p%N_SP(3) = InputFileData%StC_Z_NSP
   ! stop force stiffness
   p%K_S(1) = InputFileData%StC_X_KS
   p%K_S(2) = InputFileData%StC_Y_KS
   p%K_S(3) = InputFileData%StC_Z_KS
   ! stop force damping
   p%C_S(1) = InputFileData%StC_X_CS
   p%C_S(2) = InputFileData%StC_Y_CS
   p%C_S(3) = InputFileData%StC_Z_CS
 
   ! ground hook control damping files
   p%StC_CMODE = InputFileData%StC_CMODE
   p%StC_SA_MODE = InputFileData%StC_SA_MODE
   p%StC_X_C_HIGH = InputFileData%StC_X_C_HIGH
   p%StC_X_C_LOW = InputFileData%StC_X_C_LOW
   p%StC_Y_C_HIGH = InputFileData%StC_Y_C_HIGH
   p%StC_Y_C_LOW = InputFileData%StC_Y_C_LOW
   p%StC_Z_C_HIGH = InputFileData%StC_Z_C_HIGH
   p%StC_Z_C_LOW = InputFileData%StC_Z_C_LOW
   p%StC_X_C_BRAKE = InputFileData%StC_X_C_BRAKE
   p%StC_Y_C_BRAKE = InputFileData%StC_Y_C_BRAKE
   p%StC_Z_C_BRAKE = InputFileData%StC_Z_C_BRAKE

   ! User Defined Stiffness Table
   p%Use_F_TBL = InputFileData%Use_F_TBL
   if (allocated(InputFileData%F_TBL)) then
      call AllocAry(p%F_TBL,SIZE(InputFiledata%F_TBL,1),SIZE(InputFiledata%F_TBL,2),'F_TBL', ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName);  if (ErrStat >= ErrID_Fatal) return
      p%F_TBL = InputFileData%F_TBL
   endif

   ! Spring preload in ZDof
   p%StC_Z_PreLd = 0.0_ReKi
   TmpCh = trim(InputFileData%StC_Z_PreLdC)
   call Conv2UC(TmpCh)
   if (InputFileData%StC_DOF_MODE == DOFMode_Indept .and. InputFileData%StC_Z_DOF .and. &
         (INDEX(TmpCh, "NONE") /= 1) ) then
      if (INDEX(TmpCh, "GRAVITY") == 1 ) then  ! if gravity, then calculate
         p%StC_Z_PreLd = -p%Gravity(3)*p%M_Z
      else
         READ (InputFileData%StC_Z_PreLdC,*,IOSTAT=ErrStat2)   p%StC_Z_PreLd     ! Read a real number and store
         if (ErrStat2 /= 0) call SetErrStat(ErrID_Fatal,'StC_Z_PreLd must be "gravity", "none" or a real number', ErrStat,ErrMsg,RoutineName)
         if (ErrStat >= ErrID_Fatal) return
      endif
   endif

   ! Prescribed forces
   p%PrescribedForcesCoordSys =  InputFileData%PrescribedForcesCoordSys
   if (allocated(InputFileData%StC_PrescribedForce)) then
      call AllocAry( p%StC_PrescribedForce, size(InputFileData%StC_PrescribedForce,1), size(InputFileData%StC_PrescribedForce,2),"Array of force data", ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName);  if (ErrStat >= ErrID_Fatal) return
      p%StC_PrescribedForce = InputFileData%StC_PrescribedForce
   endif

   ! StC Control channels
   call AllocAry( p%StC_CChan, p%NumMeshPts, 'p%StC_CChan', ErrStat2, ErrMsg2 )
   if (p%StC_CMODE == CMODE_ActiveDLL ) then
      p%StC_CChan = InputFileData%StC_CChan
   else
      p%StC_CChan = 0   ! turn off regardless of input file request.
   endif

END SUBROUTINE StC_SetParameters


subroutine StC_ParseTimeSeriesFileInfo( InputFile, FileInfo_In, InputFileData, UnEcho, ErrStat, ErrMsg )

   implicit    none

      ! Passed variables
   CHARACTER(*),           intent(in   )  :: InputFile         !< Name of the file containing the primary input data
   type(StC_InputFile),    intent(inout)  :: InputFileData     !< All the data in the StrucCtrl input file
   type(FileInfoType),     intent(in   )  :: FileInfo_In       !< The derived type for holding the file information.
   integer(IntKi),         intent(inout)  :: UnEcho            !< The local unit number for this module's echo file
   integer(IntKi),         intent(  out)  :: ErrStat           !< Error status
   CHARACTER(ErrMsgLen),   intent(  out)  :: ErrMsg            !< Error message

      ! Local variables:
   integer(IntKi)                         :: i                 !< generic counter
   integer(IntKi)                         :: ErrStat2          !< Temporary Error status
   character(ErrMsgLen)                   :: ErrMsg2           !< Temporary Error message
   integer(IntKi)                         :: CurLine           !< current entry in FileInfo_In%Lines array
   real(ReKi)                             :: TmpRe7(7)         !< temporary 7 number array for reading values in
   character(*), parameter                :: RoutineName='StC_ParseTimeSeriesFileInfo'

      ! Initialization of subroutine
   ErrMsg      =  ''
   ErrMsg2     =  ''
   ErrStat     =  ErrID_None
   ErrStat2    =  ErrID_None

   !  This file should only contain a table.  Header lines etc should be commented out. Any blank lines at the
   !  end get removed by the ProcessCom
   call AllocAry( InputFileData%StC_PrescribedForce, 7, FileInfo_In%NumLines, "Array of force data", ErrStat2, ErrMsg2 )
      if (Failed())  return;

   ! Loop over all table lines.  Expecting 7 colunns
   CurLine=1
   do i=1,FileInfo_In%NumLines
      call ParseAry ( FileInfo_In, CurLine, 'Coordinates', TmpRe7, 7, ErrStat2, ErrMsg2, UnEcho )
         if (Failed())  return;
      InputFileData%StC_PrescribedForce(1:7,i) = TmpRe7
   enddo

contains
   !-------------------------------------------------------------------------------------------------
   logical function Failed()
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      Failed = ErrStat >= AbortErrLev
   end function Failed
end subroutine StC_ParseTimeSeriesFileInfo

!----------------------------------------------------------------------------------------------------------------------------------
END MODULE StrucCtrl
!**********************************************************************************************************************************
