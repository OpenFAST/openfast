module AeroDyn
    
   use NWTC_Library
   use AeroDyn_Types
   use BEMT
   use AirfoilInfo
   use BladeElement, only : SkewMod_Uncoupled, SkewMod_PittPeters, SkewMod_Coupled
   use NWTC_LAPACK
   
   implicit none

   private
   
   type(ProgDesc), parameter  :: AD_Ver = ProgDesc( 'AeroDyn', 'v15.00.00a-bjj', '29-May-2015' )
   character(*),   parameter  :: AD_Nickname = 'AD'
   integer(IntKi), parameter  :: AD_numChanPerNode = 15  ! TODO This needs to be set dynamically ?? 9/18/14 GJH
   
   
   INTEGER(IntKi), PARAMETER        :: MaxBl    =  3                                   ! Maximum number of blades allowed in simulation
   INTEGER(IntKi), PARAMETER        :: MaxOutPts = 1  !bjj: fix me!!!!!
   
   ! model identifiers
   integer(intKi), parameter        :: ModelUnknown  = -1
   
   integer(intKi), parameter        :: WakeMod_none  = 0  
   integer(intKi), parameter        :: WakeMod_BEMT  = 1
   
   integer(intKi), parameter        :: AFAeroMod_steady = 1       ! steady model
   integer(intKi), parameter        :: AFAeroMod_BL_unsteady = 2  ! Beddoes-Leishman unsteady model
   
   integer(intKi), parameter        :: TwrInflnc_none = 0  ! none
   integer(intKi), parameter        :: TwrInflnc_Potent = 1  ! potential flow
   integer(intKi), parameter        :: TwrInflnc_Shadow = 2  ! tower shadow 
   
   ! ..... Public Subroutines ...................................................................................................

   public :: AD_Init                           ! Initialization routine
   public :: AD_End                            ! Ending routine (includes clean up)

   public :: AD_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                               !   continuous states, and updating discrete states
   public :: AD_CalcOutput                     ! Routine for computing outputs

   public :: AD_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   
  
contains    
subroutine AD_SetInitOut(p, InitOut, errStat, errMsg)

   type(AD_InitOutputType),       intent(  out)  :: InitOut          ! output data
   type(AD_ParameterType),        intent(in   )  :: p                ! Parameters
   integer(IntKi),                intent(inout)  :: errStat          ! Error status of the operation
   character(*),                  intent(inout)  :: errMsg           ! Error message if ErrStat /= ErrID_None


      ! Local variables
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'AD_SetInitOut'
   
   
   
   integer(IntKi)                                :: j, k
   character(len=10) :: chanPrefix
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   
   call AllocAry( InitOut%WriteOutputHdr, p%numOuts, 'WriteOutputHdr', errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   
   call AllocAry( InitOut%WriteOutputUnt, p%numOuts, 'WriteOutputUnt', errStat2, errMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

      
   if (ErrStat >= AbortErrLev) return
   
      ! Loop over blades and nodes to populate the output channel names and units
   
   do k=1,p%numBlades
      do j=1,p%NumBlNds
         
         chanPrefix = "B"//trim(num2lstr(k))//"N"//trim(num2lstr(j))
         InitOut%WriteOutputHdr( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 1 ) = trim(chanPrefix)//"Theta"
         InitOut%WriteOutputUnt( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 1 ) = '  (deg)  '
         InitOut%WriteOutputHdr( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 2 ) = trim(chanPrefix)//"Psi"
         InitOut%WriteOutputUnt( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 2 ) = '  (deg)  '
         InitOut%WriteOutputHdr( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 3 ) = trim(chanPrefix)//"Vx"
         InitOut%WriteOutputUnt( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 3 ) = '  (m/s)  '
         InitOut%WriteOutputHdr( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 4 ) = trim(chanPrefix)//"Vy"
         InitOut%WriteOutputUnt( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 4 ) = '  (m/s)  '
         InitOut%WriteOutputHdr( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 5 ) = ' '//trim(chanPrefix)//"AxInd"
         InitOut%WriteOutputUnt( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 5 ) = '  (deg)  '
         InitOut%WriteOutputHdr( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 6 ) = ' '//trim(chanPrefix)//"TanInd"
         InitOut%WriteOutputUnt( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 6 ) = '  (deg)  '
         InitOut%WriteOutputHdr( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 7 ) = trim(chanPrefix)//"IndVel"
         InitOut%WriteOutputUnt( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 7 ) = '  (m/s)  '
         InitOut%WriteOutputHdr( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 8 ) = ' '//trim(chanPrefix)//"Phi"
         InitOut%WriteOutputUnt( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 8 ) = '  (deg)  '
         InitOut%WriteOutputHdr( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 9 ) = ' '//trim(chanPrefix)//"AOA"
         InitOut%WriteOutputUnt( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 9 ) = '  (deg)  '
         InitOut%WriteOutputHdr( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 10 ) = ' '//trim(chanPrefix)//"Cl"
         InitOut%WriteOutputUnt( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 10 ) = '   (-)   '
         InitOut%WriteOutputHdr( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 11 ) = ' '//trim(chanPrefix)//"Cd"
         InitOut%WriteOutputUnt( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 11 ) = '   (-)   '
         InitOut%WriteOutputHdr( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 12 ) = ' '//trim(chanPrefix)//"Cx"
         InitOut%WriteOutputUnt( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 12 ) = '   (-)   '
         InitOut%WriteOutputHdr( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 13 ) = ' '//trim(chanPrefix)//"Cy"
         InitOut%WriteOutputUnt( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 13 ) = '   (-)   '
         InitOut%WriteOutputHdr( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 14 ) = ' '//trim(chanPrefix)//"Fx"
         InitOut%WriteOutputUnt( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 14 ) = '  (N/m)  '
         InitOut%WriteOutputHdr( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 15 ) = ' '//trim(chanPrefix)//"Fy"
         InitOut%WriteOutputUnt( (k-1)*p%NumBlNds*AD_numChanPerNode + (j-1)*AD_numChanPerNode + 15 ) = '  (N/m)  '
      end do
   end do
   
   InitOut%Ver = AD_Ver
   
end subroutine AD_SetInitOut
!----------------------------------------------------------------------------------------------------------------------------------   
subroutine AD_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

   type(AD_InitInputType),       intent(in   ) :: InitInp       ! Input data for initialization routine
   type(AD_InputType),           intent(  out) :: u             ! An initial guess for the input; input mesh must be defined
   type(AD_ParameterType),       intent(  out) :: p             ! Parameters
   type(AD_ContinuousStateType), intent(  out) :: x             ! Initial continuous states
   type(AD_DiscreteStateType),   intent(  out) :: xd            ! Initial discrete states
   type(AD_ConstraintStateType), intent(  out) :: z             ! Initial guess of the constraint states
   type(AD_OtherStateType),      intent(  out) :: OtherState    ! Initial other/optimization states
   type(AD_OutputType),          intent(  out) :: y             ! Initial system outputs (outputs are not calculated;
                                                                !   only the output mesh is initialized)
   real(DbKi),                   intent(inout) :: interval      ! Coupling interval in seconds: the rate that
                                                                !   (1) AD_UpdateStates() is called in loose coupling &
                                                                !   (2) AD_UpdateDiscState() is called in tight coupling.
                                                                !   Input is the suggested time from the glue code;
                                                                !   Output is the actual coupling interval that will be used
                                                                !   by the glue code.
   type(AD_InitOutputType),      intent(  out) :: InitOut       ! Output for initialization routine
   integer(IntKi),               intent(  out) :: errStat       ! Error status of the operation
   character(*),                 intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   

      ! Local variables
   integer(IntKi)                              :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                        :: errMsg2       ! temporary error message 
      
   type(AD_InputFile)                          :: InputFileData ! Data stored in the module's input file
   integer(IntKi)                              :: UnEcho        ! Unit number for the echo file
   
   character(*), parameter                     :: RoutineName = 'AD_Init'
   
   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   UnEcho  = -1

      ! Initialize the NWTC Subroutine Library

   call NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information

   call DispNVD( AD_Ver )
   
   
   p%NumBlades = InitInp%NumBlades ! need this before reading the AD input file so that we know how many blade files to read
   p%RootName  = TRIM(InitInp%RootName)//'.AD'
   
      ! Read the primary AeroDyn input file
   call AD_ReadInput( InitInp%InputFile, InputFileData, interval, p%RootName, p%NumBlades, UnEcho, ErrStat2, ErrMsg2 )   
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
         
      
      ! Validate the inputs
   call ValidateInputData( InputFileData, p%NumBlades, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
      
      !............................................................................................
      ! Define parameters
      !............................................................................................
      
      ! Initialize AFI module
   call Init_AFIparams( InputFileData, p%AFI, UnEcho, p%NumBlades, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
         
      
      ! set the rest of the parameters
   call SetParameters( InitInp, InputFileData, p, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
   
      !............................................................................................
      ! Define and initialize inputs here 
      !............................................................................................
   
   call Init_u( u, p, InputFileData, InitInp, errStat2, errMsg2 ) 
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if

      ! 

      !............................................................................................
      ! Initialize the BEMT module (also sets other variables for sub module)
      !............................................................................................
      
   if (p%WakeMod == WakeMod_BEMT) then
      ! initialize BEMT after setting parameters and inputs because we are going to use the already-
      ! calculated node positions from the input meshes
      
      call Init_BEMTmodule( InputFileData, u, OtherState%BEMT_u, p, x%BEMT, xd%BEMT, z%BEMT, &
                            OtherState%BEMT, OtherState%BEMT_y, ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
         if (ErrStat >= AbortErrLev) then
            call Cleanup()
            return
         end if
         
   end if      
      
      
      !............................................................................................
      ! Define outputs here
      !............................................................................................
   call Init_y(y, u, p, errStat2, errMsg2) ! do this after input meshes have been initialized
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
   
   
      !............................................................................................
      ! Initialize states
      !............................................................................................
      
      ! many states are in the BEMT module, which were initialized in BEMT_Init()
      
   call Init_OtherStates(OtherState, p, errStat2, errMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      
      !............................................................................................
      ! Define initialization output here
      !............................................................................................
   call AD_SetInitOut(p, InitOut, errStat2, errMsg2)
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
   
            
   call Cleanup() 
      
contains
   subroutine Cleanup()

      CALL AD_DestroyInputFile( InputFileData, ErrStat2, ErrMsg2 )
      IF ( UnEcho > 0 ) CLOSE( UnEcho )
      
   end subroutine Cleanup

end subroutine AD_Init
!----------------------------------------------------------------------------------------------------------------------------------   
subroutine Init_OtherStates(OtherState, p, errStat, errMsg)
   type(AD_OtherStateType),       intent(inout)  :: OtherState       ! Otherstate data (not defined in submodules)
   type(AD_ParameterType),        intent(in   )  :: p                ! Parameters
   integer(IntKi),                intent(inout)  :: errStat          ! Error status of the operation
   character(*),                  intent(inout)  :: errMsg           ! Error message if ErrStat /= ErrID_None


      ! Local variables
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_OtherStates'

      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   
   call AllocAry( OtherState%DisturbedInflow, 3_IntKi, p%NumBlNds, p%numBlades, 'OtherState%DisturbedInflow', ErrStat2, ErrMsg2 ) ! must be same size as u%InflowOnBlade
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   call AllocAry( OtherState%WithoutSweepPitchTwist, 3_IntKi, 3_IntKi, p%NumBlNds, p%numBlades, 'OtherState%WithoutSweepPitchTwist', ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      
   if (ErrStat >= AbortErrLev) RETURN
      
   
end subroutine Init_OtherStates
!----------------------------------------------------------------------------------------------------------------------------------   
subroutine Init_y(y, u, p, errStat, errMsg)
   type(AD_OutputType),           intent(  out)  :: y               ! Module outputs
   type(AD_InputType),            intent(inout)  :: u               ! Module inputs -- intent(out) because of mesh sibling copy
   type(AD_ParameterType),        intent(in   )  :: p               ! Parameters
   integer(IntKi),                intent(inout)  :: errStat         ! Error status of the operation
   character(*),                  intent(inout)  :: errMsg          ! Error message if ErrStat /= ErrID_None


      ! Local variables
   integer(intKi)                               :: k                 ! loop counter for blades
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_y'

      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   
         
   if (p%TwrAero) then
            
      call MeshCopy ( SrcMesh  = u%TowerMotion    &
                    , DestMesh = y%TowerLoad      &
                    , CtrlCode = MESH_SIBLING     &
                    , IOS      = COMPONENT_OUTPUT &
                    , force    = .TRUE.           &
                    , moment   = .TRUE.           &
                    , ErrStat  = ErrStat2         &
                    , ErrMess  = ErrMsg2          )
   
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
         if (ErrStat >= AbortErrLev) RETURN         
         
         !y%TowerLoad%force = 0.0_ReKi  ! shouldn't have to initialize this
         !y%TowerLoad%moment= 0.0_ReKi  ! shouldn't have to initialize this
         
   end if

   
   allocate( y%BladeLoad(p%numBlades), stat=ErrStat2 )
   if (errStat2 /= 0) then
      call SetErrStat( ErrID_Fatal, 'Error allocating y%BladeLoad.', ErrStat, ErrMsg, RoutineName )      
      return
   end if
   

   do k = 1, p%numBlades
   
      call MeshCopy ( SrcMesh  = u%BladeMotion(k) &
                    , DestMesh = y%BladeLoad(k)   &
                    , CtrlCode = MESH_SIBLING     &
                    , IOS      = COMPONENT_OUTPUT &
                    , force    = .TRUE.           &
                    , moment   = .TRUE.           &
                    , ErrStat  = ErrStat2         &
                    , ErrMess  = ErrMsg2          )
   
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
                           
   end do
   
   call AllocAry( y%WriteOutput, p%numOuts, 'WriteOutput', errStat2, errMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   if (ErrStat >= AbortErrLev) RETURN      
   y%WriteOutput = 0.0_ReKi
   
   
   
end subroutine Init_y
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Init_u( u, p, InputFileData, InitInp, errStat, errMsg )
! This routine is called from AD_Init.
!  it allocates and initializes the inputs to AeroDyn
!..................................................................................................................................

   type(AD_InputType),           intent(  out)  :: u                 ! Input data
   type(AD_ParameterType),       intent(in   )  :: p                 ! Parameters
   type(AD_InputFile),           intent(in   )  :: InputFileData     ! Data stored in the module's input file
   type(AD_InitInputType),       intent(in   )  :: InitInp           ! Input data for AD initialization routine
   integer(IntKi),               intent(inout)  :: errStat           ! Error status of the operation
   character(*),                 intent(inout)  :: errMsg            ! Error message if ErrStat /= ErrID_None


      ! Local variables
   real(reKi)                                   :: position(3)       ! node reference position
   real(reKi)                                   :: positionL(3)      ! node local position
   real(reKi)                                   :: theta(3)          ! Euler angles
   real(reKi)                                   :: orientation(3,3)  ! node reference orientation
   real(reKi)                                   :: orientationL(3,3) ! node local orientation
   
   integer(intKi)                               :: j                 ! counter for nodes
   integer(intKi)                               :: k                 ! counter for blades
   
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'Init_u'

      ! Initialize variables for this routine

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Arrays for InflowWind inputs:
   
   call AllocAry( u%InflowOnBlade, 3_IntKi, p%NumBlNds, p%numBlades, 'u%InflowOnBlade', ErrStat2, ErrMsg2 )
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
   call AllocAry( u%InflowOnTower, 3_IntKi, p%NumTwrNds, 'u%InflowOnTower', ErrStat2, ErrMsg2 ) ! could be size zero
      call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
                
   if (errStat >= AbortErrLev) return      
      
   u%InflowOnBlade = 0.0_ReKi
   
      ! Meshes for motion inputs (ElastoDyn and/or BeamDyn)
         !................
         ! tower
         !................
   if (p%NumTwrNds > 0) then
      
      u%InflowOnTower = 0.0_ReKi 
      
      call MeshCreate ( BlankMesh = u%TowerMotion   &
                       ,IOS       = COMPONENT_INPUT &
                       ,Nnodes    = p%NumTwrNds     &
                       ,ErrStat   = ErrStat2        &
                       ,ErrMess   = ErrMsg2         &
                       ,Orientation     = .true.    &
                       ,TranslationDisp = .true.    &
                       ,TranslationVel  = .true.    &
                       ,RotationVel     = .true.    &
                      )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

      if (errStat >= AbortErrLev) return
            
         ! set node initial position/orientation
      position = 0.0_ReKi
      do j=1,p%NumTwrNds         
         position(3) = InputFileData%TwrElev(j)
         
         call MeshPositionNode(u%TowerMotion, j, position, errStat2, errMsg2)  ! orientation is identity by default
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      end do !j
         
         ! create line2 elements
      do j=1,p%NumTwrNds-1
         call MeshConstructElement( u%TowerMotion, ELEMENT_LINE2, errStat2, errMsg2, p1=j, p2=j+1 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
      end do !j
            
      call MeshCommit(u%TowerMotion, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
      if (errStat >= AbortErrLev) return

      
      u%TowerMotion%Orientation     = u%TowerMotion%RefOrientation
      u%TowerMotion%TranslationDisp = 0.0_ReKi
      u%TowerMotion%TranslationVel  = 0.0_ReKi
      u%TowerMotion%RotationVel     = 0.0_ReKi
      
   end if ! we compute tower loads
   
         !................
         ! hub
         !................
   
      call MeshCreate ( BlankMesh = u%HubMotion     &
                       ,IOS       = COMPONENT_INPUT &
                       ,Nnodes    = 1               &
                       ,ErrStat   = ErrStat2        &
                       ,ErrMess   = ErrMsg2         &
                       ,Orientation     = .true.    &
                       ,TranslationDisp = .true.    &
                       ,TranslationVel  = .true.    &
                       ,RotationVel     = .true.    &
                      )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

      if (errStat >= AbortErrLev) return
                     
      call MeshPositionNode(u%HubMotion, 1, InitInp%HubPosition, errStat2, errMsg2, InitInp%HubOrientation)
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         
      call MeshConstructElement( u%HubMotion, ELEMENT_POINT, errStat2, errMsg2, p1=1 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
      call MeshCommit(u%HubMotion, errStat2, errMsg2 )
         call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
      if (errStat >= AbortErrLev) return

         
      u%HubMotion%Orientation     = u%HubMotion%RefOrientation
      u%HubMotion%TranslationDisp = 0.0_ReKi
      u%HubMotion%TranslationVel  = 0.0_ReKi
      u%HubMotion%RotationVel     = 0.0_ReKi   
      
   
         !................
         ! blade roots
         !................
         
      allocate( u%BladeRootMotion(p%NumBlades), STAT = ErrStat2 )
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Fatal, 'Error allocating u%BladeRootMotion array.', ErrStat, ErrMsg, RoutineName )
         return
      end if      
      
      do k=1,p%NumBlades
         call MeshCreate ( BlankMesh = u%BladeRootMotion(k)                  &
                          ,IOS       = COMPONENT_INPUT                       &
                          ,Nnodes    = 1                                     &
                          ,ErrStat   = ErrStat2                              &
                          ,ErrMess   = ErrMsg2                               &
                          ,Orientation     = .true.                          &
                          ,TranslationDisp = .true.                          &
                          ,TranslationVel  = .true.                          &
                          ,RotationVel     = .true.                          &
                         )
               call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

         if (errStat >= AbortErrLev) return
            
         call MeshPositionNode(u%BladeRootMotion(k), 1, InitInp%BladeRootPosition(:,k), errStat2, errMsg2, InitInp%BladeRootOrientation(:,:,k))
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
                     
         call MeshConstructElement( u%BladeRootMotion(k), ELEMENT_POINT, errStat2, errMsg2, p1=1 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
         call MeshCommit(u%BladeRootMotion(k), errStat2, errMsg2 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
         if (errStat >= AbortErrLev) return

      
         u%BladeRootMotion(k)%Orientation     = u%BladeRootMotion(k)%RefOrientation
         u%BladeRootMotion(k)%TranslationDisp = 0.0_ReKi
         u%BladeRootMotion(k)%TranslationVel  = 0.0_ReKi
         u%BladeRootMotion(k)%RotationVel     = 0.0_ReKi
   
   end do !k=numBlades      
      
      
         !................
         ! blades
         !................
   
      allocate( u%BladeMotion(p%NumBlades), STAT = ErrStat2 )
      if (ErrStat2 /= 0) then
         call SetErrStat( ErrID_Fatal, 'Error allocating u%BladeMotion array.', ErrStat, ErrMsg, RoutineName )
         return
      end if
      
      do k=1,p%NumBlades
         call MeshCreate ( BlankMesh = u%BladeMotion(k)                     &
                          ,IOS       = COMPONENT_INPUT                      &
                          ,Nnodes    = InputFileData%BladeProps(k)%NumBlNds &
                          ,ErrStat   = ErrStat2                             &
                          ,ErrMess   = ErrMsg2                              &
                          ,Orientation     = .true.                         &
                          ,TranslationDisp = .true.                         &
                          ,TranslationVel  = .true.                         &
                          ,RotationVel     = .true.                         &
                         )
               call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )

         if (errStat >= AbortErrLev) return
            
                        
         do j=1,InputFileData%BladeProps(k)%NumBlNds

               ! reference position of the jth node in the kth blade, relative to the root in the local blade coordinate system:
            positionL(1) = InputFileData%BladeProps(k)%BlCrvAC(j)
            positionL(2) = InputFileData%BladeProps(k)%BlSwpAC(j)
            positionL(3) = InputFileData%BladeProps(k)%BlSpn(  j)
            
               ! reference orientation of the jth node in the kth blade, relative to the root in the local blade coordinate system:
            theta(1)     =  0.0_ReKi
            theta(2)     =  InputFileData%BladeProps(k)%BlCrvAng(j)
            theta(3)     = -InputFileData%BladeProps(k)%BlTwist( j)            
            orientationL = EulerConstruct( theta )
                     
               ! reference position of the jth node in the kth blade:
            position = InitInp%BladeRootPosition(:,k) + matmul(positionL,u%BladeRootMotion(k)%RefOrientation(:,:,1))  ! note that because positionL is a 1-D array, we're doing the transpose of matmul(transpose(u%BladeRootMotion(k)%RefOrientation),positionL)
            
               ! reference orientation of the jth node in the kth blade
            orientation = matmul( orientationL, u%BladeRootMotion(k)%RefOrientation(:,:,1) )

            call MeshPositionNode(u%BladeMotion(k), j, position, errStat2, errMsg2, orientation)
               call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
               
         end do ! j=blade nodes
         
            ! create line2 elements
         do j=1,InputFileData%BladeProps(k)%NumBlNds-1
            call MeshConstructElement( u%BladeMotion(k), ELEMENT_LINE2, errStat2, errMsg2, p1=j, p2=j+1 )
               call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
         end do !j
            
         call MeshCommit(u%BladeMotion(k), errStat2, errMsg2 )
            call SetErrStat( errStat2, errMsg2, errStat, errMsg, RoutineName )
            
         if (errStat >= AbortErrLev) return

      
         u%BladeMotion(k)%Orientation     = u%BladeMotion(k)%RefOrientation
         u%BladeMotion(k)%TranslationDisp = 0.0_ReKi
         u%BladeMotion(k)%TranslationVel  = 0.0_ReKi
         u%BladeMotion(k)%RotationVel     = 0.0_ReKi
   
   end do !k=numBlades
   
   
end subroutine Init_u
!----------------------------------------------------------------------------------------------------------------------------------
subroutine SetParameters( InitInp, InputFileData, p, ErrStat, ErrMsg )
! This routine is called from AD_Init.
! The parameters are set here and not changed during the simulation.
!..................................................................................................................................
   TYPE(AD_InitInputType),       intent(in   )  :: InitInp          ! Input data for initialization routine, out is needed because of copy below
   TYPE(AD_InputFile),           INTENT(INout)  :: InputFileData    ! Data stored in the module's input file -- intent(out) only for move_alloc statements
   TYPE(AD_ParameterType),       INTENT(INOUT)  :: p                ! Parameters
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat          ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg           ! Error message if ErrStat /= ErrID_None


      ! Local variables
   !CHARACTER(ErrMsgLen)                          :: ErrMsg2         ! temporary Error message if ErrStat /= ErrID_None
   !INTEGER(IntKi)                                :: ErrStat2        ! temporary Error status of the operation
   !INTEGER(IntKi)                                :: i, j
   character(*), parameter                       :: RoutineName = 'SetParameters'
   
      ! Initialize variables for this routine

   ErrStat  = ErrID_None
   ErrMsg   = ""

   p%DT               = InputFileData%DTAero      
   p%WakeMod          = InputFileData%WakeMod
   p%TwrInflnc        = InputFileData%TwrInflnc
   p%TwrAero          = InputFileData%TwrAero
   
 ! p%numBlades        = InitInp%numBlades    ! this was set earlier because it was necessary
   p%NumBlNds         = InputFileData%BladeProps(1)%NumBlNds
   if (p%TwrInflnc == TwrInflnc_none .and. .not. p%TwrAero) then
      p%NumTwrNds     = 0
   else
      p%NumTwrNds     = InputFileData%NumTwrNds
      
      call move_alloc( InputFileData%TwrDiam, p%TwrDiam )
      call move_alloc( InputFileData%TwrCd,   p%TwrCd )      
   end if
   
   p%AirDens          = InputFileData%AirDens          
   p%KinVisc          = InputFileData%KinVisc
   p%SpdSound         = InputFileData%SpdSound
   
  !p%AFI     ! set in call to AFI_Init() [called early because it wants to use the same echo file as AD]
  !p%BEMT    ! set in call to BEMT_Init()

  !p%numOuts        = p%numBladeNodes*p%numBlades*AD_numChanPerNode
   p%numOuts        = InputFileData%BladeProps(1)%NumBlNds*InitInp%numBlades*AD_numChanPerNode
  !p%RootName       = TRIM(InitInp%RootName)//'.AD'   ! set earlier to it could be used
  !p%OutParam STILL NEEDS TO BE SET!!! FIX ME!!!! 
   
   
      
      
   if (ErrStat >= AbortErrLev) RETURN
      
            
   
end subroutine SetParameters
!----------------------------------------------------------------------------------------------------------------------------------
subroutine AD_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(AD_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(AD_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
      TYPE(AD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(AD_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(AD_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(AD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(AD_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Place any last minute operations or calculations here:


         ! Close files here:



         ! Destroy the input data:

      CALL AD_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:

      CALL AD_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:

      CALL AD_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL AD_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL AD_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL AD_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )


         ! Destroy the output data:

      CALL AD_DestroyOutput( y, ErrStat, ErrMsg )




END SUBROUTINE AD_End
!----------------------------------------------------------------------------------------------------------------------------------
subroutine AD_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, errStat, errMsg )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input Time t; Continuous and discrete states are updated for t + Interval
!..................................................................................................................................

   real(DbKi),                     intent(in   ) :: t          ! Current simulation time in seconds
   integer(IntKi),                 intent(in   ) :: n          ! Current simulation time step n = 0,1,...
   type(AD_InputType),             intent(inout) :: u(:)       ! Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
   real(DbKi),                     intent(in   ) :: utimes(:)  ! Times associated with u(:), in seconds
   type(AD_ParameterType),         intent(in   ) :: p          ! Parameters
   type(AD_ContinuousStateType),   intent(inout) :: x          ! Input: Continuous states at t;
                                                               !   Output: Continuous states at t + Interval
   type(AD_DiscreteStateType),     intent(inout) :: xd         ! Input: Discrete states at t;
                                                               !   Output: Discrete states at t  + Interval
   type(AD_ConstraintStateType),   intent(inout) :: z          ! Input: Initial guess of constraint states at t+dt;
                                                               !   Output: Constraint states at t+dt
   type(AD_OtherStateType),        intent(inout) :: OtherState ! Other/optimization states
   integer(IntKi),                 intent(  out) :: errStat    ! Error status of the operation
   character(*),                   intent(  out) :: errMsg     ! Error message if ErrStat /= ErrID_None

   ! local variables
   type(AD_InputType)                           :: uInterp     ! Interpolated/Extrapolated input
   
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'AD_UpdateStates'
      
   ErrStat = ErrID_None
   ErrMsg  = ""
     

   call AD_CopyInput( u(1), uInterp, MESH_NEWCOPY, errStat2, errMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      if (ErrStat >= AbortErrLev) then
         call cleanup()
         return
      end if
      
   call AD_Input_ExtrapInterp(u,utimes,uInterp,t, errStat2, errMsg2)
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         
   call SetInputs(p, uInterp, OtherState, errStat2, errMsg2)      
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
   !modify for tower influence
      
   
   if ( p%WakeMod == WakeMod_BEMT ) then
                  
         ! Call into the BEMT update states    NOTE:  This is a non-standard framework interface!!!!!  GJH
      call BEMT_UpdateStates(t, n, OtherState%BEMT_u,  p%BEMT, x%BEMT, xd%BEMT, z%BEMT, OtherState%BEMT, p%AFI%AFInfo, errStat2, errMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
   end if
           
   call cleanup()
   
contains
   subroutine cleanup()
      call AD_DestroyInput( uInterp, errStat2, errMsg2)
   end subroutine cleanup
end subroutine AD_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
subroutine AD_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! Routine for computing outputs, used in both loose and tight coupling.
! This SUBROUTINE is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
! NOTE: the descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
! for a complete description of each output parameter.
! NOTE: no matter how many channels are selected for output, all of the outputs are calcalated
! All of the calculated output channels are placed into the OtherState%AllOuts(:), while the channels selected for outputs are
! placed in the y%WriteOutput(:) array.
!..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(AD_InputType),           INTENT(IN   )  :: u           ! Inputs at Time t
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(AD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(AD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(AD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(AD_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                               !   nectivity information does not have to be recalculated)
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


   INTEGER(IntKi)                               :: i,j
   character(10)                                :: chanPrefix
   real(ReKi)                                   :: q
   
   integer(intKi)                               :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2
   character(*), parameter                      :: RoutineName = 'AD_CalcOutput'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""

   
   call SetInputs(p, u, OtherState, errStat2, errMsg2)      
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
!-------------------------------------------------------
!     Start of BEMT-level calculations of outputs
!-------------------------------------------------------
   
   
   if (p%WakeMod == WakeMod_BEMT) then
         
      ! Call the BEMT module CalcOutput.  Notice that the BEMT outputs are purposely attached to AeroDyn's OtherState structure to
      ! avoid issues with the coupling code
      !
      ! NOTE: the x%BEMT, xd%BEMT, and OtherState%BEMT are simply dummy variables because the BEMT module does not use them or return them
      !
   
      call BEMT_CalcOutput(t, OtherState%BEMT_u, p%BEMT, x%BEMT, xd%BEMT, z%BEMT, OtherState%BEMT, p%AFI%AFInfo, OtherState%BEMT_y, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
                  
      call SetOutputsFromBEMT(p, u, OtherState%WithoutSweepPitchTwist, OtherState%BEMT_y, y, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         
         ! TODO Check error status
         
         
            ! Loop over blades and nodes to populate the output channel names and units
   
   do j=1,p%numBlades
      do i=1,p%NumBlNds
         
         
         chanPrefix = "B"//trim(num2lstr(j))//"N"//trim(num2lstr(i))
         y%WriteOutput( (j-1)*p%NumBlNds*AD_numChanPerNode + (i-1)*AD_numChanPerNode +  1 ) = OtherState%BEMT_u%theta(i,j)*R2D      
         y%WriteOutput( (j-1)*p%NumBlNds*AD_numChanPerNode + (i-1)*AD_numChanPerNode +  2 ) = OtherState%BEMT_u%psi(j)*R2D        
         y%WriteOutput( (j-1)*p%NumBlNds*AD_numChanPerNode + (i-1)*AD_numChanPerNode +  3 ) = OtherState%BEMT_u%Vx(i,j)         
         y%WriteOutput( (j-1)*p%NumBlNds*AD_numChanPerNode + (i-1)*AD_numChanPerNode +  4 ) = OtherState%BEMT_u%Vy(i,j)         
         y%WriteOutput( (j-1)*p%NumBlNds*AD_numChanPerNode + (i-1)*AD_numChanPerNode +  5 ) = OtherState%BEMT_y%AxInduction(i,j)         
         y%WriteOutput( (j-1)*p%NumBlNds*AD_numChanPerNode + (i-1)*AD_numChanPerNode +  6 ) = OtherState%BEMT_y%TanInduction(i,j)         
         y%WriteOutput( (j-1)*p%NumBlNds*AD_numChanPerNode + (i-1)*AD_numChanPerNode +  7 ) = OtherState%BEMT_y%inducedVel(i,j)         
         y%WriteOutput( (j-1)*p%NumBlNds*AD_numChanPerNode + (i-1)*AD_numChanPerNode +  8 ) = OtherState%BEMT_y%phi(i,j)*R2D         
         y%WriteOutput( (j-1)*p%NumBlNds*AD_numChanPerNode + (i-1)*AD_numChanPerNode +  9 ) = OtherState%BEMT_y%AOA(i,j)*R2D         
         y%WriteOutput( (j-1)*p%NumBlNds*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 10 ) = OtherState%BEMT_y%Cl(i,j)         
         y%WriteOutput( (j-1)*p%NumBlNds*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 11 ) = OtherState%BEMT_y%Cd(i,j)         
         y%WriteOutput( (j-1)*p%NumBlNds*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 12 ) = OtherState%BEMT_y%Cx(i,j)         
         y%WriteOutput( (j-1)*p%NumBlNds*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 13 ) = OtherState%BEMT_y%Cy(i,j)
         
         q = 0.5*p%AirDens*p%BEMT%chord(i,j)*OtherState%BEMT_y%inducedVel(i,j)**2
         
         y%WriteOutput( (j-1)*p%NumBlNds*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 14 ) =  q*OtherState%BEMT_y%Cx(i,j)
         y%WriteOutput( (j-1)*p%NumBlNds*AD_numChanPerNode + (i-1)*AD_numChanPerNode + 15 ) = -q*OtherState%BEMT_y%Cy(i,j)
         
      end do
   end do
!   
!-------------------------------------------------------         
         
   end if
   
!-------------------------------------------------------   
!     End of BEMT calculations  
!-------------------------------------------------------
   
!-------------------------------------------------------
!     Start AeroDyn-level calculations of outputs
!-------------------------------------------------------
   
   if ( p%TwrAero ) then
      call ADTwr_CalcOutput(p, u, OtherState, y, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
   end if
         
!-------------------------------------------------------   
!     End of AeroDyn calculations  
!-------------------------------------------------------   
                      
   
   
end subroutine AD_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
subroutine AD_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMsg )
! Tight coupling routine for solving for the residual of the constraint state equations
!..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )   :: Time        ! Current simulation time in seconds
   TYPE(AD_InputType),           INTENT(IN   )   :: u           ! Inputs at Time
   TYPE(AD_ParameterType),       INTENT(IN   )   :: p           ! Parameters
   TYPE(AD_ContinuousStateType), INTENT(IN   )   :: x           ! Continuous states at Time
   TYPE(AD_DiscreteStateType),   INTENT(IN   )   :: xd          ! Discrete states at Time
   TYPE(AD_ConstraintStateType), INTENT(INOUT)   :: z           ! Constraint states at Time (possibly a guess)
   TYPE(AD_OtherStateType),      INTENT(INOUT)   :: OtherState  ! Other/optimization states
   TYPE(AD_ConstraintStateType), INTENT(  OUT)   :: z_residual  ! Residual of the constraint state equations using
                                                                !     the input values described above
   INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                  INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


   
      ! Local variables   
   integer(intKi)                               :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2
   character(*), parameter                      :: RoutineName = 'AD_CalcConstrStateResidual'
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   call SetInputs(p, u, OtherState, errStat2, errMsg2)      
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
   
   if (p%WakeMod == WakeMod_BEMT) then
      
   
      call BEMT_CalcConstrStateResidual( Time, OtherState%BEMT_u, p%BEMT, x%BEMT, xd%BEMT, z%BEMT, OtherState%BEMT, &
                                         z_residual%BEMT, p%AFI%AFInfo, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         
   end if
   
   
END SUBROUTINE AD_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
subroutine SetInputs(p, u, OtherState, errStat, errMsg)
   type(AD_ParameterType),       intent(in   )  :: p                      ! AD parameters
   type(AD_InputType),           intent(in   )  :: u                      ! AD Inputs at Time
   type(AD_OtherStateType),      intent(inout)  :: OtherState             ! OtherStates (with inputs to submodules)
   integer(IntKi),               intent(  out)  :: ErrStat                ! Error status of the operation
   character(*),                 intent(  out)  :: ErrMsg                 ! Error message if ErrStat /= ErrID_None
                                 
   ! local variables             
   !integer(intKi)                               :: j                      ! loop counter for nodes
   !integer(intKi)                               :: k                      ! loop counter for blades
   integer(intKi)                               :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2
   character(*), parameter                      :: RoutineName = 'SetInputs'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   OtherState%DisturbedInflow = u%InflowOnBlade

   
   
   if ( p%WakeMod == WakeMod_BEMT ) then
      
         ! This needs to extract the inputs from the AD data types (mesh) and massage them for the BEMT module
      call SetInputsForBEMT(p, u, OtherState%BEMT_u, OtherState%DisturbedInflow, OtherState%WithoutSweepPitchTwist, errStat2, errMsg2)  
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   end if
   
end subroutine SetInputs
!----------------------------------------------------------------------------------------------------------------------------------
subroutine SetInputsForBEMT(p, u, BEMT_u, DisturbedInflow, WithoutSweepPitchTwist, errStat, errMsg)

   type(AD_ParameterType),  intent(in   )  :: p                               ! AD parameters
   type(AD_InputType),      intent(in   )  :: u                               ! AD Inputs at Time
   type(BEMT_InputType),    intent(inout)  :: BEMT_u                          ! BEMT Inputs at Time
   real(ReKi),              intent(in   )  :: DisturbedInflow(:,:,:)          ! inflow wind velocity disturbed by tower influence 
   real(ReKi),              intent(  out)  :: WithoutSweepPitchTwist(:,:,:,:) ! modified orientation matrix
   integer(IntKi),          intent(  out)  :: ErrStat                         ! Error status of the operation
   character(*),            intent(  out)  :: ErrMsg                          ! Error message if ErrStat /= ErrID_None
      
   ! local variables
   real(ReKi)                              :: V_diskAvg(3)
   real(ReKi)                              :: x_hat(3)
   real(ReKi)                              :: y_hat(3)
   real(ReKi)                              :: z_hat(3)
   real(ReKi)                              :: x_hat_disk(3)
   real(ReKi)                              :: y_hat_disk(3)
   real(ReKi)                              :: z_hat_disk(3)
   real(ReKi)                              :: tmp(3)
   real(ReKi)                              :: theta(3)
   real(ReKi)                              :: orientation(3,3)
   real(ReKi)                              :: orientation_nopitch(3,3)
   real(ReKi)                              :: tmp_sz, tmp_sz_y
   
   integer(intKi)                          :: j                      ! loop counter for nodes
   integer(intKi)                          :: k                      ! loop counter for blades
   integer(intKi)                          :: ErrStat2
   character(ErrMsgLen)                    :: ErrMsg2
   character(*), parameter                 :: RoutineName = 'SetInputsForBEMT'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! calculate disk-averaged relative wind speed, V_DiskAvg
   V_diskAvg = 0.0_ReKi
   do k=1,p%NumBlades
      do j=1,p%NumBlNds
         tmp = DisturbedInflow(:,j,k) - u%BladeMotion(k)%TranslationVel(:,j)
         V_diskAvg = V_diskAvg + tmp         
      end do
   end do
   
      ! orientation vectors:
   x_hat_disk = u%HubMotion%Orientation(1,:,1) !actually also x_hat_hub      
   
  tmp    = dot_product( V_diskAvg, x_hat_disk ) * x_hat_disk - V_diskAvg
  tmp_sz = TwoNorm(tmp)
  if ( EqualRealNos( tmp_sz, 0.0_ReKi ) ) then
     y_hat_disk = u%HubMotion%Orientation(2,:,1)
     z_hat_disk = u%HubMotion%Orientation(3,:,1)
  else
     y_hat_disk = tmp / tmp_sz
     z_hat_disk = cross_product( V_diskAvg, x_hat_disk ) / tmp_sz
  end if
     
      ! "Angular velocity of rotor" rad/s
   BEMT_u%omega   = dot_product( u%HubMotion%RotationVel(:,1), x_hat_disk )    
   
      ! "Angle between the vector normal to the rotor plane and the wind vector (e.g., the yaw angle in the case of no tilt)" rad 
   tmp_sz = TwoNorm( V_diskAvg )
   if ( EqualRealNos( tmp_sz, 0.0_ReKi ) ) then
      BEMT_u%chi0 = 0.0_ReKi
   else
      BEMT_u%chi0 = acos( dot_product(V_diskAvg,x_hat_disk) / tmp_sz)      
   end if
   
      ! "Azimuth angle" rad
   do k=1,p%NumBlades
      z_hat = u%BladeRootMotion(k)%Orientation(3,:,1)      
      tmp_sz_y = -1.0*dot_product(z_hat,y_hat_disk)
      tmp_sz   =      dot_product(z_hat,z_hat_disk)
      if ( EqualRealNos(tmp_sz_y,0.0_ReKi) .and. EqualRealNos(tmp_sz,0.0_ReKi) ) then
         BEMT_u%psi(k) = 0.0_ReKi
      else
         BEMT_u%psi(k) = atan2( tmp_sz_y, tmp_sz )
      end if      
   end do
   
      ! theta, "Twist angle (includes all sources of twist)" rad
      ! Vx, "Local axial velocity at node" m/s
      ! Vy, "Local tangential velocity at node" m/s
   do k=1,p%NumBlades
      
         ! construct system equivalent to u%BladeRootMotion(k)%Orientation, but without the blade-pitch angle:
      
      ! orientation = matmul( u%BladeRootMotion(k)%Orientation, transpose(u%HubMotion%Orientation(:,:,1)) )
      call LAPACK_gemm( 'n', 't', 1.0_ReKi, u%BladeRootMotion(k)%Orientation(:,:,1), u%HubMotion%Orientation(:,:,1), 0.0_ReKi, orientation, errStat2, errMsg2)
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      theta = EulerExtract( orientation ) !hub_theta_root(k)
      theta(3) = 0.0_ReKi         
      orientation_nopitch = EulerConstruct( theta ) ! withoutPitch_theta_Root(k)
      
      do j=1,p%NumBlNds         
         
            ! form coordinate system equivalent to u%BladeMotion(k)%Orientation(:,:,j) but iwthout live sweep (due to in-plane
            ! deflection), blade-pitch and twist (aerodynamic + elastic) angles:
         
         ! orientation = matmul( u%BladeMotion(k)%Orientation(:,:,j), transpose(orientation_nopitch) )
         call LAPACK_gemm( 'n', 't', 1.0_ReKi, u%BladeMotion(k)%Orientation(:,:,j), orientation_nopitch, 0.0_ReKi, orientation, errStat2, errMsg2)
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         theta = EulerExtract( orientation ) !root(k)WithoutPitch_theta(j)_blade(k)
         
         BEMT_u%theta(j,k) = -theta(3) ! local pitch + twist (aerodyanmic + elastic) angle of the jth node in the kth blade
         
         theta(1) = 0.0_ReKi
         theta(3) = 0.0_ReKi
         WithoutSweepPitchTwist(:,:,j,k) = matmul( EulerConstruct( theta ), orientation_nopitch ) ! WithoutSweepPitch+Twist_theta(j)_Blade(k)
                           
         x_hat = WithoutSweepPitchTwist(1,:,j,k)
         y_hat = WithoutSweepPitchTwist(2,:,j,k)
         tmp   = DisturbedInflow(:,j,k) - u%BladeMotion(k)%TranslationVel(:,j)
         
         BEMT_u%Vx(j,k) = dot_product( tmp, x_hat ) ! normal component (normal to the plane, not chord) of the inflow velocity of the jth node in the kth blade
         BEMT_u%Vy(j,k) = dot_product( tmp, y_hat ) ! tangential component (tangential to the plane, not chord) of the inflow velocity of the jth node in the kth blade
         
      end do !j=nodes
   end do !k=blades
   
   
      ! "Radial distance from center-of-rotation to node" m
   do k=1,p%NumBlades
      do j=1,p%NumBlNds
         
            ! displaced position of the jth node in the kth blade relative to the hub:
         tmp =  u%BladeMotion(k)%Position(:,j) + u%BladeMotion(k)%TranslationDisp(:,j) &
              - u%HubMotion%Position(:,1)      - u%HubMotion%TranslationDisp(:,1)
         
            ! local radius (normalized distance from rotor centerline)
         tmp_sz_y = dot_product( tmp, y_hat_disk )**2
         tmp_sz   = dot_product( tmp, z_hat_disk )**2
         BEMT_u%rLocal(j,k) = sqrt( tmp_sz + tmp_sz_y )
         
      end do !j=nodes
      BEMT_u%rTip(k) =  BEMT_u%rLocal(p%NumBlNds,k)  ! "Radial distance from center-of-rotation to blade tip" m
   end do !k=blades
         
   ! values for coupled model:
! FIX ME!!!!
 !???  ! "Local upstream velocity at node" m/s
   do k=1,p%NumBlades
      do j=1,p%NumBlNds
         BEMT_u%Vinf(j,k) = TwoNorm( DisturbedInflow(:,j,k) ) 
      end do
   end do
   
    !BEMT_u%lambda  =                        !-        ! "Average tip speed ratio for rotor disk" -
   
   
end subroutine SetInputsForBEMT
!----------------------------------------------------------------------------------------------------------------------------------
subroutine SetOutputsFromBEMT(p, u, WithoutSweepPitchTwist, BEMT_y, y, ErrStat, ErrMsg )

   type(AD_ParameterType),  intent(in   )  :: p                               ! AD parameters
   type(AD_OutputType),     intent(inout)  :: y                               ! AD outputs 
   type(AD_InputType),      intent(in   )  :: u                               ! AD inputs 
   type(BEMT_OutputType),   intent(in   )  :: BEMT_y                          ! BEMT outputs
   real(ReKi),              intent(in   )  :: WithoutSweepPitchTwist(:,:,:,:) ! modified orientation matrix
   integer(IntKi),          intent(  out)  :: ErrStat                         ! Error status of the operation
   character(*),            intent(  out)  :: ErrMsg                          ! Error message if ErrStat /= ErrID_None

   integer(intKi)                          :: j                      ! loop counter for nodes
   integer(intKi)                          :: k                      ! loop counter for blades
   real(reKi)                              :: orientation(3,3)
   real(reki)                              :: force(3)
   real(reki)                              :: moment(3)
   real(reki)                              :: q
   
   integer(intKi)                          :: ErrStat2
   character(ErrMsgLen)                    :: ErrMsg2
   character(*), parameter                 :: RoutineName = 'SetOutputsFromBEMT'
   
   
   force(3)    =  0.0_ReKi          
   moment(1:2) =  0.0_ReKi          
   do k=1,p%NumBlades
      do j=1,p%NumBlNds
         
         ! orientation = matmul( u%BladeMotion(k)%Orientation(:,:,j), transpose(WithoutSweepPitchTwist(:,:,j,k)) )
         call LAPACK_gemm( 'n', 't', 1.0_ReKi, u%BladeMotion(k)%Orientation(:,:,j), WithoutSweepPitchTwist(:,:,j,k), 0.0_ReKi, orientation, errStat2, errMsg2)
            call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
             
         q = 0.5 * p%airDens * BEMT_y%inducedVel(j,k)**2    ! dynamic pressure of the jth node in the kth blade
         force(1) =  BEMT_y%cx(j,k) * p%BEMT%chord(j,k)     ! X = normal force per unit length (normal to the plane, not chord) of the jth node in the kth blade
         force(2) = -BEMT_y%cy(j,k) * p%BEMT%chord(j,k)     ! Y = tangential force per unit length (tangential to the plane, not chord) of the jth node in the kth blade
         moment(3)= -BEMT_y%cm(j,k) * p%BEMT%chord(j,k)**2  ! M = pitching moment per unit length of the jth node in the kth blade
         
         y%BladeLoad(k)%Force(:,j)  = matmul( orientation, force )   ! force per unit length of the jth node in the kth blade
         y%BladeLoad(k)%Moment(:,j) = matmul( orientation, moment )  ! moment per unit length of the jth node in the kth blade
         
      end do !j=nodes
   end do !k=blades
   
   
end subroutine SetOutputsFromBEMT
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD_ReadInput( InputFileName, InputFileData, Default_DT, OutFileRoot, NumBlades, UnEcho, ErrStat, ErrMsg )
! This subroutine reads the input file and stores all the data in the AD_InputFile structure.
! It does not perform data validation.
!..................................................................................................................................

      ! Passed variables
   REAL(DbKi),              INTENT(IN)    :: Default_DT      ! The default DT (from glue code)

   CHARACTER(*),            INTENT(IN)    :: InputFileName   ! Name of the input file
   CHARACTER(*),            INTENT(IN)    :: OutFileRoot     ! The rootname of all the output files written by this routine.

   TYPE(AD_InputFile),      INTENT(OUT)   :: InputFileData   ! Data stored in the module's input file
   INTEGER(IntKi),          INTENT(OUT)   :: UnEcho          ! Unit number for the echo file

   INTEGER(IntKi),          INTENT(IN)    :: NumBlades       ! Number of blades for this model
   INTEGER(IntKi),          INTENT(OUT)   :: ErrStat         ! The error status code
   CHARACTER(*),            INTENT(OUT)   :: ErrMsg          ! The error message, if an error occurred

      ! local variables

   INTEGER(IntKi)                         :: I
   INTEGER(IntKi)                         :: ErrStat2        ! The error status code
   CHARACTER(ErrMsgLen)                   :: ErrMsg2         ! The error message, if an error occurred

   CHARACTER(1024)                        :: ADBlFile(MaxBl) ! File that contains the blade information (specified in the primary input file)
   CHARACTER(*), PARAMETER                :: RoutineName = 'AD_ReadInput'
   
   
      ! initialize values:

   ErrStat = ErrID_None
   ErrMsg  = ''
   UnEcho  = -1
   InputFileData%DTAero = Default_DT  ! the glue code's suggested DT for the module (may be overwritten in ReadPrimaryFile())

      ! get the primary/platform input-file data
      ! sets UnEcho, ADBlFile
   
   CALL ReadPrimaryFile( InputFileName, InputFileData, ADBlFile, OutFileRoot, UnEcho, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
      

      ! get the blade input-file data
      
   ALLOCATE( InputFileData%BladeProps( NumBlades ), STAT = ErrStat2 )
   IF (ErrStat2 /= 0) THEN
      CALL SetErrStat(ErrID_Fatal,"Error allocating memory for BladeProps.", ErrStat, ErrMsg, RoutineName)
      CALL Cleanup()
      RETURN
   END IF
      
   DO I=1,NumBlades
      CALL ReadBladeInputs ( ADBlFile(I), InputFileData%BladeProps(I), UnEcho, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName//TRIM(':Blade')//TRIM(Num2LStr(I)))
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
   END DO
   
      

   CALL Cleanup ( )


CONTAINS
   !...............................................................................................................................
   SUBROUTINE Cleanup()
   ! This subroutine cleans up before exiting this subroutine
   !...............................................................................................................................

      ! IF ( UnEcho > 0 ) CLOSE( UnEcho )

   END SUBROUTINE Cleanup

END SUBROUTINE AD_ReadInput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadPrimaryFile( InputFile, InputFileData, ADBlFile, OutFileRoot, UnEc, ErrStat, ErrMsg )
! This routine reads in the primary AeroDyn input file and places the values it reads in the InputFileData structure.
!   It opens and prints to an echo file if requested.
!..................................................................................................................................


   implicit                        none

      ! Passed variables
   integer(IntKi),     intent(out)     :: UnEc                                ! I/O unit for echo file. If > 0, file is open for writing.
   integer(IntKi),     intent(out)     :: ErrStat                             ! Error status

   character(*),       intent(out)     :: ADBlFile(MaxBl)                     ! name of the files containing blade inputs
   character(*),       intent(in)      :: InputFile                           ! Name of the file containing the primary input data
   character(*),       intent(out)     :: ErrMsg                              ! Error message
   character(*),       intent(in)      :: OutFileRoot                         ! The rootname of the echo file, possibly opened in this routine

   type(AD_InputFile), intent(inout)   :: InputFileData                       ! All the data in the AeroDyn input file
   
      ! Local variables:
   real(ReKi)                    :: TmpAry(3)                                 ! array to help read tower properties table
   integer(IntKi)                :: I                                         ! loop counter
   integer(IntKi)                :: UnIn                                      ! Unit number for reading file
     
   integer(IntKi)                :: ErrStat2, IOS                             ! Temporary Error status
   logical                       :: Echo                                      ! Determines if an echo file should be written
   character(ErrMsgLen)          :: ErrMsg2                                   ! Temporary Error message
   character(1024)               :: PriPath                                   ! Path name of the primary file
   character(1024)               :: FTitle                                    ! "File Title": the 2nd line of the input file, which contains a description of its contents
   character(200)                :: Line                                      ! Temporary storage of a line from the input file (to compare with "default")
   character(*), parameter       :: RoutineName = 'ReadPrimaryFile'
   
   
      ! Initialize some variables:
   ErrStat = ErrID_None
   ErrMsg  = ""
      
   UnEc = -1
   Echo = .FALSE.   
   CALL GetPath( InputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.
   

   CALL AllocAry( InputFileData%OutList, MaxOutPts, "Outlist", ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
        
      ! Get an available unit number for the file.

   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! Open the Primary input file.

   CALL OpenFInpFile ( UnIn, InputFile, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      END IF
      
                  
      
   ! Read the lines up/including to the "Echo" simulation control variable
   ! If echo is FALSE, don't write these lines to the echo file. 
   ! If Echo is TRUE, rewind and write on the second try.
   
   I = 1 !set the number of times we've read the file
   DO 
   !----------- HEADER -------------------------------------------------------------
   
      CALL ReadCom( UnIn, InputFile, 'File header: Module Version (line 1)', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
      CALL ReadStr( UnIn, InputFile, FTitle, 'FTitle', 'File Header: File Description (line 2)', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
   
   
   !----------- GENERAL OPTIONS ----------------------------------------------------
   
      CALL ReadCom( UnIn, InputFile, 'Section Header: General Options', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
         ! Echo - Echo input to "<RootName>.AD.ech".
   
      CALL ReadVar( UnIn, InputFile, Echo, 'Echo',   'Echo flag', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   
      IF (.NOT. Echo .OR. I > 1) EXIT !exit this loop
   
         ! Otherwise, open the echo file, then rewind the input file and echo everything we've read
      
      I = I + 1         ! make sure we do this only once (increment counter that says how many times we've read this file)
   
      CALL OpenEcho ( UnEc, TRIM(OutFileRoot)//'.ech', ErrStat2, ErrMsg2, AD_Ver )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
   
      IF ( UnEc > 0 )  WRITE (UnEc,'(/,A,/)')  'Data from '//TRIM(AD_Ver%Name)//' primary input file "'//TRIM( InputFile )//'":'
         
      REWIND( UnIn, IOSTAT=ErrStat2 )  
         IF (ErrStat2 /= 0_IntKi ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error rewinding file "'//TRIM(InputFile)//'".', ErrStat, ErrMsg, RoutineName )
            CALL Cleanup()
            RETURN
         END IF         
      
   END DO    

   IF (NWTC_VerboseLevel == NWTC_Verbose) THEN
      CALL WrScr( ' Heading of the '//TRIM(aD_Ver%Name)//' input file: ' )      
      CALL WrScr( '   '//TRIM( FTitle ) )
   END IF
   
   
      ! DTAero - Time interval for aerodynamic calculations {or default} (s):
   Line = ""
   CALL ReadVar( UnIn, InputFile, Line, "DTAero", "Time interval for aerodynamic calculations {or default} (s)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
         
      CALL Conv2UC( Line )
      IF ( INDEX(Line, "DEFAULT" ) /= 1 ) THEN ! If it's not "default", read this variable; otherwise use the value already stored in InputFileData%DTAero
         READ( Line, *, IOSTAT=IOS) InputFileData%DTAero
            CALL CheckIOS ( IOS, InputFile, 'DTAero', NumType, ErrStat2, ErrMsg2 )
            CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      END IF   
      
      ! WakeMod - Type of wake/induction model {0=none, 1=BEMT} (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%WakeMod, "WakeMod", "Type of wake/induction model {0=none, 1=BEMT} (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! AFAeroMod - Type of airfoil aerodynamics model {1=steady model, 2=Beddoes-Leishman unsteady model} (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%AFAeroMod, "AFAeroMod", "Type of airfoil aerodynamics model {1=steady model, 2=Beddoes-Leishman unsteady model} (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! TwrInflnc - Type of tower influence model {0=none, 1=potential flow, 2=tower shadow} (switch) :
   CALL ReadVar( UnIn, InputFile, InputFileData%TwrInflnc, "TwrInflnc", "Type of tower influence model {0=none, 1=potential flow, 2=tower shadow} (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! TwrAero - Calculate tower aerodynamic loads? (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwrAero, "TwrAero", "Calculate tower aerodynamic loads? (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
            
   !----------- ENVIRONMENTAL CONDITIONS -------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Environmental Conditions', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! AirDens - Air density (kg/m^3):
   CALL ReadVar( UnIn, InputFile, InputFileData%AirDens, "AirDens", "Air density (kg/m^3)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! KinVisc - Kinematic air viscosity (m^2/s):
   CALL ReadVar( UnIn, InputFile, InputFileData%KinVisc, "KinVisc", "Kinematic air viscosity (m^2/s)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! SpdSound - Speed of sound (m/s):
   CALL ReadVar( UnIn, InputFile, InputFileData%SpdSound, "SpdSound", "Speed of sound (m/s)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF

   !----------- BLADE-ELEMENT/MOMENTUM THEORY OPTIONS ------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Blade-Element/Momentum Theory Options', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! SkewMod - Type of skewed-wake correction model {1=uncoupled, 2=Pitt/Peters, 3=coupled} [used only when WakeMod=1] (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%SkewMod, "SkewMod", "Type of skewed-wake correction model {1=uncoupled, 2=Pitt/Peters, 3=coupled} [used only when WakeMod=1] (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! TipLoss - Use the Prandtl tip-loss model? [used only when WakeMod=1] (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TipLoss, "TipLoss", "Use the Prandtl tip-loss model? [used only when WakeMod=1] (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! HubLoss - Use the Prandtl hub-loss model? [used only when WakeMod=1] (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%HubLoss, "HubLoss", "Use the Prandtl hub-loss model? [used only when WakeMod=1] (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! TanInd - Include tangential induction in BEMT calculations? [used only when WakeMod=1] (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TanInd, "TanInd", "Include tangential induction in BEMT calculations? [used only when WakeMod=1] (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! AIDrag - Include the drag term in the axial-induction calculation? [used only when WakeMod=1] (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%AIDrag, "AIDrag", "Include the drag term in the axial-induction calculation? [used only when WakeMod=1] (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! TIDrag - Include the drag term in the tangential-induction calculation? [used only when WakeMod=1 and TanInd=TRUE] (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TIDrag, "TIDrag", "Include the drag term in the tangential-induction calculation? [used only when WakeMod=1 and TanInd=TRUE] (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! IndToler - Convergence tolerance for BEM induction factors [used only when WakeMod=1] (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%IndToler, "IndToler", "Convergence tolerance for BEM induction factors [used only when WakeMod=1] (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! MaxIter - Maximum number of iteration steps [used only when WakeMod=1] (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%MaxIter, "MaxIter", "Maximum number of iteration steps [used only when WakeMod=1] (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )      
      
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
         
   !----------- BEDDOES-LEISHMAN UNSTEADY AIRFOIL AERODYNAMICS OPTIONS -------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Beddoes-Leishman Unsteady Airfoil Aerodynamics Options', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! UAMod - Unsteady Aero Model Switch (switch) {1=Baseline model (Original), 2=Gonzalez’s variant (changes in Cn,Cc,Cm), 3=Minemma/Pierce variant (changes in Cc and Cm)} [used only when AFAreoMod=2] (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%UAMod, "UAMod", "Unsteady Aero Model Switch (switch) {1=Baseline model (Original), 2=Gonzalez’s variant (changes in Cn,Cc,Cm), 3=Minemma/Pierce variant (changes in Cc and Cm)} [used only when AFAreoMod=2] (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! FLookup - Flag to indicate whether a lookup for f’ will be calculated (TRUE) or whether best-fit exponential equations will be used (FALSE); if FALSE S1-S4 must be provided in airfoil input files [used only when AFAreoMod=2] (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%FLookup, "FLookup", "Flag to indicate whether a lookup for f’ will be calculated (TRUE) or whether best-fit exponential equations will be used (FALSE); if FALSE S1-S4 must be provided in airfoil input files [used only when AFAreoMod=2] (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )     
      
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
                  
   !----------- AIRFOIL INFORMATION ------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Airfoil Information', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )


      ! InCol_Alfa - The column in the airfoil tables that contains the angle of attack (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%InCol_Alfa, "InCol_Alfa", "The column in the airfoil tables that contains the angle of attack (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! InCol_Cl - The column in the airfoil tables that contains the lift coefficient (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%InCol_Cl, "InCol_Cl", "The column in the airfoil tables that contains the lift coefficient (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! InCol_Cd - The column in the airfoil tables that contains the drag coefficient (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%InCol_Cd, "InCol_Cd", "The column in the airfoil tables that contains the drag coefficient (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! InCol_Cm - The column in the airfoil tables that contains the pitching-moment coefficient; use zero if there is no Cm column (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%InCol_Cm, "InCol_Cm", "The column in the airfoil tables that contains the pitching-moment coefficient; use zero if there is no Cm column (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! InCol_Cpmin - The column in the airfoil tables that contains the drag coefficient; use zero if there is no Cpmin column (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%InCol_Cpmin, "InCol_Cpmin", "The column in the airfoil tables that contains the drag coefficient; use zero if there is no Cpmin column (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NumAFfiles - Number of airfoil files used (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%NumAFfiles, "NumAFfiles", "Number of airfoil files used (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

         ! Allocate space to hold AFNames
      ALLOCATE( InputFileData%AFNames(InputFileData%NumAFfiles), STAT=ErrStat2)
         IF (ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, "Error allocating AFNames.", ErrStat, ErrMsg, RoutineName)
            CALL Cleanup()
            RETURN
         END IF
               
      ! AFNames - Airfoil file names (NumAFfiles lines) (quoted strings):
   DO I = 1,InputFileData%NumAFfiles            
      CALL ReadVar ( UnIn, InputFile, InputFileData%AFNames(I), 'AFNames('//TRIM(Num2Lstr(I))//')', 'Airfoil '//TRIM(Num2Lstr(I))//' file name', ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( PathIsRelative( InputFileData%AFNames(I) ) ) InputFileData%AFNames(I) = TRIM(PriPath)//TRIM(InputFileData%AFNames(I))
   END DO      
             
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF

   !----------- ROTOR/BLADE PROPERTIES  --------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Rotor/Blade Properties', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! UseBlCm - Include aerodynamic pitching moment in calculations? (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%UseBlCm, "UseBlCm", "Include aerodynamic pitching moment in calculations? (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN
            
   !   ! NumBlNds - Number of blade nodes used in the analysis (-):
   !CALL ReadVar( UnIn, InputFile, InputFileData%NumBlNds, "NumBlNds", "Number of blade nodes used in the analysis (-)", ErrStat2, ErrMsg2, UnEc)
   !   CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   !   IF ( ErrStat >= AbortErrLev ) RETURN
      
      ! ADBlFile - Names of files containing distributed aerodynamic properties for each blade (see AD_BladeInputFile type):
   DO I = 1,MaxBl            
      CALL ReadVar ( UnIn, InputFile, ADBlFile(I), 'ADBlFile('//TRIM(Num2Lstr(I))//')', 'Name of file containing distributed aerodynamic properties for blade '//TRIM(Num2Lstr(I)), ErrStat2, ErrMsg2, UnEc )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( PathIsRelative( ADBlFile(I) ) ) ADBlFile(I) = TRIM(PriPath)//TRIM(ADBlFile(I))
   END DO      
      
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF

   !----------- TOWER INFLUENCE AND AERODYNAMICS  ----------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Tower Influence and Aerodynamics', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! TwrPFBak - Include Bak correction to potential flow? [used only when TwrInflnc=1] (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%TwrPFBak, "TwrPFBak", "Include Bak correction to potential flow? [used only when TwrInflnc=1] (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      ! NumTwrNds - Number of tower nodes used in the analysis (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%NumTwrNds, "NumTwrNds", "Number of tower nodes used in the analysis (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

      
   !....... tower properties ...................
   CALL ReadCom( UnIn, InputFile, 'Section Header: Tower Property Channels', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL ReadCom( UnIn, InputFile, 'Section Header: Tower Property Units', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! allocate space for tower inputs:
   CALL AllocAry( InputFileData%TwrElev,  InputFileData%NumTwrNds, 'TwrElev',  ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( InputFileData%TwrDiam, InputFileData%NumTwrNds, 'TwrDiam', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( InputFileData%TwrCd, InputFileData%NumTwrNds, 'TwrCd', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! Return on error if we didn't allocate space for the next inputs
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
            
   DO I=1,InputFileData%NumTwrNds
      call ReadAry ( UnIn, InputFile, TmpAry,  3, 'TwrNds',  'Properties for tower node ' &
                     //trim( Int2LStr( I ) )//'.', errStat2, errMsg2, UnEc )
         call setErrStat( errStat2, ErrMsg2 , errStat, ErrMsg , RoutineName )
            
      InputFileData%TwrElev(I) = TmpAry( 1)
      InputFileData%TwrDiam(I) = TmpAry( 2)
      InputFileData%TwrCd(I)   = TmpAry( 3)      
   END DO
               
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
                  
   !----------- OUTPUTS  -----------------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: Outputs', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! SumPrint - Generate a summary file listing input options and interpolated properties to <rootname>.AD.sum? (flag):
   CALL ReadVar( UnIn, InputFile, InputFileData%SumPrint, "SumPrint", "Generate a summary file listing input options and interpolated properties to <rootname>.AD.sum? (flag)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      ! NBlOuts - Number of blade node outputs [0 - 9] (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%NBlOuts, "NBlOuts", "Number of blade node outputs [0 - 9] (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      IF ( InputFileData%NBlOuts > SIZE(InputFileData%BlOutNd) ) THEN
         CALL SetErrStat( ErrID_Warn, ' Warning: number of blade output nodes exceeds '//&
                           TRIM(Num2LStr(SIZE(InputFileData%BlOutNd))) //'.', ErrStat, ErrMsg, RoutineName )
         InputFileData%NBlOuts = SIZE(InputFileData%BlOutNd)
      END IF
      
      ! BlOutNd - Blade nodes whose values will be output (-):
   CALL ReadAry( UnIn, InputFile, InputFileData%BlOutNd, InputFileData%NBlOuts, "BlOutNd", "Blade nodes whose values will be output (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
            
      ! NTwOuts - Number of tower node outputs [0 - 9] (-):
   CALL ReadVar( UnIn, InputFile, InputFileData%NTwOuts, "NTwOuts", "Number of tower node outputs [0 - 9] (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      IF ( InputFileData%NTwOuts > SIZE(InputFileData%TwOutNd) ) THEN
         CALL SetErrStat( ErrID_Warn, ' Warning: number of tower output nodes exceeds '//&
                           TRIM(Num2LStr(SIZE(InputFileData%TwOutNd))) //'.', ErrStat, ErrMsg, RoutineName )
         InputFileData%NTwOuts = SIZE(InputFileData%TwOutNd)
      END IF
      
      ! TwOutNd - Tower nodes whose values will be output (-):
   CALL ReadAry( UnIn, InputFile, InputFileData%TwOutNd, InputFileData%NTwOuts, "TwOutNd", "Tower nodes whose values will be output (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! Return on error at end of section
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
                  
   !----------- OUTLIST  -----------------------------------------------------------
   CALL ReadCom( UnIn, InputFile, 'Section Header: OutList', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! OutList - List of user-requested output channels (-):
   CALL ReadOutputList ( UnIn, InputFile, InputFileData%OutList, InputFileData%NumOuts, 'OutList', "List of user-requested output channels", ErrStat2, ErrMsg2, UnEc  )     ! Routine in NWTC Subroutine Library
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   !---------------------- END OF FILE -----------------------------------------
      
   CALL Cleanup( )
   RETURN


CONTAINS
   !...............................................................................................................................
   SUBROUTINE Cleanup()
   ! This subroutine cleans up any local variables and closes input files
   !...............................................................................................................................

   IF (UnIn > 0) CLOSE ( UnIn )

   END SUBROUTINE Cleanup
   !...............................................................................................................................
END SUBROUTINE ReadPrimaryFile      
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadBladeInputs ( ADBlFile, BladeKInputFileData, UnEc, ErrStat, ErrMsg )
! This routine reads a blade input file.
!..................................................................................................................................


      ! Passed variables:

   TYPE(AD_BladePropsType),  INTENT(INOUT)  :: BladeKInputFileData                 ! Data for Blade K stored in the module's input file
   CHARACTER(*),             INTENT(IN)     :: ADBlFile                            ! Name of the blade input file data
   INTEGER(IntKi),           INTENT(IN)     :: UnEc                                ! I/O unit for echo file. If present and > 0, write to UnEc

   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                             ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                              ! Error message


      ! Local variables:

   INTEGER(IntKi)               :: I                                               ! A generic DO index.
   INTEGER( IntKi )             :: UnIn                                            ! Unit number for reading file
   INTEGER(IntKi)               :: ErrStat2 , IOS                                  ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                                         ! Temporary Err msg
   CHARACTER(*), PARAMETER      :: RoutineName = 'ReadBladeInputs'

   ErrStat = ErrID_None
   ErrMsg  = ""
   UnIn = -1
      
   ! Allocate space for these variables
   
   
   
   
   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)


      ! Open the input file for blade K.

   CALL OpenFInpFile ( UnIn, ADBlFile, ErrStat2, ErrMsg2 )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      IF ( ErrStat >= AbortErrLev ) RETURN


   !  -------------- HEADER -------------------------------------------------------

      ! Skip the header.

   CALL ReadCom ( UnIn, ADBlFile, 'unused blade file header line 1', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   CALL ReadCom ( UnIn, ADBlFile, 'unused blade file header line 2', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   !  -------------- Blade properties table ------------------------------------------                                    
   CALL ReadCom ( UnIn, ADBlFile, 'Section header: Blade Properties', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      ! NumBlNds - Number of blade nodes used in the analysis (-):
   CALL ReadVar( UnIn, ADBlFile, BladeKInputFileData%NumBlNds, "NumBlNds", "Number of blade nodes used in the analysis (-)", ErrStat2, ErrMsg2, UnEc)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF ( ErrStat >= AbortErrLev ) RETURN

   CALL ReadCom ( UnIn, ADBlFile, 'Table header: names', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   CALL ReadCom ( UnIn, ADBlFile, 'Table header: units', ErrStat2, ErrMsg2, UnEc )
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
   IF ( ErrStat>= AbortErrLev ) THEN 
      CALL Cleanup()
      RETURN
   END IF
   
      
      ! allocate space for blade inputs:
   CALL AllocAry( BladeKInputFileData%BlSpn,   BladeKInputFileData%NumBlNds, 'BlSpn',   ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlCrvAC, BladeKInputFileData%NumBlNds, 'BlCrvAC', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlSwpAC, BladeKInputFileData%NumBlNds, 'BlSwpAC', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlCrvAng,BladeKInputFileData%NumBlNds, 'BlCrvAng',ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlTwist, BladeKInputFileData%NumBlNds, 'BlTwist', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlChord, BladeKInputFileData%NumBlNds, 'BlChord', ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   CALL AllocAry( BladeKInputFileData%BlAFID,  BladeKInputFileData%NumBlNds, 'BlAFID',  ErrStat2, ErrMsg2)
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      ! Return on error if we didn't allocate space for the next inputs
   IF ( ErrStat >= AbortErrLev ) THEN
      CALL Cleanup()
      RETURN
   END IF
            
   DO I=1,BladeKInputFileData%NumBlNds
      READ( UnIn, *, IOStat=IOS ) BladeKInputFileData%BlSpn(I), BladeKInputFileData%BlCrvAC(I), BladeKInputFileData%BlSwpAC(I), &
                                  BladeKInputFileData%BlCrvAng(I), BladeKInputFileData%BlTwist(I), BladeKInputFileData%BlChord(I), &
                                  BladeKInputFileData%BlAFID(I)  
         CALL CheckIOS( IOS, ADBlFile, 'Blade properties row '//TRIM(Num2LStr(I)), NumType, ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
               ! Return on error if we couldn't read this line
            IF ( ErrStat >= AbortErrLev ) THEN
               CALL Cleanup()
               RETURN
            END IF
         
         IF (UnEc > 0) THEN
            WRITE( UnEc, "(6(F9.4,1x),I9)", IOStat=IOS) BladeKInputFileData%BlSpn(I), BladeKInputFileData%BlCrvAC(I), BladeKInputFileData%BlSwpAC(I), &
                                  BladeKInputFileData%BlCrvAng(I), BladeKInputFileData%BlTwist(I), BladeKInputFileData%BlChord(I), &
                                  BladeKInputFileData%BlAFID(I)
         END IF         
   END DO
   BladeKInputFileData%BlCrvAng = BladeKInputFileData%BlCrvAng*D2R
   BladeKInputFileData%BlTwist  = BladeKInputFileData%BlTwist*D2R
                  
   !  -------------- END OF FILE --------------------------------------------

   CALL Cleanup()
   RETURN


CONTAINS
   !...............................................................................................................................
   SUBROUTINE Cleanup()
   ! This subroutine cleans up local variables and closes files
   !...............................................................................................................................

      IF (UnIn > 0) CLOSE(UnIn)

   END SUBROUTINE Cleanup

END SUBROUTINE ReadBladeInputs      
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ValidateInputData( InputFileData, NumBl, ErrStat, ErrMsg )
! This routine validates the inputs from the AeroDyn input files.
!..................................................................................................................................
      
      ! Passed variables:

   type(AD_InputFile),       intent(in)     :: InputFileData                       ! All the data in the AeroDyn input file
   integer(IntKi),           intent(in)     :: NumBl                               ! Number of blades
   integer(IntKi),           intent(out)    :: ErrStat                             ! Error status
   character(*),             intent(out)    :: ErrMsg                              ! Error message

   
      ! local variables
   integer(IntKi)                           :: k                                   ! Blade number
   integer(IntKi)                           :: j                                   ! node number
   character(*), parameter                  :: RoutineName = 'ValidateInputData'
   
   ErrStat = ErrID_None
   ErrMsg  = ""
         
   if (InputFileData%DTAero <= 0.0)  call SetErrStat ( ErrID_Fatal, 'DTAero must be greater than zero.', ErrStat, ErrMsg, RoutineName )
   if (InputFileData%WakeMod /= WakeMod_None .and. InputFileData%WakeMod /= WakeMod_BEMT) call SetErrStat ( ErrID_Fatal, &
      'WakeMod must '//trim(num2lstr(WakeMod_None))//' (none) or '//trim(num2lstr(WakeMod_BEMT))//' (BEMT).', ErrStat, ErrMsg, RoutineName ) 
   if (InputFileData%AFAeroMod /= AFAeroMod_Steady .and. InputFileData%AFAeroMod /= AFAeroMod_BL_unsteady) then
      call SetErrStat ( ErrID_Fatal, 'AFAeroMod must be '//trim(num2lstr(AFAeroMod_Steady))//' (steady) or '//&
                        trim(num2lstr(AFAeroMod_BL_unsteady))//' (Beddoes-Leishman unsteady).', ErrStat, ErrMsg, RoutineName ) 
   end if
   if (InputFileData%TwrInflnc /= TwrInflnc_none .and. InputFileData%TwrInflnc /= TwrInflnc_Potent .and. InputFileData%TwrInflnc /= TwrInflnc_Shadow) then
      call SetErrStat ( ErrID_Fatal, 'TwrInflnc must be 0 (none), 1 (potential flow), or 2 (tower shadow).', ErrStat, ErrMsg, RoutineName ) 
   end if   
   
   if (InputFileData%AirDens <= 0.0) call SetErrStat ( ErrID_Fatal, 'The air density (AirDens) must be greater than zero.', ErrStat, ErrMsg, RoutineName )
   if (InputFileData%KinVisc <= 0.0) call SetErrStat ( ErrID_Fatal, 'The kinesmatic viscosity (KinVisc) must be greater than zero.', ErrStat, ErrMsg, RoutineName )
   if (InputFileData%SpdSound <= 0.0) call SetErrStat ( ErrID_Fatal, 'The speed of sound (SpdSound) must be greater than zero.', ErrStat, ErrMsg, RoutineName )
      
   
      ! BEMT inputs
      ! these checks should probably go into BEMT where they are used...
   if (InputFileData%WakeMod == WakeMod_BEMT) then
      if ( InputFileData%MaxIter < 1 ) call SetErrStat( ErrID_Fatal, 'MaxIter must be greater than 0.', ErrStat, ErrMsg, RoutineName )
      
      if ( InputFileData%IndToler < 0.0 .or. EqualRealNos(InputFileData%IndToler, 0.0_ReKi) ) &
         call SetErrStat( ErrID_Fatal, 'IndToler must be greater than 0.', ErrStat, ErrMsg, RoutineName )
   
      if ( InputFileData%SkewMod /= SkewMod_Uncoupled .and. InputFileData%SkewMod /= SkewMod_PittPeters .and. &
           InputFileData%SkewMod /= SkewMod_Coupled ) call SetErrStat( ErrID_Fatal, 'SkewMod must be 1, 2, or 3.', ErrStat, ErrMsg, RoutineName )      
      
   end if !BEMT checks
   
      ! UA inputs
   if (InputFileData%AFAeroMod == AFAeroMod_BL_unsteady ) then
      if (InputFileData%UAMod < 1 .or. InputFileData%UAMod > 3 ) call SetErrStat( ErrID_Fatal, &
         "UAMod must be 1 (baseline/original), 2 (Gonzalez's variant), or 3 (Minemma/Pierce variant).", ErrStat, ErrMsg, RoutineName )
   end if
   
   
         ! validate the AFI input data because it doesn't appear to be done in AFI
   if (InputFileData%NumAFfiles < 1) call SetErrStat( ErrID_Fatal, 'The number of unique airfoil tables (NumAFfiles) must be greater than zero.', ErrStat, ErrMsg, RoutineName )   
   if (InputFileData%InCol_Alfa  < 0) call SetErrStat( ErrID_Fatal, 'InCol_Alfa must not be a negative number.', ErrStat, ErrMsg, RoutineName )
   if (InputFileData%InCol_Cl    < 0) call SetErrStat( ErrID_Fatal, 'InCol_Cl must not be a negative number.', ErrStat, ErrMsg, RoutineName )
   if (InputFileData%InCol_Cd    < 0) call SetErrStat( ErrID_Fatal, 'InCol_Cd must not be a negative number.', ErrStat, ErrMsg, RoutineName )
   if (InputFileData%InCol_Cm    < 0) call SetErrStat( ErrID_Fatal, 'InCol_Cm must not be a negative number.', ErrStat, ErrMsg, RoutineName )
   if (InputFileData%InCol_Cpmin < 0) call SetErrStat( ErrID_Fatal, 'InCol_Cpmin must not be a negative number.', ErrStat, ErrMsg, RoutineName )
   
      ! .............................
      ! check blade mesh data:
      ! .............................
   if ( InputFileData%BladeProps(1)%NumBlNds < 2 ) call SetErrStat( ErrID_Fatal, 'There must be at least two nodes per blade.',ErrStat, ErrMsg, RoutineName )
   do k=2,NumBl
      if ( InputFileData%BladeProps(k)%NumBlNds /= InputFileData%BladeProps(k-1)%NumBlNds ) then
         call SetErrStat( ErrID_Fatal, 'All blade property files must have the same number of blade nodes.', ErrStat, ErrMsg, RoutineName )
         exit  ! exit do loop
      end if
   end do
   
      ! Check the list of airfoil tables for blades to make sure they are all within limits.
   do k=1,NumBl
      do j=1,InputFileData%BladeProps(k)%NumBlNds
         if ( ( InputFileData%BladeProps(k)%BlAFID(j) < 1 ) .OR. ( InputFileData%BladeProps(k)%BlAFID(j) > InputFileData%NumAFfiles ) )  then
            call SetErrStat( ErrID_Fatal, 'Blade '//trim(Num2LStr(k))//' node '//trim(Num2LStr(j))//' must be a number between 1 and NumAFfiles (' &
               //TRIM(Num2LStr(InputFileData%NumAFfiles))//').', ErrStat, ErrMsg, RoutineName )
         end if
      end do ! j=nodes
   end do ! k=blades
            
      ! Check that the blade chord is > 0.
   do k=1,NumBl
      do j=1,InputFileData%BladeProps(k)%NumBlNds
         if ( InputFileData%BladeProps(k)%BlChord(j) <= 0.0_ReKi )  then
            call SetErrStat( ErrID_Fatal, 'The chord for blade '//trim(Num2LStr(k))//' node '//trim(Num2LStr(j)) &
                             //' must be greater than 0.', ErrStat, ErrMsg, RoutineName )
         endif
      end do ! j=nodes
   end do ! k=blades
   
   do k=1,NumBl
      if ( .not. EqualRealNos(InputFileData%BladeProps(k)%BlSpn(1), 0.0_ReKi) ) call SetErrStat( ErrID_Fatal, 'Blade '//trim(Num2LStr(k))//' span location must start at 0.0 m', ErrStat, ErrMsg, RoutineName)       
      do j=2,InputFileData%BladeProps(k)%NumBlNds
         if ( InputFileData%BladeProps(k)%BlSpn(j) <= InputFileData%BladeProps(k)%BlSpn(j-1) )  then
            call SetErrStat( ErrID_Fatal, 'Blade '//trim(Num2LStr(k))//' nodes must be entered in increasing elevation.', ErrStat, ErrMsg, RoutineName )
            exit
         end if
      end do ! j=nodes
   end do ! k=blades
   
      ! .............................
      ! check tower mesh data:
      ! .............................
   
   if ( InputFileData%TwrInflnc /= TwrInflnc_None .or. InputFileData%TwrAero ) then
      
      if (InputFileData%NumTwrNds < 2) call SetErrStat( ErrID_Fatal, 'There must be at least two nodes on the tower.',ErrStat, ErrMsg, RoutineName )
         
         ! Check that the tower diameter is > 0.
      do j=1,InputFileData%NumTwrNds
         if ( InputFileData%TwrDiam(j) <= 0.0_ReKi )  then
            call SetErrStat( ErrID_Fatal, 'The diameter for tower node '//trim(Num2LStr(j))//' must be greater than 0.' &
                            , ErrStat, ErrMsg, RoutineName )
         end if
      end do ! j=nodes
      
         ! check that the elevation starts at 0 and is increasing:
      if ( .not. EqualRealNos(InputFileData%TwrElev(1), 0.0_ReKi) ) call SetErrStat( ErrID_Fatal, 'The first tower elevation must be 0.0 m', ErrStat, ErrMsg, RoutineName) 
      do j=2,InputFileData%NumTwrNds
         if ( InputFileData%TwrElev(j) <= InputFileData%TwrElev(j-1) )  then
            call SetErrStat( ErrID_Fatal, 'The tower nodes must be entered in increasing elevation.', ErrStat, ErrMsg, RoutineName )
            exit
         end if
      end do ! j=nodes
            
   end if
   
         
END SUBROUTINE ValidateInputData
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Init_AFIparams( InputFileData, p_AFI, UnEc, NumBl, ErrStat, ErrMsg )


      ! Passed variables
   type(AD_InputFile),      intent(inout)   :: InputFileData      ! All the data in the AeroDyn input file (intent(out) only because of the call to MOVE_ALLOC)
   type(AFI_ParameterType), intent(  out)   :: p_AFI              ! parameters returned from the AFI (airfoil info) module
   integer(IntKi),          intent(in   )   :: UnEc               ! I/O unit for echo file. If > 0, file is open for writing.
   integer(IntKi),          intent(in   )   :: NumBl              ! number of blades (for performing check on valid airfoil data read in)
   integer(IntKi),          intent(  out)   :: ErrStat            ! Error status
   character(*),            intent(  out)   :: ErrMsg             ! Error message

      ! local variables
   type(AFI_InitInputType)                  :: AFI_InitInputs     ! initialization data for the AFI routines
   
   integer(IntKi)                           :: j                  ! loop counter for nodes
   integer(IntKi)                           :: k                  ! loop counter for blades
   integer(IntKi)                           :: File               ! loop counter for airfoil files
   integer(IntKi)                           :: Table              ! loop counter for airfoil tables in a file
   logical, allocatable                     :: fileUsed(:)
   
   integer(IntKi)                           :: ErrStat2
   character(ErrMsgLen)                     :: ErrMsg2
   character(*), parameter                  :: RoutineName = 'Init_AFIparams'

   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! Setup Airfoil InitInput data structure:
   AFI_InitInputs%NumAFfiles = InputFileData%NumAFfiles
   call MOVE_ALLOC( InputFileData%AFNames, AFI_InitInputs%FileNames ) ! move from AFNames to FileNames      
   AFI_InitInputs%InCol_Alfa  = InputFileData%InCol_Alfa
   AFI_InitInputs%InCol_Cl    = InputFileData%InCol_Cl
   AFI_InitInputs%InCol_Cd    = InputFileData%InCol_Cd
   AFI_InitInputs%InCol_Cm    = InputFileData%InCol_Cm
   AFI_InitInputs%InCol_Cpmin = InputFileData%InCol_Cpmin
               
      ! Call AFI_Init to read in and process the airfoil files.
      ! This includes creating the spline coefficients to be used for interpolation.

   call AFI_Init ( AFI_InitInputs, p_AFI, ErrStat2, ErrMsg2, UnEc )
      call SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)   
   
      
   call MOVE_ALLOC( AFI_InitInputs%FileNames, InputFileData%AFNames ) ! move from FileNames back to AFNames
   call AFI_DestroyInitInput( AFI_InitInputs, ErrStat2, ErrMsg2 )
   
   if (ErrStat >= AbortErrLev) return
   
   
      ! check that we read the correct airfoil parameters for UA:      
   if ( InputFileData%AFAeroMod == AFAeroMod_BL_unsteady ) then
      
         
         ! determine which airfoil files will be used
      call AllocAry( fileUsed, InputFileData%NumAFfiles, 'fileUsed', errStat2, errMsg2 )
         call SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)   
      if (errStat >= AbortErrLev) return
      fileUsed = .false.
            
      do k=1,NumBl
         do j=1,InputFileData%BladeProps(k)%NumBlNds
            fileUsed ( InputFileData%BladeProps(k)%BlAFID(j) ) = .true.
         end do ! j=nodes
      end do ! k=blades
      
         ! make sure all files in use have UA input parameters:
      do File = 1,InputFileData%NumAFfiles
         
         if (fileUsed(File)) then
            do Table=1,p_AFI%AFInfo(File)%NumTabs            
               if ( .not. p_AFI%AFInfo(File)%Table(Table)%InclUAdata ) then
                  call SetErrStat( ErrID_Fatal, 'Airfoil file '//trim(InputFileData%AFNames(File))//', table #'// &
                        trim(num2lstr(Table))//' does not contain parameters for UA data.', ErrStat, ErrMsg, RoutineName )
               end if
            end do
         end if
         
      end do
      
      if ( allocated(fileUsed) ) deallocate(fileUsed)
      
   end if
   
   
END SUBROUTINE Init_AFIparams
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Init_BEMTmodule( InputFileData, u_AD, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
! This routine initializes the BEMT module from within AeroDyn
!..................................................................................................................................

   type(AD_InputFile),             intent(in   ) :: InputFileData  ! All the data in the AeroDyn input file
   type(AD_InputType),             intent(in   ) :: u_AD           ! AD inputs - used for input mesh node positions
   type(BEMT_InputType),           intent(  out) :: u              ! An initial guess for the input; input mesh must be defined
   type(AD_ParameterType),         intent(inout) :: p              ! Parameters ! intent out b/c we set the BEMT parameters here
   type(BEMT_ContinuousStateType), intent(  out) :: x              ! Initial continuous states
   type(BEMT_DiscreteStateType),   intent(  out) :: xd             ! Initial discrete states
   type(BEMT_ConstraintStateType), intent(  out) :: z              ! Initial guess of the constraint states
   type(BEMT_OtherStateType),      intent(  out) :: OtherState     ! Initial other/optimization states
   type(BEMT_OutputType),          intent(  out) :: y              ! Initial system outputs (outputs are not calculated;
                                                                   !   only the output mesh is initialized)
   integer(IntKi),                 intent(  out) :: errStat        ! Error status of the operation
   character(*),                   intent(  out) :: errMsg         ! Error message if ErrStat /= ErrID_None


      ! Local variables
   real(DbKi)                                    :: Interval       ! Coupling interval in seconds: the rate that
                                                                   !   (1) BEMT_UpdateStates() is called in loose coupling &
                                                                   !   (2) BEMT_UpdateDiscState() is called in tight coupling.
                                                                   !   Input is the suggested time from the glue code;
                                                                   !   Output is the actual coupling interval that will be used
                                                                   !   by the glue code.
   type(BEMT_InitInputType)                      :: InitInp        ! Input data for initialization routine
   type(BEMT_InitOutputType)                     :: InitOut        ! Output for initialization routine
                                                 
   integer(intKi)                                :: j              ! node index
   integer(intKi)                                :: k              ! blade index
   integer(IntKi)                                :: ErrStat2
   character(ErrMsgLen)                          :: ErrMsg2
   character(*), parameter                       :: RoutineName = 'Init_BEMTmodule'

   ! note here that each blade is required to have the same number of nodes
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   
      ! set initialization data here:   
   Interval                 = p%DT   
   InitInp%numBlades        = p%NumBlades
   
   InitInp%airDens          = InputFileData%AirDens 
   InitInp%kinVisc          = InputFileData%KinVisc                  
   InitInp%skewWakeMod      = InputFileData%SkewMod      ! bjj: FIX ME: check what these variables mean.... and make sure AD input file matches BEMT code
   InitInp%aTol             = InputFileData%IndToler
   InitInp%useTipLoss       = InputFileData%TipLoss
   InitInp%useHubLoss       = InputFileData%HubLoss
   InitInp%useInduction     = InputFileData%WakeMod == 1
   InitInp%useTanInd        = InputFileData%TanInd
   InitInp%useAIDrag        = InputFileData%AIDrag        
   InitInp%useTIDrag        = InputFileData%TIDrag  
   InitInp%numBladeNodes    = p%NumBlNds
   InitInp%numReIterations  = 1                              ! This is currently not available in the input file and is only for testing  
   InitInp%maxIndIterations = InputFileData%MaxIter 
   
   call AllocAry(InitInp%chord, InitInp%numBladeNodes,InitInp%numBlades,'chord', ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
   call AllocAry(InitInp%AFindx,InitInp%numBladeNodes,InitInp%numBlades,'AFindx',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
   call AllocAry(InitInp%zHub,                        InitInp%numBlades,'zHub',  ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   call AllocAry(InitInp%zLocal,InitInp%numBladeNodes,InitInp%numBlades,'zLocal',ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)   
   call AllocAry(InitInp%zTip,                        InitInp%numBlades,'zTip',  ErrStat2,ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)

   if ( ErrStat >= AbortErrLev ) then
      call Cleanup()
      return
   end if  

   do k=1,p%numBlades
      
      InitInp%zHub(k) = TwoNorm( u_AD%BladeRootMotion(k)%Position(:,1) - u_AD%HubMotion%Position(:,1) )  
      
      InitInp%zLocal(1,k) = InitInp%zHub(k) + TwoNorm( u_AD%BladeMotion(k)%Position(:,1) - u_AD%BladeRootMotion(k)%Position(:,1) )
      do j=2,p%NumBlNds
         InitInp%zLocal(j,k) = InitInp%zLocal(j-1,k) + TwoNorm( u_AD%BladeMotion(k)%Position(:,j) - u_AD%BladeMotion(k)%Position(:,j-1) ) 
      end do !j=nodes
      
      InitInp%zTip(k) = InitInp%zLocal(p%NumBlNds,k)
      
   end do !k=blades
   
   
               
  do k=1,p%numBlades
     do j=1,p%NumBlNds
        InitInp%chord (j,k)  = InputFileData%BladeProps(k)%BlChord(j)
        InitInp%AFindx(j,k)  = InputFileData%BladeProps(k)%BlAFID(j)
     end do
  end do
   
   InitInp%UA_Flag = InputFileData%AFAeroMod == AFAeroMod_BL_unsteady
   InitInp%UAMod   = InputFileData%UAMod
   InitInp%Flookup = InputFileData%Flookup
   InitInp%a_s     = InputFileData%SpdSound
   
   
   call BEMT_Init(InitInp, u, p%BEMT,  x, xd, z, OtherState, y, Interval, InitOut, ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2,ErrMsg2, ErrStat, ErrMsg, RoutineName)   
         
   if (.not. equalRealNos(Interval, p%DT) ) &
      call SetErrStat( ErrID_Fatal, "DTAero was changed in Init_BEMTmodule(); this is not allowed.", ErrStat2, ErrMsg2, RoutineName)
   
   call Cleanup()
   return
      
contains   
   subroutine Cleanup()
      call BEMT_DestroyInitInput( InitInp, ErrStat2, ErrMsg2 )   
      call BEMT_DestroyInitOutput( InitOut, ErrStat2, ErrMsg2 )   
   end subroutine Cleanup
   
END SUBROUTINE Init_BEMTmodule
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ADTwr_CalcOutput(p, u, OtherState, y, ErrStat, ErrMsg )

   TYPE(AD_InputType),           INTENT(IN   )  :: u           ! Inputs at Time t
   TYPE(AD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(AD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(AD_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                               !   nectivity information does not have to be recalculated)
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


   INTEGER(IntKi)                               :: j
   real(ReKi)                                   :: q
   real(ReKi)                                   :: W           ! relative wind speed normal to the tower at a tower node 
   real(ReKi)                                   :: V_rel(3)    ! relative wind speed on a tower node
   real(ReKi)                                   :: VL(2)       ! relative local x- and y-components of the wind speed on a tower node
   
   !integer(intKi)                               :: ErrStat2
   !character(ErrMsgLen)                         :: ErrMsg2
   character(*), parameter                      :: RoutineName = 'ADTwr_CalcOutput'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""

   
   do j=1,p%NumTwrNds
      
      V_rel = u%InflowOnTower(:,j) - u%TowerMotion%TranslationDisp(:,j) ! relative wind speed at tower node
   
      VL(1) = dot_product( V_Rel, u%TowerMotion%Orientation(1,:,j) ) ! relative local x-component of wind speed of the jth node in the tower
      VL(2) = dot_product( V_Rel, u%TowerMotion%Orientation(2,:,j) ) ! relative local y-component of wind speed of the jth node in the tower
      
      W     =  TwoNorm( VL )  ! relative wind speed normal to the tower at node j      
      q     = 0.5 * p%TwrCd(j) * p%AirDens * p%TwrDiam(j) * W
      
         ! force per unit length of the jth node in the tower
      y%TowerLoad%force(1,j) = q * VL(1)
      y%TowerLoad%force(2,j) = q * VL(2)
      y%TowerLoad%force(3,j) = 0.0_ReKi
      
         ! moment per unit length of the jth node in the tower
      y%TowerLoad%moment(:,j) = 0.0_ReKi
      
   end do
   

END SUBROUTINE ADTwr_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE AeroDyn