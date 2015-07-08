!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2012-2015  National Renewable Energy Laboratory
!
!    This file is part of AeroDyn.
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
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE AeroDyn14

   USE AeroDyn14_Types
   USE AeroSubs
   USE NWTC_Library


   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: AD14_Ver = ProgDesc( 'AeroDyn14', 'v14.04.00a-bjj', '14-Apr-2015' )

      ! ..... Public Subroutines ............

   PUBLIC :: AD14_Init                           ! Initialization routine
   PUBLIC :: AD14_End                            ! Ending routine (includes clean up)

   PUBLIC :: AD14_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                 !   continuous states, and updating discrete states
   PUBLIC :: AD14_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: AD14_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: AD14_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: AD14_UpdateDiscState                ! Tight coupling routine for updating discrete states


CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD14_Init( InitInp, u, p, x, xd, z, O, y, Interval, InitOut, ErrStat, ErrMess )
!..................................................................................................................................
   USE               AeroGenSubs,   ONLY: ElemOpen
   USE DWM
   IMPLICIT NONE

   TYPE(AD14_InitInputType),       INTENT(INOUT)  :: InitInp     ! Input data for initialization routine
   TYPE(AD14_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
   TYPE(AD14_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
   TYPE(AD14_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
   TYPE(AD14_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
   TYPE(AD14_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
   TYPE(AD14_OtherStateType),      INTENT(  OUT)  :: O !therState  Initial other/optimization states
   TYPE(AD14_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                               !   only the output mesh is initialized)
   REAL(DbKi),                     INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                                 !   (1) AD14_UpdateStates() is called in loose coupling &
                                                                 !   (2) AD14_UpdateDiscState() is called in tight coupling.
                                                                 !   Input is the suggested time from the glue code;
                                                                 !   Output is the actual coupling interval that will be used
                                                                 !   by the glue code.
   TYPE(AD14_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMess     ! Error message if ErrStat /= ErrID_None


      ! Internal variables
   REAL(ReKi)                       :: CosPrecone
   REAL(ReKi)                       :: DTip, ElemRad, Dhub, Rhub     ! variables for calculating hub- and tip-loss constants
   REAL(ReKi)                       :: HubRadius
!   REAL(ReKi)                       :: MeanWind
   REAL(ReKi)                       :: TipRadius
   REAL(ReKi)                       :: TmpVar
   REAL(ReKi)                       :: TmpPos(3)
   REAL(ReKi)                                :: TwrNodeHt                     ! The height of the current tower node.

   INTEGER                          :: IB, IE 
   INTEGER                          :: IELM

!   CHARACTER(1024)                  :: Title

   INTEGER                                   :: Elem                          ! Index for mesh element.
   INTEGER                                   :: InterpIndx                 ! Index telling the interpolation routine where to start in the array.
   INTEGER                                   :: Node                          ! Index used to pull points out of the array of values at given node location.
   INTEGER                                   :: ErrStatLcL        ! Error status returned by called routines.

   CHARACTER(ErrMsgLen)                      :: ErrMessLcl          ! Error message returned by called routines.
   CHARACTER(*), PARAMETER                   :: RoutineName = 'AD14_Init'

         ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMess  = ""
   InterpIndx = 1
   !-------------------------------------------------------------------------------------------------
   ! Check that the module hasn't already been initialized.
   !-------------------------------------------------------------------------------------------------

   IF ( p%Initialized ) THEN
      CALL SetErrStat( ErrID_Warn,'AeroDyn has already been initialized.',ErrStat,ErrMess,RoutineName)
      RETURN
   ELSE
      p%Initialized = .TRUE.
      CALL NWTC_Init( )
   END IF

   
         ! Display the module information

   CALL DispNVD( AD14_Ver )
   
   InitOut%Ver = AD14_Ver
   O%FirstWarn = .TRUE.
   !-------------------------------------------------------------------------------------------------
   ! Set up AD variables
   !-------------------------------------------------------------------------------------------------

   p%LinearizeFlag     = .FALSE.             ! InitInp%LinearizeFlag
   p%Blade%BladeLength = InitInp%TurbineComponents%BladeLength
   p%DtAero            = Interval            ! set the default DT here; may be overwritten later, when we read the input file in AD14_GetInput()
   p%UseDWM            = InitInp%UseDWM

         ! Define parameters here:

   p%WrOptFile   = InitInp%WrSumFile

   p%NumBl   = SIZE( InitInp%TurbineComponents%Blade )
   IF ( p%NumBl < 1 ) THEN
      CALL SetErrStat( ErrID_Fatal,'AeroDyn cannot run without blades in the model.',ErrStat,ErrMess,RoutineName)
      RETURN
   END IF
!bjj: what's the difference between p%NumBl, p%Blade%NB, and InitInp%NumBl?
!MLB: Heck if I know!

         ! Define initial system states here:
   !-------------------------------------------------------------------------------------------------
   ! Read the AeroDyn14 input file and open the output file if requested
   ! bjj: these should perhaps be combined
   !-------------------------------------------------------------------------------------------------
   CALL AD14_GetInput(InitInp, P, x, xd, z, O, y, ErrStatLcl, ErrMessLcl )
      CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
      IF (ErrStat >= AbortErrLev ) RETURN 


      ! allocate variables for aerodyn forces
   p%LinearizeFlag     = .FALSE.

   Interval = p%DtAero
      

   IF ( .NOT. ALLOCATED( o%StoredForces  )) THEN
      CALL AllocAry(O%StoredForces, 3,p%Element%NELM,p%NumBl,'O%StoredForces',ErrStatLcl,ErrMessLcl  )
         CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
   END IF
   IF ( .NOT. ALLOCATED( o%StoredMoments ))  THEN
      CALL AllocAry(O%StoredMoments, 3,p%Element%NELM,p%NumBl,'O%StoredForces',ErrStatLcl,ErrMessLcl  )
         CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
   END IF
     
   IF (.NOT. ALLOCATED(O%Element%W2) ) THEN
      CALL AllocAry(O%Element%W2, p%Element%NELM, p%NumBl,'O%Element%W2',ErrStatLcl,ErrMessLcl  )
         CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
   END IF

   IF (.NOT. ALLOCATED(O%Element%Alpha) ) THEN
      CALL AllocAry(O%Element%Alpha, p%Element%NELM, p%NumBl,'O%Element%Alpha',ErrStatLcl,ErrMessLcl  )
         CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
   END IF
   IF (ErrStat >= AbortErrLev ) RETURN 

   
   P%UnWndOut = -1
   P%UnElem = -1   
   IF ( p%ElemPrn )  THEN
      CALL ElemOpen ( TRIM( InitInp%OutRootName )//'.AD.out', P, O, ErrStat, ErrMess, AD14_Ver )
         CALL SetErrStat( ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName)
         IF (ErrStat >= AbortErrLev ) RETURN 
   END IF
      

   !-------------------------------------------------------------------------------------------------
   ! Calculate the rotor and hub radaii from the input values
   !-------------------------------------------------------------------------------------------------
   HubRadius = DOT_PRODUCT( InitInp%TurbineComponents%Blade(1)%Position(:)        &
                          - InitInp%TurbineComponents%Hub%Position(:),            &
                            InitInp%TurbineComponents%Blade(1)%Orientation(3,:) )

   DO IB = 2,p%NumBl
      TmpVar    = DOT_PRODUCT( InitInp%TurbineComponents%Blade(IB)%Position(:)    &
                             - InitInp%TurbineComponents%Hub%Position(:),         &
                               InitInp%TurbineComponents%Blade(IB)%Orientation(3,:) )
      IF ( ABS( TmpVar - HubRadius ) > 0.001 ) THEN ! within 1 mm
         CALL ProgWarn( ' AeroDyn\AD14_Init() calculated HubRadius is not the same for all '// &
                           'blades. Using value from blade 1.' )
         EXIT
      END IF
   END DO !IB

   TipRadius = InitInp%TurbineComponents%BladeLength + HubRadius

   CosPrecone = ASIN( DOT_PRODUCT( InitInp%TurbineComponents%Blade(1)%Orientation(3,:), &
                                   InitInp%TurbineComponents%Hub%Orientation(1,:) ) )  ! precone angle -- do COS later

   DO IB = 2,p%NumBl
      TmpVar  = ASIN( DOT_PRODUCT( InitInp%TurbineComponents%Blade(IB)%Orientation(3,:), &
                                   InitInp%TurbineComponents%Hub%Orientation(1,:) ) )
      IF ( ABS( TmpVar - CosPrecone ) > 0.009 ) THEN     ! within ~ 1/2 degree
         CALL ProgWarn( ' AeroDyn\AD14_Init() calculated precone angle is not the same for all'// &
                           ' blades. Using value from blade 1.' )
         EXIT
      END IF
   END DO !IBld

   CosPrecone = COS( CosPrecone )

   p%Blade%R = TipRadius * CosPrecone
   RHub = HubRadius * CosPrecone
   p%HubRad = RHub

      ! Check that the AeroDyn input DR and RElm match (use the HubRadius and TipRadius to verify)
      ! before using them to calculate the tip- and hub-loss constants
   CALL CheckRComp( P, x, xd, z, O, y, ErrStat, ErrMess, &
                    InitInp%ADFileName, HubRadius, TipRadius )

   IF ( ErrStat /= ErrID_None ) RETURN

   !-------------------------------------------------------------------------------------------------
   ! Calculate tip-loss constants
   !-------------------------------------------------------------------------------------------------
   DO IElm = 1,p%Element%NElm  ! Loop through all blade elements

      ElemRad = p%Element%RELM(IElm)*CosPrecone

      IF( ElemRad == 0.0 )  THEN  !BJJ: should this be 0.001 (or another small number) instead of exactly 0.0?
         CALL SetErrStat( ErrID_Fatal,'Error calculating tip loss constant for element '//TRIM(Int2LStr(IElm))//&
                          '. Division by zero.',ErrStat,ErrMess,RoutineName)
         
         RETURN
      ELSE
         DTip         = p%Blade%R - ElemRad
         p%Element%TLCNST(IElm) = 0.5 * p%NumBl * DTip / ElemRad
      ENDIF

   ENDDO             ! IElm - all blade elements


   !-------------------------------------------------------------------------------------------------
   ! Calculate hub-loss constants
   !-------------------------------------------------------------------------------------------------
   IF ( RHub > 0.001 )  THEN

      DO Ielm = 1,p%Element%NELM  ! Loop through all blade elements

         ElemRad = p%Element%RELM(Ielm)*CosPrecone  ! Use only the precone angle of blade 1 (assumed very similar to other blades)

         DHub         = ElemRad - RHub
         p%Element%HLCNST(Ielm) = 0.5 * p%NumBl * DHub / RHub

      ENDDO             ! IELM - all blade elements

   ELSE

      p%Element%HLCNST(:) = 0.0

   ENDIF



      !-------------------------------------------------------------------------------------------------
      ! Interpolate the tower diameter at ElastoDyn's tower nodes if we will be computing tower aerodynamics.
      !-------------------------------------------------------------------------------------------------

   IF ( p%TwrProps%CalcTwrAero )  THEN

         !-------------------------------------------------------------------------------------------------
         ! IMPORTANT NOTES:
         !     o  Supposedly, the glue code will not try to do anything with the tower-aero mesh if is is
         !        not created, so the creation is inside the test for CalcTwrAero.
         !     o  The tower properties from AeroDyn's tower file are for heights from the origin (ground or
         !        MSL) to the hub height--not the top of the tower.
         !     o  For now, we are allowing only one set of Cd for the entire tower.
         !     o  InterpIndx is initialize to 1 at compile time.
         !-------------------------------------------------------------------------------------------------


         ! Create the mesh for the tower aerodynamics.

      CALL MeshCreate ( BlankMesh       = u%Twr_InputMarkers      &
                      , IOS             = COMPONENT_INPUT         &
                      , NNodes          = InitInp%NumTwrNodes     &
                      , Orientation     = .TRUE.                  &
                      , TranslationDisp = .TRUE.                  &
                      , TranslationVel  = .TRUE.                  &
                      , ErrStat         = ErrStatLcl              &
                      , ErrMess         = ErrMessLcl              )

      CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
      IF ( ErrStat >= AbortErrLev )  RETURN


         ! Set the positions of the nodes.  MeshCreate() allocated the Position array.

      DO Node = 1,u%Twr_InputMarkers%Nnodes
         CALL MeshPositionNode ( Mesh  = u%Twr_InputMarkers          &
                                ,INode = Node                        &
                                ,Pos   = InitInp%TwrNodeLocs(:,Node) &  
                                ,ErrStat   = ErrStatLcl              &
                                ,ErrMess   = ErrMessLcl              )
         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN

      END DO         


         ! Construct the tower with Line-2 elements.

      DO Elem=1,u%Twr_InputMarkers%Nnodes-1

         CALL MeshConstructElement ( Mesh     = u%Twr_InputMarkers &
                                   , Xelement = ELEMENT_LINE2      &
                                   , P1       = Elem               &
                                   , P2       = Elem+1             &
                                   , ErrStat  = ErrStatLcl         &
                                   , ErrMess  = ErrMessLcl         )

         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN

      ENDDO


         ! Commit the mesh to the funny farm.

      CALL MeshCommit ( u%Twr_InputMarkers, ErrStatLcl, ErrMessLcl )
         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN


         ! Copy the input mesh to create the output mesh.  Does

      CALL MeshCopy ( SrcMesh  = u%Twr_InputMarkers &
                    , DestMesh = y%Twr_OutputLoads  &
                    , CtrlCode = MESH_SIBLING       &
                    , Force    = .TRUE.             &
                    , ErrStat  = ErrStatLcl         &
                    , ErrMess  = ErrMessLcl         )

         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN


         ! Check to ensure that the user did not specify more than one set of Cd(Re) tables.  Temporary restriction.

      IF ( p%TwrProps%NTwrCD /= 1 )  THEN
         CALL SetErrStat(ErrID_Fatal,'You must have one and only one set of drag coefficients for the AeroDyn tower file.',ErrStat,ErrMess,RoutineName )
         RETURN
      END IF


         ! Build the TwrNodeWidth array.

      p%TwrProps%NumTwrNodes = InitInp%NumTwrNodes

      IF (.NOT. ALLOCATED( p%TwrProps%TwrNodeWidth ) ) THEN
         CALL AllocAry( p%TwrProps%TwrNodeWidth, p%TwrProps%NumTwrNodes, "array for tower widths at ED node locations", ErrStatLcl, ErrMessLcl )
         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN
      END IF

      DO Node=1,p%TwrProps%NumTwrNodes

         TwrNodeHt = InitInp%TwrNodeLocs(3,Node)/p%Rotor%HH

         p%TwrProps%TwrNodeWidth(Node) = InterpStp( TwrNodeHt, p%TwrProps%TwrHtFr, p%TwrProps%TwrWid, InterpIndx, p%TwrProps%NTwrHT )

      END DO ! Node

   ELSE
      u%Twr_InputMarkers%Nnodes = 0
      y%Twr_OutputLoads%Nnodes  = 0
   END IF ! ( p%TwrProps%CalcTwrAero )


   !-------------------------------------------------------------------------------------------------
   ! Write the summary (opt) file, then close it
   !-------------------------------------------------------------------------------------------------

   IF (p%WrOptFile) THEN

      CALL ADOut(InitInp, P, O, AD14_Ver, TRIM(InitInp%OutRootName)//'.AD.sum', ErrStatLcl, ErrMessLcl )
         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN

   ENDIF


   !-------------------------------------------------------------------------------------------------
   ! Initialize the inputs from the wind inflow module
   !-------------------------------------------------------------------------------------------------
   CALL AllocAry( u%InflowVelocity, 3, p%Element%NElm*p%NumBl + u%Twr_InputMarkers%Nnodes, 'u%InflowVelocity', ErrStatLcl, ErrMessLcl )
         CALL SetErrStat(ErrStatLcl,ErrMessLcl,ErrStat,ErrMess,RoutineName )
         IF ( ErrStat >= AbortErrLev )  RETURN
   u%InflowVelocity = 0.0_ReKi
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Calling the DWM
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   IF ( p%UseDWM ) THEN   
      ! InitInp%DWM_InitInputs%InputFileName is already set in FAST
         
         ! bjj: all this stuff should be put in DWM_Init.....>
      p%DWM_Params%RR              = p%Blade%R
      p%DWM_Params%BNum            = p%NumBl
      p%DWM_Params%ElementNum      = O%ElOut%NumElOut  !bjj: NumElOut is the number of elements to be printed in an output file. I really think you want the number of blade elements. I guess we should check that NumElOut is the same as p%Element%NElm
      p%DWM_Params%air_density     = p%Wind%Rho
   
      IF (.NOT. ALLOCATED(o%DWM_otherstates%Nforce  )) ALLOCATE ( o%DWM_otherstates%Nforce(  p%Element%NElm,p%NumBl),STAT=ErrStatLcl);CALL SetErrStat(ErrStatLcl, 'Error allocating DWM Nforce array', ErrStat,ErrMess,RoutineName )
      IF (.NOT. ALLOCATED(o%DWM_otherstates%blade_dr)) ALLOCATE ( o%DWM_otherstates%blade_dr(p%Element%NElm),        STAT=ErrStatLcl);CALL SetErrStat(ErrStatLcl, 'Error allocating DWM blade_dr array', ErrStat,ErrMess,RoutineName )
      IF (.NOT. ALLOCATED(p%DWM_Params%ElementRad   )) ALLOCATE ( p%DWM_Params%ElementRad(   p%Element%NElm),        STAT=ErrStatLcl);CALL SetErrStat(ErrStatLcl, 'Error allocating DWM ElementRad array', ErrStat,ErrMess,RoutineName )
     
      o%DWM_otherstates%blade_dr(:) = p%Blade%DR(:)
      p%DWM_Params%ElementRad(:)    = p%Element%RELM(:)   
   
      CALL DWM_Init( InitInp%DWM_InitInputs, O%DWM_Inputs, p%DWM_Params, x%DWM_ContStates, xd%DWM_DiscStates, z%DWM_ConstrStates, & 
                     O%DWM_OtherStates, y%DWM_Outputs, Interval, InitOut%DWM_InitOutput, ErrStatLcl, ErrMessLcl )
   

      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN
      
   END IF !UseDWM
   
   !-------------------------------------------------------------------------------------------------
   ! Turn off dynamic inflow for wind less than 8 m/s (per DJL: 8 m/s is really just an empirical guess)
   ! DJL: Comment out this code when using new proposed GDW check in ELEMFRC
   ! BJJ: FIX THIS!!!!
   !-------------------------------------------------------------------------------------------------

   ! BJJ: can't put this here b/c we need InitInp%MWS from InflowWind
   !IF (p%DynInfl) THEN
   !            
   !   IF ( InitInp%MWS  < 8.0 ) THEN
   !      p%DynInfl = .FALSE.
   !      CALL SetErrStat(ErrID_Info,'Estimated average inflow wind speed is less than 8 m/s. Dynamic Inflow will be turned off.',ErrStat,ErrMess,RoutineName )
   !      IF ( ErrStat >= AbortErrLev )  RETURN      
   !   END IF
   !
   !ENDIF   
   

   !-------------------------------------------------------------------------------------------------
   ! Set initial guesses for inputs:
   !-------------------------------------------------------------------------------------------------
   
   !..........
   ! u%TurbineComponents
   !..........

   CALL AD14_CopyAeroConfig( InitInp%TurbineComponents, u%TurbineComponents, MESH_NEWCOPY, ErrStatLcl, ErrMessLcl )
      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN
   
   !..........
   ! u%InputMarkers (blade meshes):
   !..........
   
   ALLOCATE( u%InputMarkers(p%NumBl), STAT=ErrStatLcl )
   IF (ErrStatLcl /= 0 ) THEN
      CALL SetErrStat ( ErrID_Fatal, 'Could not allocate u%InputMarkers (meshes)', ErrStat,ErrMess,RoutineName )
      RETURN
   END IF


   DO IB = 1, p%NumBl
      CALL MeshCreate( BlankMesh      = u%InputMarkers(IB)    &
                     ,IOS            = COMPONENT_INPUT        &
                     ,NNodes         = p%Element%NELM         &
                     ,Orientation    = .TRUE.                 &
                     ,TranslationVel = .TRUE.                 &
                     ,TranslationAcc = .TRUE.                 &  !bjj: added for MHK turbines
                     ,RotationVel    = .TRUE.                 &
                     ,nScalars       = 2                      &  ! scalar 1 is W, scalar 2 is Alpha
                     ,ErrStat        = ErrStatLcl             &
                     ,ErrMess        = ErrMessLcl             )

      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
      IF (ErrStat >= AbortErrLev) RETURN
      
      ! create the elements
      DO IE = 1, p%Element%NELM-1 ! construct the blades into Line2 elements
         CALL MeshConstructElement ( Mesh = u%InputMarkers(IB)    &
                                  ,Xelement = ELEMENT_LINE2       &
                                  ,P1       = IE                  &
                                  ,P2       = IE+1                &
                                  ,ErrStat  = ErrStatLcl          &
                                  ,ErrMess  = ErrMessLcl          )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
         
      ENDDO
     
      ! position/orient the nodes
      DO IE = 1, p%Element%NELM
         TmpPos(1) = 0.
         TmpPos(2) = 0.
         TmpPos(3) = p%Element%Relm(IE) - HubRadius
         CALL MeshPositionNode ( Mesh = u%InputMarkers(IB)              &
                                 ,INode = IE                            &
                                 ,Pos= TmpPos                           &  ! this info comes from FAST (not yet)
                                 ,ErrStat   = ErrStatLcl                &
                                 ,ErrMess   = ErrMessLcl                )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN

         ! RELATIVE ORIENTATION OF BLADE ELEMENTS
         u%InputMarkers(IB)%Orientation(1,1,IE) = COS( P%Element%TWIST(IE) )
         u%InputMarkers(IB)%Orientation(2,1,IE) = SIN( P%Element%TWIST(IE) )
         u%InputMarkers(IB)%Orientation(3,1,IE) = SIN( P%Element%TWIST(IE) )
         u%InputMarkers(IB)%Orientation(1,2,IE) = -1. * u%InputMarkers(IB)%Orientation(2,1,IE)
         u%InputMarkers(IB)%Orientation(2,2,IE) =       u%InputMarkers(IB)%Orientation(1,1,IE)
         u%InputMarkers(IB)%Orientation(3,2,IE) = 0.0
         u%InputMarkers(IB)%Orientation(1,3,IE) = 0.0
         u%InputMarkers(IB)%Orientation(2,3,IE) = 0.0
         u%InputMarkers(IB)%Orientation(3,3,IE) = 1.0
      ENDDO
     
       CALL MeshCommit ( Mesh = u%InputMarkers(IB)    &
                        ,ErrStat  = ErrStatLcl        &
                        ,ErrMess  = ErrMessLcl        )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
          
   ENDDO
         

   !..........
   ! u%Twr_InputMarkers (tower meshes):
   !..........
     
   !bjj: done above in section for IF Tower Loads is on
   

   !..........
   ! u%MulTabLoc:
   !..........
   
   IF (.NOT. ALLOCATED(u%MulTabLoc)) THEN
      ALLOCATE( u%MulTabLoc(p%Element%NELM, p%NumBl), STAT = ErrStatLcl )
      IF (ErrStatLcl /= 0) THEN
         CALL SetErrStat ( ErrID_Fatal, 'Could not allocate u%MulTabLoc', ErrStat,ErrMess,RoutineName )
         RETURN
      END IF
   END IF

   u%MulTabLoc(:,:) = 0.0
   
   
   !-------------------------------------------------------------------------------------------------
   ! Allocate space for outputs and set up output meshes:
   !-------------------------------------------------------------------------------------------------
   
   !..........
   ! y%OutputLoads (blade meshes):
   !..........
   
   
   ALLOCATE( y%OutputLoads(p%NumBl), STAT = ErrStatLcl )
      IF (ErrStatLcl /= 0) THEN
         CALL SetErrStat ( ErrID_Fatal, 'Could not allocate y%OutputLoads (meshes)', ErrStat,ErrMess,RoutineName )
         RETURN
      END IF
   
   DO IB = 1, p%NumBl

       CALL MeshCopy ( SrcMesh  = u%InputMarkers(IB)  &
                      ,DestMesh = y%OutputLoads(IB)   &
                      ,CtrlCode = MESH_SIBLING        &
                      ,Force    = .TRUE.              &
                      ,Moment   = .TRUE.              &
                      ,ErrStat  = ErrStatLcl          &
                      ,ErrMess  = ErrMessLcl          )
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,RoutineName )
         IF (ErrStat >= AbortErrLev) RETURN
   ENDDO
   
   
   !..........
   ! y%Twr_OutputLoads (tower meshes):
   !..........

   !bjj: done above in section for IF Tower Loads is on
   
     
   !-------------------------------------------------------------------------------------------------
   ! Initialize AeroDyn variables not initialized elsewhere (except in module initialization)
   ! and return
   !-------------------------------------------------------------------------------------------------
   o%InducedVel%SumInfl  = 0.0_ReKi
   o%Rotor%AvgInfl       = 0.0_ReKi
   o%OldTime             = 0.0_DbKi
   O%SuperSonic          = .FALSE.   
   o%NoLoadsCalculated   = .TRUE.

   p%TwoPiNB     = TwoPi / REAL( p%NumBl, ReKi )

   p%Initialized = .TRUE.
   
        
   
   DO ie = 1, maxInfl
      p%DynInflow%xMinv(ie) = PIBY2 / hfunc(MRvector(ie), NJvector(ie))   !bjj: this is really just a Fortran parameter, too.
   END DO !ie   

   RETURN

END SUBROUTINE AD14_Init

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD14_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMess )
! This routine is called at the end of the simulation.
!..................................................................................................................................
      USE DWM_Types
      USE DWM

      TYPE(AD14_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(AD14_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
      TYPE(AD14_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(AD14_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(AD14_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(AD14_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(AD14_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMess     ! Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMess  = ""


         ! Place any last minute operations or calculations here:
         
      IF (p%UseDWM ) THEN
         !----- Call the DWM ------- 
      
         CALL DWM_End( OtherState%DWM_Inputs, p%DWM_Params, x%DWM_ContStates, xd%DWM_DiscStates, z%DWM_ConstrStates, &
                                           OtherState%DWM_OtherStates, y%DWM_Outputs, ErrStat, ErrMess )
      END IF ! UseDWM
      
      !--------------------------


         ! Close files here:

      ! AD14_IOParams
   IF (P%UnEc > 0)    CLOSE(P%UnEc) ! not currently used

   IF (P%UnWndOut > 0) CLOSE(P%UnWndOut)
   IF (P%UnElem   > 0) CLOSE(P%UnElem)

         ! Destroy the input data:

      CALL AD14_DestroyInput( u, ErrStat, ErrMess )


         ! Destroy the parameter data:

      CALL AD14_DestroyParam( p, ErrStat, ErrMess )


         ! Destroy the state data:

      CALL AD14_DestroyContState(   x,           ErrStat, ErrMess )
      CALL AD14_DestroyDiscState(   xd,          ErrStat, ErrMess )
      CALL AD14_DestroyConstrState( z,           ErrStat, ErrMess )
      CALL AD14_DestroyOtherState(  OtherState,  ErrStat, ErrMess )


         ! Destroy the output data:

      CALL AD14_DestroyOutput( y, ErrStat, ErrMess )




END SUBROUTINE AD14_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD14_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMess )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input Time; Continuous and discrete states are updated for Time + Interval
!..................................................................................................................................

      REAL(DbKi),                           INTENT(IN   ) :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                       INTENT(IN   ) :: n           ! Current simulation time step n = 0,1,...
      TYPE(AD14_InputType),                 INTENT(INOUT) :: u(:)        ! Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
      REAL(DbKi),                           INTENT(IN   ) :: utimes(:)   ! Times associated with u(:), in seconds
      TYPE(AD14_ParameterType),             INTENT(IN   ) :: p           ! Parameters
      TYPE(AD14_ContinuousStateType),       INTENT(INOUT) :: x           ! Input: Continuous states at t;
                                                                       !   Output: Continuous states at t + Interval
      TYPE(AD14_DiscreteStateType),         INTENT(INOUT) :: xd          ! Input: Discrete states at t;
                                                                       !   Output: Discrete states at t  + Interval
      TYPE(AD14_ConstraintStateType),       INTENT(INOUT) :: z           ! Input: Initial guess of constraint states at t;
                                                                       !   Output: Constraint states at t
      TYPE(AD14_OtherStateType),            INTENT(INOUT) :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                       INTENT(  OUT) :: ErrStat     ! Error status of the operation
      CHARACTER(*),                         INTENT(  OUT) :: ErrMess     ! Error message if ErrStat /= ErrID_None

         ! Local variables

      TYPE(AD14_ContinuousStateType)                 :: dxdt        ! Continuous state derivatives at Time
      TYPE(AD14_ConstraintStateType)                 :: z_Residual  ! Residual of the constraint state equations (Z)

!      INTEGER(IntKi)                                    :: ErrStat2    ! Error status of the operation (occurs after initial error)
!      CHARACTER(ErrMsgLen)                              :: ErrMess2     ! Error message if ErrStat2 /= ErrID_None

         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMess  = ""





END SUBROUTINE AD14_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD14_CalcOutput( Time, u, p, x, xd, z, O, y, ErrStat, ErrMess )
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................
   
      USE               AeroGenSubs,   ONLY: ElemOut
      USE               DWM_Types
      USE               DWM

      REAL(DbKi),                     INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(AD14_InputType),           INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(AD14_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(AD14_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(AD14_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(AD14_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(AD14_OtherStateType),      INTENT(INOUT)  :: O!therState ! Other/optimization states
      TYPE(AD14_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at Time (Input only so that mesh con-
                                                                       !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMess     ! Error message if ErrStat /= ErrID_None


      ! Local variables
   REAL(DbKi), PARAMETER      :: OnePlusEpsilon = 1 + EPSILON(Time)

   REAL(ReKi)                 :: VNElement
   REAL(ReKi)                 :: VelNormalToRotor2
   REAL(ReKi)                 :: VNWind
   REAL(ReKi)                 :: VTTotal
   REAL(ReKi)                 :: DFN
   REAL(ReKi)                 :: DFT
   REAL(ReKi)                 :: PMA
   REAL(ReKi)                 :: SPitch                     ! sine of PitNow
   REAL(ReKi)                 :: CPitch                     ! cosine of PitNow

   REAL(ReKi)                 :: AvgVelNacelleRotorFurlYaw
   REAL(ReKi)                 :: AvgVelTowerBaseNacelleYaw
   REAL(ReKi)                 :: AvgVelTowerBaseYaw
   REAL(ReKi)                 :: AzimuthAngle
   REAL(ReKi)                 :: rNacelleHub   (2)
   REAL(ReKi)                 :: rLocal
   REAL(ReKi)                 :: rRotorFurlHub (2)
   REAL(ReKi)                 :: rTowerBaseHub (2)

   REAL(ReKi)                 :: tmpVector     (3)
   REAL(ReKi)                 :: VelocityVec   (3)

   INTEGER                    :: ErrStatLcL        ! Error status returned by called routines.
   INTEGER                    :: IBlade
   INTEGER                    :: IElement
   INTEGER                    :: Node              ! Node index.

   INTEGER                    :: I
   CHARACTER(ErrMsgLen)       :: ErrMessLcl          ! Error message returned by called routines.



   ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMess  = ""

   !-------------------------------------------------------------------------------------------------
   ! Check that the module has been initialized.
   !-------------------------------------------------------------------------------------------------
   IF ( .NOT. p%Initialized ) THEN
      ErrStat = ErrID_Fatal
      ErrMess = 'AeroDyn must be initialized before trying to calculate aerodynamic loads.'
      RETURN
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Determine if loads should be recalculated or just returned
   !-------------------------------------------------------------------------------------------------
      ! NOTE: Time is scaled by OnePlusEps to ensure that loads are calculated at every
   !       time step when DTAero = DT, even in the presence of numerical precision errors.

   IF ( o%NoLoadsCalculated .OR. ( Time*OnePlusEpsilon - o%OldTime ) >= p%DTAERO )  THEN
         ! It's time to update the aero forces

         ! First we reset the DTAERO parameters for next time
      o%DT      = Time - o%OldTime     !bjj: DT = 0 on first step,
                                       !but the subroutines that use DT check for NoLoadsCalculated (or time > 0)
      o%OldTime = Time

   ELSE IF ( .NOT. p%LinearizeFlag ) THEN

         ! Return the previously-calculated loads

!      CurrentOutputs = ADCurrentLoads

      DO IBlade=1,p%NumBl
       DO IElement=1,p%Element%Nelm
         y%OutputLoads(IBlade)%Force(:,IElement)  = o%StoredForces(:,IElement,IBlade)
         y%OutputLoads(IBlade)%Moment(:,IElement) = o%StoredMoments(:,IElement,IBlade)
       ENDDO
      ENDDO

      IF ( O%FirstWarn ) THEN
         CALL SetErrStat ( ErrID_Warn, 'AeroDyn was designed for an explicit-loose coupling scheme. '//&
            'Using last calculated values from AeroDyn on all subsequent calls until time is advanced. '//&
            'Warning will not be displayed again.', ErrStat,ErrMess,'AD14_CalcOutput' )
         O%FirstWarn = .FALSE.       
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
         END IF
      END IF
            
      RETURN

   ENDIF

   
   !-------------------------------------------------------------------------------------------------
   ! Calculate the forces and moments for the blade: SUBROUTINE AeroFrcIntrface( FirstLoop, JElemt, DFN, DFT, PMA )
   !-------------------------------------------------------------------------------------------------

      ! calculate rotor speed
      ! note: Subtracting the RotorFurl rotational velocity for REVS is needed to get the
      ! same answers as before v13.00.00. RotorFurl shouldn't be needed.

   o%Rotor%REVS = ABS( DOT_PRODUCT( u%TurbineComponents%Hub%RotationVel(:) - u%TurbineComponents%RotorFurl%RotationVel(:), &
                                    u%TurbineComponents%Hub%Orientation(1,:) ) )


      ! calculate yaw angle
      ! note: YawAng should use the Hub instead of the RotorFurl, but it is calculated this way to
      ! get the same answers as previous version.
   o%Rotor%YawAng = ATAN2( -1.*u%TurbineComponents%RotorFurl%Orientation(1,2), u%TurbineComponents%RotorFurl%Orientation(1,1) )
   o%Rotor%SYaw   = SIN( o%Rotor%YawAng )
   o%Rotor%CYaw   = COS( o%Rotor%YawAng )

      ! tilt angle
      ! note: tilt angle should use the Hub instead of RotorFurl, but it needs hub to get the same
      ! answers as the version before v13.00.00

   o%Rotor%Tilt = ATAN2( u%TurbineComponents%RotorFurl%Orientation(1,3), &
                         SQRT( u%TurbineComponents%RotorFurl%Orientation(1,1)**2 + &
                         u%TurbineComponents%RotorFurl%Orientation(1,2)**2 ) )

   o%Rotor%CTilt     = COS( o%Rotor%Tilt )
   o%Rotor%STilt     = SIN( o%Rotor%Tilt )


      ! HubVDue2Yaw - yaw velocity due solely to yaw

  AvgVelNacelleRotorFurlYaw = u%TurbineComponents%RotorFurl%RotationVel(3) - u%TurbineComponents%Nacelle%RotationVel(3)
  AvgVelTowerBaseNacelleYaw = u%TurbineComponents%Nacelle%RotationVel(3)   - u%TurbineComponents%Tower%RotationVel(3)
  AvgVelTowerBaseYaw        = u%TurbineComponents%Tower%RotationVel(3)

  rRotorFurlHub(1:2)        = u%TurbineComponents%Hub%Position(1:2) - u%TurbineComponents%RotorFurl%Position(1:2)
  rNacelleHub(1:2)          = u%TurbineComponents%Hub%Position(1:2) - u%TurbineComponents%Nacelle%Position(1:2)
  rTowerBaseHub(1:2)        = u%TurbineComponents%Hub%Position(1:2) - u%TurbineComponents%Tower%Position(1:2)

  o%Rotor%YawVel =   ( AvgVelNacelleRotorFurlYaw * rRotorFurlHub(2) + AvgVelTowerBaseNacelleYaw * rNacelleHub(2) &
                         + AvgVelTowerBaseYaw * rTowerBaseHub(2) ) * o%Rotor%SYaw &
                  - ( AvgVelNacelleRotorFurlYaw * rRotorFurlHub(1) + AvgVelTowerBaseNacelleYaw * rNacelleHub(1) &
                         + AvgVelTowerBaseYaw * rTowerBaseHub(1) ) * o%Rotor%CYaw


   !.................................................................................................
   ! start of NewTime routine
   !.................................................................................................

   o%Rotor%AvgInfl = o%InducedVel%SumInfl * 2.0 / (p%Blade%R*p%Blade%R*p%NumBl)  ! Average inflow from the previous time step
   o%InducedVel%SumInfl = 0.0   ! reset to sum for the current time step

   CALL DiskVel(Time, P, O, u%AvgInfVel, ErrStatLcl, ErrMessLcl)  ! Get a sort of "Average velocity" - sets a bunch of stored variables...
      CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput/DiskVel' )
      IF (ErrStat >= AbortErrLev) THEN
         CALL CleanUp()
         RETURN
      END IF
         
   IF ( P%DStall ) CALL BedUpdate( O )   ! that's an 'O' as in 'OtherState'

   ! Enter the dynamic inflow routines here

   IF ( p%Wake )  THEN
      CALL Inflow(Time, P, O, ErrStatLcl, ErrMessLcl)
         CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput/Inflow' )
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp()
            RETURN
         END IF   
   END IF      
   
       !bjj: perhaps we should send NoLoadsCalculated to initialize dynamic inflow [subroutine Infinit()]
       !bjj: instead of the check that time > 0...?

   !.................................................................................................
   ! end of NewTime routine
   !.................................................................................................

   Node = 0
   DO IBlade = 1,p%NumBl

         ! calculate the azimuth angle ( we add pi because AeroDyn defines 0 as pointing downward)
         ! note: the equation below should use TurbineComponents%Blade markers, but this is used to get the
         ! same answers as the previous version (before v13.00.00)

      AzimuthAngle = ATAN2( -1.*DOT_PRODUCT( u%TurbineComponents%Hub%Orientation(3,:),         &
                                             u%TurbineComponents%RotorFurl%Orientation(2,:) ), &
                                DOT_PRODUCT( u%TurbineComponents%Hub%Orientation(3,:),         &
                                             u%TurbineComponents%RotorFurl%Orientation(3,:) )  ) + pi + (IBlade - 1)*p%TwoPiNB



      DO IElement = 1,p%Element%NElm

            ! calculate element pitch

         o%Element%PitNow    = -1.*ATAN2( -1.*DOT_PRODUCT( u%TurbineComponents%Blade(IBlade)%Orientation(1,:),    &
                                                           u%InputMarkers(IBlade)%Orientation(2,:,IElement) ) , &
                                              DOT_PRODUCT( u%TurbineComponents%Blade(IBlade)%Orientation(1,:),    &
                                                           u%InputMarkers(IBlade)%Orientation(1,:,IElement) )   )

         SPitch    = SIN( o%Element%PitNow )
         CPitch    = COS( o%Element%PitNow )


            ! calculate distance between hub and element

         tmpVector = u%InputMarkers(IBlade)%Position(:,IElement) - u%TurbineComponents%Hub%Position(:)
         rLocal = SQRT(   DOT_PRODUCT( tmpVector, u%TurbineComponents%Hub%Orientation(2,:) )**2  &
                        + DOT_PRODUCT( tmpVector, u%TurbineComponents%Hub%Orientation(3,:) )**2  )

            ! determine if MulTabLoc should be set.  
     
         IF (.not. p%Reynolds) O%AirFoil%MulTabLoc = u%MulTabLoc(IElement,IBlade)
         
         !-------------------------------------------------------------------------------------------
         ! Get wind velocity components; calculate velocity normal to the rotor squared
         ! Save variables for printing in a file later;
         !-------------------------------------------------------------------------------------------
         Node = Node + 1
         VelocityVec(:)    = AD_WindVelocityWithDisturbance( Time, u, p, x, xd, z, O, y, ErrStatLcl, ErrMessLcl, &
                                                             u%InputMarkers(IBlade)%Position(:,IElement), &
                                                             u%InflowVelocity(:,Node) )
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput' )
            IF (ErrStat >= AbortErrLev) THEN
               CALL CleanUp()
               RETURN
            END IF
         
         !-------------------------------------------------------------------------------------------
         ! DWM wind input update phase 1
         !-------------------------------------------------------------------------------------------
         IF (p%UseDWM) THEN
            !bjj: FIX THIS!!!!         
            !bjj: where do p%DWM_Params%RTPD%SimulationOrder_index and p%DWM_Params%RTPD%upwindturbine_number get set?
         
            IF ( p%DWM_Params%RTPD%SimulationOrder_index > 1) THEN
               IF(  p%DWM_Params%RTPD%upwindturbine_number /= 0 ) THEN
                       
                  o%DWM_otherstates%position_y = u%InputMarkers(IBlade)%Position(2,IElement)
                                         
                  o%DWM_otherstates%position_z = u%InputMarkers(IBlade)%Position(3,IElement)
              
                  o%DWM_otherstates%velocity_wake_mean = 1
              
                  DO I = 1,p%DWM_Params%RTPD%upwindturbine_number
                     o%DWM_otherstates%DWM_tb%Aerodyn_turbine_num = I
                 
                     CALL   DWM_phase1( Time, O%DWM_Inputs, p%DWM_Params, x%DWM_contstates, xd%DWM_discstates, z%DWM_constrstates, &
                                           o%DWM_otherstates, y%DWM_outputs, ErrStatLcl, ErrMessLcl )
                 
                     CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput/DWM_phase1' )
                     IF (ErrStat >= AbortErrLev) THEN
                        CALL CleanUp()
                        RETURN
                     END IF   
                                                           
                     o%DWM_otherstates%velocity_wake_mean = (1-((1-o%DWM_otherstates%velocity_wake_mean)**2 + (1-o%DWM_otherstates%shifted_velocity_aerodyn)**2)**0.5)
                  END DO
              
                  o%DWM_otherstates%velocity_wake_mean    = o%DWM_otherstates%velocity_wake_mean * p%DWM_Params%Wind_file_Mean_u
              
                  VelocityVec(1) = (VelocityVec(1) - p%DWM_Params%Wind_file_Mean_u)*(O%DWM_Inputs%Upwind_result%upwind_small_TI(1)/p%DWM_Params%TI_amb) &
                                  + o%DWM_otherstates%velocity_wake_mean
              
               END IF
            END IF
                                     
           !------------------------DWM PHASE 2-----------------------------------------------
            IF (Time > 50.00 ) THEN
               o%DWM_otherstates%U_velocity           = VelocityVec(1)
               o%DWM_otherstates%V_velocity           = VelocityVec(2)
               o%DWM_otherstates%NacYaw               = o%Rotor%YawAng 
               o%DWM_otherstates%DWM_tb%Blade_index   = IBlade
               o%DWM_otherstates%DWM_tb%Element_index = IElement    

               CALL   DWM_phase2( Time, O%DWM_Inputs, p%DWM_Params, x%DWM_contstates, xd%DWM_discstates, z%DWM_constrstates, &
                                           o%DWM_otherstates, y%DWM_outputs, ErrStatLcl, ErrMessLcl )
            
                     CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput/DWM_phase1' )
                     IF (ErrStat >= AbortErrLev) THEN
                        CALL CleanUp()
                        RETURN
                     END IF   
            
         
               !CALL CalVelScale(VelocityVec(1),VelocityVec(2),O%o%DWM_otherstatesutputType,DWM_ConstraintStateType)
         
               !CALL turbine_average_velocity( VelocityVec(1), IBlade, IElement, O%o%DWM_otherstatesutputType,AD14_ParameterType,DWM_ConstraintStateType)
            END IF
         END IF ! UseDWM
            
         !-----------------------------------------------------------------------------------------------------------------------
         
         
         
         
         VelNormalToRotor2 = ( VelocityVec(3) * o%Rotor%STilt + (VelocityVec(1) * o%Rotor%CYaw               &
                             - VelocityVec(2) * o%Rotor%SYaw) * o%Rotor%CTilt )**2

         !-------------------------------------------------------------------------------------------
         ! reproduce GetVNVT routine:
         !-------------------------------------------------------------------------------------------
         tmpVector =  -1.*SPitch*u%InputMarkers(IBlade)%Orientation(1,:,IElement) &
                        + CPitch*u%InputMarkers(IBlade)%Orientation(2,:,IElement)
         VTTotal   =     DOT_PRODUCT( tmpVector, VelocityVec - u%InputMarkers(IBlade)%TranslationVel(:,IElement)  )

         tmpVector =     CPitch*u%InputMarkers(IBlade)%Orientation(1,:,IElement) &
                       + SPitch*u%InputMarkers(IBlade)%Orientation(2,:,IElement)
         VNWind    =     DOT_PRODUCT( tmpVector, VelocityVec )
         VNElement = -1.*DOT_PRODUCT( tmpVector, u%InputMarkers(IBlade)%TranslationVel(:,IElement ) )

         !-------------------------------------------------------------------------------------------
         ! Get blade element forces and induced velocity
         !-------------------------------------------------------------------------------------------
         CALL ELEMFRC( p, O, ErrStatLcl, ErrMessLcl,                             &
                       AzimuthAngle, rLocal, IElement, IBlade, VelNormalToRotor2, VTTotal, VNWind, &
                       VNElement, DFN, DFT, PMA, o%NoLoadsCalculated )
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput' )
            IF (ErrStat >= AbortErrLev) THEN
               CALL CleanUp()
               RETURN
            END IF
            
         !-------------------------------------------------------------------------------------------
         ! Set up dynamic inflow parameters
         !-------------------------------------------------------------------------------------------
         IF ( p%DynInfl .OR. O%DynInit ) THEN
            CALL GetRM (P, O, ErrStatLcl, ErrMessLcl, &
                        rLocal, DFN, DFT, AzimuthAngle, IElement, IBlade)
               CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput' )
               IF (ErrStat >= AbortErrLev) THEN
                  CALL CleanUp()
                  RETURN
               END IF
            
         ENDIF
         
         IF (p%UseDWM) THEN
             o%DWM_otherstates%Nforce(IElement,IBlade) = DFN              ! 12.4.2014 add by yh
         END IF ! UseDWM

         o%StoredForces(1,IElement,IBlade)  = ( DFN*CPitch + DFT*SPitch ) / p%Blade%DR(IElement)
         o%StoredForces(2,IElement,IBlade)  = ( DFN*SPitch - DFT*CPitch ) / p%Blade%DR(IElement)
         o%StoredForces(3,IElement,IBlade)  = 0.0

         o%StoredMoments(1,IElement,IBlade)  = 0.0
         o%StoredMoments(2,IElement,IBlade)  = 0.0
         o%StoredMoments(3,IElement,IBlade)  = PMA / p%Blade%DR(IElement)

!      DO IBlade=1,p%NumBl
!       DO IElement=1,p%Element%Nelm
!         y%OutputLoads(IBlade)%Force(:,IElement)  = o%StoredForces(:,IElement,IBlade)
!         y%OutputLoads(IBlade)%Moment(:,IElement) = o%StoredMoments(:,IElement,IBlade)
!       ENDDO
!!      ENDDO

            ! save velocities for output, if requested

         IF ( O%ElOut%WndElPrList(IElement) > 0 ) THEN
            O%ElOut%SaveVX( O%ElOut%WndElPrList(IElement), IBlade ) = VelocityVec(1)
            O%ElOut%SaveVY( O%ElOut%WndElPrList(IElement), IBlade ) = VelocityVec(2)
            O%ElOut%SaveVZ( O%ElOut%WndElPrList(IElement), IBlade ) = VelocityVec(3)
         ENDIF


      END DO !IElement

      IF ( IBlade == 1 .AND. p%ElemPrn ) THEN
         O%ElOut%VXSAV  = VelocityVec(1)
         O%ElOut%VYSAV  = VelocityVec(2)
         O%ElOut%VZSAV  = VelocityVec(3)
      ENDIF


   END DO !IBlade

   O%NoLoadsCalculated = .FALSE.


   DO IBlade=1,p%NumBl
     DO IElement=1,p%Element%Nelm
       y%OutputLoads(IBlade)%Force(:,IElement)  = o%StoredForces(:,IElement,IBlade)
       y%OutputLoads(IBlade)%Moment(:,IElement) = o%StoredMoments(:,IElement,IBlade)
     ENDDO
   ENDDO
   
   
   !------------------------DWM PHASE 3-----------------------------------------------
   IF (p%UseDWM) THEN
   
      IF (Time > 50.00 ) THEN !BJJ: why is 50 hard-coded here and above???
            
         !o%DWM_otherstates%Nforce(:,:)    = o%DWM_otherstates%DFN_DWM(:,:) 
         CALL   DWM_phase3( Time, O%DWM_Inputs, p%DWM_Params, x%DWM_contstates, xd%DWM_discstates, z%DWM_constrstates, &
                                o%DWM_otherstates, y%DWM_outputs, ErrStatLcl, ErrMessLcl )
    
            CALL SetErrStat ( ErrStatLcl, ErrMessLcl, ErrStat,ErrMess,'AD14_CalcOutput/DWM_phase3' )
            IF (ErrStat >= AbortErrLev) THEN
               CALL CleanUp()
               RETURN
            END IF   
      
         !CALL filter_average_induction_factor( AD14_ParameterType, DWM_ConstraintStateType, O%o%DWM_otherstatesutputType )
      END IF
   END IF !UseDWM
   
   !-----------------------------------------------------------------------------------


   
       
         

      ! Loop through all the tower nodes to calculate the aerodynamic loads on the tower if aerodynamics were requested.

   IF ( p%TwrProps%CalcTwrAero )  THEN

      DO Node=1,u%Twr_InputMarkers%Nnodes


            ! Calculate the aerodynamic load on this tower node: TwrAeroLoads ( p, Node, NodeDCMGbl, NodeVelGbl, NodeWindVelGbl, NodeFrcGbl )

         CALL TwrAeroLoads ( p, Node, u%Twr_InputMarkers%Orientation(:,:,Node), u%Twr_InputMarkers%TranslationVel(:,Node) &
                           , u%InflowVelocity(:,Node+p%NumBl*p%Element%NElm), y%Twr_OutputLoads%Force(:,Node) )

      END DO ! Node

   END IF ! ( p%TwrProps%CalcTwrAero )

   !................................................................................................
      
   
   CALL ElemOut(time, P, O )
   
   CALL CleanUp (  )

   RETURN

   !=======================================================================
   CONTAINS
   !=======================================================================
   SUBROUTINE CleanUp ( )


      ! This subroutine cleans up the parent routine before exiting.


      !   ! Deallocate the IfW_Inputs%Position array if it had been allocated.
      !
      !CALL IfW_DestroyInput( IfW_Inputs, ErrStatLcl, ErrMessLcl )


      RETURN

   END SUBROUTINE CleanUp 
   

END SUBROUTINE AD14_CalcOutput

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD14_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMess )
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(AD14_InputType),           INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(AD14_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(AD14_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(AD14_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(AD14_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(AD14_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(AD14_ContinuousStateType), INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at Time
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMess     ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMess  = ""


         ! Compute the first time derivatives of the continuous states here:

!     dxdt%DummyDiscState = 0.


END SUBROUTINE AD14_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD14_UpdateDiscState( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMess )
! Tight coupling routine for updating discrete states
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(AD14_InputType),           INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(AD14_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(AD14_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(AD14_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at Time;
                                                                    !   Output: Discrete states at Time + Interval
      TYPE(AD14_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(AD14_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMess     ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMess  = ""


         ! Update discrete states here:

      ! StateData%DiscState =

END SUBROUTINE AD14_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE AD14_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMess )
! Tight coupling routine for solving for the residual of the constraint state equations
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(AD14_InputType),           INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(AD14_ParameterType),       INTENT(IN   )  :: p           ! Parameters
      TYPE(AD14_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(AD14_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(AD14_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at Time (possibly a guess)
      TYPE(AD14_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(AD14_ConstraintStateType), INTENT(  OUT)  :: z_residual  ! Residual of the constraint state equations using
                                                                    !     the input values described above
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMess     ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMess  = ""


         ! Solve for the constraint states here:

!      z_residual%DummyConstrState = 0.

END SUBROUTINE AD14_CalcConstrStateResidual



!====================================================================================================
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! WE ARE NOT YET IMPLEMENTING THE JACOBIANS...
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------------------------------------------------------------------------------------------------------------------------------

END MODULE AeroDyn14
!**********************************************************************************************************************************

