!**********************************************************************************************************************************
! $Id: InflowWind.f90 125 2014-10-29 22:28:35Z aplatt $
!
! This module is used to read and process the (undisturbed) inflow winds.  It must be initialized
! using InflowWind_Init() with the name of the file, the file type, and possibly reference height and
! width (depending on the type of wind file being used).  This module calls appropriate routines
! in the wind modules so that the type of wind becomes seamless to the user.  InflowWind_End()
! should be called when the program has finshed.
!
! Data are assumed to be in units of meters and seconds.  Z is measured from the ground (NOT the hub!).
!
!  7 Oct 2009    Initial Release with AeroDyn 13.00.00         B. Jonkman, NREL/NWTC
! 14 Nov 2011    v1.00.01b-bjj                                 B. Jonkman
!  1 Aug 2012    v1.01.00a-bjj                                 B. Jonkman
! 10 Aug 2012    v1.01.00b-bjj                                 B. Jonkman
!    Feb 2013    v2.00.00a-adp   conversion to Framework       A. Platt
!
!..................................................................................................................................
! Files with this module:
!  InflowWind_Subs.f90
!  InflowWind.txt       -- InflowWind_Types will be auto-generated based on the descriptions found in this file.
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of InflowWind.
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
! File last committed: $Date: 2014-10-29 16:28:35 -0600 (Wed, 29 Oct 2014) $
! (File) Revision #: $Rev: 125 $
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE InflowWind


   USE                              InflowWind_Types
   USE                              NWTC_Library
   USE                              InflowWind_Subs

      !-------------------------------------------------------------------------------------------------
      ! The included wind modules
      !-------------------------------------------------------------------------------------------------

   USE                              IfW_UniformWind_Types      ! Types for IfW_UniformWind
   USE                              IfW_UniformWind            ! uniform wind files (text files)
!!!   USE                              IfW_FFWind_Types           ! Types for IfW_FFWind
!!!   USE                              IfW_FFWind                 ! full-field binary wind files
!!!   USE                              HAWCWind                   ! full-field binary wind files in HAWC format
!!!   USE                              FDWind                     ! 4-D binary wind files
!!!   USE                              CTWind                     ! coherent turbulence from KH billow - binary file superimposed on another wind type
!!!   USE                              UserWind                   ! user-defined wind module




   IMPLICIT NONE
   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: IfW_Ver = ProgDesc( 'InflowWind', 'v3.00.01a-adp', '01-Nov-2014' )



      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: InflowWind_Init                                   !< Initialization routine
   PUBLIC :: InflowWind_CalcOutput                             !< Calculate the wind velocities
   PUBLIC :: InflowWind_End                                    !< Ending routine (includes clean up)

   PUBLIC :: WindInf_ADhack_diskVel

      ! These routines satisfy the framework, but do nothing at present.
   PUBLIC :: InflowWind_UpdateStates               !< Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
   PUBLIC :: InflowWind_CalcConstrStateResidual    !< Tight coupling routine for returning the constraint state residual
   PUBLIC :: InflowWind_CalcContStateDeriv         !< Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: InflowWind_UpdateDiscState            !< Tight coupling routine for updating discrete states


      ! Not coded
   !NOTE: Jacobians have not been coded.



CONTAINS
!====================================================================================================
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
!! Since this module acts as an interface to other modules, on some things are set before initiating
!! calls to the lower modules.
!----------------------------------------------------------------------------------------------------
SUBROUTINE InflowWind_Init( InitData,   InputGuess,    ParamData,                   &
                     ContStates, DiscStates,    ConstrStateGuess,    OtherStates,   &
                     OutData,    TimeInterval,  InitOutData,                        &
                     ErrStat,    ErrMsg )



         ! Initialization data and guesses

      TYPE(InflowWind_InitInputType),        INTENT(IN   )  :: InitData          !< Input data for initialization
      TYPE(InflowWind_InputType),            INTENT(  OUT)  :: InputGuess        !< An initial guess for the input; the input mesh must be defined
      TYPE(InflowWind_ParameterType),        INTENT(  OUT)  :: ParamData         !< Parameters
      TYPE(InflowWind_ContinuousStateType),  INTENT(  OUT)  :: ContStates        !< Initial continuous states
      TYPE(InflowWind_DiscreteStateType),    INTENT(  OUT)  :: DiscStates        !< Initial discrete states
      TYPE(InflowWind_ConstraintStateType),  INTENT(  OUT)  :: ConstrStateGuess  !< Initial guess of the constraint states
      TYPE(InflowWind_OtherStateType),       INTENT(  OUT)  :: OtherStates       !< Initial other/optimization states
      TYPE(InflowWind_OutputType),           INTENT(  OUT)  :: OutData           !< Initial output (outputs are not calculated; only the output mesh is initialized)
      REAL(DbKi),                            INTENT(IN   )  :: TimeInterval      !< Coupling time interval in seconds: InflowWind does not change this.
      TYPE(InflowWind_InitOutputType),       INTENT(  OUT)  :: InitOutData       !< Initial output data -- Names, units, and version info.


         ! Error Handling

      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat           !< Error status of the operation
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None


         ! Local variables

      TYPE(InflowWind_InputFile)                            :: InputFileData     !< Data from input file

      TYPE(IfW_UniformWind_InitInputType)                   :: Uniform_InitData       !< initialization info
      TYPE(IfW_UniformWind_InputType)                       :: Uniform_InitGuess      !< input positions.
      TYPE(IfW_UniformWind_OutputType)                      :: Uniform_OutData        !< output velocities

!!!      TYPE(IfW_FFWind_InitInputType)                        :: FF_InitData       !< initialization info
!!!      TYPE(IfW_FFWind_InputType)                            :: FF_InitGuess      !< input positions.
!!!      TYPE(IfW_FFWind_ContinuousStateType)                  :: FF_ContStates     !< Unused
!!!      TYPE(IfW_FFWind_DiscreteStateType)                    :: FF_DiscStates     !< Unused
!!!      TYPE(IfW_FFWind_ConstraintStateType)                  :: FF_ConstrStates   !< Unused
!!!      TYPE(IfW_FFWind_OutputType)                           :: FF_OutData        !< output velocities


!!!     TYPE(CTTS_Backgr)                                        :: BackGrndValues


         ! Temporary variables for error handling
      INTEGER(IntKi)                                        :: TmpErrStat
      CHARACTER(LEN(ErrMsg))                                :: TmpErrMsg         !< temporary error message


         ! Local Variables
      INTEGER(IntKi)                                        :: I                 !< Generic counter



         !----------------------------------------------------------------------------------------------
         ! Initialize variables and check to see if this module has been initialized before.
         !----------------------------------------------------------------------------------------------

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Set a few variables.

      ParamData%DT            = TimeInterval             ! InflowWind does not require a specific time interval, so this is never changed.
      CALL NWTC_Init()
      CALL DispNVD( IfW_Ver )



         ! check to see if we are already initialized. Return if it has.
         ! If for some reason a different type of windfile should be used, then call InflowWind_End first, then reinitialize.

      IF ( ParamData%Initialized ) THEN
         CALL SetErrStat( ErrID_Warn, ' InflowWind has already been initialized.', ErrStat, ErrMsg, ' IfW_Init' )
         IF ( ErrStat >= AbortErrLev ) RETURN
      ENDIF


         !----------------------------------------------------------------------------------------------
         ! Read the input file
         !----------------------------------------------------------------------------------------------


         ! Set the names of the files based on the inputfilename
      ParamData%InputFileName = InitData%InputFileName
      CALL GetRoot( ParamData%InputFileName, ParamData%RootFileName )
      ParamData%EchoFileName  = TRIM(ParamData%RootFileName)//".IfW.ech"
      ParamData%SumFileName   = TRIM(ParamData%RootFileName)//".IfW.sum"


         ! Parse all the InflowWind related input files and populate the *_InitDataType derived types

      IF ( InitData%UseInputFile ) THEN
         CALL InflowWind_ReadInput( ParamData%InputFileName, ParamData%EchoFileName, InputFileData, TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF
      ELSE
            !  In the future we will make it possible to just copy the data from the InitData%PassedFileData (of derived type
            !  InflowWind_InputFile), and run things as normal using that data.  For now though, we will not allow this.
         CALL SetErrStat(ErrID_Fatal,' The UseInputFile flag cannot be set to .FALSE. at present.  This feature is not '// &
               'presently supported.',ErrStat,ErrMsg,'InflowWind_Init')
         CALL Cleanup()
         RETURN
      ENDIF





         ! Validate the InflowWind input file information.

      CALL InflowWind_ValidateInput( InputFileData, TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
      IF ( ErrStat>= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      ENDIF


         ! Set the ParamData and OtherStates for InflowWind using the input file information.

      CALL InflowWind_SetParameters( InputFileData, ParamData, OtherStates, TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
      IF ( ErrStat>= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      ENDIF



         ! Allocate arrays for the WriteOutput

      CALL AllocAry( OutData%WriteOutput, ParamData%NumOuts, 'WriteOutput', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
      IF ( ErrStat>= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      ENDIF
      OutData%WriteOutput = 0.0_ReKi
      
      CALL AllocAry( InitOutData%WriteOutputHdr, ParamData%NumOuts, 'WriteOutputHdr', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
      IF ( ErrStat>= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      ENDIF

      CALL AllocAry( InitOutData%WriteOutputUnt, ParamData%NumOuts, 'WriteOutputUnt', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
      IF ( ErrStat>= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      ENDIF
   
      InitOutData%WriteOutputHdr = ParamData%OutParam(1:ParamData%NumOuts)%Name
      InitOutData%WriteOutputUnt = ParamData%OutParam(1:ParamData%NumOuts)%Units     
 


         ! Allocate the arrays for passing points in and velocities out
      IF ( .NOT. ALLOCATED(InputGuess%PositionXYZ) ) THEN
         CALL AllocAry( InputGuess%PositionXYZ, 3, InitData%NumWindPoints, &
                     "Array of positions at which to find wind velocities", TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
         IF ( ErrStat>= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF
      ENDIF
      IF ( .NOT. ALLOCATED(OutData%VelocityUVW) ) THEN
         CALL AllocAry( OutData%VelocityUVW, 3, InitData%NumWindPoints, &
                     "Array of wind velocities returned by InflowWind", TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
         IF ( ErrStat>= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF
      ENDIF


         ! Check that the arrays for the points and velocities are the same size
      IF ( SIZE( InputGuess%PositionXYZ, DIM = 2 ) /= SIZE( OutData%VelocityUVW, DIM = 2 ) ) THEN
         CALL SetErrStat(ErrID_Fatal,' Programming error: Different number of XYZ coordinates and expected output velocities.', &
                     ErrStat,ErrMsg,'InflowWind_Init')
      ENDIF



      !-----------------------------------------------------------------
      ! Initialize the submodules based on the WindType
      !-----------------------------------------------------------------


      SELECT CASE ( ParamData%WindType )


         CASE ( Steady_WindNumber )

               !  This is a simplified case of the Uniform wind.  For this, we set the OtherStates data manually and don't
               !  call UniformWind_Init.  We do however call it for the calculations.


               ! Set InitData information -- It isn't necessary to do this since that information is only used in the Uniform_Init routine, which is not called.

               ! Set the Otherstates information
            OtherStates%UniformWind%NumDataLines   =  1_IntKi
            OtherStates%UniformWind%RefHt          =  InputFileData%Steady_RefHt
            OtherStates%UniformWind%RefLength      =  1.0_ReKi    ! This is not used since no shear gusts are used.  Set to 1.0 so calculations don't bomb. 
            OtherStates%UniformWind%TimeIndex      =  1           ! Used in UniformWind as a check if initialization was done.


            IF (.NOT. ALLOCATED(OtherStates%UniformWind%Tdata) ) THEN
               CALL AllocAry( OtherStates%UniformWind%Tdata, OtherStates%UniformWind%NumDataLines, 'Uniform wind time', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF
         
            IF (.NOT. ALLOCATED(OtherStates%UniformWind%V) ) THEN
               CALL AllocAry( OtherStates%UniformWind%V, OtherStates%UniformWind%NumDataLines, 'Uniform wind horizontal wind speed', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF
         
            IF (.NOT. ALLOCATED(OtherStates%UniformWind%Delta) ) THEN
               CALL AllocAry( OtherStates%UniformWind%Delta, OtherStates%UniformWind%NumDataLines, 'Uniform wind direction', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF
         
            IF (.NOT. ALLOCATED(OtherStates%UniformWind%VZ) ) THEN
               CALL AllocAry( OtherStates%UniformWind%VZ, OtherStates%UniformWind%NumDataLines, 'Uniform vertical wind speed', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF
         
            IF (.NOT. ALLOCATED(OtherStates%UniformWind%HShr) ) THEN
               CALL AllocAry( OtherStates%UniformWind%HShr, OtherStates%UniformWind%NumDataLines, 'Uniform horizontal linear shear', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF
         
            IF (.NOT. ALLOCATED(OtherStates%UniformWind%VShr) ) THEN
               CALL AllocAry( OtherStates%UniformWind%VShr, OtherStates%UniformWind%NumDataLines, 'Uniform vertical power-law shear exponent', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF
         
            IF (.NOT. ALLOCATED(OtherStates%UniformWind%VLinShr) ) THEN
               CALL AllocAry( OtherStates%UniformWind%VLinShr, OtherStates%UniformWind%NumDataLines, 'Uniform vertical linear shear', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF
         
            IF (.NOT. ALLOCATED(OtherStates%UniformWind%VGust) ) THEN
               CALL AllocAry( OtherStates%UniformWind%VGust, OtherStates%UniformWind%NumDataLines, 'Uniform gust velocity', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF

               ! Set the array information
            
            OtherStates%UniformWind%Tdata(  :) = 0.0_ReKi
            OtherStates%UniformWind%V(      :) = InputFileData%Steady_HWindSpeed
            OtherStates%UniformWind%Delta(  :) = 0.0_ReKi
            OtherStates%UniformWind%VZ(     :) = 0.0_ReKi
            OtherStates%UniformWind%HShr(   :) = 0.0_ReKi
            OtherStates%UniformWind%VShr(   :) = InputFileData%Steady_PLexp
            OtherStates%UniformWind%VLinShr(:) = 0.0_ReKi
            OtherStates%UniformWind%VGust(  :) = 0.0_ReKi



               ! Now we have in effect initialized the IfW_UniformWind module, so set the parameters
            ParamData%UniformWind%ReferenceHeight  =  Uniform_InitData%RefLength
            ParamData%UniformWind%RefLength        =  InputFileData%Steady_RefHt
            ParamData%UniformWind%Initialized      =  .TRUE.

               ! Output reporting data
            InitOutData%HubHeight                  =  InputFileData%Steady_RefHt
            InitOutData%WindFileConstantDT         =  .FALSE.
            InitOutData%WindFileDT                 =  0.0_ReKi
            InitOutData%WindFileTRange             =  (/ 0.0_ReKi, 0.0_ReKi /)
            InitOutData%WindFileNumTSteps          =  1_IntKi
       



         CASE ( Uniform_WindNumber )

               ! Set InitData information
            Uniform_InitData%ReferenceHeight =  InputFileData%Uniform_RefHt
            Uniform_InitData%RefLength       =  InputFileData%Uniform_RefLength 
            Uniform_InitData%WindFileName    =  InputFileData%Uniform_FileName

               ! Initialize the UniformWind module
            CALL IfW_UniformWind_Init(Uniform_InitData, InputGuess%PositionXYZ, ParamData%UniformWind, OtherStates%UniformWind, &
                        Uniform_OutData,    TimeInterval,  InitOutData%UniformWind,  TmpErrStat,          TmpErrMsg)

            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_Init' )
            IF ( ErrStat >= AbortErrLev ) RETURN


               ! Copy Relevant info over to InitOutData

               ! Allocate and copy over the WriteOutputHdr info
            CALL AllocAry( InitOutData%WriteOutputHdr, SIZE(InitOutData%UniformWind%WriteOutputHdr,1), &
                              'Empty array for names of outputable information.', TmpErrStat, TmpErrMsg )
            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_Init' )
            IF ( ErrStat >= AbortErrLev ) RETURN
            InitOutData%WriteOutputHdr    =  InitOutData%UniformWind%WriteOutputHdr

               ! Allocate and copy over the WriteOutputUnt info
            CALL AllocAry( InitOutData%WriteOutputUnt, SIZE(InitOutData%UniformWind%WriteOutputUnt,1), &
                           'Empty array for units of outputable information.', TmpErrStat, TmpErrMsg )
            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_Init' )
            IF ( ErrStat >= AbortErrLev ) RETURN

            InitOutData%WriteOutputUnt    =  InitOutData%UniformWind%WriteOutputUnt

               ! Copy the hub height info over
            InitOutData%HubHeight         =  InitOutData%UniformWind%HubHeight

               ! Timestep information
            InitOutData%WindFileConstantDT=  InitOutData%UniformWind%WindFileConstantDT
            InitOutData%WindFileDT        =  InitOutData%UniformWind%WindFileDT
            InitOutData%WindFileTRange    =  InitOutData%UniformWind%WindFileTRange
            InitOutData%WindFileNumTSteps =  InitOutData%UniformWind%WindFileNumTSteps

!FIXME: add warning if the first data point in the windfile has an angle and we are using a coordinate transform (rotation of windfield).


         CASE ( TSFF_WindNumber )
               ! Initialize the TSFFWind module
            CALL SetErrStat( ErrID_Fatal,' TSFF winds not supported yet.',ErrStat,ErrMsg,'InflowWind_Init')



         CASE ( BladedFF_WindNumber )
               ! Initialize the BladedFFWind module
            CALL SetErrStat( ErrID_Fatal,' BladedFF winds not supported yet.',ErrStat,ErrMsg,'InflowWind_Init')



         CASE ( HAWC_WindNumber )
               ! Initialize the HAWCWind module
            CALL SetErrStat( ErrID_Fatal,' HAWC winds not supported yet.',ErrStat,ErrMsg,'InflowWind_Init')



         CASE ( User_WindNumber )
               ! Initialize the User_Wind module
            CALL SetErrStat( ErrID_Fatal,' User winds not supported yet.',ErrStat,ErrMsg,'InflowWind_Init')



         CASE DEFAULT  ! keep this check to make sure that all new wind types have been accounted for
            CALL SetErrStat(ErrID_Fatal,' Undefined wind type.',ErrStat,ErrMsg,'InflowWind_Init()')




      END SELECT


      IF ( ParamData%CTTS_Flag ) THEN
         ! Initialize the CTTS_Wind module
      ENDIF



!!!         !----------------------------------------------------------------------------------------------
!!!         ! Check for coherent turbulence file (KH superimposed on a background wind file)
!!!         ! Initialize the CTWind module and initialize the module of the other wind type.
!!!         !----------------------------------------------------------------------------------------------
!!!
!!!      IF ( ParamData%WindType == CTP_WindNumber ) THEN
!!!
!!!!FIXME: remove this error message when we add CTP_Wind in
!!!            CALL SetErrStat( ErrID_Fatal, ' InflowWind cannot currently handle the CTP_Wind type.', ErrStat, ErrMsg, ' IfW_Init' )
!!!            RETURN
!!!
!!!         CALL CTTS_Init(UnWind, ParamData%WindFileName, BackGrndValues, ErrStat, ErrMsg)
!!!         IF (ErrStat /= 0) THEN
!!!   !         CALL IfW_End( ParamData, ErrStat )
!!!   !FIXME: cannot call IfW_End here -- requires InitData to be INOUT. Not allowed by framework.
!!!   !         CALL IfW_End( InitData, ParamData, ContStates, DiscStates, ConstrStateGuess, OtherStates, &
!!!   !                       OutData, ErrStat, ErrMsg )
!!!            ParamData%WindType = Undef_Wind
!!!            ErrStat  = 1
!!!            RETURN
!!!         END IF
!!!
!!!   !FIXME: check this
!!!         ParamData%WindFileName = BackGrndValues%WindFile
!!!         ParamData%WindType = BackGrndValues%WindType
!!!   !      CTTS_Flag  = BackGrndValues%CoherentStr
!!!         ParamData%CTTS_Flag  = BackGrndValues%CoherentStr    ! This might be wrong
!!!
!!!      ELSE
!!!
!!!         ParamData%CTTS_Flag  = .FALSE.
!!!
!!!      END IF
!!!
!!!         !----------------------------------------------------------------------------------------------
!!!         ! Initialize based on the wind type
!!!         !----------------------------------------------------------------------------------------------
!!!
!!!      SELECT CASE ( ParamData%WindType )
!!!
!!!         CASE (Uniform_WindNumber)
!!!
!!!            Uniform_InitData%ReferenceHeight = InitData%ReferenceHeight
!!!            Uniform_InitData%Width           = InitData%Width
!!!            Uniform_InitData%WindFileName    = ParamData%WindFileName
!!!
!!!            CALL IfW_UniformWind_Init(Uniform_InitData,   Uniform_InitGuess,  ParamData%UniformWind,                         &
!!!                                 Uniform_ContStates, Uniform_DiscStates, Uniform_ConstrStates,     OtherStates%UniformWind,  &
!!!                                 Uniform_OutData,    TimeInterval,  InitOutData%UniformWind,  TmpErrStat,          TmpErrMsg)
!!!
!!!               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_Init' )
!!!               IF ( ErrStat >= AbortErrLev ) RETURN
!!!
!!!
!!!              ! Copy Relevant info over to InitOutData
!!!
!!!                  ! Allocate and copy over the WriteOutputHdr info
!!!               CALL AllocAry( InitOutData%WriteOutputHdr, SIZE(InitOutData%UniformWind%WriteOutputHdr,1), &
!!!                              'Empty array for names of outputable information.', TmpErrStat, TmpErrMsg )
!!!                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_Init' )
!!!                  IF ( ErrStat >= AbortErrLev ) RETURN
!!!               InitOutData%WriteOutputHdr    =  InitOutData%UniformWind%WriteOutputHdr
!!!
!!!                  ! Allocate and copy over the WriteOutputUnt info
!!!               CALL AllocAry( InitOutData%WriteOutputUnt, SIZE(InitOutData%UniformWind%WriteOutputUnt,1), &
!!!                              'Empty array for units of outputable information.', TmpErrStat, TmpErrMsg )
!!!                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_Init' )
!!!                  IF ( ErrStat >= AbortErrLev ) RETURN
!!!
!!!               InitOutData%WriteOutputUnt    =  InitOutData%UniformWind%WriteOutputUnt
!!!
!!!                  ! Copy the hub height info over
!!!               InitOutData%HubHeight         =  InitOutData%UniformWind%HubHeight
!!!
!!!!           IF (CTTS_Flag) CALL CTTS_SetRefVal(FileInfo%ReferenceHeight, 0.5*FileInfo%Width, ErrStat)  !FIXME: check if this was originally used
!!!!           IF (ErrStat == ErrID_None .AND. ParamData%CTTS_Flag) &
!!!!              CALL CTTS_SetRefVal(InitData%ReferenceHeight, REAL(0.0, ReKi), ErrStat, ErrMsg)      !FIXME: will need to put this routine in the Init of CT
!!!
!!!
!!!         CASE (FF_WindNumber)
!!!
!!!            FF_InitData%ReferenceHeight = InitData%ReferenceHeight
!!!            FF_InitData%Width           = InitData%Width
!!!            FF_InitData%WindFileName    = ParamData%WindFileName
!!!
!!!            CALL IfW_FFWind_Init(FF_InitData,   FF_InitGuess,  ParamData%FFWind,                         &
!!!                                 FF_ContStates, FF_DiscStates, FF_ConstrStates,     OtherStates%FFWind,  &
!!!                                 FF_OutData,    TimeInterval,  InitOutData%FFWind,  TmpErrStat,          TmpErrMsg)
!!!
!!!               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_Init' )
!!!               IF ( ErrStat >= AbortErrLev ) RETURN
!!!
!!!
!!!              ! Copy Relevant info over to InitOutData
!!!
!!!                  ! Allocate and copy over the WriteOutputHdr info
!!!               CALL AllocAry( InitOutData%WriteOutputHdr, SIZE(InitOutData%FFWind%WriteOutputHdr,1), &
!!!                              'Empty array for names of outputable information.', TmpErrStat, TmpErrMsg )
!!!                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_Init' )
!!!                  IF ( ErrStat >= AbortErrLev ) RETURN
!!!
!!!               InitOutData%WriteOutputHdr    =  InitOutData%FFWind%WriteOutputHdr
!!!
!!!                  ! Allocate and copy over the WriteOutputUnt info
!!!               CALL AllocAry( InitOutData%WriteOutputUnt, SIZE(InitOutData%FFWind%WriteOutputUnt,1), &
!!!                              'Empty array for units of outputable information.', TmpErrStat, TmpErrMsg )
!!!                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_Init' )
!!!                  IF ( ErrStat >= AbortErrLev ) RETURN
!!!
!!!               InitOutData%WriteOutputUnt    =  InitOutData%FFWind%WriteOutputUnt
!!!
!!!                  ! Copy the hub height info over
!!!               InitOutData%HubHeight         =  InitOutData%FFWind%HubHeight
!!!
!!!            !FIXME: Fix this when CTTS_Wind is available
!!!!               ! Set CT parameters
!!!!            IF ( ErrStat == ErrID_None .AND. ParamData%CTTS_Flag ) THEN
!!!!               Height     = FF_GetValue('HubHeight', ErrStat, ErrMsg)
!!!!               IF ( ErrStat /= 0 ) Height = InitData%ReferenceHeight
!!!!
!!!!               HalfWidth  = 0.5*FF_GetValue('GridWidth', ErrStat, ErrMsg)
!!!!               IF ( ErrStat /= 0 ) HalfWidth = 0
!!!!
!!!!               CALL CTTS_SetRefVal(Height, HalfWidth, ErrStat, ErrMsg)
!!!!            END IF
!!!
!!!
!!!         CASE (User_WindNumber)
!!!
!               ! Timestep information
!            InitOutData%WindFileConstantDT=  InitOutData%UniformWind%WindFileConstantDT
!            InitOutData%WindFileDT        =  InitOutData%UniformWind%WindFileDT
!            InitOutData%WindFileTRange    =  InitOutData%UniformWind%WindFileTRange
!            InitOutData%WindFileNumTSteps =  InitOutData%UniformWind%WindFileNumTSteps
!!!               !FIXME: remove this error message when we add UD_Wind in
!!!            CALL SetErrStat( ErrID_Fatal, ' InflowWind cannot currently handle the UD_Wind type.', ErrStat, ErrMsg, ' IfW_Init' )
!!!            RETURN
!!!
!!!!            CALL UsrWnd_Init(ErrStat)
!!!
!!!
!!!         CASE (FD_WindNumber)
!!!
!!!               !FIXME: remove this error message when we add FD_Wind in
!!!            CALL SetErrStat( ErrID_Fatal, ' InflowWind cannot currently handle the FD_Wind type.', ErrStat, ErrMsg, ' IfW_Init' )
!!!            RETURN
!!!
!!!!           CALL IfW_FDWind_Init(UnWind, ParamData%WindFileName, InitData%ReferenceHeight, ErrStat)
!!!
!               ! Timestep information
!            InitOutData%WindFileConstantDT=  InitOutData%UniformWind%WindFileConstantDT
!            InitOutData%WindFileDT        =  InitOutData%UniformWind%WindFileDT
!            InitOutData%WindFileTRange    =  InitOutData%UniformWind%WindFileTRange
!            InitOutData%WindFileNumTSteps =  InitOutData%UniformWind%WindFileNumTSteps
!!!
!!!         CASE (HAWC_WindNumber)
!!!
!!!               !FIXME: remove this error message when we add HAWC_Wind in
!!!            CALL SetErrStat( ErrID_Fatal, ' InflowWind cannot currently handle the HAWC_Wind type.', ErrStat, ErrMsg, ' IfW_Init' )
!!!            RETURN
!!!
!!!!           CALL HW_Init( UnWind, ParamData%WindFileName, ErrStat )
!               ! Timestep information
!            InitOutData%WindFileConstantDT=  InitOutData%UniformWind%WindFileConstantDT
!            InitOutData%WindFileDT        =  InitOutData%UniformWind%WindFileDT
!            InitOutData%WindFileTRange    =  InitOutData%UniformWind%WindFileTRange
!            InitOutData%WindFileNumTSteps =  InitOutData%UniformWind%WindFileNumTSteps
!!!
!!!
!!!         CASE DEFAULT
!!!
!!!            CALL SetErrStat( ErrID_Fatal, ' Error: Undefined wind type in WindInflow_Init()', ErrStat, ErrMsg, ' IfW_Init' )
!!!            RETURN
!!!
!!!      END SELECT



         ! If we've arrived here, we haven't reached an AbortErrLev:
      ParamData%Initialized = .TRUE.

         ! Set the version information in InitOutData
      InitOutData%Ver   = IfW_Ver



      RETURN


   !----------------------------------------------------------------------------------------------------
CONTAINS

   SUBROUTINE CleanUp()

      ! add in stuff that we need to dispose of here
      CALL InflowWind_DestroyInputFile( InputFileData, TmpErrsTat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')

      CALL InflowWind_DestroyInputFile( InputFileData,  TmpErrStat, TmpErrMsg );  CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,'InflowWind_Init')


   END SUBROUTINE CleanUp



END SUBROUTINE InflowWind_Init



!====================================================================================================
!> This routine takes an input dataset of type InputType which contains a position array of dimensions 3*n. It then calculates
!! and returns the output dataset of type OutputType which contains a corresponding velocity array of dimensions 3*n. The input
!! array contains XYZ triplets for each position of interest (first index is X/Y/Z for values 1/2/3, second index is the point
!! number to evaluate). The returned values in the OutputData are similar with U/V/W for the first index of 1/2/3.
!!
!! _Coordinate transformation:_
!! The coordinates passed in are copied to the PositionXYZPrime array, then rotated by -(ParamData%PropogationDir) (radians) so
!! that the wind direction lies along the X-axis (all wind files are given this way).  The submodules are then called with
!! these PositionXYZPrime coordinates.
!!
!! After the calculation by the submodule, the PositionXYZPrime coordinate array is deallocated.  The returned VelocityUVW
!! array is then rotated by ParamData%PropogationDir so that it now corresponds the the global coordinate UVW values for wind
!! with that direction.
!----------------------------------------------------------------------------------------------------
SUBROUTINE InflowWind_CalcOutput( Time, InputData, ParamData, &
                              ContStates, DiscStates, ConstrStates, &   ! Framework required states -- empty in this case.
                              OtherStates, OutputData, ErrStat, ErrMsg )


         ! Inputs / Outputs

      REAL(DbKi),                               INTENT(IN   )  :: Time              !< Current simulation time in seconds
      TYPE(InflowWind_InputType),               INTENT(IN   )  :: InputData         !< Inputs at Time
      TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: ParamData         !< Parameters
      TYPE(InflowWind_ContinuousStateType),     INTENT(IN   )  :: ContStates        !< Continuous states at Time
      TYPE(InflowWind_DiscreteStateType),       INTENT(IN   )  :: DiscStates        !< Discrete states at Time
      TYPE(InflowWind_ConstraintStateType),     INTENT(IN   )  :: ConstrStates      !< Constraint states at Time
      TYPE(InflowWind_OtherStateType),          INTENT(INOUT)  :: OtherStates       !< Other/optimization states at Time
      TYPE(InflowWind_OutputType),              INTENT(INOUT)  :: OutputData        !< Outputs computed at Time (IN for mesh reasons and data allocation)

      INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< Error status of the operation
      CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None


         ! Local variables

      TYPE(IfW_UniformWind_InitInputType)                      :: Uniform_InitData       !< initialization info
      TYPE(IfW_UniformWind_InputType)                          :: Uniform_InData         !< input positions.
      TYPE(IfW_UniformWind_OutputType)                         :: Uniform_OutData        !< output velocities

!!!      TYPE(IfW_FFWind_InitInputType)                           :: FF_InitData       !< initialization info
!!!      TYPE(IfW_FFWind_InputType)                               :: FF_InData         !< input positions.
!!!      TYPE(IfW_FFWind_ContinuousStateType)                     :: FF_ContStates     !< Unused
!!!      TYPE(IfW_FFWind_DiscreteStateType)                       :: FF_DiscStates     !< Unused
!!!      TYPE(IfW_FFWind_ConstraintStateType)                     :: FF_ConstrStates   !< Unused
!!!      TYPE(IfW_FFWind_OutputType)                              :: FF_OutData        !< output velocities

!      REAL(ReKi)                                               :: RotToWind(3,3)       !< Rotation matrix for rotating from global XYZ to windfile X'Y'Z'
!      REAL(ReKi)                                               :: RotFromWind(3,3)     !< Rotation matrix for rotating from windfile X'Y'Z' to global XYZ
      REAL(ReKi), ALLOCATABLE                                  :: PositionXYZprime(:,:)   !< PositionXYZ array in the prime (wind) coordinates

      INTEGER(IntKi)                                           :: I                    !< Generic counters

      INTEGER(IntKi)                                           :: J                    !FIXME: remove after testing
      LOGICAL                                                  :: samecheck            !FIXME: remove after testing
      REAL(ReKi), ALLOCATABLE       :: TmpPositionXYZcheck(:,:)   !< PositionXYZ array returned from prime (wind) coordinates FIXME: remove after testing

         ! Temporary variables for error handling
      INTEGER(IntKi)                                           :: TmpErrStat
      CHARACTER(LEN(ErrMsg))                                   :: TmpErrMsg            ! temporary error message



         ! Initialize ErrStat
      ErrStat  = ErrID_None
      ErrMsg   = ""


         ! Allocate the velocity array to get out
      IF ( .NOT. ALLOCATED(OutputData%VelocityUVW) ) THEN
         CALL AllocAry( OutputData%VelocityUVW, 3, SIZE(InputData%PositionXYZ,DIM=2), &
                     "Velocity array returned from IfW_CalcOutput", TmpErrStat, TmpErrMsg )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_CalcOutput' )
         IF ( ErrStat >= AbortErrLev ) RETURN
      ELSEIF ( SIZE(InputData%PositionXYZ,DIM=2) /= SIZE(OutputData%VelocityUVW,DIM=2) ) THEN
         CALL SetErrStat( ErrID_Fatal," Programming error: Position and Velocity arrays are not sized the same.",  &
               ErrStat, ErrMsg, ' IfW_CalcOutput')
      ENDIF


!         !> Perform coordinate transformation using the PropogationDir value read in from the input file (the value stored in parameters
!         !! was converted from degrees to radians already).
!         !! @note    The PropogationDir is given in Meteorological \f$\Delta\phi\f$, so this is the negative of the \f$\Delta\phi\f$
!         !!          for polar coordinates
!
!         ! Create the rotation matrices -- rotate from XYZ to X'Y'Z' (wind aligned along X) coordinates
!      RotToWind(1,:) = (/    COS(-ParamData%PropogationDir),   SIN(-ParamData%PropogationDir),     0.0_ReKi  /)  
!      RotToWind(2,:) = (/   -SIN(-ParamData%PropogationDir),   COS(-ParamData%PropogationDir),     0.0_ReKi  /)  
!      RotToWind(3,:) = (/                          0.0_ReKi,                         0.0_ReKi,     1.0_ReKi  /)  
!
!         ! Create the rotation matrices -- rotate from X'Y'Z' (wind aligned along X) to global XYZ coordinates
!      RotFromWind =  TRANSPOSE(RotToWind)


         !> Make a copy of the InputData%PositionXYZ coordinates with the applied rotation matrix...
         !! This copy is made because if we translate it to the prime coordinates, then back again, we
         !! may shift the points by some small amount of machine error, and we don't want to do that.
         !!
         !! Note that we allocate this at every call to CalcOutput.  The reason is that we may call CalcOutput
         !! multiple times in each timestep with different sized InputData%PositionXYZ arrays (if calling for LIDAR
         !! data etc).  We don't really want to have extra copies of OtherStates lying around in order to be able to
         !! do this.
      CALL AllocAry( PositionXYZprime, 3, SIZE(InputData%PositionXYZ,DIM=2), &
                  "Array for holding the XYZprime position data", TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_CalcOutput' )
      IF ( ErrStat >= AbortErrLev ) RETURN


         ! Apply the coordinate transformation to the PositionXYZ coordinates to get the PositionXYZprime coordinate list
         ! If the PropogationDir is zero, we don't need to apply this and will simply copy the data.  Repeat for the WindViXYZ.
      IF ( EqualRealNos (ParamData%PropogationDir, 0.0_ReKi) ) THEN
         PositionXYZprime  =  InputData%PositionXYZ
      ELSE
         DO I  = 1,SIZE(InputData%PositionXYZ,DIM=2)
            PositionXYZprime(:,I)   =  MATMUL( ParamData%RotToWind, InputData%PositionXYZ(:,I) )
         ENDDO
      ENDIF




         ! Compute the wind velocities by stepping through all the data points and calling the appropriate GetWindSpeed routine
      SELECT CASE ( ParamData%WindType )


         CASE (Steady_WindNumber)

               ! Move the arrays for the Velocity information
            CALL MOVE_ALLOC( OutputData%VelocityUVW,  Uniform_OutData%Velocity )

               ! InputData only contains the Position array, so we can pass that directly.
            CALL  IfW_UniformWind_CalcOutput(  Time, InputData%PositionXYZ, ParamData%UniformWind, OtherStates%UniformWind, &
                                          Uniform_OutData, TmpErrStat, TmpErrMsg)

               ! Move the arrays back.  note that these are in the prime (wind file) coordinate frame still.
            CALL MOVE_ALLOC( Uniform_OutData%Velocity,   OutputData%VelocityUVW )

            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_CalcOutput' )
            IF ( ErrStat >= AbortErrLev ) RETURN

               ! Call IfW_UniforWind_CalcOutput again in order to get the values needed for the OutList -- note that we do not report errors from this
            IF ( ParamData%NWindVel >= 1_IntKi ) THEN
                  ! Move the arrays for the Velocity information
               CALL MOVE_ALLOC( OtherStates%WindViUVW,  Uniform_OutData%Velocity )

               CALL  IfW_UniformWind_CalcOutput(  Time, ParamData%WindViXYZprime, ParamData%UniformWind, &
                                             OtherStates%UniformWind, &
                                             Uniform_OutData, TmpErrStat, TmpErrMsg)
               TmpErrStat  = ErrID_None
               TmpErrMsg   = ''

                  ! Move the arrays back.  note that these are in the prime (wind file) coordinate frame still.
               CALL MOVE_ALLOC( Uniform_OutData%Velocity, OtherStates%WindViUVW )
            ENDIF




         CASE (Uniform_WindNumber)


               ! Move the arrays for the Position and Velocity information
            CALL MOVE_ALLOC( OutputData%VelocityUVW,  Uniform_OutData%Velocity )


               ! InputData only contains the Position array, so we can pass that directly.
            CALL  IfW_UniformWind_CalcOutput(  Time, InputData%PositionXYZ, ParamData%UniformWind, OtherStates%UniformWind, &
                                          Uniform_OutData, TmpErrStat, TmpErrMsg)

               ! Move the arrays back.  note that these are in the prime (wind file) coordinate frame still.
            CALL MOVE_ALLOC( Uniform_OutData%Velocity,   OutputData%VelocityUVW )

            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_CalcOutput' )
            IF ( ErrStat >= AbortErrLev ) RETURN


               ! Call IfW_UniforWind_CalcOutput again in order to get the values needed for the OutList -- note that we do not report errors from this
            IF ( ParamData%NWindVel >= 1_IntKi ) THEN
                  ! Move the arrays for the Velocity information
               CALL MOVE_ALLOC( OtherStates%WindViUVW,  Uniform_OutData%Velocity )

               CALL  IfW_UniformWind_CalcOutput(  Time, ParamData%WindViXYZprime, ParamData%UniformWind, &
                                             OtherStates%UniformWind, &
                                             Uniform_OutData, TmpErrStat, TmpErrMsg)
               TmpErrStat  = ErrID_None
               TmpErrMsg   = ''

                  ! Move the arrays back.  note that these are in the prime (wind file) coordinate frame still.
               CALL MOVE_ALLOC( Uniform_OutData%Velocity, OtherStates%WindViUVW )
            ENDIF





!!!         CASE (FF_WindNumber)
!!!
!!!               ! Allocate the position array to pass in
!!!            CALL AllocAry( FF_InData%Position, 3, SIZE(InputData%Position,2), &
!!!                           "Position grid for passing to IfW_FFWind_CalcOutput", TmpErrStat, TmpErrMsg )
!!!            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_CalcOutput' )
!!!            IF ( ErrStat >= AbortErrLev ) RETURN
!!!
!!!            ! Copy positions over
!!!            FF_InData%Position   = InputData%Position
!!!
!!!            CALL  IfW_FFWind_CalcOutput(  Time,          FF_InData,     ParamData%FFWind,                         &
!!!                                          FF_ContStates, FF_DiscStates, FF_ConstrStates,     OtherStates%FFWind,  &
!!!                                          FF_OutData,    TmpErrStat,    TmpErrMsg)
!!!
!!!               ! Copy the velocities over
!!!            OutputData%Velocity  = FF_OutData%Velocity
!!!
!!!            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_CalcOutput' )
!!!            IF ( ErrStat >= AbortErrLev ) RETURN


!!!               OutputData%Velocity(:,PointCounter) = FF_GetWindSpeed(     Time, InputData%Position(:,PointCounter), ErrStat, ErrMsg)


!!!         CASE (User_WindNumber)

!!!               OutputData%Velocity(:,PointCounter) = UsrWnd_GetWindSpeed( Time, InputData%Position(:,PointCounter), ErrStat )!, ErrMsg)


!!!         CASE (FD_WindNumber)

!!!               OutputData%Velocity(:,PointCounter) = FD_GetWindSpeed(     Time, InputData%Position(:,PointCounter), ErrStat )



!!!         CASE (HAWC_WindNumber)

!!!               OutputData%Velocity(:,PointCounter) = HW_GetWindSpeed(     Time, InputData%Position(:,PointCounter), ErrStat )



            ! If it isn't one of the above cases, we have a problem and won't be able to continue

         CASE DEFAULT

            CALL SetErrStat( ErrID_Fatal, ' Error: Undefined wind type in IfW_CalcOutput. ' &
                      //'Call WindInflow_Init() before calling this function.', ErrStat, ErrMsg, ' IfW_CalcOutput' )

            OutputData%VelocityUVW(:,:) = 0.0
            RETURN

      END SELECT


            ! Add coherent turbulence to background wind

!!!         IF (ParamData%CTTS_Flag) THEN
!!!
!!!            DO PointCounter = 1, SIZE(InputData%Position, 2)
!!!
!!!               TempWindSpeed = CTTS_GetWindSpeed(     Time, InputData%Position(:,PointCounter), ErrStat, ErrMsg )
!!!
!!!                  ! Error Handling -- move ErrMsg inside CTTS_GetWindSPeed and simplify
!!!               IF (ErrStat >= ErrID_Severe) THEN
!!!                  ErrMsg   = 'IfW_CalcOutput: Error in CTTS_GetWindSpeed for point number '//TRIM(Num2LStr(PointCounter))
!!!                  EXIT        ! Exit the loop
!!!               ENDIF
!!!
!!!               OutputData%Velocity(:,PointCounter) = OutputData%Velocity(:,PointCounter) + TempWindSpeed
!!!
!!!            ENDDO
!!!
!!!               ! If something went badly wrong, Return
!!!            IF (ErrStat >= ErrID_Severe ) RETURN
!!!
!!!         ENDIF
!!!
      !ENDIF


         ! The VelocityUVW array data that has been returned from the sub-modules is in the wind file (X'Y'Z') coordinates at
         ! this point.  These must be rotated to the global XYZ coordinates.  So now we apply the coordinate transformation
         ! to the VelocityUVW(prime) coordinates (in wind X'Y'Z' coordinate frame) returned from the submodules to the XYZ
         ! coordinate frame, but only if PropogationDir is not zero.  This is only a rotation of the returned wind field, so
         ! UVW contains the direction components of the wind at XYZ after translation from the U'V'W' wind velocity components
         ! in the X'Y'Z' (wind file) coordinate frame.
      IF ( .NOT. EqualRealNos (ParamData%PropogationDir, 0.0_ReKi) ) THEN
         DO I  = 1,SIZE(OutputData%VelocityUVW,DIM=2)
            OutputData%VelocityUVW(:,I)   =  MATMUL( ParamData%RotFromWind, OutputData%VelocityUVW(:,I) )
         ENDDO
      ENDIF

         ! We also need to rotate the reference frame for the WindViUVW array
      IF ( .NOT. EqualRealNos (ParamData%PropogationDir, 0.0_ReKi) ) THEN
         DO I  = 1,SIZE(OtherStates%WindViUVW,DIM=2)
            OtherStates%WindViUVW(:,I)   =  MATMUL( ParamData%RotFromWind, OtherStates%WindViUVW(:,I) )
         ENDDO
      ENDIF


!!!!!!!!!!!!!!!!!!!!!
!!FIXME: temp testing
!      CALL AllocAry( TmpPositionXYZcheck, 3, SIZE(InputData%PositionXYZ,DIM=2), &
!                  "Temporary array for holding the XYZprime position data", TmpErrStat, TmpErrMsg )
!      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_CalcOutput' )
!      IF ( ErrStat >= AbortErrLev ) RETURN
!
!      DO I  = 1,SIZE(InputData%PositionXYZ,DIM=2)
!         TmpPositionXYZcheck(:,I)  =  MATMUL( ParamData%RotFromWind, PositionXYZprime(:,I) )
!      ENDDO
!
!      samecheck   = .TRUE.
!      DO I  =  1,SIZE(InputData%PositionXYZ,DIM=2)
!         DO J  =  1,SIZE(InputData%PositionXYZ,DIM=1)
!            IF ( .NOT. EqualRealNos(TmpPositionXYZcheck(J,I), InputData%PositionXYZ(J,I)) )    samecheck = .FALSE.
!         ENDDO
!      ENDDO
!
!      print*,'  Matmul test: '
!      print*,'       PropogationDir = ',ParamData%PropogationDir
!      print*,'       RotToWind        ',RotToWind(1,:)
!      print*,'                        ',RotToWind(2,:)
!      print*,'                        ',RotToWind(3,:)
!      print*,''
!      DO I  = 1,SIZE(InputData%PositionXYZ,DIM=2)
!         print*,''
!         print*,'       PositionXYZ     (I) ',InputData%PositionXYZ(:,I)
!         print*,'       PositionXYZprime(I) ',PositionXYZprime(:,I)
!         print*,'    TmpPositionXYZcheck(I) ',TmpPositionXYZcheck(:,I),'     same: ',samecheck
!         print*,'       VelocityUVW     (I) ',OutputData%VelocityUVW(:,I)
!      ENDDO
!!!!!!!!!!!!!!!!!!!!!

      ! Done with the prime coordinates for the XYZ position information that was passed in.
   IF (ALLOCATED(PositionXYZprime)) DEALLOCATE(PositionXYZprime)

      !-----------------------------
      ! Outputs: OtherState%AllOuts
      !-----------------------------

   CALL SetAllOuts( ParamData, OtherStates, TmpErrStat, TmpErrMsg ) 
   
      ! Map to the outputs
   DO I = 1,ParamData%NumOuts  ! Loop through all selected output channels
      OutputData%WriteOutput(I) = ParamData%OutParam(I)%SignM * OtherStates%AllOuts( ParamData%OutParam(I)%Indx )
   ENDDO             ! I - All selected output channels




END SUBROUTINE InflowWind_CalcOutput



!====================================================================================================
SUBROUTINE InflowWind_End( InputData, ParamData, ContStates, DiscStates, ConstrStateGuess, OtherStates, &
                       OutData, ErrStat, ErrMsg )
   ! Clean up the allocated variables and close all open files.  Reset the initialization flag so
   ! that we have to reinitialize before calling the routines again.
   !----------------------------------------------------------------------------------------------------

         ! Initialization data and guesses

      TYPE(InflowWind_InputType),               INTENT(INOUT)  :: InputData         !< Input data for initialization
      TYPE(InflowWind_ParameterType),           INTENT(INOUT)  :: ParamData         !< Parameters
      TYPE(InflowWind_ContinuousStateType),     INTENT(INOUT)  :: ContStates        !< Continuous states
      TYPE(InflowWind_DiscreteStateType),       INTENT(INOUT)  :: DiscStates        !< Discrete states
      TYPE(InflowWind_ConstraintStateType),     INTENT(INOUT)  :: ConstrStateGuess  !< Guess of the constraint states
      TYPE(InflowWind_OtherStateType),          INTENT(INOUT)  :: OtherStates       !< Other/optimization states
      TYPE(InflowWind_OutputType),              INTENT(INOUT)  :: OutData           !< Output data


         ! Error Handling

      INTEGER( IntKi ),                         INTENT(  OUT)  :: ErrStat
      CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg

         ! Local variables

      TYPE(IfW_UniformWind_InputType)                               :: Uniform_InputData      !< input positions.
      TYPE(IfW_UniformWind_OutputType)                              :: Uniform_OutData        !< output velocities

!!!      TYPE(IfW_FFWind_InputType)                               :: FF_InputData      !< input positions.
!!!      TYPE(IfW_FFWind_ContinuousStateType)                     :: FF_ContStates     !< Unused
!!!      TYPE(IfW_FFWind_DiscreteStateType)                       :: FF_DiscStates     !< Unused
!!!      TYPE(IfW_FFWind_ConstraintStateType)                     :: FF_ConstrStates   !< Unused
!!!      TYPE(IfW_FFWind_OutputType)                              :: FF_OutData        !< output velocities


!!!     TYPE(CTTS_Backgr)                                           :: BackGrndValues


!NOTE: It isn't entirely clear what the purpose of Height is. Does it sometimes occur that Height  /= ParamData%ReferenceHeight???
!      REAL(ReKi)                                         :: Height      ! Retrieved from FF
!      REAL(ReKi)                                         :: HalfWidth   ! Retrieved from FF



         ! End the sub-modules (deallocates their arrays and closes their files):

      SELECT CASE ( ParamData%WindType )

         CASE (Steady_WindNumber)         ! The Steady wind is a simple wrapper for the UniformWind module.
            CALL IfW_UniformWind_End( InputData%PositionXYZ,   ParamData%UniformWind,                                        &
                                 OtherStates%UniformWind, Uniform_OutData, ErrStat, ErrMsg )

         CASE (Uniform_WindNumber)
            CALL IfW_UniformWind_End( InputData%PositionXYZ,   ParamData%UniformWind,                                        &
                                 OtherStates%UniformWind, Uniform_OutData, ErrStat, ErrMsg )

!!!         CASE (FF_WindNumber)
!!!            CALL IfW_FFWind_End( FF_InitData,   ParamData%FFWind,                                        &
!!!                                 FF_ContStates, FF_DiscStates,    FF_ConstrStates,  OtherStates%FFWind,  &
!!!                                 FF_OutData,    ErrStat,          ErrMsg )
!!!
!!!         CASE (User_WindNumber)
!!!            CALL UsrWnd_Terminate( ErrStat )
!!!
!!!         CASE (FD_WindNumber)
!!!            CALL FD_Terminate(     ErrStat )
!!!
!!!         CASE (HAWC_WindNumber)
!!!            CALL HW_Terminate(     ErrStat )

         CASE ( Undef_WindNumber )
            ! Do nothing

         CASE DEFAULT  ! keep this check to make sure that all new wind types have been accounted for
            CALL SetErrStat(ErrID_Fatal,' Undefined wind type.',ErrStat,ErrMsg,'InflowWind_End')

      END SELECT

!!!  !   IF (CTTS_Flag) CALL CTTS_Terminate( ErrStat ) !FIXME: should it be this line or the next?
!!!         CALL CTTS_Terminate( ErrStat, ErrMsg )


         ! Reset the wind type so that the initialization routine must be called
      ParamData%WindType = Undef_WindNumber
      ParamData%Initialized = .FALSE.
      ParamData%CTTS_Flag  = .FALSE.


END SUBROUTINE InflowWind_End



!====================================================================================================
! The following routines were added to satisfy the framework, but do nothing useful.
!====================================================================================================
SUBROUTINE InflowWind_UpdateStates( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input Time; Continuous and discrete states are updated for Time + Interval
!..................................................................................................................................

      REAL(DbKi),                               INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(InflowWind_InputType),               INTENT(IN   )  :: u           !< Inputs at Time
      TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: p           !< Parameters
      TYPE(InflowWind_ContinuousStateType),     INTENT(INOUT)  :: x           !< Input: Continuous states at Time;
                                                                              !! Output: Continuous states at Time + Interval
      TYPE(InflowWind_DiscreteStateType),       INTENT(INOUT)  :: xd          !< Input: Discrete states at Time;
                                                                              !! Output: Discrete states at Time  + Interval
      TYPE(InflowWind_ConstraintStateType),     INTENT(INOUT)  :: z           !< Input: Initial guess of constraint states at Time;
                                                                              !! Output: Constraint states at Time
      TYPE(InflowWind_OtherStateType),          INTENT(INOUT)  :: OtherState  !< Other/optimization states
      INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None

         ! Local variables

      TYPE(InflowWind_ContinuousStateType)                     :: dxdt        !< Continuous state derivatives at Time
      TYPE(InflowWind_ConstraintStateType)                     :: z_Residual  !< Residual of the constraint state equations (Z)

      INTEGER(IntKi)                                           :: TmpErrStat    !< Error status of the operation (occurs after initial error)
      CHARACTER(LEN(ErrMsg))                                   :: TmpErrMsg     !< Error message if TmpErrStat /= ErrID_None

         ! Initialize ErrStat


      ! BJJ: Please don't make my code end just because I called a routine that you don't use :)
      ErrStat = ErrID_None
      ErrMsg  = ""

      RETURN


      ErrStat = ErrID_Warn
      ErrMsg  = "IfW_UpdateStates was called.  That routine does nothing useful."



         ! Solve for the constraint states (z) here:

         ! Check if the z guess is correct and update z with a new guess.
         ! Iterate until the value is within a given tolerance.

      CALL IfW_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_Residual, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL IfW_DestroyConstrState( z_Residual, TmpErrStat, TmpErrMsg)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
         RETURN
      ENDIF

      ! DO WHILE ( z_Residual% > tolerance )
      !
      !  z =
      !
      !  CALL IfW_FFWind_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_Residual, ErrStat, ErrMsg )
      !  IF ( ErrStat >= AbortErrLev ) THEN
      !     CALL IfW_FFWind_DestroyConstrState( z_Residual, TmpErrStat, TmpErrMsg)
      !     ErrMsg = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
      !     RETURN
      !  ENDIF
      !
      ! END DO


         ! Destroy z_Residual because it is not necessary for the rest of the subroutine:

      CALL IfW_DestroyConstrState( z_Residual, ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN



         ! Get first time derivatives of continuous states (dxdt):

      CALL IfW_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL IfW_DestroyContState( dxdt, TmpErrStat, TmpErrMsg)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
         RETURN
      ENDIF


         ! Update discrete states:
         !   Note that xd [discrete state] is changed in IfW_FFWind_UpdateDiscState(), so IfW_FFWind_CalcOutput(),
         !   IfW_FFWind_CalcContStateDeriv(), and IfW_FFWind_CalcConstrStates() must be called first (see above).

      CALL IfW_UpdateDiscState(Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL IfW_DestroyContState( dxdt, TmpErrStat, TmpErrMsg)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(TmpErrMsg)
         RETURN
      ENDIF


         ! Integrate (update) continuous states (x) here:

      !x = function of dxdt and x


         ! Destroy dxdt because it is not necessary for the rest of the subroutine

      CALL IfW_DestroyContState( dxdt, ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN



END SUBROUTINE InflowWind_UpdateStates



!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE InflowWind_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................

      REAL(DbKi),                               INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(InflowWind_InputType),               INTENT(IN   )  :: u           !< Inputs at Time
      TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: p           !< Parameters
      TYPE(InflowWind_ContinuousStateType),     INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(InflowWind_DiscreteStateType),       INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(InflowWind_ConstraintStateType),     INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(InflowWind_OtherStateType),          INTENT(INOUT)  :: OtherState  !< Other/optimization states
      TYPE(InflowWind_ContinuousStateType),     INTENT(  OUT)  :: dxdt        !< Continuous state derivatives at Time
      INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Compute the first time derivatives of the continuous states here:

      dxdt%DummyContState = 0.0_ReKi


END SUBROUTINE InflowWind_CalcContStateDeriv



!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE InflowWind_UpdateDiscState( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Tight coupling routine for updating discrete states
!..................................................................................................................................

      REAL(DbKi),                               INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(InflowWind_InputType),               INTENT(IN   )  :: u           !< Inputs at Time
      TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: p           !< Parameters
      TYPE(InflowWind_ContinuousStateType),     INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(InflowWind_DiscreteStateType),       INTENT(INOUT)  :: xd          !< Input: Discrete states at Time;
                                                                              !! Output: Discrete states at Time + Interval
      TYPE(InflowWind_ConstraintStateType),     INTENT(IN   )  :: z           !< Constraint states at Time
      TYPE(InflowWind_OtherStateType),          INTENT(INOUT)  :: OtherState  !< Other/optimization states
      INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Update discrete states here:

      ! StateData%DiscState =

END SUBROUTINE InflowWind_UpdateDiscState



!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE InflowWind_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMsg )
! Tight coupling routine for solving for the residual of the constraint state equations
!..................................................................................................................................

      REAL(DbKi),                               INTENT(IN   )  :: Time        !< Current simulation time in seconds
      TYPE(InflowWind_InputType),               INTENT(IN   )  :: u           !< Inputs at Time
      TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: p           !< Parameters
      TYPE(InflowWind_ContinuousStateType),     INTENT(IN   )  :: x           !< Continuous states at Time
      TYPE(InflowWind_DiscreteStateType),       INTENT(IN   )  :: xd          !< Discrete states at Time
      TYPE(InflowWind_ConstraintStateType),     INTENT(IN   )  :: z           !< Constraint states at Time (possibly a guess)
      TYPE(InflowWind_OtherStateType),          INTENT(INOUT)  :: OtherState  !< Other/optimization states
      TYPE(InflowWind_ConstraintStateType),     INTENT(  OUT)  :: z_residual  !< Residual of the constraint state equations using
                                                                              !! the input values described above
      INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat     !< Error status of the operation
      CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Solve for the constraint states here:

      z_residual%DummyConstrState = 0

END SUBROUTINE InflowWind_CalcConstrStateResidual
!====================================================================================================
!====================================================================================================
FUNCTION WindInf_ADhack_diskVel( Time,ParamData, OtherStates,ErrStat, ErrMsg )
! This function should be deleted ASAP.  It's purpose is to reproduce results of AeroDyn 12.57;
! when a consensus on the definition of "average velocity" is determined, this function will be
! removed.
!----------------------------------------------------------------------------------------------------

      ! Passed variables

   REAL(DbKi),                                  INTENT(IN   )  :: Time              !< Time
   TYPE(InflowWind_ParameterType),              INTENT(IN   )  :: ParamData         !< Parameters
   TYPE(InflowWind_OtherStateType),             INTENT(INOUT)  :: OtherStates       !< Other/optimization states

   INTEGER(IntKi),                              INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                                INTENT(  OUT)  :: ErrMsg

      ! Function definition
   REAL(ReKi)                    :: WindInf_ADhack_diskVel(3)

      ! Local variables
   REAL(ReKi)                    :: Delta_tmp            ! interpolated Delta   at input TIME
   REAL(ReKi)                    :: P                    ! temporary storage for slope (in time) used in linear interpolation
   REAL(ReKi)                    :: V_tmp                ! interpolated V       at input TIME
   REAL(ReKi)                    :: VZ_tmp               ! interpolated VZ      at input TIME



   ErrStat = ErrID_None

   SELECT CASE ( ParamData%WindType )
!!!      CASE (Uniform_WindNumber)
!!!
!!!         !-------------------------------------------------------------------------------------------------
!!!         ! Linearly interpolate in time (or use nearest-neighbor to extrapolate)
!!!         ! (compare with NWTC_Num.f90\InterpStpReal)
!!!         !-------------------------------------------------------------------------------------------------
!!!
!!!
!!!            ! Let's check the limits.
!!!         IF ( Time <= OtherStates%UniformWind%Tdata(1) .OR. OtherStates%UniformWind%NumDataLines == 1 )  THEN
!!!
!!!            OtherStates%UniformWind%TimeIndex      = 1
!!!            V_tmp         = OtherStates%UniformWind%V      (1)
!!!            Delta_tmp     = OtherStates%UniformWind%Delta  (1)
!!!            VZ_tmp        = OtherStates%UniformWind%VZ     (1)
!!!
!!!         ELSE IF ( Time >= OtherStates%UniformWind%Tdata(OtherStates%UniformWind%NumDataLines) )  THEN
!!!
!!!            OtherStates%UniformWind%TimeIndex = OtherStates%UniformWind%NumDataLines - 1
!!!            V_tmp                 = OtherStates%UniformWind%V      (OtherStates%UniformWind%NumDataLines)
!!!            Delta_tmp             = OtherStates%UniformWind%Delta  (OtherStates%UniformWind%NumDataLines)
!!!            VZ_tmp                = OtherStates%UniformWind%VZ     (OtherStates%UniformWind%NumDataLines)
!!!
!!!         ELSE
!!!
!!!              ! Let's interpolate!
!!!
!!!            OtherStates%UniformWind%TimeIndex = MAX( MIN( OtherStates%UniformWind%TimeIndex, OtherStates%UniformWind%NumDataLines-1 ), 1 )
!!!
!!!            DO
!!!
!!!               IF ( Time < OtherStates%UniformWind%Tdata(OtherStates%UniformWind%TimeIndex) )  THEN
!!!
!!!                  OtherStates%UniformWind%TimeIndex = OtherStates%UniformWind%TimeIndex - 1
!!!
!!!               ELSE IF ( Time >= OtherStates%UniformWind%Tdata(OtherStates%UniformWind%TimeIndex+1) )  THEN
!!!
!!!                  OtherStates%UniformWind%TimeIndex = OtherStates%UniformWind%TimeIndex + 1
!!!
!!!               ELSE
!!!                  P           =  ( Time - OtherStates%UniformWind%Tdata(OtherStates%UniformWind%TimeIndex) )/     &
!!!                                 ( OtherStates%UniformWind%Tdata(OtherStates%UniformWind%TimeIndex+1)             &
!!!                                 - OtherStates%UniformWind%Tdata(OtherStates%UniformWind%TimeIndex) )
!!!                  V_tmp       =  ( OtherStates%UniformWind%V(      OtherStates%UniformWind%TimeIndex+1)           &
!!!                                 - OtherStates%UniformWind%V(      OtherStates%UniformWind%TimeIndex) )*P         &
!!!                                 + OtherStates%UniformWind%V(      OtherStates%UniformWind%TimeIndex)
!!!                  Delta_tmp   =  ( OtherStates%UniformWind%Delta(  OtherStates%UniformWind%TimeIndex+1)           &
!!!                                 - OtherStates%UniformWind%Delta(  OtherStates%UniformWind%TimeIndex) )*P         &
!!!                                 + OtherStates%UniformWind%Delta(  OtherStates%UniformWind%TimeIndex)
!!!                  VZ_tmp      =  ( OtherStates%UniformWind%VZ(     OtherStates%UniformWind%TimeIndex+1)           &
!!!                                 - OtherStates%UniformWind%VZ(     OtherStates%UniformWind%TimeIndex) )*P  &
!!!                                 + OtherStates%UniformWind%VZ(     OtherStates%UniformWind%TimeIndex)
!!!                  EXIT
!!!
!!!               END IF
!!!
!!!            END DO
!!!
!!!         END IF
!!!
!!!      !-------------------------------------------------------------------------------------------------
!!!      ! calculate the wind speed at this time
!!!      !-------------------------------------------------------------------------------------------------
!!!
!!!         WindInf_ADhack_diskVel(1) =  V_tmp * COS( Delta_tmp )
!!!         WindInf_ADhack_diskVel(2) = -V_tmp * SIN( Delta_tmp )
!!!         WindInf_ADhack_diskVel(3) =  VZ_tmp
!!!
!!!
!!!
!!!      CASE (FF_WindNumber)
!!!
!!!         WindInf_ADhack_diskVel(1)   = OtherStates%FFWind%MeanFFWS
!!!         WindInf_ADhack_diskVel(2:3) = 0.0

      CASE DEFAULT
         ErrStat = ErrID_Fatal
         ErrMsg = ' WindInf_ADhack_diskVel: Undefined wind type.'

   END SELECT

   RETURN

END FUNCTION WindInf_ADhack_diskVel


!====================================================================================================
END MODULE InflowWind
