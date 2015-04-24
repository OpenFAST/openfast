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
      ! The included wind modules (TYPES modules are inherited from InflowWind_Types, so not specified here again.)
      !-------------------------------------------------------------------------------------------------
   USE                              Lidar                      ! module for obtaining sensor data
      
   USE                              IfW_UniformWind            ! uniform wind files (text files)
   USE                              IfW_TSFFWind               ! TurbSim style full-field binary wind files
   USE                              IfW_BladedFFWind           ! Bladed style full-field binary wind files
   USE                              IfW_UserWind               ! User-defined wind module

!!!   USE                              HAWCWind                   ! full-field binary wind files in HAWC format
!!!   USE                              FDWind                     ! 4-D binary wind files
!!!   USE                              CTWind                     ! coherent turbulence from KH billow - binary file superimposed on another wind type

   
   IMPLICIT NONE
   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: IfW_Ver = ProgDesc( 'InflowWind', 'v3.00.01a-adp', '01-Nov-2014' )



      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: InflowWind_Init                                   !< Initialization routine
   PUBLIC :: InflowWind_CalcOutput                             !< Calculate the wind velocities
   PUBLIC :: InflowWind_End                                    !< Ending routine (includes clean up)


      ! These routines satisfy the framework, but do nothing at present.
   PUBLIC :: InflowWind_UpdateStates               !< Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
   PUBLIC :: InflowWind_CalcConstrStateResidual    !< Tight coupling routine for returning the constraint state residual
   PUBLIC :: InflowWind_CalcContStateDeriv         !< Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: InflowWind_UpdateDiscState            !< Tight coupling routine for updating discrete states




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

      IMPLICIT                                                 NONE

      CHARACTER(*),              PARAMETER                  :: RoutineName="InflowWind_Init"

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

      TYPE(IfW_UniformWind_InitInputType)                   :: Uniform_InitData  !< initialization info
      TYPE(IfW_UniformWind_OutputType)                      :: Uniform_OutData   !< output velocities

      TYPE(IfW_TSFFWind_InitInputType)                      :: TSFF_InitData     !< initialization info
      TYPE(IfW_TSFFWind_OutputType)                         :: TSFF_OutData      !< output velocities


      TYPE(IfW_BladedFFWind_InitInputType)                  :: BladedFF_InitData       !< initialization info
      TYPE(IfW_BladedFFWind_OutputType)                     :: BladedFF_OutData        !< output velocities

      TYPE(IfW_UserWind_InitInputType)                      :: User_InitData     !< initialization info
      TYPE(IfW_UserWind_OutputType)                         :: User_OutData      !< output velocities


!!!     TYPE(CTBladed_Backgr)                                        :: BackGrndValues


         ! Temporary variables for error handling
      INTEGER(IntKi)                                        :: TmpErrStat
      CHARACTER(LEN(ErrMsg))                                :: TmpErrMsg         !< temporary error message


         ! Local Variables
      INTEGER(IntKi)                                        :: I                 !< Generic counter
      INTEGER(IntKi)                                        :: SumFileUnit       !< Unit number for the summary file
      CHARACTER(256)                                        :: SumFileName       !< Name of the summary file
      CHARACTER(256)                                        :: EchoFileName      !< Name of the summary file


         !----------------------------------------------------------------------------------------------
         ! Initialize variables and check to see if this module has been initialized before.
         !----------------------------------------------------------------------------------------------

      ErrStat        =  ErrID_None
      ErrMsg         =  ""


         ! Set a few variables.

      ParamData%DT   = TimeInterval             ! InflowWind does not require a specific time interval, so this is never changed.
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
      EchoFileName  = TRIM(ParamData%RootFileName)//".IfW.ech"
      SumFileName   = TRIM(ParamData%RootFileName)//".IfW.sum"


         ! Parse all the InflowWind related input files and populate the *_InitDataType derived types

      IF ( InitData%UseInputFile ) THEN
         CALL InflowWind_ReadInput( ParamData%InputFileName, EchoFileName, InputFileData, TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF
         
         ! these values (and others hard-coded in lidar_init) should be set in the input file, too
        InputFileData%SensorType = InitData%lidar%SensorType
        InputFileData%NumPulseGate = InitData%lidar%NumPulseGate
        InputFileData%RotorApexOffsetPos = InitData%lidar%RotorApexOffsetPos
        InputFileData%LidRadialVel = InitData%lidar%LidRadialVel
                        
      ELSE
            !  In the future we will make it possible to just copy the data from the InitData%PassedFileData (of derived type
            !  InflowWind_InputFile), and run things as normal using that data.  For now though, we will not allow this.
         CALL SetErrStat(ErrID_Fatal,' The UseInputFile flag cannot be set to .FALSE. at present.  This feature is not '// &
               'presently supported.',ErrStat,ErrMsg,RoutineName)
         CALL Cleanup()
         RETURN
      ENDIF


         ! initialize sensor data:   
      CALL Lidar_Init( InitData,   InputGuess,    ParamData,                          &
                       ContStates, DiscStates,    ConstrStateGuess,    OtherStates,   &
                       OutData,    TimeInterval,  InitOutData,  TmpErrStat, TmpErrMsg )
                  CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )      
      

         ! Validate the InflowWind input file information.

      CALL InflowWind_ValidateInput( InputFileData, TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat>= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      ENDIF


         ! Set the ParamData and OtherStates for InflowWind using the input file information.

      CALL InflowWind_SetParameters( InputFileData, ParamData, OtherStates, TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat>= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      ENDIF

 

         ! Allocate arrays for the WriteOutput

      CALL AllocAry( OutData%WriteOutput, ParamData%NumOuts, 'WriteOutput', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat>= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      ENDIF
      OutData%WriteOutput = 0.0_ReKi
      
      CALL AllocAry( InitOutData%WriteOutputHdr, ParamData%NumOuts, 'WriteOutputHdr', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat>= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      ENDIF

      CALL AllocAry( InitOutData%WriteOutputUnt, ParamData%NumOuts, 'WriteOutputUnt', TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat>= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      ENDIF
   
      InitOutData%WriteOutputHdr = ParamData%OutParam(1:ParamData%NumOuts)%Name
      InitOutData%WriteOutputUnt = ParamData%OutParam(1:ParamData%NumOuts)%Units     
 


         ! If a summary file was requested, open it.
      IF ( InputFileData%SumPrint ) THEN

            ! Open the summary file and write some preliminary info to it
         CALL InflowWind_OpenSumFile( SumFileUnit, SumFileName, IfW_Ver, ParamData%WindType, TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         ENDIF
      ELSE
         SumFileUnit =  -1_IntKi       ! So that we don't try to write to something.  Used as indicator in submodules.
      ENDIF


         ! Allocate the arrays for passing points in and velocities out
      IF ( .NOT. ALLOCATED(InputGuess%PositionXYZ) ) THEN
         CALL AllocAry( InputGuess%PositionXYZ, 3, InitData%NumWindPoints, &
                     "Array of positions at which to find wind velocities", TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat>= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF
      ENDIF
      IF ( .NOT. ALLOCATED(OutData%VelocityUVW) ) THEN
         CALL AllocAry( OutData%VelocityUVW, 3, InitData%NumWindPoints, &
                     "Array of wind velocities returned by InflowWind", TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat>= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF
      ENDIF


         ! Check that the arrays for the points and velocities are the same size
      IF ( SIZE( InputGuess%PositionXYZ, DIM = 2 ) /= SIZE( OutData%VelocityUVW, DIM = 2 ) ) THEN
         CALL SetErrStat(ErrID_Fatal,' Programming error: Different number of XYZ coordinates and expected output velocities.', &
                     ErrStat,ErrMsg,RoutineName)
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
            ParamData%UniformWind%NumDataLines     =  1_IntKi
            ParamData%UniformWind%RefHt            =  InputFileData%Steady_RefHt
            ParamData%UniformWind%RefLength        =  1.0_ReKi    ! This is not used since no shear gusts are used.  Set to 1.0 so calculations don't bomb. 
            OtherStates%UniformWind%TimeIndex      =  1           ! Used in UniformWind as a check if initialization was done.


            IF (.NOT. ALLOCATED(ParamData%UniformWind%Tdata) ) THEN
               CALL AllocAry( ParamData%UniformWind%Tdata, ParamData%UniformWind%NumDataLines, 'Uniform wind time', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF
         
            IF (.NOT. ALLOCATED(ParamData%UniformWind%V) ) THEN
               CALL AllocAry( ParamData%UniformWind%V, ParamData%UniformWind%NumDataLines, 'Uniform wind horizontal wind speed', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF
         
            IF (.NOT. ALLOCATED(ParamData%UniformWind%Delta) ) THEN
               CALL AllocAry( ParamData%UniformWind%Delta, ParamData%UniformWind%NumDataLines, 'Uniform wind direction', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF
         
            IF (.NOT. ALLOCATED(ParamData%UniformWind%VZ) ) THEN
               CALL AllocAry( ParamData%UniformWind%VZ, ParamData%UniformWind%NumDataLines, 'Uniform vertical wind speed', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF
         
            IF (.NOT. ALLOCATED(ParamData%UniformWind%HShr) ) THEN
               CALL AllocAry( ParamData%UniformWind%HShr, ParamData%UniformWind%NumDataLines, 'Uniform horizontal linear shear', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF
         
            IF (.NOT. ALLOCATED(ParamData%UniformWind%VShr) ) THEN
               CALL AllocAry( ParamData%UniformWind%VShr, ParamData%UniformWind%NumDataLines, 'Uniform vertical power-law shear exponent', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF
         
            IF (.NOT. ALLOCATED(ParamData%UniformWind%VLinShr) ) THEN
               CALL AllocAry( ParamData%UniformWind%VLinShr, ParamData%UniformWind%NumDataLines, 'Uniform vertical linear shear', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF
         
            IF (.NOT. ALLOCATED(ParamData%UniformWind%VGust) ) THEN
               CALL AllocAry( ParamData%UniformWind%VGust, ParamData%UniformWind%NumDataLines, 'Uniform gust velocity', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
               IF ( ErrStat >= AbortErrLev ) RETURN
            END IF


               ! Set the array information
            
            ParamData%UniformWind%Tdata(  :)       = 0.0_ReKi
            ParamData%UniformWind%V(      :)       = InputFileData%Steady_HWindSpeed
            ParamData%UniformWind%Delta(  :)       = 0.0_ReKi
            ParamData%UniformWind%VZ(     :)       = 0.0_ReKi
            ParamData%UniformWind%HShr(   :)       = 0.0_ReKi
            ParamData%UniformWind%VShr(   :)       = InputFileData%Steady_PLexp
            ParamData%UniformWind%VLinShr(:)       = 0.0_ReKi
            ParamData%UniformWind%VGust(  :)       = 0.0_ReKi



               ! Now we have in effect initialized the IfW_UniformWind module, so set the parameters
            ParamData%UniformWind%RefLength        =  1.0_ReKi       ! This is done so that 0.0*(some stuff)/RefLength doesn't blow up.
            ParamData%UniformWind%RefHt            =  InputFileData%Steady_RefHt
            OtherStates%UniformWind%TimeIndex      =  1_IntKi


               ! Store wind file metadata
            InitOutData%WindFileInfo%FileName         =  ""
            InitOutData%WindFileInfo%WindType         =  Steady_WindNumber
            InitOutData%WindFileInfo%RefHt            =  InputFileData%Steady_RefHt
            InitOutData%WindFileInfo%RefHt_Set        =  .FALSE.                             ! The wind file does not set this
            InitOutData%WindFileInfo%DT               =  0.0_ReKi
            InitOutData%WindFileInfo%NumTSteps        =  1_IntKi
            InitOutData%WindFileInfo%ConstantDT       =  .FALSE.
            InitOutData%WindFileInfo%TRange           =  (/ 0.0_ReKi, 0.0_ReKi /)
            InitOutData%WindFileInfo%TRange_Limited   =  .FALSE.                             ! This is constant
            InitOutData%WindFileInfo%YRange           =  (/ 0.0_ReKi, 0.0_ReKi /)
            InitOutData%WindFileInfo%YRange_Limited   =  .FALSE.                             ! Hard boundaries not enforced in y-direction
            InitOutData%WindFileInfo%ZRange           =  (/ 0.0_ReKi, 0.0_ReKi /)
            InitOutData%WindFileInfo%ZRange_Limited   =  .FALSE.                             ! Hard boundaries not enforced in z-direction
            InitOutData%WindFileInfo%BinaryFormat     =  0_IntKi
            InitOutData%WindFileInfo%IsBinary         =  .FALSE.
            InitOutData%WindFileInfo%TI               =  0.0_ReKi
            InitOutData%WindFileInfo%TI_listed        =  .FALSE.
      

               ! Write summary file information
            IF ( SumFileUnit > 0 ) THEN
               WRITE(SumFileUnit,'(A)',IOSTAT=TmpErrStat)
               WRITE(SumFileUnit,'(A80)',IOSTAT=TmpErrStat)          'Steady wind -- Constant wind profile for entire simulation. No windfile read in.'
               WRITE(SumFileUnit,'(A40,G12.4)',IOSTAT=TmpErrStat)    '     Reference height:                  ',ParamData%UniformWind%RefHt
               WRITE(SumFileUnit,'(A40,G12.4)',IOSTAT=TmpErrStat)    '     Horizontal velocity:               ',ParamData%UniformWind%V
               WRITE(SumFileUnit,'(A40,G12.4)',IOSTAT=TmpErrStat)    '     Vertical sheer power law exponent: ',ParamData%UniformWind%VShr

                  ! We are assuming that if the last line was written ok, then all of them were.
               IF (TmpErrStat /= 0_IntKi) THEN
                  CALL SetErrStat(ErrID_Fatal,'Error writing to summary file.',ErrStat,ErrMsg,RoutineName)
                  CALL Cleanup
                  RETURN
               ENDIF   
            ENDIF 


         CASE ( Uniform_WindNumber )

               ! Set InitData information
            Uniform_InitData%ReferenceHeight          =  InputFileData%Uniform_RefHt
            Uniform_InitData%RefLength                =  InputFileData%Uniform_RefLength 
            Uniform_InitData%WindFileName             =  InputFileData%Uniform_FileName
            Uniform_InitData%SumFileUnit              =  SumFileUnit

               ! Initialize the UniformWind module
            CALL IfW_UniformWind_Init(Uniform_InitData, InputGuess%PositionXYZ, ParamData%UniformWind, OtherStates%UniformWind, &
                        Uniform_OutData,    TimeInterval,  InitOutData%UniformWind,  TmpErrStat,          TmpErrMsg)

            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_Init' )
            IF ( ErrStat >= AbortErrLev ) RETURN


               ! Store wind file metadata
            InitOutData%WindFileInfo%FileName         =  InputFileData%Uniform_FileName
            InitOutData%WindFileInfo%WindType         =  Uniform_WindNumber
            InitOutData%WindFileInfo%RefHt            =  ParamData%UniformWind%RefHt
            InitOutData%WindFileInfo%RefHt_Set        =  .FALSE.                             ! The wind file does not set this
            InitOutData%WindFileInfo%DT               =  InitOutData%UniformWind%WindFileDT
            InitOutData%WindFileInfo%NumTSteps        =  InitOutData%UniformWind%WindFileNumTSteps
            InitOutData%WindFileInfo%ConstantDT       =  InitOutData%UniformWind%WindFileConstantDT
            InitOutData%WindFileInfo%TRange           =  InitOutData%UniformWind%WindFileTRange
            InitOutData%WindFileInfo%TRange_Limited   = .FALSE.                              ! UniformWind sets to limit of file if outside time bounds
            InitOutData%WindFileInfo%YRange           =  (/ 0.0_ReKi, 0.0_ReKi /)
            InitOutData%WindFileInfo%YRange_Limited   =  .FALSE.                             ! Hard boundaries not enforced in y-direction
            InitOutData%WindFileInfo%ZRange           =  (/ 0.0_ReKi, 0.0_ReKi /)
            InitOutData%WindFileInfo%ZRange_Limited   =  .FALSE.                             ! Hard boundaries not enforced in z-direction
            InitOutData%WindFileInfo%BinaryFormat     =  0_IntKi
            InitOutData%WindFileInfo%IsBinary         =  .FALSE.
            InitOutData%WindFileInfo%TI               =  0.0_ReKi
            InitOutData%WindFileInfo%TI_listed        =  .FALSE.


               ! Check if the fist data point from the file is not along the X-axis while applying the windfield rotation
            IF ( ( .NOT. EqualRealNos (ParamData%UniformWind%Delta(1), 0.0_ReKi) ) .AND.  &
                 ( .NOT. EqualRealNos (ParamData%PropogationDir, 0.0_ReKi)       ) ) THEN
               CALL SetErrStat( ErrID_Warn,' Possible double rotation of wind field! Uniform wind file starts with a wind direction of '// &
                        TRIM(Num2LStr(ParamData%UniformWind%Delta(1)*R2D))//                       &
                        ' degrees and the InflowWind input file specifies a PropogationDir of '//  &
                        TRIM(Num2LStr(ParamData%PropogationDir*R2D))//' degrees.',                 &
                        ErrStat,ErrMsg,RoutineName )
            ENDIF



         CASE ( TSFF_WindNumber )

               ! Set InitData information
            TSFF_InitData%WindFileName                   =  InputFileData%TSFF_FileName
            TSFF_InitData%SumFileUnit                    =  SumFileUnit

               ! Initialize the TSFFWind module
            CALL IfW_TSFFWind_Init(TSFF_InitData, InputGuess%PositionXYZ, ParamData%TSFFWind, OtherStates%TSFFWind, &
                        TSFF_OutData,    TimeInterval,  InitOutData%TSFFWind,  TmpErrStat,          TmpErrMsg)
            CALL SetErrSTat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN

               ! Store wind file metadata
            InitOutData%WindFileInfo%FileName            =  InputFileData%TSFF_FileName
            InitOutData%WindFileInfo%WindType            =  TSFF_WindNumber
            InitOutData%WindFileInfo%RefHt               =  ParamData%TSFFWind%RefHt
            InitOutData%WindFileInfo%RefHt_Set           =  .TRUE.
            InitOutData%WindFileInfo%DT                  =  ParamData%TSFFWind%FFDTime
            InitOutData%WindFileInfo%NumTSteps           =  ParamData%TSFFWind%NFFSteps
            InitOutData%WindFileInfo%ConstantDT          =  .TRUE.
            IF ( ParamData%TSFFWind%Periodic ) THEN
               InitOutData%WindFileInfo%TRange           =  (/ 0.0_ReKi, ParamData%TSFFWind%TotalTime /)
               InitOutData%WindFileInfo%TRange_Limited   =  .FALSE.
            ELSE  ! Shift the time range to compensate for the shifting of the wind grid
               InitOutData%WindFileInfo%TRange           =  (/ 0.0_ReKi, ParamData%TSFFWind%TotalTime /)  &
                  -  ParamData%TSFFWind%InitXPosition*ParamData%TSFFWind%InvMFFWS
               InitOutData%WindFileInfo%TRange_Limited   =  .TRUE.
            ENDIF
            InitOutData%WindFileInfo%YRange              =  (/ -ParamData%TSFFWind%FFYHWid, ParamData%TSFFWind%FFYHWid /)
            InitOutData%WindFileInfo%YRange_Limited      =  .TRUE.      ! Hard boundaries enforced in y-direction
            IF ( ParamData%TSFFWind%TowerDataExist ) THEN        ! have tower data
               InitOutData%WindFileInfo%ZRange           =  (/ 0.0_Reki,                                                         & 
                                                               ParamData%TSFFWind%RefHt + ParamData%TSFFWind%FFZHWid /)
            ELSE
               InitOutData%WindFileInfo%ZRange           =  (/ ParamData%TSFFWind%RefHt - ParamData%TSFFWind%FFZHWid,    &
                                                               ParamData%TSFFWind%RefHt + ParamData%TSFFWind%FFZHWid /)
            ENDIF
            InitOutData%WindFileInfo%ZRange_Limited      =  .TRUE.
            InitOutData%WindFileInfo%BinaryFormat        =  ParamData%TSFFWind%WindFileFormat
            InitOutData%WindFileInfo%IsBinary            =  .TRUE.
            InitOutData%WindFileInfo%TI                  =  0.0_ReKi
            InitOutData%WindFileInfo%TI_listed           =  .FALSE.



         CASE ( BladedFF_WindNumber )

               ! Set InitData information
            BladedFF_InitData%WindFileName               =  InputFileData%BladedFF_FileName
            BladedFF_InitData%TowerFileExist             =  InputFileData%BladedFF_TowerFile
            BladedFF_InitData%SumFileUnit                =  SumFileUnit

               ! Initialize the BladedFFWind module
            CALL IfW_BladedFFWind_Init(BladedFF_InitData, InputGuess%PositionXYZ, ParamData%BladedFFWind, OtherStates%BladedFFWind, &
                        BladedFF_OutData,    TimeInterval,  InitOutData%BladedFFWind,  TmpErrStat,          TmpErrMsg)
            CALL SetErrSTat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN

               ! Store wind file metadata
            InitOutData%WindFileInfo%FileName            =  InputFileData%BladedFF_FileName
            InitOutData%WindFileInfo%WindType            =  BladedFF_WindNumber
            InitOutData%WindFileInfo%RefHt               =  ParamData%BladedFFWind%RefHt
            InitOutData%WindFileInfo%RefHt_Set           =  .TRUE.
            InitOutData%WindFileInfo%DT                  =  ParamData%BladedFFWind%FFDTime
            InitOutData%WindFileInfo%NumTSteps           =  ParamData%BladedFFWind%NFFSteps
            InitOutData%WindFileInfo%ConstantDT          =  .TRUE.
            IF ( ParamData%BladedFFWind%Periodic ) THEN
               InitOutData%WindFileInfo%TRange           =  (/ 0.0_ReKi, ParamData%BladedFFWind%TotalTime /)
               InitOutData%WindFileInfo%TRange_Limited   =  .FALSE.
            ELSE  ! Shift the time range to compensate for the shifting of the wind grid
               InitOutData%WindFileInfo%TRange           =  (/ 0.0_ReKi, ParamData%BladedFFWind%TotalTime /)  &
                  -  ParamData%BladedFFWind%InitXPosition*ParamData%BladedFFWind%InvMFFWS
               InitOutData%WindFileInfo%TRange_Limited   =  .TRUE.
            ENDIF
            InitOutData%WindFileInfo%YRange              =  (/ -ParamData%BladedFFWind%FFYHWid, ParamData%BladedFFWind%FFYHWid /)
            InitOutData%WindFileInfo%YRange_Limited      =  .TRUE.      ! Hard boundaries enforced in y-direction
            IF ( ParamData%BladedFFWind%TowerDataExist ) THEN        ! have tower data
               InitOutData%WindFileInfo%ZRange           =  (/ 0.0_Reki,                                                         & 
                                                            ParamData%BladedFFWind%RefHt + ParamData%BladedFFWind%FFZHWid /)
            ELSE
               InitOutData%WindFileInfo%ZRange           =  (/ ParamData%BladedFFWind%RefHt - ParamData%BladedFFWind%FFZHWid,    &
                                                            ParamData%BladedFFWind%RefHt + ParamData%BladedFFWind%FFZHWid /)
            ENDIF
            InitOutData%WindFileInfo%ZRange_Limited      =  .TRUE.
            InitOutData%WindFileInfo%BinaryFormat        =  ParamData%BladedFFWind%WindFileFormat
            InitOutData%WindFileInfo%IsBinary            =  .TRUE.
            InitOutData%WindFileInfo%TI                  =  InitOutData%BladedFFWind%TI
            InitOutData%WindFileInfo%TI_listed           =  .TRUE.   ! This must be listed in the file someplace



         CASE ( HAWC_WindNumber )
               ! Initialize the HAWCWind module
            CALL SetErrStat( ErrID_Fatal,' HAWC winds not supported yet.',ErrStat,ErrMsg,RoutineName )




         CASE (User_WindNumber)

               ! Initialize the UserWind module
            CALL IfW_UserWind_Init(User_InitData, InputGuess%PositionXYZ, ParamData%UserWind, OtherStates%UserWind, &
                        User_OutData,    TimeInterval,  InitOutData%UserWind,  TmpErrStat,          TmpErrMsg)
            CALL SetErrSTat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
            IF ( ErrStat >= AbortErrLev ) RETURN



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
!!!NOTE: is the CTTS used on uniform wind?
!!!!           IF (CTTS_Flag) CALL CTTS_SetRefVal(FileInfo%ReferenceHeight, 0.5*FileInfo%Width, ErrStat)  !FIXME: check if this was originally used
!!!!           IF (ErrStat == ErrID_None .AND. ParamData%CTTS_Flag) &
!!!!              CALL CTTS_SetRefVal(InitData%ReferenceHeight, REAL(0.0, ReKi), ErrStat, ErrMsg)      !FIXME: will need to put this routine in the Init of CT
!!!
!!!
!!!



!!!         CASE (HAWC_WindNumber)
!!!
!!!               !FIXME: remove this error message when we add HAWC_Wind in
!!!            CALL SetErrStat( ErrID_Fatal, ' InflowWind cannot currently handle the HAWC_Wind type.', ErrStat, ErrMsg, ' IfW_Init' )
!!!            RETURN
!!!
!!!!           CALL HW_Init( UnWind, ParamData%WindFileName, ErrStat )
!!!



         ! If we've arrived here, we haven't reached an AbortErrLev:
      ParamData%Initialized = .TRUE.

         ! Set the version information in InitOutData
      InitOutData%Ver   = IfW_Ver


         ! Close the summary file
      CALL InflowWind_CloseSumFile( SumFileUnit, TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)


      RETURN


   !----------------------------------------------------------------------------------------------------
CONTAINS

   SUBROUTINE CleanUp()

      ! add in stuff that we need to dispose of here
      CALL InflowWind_DestroyInputFile( InputFileData, TmpErrsTat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)

      CALL InflowWind_DestroyInputFile( InputFileData,  TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)

         ! Close the summary file if we were writing one
      IF ( SumFileUnit > 0 ) THEN
         CALL InflowWind_CloseSumFile( SumFileUnit, TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      ENDIF


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


      IMPLICIT                                                    NONE

      CHARACTER(*),              PARAMETER                     :: RoutineName="InflowWind_CalcOutput"


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

      TYPE(IfW_UniformWind_OutputType)                         :: Uniform_OutData     !< output velocities
      TYPE(IfW_TSFFWind_OutputType)                            :: TSFF_OutData        !< output velocities
      TYPE(IfW_BladedFFWind_OutputType)                        :: BladedFF_OutData    !< output velocities
      TYPE(IfW_UserWind_OutputType)                            :: User_OutData        !< output velocities

      REAL(ReKi), ALLOCATABLE                                  :: PositionXYZprime(:,:)   !< PositionXYZ array in the prime (wind) coordinates

      INTEGER(IntKi)                                           :: I                   !< Generic counters


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
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
         IF ( ErrStat >= AbortErrLev ) RETURN
      ELSEIF ( SIZE(InputData%PositionXYZ,DIM=2) /= SIZE(OutputData%VelocityUVW,DIM=2) ) THEN
         CALL SetErrStat( ErrID_Fatal," Programming error: Position and Velocity arrays are not sized the same.",  &
               ErrStat, ErrMsg, RoutineName)
      ENDIF

      !-----------------------------
      ! Outputs: OutputData%VelocityUVW and OutputData%DiskVel
      !-----------------------------
      

      CALL CalculateOutput( Time, InputData, ParamData, &
                       ContStates, DiscStates, ConstrStates, & 
                       OtherStates, OutputData, .TRUE., TmpErrStat, TmpErrMsg )      
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )

      !-----------------------------
      ! Outputs: OutputData%lidar%LidSpeed and OutputData%lidar%WtTrunc
      !-----------------------------
      
         ! return sensor values
      IF (ParamData%lidar%SensorType /= SensorType_None) THEN
         
         CALL Lidar_CalcOutput(Time, InputData, ParamData, &
                              ContStates, DiscStates, ConstrStates, OtherStates, &  
                              OutputData, TmpErrStat, TmpErrMsg )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
         
      END IF      
       
      
      !-----------------------------
      ! Outputs: OutputData%WriteOutput from OtherState%AllOuts and OutputData%lidar%LidSpeed
      !-----------------------------

      CALL SetAllOuts( ParamData, OutputData, OtherStates, TmpErrStat, TmpErrMsg ) 
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
   
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

      IMPLICIT                                                    NONE

      CHARACTER(*),              PARAMETER                     :: RoutineName="InflowWind_End"

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

      TYPE(IfW_UniformWind_OutputType)                         :: Uniform_OutData        !< output velocities

      TYPE(IfW_TSFFWind_OutputType)                            :: TSFF_OutData        !< output velocities

      TYPE(IfW_BladedFFWind_OutputType)                        :: BladedFF_OutData        !< output velocities

      TYPE(IfW_UserWind_OutputType)                            :: User_OutData        !< output velocities


     ErrStat = ErrID_None
     ErrMsg = ""

         ! End the sub-modules (deallocates their arrays and closes their files):

      SELECT CASE ( ParamData%WindType )

         CASE (Steady_WindNumber)         ! The Steady wind is a simple wrapper for the UniformWind module.
            CALL IfW_UniformWind_End( InputData%PositionXYZ, ParamData%UniformWind, &
                                 OtherStates%UniformWind, Uniform_OutData, ErrStat, ErrMsg )

         CASE (Uniform_WindNumber)
            CALL IfW_UniformWind_End( InputData%PositionXYZ, ParamData%UniformWind, &
                                 OtherStates%UniformWind, Uniform_OutData, ErrStat, ErrMsg )

         CASE (TSFF_WindNumber)
            CALL IfW_TSFFWind_End( InputData%PositionXYZ, ParamData%TSFFWind, &
                                 OtherStates%TSFFWind, TSFF_OutData, ErrStat, ErrMsg )

         CASE (BladedFF_WindNumber)
            CALL IfW_BladedFFWind_End( InputData%PositionXYZ, ParamData%BladedFFWind, &
                                 OtherStates%BladedFFWind, BladedFF_OutData, ErrStat, ErrMsg )

!!!         CASE (HAWC_WindNumber)
!!!            CALL HW_Terminate(     ErrStat )

         CASE (User_WindNumber)
            CALL IfW_UserWind_End( InputData%PositionXYZ, ParamData%UserWind, &
                                 OtherStates%UserWind, User_OutData, ErrStat, ErrMsg )

         CASE ( Undef_WindNumber )
            ! Do nothing

         CASE DEFAULT  ! keep this check to make sure that all new wind types have been accounted for
            CALL SetErrStat(ErrID_Fatal,' Undefined wind type.',ErrStat,ErrMsg,RoutineName)

      END SELECT

!!!  !   IF (CTTS_Flag) CALL CTTS_Terminate( ErrStat ) !FIXME: should it be this line or the next?
!!!         CALL CTTS_Terminate( ErrStat, ErrMsg )


      CALL InflowWind_DestroyInput( InputData, ErrStat, ErrMsg )         
      CALL InflowWind_DestroyParam( ParamData, ErrStat, ErrMsg )         
      CALL InflowWind_DestroyContState( ContStates, ErrStat, ErrMsg )         
      CALL InflowWind_DestroyDiscState( DiscStates, ErrStat, ErrMsg )         
      CALL InflowWind_DestroyConstrState( ConstrStateGuess, ErrStat, ErrMsg )         
      CALL InflowWind_DestroyOtherState( OtherStates, ErrStat, ErrMsg )         
      CALL InflowWind_DestroyOutput( OutData, ErrStat, ErrMsg )                     
      
      
         ! Reset the wind type so that the initialization routine must be called
      ParamData%WindType      = Undef_WindNumber
      ParamData%Initialized   = .FALSE.
      ParamData%CTTS_Flag     = .FALSE.


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


         ! Initialize ErrStat


      ! BJJ: Please don't make my code end just because I called a routine that you don't use :)
      ErrStat = ErrID_None
      ErrMsg  = ""

      x%DummyContState     = 0.0_ReKi
      xd%DummyDiscState    = 0.0_ReKi
      z%DummyConstrState   = 0.0_ReKi
      
      RETURN


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
END MODULE InflowWind
