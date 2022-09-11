!**********************************************************************************************************************************
!! This module is used to read and process the (undisturbed) inflow winds.  It must be initialized
!! using InflowWind_Init() with the name of the file, the file type, and possibly reference height and
!! width (depending on the type of wind file being used).  This module calls appropriate routines
!! in the wind modules so that the type of wind becomes seamless to the user.  InflowWind_End()
!! should be called when the program has finshed.
!!
!! Data are assumed to be in units of meters and seconds.  Z is measured from the ground (NOT the hub!).
!!
!!  7 Oct 2009    Initial Release with AeroDyn 13.00.00         B. Jonkman, NREL/NWTC
!! 14 Nov 2011    v1.00.01b-bjj                                 B. Jonkman
!!  1 Aug 2012    v1.01.00a-bjj                                 B. Jonkman
!! 10 Aug 2012    v1.01.00b-bjj                                 B. Jonkman
!!    Feb 2013    v2.00.00a-adp   conversion to Framework       A. Platt
!!    Sep 2015    v3.00.00a-adb   added separate input file     A. Platt
!
!..................................................................................................................................
! Files with this module:
!  InflowWind_Subs.f90
!  InflowWind.txt       -- InflowWind_Types will be auto-generated based on the descriptions found in this file.
!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2017  Envision Energy
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
MODULE InflowWind


   USE                              InflowWind_Types
   USE                              NWTC_Library
   USE                              InflowWind_Subs

   USE                              Lidar                      ! module for obtaining sensor data
   
   IMPLICIT NONE
   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: IfW_Ver = ProgDesc( 'InflowWind', '', '' )




      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: InflowWind_Init                                   !< Initialization routine
   PUBLIC :: InflowWind_CalcOutput                             !< Calculate the wind velocities
   PUBLIC :: InflowWind_End                                    !< Ending routine (includes clean up)

   PUBLIC :: InflowWind_Convert2HAWC               !< An extension of the FAST framework, this routine converts an InflowWind data structure to HAWC format wind files
   PUBLIC :: InflowWind_Convert2Bladed             !< An extension of the FAST framework, this routine converts an InflowWind data structure to Bladed format wind files (with shear already included)
   PUBLIC :: InflowWind_Convert2VTK                !< An extension of the FAST framework, this routine converts an InflowWind data structure to VTK format wind files


      ! These routines satisfy the framework, but do nothing at present.
   PUBLIC :: InflowWind_UpdateStates               !< Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
   PUBLIC :: InflowWind_CalcConstrStateResidual    !< Tight coupling routine for returning the constraint state residual
   PUBLIC :: InflowWind_CalcContStateDeriv         !< Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: InflowWind_UpdateDiscState            !< Tight coupling routine for updating discrete states


      ! These routines compute Jacobians; only dYdu is defined.
   PUBLIC :: InflowWind_JacobianPInput             !< Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                                   !!   (Xd), and constraint - state(Z) functions all with respect to the inputs(u)
   PUBLIC :: InflowWind_JacobianPContState         !< Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                                   !!   (Xd), and constraint - state(Z) functions all with respect to the continuous
                                                   !!   states(x)
   PUBLIC :: InflowWind_JacobianPDiscState         !< Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                                   !!   (Xd), and constraint - state(Z) functions all with respect to the discrete
                                                   !!   states(xd)
   PUBLIC :: InflowWind_JacobianPConstrState       !< Routine to compute the Jacobians of the output(Y), continuous - (X), discrete -
                                                   !!   (Xd), and constraint - state(Z) functions all with respect to the constraint
                                                   !!   states(z)
   PUBLIC :: InflowWind_GetOP                      !< Routine to pack the operating point values (for linearization) into arrays
   


CONTAINS
!====================================================================================================
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
!! Since this module acts as an interface to other modules, on some things are set before initiating
!! calls to the lower modules.
!----------------------------------------------------------------------------------------------------
SUBROUTINE InflowWind_Init( InitInp,   InputGuess,    p, ContStates, DiscStates,    ConstrStateGuess,    OtherStates,   &
                     y,    m,   TimeInterval,  InitOutData, ErrStat, ErrMsg )

      IMPLICIT                                                 NONE

      CHARACTER(*),              PARAMETER                  :: RoutineName="InflowWind_Init"

         ! Initialization data and guesses

      TYPE(InflowWind_InitInputType),        INTENT(IN   )  :: InitInp           !< Input data for initialization
      TYPE(InflowWind_InputType),            INTENT(  OUT)  :: InputGuess        !< An initial guess for the input; the input mesh must be defined
      TYPE(InflowWind_ParameterType),        INTENT(  OUT)  :: p                 !< Parameters
      TYPE(InflowWind_ContinuousStateType),  INTENT(  OUT)  :: ContStates        !< Initial continuous states
      TYPE(InflowWind_DiscreteStateType),    INTENT(  OUT)  :: DiscStates        !< Initial discrete states
      TYPE(InflowWind_ConstraintStateType),  INTENT(  OUT)  :: ConstrStateGuess  !< Initial guess of the constraint states
      TYPE(InflowWind_OtherStateType),       INTENT(  OUT)  :: OtherStates       !< Initial other/optimization states
      TYPE(InflowWind_OutputType),           INTENT(  OUT)  :: y                 !< Initial output (outputs are not calculated; only the output mesh is initialized)
      TYPE(InflowWind_MiscVarType),          INTENT(  OUT)  :: m                 !< Misc variables for optimization (not copied in glue code)
      REAL(DbKi),                            INTENT(IN   )  :: TimeInterval      !< Coupling time interval in seconds: InflowWind does not change this.
      TYPE(InflowWind_InitOutputType),       INTENT(  OUT)  :: InitOutData       !< Initial output data -- Names, units, and version info.


         ! Error Handling

      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat           !< Error status of the operation
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None


         ! Local variables

      TYPE(InflowWind_InputFile)                            :: InputFileData        !< Data from input file

      TYPE(IfW_UniformWind_InitInputType)                   :: Uniform_InitData     !< initialization info
      TYPE(IfW_UniformWind_InitOutputType)                  :: Uniform_InitOutData  !< initialization output info

      TYPE(IfW_TSFFWind_InitInputType)                      :: TSFF_InitData        !< initialization info
      TYPE(IfW_TSFFWind_InitOutputType)                     :: TSFF_InitOutData     !< initialization output info

      TYPE(IfW_HAWCWind_InitInputType)                      :: HAWC_InitData        !< initialization info
      TYPE(IfW_HAWCWind_InitOutputType)                     :: HAWC_InitOutData     !< initialization output info

      TYPE(IfW_BladedFFWind_InitInputType)                  :: BladedFF_InitData    !< initialization info
      TYPE(IfW_BladedFFWind_InitOutputType)                 :: BladedFF_InitOutData !< initialization output info

      TYPE(IfW_UserWind_InitInputType)                      :: User_InitData        !< initialization info
      TYPE(IfW_UserWind_InitOutputType)                     :: User_InitOutData     !< initialization info

      TYPE(IfW_4Dext_InitOutputType)                        :: FDext_InitOutData    !< initialization info

      TYPE(FileInfoType)                                    :: InFileInfo    !< The derived type for holding the full input file for parsing -- we may pass this in the future
!!!     TYPE(CTBladed_Backgr)                                        :: BackGrndValues


         ! Temporary variables for error handling
      INTEGER(IntKi)                                        :: TmpErrStat
      CHARACTER(ErrMsgLen)                                  :: TmpErrMsg         !< temporary error message
      CHARACTER(1024)                                       :: PriPath

         ! Local Variables
      INTEGER(IntKi)                                        :: I, j              !< Generic counter
      INTEGER(IntKi)                                        :: SumFileUnit       !< Unit number for the summary file
      CHARACTER(256)                                        :: SumFileName       !< Name of the summary file
      CHARACTER(256)                                        :: EchoFileName      !< Name of the summary file
      CHARACTER(1), PARAMETER                               :: UVW(3) = (/'U','V','W'/)
      CHARACTER(1), PARAMETER                               :: XYZ(3) = (/'X','Y','Z'/)

         !----------------------------------------------------------------------------------------------
         ! Initialize variables and check to see if this module has been initialized before.
         !----------------------------------------------------------------------------------------------

      ErrStat        =  ErrID_None
      ErrMsg         =  ""
      SumFileUnit    =  -1_IntKi ! set at beginning in case of error

         ! Set a few variables.

      p%DT   = TimeInterval             ! InflowWind does not require a specific time interval, so this is never changed.
      CALL NWTC_Init()
      CALL DispNVD( IfW_Ver )



         !----------------------------------------------------------------------------------------------
         ! Read the input file
         !----------------------------------------------------------------------------------------------


         ! Set the names of the files based on the inputfilename
      p%RootFileName  = InitInp%RootName
      IF (LEN_TRIM(p%RootFileName) == 0) CALL GetRoot( InitInp%InputFileName, p%RootFileName )
      EchoFileName  = TRIM(p%RootFileName)//".ech"
      SumFileName   = TRIM(p%RootFileName)//".sum"

         ! these values (and others hard-coded in lidar_init) should be set in the input file, too
      InputFileData%SensorType = InitInp%lidar%SensorType
      InputFileData%NumPulseGate = InitInp%lidar%NumPulseGate
      InputFileData%RotorApexOffsetPos = InitInp%lidar%RotorApexOffsetPos
      InputFileData%LidRadialVel = InitInp%lidar%LidRadialVel

         ! Parse all the InflowWind related input files and populate the *_InitDataType derived types
      CALL GetPath( InitInp%InputFileName, PriPath )

      IF ( InitInp%UseInputFile ) THEN
         CALL ProcessComFile( InitInp%InputFileName, InFileInfo, TmpErrStat, TmpErrMsg )
         ! For diagnostic purposes, the following can be used to display the contents
         ! of the InFileInfo data structure.
         ! call Print_FileInfo_Struct( CU, InFileInfo ) ! CU is the screen -- different number on different systems.

         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF
                        
      ELSE
         CALL NWTC_Library_CopyFileInfoType( InitInp%PassedFileData, InFileInfo, MESH_NEWCOPY, TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF          
         
      ENDIF

      CALL InflowWind_ParseInputFileInfo( InputFileData,  InFileInfo, PriPath, InitInp%InputFileName, EchoFileName, InitInp%FixedWindFileRootName, InitInp%TurbineID, TmpErrStat, TmpErrMsg )

      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL Cleanup()
         CALL InflowWind_DestroyInputFile( InputFileData, TmpErrStat, TmpErrMsg )
         RETURN
      ENDIF
         ! let's tell InflowWind if an external module (e.g., FAST.Farm) is going to set the velocity grids.
      
      IF ( InitInp%Use4Dext) then
         InputFileData%WindType = FDext_WindNumber      
         InputFileData%PropagationDir = 0.0_ReKi ! wind is in XYZ coordinates (already rotated if necessary), so don't rotate it again
      END IF

         ! initialize sensor data:   
      CALL Lidar_Init( InitInp, InputGuess, p, ContStates, DiscStates, ConstrStateGuess, OtherStates,   &
                       y, m, TimeInterval, InitOutData, TmpErrStat, TmpErrMsg )
         CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )      
      

         ! Validate the InflowWind input file information.

      CALL InflowWind_ValidateInput( InitInp, InputFileData, TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         
      IF ( ErrStat>= AbortErrLev ) THEN
         CALL Cleanup()
         RETURN
      ENDIF


      
         ! If a summary file was requested, open it.
      IF ( InputFileData%SumPrint ) THEN

            ! Open the summary file and write some preliminary info to it
         CALL InflowWind_OpenSumFile( SumFileUnit, SumFileName, IfW_Ver, InputFileData%WindType, TmpErrStat, TmpErrMsg )
            CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
            IF (ErrStat >= AbortErrLev) THEN
               CALL Cleanup()
               RETURN
            ENDIF
      ELSE
         SumFileUnit =  -1_IntKi       ! So that we don't try to write to something.  Used as indicator in submodules.
      ENDIF


      ! Allocate the arrays for passing points in and velocities out
      CALL AllocAry( InputGuess%PositionXYZ, 3, InitInp%NumWindPoints, &
                  "Array of positions at which to find wind velocities", TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         
      CALL AllocAry( y%VelocityUVW, 3, InitInp%NumWindPoints, &
                  "Array of wind velocities returned by InflowWind", TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat>= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF
      InputGuess%PositionXYZ = 0.0_ReKi
      y%VelocityUVW = 0.0_ReKi


      !-----------------------------------------------------------------
      ! Initialize the submodules based on the WindType
      !-----------------------------------------------------------------

      
      InitOutData%WindFileInfo%MWS = HUGE(InitOutData%WindFileInfo%MWS)

      SELECT CASE ( InputFileData%WindType )


         CASE ( Steady_WindNumber )

               !  This is a simplified case of the Uniform wind.  For this, we set the OtherStates data manually and don't
               !  call UniformWind_Init.  We do however call it for the calculations.


               ! Set InitInp information -- It isn't necessary to do this since that information is only used in the Uniform_Init routine, which is not called.

               ! Set the Otherstates information
            p%UniformWind%NumDataLines     =  1_IntKi
            p%UniformWind%RefHt            =  InputFileData%Steady_RefHt
            p%UniformWind%RefLength        =  1.0_ReKi    ! This is not used since no shear gusts are used.  Set to 1.0 so calculations don't bomb. 
            m%UniformWind%TimeIndex =  1           ! Used in UniformWind as a check if initialization was done.


            CALL AllocAry( p%UniformWind%Tdata, p%UniformWind%NumDataLines, 'Uniform wind time', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         
            CALL AllocAry( p%UniformWind%V, p%UniformWind%NumDataLines, 'Uniform wind horizontal wind speed', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         
            CALL AllocAry( p%UniformWind%Delta, p%UniformWind%NumDataLines, 'Uniform wind direction', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         
            CALL AllocAry( p%UniformWind%Upflow, p%UniformWind%NumDataLines, 'Uniform wind upflow angle', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         
            CALL AllocAry( p%UniformWind%VZ, p%UniformWind%NumDataLines, 'Uniform vertical wind speed', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         
            CALL AllocAry( p%UniformWind%HShr, p%UniformWind%NumDataLines, 'Uniform horizontal linear shear', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         
            CALL AllocAry( p%UniformWind%VShr, p%UniformWind%NumDataLines, 'Uniform vertical power-law shear exponent', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         
            CALL AllocAry( p%UniformWind%VLinShr, p%UniformWind%NumDataLines, 'Uniform vertical linear shear', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
               
            CALL AllocAry( p%UniformWind%VGust, p%UniformWind%NumDataLines, 'Uniform gust velocity', TmpErrStat, TmpErrMsg )
               CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               ENDIF

               ! Set the array information
            
            p%UniformWind%Tdata(  :)       = 0.0_ReKi
            p%UniformWind%V(      :)       = InputFileData%Steady_HWindSpeed
            p%UniformWind%Delta(  :)       = 0.0_ReKi
            p%UniformWind%Upflow( :)       = 0.0_ReKi
            p%UniformWind%VZ(     :)       = 0.0_ReKi
            p%UniformWind%HShr(   :)       = 0.0_ReKi
            p%UniformWind%VShr(   :)       = InputFileData%Steady_PLexp
            p%UniformWind%VLinShr(:)       = 0.0_ReKi
            p%UniformWind%VGust(  :)       = 0.0_ReKi



               ! Now we have in effect initialized the IfW_UniformWind module, so set the parameters
            p%UniformWind%RefLength        =  1.0_ReKi       ! This is done so that 0.0*(some stuff)/RefLength doesn't blow up.
            p%UniformWind%RefHt            =  InputFileData%Steady_RefHt
            m%UniformWind%TimeIndex        =  1_IntKi

            p%ReferenceHeight = p%UniformWind%RefHt
            
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
            InitOutData%WindFileInfo%MWS              = InputFileData%Steady_HWindSpeed

               ! Write summary file information
            IF ( SumFileUnit > 0 ) THEN
               WRITE(SumFileUnit,'(A)',IOSTAT=TmpErrStat)
               WRITE(SumFileUnit,'(A80)',IOSTAT=TmpErrStat)          'Steady wind -- Constant wind profile for entire simulation. No windfile read in.'
               WRITE(SumFileUnit,'(A40,G12.4)',IOSTAT=TmpErrStat)    '     Reference height:                  ',p%UniformWind%RefHt
               WRITE(SumFileUnit,'(A40,G12.4)',IOSTAT=TmpErrStat)    '     Horizontal velocity:               ',p%UniformWind%V
               WRITE(SumFileUnit,'(A40,G12.4)',IOSTAT=TmpErrStat)    '     Vertical sheer power law exponent: ',p%UniformWind%VShr

                  ! We are assuming that if the last line was written ok, then all of them were.
               IF (TmpErrStat /= 0_IntKi) THEN
                  CALL SetErrStat(ErrID_Fatal,'Error writing to summary file.',ErrStat,ErrMsg,RoutineName)
                  CALL Cleanup
                  RETURN
               ENDIF   
            ENDIF 


         CASE ( Uniform_WindNumber )

               ! Set InitInp information
            Uniform_InitData%ReferenceHeight          =  InputFileData%Uniform_RefHt
            Uniform_InitData%RefLength                =  InputFileData%Uniform_RefLength 
            Uniform_InitData%WindFileName             =  InputFileData%Uniform_FileName
            Uniform_InitData%SumFileUnit              =  SumFileUnit

            Uniform_InitData%UseInputFile             =  InitInp%WindType2UseInputFile
            Uniform_InitData%PassedFileData           =  InitInp%WindType2Data

               ! Initialize the UniformWind module
            CALL IfW_UniformWind_Init(Uniform_InitData, p%UniformWind, &
                        m%UniformWind, Uniform_InitOutData,  TmpErrStat, TmpErrMsg)

               CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, ' IfW_Init' )
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               ENDIF
               
            if (InitInp%Linearize) then
                  ! we'd have to redo the math to get this correct, so for now we are disabling upflow for linearization:
               if (any(p%UniformWind%upflow /= 0.0_ReKi) ) then
                  call SetErrStat(ErrID_Fatal, 'Upflow in uniform wind files must be 0 for linearization analysis in InflowWind.', ErrStat, ErrMsg, RoutineName)
                  CALL Cleanup()
                  return
               end if
            end if


            p%ReferenceHeight = p%UniformWind%RefHt

               ! Store wind file metadata
            InitOutData%WindFileInfo%FileName         =  InputFileData%Uniform_FileName
            InitOutData%WindFileInfo%WindType         =  Uniform_WindNumber
            InitOutData%WindFileInfo%RefHt            =  p%UniformWind%RefHt
            InitOutData%WindFileInfo%RefHt_Set        =  .FALSE.                             ! The wind file does not set this
            InitOutData%WindFileInfo%DT               =  Uniform_InitOutData%WindFileDT
            InitOutData%WindFileInfo%NumTSteps        =  Uniform_InitOutData%WindFileNumTSteps
            InitOutData%WindFileInfo%ConstantDT       =  Uniform_InitOutData%WindFileConstantDT
            InitOutData%WindFileInfo%TRange           =  Uniform_InitOutData%WindFileTRange
            InitOutData%WindFileInfo%TRange_Limited   = .FALSE.                              ! UniformWind sets to limit of file if outside time bounds
            InitOutData%WindFileInfo%YRange           =  (/ 0.0_ReKi, 0.0_ReKi /)
            InitOutData%WindFileInfo%YRange_Limited   =  .FALSE.                             ! Hard boundaries not enforced in y-direction
            InitOutData%WindFileInfo%ZRange           =  (/ 0.0_ReKi, 0.0_ReKi /)
            InitOutData%WindFileInfo%ZRange_Limited   =  .FALSE.                             ! Hard boundaries not enforced in z-direction
            InitOutData%WindFileInfo%BinaryFormat     =  0_IntKi
            InitOutData%WindFileInfo%IsBinary         =  .FALSE.
            InitOutData%WindFileInfo%TI               =  0.0_ReKi
            InitOutData%WindFileInfo%TI_listed        =  .FALSE.

            if (p%UniformWind%NumDataLines == 1) then
               InitOutData%WindFileInfo%MWS = p%UniformWind%V(1)
            else
               InitOutData%WindFileInfo%MWS = 0.0_ReKi               
               do i=2,p%UniformWind%NumDataLines
                  InitOutData%WindFileInfo%MWS = InitOutData%WindFileInfo%MWS + &
                                                 0.5_ReKi*(p%UniformWind%V(i)+p%UniformWind%V(i-1))*&
                                                          (p%UniformWind%Tdata(i)-p%UniformWind%Tdata(i-1))
               end do
               InitOutData%WindFileInfo%MWS = InitOutData%WindFileInfo%MWS / &
                           ( p%UniformWind%Tdata(p%UniformWind%NumDataLines) - p%UniformWind%Tdata(1) )
            end if
            
            
               ! Check if the fist data point from the file is not along the X-axis while applying the windfield rotation
            IF ( ( .NOT. EqualRealNos (p%UniformWind%Delta(1), 0.0_ReKi) ) .AND.  &
                 ( .NOT. EqualRealNos (p%PropagationDir, 0.0_ReKi)       ) ) THEN
               CALL SetErrStat( ErrID_Warn,' Possible double rotation of wind field! Uniform wind file starts with a wind direction of '// &
                        TRIM(Num2LStr(p%UniformWind%Delta(1)*R2D))//                       &
                        ' degrees and the InflowWind input file specifies a PropagationDir of '//  &
                        TRIM(Num2LStr(p%PropagationDir*R2D))//' degrees.',                 &
                        ErrStat,ErrMsg,RoutineName )
            ENDIF



         CASE ( TSFF_WindNumber )

               ! Set InitInp information
            TSFF_InitData%WindFileName                   =  InputFileData%TSFF_FileName
            TSFF_InitData%SumFileUnit                    =  SumFileUnit

               ! Initialize the TSFFWind module
            CALL IfW_TSFFWind_Init(TSFF_InitData, p%TSFFWind, &
                           m%TSFFWind, TSFF_InitOutData,  TmpErrStat, TmpErrMsg)
            CALL SetErrSTat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               ENDIF

               ! Store wind file metadata
            InitOutData%WindFileInfo%FileName            =  InputFileData%TSFF_FileName
            
            CALL SetFFInitOutData(p%TSFFWind%FF)
            p%ReferenceHeight = InitOutData%WindFileInfo%RefHt


         CASE ( BladedFF_WindNumber, BladedFF_Shr_WindNumber )

               ! Set InitInp information
            BladedFF_InitData%SumFileUnit                =  SumFileUnit
            BladedFF_InitData%FixedWindFileRootName      = InitInp%FixedWindFileRootName
            BladedFF_InitData%TurbineID                  = InitInp%TurbineID
            
            if (InputFileData%WindType /= BladedFF_Shr_WindNumber) then  
               IF ( InitInp%FixedWindFileRootName ) THEN ! .TRUE. when FAST.Farm uses multiple instances of InflowWind for ambient wind data
                  IF ( InitInp%TurbineID == 0 ) THEN     ! .TRUE. for the FAST.Farm low-resolution domain
                     InputFileData%BladedFF_FileName = TRIM(InputFileData%BladedFF_FileName)//TRIM(PathSep)//'Low'
                  ELSE                                   ! FAST.Farm high-resolution domain(s)
                     InputFileData%BladedFF_FileName = TRIM(InputFileData%BladedFF_FileName)//TRIM(PathSep)//'HighT'//TRIM(Num2Lstr(InitInp%TurbineID))
                  ENDIF
               ENDIF
               
               BladedFF_InitData%WindFileName            = TRIM(InputFileData%BladedFF_FileName)//'.wnd'
               BladedFF_InitData%TowerFileExist          =  InputFileData%BladedFF_TowerFile
               BladedFF_InitData%NativeBladedFmt         = .false.
            else
               BladedFF_InitData%WindFileName            =  InputFileData%BladedFF_FileName
               BladedFF_InitData%TowerFileExist          = .false.
               BladedFF_InitData%NativeBladedFmt         = .true.
               !call IfW_FFWind_CopyInitInput( InputFileData%FF, BladedFF_InitData%FF, MESH_NEWCOPY, TmpErrStat, TmpErrMsg)
            end if
                        
               ! Initialize the BladedFFWind module
            CALL IfW_BladedFFWind_Init(BladedFF_InitData, p%BladedFFWind, m%BladedFFWind, &
                                        BladedFF_InitOutData,  TmpErrStat, TmpErrMsg)
            CALL SetErrSTat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               ENDIF

               ! Store wind file metadata
            InitOutData%WindFileInfo%FileName            =  InputFileData%BladedFF_FileName
            
            CALL SetFFInitOutData(p%BladedFFWind%FF)

            InitOutData%WindFileInfo%TI               =  BladedFF_InitOutData%TI
            InitOutData%WindFileInfo%TI_listed        =  .TRUE.   ! This must be listed in the file someplace
            
            if (InputFileData%WindType == BladedFF_Shr_WindNumber) then
               InputFileData%WindType = BladedFF_WindNumber
               ! this overwrites the values of PropagationDir and VFlowAngle with values from the native Bladed file
               InputFileData%PropagationDir = BladedFF_InitOutData%PropagationDir
               InputFileData%VFlowAngle     = BladedFF_InitOutData%VFlowAngle 
            end if
            p%ReferenceHeight = InitOutData%WindFileInfo%RefHt

            
         CASE ( HAWC_WindNumber )
            
               ! Set InitInp information
            HAWC_InitData%WindFileName(1)    = InputFileData%HAWC_FileName_u
            HAWC_InitData%WindFileName(2)    = InputFileData%HAWC_FileName_v
            HAWC_InitData%WindFileName(3)    = InputFileData%HAWC_FileName_w
            HAWC_InitData%SumFileUnit        = SumFileUnit
            HAWC_InitData%nx                 = InputFileData%HAWC_nx
            HAWC_InitData%ny                 = InputFileData%HAWC_ny
            HAWC_InitData%nz                 = InputFileData%HAWC_nz

            HAWC_InitData%dx                 = InputFileData%HAWC_dx
            HAWC_InitData%dy                 = InputFileData%HAWC_dy
            HAWC_InitData%dz                 = InputFileData%HAWC_dz

            call IfW_FFWind_CopyInitInput( InputFileData%FF, HAWC_InitData%FF, MESH_NEWCOPY, TmpErrStat, TmpErrMsg)
                   
            
               ! Initialize the HAWCWind module
            CALL IfW_HAWCWind_Init(HAWC_InitData, p%HAWCWind, m%HAWCWind, &
                                   TimeInterval,  HAWC_InitOutData,  TmpErrStat, TmpErrMsg)
            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               ENDIF


               ! Store wind file metadata
            CALL SetFFInitOutData(p%HAWCWind%FF)
            InitOutData%WindFileInfo%FileName            =  InputFileData%HAWC_FileName_u
            p%ReferenceHeight = InitOutData%WindFileInfo%RefHt

            
         CASE (User_WindNumber)

               ! Initialize the UserWind module
            CALL IfW_UserWind_Init(User_InitData, p%UserWind, m%UserWind, &
                        TimeInterval,  User_InitOutData,  TmpErrStat, TmpErrMsg)
            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               ENDIF
            
            p%ReferenceHeight = InputFileData%Steady_RefHt ! FIXME!!!!
            
         CASE ( FDext_WindNumber )
            
               ! Initialize the UserWind module
            CALL IfW_4Dext_Init(InitInp%FDext, p%FDext, m%FDext, TimeInterval, FDext_InitOutData, TmpErrStat, TmpErrMsg)
            CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName)
               IF (ErrStat >= AbortErrLev) THEN
                  CALL Cleanup()
                  RETURN
               ENDIF
            p%ReferenceHeight = p%FDext%pZero(3) + (p%FDext%n(3)/2) * p%FDext%delta(3) ! should be middle of grid, right???? FIXME
            
         CASE DEFAULT  ! keep this check to make sure that all new wind types have been accounted for
            CALL SetErrStat(ErrID_Fatal,' Undefined wind type.',ErrStat,ErrMsg,'InflowWind_Init()')




         END SELECT
                  

      !IF ( InputFileData%CTTS_Flag ) THEN
      !   ! Initialize the CTTS_Wind module
      !ENDIF


      !............................................
      ! Set the p and OtherStates for InflowWind using the input file information.
      ! (set this after initializing modules so that we can use propagationDir and VFlowAng from native-Bladed files
      !............................................

      CALL InflowWind_SetParameters( InitInp, InputFileData, p, m, TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat>= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF

            

!!!         !----------------------------------------------------------------------------------------------
!!!         ! Check for coherent turbulence file (KH superimposed on a background wind file)
!!!         ! Initialize the CTWind module and initialize the module of the other wind type.
!!!         !----------------------------------------------------------------------------------------------
!!!
!!!      IF ( p%WindType == CTP_WindNumber ) THEN
!!!
!!!!FIXME: remove this error message when we add CTP_Wind in
!!!            CALL SetErrStat( ErrID_Fatal, ' InflowWind cannot currently handle the CTP_Wind type.', ErrStat, ErrMsg, ' IfW_Init' )
!!!            RETURN
!!!
!!!         CALL CTTS_Init(UnWind, p%WindFileName, BackGrndValues, ErrStat, ErrMsg)
!!!         IF (ErrStat /= 0) THEN
!!!            p%WindType = Undef_Wind
!!!            ErrStat  = ErrID_Fatal
!!!            RETURN
!!!         END IF
!!!
!!!   !FIXME: check this
!!!         p%WindFileName = BackGrndValues%WindFile
!!!         p%WindType = BackGrndValues%WindType
!!!   !      CTTS_Flag  = BackGrndValues%CoherentStr
!!!         p%CTTS_Flag  = BackGrndValues%CoherentStr    ! This might be wrong
!!!
!!!      ELSE
!!!
!!!         p%CTTS_Flag  = .FALSE.
!!!
!!!      END IF
!!!



      
         ! Allocate arrays for the WriteOutput

      CALL AllocAry( y%WriteOutput, p%NumOuts, 'WriteOutput', TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat>= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF
      y%WriteOutput = 0.0_ReKi
      
      CALL AllocAry( InitOutData%WriteOutputHdr, p%NumOuts, 'WriteOutputHdr', TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      CALL AllocAry( InitOutData%WriteOutputUnt, p%NumOuts, 'WriteOutputUnt', TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         IF ( ErrStat>= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         ENDIF
   
      InitOutData%WriteOutputHdr = p%OutParam(1:p%NumOuts)%Name
      InitOutData%WriteOutputUnt = p%OutParam(1:p%NumOuts)%Units     
       

      ! allocate and fill variables for linearization:
      if (InitInp%Linearize) then
         
         CALL AllocAry(InitOutData%LinNames_u, InitInp%NumWindPoints*3 + 3, 'LinNames_u', TmpErrStat, TmpErrMsg)
            CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         CALL AllocAry(InitOutData%RotFrame_u, InitInp%NumWindPoints*3 + 3, 'RotFrame_u', TmpErrStat, TmpErrMsg)
            CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         CALL AllocAry(InitOutData%IsLoad_u, InitInp%NumWindPoints*3 + 3, 'IsLoad_u', TmpErrStat, TmpErrMsg)
            CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         CALL AllocAry(InitOutData%LinNames_y, InitInp%NumWindPoints*3+p%NumOuts, 'LinNames_y', TmpErrStat, TmpErrMsg)
            CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         CALL AllocAry(InitOutData%RotFrame_y, InitInp%NumWindPoints*3+p%NumOuts, 'RotFrame_y', TmpErrStat, TmpErrMsg)
            CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
         IF (ErrStat >= AbortErrLev) THEN
            CALL Cleanup()
            RETURN
         ENDIF
         
         do i=1,InitInp%NumWindPoints
            do j=1,3
               InitOutData%LinNames_y((i-1)*3+j) = UVW(j)//'-component inflow velocity at node '//trim(num2lstr(i))//', m/s'
               InitOutData%LinNames_u((i-1)*3+j) = XYZ(j)//'-component position of node '//trim(num2lstr(i))//', m'
            end do            
         end do

         InitOutData%LinNames_u(InitInp%NumWindPoints*3 + 1) = 'Extended input: horizontal wind speed (steady/uniform wind), m/s'
         InitOutData%LinNames_u(InitInp%NumWindPoints*3 + 2) = 'Extended input: vertical power-law shear exponent, -'
         InitOutData%LinNames_u(InitInp%NumWindPoints*3 + 3) = 'Extended input: propagation direction, rad'         
         
         do i=1,p%NumOuts
            InitOutData%LinNames_y(i+3*InitInp%NumWindPoints) = trim(p%OutParam(i)%Name)//', '//p%OutParam(i)%Units
         end do

         ! IfW inputs and outputs are in the global, not rotating frame
         InitOutData%RotFrame_u = .false. 
         InitOutData%RotFrame_y = .false. 

         InitOutData%IsLoad_u = .false. ! IfW inputs for linearization are not loads
         
         !InitOutData%PropagationDir = -p%PropagationDir
         !InitOutData%RefHt = p%UniformWind%RefHt
         !InitOutData%RefLength = p%UniformWind%RefLength
         
      end if
                  

         ! Set the version information in InitOutData
      InitOutData%Ver   = IfW_Ver


      CALL CleanUp()


      RETURN


   !----------------------------------------------------------------------------------------------------
CONTAINS

   SUBROUTINE CleanUp()

      ! add in stuff that we need to dispose of here
      CALL InflowWind_DestroyInputFile( InputFileData,  TmpErrStat, TmpErrMsg )
      CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)

      ! Ignore error messages from InFileInfo destruction
      call NWTC_Library_DestroyFileInfoType( InFileInfo, TmpErrStat, TmpErrMsg )

         ! Close the summary file if we were writing one
      IF ( SumFileUnit > 0 ) THEN
         CALL InflowWind_CloseSumFile( SumFileUnit, TmpErrStat, TmpErrMsg )
         CALL SetErrStat(TmpErrStat,TmpErrMsg,ErrStat,ErrMsg,RoutineName)
      ENDIF


   END SUBROUTINE CleanUp


   SUBROUTINE SetFFInitOutData(FFp)
   
      TYPE(IfW_FFWind_ParameterType),              INTENT(IN   )  :: FFp                !< Parameters

   
            InitOutData%WindFileInfo%WindType            =  p%WindType
            InitOutData%WindFileInfo%RefHt               =  FFp%RefHt
            InitOutData%WindFileInfo%RefHt_Set           =  .TRUE.
            InitOutData%WindFileInfo%DT                  =  FFp%FFDTime
            InitOutData%WindFileInfo%NumTSteps           =  FFp%NFFSteps
            InitOutData%WindFileInfo%ConstantDT          =  .TRUE.
            IF ( FFp%Periodic ) THEN
               InitOutData%WindFileInfo%TRange           =  (/ 0.0_ReKi, FFp%TotalTime /)
               InitOutData%WindFileInfo%TRange_Limited   =  .FALSE.
            ELSE  ! Shift the time range to compensate for the shifting of the wind grid
               InitOutData%WindFileInfo%TRange           =  (/ 0.0_ReKi, FFp%TotalTime /)  -  FFp%InitXPosition*FFp%InvMFFWS
               InitOutData%WindFileInfo%TRange_Limited   =  .TRUE.
            ENDIF
            InitOutData%WindFileInfo%YRange              =  (/ -FFp%FFYHWid, FFp%FFYHWid /)
            InitOutData%WindFileInfo%YRange_Limited      =  .TRUE.      ! Hard boundaries enforced in y-direction
            IF ( p%TSFFWind%FF%NTGrids > 0 ) THEN        ! have tower data
               InitOutData%WindFileInfo%ZRange           =  (/ 0.0_Reki,  FFp%RefHt + FFp%FFZHWid /)
            ELSE
               InitOutData%WindFileInfo%ZRange           =  (/ FFp%GridBase,    &
                                                               FFp%GridBase + FFp%FFZHWid*2.0 /)
            ENDIF
            InitOutData%WindFileInfo%ZRange_Limited      =  .TRUE.
            InitOutData%WindFileInfo%BinaryFormat        =  FFp%WindFileFormat
            InitOutData%WindFileInfo%IsBinary            =  .TRUE.
            InitOutData%WindFileInfo%MWS                 = FFp%MeanFFWS

            InitOutData%WindFileInfo%TI                  =  0.0_ReKi
            InitOutData%WindFileInfo%TI_listed           =  .FALSE.
                  
   END SUBROUTINE SetFFInitOutData
   
END SUBROUTINE InflowWind_Init


!====================================================================================================
!> This routine takes an input dataset of type InputType which contains a position array of dimensions 3*n. It then calculates
!! and returns the output dataset of type OutputType which contains a corresponding velocity array of dimensions 3*n. The input
!! array contains XYZ triplets for each position of interest (first index is X/Y/Z for values 1/2/3, second index is the point
!! number to evaluate). The returned values in the OutputData are similar with U/V/W for the first index of 1/2/3.
!!
!! _Coordinate transformation:_
!! The coordinates passed in are copied to the PositionXYZPrime array, then rotated by -(p%PropagationDir) (radians) so
!! that the wind direction lies along the X-axis (all wind files are given this way).  The submodules are then called with
!! these PositionXYZPrime coordinates.
!!
!! After the calculation by the submodule, the PositionXYZPrime coordinate array is deallocated.  The returned VelocityUVW
!! array is then rotated by p%PropagationDir so that it now corresponds the the global coordinate UVW values for wind
!! with that direction.
!----------------------------------------------------------------------------------------------------
SUBROUTINE InflowWind_CalcOutput( Time, InputData, p, &
                              ContStates, DiscStates, ConstrStates, &   ! Framework required states -- empty in this case.
                              OtherStates, OutputData, m, ErrStat, ErrMsg )


   IMPLICIT                                                    NONE

   CHARACTER(*),              PARAMETER                     :: RoutineName="InflowWind_CalcOutput"


      ! Inputs / Outputs

   REAL(DbKi),                               INTENT(IN   )  :: Time              !< Current simulation time in seconds
   TYPE(InflowWind_InputType),               INTENT(IN   )  :: InputData         !< Inputs at Time
   TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: p                 !< Parameters
   TYPE(InflowWind_ContinuousStateType),     INTENT(IN   )  :: ContStates        !< Continuous states at Time
   TYPE(InflowWind_DiscreteStateType),       INTENT(IN   )  :: DiscStates        !< Discrete states at Time
   TYPE(InflowWind_ConstraintStateType),     INTENT(IN   )  :: ConstrStates      !< Constraint states at Time
   TYPE(InflowWind_OtherStateType),          INTENT(IN   )  :: OtherStates       !< Other/optimization states at Time
   TYPE(InflowWind_OutputType),              INTENT(INOUT)  :: OutputData        !< Outputs computed at Time (IN for mesh reasons and data allocation)
   TYPE(InflowWind_MiscVarType),             INTENT(INOUT)  :: m                 !< Misc variables for optimization (not copied in glue code)

   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< Error status of the operation
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None


      ! Local variables
   INTEGER(IntKi)                                           :: i
      ! Temporary variables for error handling
   INTEGER(IntKi)                                           :: TmpErrStat
   CHARACTER(ErrMsgLen)                                     :: TmpErrMsg            ! temporary error message



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
      

   CALL CalculateOutput( Time, InputData, p, &
                     ContStates, DiscStates, ConstrStates, & 
                     OtherStates, OutputData, m, .TRUE., TmpErrStat, TmpErrMsg )      
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )

   !-----------------------------
   ! Outputs: OutputData%lidar%LidSpeed and OutputData%lidar%WtTrunc
   !-----------------------------
      
      ! return sensor values
   IF (p%lidar%SensorType /= SensorType_None) THEN
         
      CALL Lidar_CalcOutput(Time, InputData, p, &
                           ContStates, DiscStates, ConstrStates, OtherStates, &  
                           OutputData, m, TmpErrStat, TmpErrMsg )
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
         
   END IF      
       
      
   !-----------------------------
   ! Outputs: OutputData%WriteOutput from m%AllOuts and OutputData%lidar%LidSpeed
   !-----------------------------

   CALL SetAllOuts( p, OutputData, m, TmpErrStat, TmpErrMsg ) 
      CALL SetErrStat( TmpErrStat, TmpErrMsg, ErrStat, ErrMsg, RoutineName )
   
      ! Map to the outputs
   DO I = 1,p%NumOuts  ! Loop through all selected output channels
      OutputData%WriteOutput(I) = p%OutParam(I)%SignM * m%AllOuts( p%OutParam(I)%Indx )
   ENDDO             ! I - All selected output channels


END SUBROUTINE InflowWind_CalcOutput



!====================================================================================================
!> Clean up the allocated variables and close all open files.  Reset the initialization flag so
!! that we have to reinitialize before calling the routines again.
SUBROUTINE InflowWind_End( InputData, p, ContStates, DiscStates, ConstrStateGuess, OtherStates, &
                       y, m, ErrStat, ErrMsg )

   IMPLICIT                                                    NONE

   CHARACTER(*),              PARAMETER                     :: RoutineName="InflowWind_End"

      ! Initialization data and guesses

   TYPE(InflowWind_InputType),               INTENT(INOUT)  :: InputData         !< Input data for initialization
   TYPE(InflowWind_ParameterType),           INTENT(INOUT)  :: p         !< Parameters
   TYPE(InflowWind_ContinuousStateType),     INTENT(INOUT)  :: ContStates        !< Continuous states
   TYPE(InflowWind_DiscreteStateType),       INTENT(INOUT)  :: DiscStates        !< Discrete states
   TYPE(InflowWind_ConstraintStateType),     INTENT(INOUT)  :: ConstrStateGuess  !< Guess of the constraint states
   TYPE(InflowWind_OtherStateType),          INTENT(INOUT)  :: OtherStates       !< Other/optimization states
   TYPE(InflowWind_OutputType),              INTENT(INOUT)  :: y           !< Output data
   TYPE(InflowWind_MiscVarType),             INTENT(INOUT)  :: m          !< Misc variables for optimization (not copied in glue code)


      ! Error Handling

   INTEGER( IntKi ),                         INTENT(  OUT)  :: ErrStat           !< error status
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< error message


   ErrStat = ErrID_None
   ErrMsg = ""

      ! End the sub-modules (deallocates their arrays and closes their files):

   SELECT CASE ( p%WindType )

      CASE (Steady_WindNumber, Uniform_WindNumber)         ! The Steady wind is a simple wrapper for the UniformWind module.
         CALL IfW_UniformWind_End( p%UniformWind, m%UniformWind, ErrStat, ErrMsg )

      CASE (TSFF_WindNumber)
         CALL IfW_TSFFWind_End( p%TSFFWind, m%TSFFWind, ErrStat, ErrMsg )

      CASE (BladedFF_WindNumber)
         CALL IfW_BladedFFWind_End( p%BladedFFWind, m%BladedFFWind, ErrStat, ErrMsg )

      CASE (HAWC_WindNumber)
         CALL IfW_HAWCWind_End( p%HAWCWind, m%HAWCWind, ErrStat, ErrMsg )

      CASE (User_WindNumber)
         CALL IfW_UserWind_End( p%UserWind, m%UserWind, ErrStat, ErrMsg )

      CASE (FDext_WindNumber)
         CALL IfW_4Dext_End( p%FDext, m%FDext, ErrStat, ErrMsg )
         
      CASE ( Undef_WindNumber )
         ! Do nothing

      CASE DEFAULT  ! keep this check to make sure that all new wind types have been accounted for
         CALL SetErrStat(ErrID_Fatal,' Undefined wind type.',ErrStat,ErrMsg,RoutineName)

   END SELECT

!!!  !   IF (CTTS_Flag) CALL CTTS_Terminate( ErrStat ) !FIXME: should it be this line or the next?
!!!         CALL CTTS_Terminate( ErrStat, ErrMsg )


   CALL InflowWind_DestroyInput( InputData, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyParam( p, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyContState( ContStates, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyDiscState( DiscStates, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyConstrState( ConstrStateGuess, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyOtherState( OtherStates, ErrStat, ErrMsg )         
   CALL InflowWind_DestroyOutput( y, ErrStat, ErrMsg )                     
   CALL InflowWind_DestroyMisc( m, ErrStat, ErrMsg )                     
      
      
      ! Reset the wind type so that the initialization routine must be called
   p%WindType      = Undef_WindNumber
   p%CTTS_Flag     = .FALSE.


END SUBROUTINE InflowWind_End


!====================================================================================================
! The following routines were added to satisfy the framework, but do nothing useful.
!====================================================================================================
!> This is a loose coupling routine for solving constraint states, integrating continuous states, and updating discrete and other 
!! states. Continuous, constraint, discrete, and other states are updated to values at t + Interval.
SUBROUTINE InflowWind_UpdateStates( t, n, Inputs, InputTimes, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

   REAL(DbKi),                            INTENT(IN   ) :: t               !< Current simulation time in seconds
   INTEGER(IntKi),                        INTENT(IN   ) :: n               !< Current step of the simulation: t = n*Interval
   TYPE(InflowWind_InputType),            INTENT(INOUT) :: Inputs(:)       !< Inputs at InputTimes (output only for mesh record-keeping in ExtrapInterp routine)
   REAL(DbKi),                            INTENT(IN   ) :: InputTimes(:)   !< Times in seconds associated with Inputs
   TYPE(InflowWind_ParameterType),        INTENT(IN   ) :: p               !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(INOUT) :: x               !< Input: Continuous states at t;
                                                                           !!    Output: Continuous states at t + Interval
   TYPE(InflowWind_DiscreteStateType),    INTENT(INOUT) :: xd              !< Input: Discrete states at t;
                                                                           !!    Output: Discrete states at t  + Interval
   TYPE(InflowWind_ConstraintStateType),  INTENT(INOUT) :: z               !< Input: Constraint states at t;
                                                                           !!   Output: Constraint states at t + Interval
   TYPE(InflowWind_OtherStateType),       INTENT(INOUT) :: OtherState      !< Other states: Other states at t;
                                                                           !!   Output: Other states at t + Interval
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT) :: m               !< Misc variables for optimization (not copied in glue code)
   INTEGER(IntKi),                        INTENT(  OUT) :: ErrStat         !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT) :: ErrMsg          !< Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""

   x%DummyContState     = 0.0_ReKi
   xd%DummyDiscState    = 0.0_ReKi
   z%DummyConstrState   = 0.0_ReKi
      
   RETURN


END SUBROUTINE InflowWind_UpdateStates

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for computing derivatives of continuous states
SUBROUTINE InflowWind_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, m, dxdt, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                               INTENT(IN   )  :: Time        !< Current simulation time in seconds
   TYPE(InflowWind_InputType),               INTENT(IN   )  :: u           !< Inputs at Time
   TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: p           !< Parameters
   TYPE(InflowWind_ContinuousStateType),     INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(InflowWind_DiscreteStateType),       INTENT(IN   )  :: xd          !< Discrete states at Time
   TYPE(InflowWind_ConstraintStateType),     INTENT(IN   )  :: z           !< Constraint states at Time
   TYPE(InflowWind_OtherStateType),          INTENT(IN   )  :: OtherState  !< Other states at Time
   TYPE(InflowWind_MiscVarType),             INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
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
!> Tight coupling routine for updating discrete states
SUBROUTINE InflowWind_UpdateDiscState( Time, u, p, x, xd, z, OtherState, m, ErrStat, ErrMsg )

   REAL(DbKi),                               INTENT(IN   )  :: Time        !< Current simulation time in seconds
   TYPE(InflowWind_InputType),               INTENT(IN   )  :: u           !< Inputs at Time
   TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: p           !< Parameters
   TYPE(InflowWind_ContinuousStateType),     INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(InflowWind_DiscreteStateType),       INTENT(INOUT)  :: xd          !< Input: Discrete states at Time;
                                                                           !! Output: Discrete states at Time + Interval
   TYPE(InflowWind_ConstraintStateType),     INTENT(IN   )  :: z           !< Constraint states at Time
   TYPE(InflowWind_OtherStateType),          INTENT(IN   )  :: OtherState  !< Other states at Time
   TYPE(InflowWind_MiscVarType),             INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Update discrete states here:

   ! StateData%DiscState =

END SUBROUTINE InflowWind_UpdateDiscState

!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
SUBROUTINE InflowWind_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )

   REAL(DbKi),                               INTENT(IN   )  :: Time        !< Current simulation time in seconds
   TYPE(InflowWind_InputType),               INTENT(IN   )  :: u           !< Inputs at Time
   TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: p           !< Parameters
   TYPE(InflowWind_ContinuousStateType),     INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(InflowWind_DiscreteStateType),       INTENT(IN   )  :: xd          !< Discrete states at Time
   TYPE(InflowWind_ConstraintStateType),     INTENT(IN   )  :: z           !< Constraint states at Time (possibly a guess)
   TYPE(InflowWind_OtherStateType),          INTENT(IN   )  :: OtherState  !< Other states at Time
   TYPE(InflowWind_MiscVarType),             INTENT(INOUT)  :: m           !< Misc variables for optimization (not copied in glue code)
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

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ###### The following four routines are Jacobian routines for linearization capabilities #######
! If the module does not implement them, set ErrStat = ErrID_Fatal in IfW_Init() when InitInp%Linearize is .true.
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and dZ/du are returned.
SUBROUTINE InflowWind_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu )
!..................................................................................................................................

   REAL(DbKi),                            INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(InflowWind_InputType),            INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(InflowWind_ParameterType),        INTENT(IN   )           :: p          !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(InflowWind_DiscreteStateType),    INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(InflowWind_ConstraintStateType),  INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(InflowWind_OtherStateType),       INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(InflowWind_OutputType),           INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required); 
                                                                                !!   Output fields are not used by this routine, but type is   
                                                                                !!   available here so that mesh parameter information (i.e.,  
                                                                                !!   connectivity) does not have to be recalculated for dYdu.
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                        INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dYdu(:,:)  !< Partial derivatives of output functions (Y) 
                                                                                !!   with respect to inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXdu(:,:)  !< Partial derivatives of continuous state functions (X) 
                                                                                !!   with respect to inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXddu(:,:) !< Partial derivatives of discrete state functions (Xd) 
                                                                                !!   with respect to inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dZdu(:,:)  !< Partial derivatives of constraint state functions (Z)  
                                                                                !!   with respect to inputs (u) [intent in to avoid deallocation]
 
      ! local variables: 
   INTEGER(IntKi)                                                 :: ErrStat2
   CHARACTER(ErrMsgLen)                                           :: ErrMsg2            ! temporary error message
   CHARACTER(*), PARAMETER                                        :: RoutineName = 'InflowWind_JacobianPInput'
      
   REAL(R8Ki)                                                     :: local_dYdu(3,6)
   integer                                                        :: i, n
   integer                                                        :: i_start, i_end  ! indices for input/output start and end
   integer                                                        :: node, comp
   
      
      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''


   IF ( PRESENT( dYdu ) ) THEN

      ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:
         
         ! outputs are all velocities at all positions plus the WriteOutput values
         !
      if (.not. ALLOCATED(dYdu)) then
         CALL AllocAry( dYdu, SIZE(u%PositionXYZ)+p%NumOuts, SIZE(u%PositionXYZ)+3, 'dYdu', ErrStat2, ErrMsg2 )
         call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      end if
         
         
      SELECT CASE ( p%WindType )
      CASE (Steady_WindNumber, Uniform_WindNumber)

            ! note that we are including the propagation direction in the analytical derivative calculated
            ! inside IfW_UniformWind_JacobianPInput, so no need to transform input position vectors first
         
         dYdu = 0.0_R8Ki ! initialize all non-diagonal entries to zero (position of node effects the output of only that node) 
         
         n = SIZE(u%PositionXYZ,2)
            ! these are the positions used in the module coupling
         do i=1,n
            ! note that p%RotToWind(1,1) = cos(p%PropagationDir) and p%RotToWind(2,1) = sin(p%PropagationDir), which are the
            ! values we need to compute the jacobian.
!!!FIX ME with the propagation values!!!!         
            call IfW_UniformWind_JacobianPInput( t, u%PositionXYZ(:,i), p%RotToWind(1,1), p%RotToWind(2,1), p%UniformWind, m%UniformWind, local_dYdu )
            
            i_end  = 3*i
            i_start= i_end - 2
            
            dYdu(i_start:i_end,i_start:i_end) = local_dYdu(:,1:3)
            
            dYdu(i_start:i_end,n*3+1:) = local_dYdu(:,4:6) ! extended inputs
            
         end do            

            ! these are the InflowWind WriteOutput velocities (and note that we may not have all of the components of each point) 
         ! they do not depend on the inputs, so the derivatives w.r.t. X, Y, Z are all zero
         do i=1, p%NumOuts               
            node  = p%OutParamLinIndx(1,i) ! output node
            comp  = p%OutParamLinIndx(2,i) ! component of output node

            if (node > 0) then
!!!FIX ME with the propagation values!!!!         
               call IfW_UniformWind_JacobianPInput( t, p%WindViXYZ(:,node), p%RotToWind(1,1), p%RotToWind(2,1), p%UniformWind, m%UniformWind, local_dYdu )                                                                                       
            else
               local_dYdu = 0.0_R8Ki
            end if            
            
            dYdu(3*n+i, 3*n+1:) = p%OutParam(i)%SignM * local_dYdu( comp , 4:6)
         end do            

      CASE DEFAULT

      END SELECT         
                  
   END IF

   IF ( PRESENT( dXdu ) ) THEN
      if (allocated(dXdu)) deallocate(dXdu) 
   END IF

   IF ( PRESENT( dXddu ) ) THEN
      if (allocated(dXddu)) deallocate(dXddu) 
   END IF

   IF ( PRESENT( dZdu ) ) THEN
      if (allocated(dZdu)) deallocate(dZdu) 
   END IF


END SUBROUTINE InflowWind_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and dZ/dx are returned.
SUBROUTINE InflowWind_JacobianPContState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx, dXdx, dXddx, dZdx )
!..................................................................................................................................

   REAL(DbKi),                            INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(InflowWind_InputType),            INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(InflowWind_ParameterType),        INTENT(IN   )           :: p          !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(InflowWind_DiscreteStateType),    INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(InflowWind_ConstraintStateType),  INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(InflowWind_OtherStateType),       INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(InflowWind_OutputType),           INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required); 
                                                                                !!   Output fields are not used by this routine, but type is   
                                                                                !!   available here so that mesh parameter information (i.e.,  
                                                                                !!   connectivity) does not have to be recalculated for dYdx.
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                        INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dYdx(:,:)  !< Partial derivatives of output functions
                                                                                !!   (Y) with respect to the continuous
                                                                                !!   states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXdx(:,:)  !< Partial derivatives of continuous state
                                                                                !!   functions (X) with respect to
                                                                                !!   the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXddx(:,:) !< Partial derivatives of discrete state
                                                                                !!   functions (Xd) with respect to
                                                                                !!   the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dZdx(:,:)  !< Partial derivatives of constraint state
                                                                                !!   functions (Z) with respect to
                                                                                !!   the continuous states (x) [intent in to avoid deallocation]


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''



   IF ( PRESENT( dYdx ) ) THEN

      ! Calculate the partial derivative of the output functions (Y) with respect to the continuous states (x) here:

      ! allocate and set dYdx

   END IF

   IF ( PRESENT( dXdx ) ) THEN

      ! Calculate the partial derivative of the continuous state functions (X) with respect to the continuous states (x) here:

      ! allocate and set dXdx

   END IF

   IF ( PRESENT( dXddx ) ) THEN

      ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the continuous states (x) here:

      ! allocate and set dXddx

   END IF

   IF ( PRESENT( dZdx ) ) THEN


      ! Calculate the partial derivative of the constraint state functions (Z) with respect to the continuous states (x) here:

      ! allocate and set dZdx

   END IF


END SUBROUTINE InflowWind_JacobianPContState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and dZ/dxd are returned.
SUBROUTINE InflowWind_JacobianPDiscState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdxd, dXdxd, dXddxd, dZdxd )
!..................................................................................................................................

   REAL(DbKi),                            INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(InflowWind_InputType),            INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(InflowWind_ParameterType),        INTENT(IN   )           :: p          !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(InflowWind_DiscreteStateType),    INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(InflowWind_ConstraintStateType),  INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(InflowWind_OtherStateType),       INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(InflowWind_OutputType),           INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required); 
                                                                                !!   Output fields are not used by this routine, but type is   
                                                                                !!   available here so that mesh parameter information (i.e.,  
                                                                                !!   connectivity) does not have to be recalculated for dYdxd.
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                        INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dYdxd(:,:) !< Partial derivatives of output functions
                                                                                !!  (Y) with respect to the discrete
                                                                                !!  states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXdxd(:,:) !< Partial derivatives of continuous state
                                                                                !!   functions (X) with respect to the
                                                                                !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXddxd(:,:) !< Partial derivatives of discrete state
                                                                                !!   functions (Xd) with respect to the
                                                                                !!   discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dZdxd(:,:) !< Partial derivatives of constraint state
                                                                                !!   functions (Z) with respect to the
                                                                                !!   discrete states (xd) [intent in to avoid deallocation]


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''


   IF ( PRESENT( dYdxd ) ) THEN

      ! Calculate the partial derivative of the output functions (Y) with respect to the discrete states (xd) here:

      ! allocate and set dYdxd

   END IF

   IF ( PRESENT( dXdxd ) ) THEN

      ! Calculate the partial derivative of the continuous state functions (X) with respect to the discrete states (xd) here:

      ! allocate and set dXdxd

   END IF

   IF ( PRESENT( dXddxd ) ) THEN

      ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the discrete states (xd) here:

      ! allocate and set dXddxd

   END IF

   IF ( PRESENT( dZdxd ) ) THEN

      ! Calculate the partial derivative of the constraint state functions (Z) with respect to the discrete states (xd) here:

      ! allocate and set dZdxd

   END IF


END SUBROUTINE InflowWind_JacobianPDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and dZ/dz are returned.
SUBROUTINE InflowWind_JacobianPConstrState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdz, dXdz, dXddz, dZdz )
!..................................................................................................................................

   REAL(DbKi),                            INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(InflowWind_InputType),            INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(InflowWind_ParameterType),        INTENT(IN   )           :: p          !< Parameters
   TYPE(InflowWind_ContinuousStateType),  INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(InflowWind_DiscreteStateType),    INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(InflowWind_ConstraintStateType),  INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(InflowWind_OtherStateType),       INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(InflowWind_OutputType),           INTENT(IN   )           :: y          !< Output (change to inout if a mesh copy is required); 
                                                                                !!   Output fields are not used by this routine, but type is   
                                                                                !!   available here so that mesh parameter information (i.e.,  
                                                                                !!   connectivity) does not have to be recalculated for dYdz.
   TYPE(InflowWind_MiscVarType),          INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                        INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                          INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dYdz(:,:)  !< Partial derivatives of output
                                                                                !!  functions (Y) with respect to the
                                                                                !!  constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXdz(:,:)  !< Partial derivatives of continuous
                                                                                !!  state functions (X) with respect to
                                                                                !!  the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dXddz(:,:) !< Partial derivatives of discrete state
                                                                                !!  functions (Xd) with respect to the
                                                                                !!  constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,     INTENT(INOUT)           :: dZdz(:,:)  !< Partial derivatives of constraint
                                                                                !! state functions (Z) with respect to
                                                                                !!  the constraint states (z) [intent in to avoid deallocation]


      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

   IF ( PRESENT( dYdz ) ) THEN

         ! Calculate the partial derivative of the output functions (Y) with respect to the constraint states (z) here:

      ! allocate and set dYdz

   END IF

   IF ( PRESENT( dXdz ) ) THEN

         ! Calculate the partial derivative of the continuous state functions (X) with respect to the constraint states (z) here:

      ! allocate and set dXdz

   END IF

   IF ( PRESENT( dXddz ) ) THEN

         ! Calculate the partial derivative of the discrete state functions (Xd) with respect to the constraint states (z) here:

      ! allocate and set dXddz

   END IF

   IF ( PRESENT( dZdz ) ) THEN

         ! Calculate the partial derivative of the constraint state functions (Z) with respect to the constraint states (z) here:

      ! allocate and set dZdz

   END IF


END SUBROUTINE InflowWind_JacobianPConstrState
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Routine to pack the data structures representing the operating points into arrays for linearization.
SUBROUTINE InflowWind_GetOP( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, u_op, y_op, x_op, dx_op, xd_op, z_op )

   REAL(DbKi),                           INTENT(IN   )           :: t          !< Time in seconds at operating point
   TYPE(InflowWind_InputType),           INTENT(IN   )           :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(InflowWind_ParameterType),       INTENT(IN   )           :: p          !< Parameters
   TYPE(InflowWind_ContinuousStateType), INTENT(IN   )           :: x          !< Continuous states at operating point
   TYPE(InflowWind_DiscreteStateType),   INTENT(IN   )           :: xd         !< Discrete states at operating point
   TYPE(InflowWind_ConstraintStateType), INTENT(IN   )           :: z          !< Constraint states at operating point
   TYPE(InflowWind_OtherStateType),      INTENT(IN   )           :: OtherState !< Other states at operating point
   TYPE(InflowWind_OutputType),          INTENT(IN   )           :: y          !< Output at operating point
   TYPE(InflowWind_MiscVarType),         INTENT(INOUT)           :: m          !< Misc/optimization variables
   INTEGER(IntKi),                       INTENT(  OUT)           :: ErrStat    !< Error status of the operation
   CHARACTER(*),                         INTENT(  OUT)           :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: u_op(:)    !< values of linearized inputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: y_op(:)    !< values of linearized outputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: x_op(:)    !< values of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: dx_op(:)   !< values of first time derivatives of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: xd_op(:)   !< values of linearized discrete states
   REAL(ReKi), ALLOCATABLE, OPTIONAL,    INTENT(INOUT)           :: z_op(:)    !< values of linearized constraint states

   
   INTEGER(IntKi)                                    :: index, i, j
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'InflowWind_GetOP'
   

      ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = ''

   IF ( PRESENT( u_op ) ) THEN
      if (.not. allocated(u_op)) then
         call AllocAry(u_op, size(u%PositionXYZ)+3, 'u_op', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >= AbortErrLev) return
      end if
      
         
      index = 0
      do i=1,size(u%PositionXYZ,2)
         do j=1,size(u%PositionXYZ,1)
            index = index + 1 !(i-1)*size(u%PositionXYZ,1)+j
            u_op(index) = u%PositionXYZ(j,i)
         end do            
      end do  
      
      call IfW_UniformWind_GetOP( t, p%UniformWind, m%UniformWind, u_op(index+1:index+2) )
      u_op(index + 3) = p%PropagationDir
      
   END IF

   IF ( PRESENT( y_op ) ) THEN
      if (.not. allocated(y_op)) then
         call AllocAry(y_op, size(u%PositionXYZ)+p%NumOuts, 'y_op', ErrStat2, ErrMsg2)
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            if (ErrStat >= AbortErrLev) return
      end if
      
      index = 0
      do i=1,size(u%PositionXYZ,2)
         do j=1,size(u%PositionXYZ,1)
            index = index + 1 !(i-1)*size(u%PositionXYZ,1)+j
            y_op(index) = y%VelocityUVW(j,i)
         end do            
      end do
         
      do i=1,p%NumOuts         
         y_op(i+index) = y%WriteOutput( i )
      end do      
         
   END IF

   IF ( PRESENT( x_op ) ) THEN

   END IF

   IF ( PRESENT( dx_op ) ) THEN

   END IF

   IF ( PRESENT( xd_op ) ) THEN

   END IF
   
   IF ( PRESENT( z_op ) ) THEN

   END IF

END SUBROUTINE InflowWind_GetOP
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!====================================================================================================
SUBROUTINE InflowWind_Convert2HAWC( FileRootName, p, m, ErrStat, ErrMsg )

   USE IfW_FFWind_Base
   IMPLICIT NONE

   CHARACTER(*),              PARAMETER                     :: RoutineName="InflowWind_Convert2HAWC"

      ! Subroutine arguments

   TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: p                 !< Parameters
   TYPE(InflowWind_MiscVarType),             INTENT(INOUT)  :: m                 !< Misc/optimization variables
   CHARACTER(*),                             INTENT(IN   )  :: FileRootName      !< RootName for output files

   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< Error status of the operation
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None


      ! Local variables
   TYPE(IfW_FFWind_ParameterType)                           :: p_ff              !< FF Parameters
   INTEGER(IntKi)                                           :: ErrStat2
   CHARACTER(ErrMsgLen)                                     :: ErrMsg2


      ! Compute the wind velocities by stepping through all the data points and calling the appropriate GetWindSpeed routine
   SELECT CASE ( p%WindType )
         
   CASE (Steady_WindNumber, Uniform_WindNumber)

      CALL Uniform_to_FF(p%UniformWind, m%UniformWind, p_ff, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
      IF (ErrStat < AbortErrLev) THEN
         CALL ConvertFFWind_to_HAWC2(FileRootName, p_ff, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      END IF
            
      CALL IfW_FFWind_DestroyParam(p_ff,ErrStat2,ErrMsg2)
   
   CASE (TSFF_WindNumber)

      CALL ConvertFFWind_to_HAWC2(FileRootName, p%TSFFWind%FF, ErrStat, ErrMsg)

   CASE (BladedFF_WindNumber)

      CALL ConvertFFWind_to_HAWC2(FileRootName, p%BladedFFWind%FF, ErrStat, ErrMsg)

   CASE ( HAWC_WindNumber )

      CALL ConvertFFWind_to_HAWC2(FileRootName, p%HAWCWind%FF, ErrStat, ErrMsg)

   CASE DEFAULT ! User_WindNumber

      ErrStat = ErrID_Warn
      ErrMsg  = 'Wind type '//TRIM(Num2LStr(p%WindType))//' cannot be converted to HAWC format.'

   END SELECT
   
END SUBROUTINE InflowWind_Convert2HAWC

!====================================================================================================
SUBROUTINE InflowWind_Convert2Bladed( FileRootName, p, m, ErrStat, ErrMsg )

   USE IfW_FFWind_Base
   IMPLICIT NONE

   CHARACTER(*),              PARAMETER                     :: RoutineName="InflowWind_Convert2Bladed"

      ! Subroutine arguments

   TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: p                 !< Parameters
   TYPE(InflowWind_MiscVarType),             INTENT(INOUT)  :: m                 !< Misc/optimization variables
   CHARACTER(*),                             INTENT(IN   )  :: FileRootName      !< RootName for output files

   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< Error status of the operation
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None
   
      ! Local variables
   TYPE(IfW_FFWind_ParameterType)                           :: p_ff              !< FF Parameters
   INTEGER(IntKi)                                           :: ErrStat2
   CHARACTER(ErrMsgLen)                                     :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ""

      ! Local variables

      ! Compute the wind velocities by stepping through all the data points and calling the appropriate GetWindSpeed routine
   SELECT CASE ( p%WindType )
         
   CASE (Steady_WindNumber, Uniform_WindNumber)

      CALL Uniform_to_FF(p%UniformWind, m%UniformWind, p_ff, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
      IF (ErrStat < AbortErrLev) THEN
         CALL ConvertFFWind_to_Bladed(FileRootName, p_ff, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      END IF
            
      CALL IfW_FFWind_DestroyParam(p_ff,ErrStat2,ErrMsg2)

   CASE (TSFF_WindNumber)

      CALL ConvertFFWind_to_Bladed(FileRootName, p%TSFFWind%FF, ErrStat, ErrMsg)

   CASE (BladedFF_WindNumber)

      CALL ConvertFFWind_to_Bladed(FileRootName, p%BladedFFWind%FF, ErrStat, ErrMsg)

   CASE ( HAWC_WindNumber )

      CALL ConvertFFWind_to_Bladed(FileRootName, p%HAWCWind%FF, ErrStat, ErrMsg)

   CASE DEFAULT ! User_WindNumber

      ErrStat = ErrID_Warn
      ErrMsg  = 'Wind type '//TRIM(Num2LStr(p%WindType))//' cannot be converted to Bladed format.'

   END SELECT
   
END SUBROUTINE InflowWind_Convert2Bladed

!====================================================================================================
SUBROUTINE InflowWind_Convert2VTK( FileRootName, p, m, ErrStat, ErrMsg )

   USE IfW_FFWind_Base
   IMPLICIT NONE

   CHARACTER(*),              PARAMETER                     :: RoutineName="InflowWind_Convert2VTK"

      ! Subroutine arguments

   TYPE(InflowWind_ParameterType),           INTENT(IN   )  :: p                 !< Parameters
   TYPE(InflowWind_MiscVarType),             INTENT(INOUT)  :: m                 !< Misc/optimization variables
   CHARACTER(*),                             INTENT(IN   )  :: FileRootName      !< RootName for output files

   INTEGER(IntKi),                           INTENT(  OUT)  :: ErrStat           !< Error status of the operation
   CHARACTER(*),                             INTENT(  OUT)  :: ErrMsg            !< Error message if ErrStat /= ErrID_None
   
      ! Local variables
   TYPE(IfW_FFWind_ParameterType)                           :: p_ff              !< FF Parameters
   INTEGER(IntKi)                                           :: ErrStat2
   CHARACTER(ErrMsgLen)                                     :: ErrMsg2

   ErrStat = ErrID_None
   ErrMsg = ""

      ! Local variables

      ! Compute the wind velocities by stepping through all the data points and calling the appropriate GetWindSpeed routine
   SELECT CASE ( p%WindType )
         
   CASE (Steady_WindNumber, Uniform_WindNumber)

      CALL Uniform_to_FF(p%UniformWind, m%UniformWind, p_ff, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
            
      IF (ErrStat < AbortErrLev) THEN
         CALL ConvertFFWind_toVTK(FileRootName, p_ff, ErrStat2, ErrMsg2)
         CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      END IF
            
      CALL IfW_FFWind_DestroyParam(p_ff,ErrStat2,ErrMsg2)

   CASE (TSFF_WindNumber)

      CALL ConvertFFWind_toVTK(FileRootName, p%TSFFWind%FF, ErrStat, ErrMsg)

   CASE (BladedFF_WindNumber)

      CALL ConvertFFWind_toVTK(FileRootName, p%BladedFFWind%FF, ErrStat, ErrMsg)

   CASE ( HAWC_WindNumber )

      CALL ConvertFFWind_toVTK(FileRootName, p%HAWCWind%FF, ErrStat, ErrMsg)

   CASE DEFAULT ! User_WindNumber

      ErrStat = ErrID_Warn
      ErrMsg  = 'Wind type '//TRIM(Num2LStr(p%WindType))//' cannot be converted to VTK format.'

   END SELECT
   
END SUBROUTINE InflowWind_Convert2VTK
!====================================================================================================
END MODULE InflowWind
