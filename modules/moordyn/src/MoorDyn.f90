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
MODULE MoorDyn

   USE MoorDyn_Types
   USE MoorDyn_IO
   USE NWTC_Library
   USE MoorDyn_Line
   USE MoorDyn_Point
   USE MoorDyn_Rod
   USE MoorDyn_Body
   USE MoorDyn_Misc
   
   !USE WAVES, only: WaveGrid_n, WaveGrid_x0, WaveGrid_dx, WaveGrid_nx, WaveGrid_y0, WaveGrid_dy, WaveGrid_ny, WaveGrid_nz  ! seeing if I can get waves data here directly...

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: MD_ProgDesc = ProgDesc( 'MoorDyn', 'v2.0.0', '2022-12-08' )

   INTEGER(IntKi), PARAMETER            :: wordy = 0   ! verbosity level. >1 = more console output

   PUBLIC :: MD_Init
   PUBLIC :: MD_UpdateStates
   PUBLIC :: MD_CalcOutput
   PUBLIC :: MD_CalcContStateDeriv
   PUBLIC :: MD_End
   PUBLIC :: MD_JacobianPContState 
   PUBLIC :: MD_JacobianPInput 
   PUBLIC :: MD_JacobianPDiscState 
   PUBLIC :: MD_JacobianPConstrState 
   PUBLIC :: MD_GetOP 

CONTAINS

   !=========================================   MD_Init   ===================================
   SUBROUTINE MD_Init(InitInp, u, p, x, xd, z, other, y, m, DTcoupling, InitOut, ErrStat, ErrMsg)

      IMPLICIT NONE

      TYPE(MD_InitInputType),       INTENT(IN   )  :: InitInp     ! INTENT(INOUT) : Input data for initialization routine
      TYPE(MD_InputType),           INTENT(  OUT)  :: u           ! INTENT( OUT) : An initial guess for the input; input mesh must be defined
      TYPE(MD_ParameterType),       INTENT(  OUT)  :: p           ! INTENT( OUT) : Parameters
      TYPE(MD_ContinuousStateType), INTENT(  OUT)  :: x           ! INTENT( OUT) : Initial continuous states
      TYPE(MD_DiscreteStateType),   INTENT(  OUT)  :: xd          ! INTENT( OUT) : Initial discrete states
      TYPE(MD_ConstraintStateType), INTENT(  OUT)  :: z           ! INTENT( OUT) : Initial guess of the constraint states
      TYPE(MD_OtherStateType),      INTENT(  OUT)  :: other       ! INTENT( OUT) : Initial other states
      TYPE(MD_OutputType),          INTENT(  OUT)  :: y           ! INTENT( OUT) : Initial system outputs (outputs are not calculated; only the output mesh is initialized)
      TYPE(MD_MiscVarType),         INTENT(  OUT)  :: m           ! INTENT( OUT) : Initial misc/optimization variables
      REAL(DbKi),                   INTENT(INOUT)  :: DTcoupling  ! Coupling interval in seconds: the rate that Output is the actual coupling interval
      TYPE(MD_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
      TYPE(MD_InputFileType)                       :: InputFileDat   ! Data read from input file for setup, but not stored after Init
      type(FileInfoType)                           :: FileInfo_In    !< The derived type for holding the full input file for parsing -- we may pass this in the future
  !    CHARACTER(1024)                              :: priPath        ! The path to the primary MoorDyn input file
      REAL(DbKi)                                   :: t              ! instantaneous time, to be used during IC generation
      INTEGER(IntKi)                               :: l              ! index
      INTEGER(IntKi)                               :: I              ! Current line number of input file 
      INTEGER(IntKi)                               :: J              ! index
      INTEGER(IntKi)                               :: K              ! index
      INTEGER(IntKi)                               :: Itemp          ! index
      INTEGER(IntKi)                               :: iTurb          ! index for turbine in FAST.Farm applications
      INTEGER(IntKi)                               :: Converged      ! flag indicating whether the dynamic relaxation has converged
      INTEGER(IntKi)                               :: N              ! convenience integer for readability: number of segments in the line
      REAL(ReKi)                                   :: rPos(3)        ! array for setting fairlead reference positions in mesh
      REAL(ReKi)                                   :: OrMat(3,3)     ! rotation matrix for setting fairlead positions correctly if there is initial platform rotation
      REAL(ReKi)                                   :: OrMat2(3,3)
      REAL(R8Ki)                                   :: OrMatRef(3,3)
      REAL(DbKi), ALLOCATABLE                      :: FairTensIC(:,:)! array of size nCpldCons, 3 to store three latest fairlead tensions of each line
      CHARACTER(20)                                :: TempString     ! temporary string for incidental use
      INTEGER(IntKi)                               :: ErrStat2       ! Error status of the operation
      CHARACTER(ErrMsgLen)                         :: ErrMsg2        ! Error message if ErrStat2 /= ErrID_None
      
      REAL(DbKi)                                      :: dtM         ! actual mooring dynamics time step
      INTEGER(IntKi)                                  :: NdtM        ! number of time steps to integrate through with RK2
      INTEGER(IntKi)                                  :: ntWave      ! number of time steps of wave data
      
      TYPE(MD_InputType)    :: u_array(1)    ! a size-one array for u to make call to TimeStep happy
      REAL(DbKi)            :: t_array(1)    ! a size-one array saying time is 0 to make call to TimeStep happy  
      TYPE(MD_InputType)                                  :: u_interp   ! interpolated instantaneous input values to be calculated for each mooring time step



      CHARACTER(MaxWrScrLen)                       :: Message
      
      ! Local variables for reading file input (Previously in MDIO_ReadInput)
      INTEGER(IntKi)               :: UnEc                 ! The local unit number for this module's echo file
      INTEGER(IntKi)               :: UnOut    ! for outputing wave kinematics data
      CHARACTER(200)               :: Frmt     ! a string to hold a format statement

      CHARACTER(1024)              :: EchoFile             ! Name of MoorDyn echo file
      CHARACTER(1024)              :: Line                 ! String to temporarially hold value of read line
      CHARACTER(20)                :: LineOutString        ! String to temporarially hold characters specifying line output options
      CHARACTER(20)                :: OptString            ! String to temporarially hold name of option variable
      CHARACTER(40)                :: OptValue             ! String to temporarially hold value of options variable input
      CHARACTER(40)                :: DepthValue           ! Temporarily stores the optional WtrDpth setting for MD, which could be a number or a filename
      CHARACTER(40)                :: WaterKinValue        ! Temporarily stores the optional WaterKin setting for MD, which is typically a filename
      INTEGER(IntKi)               :: nOpts                ! number of options lines in input file
      CHARACTER(40)                :: TempString1          !
      CHARACTER(40)                :: TempString2          !
      CHARACTER(40)                :: TempString3          !
      CHARACTER(40)                :: TempString4          !
      CHARACTER(40)                :: TempString5          !
      CHARACTER(40)                :: TempStrings(6)       ! Array of 6 strings used when parsing comma-separated items
      CHARACTER(1024)              :: FileName             !

      REAL(DbKi)                   :: depth                ! local water depth interpolated from bathymetry grid [m]
      Real(DbKi)                   :: nvec(3)              ! local seabed surface normal vector (positive out)
      
      
      CHARACTER(25)                 :: let1                ! strings used for splitting and parsing identifiers
      CHARACTER(25)                 :: num1
      CHARACTER(25)                 :: let2
      CHARACTER(25)                 :: num2
      CHARACTER(25)                 :: let3
      
      REAL(DbKi)                    :: tempArray(6)
      REAL(ReKi)                    :: rRef(6)             ! used to pass positions to mesh (real type precision)
      REAL(DbKi)                    :: rRefDub(3)
      
      INTEGER(IntKi)               :: TempIDnums(100)      ! array to hold IdNums of controlled lines for each CtrlChan
      
      ! for reading output channels
      CHARACTER(ChanLen),ALLOCATABLE :: OutList(:)          ! array of output channel request (moved here from InitInput)
      INTEGER                       :: MaxAryLen = 1000    ! Maximum length of the array being read
      INTEGER                       :: NumWords            ! Number of words contained on a line
      INTEGER                       :: Nx
      INTEGER                       :: QuoteCh                                     ! Character position.
      CHARACTER(*), PARAMETER       :: RoutineName = 'MD_Init'

      

      ErrStat = ErrID_None
      ErrMsg  = ""
      m%zeros6 = 0.0_DbKi

      ! Initialize the NWTC Subroutine Library
      CALL NWTC_Init( )

      ! Display the module information
      CALL DispNVD( MD_ProgDesc )
      InitOut%Ver = MD_ProgDesc

      CALL WrScr('   This is MoorDyn v2, with significant input file changes from v1.')
      CALL WrScr('   Copyright: (C) 2022 National Renewable Energy Laboratory, (C) 2019 Matt Hall')


      !---------------------------------------------------------------------------------------------
      !                   Get all the inputs taken care of
      !---------------------------------------------------------------------------------------------

      p%RootName = TRIM(InitInp%RootName)//'.MD'  ! all files written from this module will have this root name

      ! set default values for the simulation settings
      ! these defaults are based on the glue code
      p%dtM0                 = DTcoupling      ! default to the coupling interval (but will likely need to be smaller)
      p%Tmax                 = InitInp%Tmax
      p%g                    = InitInp%g
      p%rhoW                 = InitInp%rhoW
      ! TODO:   add MSL2SWL from OpenFAST <<<<
      ! set the following to some defaults
      p%kBot                 = 3.0E6
      p%cBot                 = 3.0E5
      InputFileDat%dtIC      = 2.0_DbKi
      InputFileDat%TMaxIC    = 60.0_DbKi
      InputFileDat%CdScaleIC = 4.0_ReKi
      InputFileDat%threshIC  = 0.01_ReKi
      p%WaveKin              = 0_IntKi
      p%Current              = 0_IntKi
      p%dtOut                = 0.0_DbKi
      p%mu_kT                = 0.0_DbKi
      p%mu_kA                = 0.0_DbKi
      p%mc                   = 1.0_DbKi
      p%cv                   = 200.0_DbKi
      DepthValue = ""  ! Start off as empty string, to only be filled if MD setting is specified (otherwise InitInp%WtrDepth is used)
                       ! DepthValue and InitInp%WtrDepth are processed later by setupBathymetry.
      WaterKinValue = ""
      
      m%PtfmInit = InitInp%PtfmInit(:,1)   ! is this copying necssary in case this is an individual instance in FAST.Farm?



      ! Check if this MoorDyn instance is being run from FAST.Farm (indicated by FarmSize > 0)
      if (InitInp%FarmSize > 0) then
         CALL WrScr('   >>> MoorDyn is running in array mode <<< ')
         ! could make sure the size of this is right: SIZE(InitInp%FarmCoupledKinematics)  
         p%nTurbines = InitInp%FarmSize
      else    ! FarmSize==0 indicates normal, FAST module mode
         p%nTurbines = 1  ! if a regular FAST module mode, we treat it like a nTurbine=1 farm case
      END IF

      ! allocate some parameter arrays that are for each turbine (size 1 if regular OpenFAST use)
      allocate( p%nCpldBodies(     p%nTurbines)) 
      allocate( p%nCpldRods  (     p%nTurbines)) 
      allocate( p%nCpldCons  (     p%nTurbines)) 
      allocate( p%TurbineRefPos(3, p%nTurbines))
      
      ! initialize the arrays (to zero, except for passed in farm turbine reference positions)
      p%nCpldBodies = 0
      p%nCpldRods   = 0
      p%nCpldCons   = 0
      
      if (InitInp%FarmSize > 0) then
         p%TurbineRefPos = InitInp%TurbineRefPos  ! copy over turbine reference positions for later use
      else    
         p%TurbineRefPos = 0.0_DbKi               ! for now assuming this is zero for FAST use
      end if

      
      !---------------------------------------------------------------------------------------------
      !            read input file and create cross-referenced mooring system objects
      !---------------------------------------------------------------------------------------------
      
      
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""


      CALL WrScr( '   Parsing MoorDyn input file: '//trim(InitInp%FileName) )


      ! -----------------------------------------------------------------
      ! Read the primary MoorDyn input file, or copy from passed input
   if (InitInp%UsePrimaryInputFile) then
      ! Read the entire input file, minus any comment lines, into the FileInfo_In
      ! data structure in memory for further processing.
      call ProcessComFile( InitInp%FileName, FileInfo_In, ErrStat2, ErrMsg2 )
      CALL GetPath( InitInp%FileName, p%PriPath )    ! Input files will be relative to the path where the primary input file is located.
   else
      call NWTC_Library_CopyFileInfoType( InitInp%PassedPrimaryInputData, FileInfo_In, MESH_NEWCOPY, ErrStat2, ErrMsg2 )
      p%PriPath = ""
   endif
   if (Failed()) return;

   ! For diagnostic purposes, the following can be used to display the contents
   ! of the FileInfo_In data structure.
   !call Print_FileInfo_Struct( CU, FileInfo_In ) ! CU is the screen -- different number on different systems.

      !  Parse the FileInfo_In structure of data from the inputfile into the InitInp%InputFile structure
!   CALL ParsePrimaryFileInfo_BuildModel( PriPath, InitInp, FileInfo_In, InputFileDat, p, m, UnEc, ErrStat2, ErrMsg2 )
!   if (Failed()) return;




!NOTE: This could be split into a separate routine for easier to read code
      !-------------------------------------------------------------------------------------------------
      ! Parsing of input file from the FileInfo_In data structure
      !     -  FileInfo_Type is essentially a string array with some metadata.
      !-------------------------------------------------------------------------------------------------

      UnEc = -1
      nOpts = 0         ! Setting here rather than implied save


      ! ----------------- go through file contents a first time, counting each entry -----------------------

      i  = 0  ! set line number counter to before first line
      Line = NextLine(i);     ! Get the line and increment counter.  See description of routine. 
      
      do while ( i <= FileInfo_In%NumLines )

         if (INDEX(Line, "---") > 0) then ! look for a header line

            if ( ( INDEX(Line, "LINE DICTIONARY") > 0) .or. ( INDEX(Line, "LINE TYPES") > 0) ) then ! if line dictionary header

               ! skip following two lines (label line and unit line)
               Line = NextLine(i)
               Line = NextLine(i)
               
               ! find how many elements of this type there are
               Line = NextLine(i)
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  p%nLineTypes = p%nLineTypes + 1
                  Line = NextLine(i)
               END DO

            else if ( (INDEX(Line, "ROD DICTIONARY") > 0) .or. ( INDEX(Line, "ROD TYPES") > 0) ) then ! if rod dictionary header

               ! skip following two lines (label line and unit line)
               Line = NextLine(i)
               Line = NextLine(i)
               
               ! find how many elements of this type there are
               Line = NextLine(i)
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  p%nRodTypes = p%nRodTypes + 1
                  Line = NextLine(i)
               END DO

            else if ((INDEX(Line, "BODIES") > 0 ) .or. (INDEX(Line, "BODY LIST") > 0 ) .or. (INDEX(Line, "BODY PROPERTIES") > 0 )) then

               ! skip following two lines (label line and unit line)
               Line = NextLine(i)
               Line = NextLine(i)
               
               ! find how many elements of this type there are
               Line = NextLine(i)
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  p%nBodies = p%nBodies + 1
                  Line = NextLine(i)
               END DO

            else if ((INDEX(Line, "RODS") > 0 ) .or. (INDEX(Line, "ROD LIST") > 0) .or. (INDEX(Line, "ROD PROPERTIES") > 0)) then ! if rod properties header

               ! skip following two lines (label line and unit line)
               Line = NextLine(i)
               Line = NextLine(i)
               
               ! find how many elements of this type there are
               Line = NextLine(i)
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  p%nRods = p%nRods + 1
                  Line = NextLine(i)
               END DO

            else if ((INDEX(Line, "POINTS") > 0 ) .or. (INDEX(Line, "CONNECTION PROPERTIES") > 0) .or. (INDEX(Line, "NODE PROPERTIES") > 0) .or. (INDEX(Line, "POINT PROPERTIES") > 0) .or. (INDEX(Line, "POINT LIST") > 0) ) then ! if node properties header

               ! skip following two lines (label line and unit line)
               Line = NextLine(i)
               Line = NextLine(i)
               
               ! find how many elements of this type there are
               Line = NextLine(i)
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  p%nConnects = p%nConnects + 1
                  Line = NextLine(i)
               END DO

            else if ((INDEX(Line, "LINES") > 0 ) .or. (INDEX(Line, "LINE PROPERTIES") > 0) .or. (INDEX(Line, "LINE LIST") > 0) ) then ! if line properties header

               ! skip following two lines (label line and unit line)
               i=i+2
               
               ! find how many elements of this type there are
               Line = NextLine(i)
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  p%nLines = p%nLines + 1
                  Line = NextLine(i)
               END DO

            else if (INDEX(Line, "CONTROL") > 0) then ! if failure conditions header

               IF (wordy > 1) print *, "   Reading control channels: ";
               
               ! skip following two lines (label line and unit line)
               Line = NextLine(i)
               Line = NextLine(i)
               
               ! find how many elements of this type there are
               Line = NextLine(i)
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  p%nCtrlChans = p%nCtrlChans + 1
                  Line = NextLine(i)
               END DO
               
            else if (INDEX(Line, "FAILURE") > 0) then ! if failure conditions header

               IF (wordy > 1) print *, "   Reading failure conditions: ";
               
               ! skip following two lines (label line and unit line)
               Line = NextLine(i)
               Line = NextLine(i)
               
               ! find how many elements of this type there are
               Line = NextLine(i)
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  p%nFails = p%nFails + 1
                  Line = NextLine(i)
               END DO
               
               
            else if (INDEX(Line, "OPTIONS") > 0) then ! if options header

               IF (wordy > 0) print *, "Reading Options"
               
               ! don't skip any lines (no column headers for the options section)               
               ! process each line in this section
               Line = NextLine(i)
               DO while (INDEX(Line, "---") == 0) ! while we DON'T find another header line
                  
                  ! parse out entries:  value, option keyword
                  READ(Line,*,IOSTAT=ErrStat2) OptValue, OptString  ! look at first two entries, ignore remaining words in line, which should be comments
                  IF ( ErrStat2 /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Failed to read options.', ErrStat, ErrMsg, RoutineName ) ! would be nice to specify which line had the error
                     CALL CleanUp()
                     RETURN
                  END IF
                              
                  CALL Conv2UC(OptString)

                  ! check all possible options types and see if OptString is one of them, in which case set the variable.
                  if ( OptString == 'WRITELOG') THEN
                     read (OptValue,*) p%writeLog
                     if (p%writeLog > 0) then   ! if not zero, open a log file for output
                        CALL GetNewUnit( p%UnLog )
                        CALL OpenFOutFile ( p%UnLog, TRIM(p%RootName)//'.log', ErrStat, ErrMsg )
                        IF ( ErrStat > AbortErrLev ) THEN
                           ErrMsg = ' Failed to open MoorDyn log file: '//TRIM(ErrMsg)
                           RETURN
                        END IF
                        write(p%UnLog,'(A)', IOSTAT=ErrStat2) "MoorDyn v2 log file with output level "//TRIM(Num2LStr(p%writeLog))
                        write(p%UnLog,'(A)', IOSTAT=ErrStat2) "Note: options above the writeLog line in the input file will not be recorded."
                     end if
                  else if ( OptString == 'DTM') THEN
                     read (OptValue,*) p%dtM0 
                  else if ( OptString == 'G') then
                     read (OptValue,*) p%g
                  else if (( OptString == 'RHOW') .or. ( OptString == 'RHO')) then
                     read (OptValue,*) p%rhoW
                  else if (( OptString == 'WTRDPTH') .or. ( OptString == 'DEPTH') .or. ( OptString == 'WATERDEPTH')) then
                     read (OptValue,*) DepthValue    ! water depth input read in as a string to be processed by setupBathymetry
                  else if (( OptString == 'KBOT') .or. ( OptString == 'KB'))  then
                     read (OptValue,*) p%kBot
                  else if (( OptString == 'CBOT') .or. ( OptString == 'CB'))  then
                     read (OptValue,*) p%cBot
                  else if ( OptString == 'DTIC')  then
                     read (OptValue,*) InputFileDat%dtIC
                  else if ( OptString == 'TMAXIC')  then
                     read (OptValue,*) InputFileDat%TMaxIC
                  else if ( OptString == 'CDSCALEIC')  then
                     read (OptValue,*) InputFileDat%CdScaleIC
                  else if ( OptString == 'THRESHIC')  then
                     read (OptValue,*) InputFileDat%threshIC
                  else if ( OptString == 'WATERKIN')  then
                     read (OptValue,*) WaterKinValue
                  else if ( OptString == 'DTOUT')  then
                     read (OptValue,*) p%dtOut
                  else if ( OptString == 'MU_KT')  then
                     read (OptValue,*) p%mu_kT
                  else if ( OptString == 'MU_KA')  then
                     read (OptValue,*) p%mu_kA
                  else if ( OptString == 'MC')  then
                     read (OptValue,*) p%mc
                  else if ( OptString == 'CV')  then
                     read (OptValue,*) p%cv
                  else
                     CALL SetErrStat( ErrID_Warn, 'Unable to interpret input '//trim(OptString)//' in OPTIONS section.', ErrStat, ErrMsg, RoutineName )
                  end if

                  nOpts = nOpts + 1
                  Line = NextLine(i)
               END DO
               

            else if (INDEX(Line, "OUTPUT") > 0) then ! if output header

               ! we don't need to count this section...

               Line = NextLine(i)


            else  ! otherwise ignore this line that isn't a recognized header line and read the next line
               Line = NextLine(i)
            end if

         else ! otherwise ignore this line, which doesn't have the "---" or header line and read the next line
            Line = NextLine(i)
         end if
     
      end do

      p%nConnectsExtra = p%nConnects + 2*p%nLines    ! set maximum number of connections, accounting for possible detachment of each line end and a connection for that

      IF (wordy > 0) print *, "  Identified ", p%nLineTypes  , "LineTypes in input file."
      IF (wordy > 0) print *, "  Identified ", p%nRodTypes   , "RodTypes in input file."
      IF (wordy > 0) print *, "  Identified ", p%nBodies     , "Bodies in input file."
      IF (wordy > 0) print *, "  Identified ", p%nRods       , "Rods in input file."
      IF (wordy > 0) print *, "  Identified ", p%nConnects   , "Connections in input file."
      IF (wordy > 0) print *, "  Identified ", p%nLines      , "Lines in input file."
      IF (wordy > 0) print *, "  Identified ", nOpts         , "Options in input file."


      ! set up seabed bathymetry
      CALL setupBathymetry(DepthValue, InitInp%WtrDepth, m%BathymetryGrid, m%BathGrid_Xs, m%BathGrid_Ys, ErrStat2, ErrMsg2)
      CALL getDepthFromBathymetry(m%BathymetryGrid, m%BathGrid_Xs, m%BathGrid_Ys, 0.0_DbKi, 0.0_DbKi, p%WtrDpth, nvec)  ! set depth at 0,0 as nominal for waves etc
      
      
      ! set up wave and current kinematics 
      CALL setupWaterKin(WaterKinValue, p, InitInp%Tmax, ErrStat2, ErrMsg2); if(Failed()) return



      ! ----------------------------- misc checks to be sorted -----------------------------


    ! make sure nLineTypes isn't zero
    IF ( p%nLineTypes < 1 ) THEN
       CALL SetErrStat( ErrID_Fatal, 'nLineTypes parameter must be greater than zero.', ErrStat, ErrMsg, RoutineName )
       CALL CleanUp()
       RETURN
    END IF
    
    ! make sure NLines is at least one
    IF ( p%NLines < 1 ) THEN
       CALL SetErrStat( ErrID_Fatal, 'NLines parameter must be at least 1.', ErrStat, ErrMsg, RoutineName )
       CALL CleanUp()
       RETURN
    END IF


    
    

      ! ----------------------------- allocate necessary arrays ----------------------------


      ! Allocate object arrays

      ALLOCATE(m%LineTypeList(p%nLineTypes), STAT = ErrStat2 ); if(AllocateFailed("LineTypeList")) return
      ALLOCATE(m%RodTypeList( p%nRodTypes ), STAT = ErrStat2 ); if(AllocateFailed("LineTypeList")) return
               
      ALLOCATE(m%BodyList(    p%nBodies   ), STAT = ErrStat2 ); if(AllocateFailed("BodyList"    )) return
      ALLOCATE(m%RodList(     p%nRods     ), STAT = ErrStat2 ); if(AllocateFailed("RodList"     )) return
      ALLOCATE(m%ConnectList( p%nConnects ), STAT = ErrStat2 ); if(AllocateFailed("ConnectList" )) return
      ALLOCATE(m%LineList(    p%nLines    ), STAT = ErrStat2 ); if(AllocateFailed("LineList"    )) return
      
      ALLOCATE(m%FailList(    p%nFails    ), STAT = ErrStat2 ); if(AllocateFailed("FailList"    )) return
    
      
      ! Allocate associated index arrays (note: some are allocated larger than will be used, for simplicity)
      ALLOCATE(m%BodyStateIs1(p%nBodies ), m%BodyStateIsN(p%nBodies ), STAT=ErrStat2); if(AllocateFailed("BodyStateIs1/N")) return
      ALLOCATE(m%RodStateIs1(p%nRods    ), m%RodStateIsN(p%nRods    ), STAT=ErrStat2); if(AllocateFailed("RodStateIs1/N" )) return
      ALLOCATE(m%ConStateIs1(p%nConnects), m%ConStateIsN(p%nConnects), STAT=ErrStat2); if(AllocateFailed("ConStateIs1/N" )) return
      ALLOCATE(m%LineStateIs1(p%nLines)  , m%LineStateIsN(p%nLines)  , STAT=ErrStat2); if(AllocateFailed("LineStateIs1/N")) return

      ALLOCATE(m%FreeBodyIs(   p%nBodies ),  STAT=ErrStat2); if(AllocateFailed("FreeBodyIs")) return
      ALLOCATE(m%FreeRodIs(    p%nRods   ),  STAT=ErrStat2); if(AllocateFailed("FreeRodIs")) return
      ALLOCATE(m%FreeConIs(   p%nConnects),  STAT=ErrStat2); if(AllocateFailed("FreeConnectIs")) return

      ALLOCATE(m%CpldBodyIs(p%nBodies , p%nTurbines), STAT=ErrStat2); if(AllocateFailed("CpldBodyIs")) return
      ALLOCATE(m%CpldRodIs( p%nRods   , p%nTurbines), STAT=ErrStat2); if(AllocateFailed("CpldRodIs")) return
      ALLOCATE(m%CpldConIs(p%nConnects, p%nTurbines), STAT=ErrStat2); if(AllocateFailed("CpldConnectIs")) return


      ! ---------------------- now go through again and process file contents --------------------

      call Body_Setup( m%GroundBody, m%zeros6, p, ErrStat2, ErrMsg2)
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

      ! note: no longer worrying about "Echo" option
      
      Nx = 0  ! set state counter to zero
      i  = 0  ! set line number counter to before first line 
      Line = NextLine(i)
      
      do while ( i <= FileInfo_In%NumLines )
      
         if (INDEX(Line, "---") > 0) then ! look for a header line
             
            CALL Conv2UC(Line)  ! allow lowercase section header names as well

            !-------------------------------------------------------------------------------------------
            if ( ( INDEX(Line, "LINE DICTIONARY") > 0) .or. ( INDEX(Line, "LINE TYPES") > 0) ) then ! if line dictionary header
               
               IF (wordy > 0) print *, "Reading line types"
               
               ! skip following two lines (label line and unit line)
               Line = NextLine(i)
               Line = NextLine(i)
               
                ! process each line
                DO l = 1,p%nLineTypes
                   
                   !read into a line
                   Line = NextLine(i)
                   
                   ! check for correct number of columns in current line
                   IF ( CountWords( Line ) /= 10 ) THEN
                       CALL SetErrStat( ErrID_Fatal, ' Unable to parse Line type '//trim(Num2LStr(l))//' on row '//trim(Num2LStr(i))//' in input file. Row has wrong number of columns. Must be 10 columns.', ErrStat, ErrMsg, RoutineName )
                       CALL CleanUp()
                       RETURN
                   END IF

                   ! parse out entries: Name  Diam MassDenInAir EA cIntDamp EI    Cd  Ca  CdAx  CaAx 
                   READ(Line,*,IOSTAT=ErrStat2) m%LineTypeList(l)%name, m%LineTypeList(l)%d,  &
                      m%LineTypeList(l)%w, tempString1, tempString2, tempString3, &
                      m%LineTypeList(l)%Cdn, m%LineTypeList(l)%Can, m%LineTypeList(l)%Cdt, m%LineTypeList(l)%Cat
                   
                    IF ( ErrStat2 /= ErrID_None ) THEN
                      CALL SetErrStat( ErrID_Fatal, 'Failed to process line type inputs of entry '//trim(Num2LStr(l))//'. Check formatting and correct number of columns.', ErrStat, ErrMsg, RoutineName )
                      CALL CleanUp()
                      RETURN
                   END IF
                   
                   !TODO: add check if %name is maximum length, which might indicate the full name was too long <<<
                   
                   ! process stiffness coefficients
                   CALL SplitByBars(tempString1, N, tempStrings)
                   if (N > 2) then
                      CALL SetErrStat( ErrID_Fatal, 'A line type EA entry can have at most 2 (comma-separated) values.', ErrStat, ErrMsg, RoutineName )
                      CALL CleanUp()
                   else if (N==2) then                               ! visco-elastic case!
                      m%LineTypeList(l)%ElasticMod = 2
                      read(tempStrings(2), *) m%LineTypeList(l)%EA_D 
                   else
                     m%LineTypeList(l)%ElasticMod = 1                ! normal case
                   end if
                   ! get the regular/static coefficient or relation in all cases (can be from a lookup table)
                   CALL getCoefficientOrCurve(tempStrings(1), m%LineTypeList(l)%EA,     &
                                                           m%LineTypeList(l)%nEApoints, &
                                                           m%LineTypeList(l)%stiffXs,   &
                                                           m%LineTypeList(l)%stiffYs,  ErrStat2, ErrMsg2)

                   
                   ! process damping coefficients 
                   CALL SplitByBars(tempString2, N, tempStrings)
                   if (N > m%LineTypeList(l)%ElasticMod) then
                      CALL SetErrStat( ErrID_Fatal, 'A line type BA entry cannot have more (comma-separated) values its EA entry.', ErrStat, ErrMsg, RoutineName )
                      CALL CleanUp()
                   else if (N==2) then                               ! visco-elastic case when two BA values provided
                      read(tempStrings(2), *) m%LineTypeList(l)%BA_D 
                   else if (m%LineTypeList(l)%ElasticMod == 2) then  ! case where there is no dynamic damping for viscoelastic model (will it work)?
                      CALL WrScr("Warning, viscoelastic model being used with zero damping on the dynamic stiffness.")
                   end if
                   ! get the regular/static coefficient or relation in all cases (can be from a lookup table?)
                   CALL getCoefficientOrCurve(tempStrings(1), m%LineTypeList(l)%BA,     &
                                                           m%LineTypeList(l)%nBApoints,  &
                                                           m%LineTypeList(l)%dampXs,    &
                                                           m%LineTypeList(l)%dampYs,   ErrStat2, ErrMsg2)
                                                           
                   ! process bending stiffness coefficients (which might use lookup tables)
                   CALL getCoefficientOrCurve(tempString3, m%LineTypeList(l)%EI,        &
                                                           m%LineTypeList(l)%nEIpoints, &
                                                           m%LineTypeList(l)%bstiffXs,  &
                                                           m%LineTypeList(l)%bstiffYs, ErrStat2, ErrMsg2)

                   ! specify IdNum of line type for error checking
                   m%LineTypeList(l)%IdNum = l  
                   
                  ! write lineType information to log file
                  if (p%writeLog > 1) then
                      write(p%UnLog, '(A12,A20)'  ) " LineType"//trim(num2lstr(l))//":"
                      write(p%UnLog, '(A12,A20)'  ) " name: ", m%LineTypeList(l)%name
                      write(p%UnLog, '(A12,f12.4)') " d   : ", m%LineTypeList(l)%d  
                      write(p%UnLog, '(A12,f12.4)') " w   : ", m%LineTypeList(l)%w  
                      write(p%UnLog, '(A12,f12.4)') " Cdn : ", m%LineTypeList(l)%Cdn
                      write(p%UnLog, '(A12,f12.4)') " Can : ", m%LineTypeList(l)%Can
                      write(p%UnLog, '(A12,f12.4)') " Cdt : ", m%LineTypeList(l)%Cdt
                      write(p%UnLog, '(A12,f12.4)') " Cat : ", m%LineTypeList(l)%Cat
                  end if

                  IF ( ErrStat2 /= ErrID_None ) THEN
                     CALL SetErrStat( ErrID_Fatal, ErrMsg2, ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF

               END DO


            !-------------------------------------------------------------------------------------------
            else if ( (INDEX(Line, "ROD DICTIONARY") > 0) .or. ( INDEX(Line, "ROD TYPES") > 0) ) then ! if rod dictionary header
               
               IF (wordy > 0) print *, "Reading rod types"
               
               ! skip following two lines (label line and unit line)
               Line = NextLine(i)
               Line = NextLine(i)
               
                ! process each line
                DO l = 1,p%nRodTypes
                   
                   !read into a line
                   Line = NextLine(i)

                   ! check for correct number of columns in current line
                   IF ( CountWords( Line ) /= 7 ) THEN
                       CALL SetErrStat( ErrID_Fatal, ' Unable to parse Rod Type '//trim(Num2LStr(l))//' on row '//trim(Num2LStr(i))//' in input file. Row has wrong number of columns. Must be 7 columns.', ErrStat, ErrMsg, RoutineName )
                       CALL CleanUp()
                       RETURN
                   END IF
                   
                   ! parse out entries: Name  Diam MassDen Cd  Ca  CdEnd  CaEnd
                   IF (ErrStat2 == 0) THEN
                      READ(Line,*,IOSTAT=ErrStat2) m%RodTypeList(l)%name, m%RodTypeList(l)%d, m%RodTypeList(l)%w, &
                         m%RodTypeList(l)%Cdn, m%RodTypeList(l)%Can, m%RodTypeList(l)%CdEnd, m%RodTypeList(l)%CaEnd   
                         
                      m%RodTypeList(l)%Cdt = 0.0_DbKi ! not used
                      m%RodTypeList(l)%Cat = 0.0_DbKi ! not used
                   END IF

                   ! specify IdNum of rod type for error checking
                   m%RodTypeList(l)%IdNum = l  
                   
                  ! write lineType information to log file
                  if (p%writeLog > 1) then
                      write(p%UnLog, '(A12,A20)'  ) " RodType"//trim(num2lstr(l))//":"
                      write(p%UnLog, '(A12,A20)'  ) " name: ", m%RodTypeList(l)%name
                      write(p%UnLog, '(A12,f12.4)') " d   : ", m%RodTypeList(l)%d  
                      write(p%UnLog, '(A12,f12.4)') " w   : ", m%RodTypeList(l)%w  
                      write(p%UnLog, '(A12,f12.4)') " Cdn : ", m%RodTypeList(l)%Cdn
                      write(p%UnLog, '(A12,f12.4)') " Can : ", m%RodTypeList(l)%Can
                      write(p%UnLog, '(A12,f12.4)') " Cdt : ", m%RodTypeList(l)%CdEnd
                      write(p%UnLog, '(A12,f12.4)') " Cat : ", m%RodTypeList(l)%CaEnd
                  end if

                   IF ( ErrStat2 /= ErrID_None ) THEN
                      CALL SetErrStat( ErrID_Fatal, 'Failed to process rod type properties for rod '//trim(Num2LStr(l))//'. Check formatting and correct number of columns.', ErrStat, ErrMsg, RoutineName )
                      CALL CleanUp()
                      RETURN
                   END IF

                END DO


            !-------------------------------------------------------------------------------------------
            else if ((INDEX(Line, "BODIES") > 0 ) .or. (INDEX(Line, "BODY LIST") > 0 ) .or. (INDEX(Line, "BODY PROPERTIES") > 0 )) then

               ! skip following two lines (label line and unit line)
               Line = NextLine(i)
               Line = NextLine(i)
               
               ! process each body
               DO l = 1,p%nBodies
                  
                  !read into a line
                  Line = NextLine(i)
                   
                  ! check for correct number of columns in current line
                  IF ( CountWords( Line ) /= 14 ) THEN
                      CALL SetErrStat( ErrID_Fatal, ' Unable to parse Body '//trim(Num2LStr(l))//' on row '//trim(Num2LStr(i))//' in input file. Row has wrong number of columns. Must be 14 columns.', ErrStat, ErrMsg, RoutineName )
                      CALL CleanUp()
                      RETURN
                  END IF

                  ! parse out entries: ID   Attachment  X0  Y0  Z0  r0  p0  y0    M  CG*  I*    V  CdA*  Ca*
                  IF (ErrStat2 == 0) THEN
                     READ(Line,*,IOSTAT=ErrStat2) m%BodyList(l)%IdNum, tempString1, &
                        tempArray(1), tempArray(2), tempArray(3), tempArray(4), tempArray(5), tempArray(6), &
                        m%BodyList(l)%bodyM, tempString2, tempString3, m%BodyList(l)%bodyV, tempString4, tempString5
                  END IF
                  
                  ! process CG
                  CALL SplitByBars(tempString2, N, tempStrings)
                  if (N == 1) then                                   ! if only one entry, it is the z coordinate
                     m%BodyList(l)%rCG(1) = 0.0_DbKi
                     m%BodyList(l)%rCG(2) = 0.0_DbKi
                     READ(tempString2, *) m%BodyList(l)%rCG(3)
                  else if (N==3) then                                ! all three coordinates provided
                     READ(tempStrings(1), *) m%BodyList(l)%rCG(1)
                     READ(tempStrings(2), *) m%BodyList(l)%rCG(2)
                     READ(tempStrings(3), *) m%BodyList(l)%rCG(3)
                  else
                     CALL SetErrStat( ErrID_Fatal, 'Body '//trim(Num2LStr(l))//' CG entry (col 10) must have 1 or 3 numbers.' , ErrStat, ErrMsg, RoutineName )
                  end if
                  ! process mements of inertia
                  CALL SplitByBars(tempString3, N, tempStrings)
                  if (N == 1) then                                   ! if only one entry, use it for all directions
                     READ(tempString3, *) m%BodyList(l)%BodyI(1)
                     m%BodyList(l)%BodyI(2) = m%BodyList(l)%BodyI(1)
                     m%BodyList(l)%BodyI(3) = m%BodyList(l)%BodyI(1)
                  else if (N==3) then                                ! all three directions provided separately
                     READ(tempStrings(1), *) m%BodyList(l)%BodyI(1)
                     READ(tempStrings(2), *) m%BodyList(l)%BodyI(2)
                     READ(tempStrings(3), *) m%BodyList(l)%BodyI(3)
                  else
                     CALL SetErrStat( ErrID_Fatal, 'Body '//trim(Num2LStr(l))//' inertia entry (col 11) must have 1 or 3 numbers.' , ErrStat, ErrMsg, RoutineName )
                  end if
                  ! process drag ceofficient by area product
                  CALL SplitByBars(tempString4, N, tempStrings)
                  if (N == 1) then                                   ! if only one entry, use it for all directions
                     READ(tempString4, *) m%BodyList(l)%BodyCdA(1)
                     m%BodyList(l)%BodyCdA(2) = m%BodyList(l)%BodyCdA(1)
                     m%BodyList(l)%BodyCdA(3) = m%BodyList(l)%BodyCdA(1)
                  else if (N==3) then                                ! all three coordinates provided
                     READ(tempStrings(1), *) m%BodyList(l)%BodyCdA(1)
                     READ(tempStrings(2), *) m%BodyList(l)%BodyCdA(2)
                     READ(tempStrings(3), *) m%BodyList(l)%BodyCdA(3)
                  else
                     CALL SetErrStat( ErrID_Fatal, 'Body '//trim(Num2LStr(l))//' CdA entry (col 13) must have 1 or 3 numbers.' , ErrStat, ErrMsg, RoutineName )
                  end if
                  ! process added mass coefficient
                  CALL SplitByBars(tempString5, N, tempStrings)
                  if (N == 1) then                                   ! if only one entry, use it for all directions
                     READ(tempString5, *) m%BodyList(l)%BodyCa(1)
                     m%BodyList(l)%BodyCa(2) = m%BodyList(l)%BodyCa(1)
                     m%BodyList(l)%BodyCa(3) = m%BodyList(l)%BodyCa(1)
                  else if (N==3) then                                ! all three coordinates provided
                     READ(tempStrings(1), *) m%BodyList(l)%BodyCa(1)
                     READ(tempStrings(2), *) m%BodyList(l)%BodyCa(2)
                     READ(tempStrings(3), *) m%BodyList(l)%BodyCa(3)
                  else
                     CALL SetErrStat( ErrID_Fatal, 'Body '//trim(Num2LStr(l))//' Ca entry (col 14) must have 1 or 3 numbers.' , ErrStat, ErrMsg, RoutineName )
                  end if


                  IF ( ErrStat2 /= 0 ) THEN
                     CALL WrScr('   Unable to parse Body '//trim(Num2LStr(l))//' on row '//trim(Num2LStr(i))//' in input file.')  ! Specific screen output because errors likely
                     CALL WrScr('   Ensure row has all 13 columns needed in MDv2 input file (13th Dec 2021).')  
                        CALL SetErrStat( ErrID_Fatal, 'Failed to read bodies.' , ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF
                  

                  !----------- process body type -----------------

                  call DecomposeString(tempString1, let1, num1, let2, num2, let3)   ! note: this call is overkill (it's just a string) but leaving it here for potential future expansions
                  
                  if ((let1 == "ANCHOR") .or. (let1 == "FIXED") .or. (let1 == "FIX")) then   ! if a fixed body (this would just be used if someone wanted to temporarly fix a body that things were attached to)
                  
                     m%BodyList(l)%typeNum = 1
                     
                  else if ((let1 == "COUPLED") .or. (let1 == "VESSEL") .or. (let1 == "CPLD") .or. (let1 == "VES")) then    ! if a coupled body
                     
                     m%BodyList(l)%typeNum = -1
                     p%nCpldBodies(1)=p%nCpldBodies(1)+1  ! add this body to coupled list                          
                     m%CpldBodyIs(p%nCpldBodies(1),1) = l

                     ! body initial position due to coupling will be adjusted later
                     
                  ! TODO: add option for body coupling to different turbines in FAST.Farm  <<<
                     
                  else if (let1 == "FREE") then    ! if a free body
                     m%BodyList(l)%typeNum = 0
                     
                     p%nFreeBodies=p%nFreeBodies+1
                     
                     m%BodyStateIs1(p%nFreeBodies) = Nx+1
                     m%BodyStateIsN(p%nFreeBodies) = Nx+12
                     Nx = Nx + 12                           ! add 12 state variables for free Body
                     
                     m%FreeBodyIs(p%nFreeBodies) = l
                     
                     m%BodyList(l)%r6 = tempArray     ! set initial body position and orientation
                     
                  else 
                     CALL SetErrStat( ErrID_Fatal,  "Unidentified Body type string for Body "//trim(Num2LStr(l))//": "//trim(tempString1), ErrStat, ErrMsg, RoutineName )
                     return
                  end if
                  

                  ! check for sequential IdNums
                  IF ( m%BodyList(l)%IdNum .NE. l ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Body numbers must be sequential starting from 1.', ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF
                  
                  
                  ! set up body
                  CALL Body_Setup( m%BodyList(l), tempArray, p, ErrStat2, ErrMsg2)
                     CALL CheckError( ErrStat2, ErrMsg2 )
                     IF (ErrStat >= AbortErrLev) RETURN

                  IF ( ErrStat2 /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Failed to read data for body '//trim(Num2LStr(l)), ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF
                  
                  IF (wordy > 1) print *, "Set up body ", l, " of type ",  m%BodyList(l)%typeNum

               END DO
               
               
            !-------------------------------------------------------------------------------------------
            else if ((INDEX(Line, "RODS") > 0 ) .or. (INDEX(Line, "ROD LIST") > 0) .or. (INDEX(Line, "ROD PROPERTIES") > 0)) then ! if rod properties header

               IF (wordy > 0) print *, "Reading Rods"
               
               ! skip following two lines (label line and unit line)
               Line = NextLine(i)
               Line = NextLine(i)
               
               ! process each rod
               DO l = 1,p%nRods
                  
                  !read into a line
                  Line = NextLine(i)

                  ! check for correct number of columns in current line
                  IF ( CountWords( Line ) /= 11 ) THEN
                      CALL SetErrStat( ErrID_Fatal, ' Unable to parse Rod '//trim(Num2LStr(l))//' on row '//trim(Num2LStr(i))//' in input file. Row has wrong number of columns. Must be 11 columns.', ErrStat, ErrMsg, RoutineName )
                      CALL CleanUp()
                      RETURN
                  END IF
                  
                  ! parse out entries: RodID  RodType  Attachment  Xa   Ya   Za   Xb   Yb   Zb  NumSegs  Flags/Outputs
                  IF (ErrStat2 == 0) THEN
                     READ(Line,*,IOSTAT=ErrStat2) m%RodList(l)%IdNum, tempString1, tempString2, &
                         tempArray(1), tempArray(2), tempArray(3), tempArray(4), tempArray(5), tempArray(6), &
                         m%RodList(l)%N, LineOutString
                  END IF

                  ! find Rod properties index
                  DO J = 1,p%nRodTypes
                     IF (trim(tempString1) == trim(m%RodTypeList(J)%name)) THEN
                       m%RodList(l)%PropsIdNum = J
                       EXIT
                     END IF
                     IF (J == p%nRodTypes) THEN   ! call an error if there is no match
                         CALL SetErrStat( ErrID_Fatal, 'Unable to find matching rod type name for Rod '//trim(Num2LStr(l))//": "//trim(tempString1), ErrStat, ErrMsg, RoutineName )
                         RETURN
                     END IF
                  END DO


                  !----------- process rod type -----------------

                  call DecomposeString(tempString2, let1, num1, let2, num2, let3)
                  
                  if ((let1 == "ANCHOR") .or. (let1 == "FIXED") .or. (let1 == "FIX")) then
                     
                     m%RodList(l)%typeNum = 2
                     CALL Body_AddRod(m%GroundBody, l, tempArray)   ! add rod l to Ground body
                           

                     else if ((let1 == "PINNED") .or. (let1 == "PIN")) then
                        m%RodList(l)%typeNum = 1
                     CALL Body_AddRod(m%GroundBody, l, tempArray)   ! add rod l to Ground body
                     
                     p%nFreeRods=p%nFreeRods+1  ! add this pinned rod to the free list because it is half free
                     
                     m%RodStateIs1(p%nFreeRods) = Nx+1
                     m%RodStateIsN(p%nFreeRods) = Nx+6
                     Nx = Nx + 6                                               ! add 6 state variables for each pinned rod
                     
                     m%FreeRodIs(p%nFreeRods) = l
                     
                  else if (let1 == "BODY") then ! attached to a body (either rididly or pinned)
                  
                     if (len_trim(num1) > 0) then
                     
                        READ(num1,*) J   ! convert to int, representing parent body index
                        
                        if ((J <= p%nBodies) .and. (J > 0)) then 
                        
                           CALL Body_AddRod(m%BodyList(J), l, tempArray)   ! add rod l to the body
                           
                           if ( (let2 == "PINNED") .or. (let2 == "PIN") ) then
                              m%RodList(l)%typeNum = 1
                              
                              p%nFreeRods=p%nFreeRods+1  ! add this pinned rod to the free list because it is half free
                              
                              m%RodStateIs1(p%nFreeRods) = Nx+1
                              m%RodStateIsN(p%nFreeRods) = Nx+6
                              Nx = Nx + 6                                               ! add 6 state variables for each pinned rod
                              
                              m%FreeRodIs(p%nFreeRods) = l
                     
                           else if (let2 == " ") then ! rod is not requested to be pinned, so add this rod as a fixed one
                              m%RodList(l)%typeNum = 2
                     
                           else
                               CALL SetErrStat( ErrID_Fatal,  "Unidentified Type/BodyID for Rod "//trim(Num2LStr(l))//": "//trim(tempString2), ErrStat, ErrMsg, RoutineName )
                               return
                           end if
                        
                        else
                           CALL SetErrStat( ErrID_Fatal,  "Body ID out of bounds for Rod "//trim(Num2LStr(l))//".", ErrStat, ErrMsg, RoutineName )  
                           return
                        end if
                     
                  else
                     CALL SetErrStat( ErrID_Fatal,  "No number provided for Rod "//trim(Num2LStr(l))//" Body attachment.", ErrStat, ErrMsg, RoutineName )   
                         return
                  end if
                  
                  else if ((let1 == "VESSEL") .or. (let1 == "VES") .or. (let1 == "COUPLED") .or. (let1 == "CPLD")) then    ! if a rigidly coupled rod, add to list and add 
                     m%RodList(l)%typeNum = -2            

                     p%nCpldRods(1)=p%nCpldRods(1)+1     ! add this rod to coupled list
                     
                     m%CpldRodIs(p%nCpldRods(1),1) = l
                 
                  else if ((let1 == "VESSELPINNED") .or. (let1 == "VESPIN") .or. (let1 == "COUPLEDPINNED") .or. (let1 == "CPLDPIN")) then  ! if a pinned coupled rod, add to list and add 
                     m%RodList(l)%typeNum = -1
                     
                     p%nCpldRods(1)=p%nCpldRods(1)+1  ! add
                     p%nFreeRods   =p%nFreeRods+1     ! add this pinned rod to the free list because it is half free
                     
                     m%RodStateIs1(p%nFreeRods) = Nx+1
                     m%RodStateIsN(p%nFreeRods) = Nx+6
                     Nx = Nx + 6                                               ! add 6 state variables for each pinned rod
                     
                     m%CpldRodIs(p%nCpldRods(1),1) = l
                     m%FreeRodIs(p%nFreeRods) = l
                   
                  ! TODO: add option for body coupling to different turbines in FAST.Farm <<<
                   
                  else if ((let1 == "CONNECT") .or. (let1 == "CON") .or. (let1 == "FREE")) then
                     m%RodList(l)%typeNum = 0
                     
                     p%nFreeRods=p%nFreeRods+1  ! add this pinned rod to the free list because it is half free
                     
                     m%RodStateIs1(p%nFreeRods) = Nx+1
                     m%RodStateIsN(p%nFreeRods) = Nx+12
                     Nx = Nx + 12                                              ! add 12 state variables for free Rod
                     
                     m%FreeRodIs(p%nFreeRods) = l
                     
                  else 
                  
                     CALL SetErrStat( ErrID_Fatal,  "Unidentified Type/BodyID for Rod "//trim(Num2LStr(l))//": "//trim(tempString2), ErrStat, ErrMsg, RoutineName )   
                     return
                  end if
                  
                
                  ! process output flag characters (LineOutString) and set line output flag array (OutFlagList)
                  m%RodList(l)%OutFlagList = 0  ! first set array all to zero
                  ! per node, 3 component
                  IF ( scan( LineOutString, 'p') > 0 )  m%RodList(l)%OutFlagList(2 ) = 1   ! node position
                  IF ( scan( LineOutString, 'v') > 0 )  m%RodList(l)%OutFlagList(3 ) = 1   ! node velocity
                  IF ( scan( LineOutString, 'U') > 0 )  m%RodList(l)%OutFlagList(4 ) = 1   ! water velocity
                  IF ( scan( LineOutString, 'B') > 0 )  m%RodList(l)%OutFlagList(5 ) = 1   ! node buoyancy force
                  IF ( scan( LineOutString, 'D') > 0 )  m%RodList(l)%OutFlagList(6 ) = 1   ! drag force
                  IF ( scan( LineOutString, 'I') > 0 )  m%RodList(l)%OutFlagList(7 ) = 1   ! inertia force
                  IF ( scan( LineOutString, 'P') > 0 )  m%RodList(l)%OutFlagList(8 ) = 1   ! dynamic pressure force
                  IF ( scan( LineOutString, 'b') > 0 )  m%RodList(l)%OutFlagList(9 ) = 1   ! seabed contact forces
                  ! per node, 1 component
                  IF ( scan( LineOutString, 'W') > 0 )  m%RodList(l)%OutFlagList(10) = 1   ! node weight/buoyancy (positive up)
                  IF ( scan( LineOutString, 'K') > 0 )  m%RodList(l)%OutFlagList(11) = 1   ! curvature at node
                  ! per element, 1 component >>> these don't apply to a rod!! <<<
                  IF ( scan( LineOutString, 't') > 0 )  m%RodList(l)%OutFlagList(12) = 1  ! segment tension force (just EA)
                  IF ( scan( LineOutString, 'c') > 0 )  m%RodList(l)%OutFlagList(13) = 1  ! segment internal damping force
                  IF ( scan( LineOutString, 's') > 0 )  m%RodList(l)%OutFlagList(14) = 1  ! Segment strain
                  IF ( scan( LineOutString, 'd') > 0 )  m%RodList(l)%OutFlagList(15) = 1  ! Segment strain rate

                  IF (SUM(m%RodList(l)%OutFlagList) > 0)   m%RodList(l)%OutFlagList(1) = 1  ! this first entry signals whether to create any output file at all
                  ! the above letter-index combinations define which OutFlagList entry corresponds to which output type


                  ! specify IdNum of line for error checking
                  m%RodList(l)%IdNum = l  

                  ! check for sequential IdNums
                  IF ( m%RodList(l)%IdNum .NE. l ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Line numbers must be sequential starting from 1.', ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF
                  
                  ! set up rod
                  CALL Rod_Setup( m%RodList(l), m%RodTypeList(m%RodList(l)%PropsIdNum), tempArray, p, ErrStat2, ErrMsg2)
                     CALL CheckError( ErrStat2, ErrMsg2 )
                     IF (ErrStat >= AbortErrLev) RETURN
                     
                  ! note: Rod was already added to its respective parent body if type > 0
                  
                  IF ( ErrStat2 /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Failed to read rod data for Rod '//trim(Num2LStr(l)), ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF

               END DO   ! l = 1,p%nRods


            !-------------------------------------------------------------------------------------------
            else if ((INDEX(Line, "POINTS") > 0 ) .or. (INDEX(Line, "CONNECTION PROPERTIES") > 0) .or. (INDEX(Line, "NODE PROPERTIES") > 0) .or. (INDEX(Line, "POINT PROPERTIES") > 0) .or. (INDEX(Line, "POINT LIST") > 0) ) then ! if node properties header
               
               IF (wordy > 0) print *, "Reading Points"
               
               ! skip following two lines (label line and unit line)
               Line = NextLine(i)
               Line = NextLine(i)
               
               ! process each point
               DO l = 1,p%nConnects
                  
                  !read into a line
                  Line = NextLine(i)

                  ! check for correct number of columns in current line
                  IF ( CountWords( Line ) /= 9 ) THEN
                      CALL SetErrStat( ErrID_Fatal, ' Unable to parse Point '//trim(Num2LStr(l))//' on row '//trim(Num2LStr(i))//' in input file. Row has wrong number of columns. Must be 9 columns.', ErrStat, ErrMsg, RoutineName )
                      CALL CleanUp()
                      RETURN
                  END IF
                  
                  ! parse out entries: PointID Attachment  X  Y  Z  M  V  CdA Ca 
                  IF (ErrStat2 == 0) THEN
                     READ(Line,*,IOSTAT=ErrStat2) m%ConnectList(l)%IdNum, tempString1, tempArray(1), &
                        tempArray(2), tempString4, m%ConnectList(l)%conM, &
                        m%ConnectList(l)%conV, m%ConnectList(l)%conCdA, m%ConnectList(l)%conCa
                                          
                     CALL Conv2UC(tempString4) ! convert to uppercase so that matching is not case-sensitive
                     
                     if ((INDEX(tempString4, "SEABED") > 0 ) .or. (INDEX(tempString4, "GROUND") > 0 ) .or. (INDEX(tempString4, "FLOOR") > 0 )) then  ! if keyword used
                        CALL WrScr('Point '//trim(Num2LStr(l))//' depth set to be on the seabed; finding z location based on depth/bathymetry')      ! interpret the anchor depth value as a 'seabed' input
                        CALL getDepthFromBathymetry(m%BathymetryGrid, m%BathGrid_Xs, m%BathGrid_Ys, tempArray(1), tempArray(2), depth, nvec)         ! meaning the anchor should be at the depth of the local bathymetry
                        tempArray(3) = -depth
                     else                                                       ! if the anchor depth input isn't one of the supported keywords, 
                        READ(tempString4, *, IOSTAT=ErrStat2) tempArray(3)      ! assume it's a scalar depth value
                        !TODO: add error check for if the above read fails
                     end if
                        
                     ! not used
                     m%ConnectList(l)%conFX = 0.0_DbKi 
                     m%ConnectList(l)%conFY = 0.0_DbKi
                     m%ConnectList(l)%conFZ = 0.0_DbKi
                        
                  END IF
                  

                  IF ( ErrStat2 /= 0 ) THEN
                     CALL WrScr('   Unable to parse Point '//trim(Num2LStr(l))//' row in input file.')  ! Specific screen output because errors likely
                     CALL WrScr('   Ensure row has all 9 columns, including CdA and Ca.')           ! to be caused by non-updated input file formats.
                        CALL SetErrStat( ErrID_Fatal, 'Failed to read connects.' , ErrStat, ErrMsg, RoutineName ) ! would be nice to specify which line <<<<<<<<<
                     CALL CleanUp()
                     RETURN
                  END IF
                  
                  m%ConnectList(l)%r = tempArray(1:3)   ! set initial, or reference, node position (for coupled or child objects, this will be the local reference location about the parent)

                  !----------- process connection type -----------------

                  call DecomposeString(tempString1, let1, num1, let2, num2, let3)
                  
                     if ((let1 == "ANCHOR") .or. (let1 == "FIXED") .or. (let1 == "FIX")) then
                         m%ConnectList(l)%typeNum = 1
                     
                     !m%ConnectList(l)%r = tempArray(1:3)   ! set initial node position
                     
                     CALL Body_AddConnect(m%GroundBody, l, tempArray(1:3))   ! add connection l to Ground body                     

                     else if (let1 == "BODY") then ! attached to a body
                     if (len_trim(num1) > 0) then                     
                        READ(num1, *) J   ! convert to int, representing parent body index
                        
                        if ((J <= p%nBodies) .and. (J > 0)) then
                           m%ConnectList(l)%typeNum = 1    

                           CALL Body_AddConnect(m%BodyList(J), l, tempArray(1:3))   ! add connection l to Ground body
                           
                        else
                           CALL SetErrStat( ErrID_Fatal,  "Body ID out of bounds for Connection "//trim(Num2LStr(l))//".", ErrStat, ErrMsg, RoutineName )  
                           return
                        end if                     
                     else
                        CALL SetErrStat( ErrID_Fatal,  "No number provided for Connection "//trim(Num2LStr(l))//" Body attachment.", ErrStat, ErrMsg, RoutineName )   
                            return
                     end if
                  
                  else if ((let1 == "VESSEL") .or. (let1 == "VES") .or. (let1 == "COUPLED") .or. (let1 == "CPLD")) then    ! if a fairlead, add to list and add 
                     m%ConnectList(l)%typeNum = -1
                     p%nCpldCons(1)=p%nCpldCons(1)+1                       
                     m%CpldConIs(p%nCpldCons(1),1) = l
                 
                  else if ((let1 == "CONNECT") .or. (let1 == "CON") .or. (let1 == "FREE")) then
                     m%ConnectList(l)%typeNum = 0
                     
                     p%nFreeCons=p%nFreeCons+1             ! add this pinned rod to the free list because it is half free
                     
                     m%ConStateIs1(p%nFreeCons) = Nx+1
                     m%ConStateIsN(p%nFreeCons) = Nx+6
                     Nx = Nx + 6                           ! add 12 state variables for free Connection
                     
                     m%FreeConIs(p%nFreeCons) = l
                     
                     !m%ConnectList(l)%r = tempArray(1:3)   ! set initial node position
                     
                  else if ((let1 == "TURBINE") .or. (let1 == "T")) then  ! turbine-coupled in FAST.Farm case
                  
                     if (len_trim(num1) > 0) then                     
                        READ(num1, *) J   ! convert to int, representing turbine index
                        
                        if ((J <= p%nTurbines) .and. (J > 0)) then
                           
                           m%ConnectList(l)%TypeNum = -1            ! set as coupled type   
                           p%nCpldCons(J) = p%nCpldCons(J) + 1      ! increment counter for the appropriate turbine                   
                           m%CpldConIs(p%nCpldCons(J),J) = l
                           CALL WrScr(' added connection '//TRIM(int2lstr(l))//' as fairlead for turbine '//trim(int2lstr(J)))
                           
                           
                        else
                           CALL SetErrStat( ErrID_Fatal,  "Turbine ID out of bounds for Connection "//trim(Num2LStr(l))//".", ErrStat, ErrMsg, RoutineName )  
                           return
                        end if   
                     else
                        CALL SetErrStat( ErrID_Fatal,  "No number provided for Connection "//trim(Num2LStr(l))//" Turbine attachment.", ErrStat, ErrMsg, RoutineName )   
                            return
                     end if
                  
                  else 
                     CALL SetErrStat( ErrID_Fatal,  "Unidentified Type/BodyID for Connection "//trim(Num2LStr(l))//": "//trim(tempString1), ErrStat, ErrMsg, RoutineName )
                     return
                  end if
                  
                  ! set initial velocity to zero
                  m%ConnectList(l)%rd(1) = 0.0_DbKi
                  m%ConnectList(l)%rd(2) = 0.0_DbKi
                  m%ConnectList(l)%rd(3) = 0.0_DbKi
                           
                  !also set number of attached lines to zero initially
                  m%ConnectList(l)%nAttached = 0


                  ! check for sequential IdNums
                  IF ( m%ConnectList(l)%IdNum .NE. l ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Connection numbers must be sequential starting from 1.', ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF
                  

                  IF ( ErrStat2 /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Failed to read data for Connection '//trim(Num2LStr(l)), ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF
                  
                  IF (wordy > 0) print *, "Set up Point ", l, " of type ",  m%ConnectList(l)%typeNum

               END DO   ! l = 1,p%nRods

            !-------------------------------------------------------------------------------------------
            else if ((INDEX(Line, "LINES") > 0 ) .or. (INDEX(Line, "LINE PROPERTIES") > 0) .or. (INDEX(Line, "LINE LIST") > 0) ) then ! if line properties header

               IF (wordy > 0) print *, "Reading Lines"
               
               ! skip following two lines (label line and unit line)
               Line = NextLine(i)
               Line = NextLine(i)
               
               ! process each line
               DO l = 1,p%nLines
                  
                  !read into a line
                  Line = NextLine(i)

                  ! check for correct number of columns in current line
                  IF ( CountWords( Line ) /= 7 ) THEN
                      CALL SetErrStat( ErrID_Fatal, ' Unable to parse Line '//trim(Num2LStr(l))//' on row '//trim(Num2LStr(i))//' in input file. Row has wrong number of columns. Must be 7 columns.', ErrStat, ErrMsg, RoutineName )
                      CALL CleanUp()
                      RETURN
                  END IF
                  
                   ! parse out entries: ID  LineType  AttachA  AttachB     UnstrLen  NumSegs   Outputs  (note: order changed Dec 13, 2021 before MDv2 release) 
                  IF (ErrStat2 == 0) THEN
                     READ(Line,*,IOSTAT=ErrStat2) m%LineList(l)%IdNum, tempString1, tempString2, tempString3, &
                        m%LineList(l)%UnstrLen, m%LineList(l)%N,  LineOutString
                  END IF
                                    
                  ! identify index of line type
                  DO J = 1,p%nLineTypes
                     IF (trim(tempString1) == trim(m%LineTypeList(J)%name)) THEN
                       m%LineList(l)%PropsIdNum = J
                       EXIT
                       END IF
                       IF (J == p%nLineTypes) THEN   ! call an error if there is no match
                           CALL SetErrStat( ErrID_Fatal, 'Unable to find matching line type name for Line '//trim(Num2LStr(l)), ErrStat, ErrMsg, RoutineName )
                           RETURN
                       END IF
                  END DO
                  
                  ! account for states of line
                  m%LineStateIs1(l) = Nx + 1
                  if (m%LineTypeList(m%LineList(l)%PropsIdNum)%ElasticMod == 2) then
                     Nx = Nx + 7*m%LineList(l)%N - 6       ! if using viscoelastic model, need one more state per segment
                     m%LineStateIsN(l) = Nx          
                  else
                     Nx = Nx + 6*m%LineList(l)%N - 6       ! normal case, just 6 states per internal node   
                     m%LineStateIsN(l) = Nx          
                  end if
                  
                  ! Process attachment identfiers and attach line ends 
                  
                  ! First for the anchor (or end A)...
                  
                  call DecomposeString(tempString2, let1, num1, let2, num2, let3)
                  
                  if (len_trim(num1)<1) then
                     CALL SetErrStat( ErrID_Fatal,  "Error: no number provided for line "//trim(Num2LStr(l))//" end A attachment.", ErrStat, ErrMsg, RoutineName )  
                     return
                  end if 

                  READ(num1, *) J   ! convert to int
                     
                  ! if id starts with an "R" or "Rod"
                  if ((let1 == "R") .or. (let1 == "ROD")) then   
                  
                     if ((J <= p%nRods) .and. (J > 0)) then                  
                        if (let2 == "A") then
                           CALL Rod_AddLine(m%RodList(J), l, 0, 0)   ! add line l (end A, denoted by 0) to rod J (end A, denoted by 0)
                        else if (let2 == "B") then 
                           CALL Rod_AddLine(m%RodList(J), l, 0, 1)   ! add line l (end A, denoted by 0) to rod J (end B, denoted by 1)
                        else
                           CALL SetErrStat( ErrID_Fatal,  "Error: rod end (A or B) must be specified for line "//trim(Num2LStr(l))//" end A attachment. Instead seeing "//let2, ErrStat, ErrMsg, RoutineName )  
                            return
                        end if
                     else
                        CALL SetErrStat( ErrID_Fatal,  "Error: rod connection ID out of bounds for line "//trim(Num2LStr(l))//" end A attachment.", ErrStat, ErrMsg, RoutineName )  
                        return
                     end if
                  
                     ! if J starts with a "C" or "Con" or goes straight ot the number then it's attached to a Connection
                  else if ((len_trim(let1)==0) .or. (let1 == "C") .or. (let1 == "CON")) then 

                     if ((J <= p%nConnects) .and. (J > 0)) then                  
                        CALL Connect_AddLine(m%ConnectList(J), l, 0)   ! add line l (end A, denoted by 0) to connection J
                     else
                        CALL SetErrStat( ErrID_Fatal,  "Error: connection out of bounds for line "//trim(Num2LStr(l))//" end A attachment.", ErrStat, ErrMsg, RoutineName )  
                        return
                     end if
                        
                  end if


                  ! Then again for the fairlead (or end B)...

                  call DecomposeString(tempString3, let1, num1, let2, num2, let3)

                  if (len_trim(num1)<1) then
                     CALL SetErrStat( ErrID_Fatal,  "Error: no number provided for line "//trim(Num2LStr(l))//" end B attachment.", ErrStat, ErrMsg, RoutineName )  
                     return
                  end if 

                  READ(num1, *) J   ! convert to int
                     
                  ! if id starts with an "R" or "Rod"
                  if ((let1 == "R") .or. (let1 == "ROD")) then   

                     if ((J <= p%nRods) .and. (J > 0)) then                  
                        if (let2 == "A") then
                           CALL Rod_AddLine(m%RodList(J), l, 1, 0)   ! add line l (end B, denoted by 1) to rod J (end A, denoted by 0)
                        else if (let2 == "B") then 
                           CALL Rod_AddLine(m%RodList(J), l, 1, 1)   ! add line l (end B, denoted by 1) to rod J (end B, denoted by 1)
                        else
                           CALL SetErrStat( ErrID_Fatal,  "Error: rod end (A or B) must be specified for line "//trim(Num2LStr(l))//" end B attachment. Instead seeing "//let2, ErrStat, ErrMsg, RoutineName )  
                            return
                        end if
                     else
                        CALL SetErrStat( ErrID_Fatal,  "Error: rod connection ID out of bounds for line "//trim(Num2LStr(l))//" end B attachment.", ErrStat, ErrMsg, RoutineName )  
                        return
                     end if

                  ! if J starts with a "C" or "Con" or goes straight ot the number then it's attached to a Connection
                  else if ((len_trim(let1)==0) .or. (let1 == "C") .or. (let1 == "CON")) then 

                     if ((J <= p%nConnects) .and. (J > 0)) then                  
                        CALL Connect_AddLine(m%ConnectList(J), l, 1)   ! add line l (end B, denoted by 1) to connection J
                     else
                        CALL SetErrStat( ErrID_Fatal,  "Error: connection out of bounds for line "//trim(Num2LStr(l))//" end B attachment.", ErrStat, ErrMsg, RoutineName )  
                        return
                     end if
                        
                  end if
                
                
                  ! process output flag characters (LineOutString) and set line output flag array (OutFlagList)
                  m%LineList(l)%OutFlagList = 0  ! first set array all to zero
                  ! per node 3 component
                  IF ( scan( LineOutString, 'p') > 0 )  m%LineList(l)%OutFlagList(2) = 1 
                  IF ( scan( LineOutString, 'v') > 0 )  m%LineList(l)%OutFlagList(3) = 1
                  IF ( scan( LineOutString, 'U') > 0 )  m%LineList(l)%OutFlagList(4) = 1
                  IF ( scan( LineOutString, 'D') > 0 )  m%LineList(l)%OutFlagList(5) = 1
                  IF ( scan( LineOutString, 'b') > 0 )  m%LineList(l)%OutFlagList(6) = 1   ! seabed contact forces
                  ! per node 1 component
                  IF ( scan( LineOutString, 'W') > 0 )  m%LineList(l)%OutFlagList(7) = 1  ! node weight/buoyancy (positive up)
                  IF ( scan( LineOutString, 'K') > 0 )  m%LineList(l)%OutFlagList(8) = 1  ! curvature at node
                  ! per element 1 component
                  IF ( scan( LineOutString, 't') > 0 )  m%LineList(l)%OutFlagList(10) = 1  ! segment tension force (just EA)
                  IF ( scan( LineOutString, 'c') > 0 )  m%LineList(l)%OutFlagList(11) = 1  ! segment internal damping force
                  IF ( scan( LineOutString, 's') > 0 )  m%LineList(l)%OutFlagList(12) = 1  ! Segment strain
                  IF ( scan( LineOutString, 'd') > 0 )  m%LineList(l)%OutFlagList(13) = 1  ! Segment strain rate
                  IF ( scan( LineOutString, 'l') > 0 )  m%LineList(l)%OutFlagList(14) = 1  ! Segment stretched length

                  IF (SUM(m%LineList(l)%OutFlagList) > 0)   m%LineList(l)%OutFlagList(1) = 1  ! this first entry signals whether to create any output file at all
                  ! the above letter-index combinations define which OutFlagList entry corresponds to which output type


                  ! specify IdNum of line for error checking
                  m%LineList(l)%IdNum = l  


                  ! check for sequential IdNums
                  IF ( m%LineList(l)%IdNum .NE. l ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Line numbers must be sequential starting from 1.', ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF
                  
                  
                  ! setup line
                  CALL SetupLine( m%LineList(l), m%LineTypeList(m%LineList(l)%PropsIdNum), p, ErrStat2, ErrMsg2)
                     CALL CheckError( ErrStat2, ErrMsg2 )
                     IF (ErrStat >= AbortErrLev) RETURN
                     

                  IF ( ErrStat2 /= 0 ) THEN
                     CALL SetErrStat( ErrID_Fatal, 'Failed to read line data for Line '//trim(Num2LStr(l)), ErrStat, ErrMsg, RoutineName )
                     CALL CleanUp()
                     RETURN
                  END IF

               END DO   ! l = 1,p%nLines



            !-------------------------------------------------------------------------------------------
            else if (INDEX(Line, "CONTROL") > 0) then ! if control inputs header

               IF (wordy > 0) print *, "   Reading control inputs";
               
               ! TODO: add stuff <<<<<<<<

               ! skip following two lines (label line and unit line)
               Line = NextLine(i)
               Line = NextLine(i)
               
               ! process each line
               DO l = 1,p%nCtrlChans
                  
                  !read into a line
                  Line = NextLine(i)
                  
                  ! count commas to determine how many line IDs specified for this channel
                  N = count(transfer(Line, 'a', len(Line)) == ",") + 1   ! number of line IDs given
                  
                  ! parse out entries:        CtrlChan, LineIdNums
                  read(Line, *) Itemp, TempIDnums(1:N)                                   ! parse out each line ID
                  
                  DO J = 1,N
                     if (TempIDnums(J) <= p%nLines) then      ! ensure line ID is in range
                        if (m%LineList( TempIDnums(J) )%CtrlChan == 0) then      ! ensure line doesn't already have a CtrlChan assigned 
                           m%LineList( TempIDnums(J) )%CtrlChan = Itemp
                           CALL WrScr('Assigned Line '//TRIM(Int2LStr(TempIDnums(J)))//' to control channel '//TRIM(Int2LStr(Itemp)))
                        else
                           CALL WrScr('Error: Line '//TRIM(Int2LStr(TempIDnums(J)))//' already is assigned to control channel '//TRIM(Int2LStr(m%LineList( TempIDnums(J) )%CtrlChan))//' so cannot also be assigned to channel '//TRIM(Int2LStr(Itemp)))
                        end if                     
                     else
                        CALL WrScr('Error: Line ID '//TRIM(Int2LStr(TempIDnums(J)))//' of CtrlChan '//TRIM(Int2LStr(Itemp))//' is out of range')
                     end if
                  
                  END DO
                  
               END DO
               

            !-------------------------------------------------------------------------------------------
            else if (INDEX(Line, "FAILURE") > 0) then ! if failure conditions header

               CALL WrScr("   Warning: Failure capabilities are not yet implemented in MoorDyn.")
               
               ! skip following two lines (label line and unit line)
               Line = NextLine(i)
               Line = NextLine(i)
               
               ! process each line
               DO l = 1,p%nFails
                  
                  !read into a line
                  Line = NextLine(i)
                  
                  ! TODO: Failure capabilities still need to be completed
                  READ(Line,*,IOSTAT=ErrStat2) m%LineList(l)%IdNum, tempString1, m%LineList(l)%UnstrLen, &
                       m%LineList(l)%N, tempString2, tempString3, LineOutString
                  
               END DO
               
               
            !-------------------------------------------------------------------------------------------
            else if (INDEX(Line, "OUTPUT") > 0) then ! if output header

               IF (wordy > 0) print *, "Reading Outputs"
               
               ! (don't skip any lines)
                        
               ! allocate InitInp%Outliest (to a really big number for now...)
               CALL AllocAry( OutList, MaxAryLen, "MoorDyn Input File's Outlist", ErrStat2, ErrMsg2 ); if(Failed()) return

 
               ! Initialize some values
               p%NumOuts = 0    ! start counter at zero
               OutList = ''


               ! Read in all of the lines containing output parameters and store them in OutList(:)
               ! customm implementation to avoid need for "END" keyword line
               DO
                  ! read a line
                  Line = NextLine(i)
                  Line = adjustl(trim(Line))   ! remove leading whitespace

                  CALL Conv2UC(Line)   ! convert to uppercase for easy string matching

                  if ((INDEX(Line, "---") > 0) .or. (INDEX(Line, "END") > 0)) EXIT ! stop if we hit a header line or the keyword "END"
                     
                  ! Check if we have a quoted string at the beginning.  Ignore anything outside the quotes if so (this is the ReadVar behaviour for quoted strings).
                  IF (SCAN(Line(1:1), '''"' ) == 1_IntKi ) THEN
                     QuoteCh = SCAN( Line(2:), '''"' )            ! last quote
                     IF (QuoteCh < 1)  QuoteCh = LEN_TRIM(Line)   ! in case no end quote
                     Line(QuoteCh+2:) = ' '    ! blank out everything after last quote
                  END IF

                  NumWords = CountWords( Line )    ! The number of words in Line.

                  p%NumOuts = p%NumOuts + NumWords  ! The total number of output channels read in so far.


                  IF ( p%NumOuts > MaxAryLen )  THEN  ! Check to see if the maximum # allowable in the array has been reached.

                     ErrStat = ErrID_Fatal
                     ErrMsg = 'Error while reading output channels: The maximum number of output channels allowed is '//TRIM( Int2LStr(MaxAryLen) )//'.'
                     EXIT

                  ELSE
                     CALL GetWords ( Line, OutList((p%NumOuts - NumWords + 1):p%NumOuts), NumWords )

                  END IF

               END DO

               ! process the OutList array and set up the index arrays for the requested output quantities
               CALL MDIO_ProcessOutList(OutList, p, m, y, InitOut, ErrStat2, ErrMsg2 )
                  CALL CheckError( ErrStat2, ErrMsg2 )
                  IF (ErrStat >= AbortErrLev) RETURN


            !-------------------------------------------------------------------------------------------
            else  ! otherwise ignore this line that isn't a recognized header line and read the next line
               Line = NextLine(i)
            end if

            !-------------------------------------------------------------------------------------------
         
         else ! otherwise ignore this line, which doesn't have the "---" or header line and read the next line
            Line = NextLine(i)
         end if
     
      end do


      ! this is the end of parsing the input file, so cleanup anything we don't need anymore
      CALL CleanUp()
      
      ! End of input file parsing from the FileInfo_In data structure
      !-------------------------------------------------------------------------------------------------
      
      
      
      
      
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN



      !-------------------------------------------------------------------------------------------------
      !          Connect mooring system together and make necessary allocations
      !-------------------------------------------------------------------------------------------------

      CALL WrNr('   Created mooring system: ' )

!     p%NAnchs = 0   ! this is the number of "fixed" type Connections. <<<<<<<<<<<<<<

      CALL WrScr(trim(Num2LStr(p%nLines))//' lines, '//trim(Num2LStr(p%NConnects))//' points, '//trim(Num2LStr(p%nRods))//' rods, '//trim(Num2LStr(p%nBodies))//' bodies.')




 !     ! now go back through and record the fairlead Id numbers (this >>>WAS<<< all the "connecting" that's required) <<<<
 !     J = 1  ! counter for fairlead number
 !     K = 1  ! counter for connect number
 !     DO I = 1,p%NConnects
 !        IF (m%ConnectList(I)%typeNum == 1) THEN
 !          m%CpldConIs(J) = I             ! if a vessel connection, add ID to list
 !          J = J + 1
 !        ELSE IF (m%ConnectList(I)%typeNum == 2) THEN
 !          m%FreeConIs(K) = I             ! if a connect connection, add ID to list
 !          K = K + 1
 !        END IF
 !     END DO

   IF (wordy > 1) print *, "nLineTypes     = ",p%nLineTypes    
   IF (wordy > 1) print *, "nRodTypes      = ",p%nRodTypes     
   IF (wordy > 1) print *, "nConnects      = ",p%nConnects     
   IF (wordy > 1) print *, "nConnectsExtra = ",p%nConnectsExtra
   IF (wordy > 1) print *, "nBodies        = ",p%nBodies       
   IF (wordy > 1) print *, "nRods          = ",p%nRods         
   IF (wordy > 1) print *, "nLines         = ",p%nLines        
   IF (wordy > 1) print *, "nCtrlChans     = ",p%nCtrlChans        
   IF (wordy > 1) print *, "nFails         = ",p%nFails        
   IF (wordy > 1) print *, "nFreeBodies    = ",p%nFreeBodies   
   IF (wordy > 1) print *, "nFreeRods      = ",p%nFreeRods     
   IF (wordy > 1) print *, "nFreeCons      = ",p%nFreeCons     
   IF (wordy > 1) print *, "nCpldBodies    = ",p%nCpldBodies   
   IF (wordy > 1) print *, "nCpldRods      = ",p%nCpldRods     
   IF (wordy > 1) print *, "nCpldCons      = ",p%nCpldCons     
   IF (wordy > 1) print *, "NConns         = ",p%NConns        
   IF (wordy > 1) print *, "NAnchs         = ",p%NAnchs              
      
   IF (wordy > 2) print *, "FreeConIs are ", m%FreeConIs
   IF (wordy > 2) print *, "CpldConIs are ", m%CpldConIs


   ! write system description to log file
   if (p%writeLog > 1) then
      write(p%UnLog, '(A)') "----- MoorDyn Model Summary (to be written) -----"
   end if



      !------------------------------------------------------------------------------------
      !                          fill in state vector index record holders
      !------------------------------------------------------------------------------------

      ! allocate state vector index record holders...

      

 !    ! allocate list of starting and ending state vector indices for each free connection
 !    ALLOCATE ( m%ConStateIs1(p%nFreeCons), m%ConStateIsN(p%nFreeCons), STAT = ErrStat )
 !    IF ( ErrStat /= ErrID_None ) THEN
 !      CALL CheckError(ErrID_Fatal, ' Error allocating ConStateIs array.')
 !      RETURN
 !    END IF
 !    
 !    ! allocate list of starting and ending state vector indices for each line  - does this belong elsewhere?
 !    ALLOCATE ( m%LineStateIs1(p%nLines), m%LineStateIsN(p%nLines), STAT = ErrStat )
 !    IF ( ErrStat /= ErrID_None ) THEN
 !      CALL CheckError(ErrID_Fatal, ' Error allocating LineStateIs arrays.')
 !      RETURN
 !    END IF
 !          
 !    
 !    ! fill in values for state vector index record holders...
 !    
 !    J=0 ! start off index counter at zero
 !    
 !    ! Free Bodies...
 !    ! Free Rods...
 !    
 !    ! Free Connections...
 !    DO l = 1, p%nFreeCons
 !       J = J + 1            ! assign start index
 !       m%ConStateIs1(l) = J
 !       
 !       J = J + 5            ! assign end index (5 entries further, since nodes have 2*3 states)
 !       m%ConStateIsN(l) = J
 !    END DO
 !    
 !    ! Lines
 !    DO l = 1, p%nLines
 !       J = J + 1            ! assign start index
 !       m%LineStateIs1(l) = J
 !       
 !       J = J + 6*(m%LineList(l)%N - 1) - 1  ! !add 6 state variables for each internal node
 !       m%LineStateIsN(l) = J
 !    END DO
 !
 !
 !    ! record number of states
 !    m%Nx = J
 
 
      !------------------------------------------------------------------------------------
      !                               prepare state vector etc.
      !------------------------------------------------------------------------------------

      ! the number of states is Nx
      m%Nx = Nx
      
      IF (wordy > 0) print *, "allocating state vectors to size ", Nx

      ! allocate state vector and temporary state vectors based on size just calculated
      ALLOCATE ( x%states(m%Nx), m%xTemp%states(m%Nx), m%xdTemp%states(m%Nx), STAT = ErrStat2 )
      IF ( ErrStat2 /= ErrID_None ) THEN
        ErrMsg  = ' Error allocating state vectors.'
        !CALL CleanUp()
        RETURN
      END IF
      x%states        = 0.0_DbKi
      m%xTemp%states  = 0.0_DbKi
      m%xdTemp%states = 0.0_DbKi



      ! ================================ initialize system ================================
      ! This will also set the initial positions of any dependent (child) objects

      ! call ground body to update all the fixed things...
      m%GroundBody%r6(4:6) = 0.0_DbKi
      CALL Body_SetDependentKin(m%GroundBody, 0.0_DbKi, m)

    ! m%GroundBody%OrMat = EulerConstruct( m%GroundBody%r6(4:6) ) ! make sure it's OrMat is set up  <<< need to check this approach
      
   !   ! first set/update the kinematics of all the fixed things (>>>> eventually do this by using a ground body <<<<)
   !   ! only doing connections so far
   !   DO J = 1,p%nConnects 
   !      if (m%ConnectList(J)%typeNum == 1) then
   !         ! set the attached line endpoint positions:
   !         CALL Connect_SetKinematics(m%ConnectList(J), m%ConnectList(J)%r, (/0.0_DbKi,0.0_DbKi,0.0_DbKi/), 0.0_DbKi, m%LineList) 
   !      end if
   !   END DO 
      
      
      ! Initialize coupled objects based on passed kinematics
      ! (set up initial condition of each coupled object based on values specified by glue code) 
      ! Also create i/o meshes       
      
      ALLOCATE ( u%CoupledKinematics(p%nTurbines), STAT = ErrStat2 )
      IF ( ErrStat2 /= ErrID_None ) THEN
        CALL CheckError(ErrID_Fatal, ' Error allocating CoupledKinematics input array.')
        RETURN
      END IF
      ALLOCATE ( y%CoupledLoads(p%nTurbines), STAT = ErrStat2 )
      IF ( ErrStat2 /= ErrID_None ) THEN
        CALL CheckError(ErrID_Fatal, ' Error allocating CoupledLoads output array.')
        RETURN
      END IF
      
      ! Go through each turbine and set up its mesh and initial positions of coupled objects
      DO iTurb = 1,p%nTurbines

         ! calculate rotation matrix OrMat for the initial orientation provided for this turbine
         CALL SmllRotTrans('PtfmInit', InitInp%PtfmInit(4,iTurb),InitInp%PtfmInit(5,iTurb),InitInp%PtfmInit(6,iTurb), OrMat, '', ErrStat2, ErrMsg2)
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN
         
         ! count number of coupling nodes needed for the mesh of this turbine
         K = p%nCpldBodies(iTurb) + p%nCpldRods(iTurb) + p%nCpldCons(iTurb)
         if (K == 0) K = 1         ! Always have at least one node (it will be a dummy node if no fairleads are attached)

         ! create input mesh for fairlead kinematics
         CALL MeshCreate(BlankMesh=u%CoupledKinematics(iTurb) , &
                       IOS= COMPONENT_INPUT, Nnodes = K, &
                       TranslationDisp=.TRUE., TranslationVel=.TRUE., &
                       Orientation=.TRUE., RotationVel=.TRUE., &
                       TranslationAcc=.TRUE., RotationAcc= .TRUE., &
                       ErrStat=ErrStat2, ErrMess=ErrMsg2)

         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN
      
         ! note: in MoorDyn-F v2, the points in the mesh correspond in order to all the coupled bodies, then rods, then connections
         ! >>> make sure all coupled objects have been offset correctly by the PtfmInit values, including if it's a farm situation -- below or where the objects are first created <<<<
         
         
         J = 0 ! this is the counter through the mesh points for each turbine
         
         DO l = 1,p%nCpldBodies(iTurb)
            J = J + 1
         
            rRef = m%BodyList(m%CpldBodyIs(l,iTurb))%r6              ! for now set reference position as per input file <<< 
            !OrMatRef = 
            
            CALL MeshPositionNode(u%CoupledKinematics(iTurb), J, rRef(1:3), ErrStat2, ErrMsg2) ! defaults to identity orientation matrix
            !TODO: >>> should also maybe set reference orientation (which might make part of a couple lines down redundant) <<<
            
            ! calculate initial point relative position, adjusted due to initial platform translations
            u%CoupledKinematics(iTurb)%TranslationDisp(:,J) = InitInp%PtfmInit(1:3,iTurb) - rRef(1:3)
            
            OrMat2 = MATMUL(OrMat, TRANSPOSE( EulerConstruct( rRef(4:6))))  ! combine the Body's relative orientation with the turbine's initial orientation
            u%CoupledKinematics(iTurb)%Orientation(:,:,J) = OrMat2          ! set the result as the current orientation of the body <<<

            ! set absolute initial positions in MoorDyn
            m%BodyList(m%CpldBodyIs(l,iTurb))%r6(1:3) = u%CoupledKinematics(iTurb)%Position(:,J) + u%CoupledKinematics(iTurb)%TranslationDisp(:,J) + p%TurbineRefPos(:,iTurb)
            m%BodyList(m%CpldBodyIs(l,iTurb))%r6(4:6) = EulerExtract(MATMUL(OrMat, TRANSPOSE( EulerConstruct( rRef(4:6)))))     ! apply rotation from PtfmInit onto input file's body orientation to get its true initial orientation

            CALL MeshConstructElement(u%CoupledKinematics(iTurb), ELEMENT_POINT, ErrStat2, ErrMsg2, J)      ! set node as point element
            
            ! lastly, do this to initialize any attached Rods or Points and set their positions
            CALL Body_InitializeUnfree( m%BodyList(m%CpldBodyIs(l,iTurb)), m )

         END DO 
         
         DO l = 1,p%nCpldRods(iTurb)   ! keeping this one simple for now, positioning at whatever is specified in input file <<<<< should change to glue code!
            J = J + 1
            
            rRef = m%RodList(m%CpldRodIs(l,iTurb))%r6          ! for now set reference position as per input file <<< 
            OrMatRef = TRANSPOSE( m%RodList(m%CpldRodIs(l,iTurb))%OrMat )  ! for now set reference orientation as per input file <<< 
            CALL MeshPositionNode(u%CoupledKinematics(iTurb), J, rRef(1:3), ErrStat2, ErrMsg2, OrMatRef)  ! assign the reference position and orientation
            
            ! calculate initial point relative position, adjusted due to initial platform rotations and translations  <<< could convert to array math
            u%CoupledKinematics(iTurb)%TranslationDisp(1,J) = InitInp%PtfmInit(1,iTurb) + OrMat(1,1)*rRef(1) + OrMat(2,1)*rRef(2) + OrMat(3,1)*rRef(3) - rRef(1)
            u%CoupledKinematics(iTurb)%TranslationDisp(2,J) = InitInp%PtfmInit(2,iTurb) + OrMat(1,2)*rRef(1) + OrMat(2,2)*rRef(2) + OrMat(3,2)*rRef(3) - rRef(2)
            u%CoupledKinematics(iTurb)%TranslationDisp(3,J) = InitInp%PtfmInit(3,iTurb) + OrMat(1,3)*rRef(1) + OrMat(2,3)*rRef(2) + OrMat(3,3)*rRef(3) - rRef(3)
            
            OrMat2 = MATMUL(OrMat, TRANSPOSE( EulerConstruct( rRef(4:6))))  ! combine the Rod's relative orientation with the turbine's initial orientation
            u%CoupledKinematics(iTurb)%Orientation(:,:,J) = OrMat2          ! set the result as the current orientation of the rod <<<
            
            ! set absolute initial positions in MoorDyn
            m%RodList(m%CpldRodIs(l,iTurb))%r6(1:3) = u%CoupledKinematics(iTurb)%Position(:,J) + u%CoupledKinematics(iTurb)%TranslationDisp(:,J) + p%TurbineRefPos(:,iTurb)
            m%RodList(m%CpldRodIs(l,iTurb))%r6(4:6) = EulerExtract(MATMUL(OrMat, OrMatRef))     ! apply rotation from PtfmInit onto input file's rod orientation to get its true initial orientation
            
            ! >>> still need to set Rod initial orientations accounting for PtfmInit rotation <<<
            
            CALL MeshConstructElement(u%CoupledKinematics(iTurb), ELEMENT_POINT, ErrStat2, ErrMsg2, J)
            
            ! lastly, do this to set the attached line endpoint positions:
            CALL Rod_SetKinematics(m%RodList(m%CpldRodIs(l,iTurb)), REAL(rRef,R8Ki), m%zeros6, m%zeros6, 0.0_DbKi, m)
         END DO 

         DO l = 1,p%nCpldCons(iTurb)   ! keeping this one simple for now, positioning at whatever is specified by glue code <<<
            J = J + 1
            
            ! set reference position as per input file  <<< what about turbine positions in array?
            rRef(1:3) = m%ConnectList(m%CpldConIs(l,iTurb))%r                           
            CALL MeshPositionNode(u%CoupledKinematics(iTurb), J, rRef(1:3), ErrStat2, ErrMsg2)  
            
            ! calculate initial point relative position, adjusted due to initial platform rotations and translations  <<< could convert to array math
            u%CoupledKinematics(iTurb)%TranslationDisp(1,J) = InitInp%PtfmInit(1,iTurb) + OrMat(1,1)*rRef(1) + OrMat(2,1)*rRef(2) + OrMat(3,1)*rRef(3) - rRef(1)
            u%CoupledKinematics(iTurb)%TranslationDisp(2,J) = InitInp%PtfmInit(2,iTurb) + OrMat(1,2)*rRef(1) + OrMat(2,2)*rRef(2) + OrMat(3,2)*rRef(3) - rRef(2)
            u%CoupledKinematics(iTurb)%TranslationDisp(3,J) = InitInp%PtfmInit(3,iTurb) + OrMat(1,3)*rRef(1) + OrMat(2,3)*rRef(2) + OrMat(3,3)*rRef(3) - rRef(3)
                 
            ! set absolute initial positions in MoorDyn
            m%ConnectList(m%CpldConIs(l,iTurb))%r = u%CoupledKinematics(iTurb)%Position(:,J) + u%CoupledKinematics(iTurb)%TranslationDisp(:,J) + p%TurbineRefPos(:,iTurb)
            
            CALL MeshConstructElement(u%CoupledKinematics(iTurb), ELEMENT_POINT, ErrStat2, ErrMsg2, J)

            ! lastly, do this to set the attached line endpoint positions:
            rRefDub = rRef(1:3)
            CALL Connect_SetKinematics(m%ConnectList(m%CpldConIs(l,iTurb)), rRefDub, m%zeros6(1:3), m%zeros6(1:3), 0.0_DbKi, m)
         END DO 

         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN
      
         ! if no coupled objects exist for this turbine, add a single dummy element to keep I/O interp/extrap routines happy
         if (J == 0) then
            rRef = 0.0_DbKi       ! position at PRP
            CALL MeshPositionNode(u%CoupledKinematics(iTurb), 1, rRef, ErrStat2, ErrMsg2)
            CALL MeshConstructElement(u%CoupledKinematics(iTurb), ELEMENT_POINT, ErrStat2, ErrMsg2, 1)
            CALL CheckError( ErrStat2, ErrMsg2 )
            IF (ErrStat >= AbortErrLev) RETURN
         end if
         
         ! set velocities/accelerations of all mesh nodes to zero
         u%CoupledKinematics(iTurb)%TranslationVel = 0.0_ReKi
         u%CoupledKinematics(iTurb)%TranslationAcc = 0.0_ReKi
         u%CoupledKinematics(iTurb)%RotationVel    = 0.0_ReKi
         u%CoupledKinematics(iTurb)%RotationAcc    = 0.0_ReKi

         CALL MeshCommit ( u%CoupledKinematics(iTurb), ErrStat2, ErrMsg )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN

         ! copy the input fairlead kinematics mesh to make the output mesh for fairlead loads, PtFairleadLoad
         CALL MeshCopy ( SrcMesh  = u%CoupledKinematics(iTurb),   DestMesh = y%CoupledLoads(iTurb), &
                         CtrlCode = MESH_SIBLING,  IOS = COMPONENT_OUTPUT, &
                         Force = .TRUE.,  Moment = .TRUE.,  ErrStat  = ErrStat2, ErrMess=ErrMsg2 )

         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN

      end do  ! iTurb
   
      ! >>>>>> ensure the output mesh includes all elements from u%(Farm)CoupledKinematics, OR make a seperate array of output meshes for each turbine <<<<<<<<<
      

      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

      
      ! ----------------------------- Arrays for active tensioning ---------------------------
      
      ! size active tensioning inputs arrays based on highest channel number read from input file for now <<<<<<<
      
      ! find the highest channel number
      N = 0
      DO I = 1, p%NLines
         IF ( m%LineList(I)%CtrlChan > N ) then
            N = m%LineList(I)%CtrlChan       
         END IF
      END DO   
      
      ! note: it would be nice to just have input arrays of the number of control channels used, rather than from 1 up to N (the highest CtrlChan)
      
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
            if (m%LineList(J)%CtrlChan > 0)  InitOut%CableCChanRqst(m%LineList(J)%CtrlChan) = .TRUE.  ! set the flag of the corresponding channel to true
         enddo
      endif
      
      
      ! >>> set up wave stuff here??? <<<
      
      
      m%WaveTi = 1   ! set initial wave grid time interpolation index to 1 to start with
      
      
   !   Frmt = '(A10,'//TRIM(Int2LStr(p%NumOuts))//'(A1,A12))'
   !
   !   WRITE(p%MDUnOut,Frmt, IOSTAT=ErrStat2)  TRIM( 'Time' ), ( p%Delim, TRIM( p%OutParam(I)%Name), I=1,p%NumOuts )
   !
   !   WRITE(p%MDUnOut,Frmt)  TRIM( '(s)' ), ( p%Delim, TRIM( p%OutParam(I)%Units ), I=1,p%NumOuts )
   !
   !
   !
   !   ! Write the output parameters to the file
   !
   !   Frmt = '(F10.4,'//TRIM(Int2LStr(p%NumOuts))//'(A1,e10.4))' 
   !
   !   WRITE(p%MDUnOut,Frmt)  Time, ( p%Delim, y%WriteOutput(I), I=1,p%NumOuts )
   
      
      
      ! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::      


      ! if any of the coupled objects need initialization steps, that should have been taken care of already <<<<
      
      
      ! initialize objects with states, writing their initial states to the master state vector (x%states)
      
      
      !TODO: apply any initial adjustment of line length from active tensioning <<<<<<<<<<<<
      ! >>> maybe this should be skipped <<<<

      
       ! Go through Bodys and write the coordinates to the state vector
      DO l = 1,p%nFreeBodies
         CALL Body_Initialize(m%BodyList(m%FreeBodyIs(l)), x%states(m%BodyStateIs1(l) : m%BodyStateIsN(l)), m)
      END DO
      
      ! Go through independent (including pinned) Rods and write the coordinates to the state vector
      DO l = 1,p%nFreeRods
         CALL Rod_Initialize(m%RodList(m%FreeRodIs(l)), x%states(m%RodStateIs1(l):m%RodStateIsN(l)), m)
      END DO

      ! Go through independent connections (Connects) and write the coordinates to the state vector and set positions of attached line ends
      DO l = 1, p%nFreeCons
         CALL Connect_Initialize(m%ConnectList(m%FreeConIs(l)), x%states(m%ConStateIs1(l) : m%conStateIsN(l)), m)
      END DO


      ! Lastly, go through lines and initialize internal node positions using quasi-static model
      DO l = 1, p%NLines

         N = m%LineList(l)%N ! for convenience

   !      ! set end node positions and velocities from connect objects
   !      m%LineList(l)%r(:,N) = m%ConnectList(m%LineList(l)%FairConnect)%r
   !      m%LineList(l)%r(:,0) = m%ConnectList(m%LineList(l)%AnchConnect)%r
   !      m%LineList(l)%rd(:,N) = (/ 0.0, 0.0, 0.0 /)  ! set anchor end velocities to zero
   !      m%LineList(l)%rd(:,0) = (/ 0.0, 0.0, 0.0 /)  ! set fairlead end velocities to zero

         ! set initial line internal node positions using quasi-static model or straight-line interpolation from anchor to fairlead
         CALL Line_Initialize( m%LineList(l), m%LineTypeList(m%LineList(l)%PropsIdNum), p%rhoW ,  ErrStat2, ErrMsg2)
            CALL CheckError( ErrStat2, ErrMsg2 )
            IF (ErrStat >= AbortErrLev) RETURN
            !IF (ErrStat >= ErrId_Warn) CALL WrScr("   Note: Catenary pre-solver was unsuccessful for one or more lines so started with linear node spacing instead.")  ! make this statement more accurate

         IF (wordy > 2) print *, "Line ", l, " with NumSegs =", N
         IF (wordy > 2) print *, "its states range from index ", m%LineStateIs1(l), " to ", m%LineStateIsN(l)

         ! assign the resulting internal node positions to the integrator initial state vector! (velocities leave at 0)
         DO I = 1, N-1
!            print *, "I=", I
            DO J = 1, 3
!               print*, J, " ... writing position state to index ", 1*(m%LineStateIs1(l) + 3*N-3 + 3*I-3 + J-1)
               x%states(m%LineStateIs1(l) + 3*N-3 + 3*I-3 + J-1 ) = m%LineList(l)%r(J,I) ! assign position
               x%states(m%LineStateIs1(l)         + 3*I-3 + J-1 ) = 0.0_DbKi ! assign velocities (of zero)
            END DO
!            print *, m%LineList(l)%r(:,I)
         END DO
               
         ! if using viscoelastic model, initialize the internal states
         if (m%LineList(l)%ElasticMod == 2) then
            do I = 1,N
               x%states(m%LineStateIs1(l) + 6*N-6 + I-1) = m%LineList(l)%dl_1(I)   ! should be zero
            end do
         end if
         

      END DO    !l = 1, p%NLines



      ! --------------------------------------------------------------------
      !          open output file(s) and write header lines
      CALL MDIO_OpenOutput( MD_ProgDesc, p, m, InitOut, ErrStat2, ErrMsg2 )
         CALL CheckError( ErrStat2, ErrMsg2 )
         IF (ErrStat >= AbortErrLev) RETURN
      ! --------------------------------------------------------------------


      IF (wordy > 2) THEN      
         print *,"Done setup of the system (before any dynamic relaxation. State vector is as follows:"
         
         DO I = 1, m%Nx
            print *, x%states(I)
         END DO
      END IF

!      ! try writing output for troubleshooting purposes (TEMPORARY)
!      CALL MDIO_WriteOutputs(-1.0_DbKi, p, m, y, ErrStat, ErrMsg)
!      IF ( ErrStat >= AbortErrLev ) THEN
!         ErrMsg = ' Error in MDIO_WriteOutputs: '//TRIM(ErrMsg)
!         RETURN
!      END IF
!      END DO

      ! -------------------------------------------------------------------
      !        if log file, compute and write some object properties
      ! -------------------------------------------------------------------
      if (p%writeLog > 1) then
      
         write(p%UnLog, '(A)'  ) "  Bodies:"         
         DO l = 1,p%nBodies
            write(p%UnLog, '(A)'  )         "    Body"//trim(num2lstr(l))//":"            
            write(p%UnLog, '(A12, f12.4)')  "      mass: ", m%BodyList(l)%M(1,1)
         END DO
         
         write(p%UnLog, '(A)'  ) "  Rods:"
         DO l = 1,p%nRods
            write(p%UnLog, '(A)'  )         "    Rod"//trim(num2lstr(l))//":"  
            ! m%RodList(l) 
         END DO
         
         write(p%UnLog, '(A)'  ) "  Points:"
         DO l = 1,p%nFreeCons
            write(p%UnLog, '(A)'  )         "    Point"//trim(num2lstr(l))//":"  
            ! m%ConnectList(l)
         END DO
         
         write(p%UnLog, '(A)'  ) "  Lines:"
         DO l = 1,p%nLines
            write(p%UnLog, '(A)'  )         "    Line"//trim(num2lstr(l))//":"  
            ! m%LineList(l)
         END DO
      
      end if


      ! --------------------------------------------------------------------
      !           do dynamic relaxation to get ICs
      ! --------------------------------------------------------------------

      ! only do this if TMaxIC > 0
      if (InputFileDat%TMaxIC > 0.0_DbKi) then

         CALL WrScr("   Finalizing initial conditions using dynamic relaxation."//NewLine)  ! newline because next line writes over itself

         ! boost drag coefficient of each line type  <<<<<<<< does this actually do anything or do lines hold these coefficients???
         DO I = 1, p%nLineTypes
            m%LineTypeList(I)%Cdn = m%LineTypeList(I)%Cdn * InputFileDat%CdScaleIC
            m%LineTypeList(I)%Cdt = m%LineTypeList(I)%Cdt * InputFileDat%CdScaleIC   ! <<<<< need to update this to apply to all objects' drag
         END DO

         ! allocate array holding 10 latest fairlead tensions
         ALLOCATE ( FairTensIC(p%nLines, 10), STAT = ErrStat2 )
         IF ( ErrStat2 /= ErrID_None ) THEN
            CALL CheckError( ErrID_Fatal, ErrMsg2 )
            RETURN
         END IF

         ! initialize fairlead tension memory at changing values so things start unconverged
         DO J = 1,p%nLines
            DO I = 1, 10
               FairTensIC(J,I) = I
            END DO
         END DO


         ! round dt to integer number of time steps  
         NdtM = ceiling(InputFileDat%dtIC/p%dtM0)            ! get number of mooring time steps to do based on desired time step size
         dtM = InputFileDat%dtIC/real(NdtM, DbKi)            ! adjust desired time step to satisfy dt with an integer number of time steps

         t = 0.0_DbKi     ! start time at zero

         ! because TimeStep wants an array...
         call MD_CopyInput( u, u_array(1), MESH_NEWCOPY, ErrStat2, ErrMsg2 )  ! make a size=1 array of inputs (since MD_RK2 expects an array to InterpExtrap)
         call MD_CopyInput( u,  u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2 )  ! also make an inputs object to interpExtrap to
         t_array(1) = t                                                       ! fill in the times "array" for u_array

         DO I = 1, ceiling(InputFileDat%TMaxIC/InputFileDat%dtIC)   ! loop through IC gen time steps, up to maximum


            !loop through line integration time steps
            DO J = 1, NdtM                                 ! for (double ts=t; ts<=t+ICdt-dts; ts+=dts)

               CALL MD_RK2(t, dtM, u_interp, u_array, t_array, p, x, xd, z, other, m, ErrStat2, ErrMsg2)
                              
               ! check for NaNs - is this a good place/way to do it?
               DO K = 1, m%Nx
                  IF (Is_NaN(x%states(K))) THEN
                     ErrStat = ErrID_Fatal
                     ErrMsg = ' NaN state detected.'
                     EXIT
                  END IF
               END DO
               
               IF (ErrStat == ErrID_Fatal) THEN
                  CALL WrScr("NaN detected at time "//TRIM(Num2LStr(t))//" during MoorDyn's dynamic relaxation process.")
                  IF (wordy > 1) THEN
                     print *, "Here is the state vector: "
                     print *, x%states
                  END IF
                  EXIT
               END IF

            END DO  ! J  time steps

         !   ! integrate the EOMs one DTIC s time step
         !   CALL TimeStep ( t, InputFileDat%dtIC, u_array, t_array, p, x, xd, z, other, m, ErrStat, ErrMsg )
         !      CALL CheckError( ErrStat2, ErrMsg2 )
         !      IF (ErrStat >= AbortErrLev) RETURN

            ! store new fairlead tension (and previous fairlead tensions for comparison)
            DO l = 1, p%nLines
            
               DO K=0,8   ! we want to count down from 10 to 2 .
                  FairTensIC(l, 10-K) = FairTensIC(l, 9-K)   ! this pushes stored values up in the array
               END DO
                  
               ! now store latest value of each line's fairlead (end B) tension   
               FairTensIC(l,1) = TwoNorm(m%LineList(l)%Fnet(:, m%LineList(l)%N))
            END DO


            ! provide status message
            ! bjj: putting this in a string so we get blanks to cover up previous values (if current string is shorter than previous one)
            Message = '   t='//trim(Num2LStr(t))//'  FairTen 1: '//trim(Num2LStr(FairTensIC(1,1)))// &
                           ', '//trim(Num2LStr(FairTensIC(1,2)))//', '//trim(Num2LStr(FairTensIC(1,3))) 
            CALL WrOver( Message )

            ! check for convergence (compare current tension at each fairlead with previous 9 values)
            IF (I > 9) THEN
               
               Converged = 1
               
               ! check for non-convergence
               
               DO l = 1, p%nLines   
                  DO K = 1,9
                     IF ( abs( FairTensIC(l,K)/FairTensIC(l,K+1) - 1.0 ) > InputFileDat%threshIC ) THEN
                        Converged = 0
                        EXIT
                     END IF
                  END DO
                  
                  IF (Converged == 0) EXIT   ! make sure we exit this loop too
               END DO

               IF (Converged == 1)  THEN  ! if we made it with all cases satisfying the threshold
                  CALL WrScr('   Fairlead tensions converged to '//trim(Num2LStr(100.0*InputFileDat%threshIC))//'% after '//trim(Num2LStr(t))//' seconds.')
                  EXIT  ! break out of the time stepping loop
               END IF
            END IF

            IF (I == ceiling(InputFileDat%TMaxIC/InputFileDat%dtIC) ) THEN
               CALL WrScr('   Fairlead tensions did not converge within TMaxIC='//trim(Num2LStr(InputFileDat%TMaxIC))//' seconds.')
               !ErrStat = ErrID_Warn
               !ErrMsg = '  MD_Init: ran dynamic convergence to TMaxIC without convergence'
            END IF

         END DO ! I ... looping through time steps



         CALL MD_DestroyInput( u_array(1), ErrStat2, ErrMsg2 )

         ! UNboost drag coefficient of each line type   <<<
         DO I = 1, p%nLineTypes
            m%LineTypeList(I)%Cdn = m%LineTypeList(I)%Cdn / InputFileDat%CdScaleIC
            m%LineTypeList(I)%Cdt = m%LineTypeList(I)%Cdt / InputFileDat%CdScaleIC
         END DO

      end if ! InputFileDat%TMaxIC > 0
      

      p%dtCoupling = DTcoupling  ! store coupling time step for use in updatestates

      other%dummy = 0
      xd%dummy    = 0
      z%dummy     = 0      
      
      if (InitInp%Linearize) then
         call MD_Init_Jacobian(InitInp, p, u, y, m, InitOut, ErrStat2, ErrMsg2); if(Failed()) return
      endif
      
      CALL WrScr('   MoorDyn initialization completed.')
      
      m%LastOutTime = -1.0_DbKi    ! set to nonzero to ensure that output happens at the start of simulation at t=0
      
      ! TODO: add feature for automatic water depth increase based on max anchor depth!

   CONTAINS


      LOGICAL FUNCTION AllocateFailed(arrayName)

         CHARACTER(*), INTENT(IN   )      :: arrayName     ! The array name
         
         call SetErrStat(ErrStat2, "Error allocating space for "//trim(arrayName)//" array.", ErrStat, ErrMsg, 'MD_Init') 
         AllocateFailed = ErrStat2 >= AbortErrLev
         if (AllocateFailed) call CleanUp() !<<<<<<<<<< need to fix this up
      END FUNCTION AllocateFailed
      
      
      LOGICAL FUNCTION Failed()

         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MD_Init') 
         Failed = ErrStat >= AbortErrLev
         if (Failed) call CleanUp()
      END FUNCTION Failed


      SUBROUTINE CheckError(ErrID,Msg)
         ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev

            ! Passed arguments
         INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
         CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

         INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
         CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)

         ! Set error status/message;
         IF ( ErrID /= ErrID_None ) THEN

            IF (ErrStat /= ErrID_None) ErrMsg = TRIM(ErrMsg)//NewLine   ! if there's a pre-existing warning/error, retain the message and start a new line

            ErrMsg = TRIM(ErrMsg)//' MD_Init:'//TRIM(Msg)
            ErrStat = MAX(ErrStat, ErrID)

            ! Clean up if we're going to return on error: close files, deallocate local arrays


            IF ( ErrStat >= AbortErrLev ) THEN                
               IF (ALLOCATED(m%CpldConIs        ))  DEALLOCATE(m%CpldConIs       )
               IF (ALLOCATED(m%FreeConIs       ))  DEALLOCATE(m%FreeConIs       )
               IF (ALLOCATED(m%LineStateIs1     ))  DEALLOCATE(m%LineStateIs1     )
               IF (ALLOCATED(m%LineStateIsN     ))  DEALLOCATE(m%LineStateIsN     )
               IF (ALLOCATED(m%ConStateIs1      ))  DEALLOCATE(m%ConStateIs1     )
               IF (ALLOCATED(m%ConStateIsN      ))  DEALLOCATE(m%ConStateIsN     )
               IF (ALLOCATED(x%states           ))  DEALLOCATE(x%states           )
               IF (ALLOCATED(FairTensIC         ))  DEALLOCATE(FairTensIC         )

               call CleanUp()    ! make sure to close files 
            END IF
         END IF

      END SUBROUTINE CheckError

      SUBROUTINE CleanUp()
        ! ErrStat = ErrID_Fatal  
        call MD_DestroyInputFileType( InputFileDat, ErrStat2, ErrMsg2 )    ! Ignore any error messages from this
        IF (p%UnLog > 0_IntKi) CLOSE( p%UnLog )       ! Remove this when the log file is kept open during the full simulation
      END SUBROUTINE

      !> If for some reason the file is truncated, it is possible to get into an infinite loop
      !! in a while looking for the next section and accidentally overstep the end of the array
      !! resulting in a segfault.  This function will trap that issue and return a section break
      CHARACTER(1024) function NextLine(i)
         integer, intent(inout) :: i      ! Current line number corresponding to contents of NextLine
         i=i+1             ! Increment to line next line.
         if (i>FileInfo_In%NumLines) then
            NextLine="---"       ! Set as a separator so we can escape some of the while loops
         else
            NextLine=trim(FileInfo_In%Lines(i))
            !TODO: add comment character recognition here? (discard any characters past a #)
         endif
      end function NextLine

   END SUBROUTINE MD_Init
   !----------------------------------------------------------------------------------------======




   !----------------------------------------------------------------------------------------======
   SUBROUTINE MD_UpdateStates( t, n, u, t_array, p, x, xd, z, other, m, ErrStat, ErrMsg)

      REAL(DbKi)                      , INTENT(IN   ) :: t
      INTEGER(IntKi)                  , INTENT(IN   ) :: n
      TYPE(MD_InputType)              , INTENT(INOUT) :: u(:)       ! INTENT(INOUT) ! had to change this to INOUT
      REAL(DbKi)                      , INTENT(IN   ) :: t_array(:)
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
  
      TYPE(MD_InputType)                                  :: u_interp   ! interpolated instantaneous input values to be calculated for each mooring time step

      REAL(DbKi)                                      :: t2          ! copy of time variable that will get advanced by the integrator (not sure this is necessary<<<)
      REAL(DbKi)                                      :: dtM         ! actual mooring dynamics time step
      INTEGER(IntKi)                                  :: NdtM        ! number of time steps to integrate through with RK2
      INTEGER(IntKi)                                  :: I
      INTEGER(IntKi)                                  :: J

      nTime = size(u) ! the number of times of input data provided? <<<<<<< not used

      t2 = t

! >>> removing this section and putting it inside loop of TimeStep (to be done at every time step) <<<
!      ! create space for arrays/meshes in u_interp
!      CALL MD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat2, ErrMsg2)
!         CALL CheckError( ErrStat2, ErrMsg2 )
!         IF (ErrStat >= AbortErrLev) RETURN
!
!      ! interpolate input mesh to correct time
!      CALL MD_Input_ExtrapInterp(u, t_array, u_interp, t, ErrStat2, ErrMsg2)
!         CALL CheckError( ErrStat2, ErrMsg2 )
!         IF (ErrStat >= AbortErrLev) RETURN
!
!
!      ! go through fairleads and apply motions from driver
!      DO I = 1, p%nCpldCons
!         DO J = 1,3
!            m%ConnectList(m%CpldConIs(I))%r(J)  = u_interp%PtFairleadDisplacement%Position(J,I) + u_interp%PtFairleadDisplacement%TranslationDisp(J,I)
!            m%ConnectList(m%CpldConIs(I))%rd(J) = u_interp%PtFairleadDisplacement%TranslationVel(J,I)  ! is this right? <<<
!         END DO
!      END DO
!



!      ! call function that loops through mooring model time steps
!      CALL TimeStep ( t2, p%dtCoupling, u, t_array, p, x, xd, z, other, m, ErrStat2, ErrMsg2 )
!         CALL CheckError( ErrStat2, ErrMsg2 )
!         IF (ErrStat >= AbortErrLev) RETURN

         
      ! create space for arrays/meshes in u_interp   ... is it efficient to do this every time step???
      CALL MD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg)


      ! round dt to integer number of time steps   <<<< should this be calculated only once, up front?
      NdtM = ceiling(p%dtCoupling/p%dtM0)            ! get number of mooring time steps to do based on desired time step size
      dtM = p%dtCoupling/REAL(NdtM,DbKi)             ! adjust desired time step to satisfy dt with an integer number of time steps


      !loop through line integration time steps
      DO I = 1, NdtM                                 ! for (double ts=t; ts<=t+ICdt-dts; ts+=dts)

         CALL MD_RK2(t2, dtM, u_interp, u, t_array, p, x, xd, z, other, m, ErrStat2, ErrMsg2)
         
         
         ! check for NaNs - is this a good place/way to do it?
         DO J = 1, m%Nx
            IF (Is_NaN(x%states(J))) THEN
               ErrStat = ErrID_Fatal
               ErrMsg = ' NaN state detected.'
               EXIT
            END IF
         END DO
         
         IF (ErrStat == ErrID_Fatal) THEN
            CALL WrScr("NaN detected at time "//TRIM(Num2LStr(t2))//" in MoorDyn.")
            IF (wordy > 1) THEN
               print *, ". Here is the state vector: "
               print *, x%states
            END IF
            EXIT
         END IF
      
      END DO  ! I  time steps


      ! destroy dxdt and x2, and u_interp
      !CALL MD_DestroyContState( dxdt, ErrStat, ErrMsg)
      !CALL MD_DestroyContState( x2, ErrStat, ErrMsg)
      CALL MD_DestroyInput(u_interp, ErrStat, ErrMsg)
      IF ( ErrStat /= ErrID_None ) THEN
         ErrMsg  = ' Error destroying dxdt or x2.'
      END IF
      !   CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MD_UpdateStates')


      ! check for NaNs - is this a good place/way to do it?
      DO J = 1, m%Nx
         IF (Is_NaN(x%states(J))) THEN
            ErrStat = ErrID_Fatal
            ErrMsg = ' NaN state detected.'
            EXIT
         END IF
      END DO
      
      IF (ErrStat == ErrID_Fatal) THEN
         CALL WrScr("NaN detected at time "//TRIM(Num2LStr(t2))//" in MoorDyn.")
         IF (wordy > 1) THEN
            print *, ". Here is the state vector: "
            print *, x%states
         END IF
      END IF

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
   !----------------------------------------------------------------------------------------



   !----------------------------------------------------------------------------------------
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

   !   TYPE(MD_ContinuousStateType)                   :: dxdt    ! time derivatives of continuous states (initialized in CalcContStateDeriv)
      INTEGER(IntKi)                                 :: I         ! counter
      INTEGER(IntKi)                                 :: J         ! counter
      INTEGER(IntKi)                                 :: K         ! counter
      INTEGER(IntKi)                                 :: l         ! index used for objects
      INTEGER(IntKi)                                 :: iTurb   ! counter
      
      Real(DbKi)                                     :: F6net(6)  ! net force and moment calculated on coupled objects
    
      INTEGER(IntKi)                                 :: ErrStat2  ! Error status of the operation
      CHARACTER(ErrMsgLen)                           :: ErrMsg2   ! Error message if ErrStat2 /= ErrID_None


      ! below updated to make sure outputs are current (based on provided x and u)  - similar to what's in UpdateStates

    !  ! go through fairleads and apply motions from driver
    !  DO I = 1, p%nCpldCons
    !     DO J = 1,3
    !        m%ConnectList(m%CpldConIs(I))%r(J)  = u%CoupledKinematics%Position(J,I) + u%CoupledKinematics%TranslationDisp(J,I)
    !        m%ConnectList(m%CpldConIs(I))%rd(J) = u%CoupledKinematics%TranslationVel(J,I)  ! is this right? <<<
    !     END DO
    !  END DO
      
      
    ! ! go through nodes and apply wave kinematics from driver (if water kinematics were passed in at each node in future)
    ! IF (p%WaterKin > 0) THEN
    ! 
    !    J=0
    !    ! Body reference point coordinates
    !    DO I = 1, p%nBodies
    !       J = J + 1                     
    !       m%BodyList(I)%U    = u%U(:,J)
    !       m%BodyList(I)%Ud   = u%Ud(:,J)
    !       m%BodyList(I)%zeta = u%zeta(J)
    !    END DO
    !    ! Rod node coordinates
    !    DO I = 1, p%nRods
    !       DO K = 0,m%RodList(I)%N  
    !          J = J + 1             
    !          m%RodList(I)%U (:,K) = u%U(:,J)
    !          m%RodList(I)%Ud(:,K) = u%Ud(:,J)
    !          m%RodList(I)%zeta(K) = u%zeta(J)
    !          m%RodList(I)%PDyn(K) = u%PDyn(J)
    !       END DO
    !    END DO
    !    ! Point reference point coordinates
    !    DO I = 1, p%nConnects
    !       J = J + 1
    !       m%ConnectList(I)%U    = u%U(:,J)
    !       m%ConnectList(I)%Ud   = u%Ud(:,J)
    !       m%ConnectList(I)%zeta = u%zeta(J)
    !    END DO      
    !    ! Line internal node coordinates
    !    DO I = 1, p%nLines
    !       DO K = 1, m%LineList(I)%N-1
    !          J = J + 1               
    !          m%LineList(I)%U (:,K) = u%U(:,J)
    !          m%LineList(I)%Ud(:,K) = u%Ud(:,J)
    !          m%LineList(I)%zeta(K) = u%zeta(J)
    !       END DO
    !    END DO   
    ! 
    ! END IF
      
      

      ! call CalcContStateDeriv in order to run model and calculate dynamics with provided x and u
      CALL MD_CalcContStateDeriv( t, u, p, x, xd, z, other, m, m%xdTemp, ErrStat, ErrMsg )

    !  ! assign net force on fairlead Connects to the fairlead force output mesh
    !  DO i = 1, p%nCpldCons
    !     DO J=1,3
    !        y%PtFairleadLoad%Force(J,I) = m%ConnectList(m%CpldConIs(I))%Fnet(J)
    !     END DO
    !  END DO
      
      ! now that forces have been updated, write them to the output mesh
      
      do iTurb = 1,p%nTurbines
      
         J = 0    ! mesh index
         DO l = 1,p%nCpldBodies(iTurb)
            J = J + 1
            CALL Body_GetCoupledForce(m%BodyList(m%CpldBodyIs(l,iTurb)), F6net, m, p)
            y%CoupledLoads(iTurb)%Force( :,J) = F6net(1:3)
            y%CoupledLoads(iTurb)%Moment(:,J) = F6net(4:6)
         END DO
               
         DO l = 1,p%nCpldRods(iTurb)
            J = J + 1
            CALL Rod_GetCoupledForce(m%RodList(m%CpldRodIs(l,iTurb)), F6net, m, p)
            y%CoupledLoads(iTurb)%Force( :,J) = F6net(1:3)
            y%CoupledLoads(iTurb)%Moment(:,J) = F6net(4:6)
         END DO
         
         DO l = 1,p%nCpldCons(iTurb)
            J = J + 1
            CALL Connect_GetCoupledForce(m%ConnectList(m%CpldConIs(l,iTurb)), F6net(1:3), m, p)
            y%CoupledLoads(iTurb)%Force(:,J) = F6net(1:3)
         END DO
         
      end do
      
   !  ! write all node positions to the node positons output array (if water kinematics were passed in at each node in future)
   !  ! go through the nodes and fill in the data (this should maybe be turned into a global function)
   !  J=0
   !  ! Body reference point coordinates
   !  DO I = 1, p%nBodies
   !     J = J + 1                     
   !     y%rAll(:,J) = m%BodyList(I)%r6(1:3)         
   !  END DO
   !  ! Rod node coordinates
   !  DO I = 1, p%nRods
   !     DO K = 0,m%RodList(I)%N  
   !        J = J + 1             
   !        y%rAll(:,J) = m%RodList(I)%r(:,K)
   !     END DO
   !  END DO
   !  ! Point reference point coordinates
   !  DO I = 1, p%nConnects
   !     J = J + 1
   !     y%rAll(:,J) = m%ConnectList(I)%r
   !  END DO      
   !  ! Line internal node coordinates
   !  DO I = 1, p%nLines
   !     DO K = 1, m%LineList(I)%N-1
   !        J = J + 1               
   !        y%rAll(:,J) = m%LineList(I)%r(:,K)
   !     END DO
   !  END DO   
   

      ! calculate outputs (y%WriteOutput) for glue code and write any m outputs to MoorDyn output files
      CALL MDIO_WriteOutputs(REAL(t,DbKi) , p, m, y, ErrStat2, ErrMsg2)
      CALL CheckError(ErrStat2, 'In MDIO_WriteOutputs: '//trim(ErrMsg2))
      IF ( ErrStat >= AbortErrLev ) RETURN


  !    ! destroy dxdt
  !    CALL MD_DestroyContState( dxdt, ErrStat2, ErrMsg2)
  !    CALL CheckError(ErrStat2, 'When destroying dxdt: '//trim(ErrMsg2))
  !    IF ( ErrStat >= AbortErrLev ) RETURN



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

      !      IF( ErrStat > ErrID_Warn ) THEN
      !          CALL MD_DestroyContState( dxdt, ErrStat2, ErrMsg2)
      !      END IF
         END IF

      END SUBROUTINE CheckError

   END SUBROUTINE MD_CalcOutput
   !----------------------------------------------------------------------------------------


   !----------------------------------------------------------------------------------------
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
      TYPE(MD_ContinuousStateType),       INTENT(INOUT)  :: dxdt    ! Continuous state derivatives at t
      INTEGER(IntKi),                     INTENT( OUT)   :: ErrStat ! Error status of the operation
      CHARACTER(*),                       INTENT( OUT)   :: ErrMsg  ! Error message if ErrStat /= ErrID_None


      INTEGER(IntKi)                                     :: L       ! index
      INTEGER(IntKi)                                     :: I       ! index
      INTEGER(IntKi)                                     :: J       ! index
      INTEGER(IntKi)                                     :: K       ! index
      INTEGER(IntKi)                                     :: iTurb   ! index
      INTEGER(IntKi)                                     :: Istart  ! start index of line/connect in state vector
      INTEGER(IntKi)                                     :: Iend    ! end index of line/connect in state vector

      REAL(DbKi)                                         :: temp(3) ! temporary for passing kinematics
      
      REAL(DbKi)                                         :: r6_in(6) ! temporary for passing kinematics
      REAL(DbKi)                                         :: v6_in(6) ! temporary for passing kinematics
      REAL(DbKi)                                         :: a6_in(6) ! temporary for passing kinematics
      REAL(DbKi)                                         :: r_in(3)  ! temporary for passing kinematics
      REAL(DbKi)                                         :: rd_in(3) ! temporary for passing kinematics
      REAL(DbKi)                                         :: a_in(3)  ! temporary for passing kinematics

      INTEGER(IntKi)                                     :: ErrStat2 ! Error status of the operation
      CHARACTER(ErrMsgLen)                               :: ErrMsg2  ! Error message if ErrStat2 /= ErrID_None
      character(*), parameter                            :: RoutineName = 'MD_CalcContStateDeriv'
      
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg = ""

      ! allocate dxdt if not already allocated (e.g. if called for linearization)
      IF (.NOT. ALLOCATED(dxdt%states) ) THEN
         CALL AllocAry( dxdt%states,  SIZE(x%states),  'dxdt%states',  ErrStat2, ErrMsg2 )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         IF ( ErrStat >= AbortErrLev ) RETURN
      END IF

      ! clear connection force and mass values <M<<<<<<<<<<<<<<<<<<<<<<<
      DO L = 1, p%NConnects
        DO J = 1,3
          m%ConnectList(L)%Fnet(J) = 0.0_DbKi
          m%ConnectList(L)%Fnet(J) = 0.0_DbKi
          DO K = 1,3
            m%ConnectList(L)%M   (K,J) = 0.0_DbKi
            m%ConnectList(L)%M   (K,J) = 0.0_DbKi
          END DO
        END DO
      END DO


      ! call ground body to update all the fixed things...
      !GroundBody->updateFairlead( t );                  <<<< manually set anchored connection stuff for now here
      r6_in = 0.0_DbKi
      v6_in = 0.0_DbKi
      CALL Body_SetKinematics(m%GroundBody, r6_in, v6_in, m%zeros6, t, m)
      
      ! ---------------------------------- coupled things ---------------------------------
      ! Apply displacement and velocity terms here. Accelerations will be considered to calculate inertial loads at the end.      
      ! Note: TurbineRefPos is to offset into farm's true global reference based on turbine X and Y reference positions (these should be 0 for regular FAST use)
         
      
      DO iTurb = 1, p%nTurbines
         
         J = 0  ! J is the index of the coupling points in the input mesh CoupledKinematics
         ! any coupled bodies (type -1)
         DO l = 1,p%nCpldBodies(iTurb)
            J = J + 1
            r6_in(1:3) = u%CoupledKinematics(iTurb)%Position(:,J) + u%CoupledKinematics(iTurb)%TranslationDisp(:,J) + p%TurbineRefPos(:,iTurb)
            !r6_in(4:6) = EulerExtract( TRANSPOSE( u%CoupledKinematics(iTurb)%Orientation(:,:,J) ) )
            r6_in(4:6) = EulerExtract( u%CoupledKinematics(iTurb)%Orientation(:,:,J) )   ! <<< changing back
            v6_in(1:3) = u%CoupledKinematics(iTurb)%TranslationVel(:,J)
            v6_in(4:6) = u%CoupledKinematics(iTurb)%RotationVel(:,J)
            a6_in(1:3) = u%CoupledKinematics(iTurb)%TranslationAcc(:,J)
            a6_in(4:6) = u%CoupledKinematics(iTurb)%RotationAcc(:,J)
         
            CALL Body_SetKinematics(m%BodyList(m%CpldBodyIs(l,iTurb)), r6_in, v6_in, a6_in, t, m)
         END DO
         
         ! any coupled rods (type -1 or -2)    note, rotations ignored if it's a pinned rod
         DO l = 1,p%nCpldRods(iTurb)
            J = J + 1

            r6_in(1:3) = u%CoupledKinematics(iTurb)%Position(:,J) + u%CoupledKinematics(iTurb)%TranslationDisp(:,J) + p%TurbineRefPos(:,iTurb)
            r6_in(4:6) = MATMUL( u%CoupledKinematics(iTurb)%Orientation(:,:,J) , (/0.0, 0.0, 1.0/) ) ! <<<< CHECK ! adjustment because rod's rotational entries are a unit vector, q
            v6_in(1:3) = u%CoupledKinematics(iTurb)%TranslationVel(:,J)
            v6_in(4:6) = u%CoupledKinematics(iTurb)%RotationVel(:,J)
            a6_in(1:3) = u%CoupledKinematics(iTurb)%TranslationAcc(:,J)
            a6_in(4:6) = u%CoupledKinematics(iTurb)%RotationAcc(:,J)
         
            CALL Rod_SetKinematics(m%RodList(m%CpldRodIs(l,iTurb)), r6_in, v6_in, a6_in, t, m)
    
         END DO
         
         ! any coupled points (type -1)
         DO l = 1, p%nCpldCons(iTurb)
            J = J + 1
            
            r_in  = u%CoupledKinematics(iTurb)%Position(:,J) + u%CoupledKinematics(iTurb)%TranslationDisp(:,J) + p%TurbineRefPos(:,iTurb)
            rd_in = u%CoupledKinematics(iTurb)%TranslationVel(:,J)
            a_in(1:3) = u%CoupledKinematics(iTurb)%TranslationAcc(:,J)
            CALL Connect_SetKinematics(m%ConnectList(m%CpldConIs(l,iTurb)), r_in, rd_in, a_in, t, m)
            
            !print "(f8.5, f12.6, f12.6, f8.4, f8.4, f8.4, f8.4)", t, r_in(1), r_in(3), rd_in(1), rd_in(3), a_in(1), a_in(3)
            
         END DO
         
      end do  ! iTurb
      
      
      ! >>>>> in theory I would repeat the above but for each turbine in the case of array use here <<<<<
      !  DO I = 1,p%nTurbines
      !     J = 0?
      !     other logic?
      !     nvm: need to get kinematics from entries in u%FarmCoupledKinematics(I)%Position etc.
      !     nvm: using knowledge of p%meshIndex or something
      !   in theory might also support individual line tensioning control commands from turbines this way too, or maybe it's supercontroller level (not a short term problem though)
      
      
      ! apply line length changes from active tensioning if applicable
      DO L = 1, p%NLines
         IF (m%LineList(L)%CtrlChan > 0) then

            ! do a bounds check to prohibit excessive segment length changes (until a method to add/remove segments is created)
            IF ( u%DeltaL(m%LineList(L)%CtrlChan) > m%LineList(L)%UnstrLen / m%LineList(L)%N ) then
                ErrStat = ErrID_Fatal
                ErrMsg  = ' Active tension command will make a segment longer than the limit of twice its original length.'
                call WrScr(trim(Num2LStr(u%DeltaL(m%LineList(L)%CtrlChan)))//" is an increase of more than "//trim(Num2LStr(m%LineList(L)%UnstrLen / m%LineList(L)%N)))
                IF (wordy > 0) print *, u%DeltaL
                IF (wordy > 0) print*, m%LineList(L)%CtrlChan
                RETURN
            END IF
            IF ( u%DeltaL(m%LineList(L)%CtrlChan) < -0.5 * m%LineList(L)%UnstrLen / m%LineList(L)%N ) then
             ErrStat = ErrID_Fatal
                ErrMsg  = ' Active tension command will make a segment shorter than the limit of half its original length.'
                call WrScr(trim(Num2LStr(u%DeltaL(m%LineList(L)%CtrlChan)))//" is a reduction of more than half of "//trim(Num2LStr(m%LineList(L)%UnstrLen / m%LineList(L)%N)))
                IF (wordy > 0) print *, u%DeltaL
                IF (wordy > 0) print*, m%LineList(L)%CtrlChan
                RETURN
            END IF                

            ! for now this approach only acts on the fairlead end segment, and assumes all segment lengths are otherwise equal size
            m%LineList(L)%l( m%LineList(L)%N) = m%LineList(L)%UnstrLen/m%LineList(L)%N + u%DeltaL(m%LineList(L)%CtrlChan)       
            m%LineList(L)%ld(m%LineList(L)%N) =                                       u%DeltaLdot(m%LineList(L)%CtrlChan)       
         END IF
      END DO      
      
      
   !  ! go through nodes and apply wave kinematics from driver (if water kinematics were passed in at each node in future)
   !  IF (p%WaterKin > 0) THEN
   !  
   !     J=0
   !     ! Body reference point coordinates
   !     DO I = 1, p%nBodies
   !        J = J + 1                     
   !        m%BodyList(I)%U    = u%U(:,J)
   !        m%BodyList(I)%Ud   = u%Ud(:,J)
   !        m%BodyList(I)%zeta = u%zeta(J)
   !     END DO
   !     ! Rod node coordinates
   !     DO I = 1, p%nRods
   !        DO K = 0,m%RodList(I)%N  
   !           J = J + 1             
   !           m%RodList(I)%U (:,K) = u%U(:,J)
   !           m%RodList(I)%Ud(:,K) = u%Ud(:,J)
   !           m%RodList(I)%zeta(K) = u%zeta(J)
   !           m%RodList(I)%PDyn(K) = u%PDyn(J)
   !        END DO
   !     END DO
   !     ! Point reference point coordinates
   !     DO I = 1, p%nConnects
   !        J = J + 1
   !        m%ConnectList(I)%U    = u%U(:,J)
   !        m%ConnectList(I)%Ud   = u%Ud(:,J)
   !        m%ConnectList(I)%zeta = u%zeta(J)
   !     END DO      
   !     ! Line internal node coordinates
   !     DO I = 1, p%nLines
   !        DO K = 1, m%LineList(I)%N-1
   !           J = J + 1               
   !           m%LineList(I)%U (:,K) = u%U(:,J)
   !           m%LineList(I)%Ud(:,K) = u%Ud(:,J)
   !           m%LineList(I)%zeta(K) = u%zeta(J)
   !        END DO
   !     END DO   
   !  
   !  END IF
      
      
      ! independent or semi-independent things with their own states...
      
      ! give Bodies latest state variables (kinematics will also be assigned to dependent connections and rods, and thus line ends)
      DO l = 1,p%nFreeBodies
         CALL Body_SetState(m%BodyList(m%FreeBodyIs(l)), x%states(m%BodyStateIs1(l):m%BodyStateIsN(l)), t, m)
      END DO
      
      ! give independent or pinned rods' latest state variables (kinematics will also be assigned to attached line ends)
      DO l = 1,p%nFreeRods
         CALL Rod_SetState(m%RodList(m%FreeRodIs(l)), x%states(m%RodStateIs1(l):m%RodStateIsN(l)), t, m)
      END DO
      
      ! give Connects (independent connections) latest state variable values (kinematics will also be assigned to attached line ends)
      DO l = 1,p%nFreeCons
  !       Print *, "calling SetState for free connection, con#", m%FreeConIs(l), " with state range: ", m%ConStateIs1(l), "-", m%ConStateIsN(l)
         !K=K+1
         CALL Connect_SetState(m%ConnectList(m%FreeConIs(l)), x%states(m%ConStateIs1(l):m%ConStateIsN(l)), t, m)
      END DO
      
      ! give Lines latest state variable values for internal nodes
      DO l = 1,p%nLines
         CALL Line_SetState(m%LineList(l), x%states(m%LineStateIs1(l):m%LineStateIsN(l)), t)
      END DO

      ! calculate dynamics of free objects (will also calculate forces (doRHS()) from any child/dependent objects)...
         
      ! calculate line dynamics (and calculate line forces and masses attributed to connections)
      DO l = 1,p%nLines
         CALL Line_GetStateDeriv(m%LineList(l), dxdt%states(m%LineStateIs1(l):m%LineStateIsN(l)), m, p)  !dt might also be passed for fancy friction models
      END DO
      
      ! calculate connect dynamics (including contributions from attached lines
      ! as well as hydrodynamic forces etc. on connect object itself if applicable)
      DO l = 1,p%nFreeCons
         CALL Connect_GetStateDeriv(m%ConnectList(m%FreeConIs(l)), dxdt%states(m%ConStateIs1(l):m%ConStateIsN(l)), m, p)
      END DO
      
      ! calculate dynamics of independent Rods 
      DO l = 1,p%nFreeRods
         CALL Rod_GetStateDeriv(m%RodList(m%FreeRodIs(l)), dxdt%states(m%RodStateIs1(l):m%RodStateIsN(l)), m, p)
      END DO
      
      ! calculate dynamics of Bodies
      DO l = 1,p%nFreeBodies
         CALL Body_GetStateDeriv(m%BodyList(m%FreeBodyIs(l)), dxdt%states(m%BodyStateIs1(l):m%BodyStateIsN(l)), m, p)
      END DO
      
      
      
      ! get dynamics/forces (doRHS()) of coupled objects, which weren't addressed in above calls (this includes inertial loads)
      ! note: can do this in any order since there are no dependencies among coupled objects
      
      DO iTurb = 1,p%nTurbines
         DO l = 1,p%nCpldCons(iTurb)
         
    !        >>>>>>>> here we should pass along accelerations and include inertial loads in the calculation!!! <<<??
    !               in other words are the below good enough or do I need to call _getCoupledFOrce??
         
            CALL Connect_DoRHS(m%ConnectList(m%CpldConIs(l,iTurb)), m, p)
         END DO
         
         DO l = 1,p%nCpldRods(iTurb)
            CALL Rod_DoRHS(m%RodList(m%CpldRodIs(l,iTurb)), m, p)
            ! NOTE: this won't compute net loads on Rod. Need Rod_GetNetForceAndMass for that. Change? <<<<
         END DO
         
         DO l = 1,p%nCpldBodies(iTurb)
            CALL Body_DoRHS(m%BodyList(m%CpldBodyIs(l,iTurb)), m, p)
         END DO
      end do

      ! call ground body to update all the fixed things
      CALL Body_DoRHS(m%GroundBody, m, p)
      
      
      !print *, t, m%LineList(1)%T(1,9), m%LineList(1)%T(2,9), m%LineList(1)%T(3,9), m%LineList(3)%T(1,9), m%LineList(3)%T(2,9), m%LineList(3)%T(3,9)

   
   END SUBROUTINE MD_CalcContStateDeriv
   !----------------------------------------------------------------------------------------=====


   !----------------------------------------------------------------------------------------=======
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
         
      IF (p%UnLog > 0_IntKi) CLOSE( p%UnLog )  ! close log file if it's open
         !TODO: any need to specifically deallocate things like m%xTemp%states in the above? <<<<

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
   !----------------------------------------------------------------------------------------==================


!!==========   MD_CheckError   =======     <---------------------------------------------------------------+
! SUBROUTINE MD_CheckError(InMsg,OutMsg)
!   ! Passed arguments
!!   CHARACTER(*), INTENT(IN   ) :: InMsg       ! The input string
!   CHARACTER(*), INTENT(INOUT) :: OutMsg      ! The error message (ErrMsg)!
!
 !  OutMsg = InMsg
 !  RETURN
 !END SUBROUTINE MD_CheckError                                                                  !   -------+
 !----------------------------------------------------------------------------------------==================


   ! RK2 integrater (part of what was in TimeStep)
   !--------------------------------------------------------------
   SUBROUTINE MD_RK2 ( t, dtM, u_interp, u, t_array, p, x, xd, z, other, m, ErrStat, ErrMsg )
   
      REAL(DbKi)                     , INTENT(INOUT)      :: t          ! intial time (s) for this integration step
      REAL(DbKi)                     , INTENT(IN   )      :: dtM        ! single time step  size (s) for this integration step
      TYPE( MD_InputType )           , INTENT(INOUT)      :: u_interp   ! interpolated instantaneous input values to be calculated for each mooring time step
      TYPE( MD_InputType )           , INTENT(INOUT)      :: u(:)       ! INTENT(IN   )
      REAL(DbKi)                     , INTENT(IN   )      :: t_array(:)  ! times corresponding to elements of u(:)?
      TYPE( MD_ParameterType )       , INTENT(IN   )      :: p          ! INTENT(IN   )
      TYPE( MD_ContinuousStateType ) , INTENT(INOUT)      :: x
      TYPE( MD_DiscreteStateType )   , INTENT(IN   )      :: xd         ! INTENT(IN   )
      TYPE( MD_ConstraintStateType ) , INTENT(IN   )      :: z          ! INTENT(IN   )
      TYPE( MD_OtherStateType )      , INTENT(IN   )      :: other      ! INTENT(INOUT)
      TYPE(MD_MiscVarType)           , INTENT(INOUT)      :: m          ! INTENT(INOUT)
      INTEGER(IntKi)                 , INTENT(  OUT)      :: ErrStat
      CHARACTER(*)                   , INTENT(  OUT)      :: ErrMsg


      INTEGER(IntKi)                                      :: I          ! counter
      INTEGER(IntKi)                                      :: J          ! counter
         
   
      ! -------------------------------------------------------------------------------
      !       RK2 integrator written here, now calling CalcContStateDeriv
      !--------------------------------------------------------------------------------

      ! step 1

      CALL MD_Input_ExtrapInterp(u, t_array, u_interp, t          , ErrStat, ErrMsg)   ! interpolate input mesh to correct time (t)
   
      CALL MD_CalcContStateDeriv( t, u_interp, p, x, xd, z, other, m, m%xdTemp, ErrStat, ErrMsg )
      DO J = 1, m%Nx
         m%xTemp%states(J) = x%states(J) + 0.5*dtM*m%xdTemp%states(J)                                           !x1 = x0 + dt*f0/2.0;
      END DO

      ! step 2

      CALL MD_Input_ExtrapInterp(u, t_array, u_interp, t + 0.5_DbKi*dtM, ErrStat, ErrMsg)   ! interpolate input mesh to correct time (t+0.5*dtM)
         
      CALL MD_CalcContStateDeriv( (t + 0.5_DbKi*dtM), u_interp, p, m%xTemp, xd, z, other, m, m%xdTemp, ErrStat, ErrMsg )       !called with updated states x2 and time = t + dt/2.0
      DO J = 1, m%Nx
         x%states(J) = x%states(J) + dtM*m%xdTemp%states(J)
      END DO

      t = t + dtM  ! update time
      
      !TODO error check? <<<<

   END SUBROUTINE MD_RK2
   !--------------------------------------------------------------


   !----------------------------------------------------------------------------------------================
   ! this would do a full (coupling) time step and is no longer used
   SUBROUTINE TimeStep ( t, dtStep, u, t_array, p, x, xd, z, other, m, ErrStat, ErrMsg )
   
      REAL(DbKi)                     , INTENT(INOUT)      :: t
      REAL(DbKi)                     , INTENT(IN   )      :: dtStep     ! how long to advance the time for
      TYPE( MD_InputType )           , INTENT(INOUT)      :: u(:)       ! INTENT(IN   )
      REAL(DbKi)                     , INTENT(IN   )      :: t_array(:)  ! times corresponding to elements of u(:)?
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

  !    Real(DbKi)                                          :: tDbKi   ! double version because that's what MD_Input_ExtrapInterp needs.
      
      
      ! allocate space for x2
      CALL MD_CopyContState( x, x2, 0, ErrStat, ErrMsg)
         
      ! create space for arrays/meshes in u_interp   ... is it efficient to do this every time step???
      CALL MD_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg)
         

      Nx = size(x%states)   ! <<<< should this be the m%Nx parameter instead?


      ! round dt to integer number of time steps
      NdtM = ceiling(dtStep/p%dtM0)                  ! get number of mooring time steps to do based on desired time step size
      dtM = dtStep/REAL(NdtM,DbKi)                   ! adjust desired time step to satisfy dt with an integer number of time steps


      !loop through line integration time steps
      DO I = 1, NdtM                                 ! for (double ts=t; ts<=t+ICdt-dts; ts+=dts)
      
      
   !      tDbKi = t        ! get DbKi version of current time (why does ExtrapInterp except different time type than UpdateStates?)
         
      
         ! -------------------------------------------------------------------------------
         !       RK2 integrator written here, now calling CalcContStateDeriv
         !--------------------------------------------------------------------------------

         ! step 1

         CALL MD_Input_ExtrapInterp(u, t_array, u_interp, t          , ErrStat, ErrMsg)   ! interpolate input mesh to correct time (t)
      
         CALL MD_CalcContStateDeriv( t, u_interp, p, x, xd, z, other, m, dxdt, ErrStat, ErrMsg )
         DO J = 1, Nx
            x2%states(J) = x%states(J) + 0.5*dtM*dxdt%states(J)                                           !x1 = x0 + dt*f0/2.0;
         END DO

         ! step 2
   
         CALL MD_Input_ExtrapInterp(u, t_array, u_interp, t + 0.5_DbKi*dtM, ErrStat, ErrMsg)   ! interpolate input mesh to correct time (t+0.5*dtM)
            
         CALL MD_CalcContStateDeriv( (t + 0.5_DbKi*dtM), u_interp, p, x2, xd, z, other, m, dxdt, ErrStat, ErrMsg )       !called with updated states x2 and time = t + dt/2.0
         DO J = 1, Nx
            x%states(J) = x%states(J) + dtM*dxdt%states(J)
         END DO

         t = t + dtM  ! update time

         !----------------------------------------------------------------------------------

   ! >>> below should no longer be necessary thanks to using ExtrapInterp of u(:) within the mooring time stepping loop.. <<<
   !      ! update Fairlead positions by integrating velocity and last position (do this AFTER the processing of the time step rather than before)
   !      DO J = 1, p%nCpldCons
   !         DO K = 1, 3
   !          m%ConnectList(m%CpldConIs(J))%r(K) = m%ConnectList(m%CpldConIs(J))%r(K) + m%ConnectList(m%CpldConIs(J))%rd(K)*dtM
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
         IF (Is_NaN(x%states(J))) THEN
            ErrStat = ErrID_Fatal
            ErrMsg = ' NaN state detected.'
         END IF
      END DO
 

   END SUBROUTINE TimeStep
   !--------------------------------------------------------------



!--------------------------------------------------------------
!            Connection-Specific Subroutines
!--------------------------------------------------------------




!--------------------------------------------------------------
!            Rod-Specific Subroutines
!--------------------------------------------------------------








!--------------------------------------------------------------
!            Body-Specific Subroutines
!--------------------------------------------------------------



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ###### The following four routines are Jacobian routines for linearization capabilities #######
! If the module does not implement them, set ErrStat = ErrID_Fatal in SD_Init() when InitInp%Linearize is .true.
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the inputs (u). The partial derivatives dY/du, dX/du, dXd/du, and DZ/du are returned.
SUBROUTINE MD_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdu, dXdu, dXddu, dZdu)
   REAL(DbKi),                        INTENT(IN   ) :: t                  !< Time in seconds at operating point
   TYPE(MD_InputType),                INTENT(INOUT) :: u                  !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(MD_ParameterType),            INTENT(IN   ) :: p                  !< Parameters
   TYPE(MD_ContinuousStateType),      INTENT(IN   ) :: x                  !< Continuous states at operating point
   TYPE(MD_DiscreteStateType),        INTENT(IN   ) :: xd                 !< Discrete states at operating point
   TYPE(MD_ConstraintStateType),      INTENT(IN   ) :: z                  !< Constraint states at operating point
   TYPE(MD_OtherStateType),           INTENT(IN   ) :: OtherState         !< Other states at operating point
   TYPE(MD_OutputType),               INTENT(INOUT) :: y                  !< Output (change to inout if a mesh copy is required); Output fields are not used by this routine, but type is available here so that mesh parameter information (i.e., connectivity) does not have to be recalculated for dYdu.
   TYPE(MD_MiscVarType),              INTENT(INOUT) :: m                  !< Misc/optimization variables
   INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat            !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT) :: ErrMsg             !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dYdu(:,:)          !< Partial derivatives of output functions (Y) wrt the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXdu(:,:)          !< Partial derivatives of continuous state functions (X) wrt the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXddu(:,:)         !< Partial derivatives of discrete state functions (Xd) wrt the inputs (u) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dZdu(:,:)          !< Partial derivatives of constraint state functions (Z) wrt the inputs (u) [intent in to avoid deallocation]
   
   ! local variables
   TYPE(MD_OutputType)          :: y_m, y_p
   TYPE(MD_ContinuousStateType) :: x_m, x_p
   TYPE(MD_InputType)           :: u_perturb
   REAL(R8Ki)                   :: delta_p, delta_m   ! delta change in input (plus, minus)
   INTEGER(IntKi)               :: i
   integer(intKi)               :: ErrStat2
   character(ErrMsgLen)         :: ErrMsg2
   character(*), parameter      :: RoutineName = 'MD_JacobianPInput'
   
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''
   
   ! get OP values here:
   call MD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 ); if(Failed()) return
   
   ! make a copy of the inputs to perturb
   call MD_CopyInput( u, u_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
   
   IF ( PRESENT( dYdu ) ) THEN
      ! Calculate the partial derivative of the output functions (Y) with respect to the inputs (u) here:
      if (.not. allocated(dYdu) ) then
         call AllocAry(dYdu, p%Jac_ny, size(p%Jac_u_indx,1),'dYdu', ErrStat2, ErrMsg2); if(Failed()) return
      end if
      ! make a copy of outputs because we will need two for the central difference computations (with orientations)
      call MD_CopyOutput( y, y_p, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
      call MD_CopyOutput( y, y_m, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
      do i=1,size(p%Jac_u_indx,1)
         ! get u_op + delta_p u
         call MD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         call MD_Perturb_u( p, i, 1, u_perturb, delta_p )
         ! compute y at u_op + delta_p u
         call MD_CalcOutput( t, u_perturb, p, x, xd, z, OtherState, y_p, m, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         ! get u_op - delta_m u
         call MD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         call MD_Perturb_u( p, i, -1, u_perturb, delta_m )
         ! compute y at u_op - delta_m u
         call MD_CalcOutput( t, u_perturb, p, x, xd, z, OtherState, y_m, m, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         ! get central difference:
         call MD_Compute_dY( p, y_p, y_m, delta_p, dYdu(:,i) )
      end do
      if(Failed()) return
   END IF
   IF ( PRESENT( dXdu ) ) THEN
      if (.not. allocated(dXdu)) then
         call AllocAry(dXdu, p%Jac_nx, size(p%Jac_u_indx,1), 'dXdu', ErrStat2, ErrMsg2); if (Failed()) return
      endif
      do i=1,size(p%Jac_u_indx,1)
         ! get u_op + delta u
         call MD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         call MD_Perturb_u( p, i, 1, u_perturb, delta_p )
         ! compute x at u_op + delta u
         call MD_CalcContStateDeriv( t, u_perturb, p, x, xd, z, OtherState, m, x_p, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         ! get u_op - delta u
         call MD_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
         call MD_Perturb_u( p, i, -1, u_perturb, delta_m )
         ! compute x at u_op - delta u
         call MD_CalcContStateDeriv( t, u_perturb, p, x, xd, z, OtherState, m, x_m, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
         ! get central difference:
         ! we may have had an error allocating memory, so we'll check
         if(Failed()) return
         ! get central difference (state entries are mapped the the dXdu column in routine):
         call MD_Compute_dX( p, x_p, x_m, delta_p, dXdu(:,i) )
      end do
   END IF ! dXdu
   IF ( PRESENT( dXddu ) ) THEN
      if (allocated(dXddu)) deallocate(dXddu)
   END IF
   IF ( PRESENT( dZdu ) ) THEN
      if (allocated(dZdu)) deallocate(dZdu)
   END IF
   call CleanUp()
contains

   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed

   subroutine CleanUp()
      call MD_DestroyContState(  x_p, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call MD_DestroyContState(  x_m, ErrStat2, ErrMsg2 ) ! we don't need this any more
      call MD_DestroyOutput(     y_p, ErrStat2, ErrMsg2 )
      call MD_DestroyOutput(     y_m, ErrStat2, ErrMsg2 )
      call MD_DestroyInput(u_perturb, ErrStat2, ErrMsg2 )
   end subroutine cleanup

END SUBROUTINE MD_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the continuous states (x). The partial derivatives dY/dx, dX/dx, dXd/dx, and dZ/dx are returned.
SUBROUTINE MD_JacobianPContState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdx, dXdx, dXddx, dZdx)
   REAL(DbKi),                        INTENT(IN   ) :: t                  !< Time in seconds at operating point
   TYPE(MD_InputType),                INTENT(INOUT) :: u                  !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(MD_ParameterType),            INTENT(IN   ) :: p                  !< Parameters
   TYPE(MD_ContinuousStateType),      INTENT(IN   ) :: x                  !< Continuous states at operating point
   TYPE(MD_DiscreteStateType),        INTENT(IN   ) :: xd                 !< Discrete states at operating point
   TYPE(MD_ConstraintStateType),      INTENT(IN   ) :: z                  !< Constraint states at operating point
   TYPE(MD_OtherStateType),           INTENT(IN   ) :: OtherState         !< Other states at operating point
   TYPE(MD_OutputType),               INTENT(INOUT) :: y                  !< Output (change to inout if a mesh copy is required); Output fields are not used by this routine, but type is available here so that mesh parameter information (i.e., connectivity) does not have to be recalculated for dYdx.
   TYPE(MD_MiscVarType),              INTENT(INOUT) :: m                  !< Misc/optimization variables
   INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat            !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT) :: ErrMsg             !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dYdx(:,:)          !< Partial derivatives of output functions wrt the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXdx(:,:)          !< Partial derivatives of continuous state functions (X) wrt the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXddx(:,:)         !< Partial derivatives of discrete state functions (Xd) wrt the continuous states (x) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dZdx(:,:)          !< Partial derivatives of constraint state functions (Z) wrt the continuous states (x) [intent in to avoid deallocation]
   ! local variables
   TYPE(MD_OutputType)          :: y_p, y_m
   TYPE(MD_ContinuousStateType) :: x_p, x_m
   TYPE(MD_ContinuousStateType) :: x_perturb
   REAL(R8Ki)                   :: delta        ! delta change in input or state
   INTEGER(IntKi)               :: i, k
   INTEGER(IntKi)               :: ErrStat2
   CHARACTER(ErrMsgLen)         :: ErrMsg2
   CHARACTER(*), PARAMETER      :: RoutineName = 'MD_JacobianPContState'
   
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''
   
   ! make a copy of the continuous states to perturb NOTE: MESH_NEWCOPY
   call MD_CopyContState( x, x_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
   
   IF ( PRESENT( dYdx ) ) THEN
      ! Calculate the partial derivative of the output functions (Y) with respect to the continuous states (x) here:
      if (.not. allocated(dYdx)) then
         call AllocAry(dYdx, p%Jac_ny, p%Jac_nx, 'dYdx', ErrStat2, ErrMsg2); if(Failed()) return
      end if
      ! make a copy of outputs because we will need two for the central difference computations (with orientations)
      call MD_CopyOutput( y, y_p, MESH_NEWCOPY, ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
      call MD_CopyOutput( y, y_m, MESH_NEWCOPY, ErrStat2, ErrMsg2); if(Failed()) return
      !  Loop over the dx dimension of the dYdx array.  Perturb the corresponding state (note difference in ordering of dYdx and x%states).
      !  The p%dxIdx_map2_xStateIdx(i) is the index to the state array for the given dx index
      do i=1,p%Jac_nx      ! index into dx dimension
         ! get x_op + delta x
         call MD_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         call MD_perturb_x(p, p%dxIdx_map2_xStateIdx(i), 1, x_perturb, delta )
         ! compute y at x_op + delta x
         call MD_CalcOutput( t, u, p, x_perturb, xd, z, OtherState, y_p, m, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         ! get x_op - delta x
         call MD_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         call MD_perturb_x(p, p%dxIdx_map2_xStateIdx(i), -1, x_perturb, delta )
         ! compute y at x_op - delta x
         call MD_CalcOutput( t, u, p, x_perturb, xd, z, OtherState, y_m, m, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         ! get central difference:
         call MD_Compute_dY( p, y_p, y_m, delta, dYdx(:,i) )
      end do
      if(Failed()) return
   END IF
   
   IF ( PRESENT( dXdx ) ) THEN
      ! Calculate the partial derivative of the continuous state functions (X) with respect to the continuous states (x) here:
      if (.not. allocated(dXdx)) then
         call AllocAry(dXdx, p%Jac_nx, p%Jac_nx, 'dXdx', ErrStat2, ErrMsg2); if(Failed()) return
      end if
      !  Loop over the dx dimension of the array.  Perturb the corresponding state (note difference in ordering of dXdx and x%states).
      !  The resulting x_p and x_m are used to calculate the column for dXdx (mapping of state entry to dXdx row entry occurs in MD_Compute_dX)
      !  The p%dxIdx_map2_xStateIdx(i) is the index to the state array for the given dx index
      do i=1,p%Jac_nx      ! index into dx dimension
         ! get x_op + delta x
         call MD_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         call MD_perturb_x(p, p%dxIdx_map2_xStateIdx(i), 1, x_perturb, delta )
         ! compute x at x_op + delta x
         call MD_CalcContStateDeriv( t, u, p, x_perturb, xd, z, OtherState, m, x_p, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         ! get x_op - delta x
         call MD_CopyContState( x, x_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
         call MD_perturb_x(p, p%dxIdx_map2_xStateIdx(i), -1, x_perturb, delta )
         ! compute x at x_op - delta x
         call MD_CalcContStateDeriv( t, u, p, x_perturb, xd, z, OtherState, m, x_m, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName) 
         if(Failed()) return
         ! get central difference:
         call MD_Compute_dX( p, x_p, x_m, delta, dXdx(:,i) )
      end do
   END IF
   IF ( PRESENT( dXddx ) ) THEN
      if (allocated(dXddx)) deallocate(dXddx)
   END IF
   IF ( PRESENT( dZdx ) ) THEN
      if (allocated(dZdx)) deallocate(dZdx)
   END IF
   call CleanUp()
   
contains

   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MD_JacobianPContState') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed

   subroutine CleanUp()
      call MD_DestroyOutput(         y_p, ErrStat2, ErrMsg2 )
      call MD_DestroyOutput(         y_m, ErrStat2, ErrMsg2 )
      call MD_DestroyContState(      x_p, ErrStat2, ErrMsg2 )
      call MD_DestroyContState(      x_m, ErrStat2, ErrMsg2 )
      call MD_DestroyContState(x_perturb, ErrStat2, ErrMsg2 )
   end subroutine cleanup

END SUBROUTINE MD_JacobianPContState

!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the discrete states (xd). The partial derivatives dY/dxd, dX/dxd, dXd/dxd, and DZ/dxd are returned.
SUBROUTINE MD_JacobianPDiscState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdxd, dXdxd, dXddxd, dZdxd )
   REAL(DbKi),                        INTENT(IN   ) :: t                  !< Time in seconds at operating point
   TYPE(MD_InputType),                INTENT(INOUT) :: u                  !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(MD_ParameterType),            INTENT(IN   ) :: p                  !< Parameters
   TYPE(MD_ContinuousStateType),      INTENT(IN   ) :: x                  !< Continuous states at operating point
   TYPE(MD_DiscreteStateType),        INTENT(IN   ) :: xd                 !< Discrete states at operating point
   TYPE(MD_ConstraintStateType),      INTENT(IN   ) :: z                  !< Constraint states at operating point
   TYPE(MD_OtherStateType),           INTENT(IN   ) :: OtherState         !< Other states at operating point
   TYPE(MD_OutputType),               INTENT(INOUT) :: y                  !< Output (change to inout if a mesh copy is required); Output fields are not used by this routine, but type is available here so that mesh parameter information (i.e., connectivity) does not have to be recalculated for dYdx.
   TYPE(MD_MiscVarType),              INTENT(INOUT) :: m                  !< Misc/optimization variables
   INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat    !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dYdxd(:,:) !< Partial derivatives of output functions (Y) wrt the discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXdxd(:,:) !< Partial derivatives of continuous state functions (X) wrt the  discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXddxd(:,:)!< Partial derivatives of discrete state functions (Xd) wrt the discrete states (xd) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dZdxd(:,:) !< Partial derivatives of constraint state functions (Z) wrt discrete states (xd) [intent in to avoid deallocation]
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''
   IF ( PRESENT( dYdxd ) ) THEN
   END IF
   IF ( PRESENT( dXdxd ) ) THEN
   END IF
   IF ( PRESENT( dXddxd ) ) THEN
   END IF
   IF ( PRESENT( dZdxd ) ) THEN
   END IF
END SUBROUTINE MD_JacobianPDiscState
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) functions
!! with respect to the constraint states (z). The partial derivatives dY/dz, dX/dz, dXd/dz, and DZ/dz are returned.
SUBROUTINE MD_JacobianPConstrState( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, dYdz, dXdz, dXddz, dZdz )
   REAL(DbKi),                        INTENT(IN   ) :: t                  !< Time in seconds at operating point
   TYPE(MD_InputType),                INTENT(INOUT) :: u                  !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(MD_ParameterType),            INTENT(IN   ) :: p                  !< Parameters
   TYPE(MD_ContinuousStateType),      INTENT(IN   ) :: x                  !< Continuous states at operating point
   TYPE(MD_DiscreteStateType),        INTENT(IN   ) :: xd                 !< Discrete states at operating point
   TYPE(MD_ConstraintStateType),      INTENT(IN   ) :: z                  !< Constraint states at operating point
   TYPE(MD_OtherStateType),           INTENT(IN   ) :: OtherState         !< Other states at operating point
   TYPE(MD_OutputType),               INTENT(INOUT) :: y                  !< Output (change to inout if a mesh copy is required); Output fields are not used by this routine, but type is available here so that mesh parameter information (i.e., connectivity) does not have to be recalculated for dYdx.
   TYPE(MD_MiscVarType),              INTENT(INOUT) :: m                  !< Misc/optimization variables
   INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat    !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dYdz(:,:)  !< Partial derivatives of output functions (Y) with respect to the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXdz(:,:)  !< Partial derivatives of continuous state functions (X) with respect to the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dXddz(:,:) !< Partial derivatives of discrete state functions (Xd) with respect to the constraint states (z) [intent in to avoid deallocation]
   REAL(R8Ki), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dZdz(:,:)  !< Partial derivatives of constraint state functions (Z) with respect to the constraint states (z) [intent in to avoid deallocation]
   ! local variables
   character(*), parameter                                       :: RoutineName = 'MD_JacobianPConstrState'
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ''
   IF ( PRESENT( dYdz ) ) THEN
   END IF
   IF ( PRESENT( dXdz ) ) THEN
      if (allocated(dXdz)) deallocate(dXdz)
   END IF
   IF ( PRESENT( dXddz ) ) THEN
      if (allocated(dXddz)) deallocate(dXddz)
   END IF
   IF ( PRESENT(dZdz) ) THEN
   END IF
END SUBROUTINE MD_JacobianPConstrState
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Routine to pack the data structures representing the operating points into arrays for linearization.
SUBROUTINE MD_GetOP( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg, u_op, y_op, x_op, dx_op, xd_op, z_op )
   REAL(DbKi),                        INTENT(IN   ) :: t          !< Time in seconds at operating point
   TYPE(MD_InputType),                INTENT(INOUT) :: u          !< Inputs at operating point (may change to inout if a mesh copy is required)
   TYPE(MD_ParameterType),            INTENT(IN   ) :: p          !< Parameters
   TYPE(MD_ContinuousStateType),      INTENT(IN   ) :: x          !< Continuous states at operating point
   TYPE(MD_DiscreteStateType),        INTENT(IN   ) :: xd         !< Discrete states at operating point
   TYPE(MD_ConstraintStateType),      INTENT(IN   ) :: z          !< Constraint states at operating point
   TYPE(MD_OtherStateType),           INTENT(IN   ) :: OtherState !< Other states at operating point
   TYPE(MD_OutputType),               INTENT(IN   ) :: y          !< Output at operating point
   TYPE(MD_MiscVarType),              INTENT(INOUT) :: m          !< Misc/optimization variables
   INTEGER(IntKi),                    INTENT(  OUT) :: ErrStat    !< Error status of the operation
   CHARACTER(*),                      INTENT(  OUT) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(ReKi), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: u_op(:)    !< values of linearized inputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: y_op(:)    !< values of linearized outputs
   REAL(ReKi), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: x_op(:)    !< values of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: dx_op(:)   !< values of first time derivatives of linearized continuous states
   REAL(ReKi), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: xd_op(:)   !< values of linearized discrete states
   REAL(ReKi), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: z_op(:)    !< values of linearized constraint states
   ! Local
   INTEGER(IntKi)                                                :: idx, i
   INTEGER(IntKi)                                                :: nu
   INTEGER(IntKi)                                                :: ny
   INTEGER(IntKi)                                                :: ErrStat2
   CHARACTER(ErrMsgLen)                                          :: ErrMsg2
   CHARACTER(*), PARAMETER                                       :: RoutineName = 'MD_GetOP'
   LOGICAL                                                       :: FieldMask(FIELDMASK_SIZE)
   TYPE(MD_ContinuousStateType)                                  :: dx          ! derivative of continuous states at operating point
   ErrStat = ErrID_None
   ErrMsg  = ''
   ! inputs
   IF ( PRESENT( u_op ) ) THEN
      nu = size(p%Jac_u_indx,1) + u%CoupledKinematics(1)%NNodes * 6  ! Jac_u_indx has 3 orientation angles, but the OP needs the full 9 elements of the DCM (thus 6 more per node)
      if (.not. allocated(u_op)) then
         call AllocAry(u_op, nu, 'u_op', ErrStat2, ErrMsg2); if(Failed()) return
      end if
      idx = 1
      FieldMask = .false.
      FieldMask(MASKID_TranslationDisp) = .true.
      FieldMask(MASKID_Orientation)     = .true.
      FieldMask(MASKID_TranslationVel)  = .true.
      FieldMask(MASKID_RotationVel)     = .true.
      FieldMask(MASKID_TranslationAcc)  = .true.
      FieldMask(MASKID_RotationAcc)     = .true.
      ! fill in the u_op values from the input mesh
      call PackMotionMesh(u%CoupledKinematics(1), u_op, idx, FieldMask=FieldMask)
      
      ! now do the active tensioning commands if there are any
      if (allocated(u%DeltaL)) then
         do i=1,size(u%DeltaL)
            u_op(idx) = u%DeltaL(i)
            idx = idx + 1
            u_op(idx) = u%DeltaLdot(i)
            idx = idx + 1
         end do
      endif
   END IF
   ! outputs
   IF ( PRESENT( y_op ) ) THEN
      ny = p%Jac_ny + y%CoupledLoads(1)%NNodes * 6  ! Jac_ny has 3 orientation angles, but the OP needs the full 9 elements of the DCM (thus 6 more per node)
      if (.not. allocated(y_op)) then
         call AllocAry(y_op, ny, 'y_op', ErrStat2, ErrMsg2); if(Failed()) return
      end if
      idx = 1
      call PackLoadMesh(y%CoupledLoads(1), y_op, idx)
      do i=1,p%NumOuts
         y_op(idx) = y%WriteOutput(i)
         idx = idx + 1
      end do
   END IF
   ! states
   IF ( PRESENT( x_op ) ) THEN
      if (.not. allocated(x_op)) then
         call AllocAry(x_op, p%Jac_nx,'x_op',ErrStat2,ErrMsg2); if (Failed()) return
      end if
      do i=1, p%Jac_nx
         x_op(i) = x%states(p%dxIdx_map2_xStateIdx(i))      ! x for lin is different order, so use mapping
      end do
   END IF
   ! state derivatives?
   IF ( PRESENT( dx_op ) ) THEN
      if (.not. allocated(dx_op)) then
         call AllocAry(dx_op, p%Jac_nx,'dx_op',ErrStat2,ErrMsg2); if(failed()) return
      end if
      call MD_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, m, dx, ErrStat2, ErrMsg2 ) ; if(Failed()) return
      do i=1, p%Jac_nx
         dx_op(i) = dx%states(p%dxIdx_map2_xStateIdx(i))    ! x for lin is different order, so use mapping
      end do
   END IF
   IF ( PRESENT( xd_op ) ) THEN
      ! pass
   END IF
   IF ( PRESENT( z_op ) ) THEN
      ! pass
   END IF
   call CleanUp()
contains
   logical function Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'MD_GetOP') 
        Failed =  ErrStat >= AbortErrLev
        if (Failed) call CleanUp()
   end function Failed

   subroutine CleanUp()
      call MD_DestroyContState(dx, ErrStat2, ErrMsg2);
   end subroutine
END SUBROUTINE MD_GetOP



!====================================================================================================
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> This routine initializes the array that maps rows/columns of the Jacobian to specific mesh fields.
!! Do not change the order of this packing without changing subroutines calculating dXdx etc (MD_Compute_dX)
SUBROUTINE MD_Init_Jacobian(Init, p, u, y, m, InitOut, ErrStat, ErrMsg)
   TYPE(MD_InitInputType)            , INTENT(IN   ) :: Init                  !< Init
   TYPE(MD_ParameterType)            , INTENT(INOUT) :: p                     !< parameters
   TYPE(MD_InputType)                , INTENT(IN   ) :: u                     !< inputs
   TYPE(MD_OutputType)               , INTENT(IN   ) :: y                     !< outputs
   TYPE(MD_MiscVarType)              , INTENT(INOUT) :: m                     !< misc variables <<<<<<<<
   TYPE(MD_InitOutputType)           , INTENT(INOUT) :: InitOut               !< Initialization output data (for Jacobian row/column names)
   INTEGER(IntKi)                    , INTENT(  OUT) :: ErrStat               !< Error status of the operation
   CHARACTER(*)                      , INTENT(  OUT) :: ErrMsg                !< Error message if ErrStat /= ErrID_None
   
   INTEGER(IntKi)                                    :: ErrStat2
   CHARACTER(ErrMsgLen)                              :: ErrMsg2
   CHARACTER(*), PARAMETER                           :: RoutineName = 'SD_Init_Jacobian'
   real(ReKi) :: dx, dy, dz, maxDim
   
   INTEGER(IntKi)                                    :: l, I
   real(ReKi)                                        :: dl_slack     ! how much a given line segment is stretched [m] 
   real(ReKi)                                        :: dl_slack_min ! minimum change in a node position for the least-strained segment in the simulation to go slack [m]
   
   
   ! local variables:
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   !! --- System dimension
   !dx = maxval(Init%Nodes(:,2))- minval(Init%Nodes(:,2))
   !dy = maxval(Init%Nodes(:,3))- minval(Init%Nodes(:,3))
   !dz = maxval(Init%Nodes(:,4))- minval(Init%Nodes(:,4))
   !maxDim = max(dx, dy, dz)
   
   
   ! Figure out appropriate transverse perturbation size to avoid slack segments
   dl_slack_min = 0.1_ReKi  ! start at 0.1 m
   
   do l = 1,p%nLines
      do I = 1, m%LineList(l)%N
         dl_slack = m%LineList(l)%lstr(I) - m%LineList(l)%l(I)
      
         ! store the smallest positive length margin to a segment going slack
         if (( dl_slack > 0.0_ReKi) .and. (dl_slack < dl_slack_min)) then
            dl_slack_min = dl_slack  
         end if
      end do
   end do
   
   dl_slack_min = 0.5*dl_slack_min  ! apply 0.5 safety factor
   
   !TODO: consider attachment radii to also produce a rotational perturbation size from the above
   
   
   ! --- System dimension
   call Init_Jacobian_y(); if (Failed()) return
   call Init_Jacobian_x(); if (Failed()) return
   call Init_Jacobian_u(); if (Failed()) return

contains
   LOGICAL FUNCTION Failed()
        call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'SD_Init_Jacobian') 
        Failed =  ErrStat >= AbortErrLev
   END FUNCTION Failed
   
   !> This routine initializes the Jacobian parameters and initialization outputs for the linearized outputs.
   SUBROUTINE Init_Jacobian_y()
      INTEGER(IntKi) :: index_next, i
      
      ! Number of outputs
      p%Jac_ny = y%CoupledLoads(1)%nNodes * 6     & ! 3 forces + 3 moments at each node (moments may be zero)
               + p%NumOuts                       ! WriteOutput values 
      ! Storage info for each output (names, rotframe)
      call AllocAry(InitOut%LinNames_y, p%Jac_ny, 'LinNames_y',ErrStat2,ErrMsg2); if(ErrStat2/=ErrID_None) return
      call AllocAry(InitOut%RotFrame_y, p%Jac_ny, 'RotFrame_y',ErrStat2,ErrMsg2); if(ErrStat2/=ErrID_None) return
      ! Names
      index_next = 1
      call PackLoadMesh_Names(  y%CoupledLoads(1), 'LinNames_y', InitOut%LinNames_y, index_next)  ! <<< should a specific name be provided here?
      do i=1,p%NumOuts
         InitOut%LinNames_y(i+index_next-1) = trim(InitOut%WriteOutputHdr(i))//', '//trim(InitOut%WriteOutputUnt(i))
      end do
      
      InitOut%RotFrame_y(:) = .false.
   END SUBROUTINE Init_Jacobian_y

   !> This routine initializes the Jacobian parameters and initialization outputs for the linearized continuous states.
   SUBROUTINE Init_Jacobian_x()
      INTEGER(IntKi) :: idx      ! index into the LinNames_x array
      INTEGER(IntKi) :: i
      INTEGER(IntKi) :: l
      INTEGER(IntKi) :: N

      
      p%Jac_nx = m%Nx ! size of (continuous) state vector (includes the first derivatives)
      
      ! allocate space for the row/column names and for perturbation sizes
      CALL AllocAry(InitOut%LinNames_x    , p%Jac_nx, 'LinNames_x'            , ErrStat2, ErrMsg2); if(ErrStat/=ErrID_None) return
      CALL AllocAry(InitOut%RotFrame_x    , p%Jac_nx, 'RotFrame_x'            , ErrStat2, ErrMsg2); if(ErrStat/=ErrID_None) return
      CALL AllocAry(InitOut%DerivOrder_x  , p%Jac_nx, 'DerivOrder_x'          , ErrStat2, ErrMsg2); if(ErrStat/=ErrID_None) return
      CALL AllocAry(p%dx                  , p%Jac_nx, 'p%dx'                  , ErrStat2, ErrMsg2); if(ErrStat/=ErrID_None) return
      CALL AllocAry(p%dxIdx_map2_xStateIdx, p%Jac_nx, 'p%dxIdx_map2_xStateIdx', ErrStat2, ErrMsg2); if(ErrStat/=ErrID_None) return

      p%dxIdx_map2_xStateIdx = 0_IntKi ! all values should be overwritten by logic below

      ! set linearization output names and default perturbations, p%dx:
      !  NOTE: the order is different than the order of the internal states.  This is to
      !        match what the OpenFAST framework is expecting: all positions first, then all
      !        derviatives of positions (velocity terms) second.  This adds slight complexity
      !        here, but considerably simplifies post processing of the full OpenFAST results
      !        for linearization.
      !        The p%dxIdx_map2_xStateIdx array holds the index for the x%states array
      !        corresponding to the current jacobian index.

      !-----------------
      ! position states
      !-----------------
      idx = 0
      ! Free bodies
      DO l = 1,p%nFreeBodies                 ! Body m%BodyList(m%FreeBodyIs(l))
         p%dx(idx+1:idx+3) = dl_slack_min    ! body displacement [m]
         p%dx(idx+4:idx+6) = 0.02            ! body rotation [rad]
         ! corresponds to state indices: (m%BodyStateIs1(l)+6:m%BodyStateIs1(l)+11)
         InitOut%LinNames_x(idx+1) = 'Body '//trim(num2lstr(m%FreeBodyIs(l)))//' Px, m'
         InitOut%LinNames_x(idx+2) = 'Body '//trim(num2lstr(m%FreeBodyIs(l)))//' Py, m'
         InitOut%LinNames_x(idx+3) = 'Body '//trim(num2lstr(m%FreeBodyIs(l)))//' Pz, m'
         InitOut%LinNames_x(idx+4) = 'Body '//trim(num2lstr(m%FreeBodyIs(l)))//' rot_x, rad'
         InitOut%LinNames_x(idx+5) = 'Body '//trim(num2lstr(m%FreeBodyIs(l)))//' rot_y, rad'
         InitOut%LinNames_x(idx+6) = 'Body '//trim(num2lstr(m%FreeBodyIs(l)))//' rot_z, rad'
         p%dxIdx_map2_xStateIdx(idx+1) = m%BodyStateIs1(l)+6         ! x%state index for Px
         p%dxIdx_map2_xStateIdx(idx+2) = m%BodyStateIs1(l)+7         ! x%state index for Py
         p%dxIdx_map2_xStateIdx(idx+3) = m%BodyStateIs1(l)+8         ! x%state index for Pz
         p%dxIdx_map2_xStateIdx(idx+4) = m%BodyStateIs1(l)+9         ! x%state index for rot_x
         p%dxIdx_map2_xStateIdx(idx+5) = m%BodyStateIs1(l)+10        ! x%state index for rot_y
         p%dxIdx_map2_xStateIdx(idx+6) = m%BodyStateIs1(l)+11        ! x%state index for rot_z
         idx = idx + 6
      END DO      

      ! Rods
      DO l = 1,p%nFreeRods                   ! Rod m%RodList(m%FreeRodIs(l))
         if (m%RodList(m%FreeRodIs(l))%typeNum == 1) then  ! pinned rod
            p%dx(idx+1:idx+3) = 0.02         ! rod rotation [rad]
            ! corresponds to state indices: (m%RodStateIs1(l)+3:m%RodStateIs1(l)+5)
            InitOut%LinNames_x(idx+1) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' rot_x, rad'
            InitOut%LinNames_x(idx+2) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' rot_y, rad'
            InitOut%LinNames_x(idx+3) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' rot_z, rad'   
            p%dxIdx_map2_xStateIdx(idx+4) = m%RodStateIs1(l)+3          ! x%state index for rot_x
            p%dxIdx_map2_xStateIdx(idx+5) = m%RodStateIs1(l)+4          ! x%state index for rot_y
            p%dxIdx_map2_xStateIdx(idx+6) = m%RodStateIs1(l)+5          ! x%state index for rot_z
            idx = idx + 3
         else                                ! free rod
            p%dx(idx+1:idx+3) = dl_slack_min ! rod displacement [m]
            p%dx(idx+4:idx+6) = 0.02         ! rod rotation [rad]
            ! corresponds to state indices: (m%RodStateIs1(l)+6:m%RodStateIs1(l)+11)
            InitOut%LinNames_x(idx+1) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' Px, m'
            InitOut%LinNames_x(idx+2) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' Py, m'
            InitOut%LinNames_x(idx+3) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' Pz, m'
            InitOut%LinNames_x(idx+4) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' rot_x, rad'
            InitOut%LinNames_x(idx+5) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' rot_y, rad'
            InitOut%LinNames_x(idx+6) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' rot_z, rad'   
            p%dxIdx_map2_xStateIdx(idx+1) = m%RodStateIs1(l)+6          ! x%state index for Px
            p%dxIdx_map2_xStateIdx(idx+2) = m%RodStateIs1(l)+7          ! x%state index for Py
            p%dxIdx_map2_xStateIdx(idx+3) = m%RodStateIs1(l)+8          ! x%state index for Pz
            p%dxIdx_map2_xStateIdx(idx+4) = m%RodStateIs1(l)+9          ! x%state index for rot_x
            p%dxIdx_map2_xStateIdx(idx+5) = m%RodStateIs1(l)+10         ! x%state index for rot_y
            p%dxIdx_map2_xStateIdx(idx+6) = m%RodStateIs1(l)+11         ! x%state index for rot_z
            idx = idx + 6
         end if
      END DO      

      ! Free Connnections
      DO l = 1,p%nFreeCons                   ! Point m%ConnectList(m%FreeConIs(l))
         ! corresponds to state indices: (m%ConStateIs1(l)+3:m%ConStateIs1(l)+5)
         p%dx(idx+1:idx+3) = dl_slack_min    ! point displacement [m]
         InitOut%LinNames_x(idx+1) = 'Point '//trim(num2lstr(m%FreeConIs(l)))//' Px, m'
         InitOut%LinNames_x(idx+2) = 'Point '//trim(num2lstr(m%FreeConIs(l)))//' Py, m'
         InitOut%LinNames_x(idx+3) = 'Point '//trim(num2lstr(m%FreeConIs(l)))//' Pz, m'
         p%dxIdx_map2_xStateIdx(idx+1) = m%ConStateIs1(l)+3          ! x%state index for Px
         p%dxIdx_map2_xStateIdx(idx+2) = m%ConStateIs1(l)+4          ! x%state index for Py
         p%dxIdx_map2_xStateIdx(idx+3) = m%ConStateIs1(l)+5          ! x%state index for Pz
         idx = idx + 3
      END DO

      ! Lines
      DO l = 1,p%nLines                      ! Line m%LineList(l)         
         ! corresponds to state indices: (m%LineStateIs1(l)+3*N-3:m%LineStateIs1(l)+6*N-7) -- NOTE: end nodes not included 
         N = m%LineList(l)%N                 ! number of segments in the line
         DO i = 0,N-2
            p%dx(idx+1:idx+3) = dl_slack_min ! line internal node displacement [m]
            InitOut%LinNames_x(idx+1) = 'Line '//trim(num2lstr(l))//' node '//trim(num2lstr(i+1))//' Px, m'
            InitOut%LinNames_x(idx+2) = 'Line '//trim(num2lstr(l))//' node '//trim(num2lstr(i+1))//' Py, m'
            InitOut%LinNames_x(idx+3) = 'Line '//trim(num2lstr(l))//' node '//trim(num2lstr(i+1))//' Pz, m'
            p%dxIdx_map2_xStateIdx(idx+1) = m%LineStateIs1(l)+3*N+3*i-3 ! x%state index for Px
            p%dxIdx_map2_xStateIdx(idx+2) = m%LineStateIs1(l)+3*N+3*i-2 ! x%state index for Py
            p%dxIdx_map2_xStateIdx(idx+3) = m%LineStateIs1(l)+3*N+3*i-1 ! x%state index for Pz
            idx = idx + 3
         END DO         
      END DO

      !-----------------
      ! velocity states
      !-----------------
      ! Free bodies
      DO l = 1,p%nFreeBodies                 ! Body m%BodyList(m%FreeBodyIs(l))
         ! corresponds to state indices: (m%BodyStateIs1(l):m%BodyStateIs1(l)+5)
         p%dx(idx+1:idx+3) = 0.1             ! body translational velocity [m/s]
         p%dx(idx+4:idx+6) = 0.1             ! body rotational velocity [rad/s]
         InitOut%LinNames_x(idx+1) = 'Body '//trim(num2lstr(m%FreeBodyIs(l)))//' Vx, m/s'
         InitOut%LinNames_x(idx+2) = 'Body '//trim(num2lstr(m%FreeBodyIs(l)))//' Vy, m/s'
         InitOut%LinNames_x(idx+3) = 'Body '//trim(num2lstr(m%FreeBodyIs(l)))//' Vz, m/s'
         InitOut%LinNames_x(idx+4) = 'Body '//trim(num2lstr(m%FreeBodyIs(l)))//' omega_x, rad/s'
         InitOut%LinNames_x(idx+5) = 'Body '//trim(num2lstr(m%FreeBodyIs(l)))//' omega_y, rad/s'
         InitOut%LinNames_x(idx+6) = 'Body '//trim(num2lstr(m%FreeBodyIs(l)))//' omega_z, rad/s'
         p%dxIdx_map2_xStateIdx(idx+1) = m%BodyStateIs1(l)+0         ! x%state index for Rx
         p%dxIdx_map2_xStateIdx(idx+2) = m%BodyStateIs1(l)+1         ! x%state index for Ry
         p%dxIdx_map2_xStateIdx(idx+3) = m%BodyStateIs1(l)+2         ! x%state index for Rz
         p%dxIdx_map2_xStateIdx(idx+4) = m%BodyStateIs1(l)+3         ! x%state index for omega_x
         p%dxIdx_map2_xStateIdx(idx+5) = m%BodyStateIs1(l)+4         ! x%state index for omega_y
         p%dxIdx_map2_xStateIdx(idx+6) = m%BodyStateIs1(l)+5         ! x%state index for omega_z
         idx = idx + 6
      END DO      

      ! Rods
      DO l = 1,p%nFreeRods                   ! Rod m%RodList(m%FreeRodIs(l))
         if (m%RodList(m%FreeRodIs(l))%typeNum == 1) then ! pinned rod
            ! corresponds to state indices: (m%RodStateIs1(l):m%RodStateIs1(l)+2)
            p%dx(idx+1:idx+3) = 0.1          ! body rotational velocity [rad/s]
            InitOut%LinNames_x(idx+1) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' omega_x, rad/s'
            InitOut%LinNames_x(idx+2) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' omega_y, rad/s'
            InitOut%LinNames_x(idx+3) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' omega_z, rad/s'
            p%dxIdx_map2_xStateIdx(idx+1) = m%RodStateIs1(l)+0          ! x%state index for Vx
            p%dxIdx_map2_xStateIdx(idx+2) = m%RodStateIs1(l)+1          ! x%state index for Vy
            p%dxIdx_map2_xStateIdx(idx+3) = m%RodStateIs1(l)+2          ! x%state index for Vz
            idx = idx + 3
         else                                ! free rod
            ! corresponds to state indices: (m%RodStateIs1(l):m%RodStateIs1(l)+5)
            p%dx(idx+1:idx+3) = 0.1          ! body translational velocity [m/s]
            p%dx(idx+4:idx+6) = 0.02         ! body rotational velocity [rad/s]
            InitOut%LinNames_x(idx+1) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' Vx, m/s'
            InitOut%LinNames_x(idx+2) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' Vy, m/s'
            InitOut%LinNames_x(idx+3) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' Vz, m/s'
            InitOut%LinNames_x(idx+4) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' omega_x, rad/s'
            InitOut%LinNames_x(idx+5) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' omega_y, rad/s'
            InitOut%LinNames_x(idx+6) = 'Rod '//trim(num2lstr(m%FreeRodIs(l)))//' omega_z, rad/s'
            p%dxIdx_map2_xStateIdx(idx+1) = m%RodStateIs1(l)+0          ! x%state index for Vx
            p%dxIdx_map2_xStateIdx(idx+2) = m%RodStateIs1(l)+1          ! x%state index for Vy
            p%dxIdx_map2_xStateIdx(idx+3) = m%RodStateIs1(l)+2          ! x%state index for Vz
            p%dxIdx_map2_xStateIdx(idx+4) = m%RodStateIs1(l)+3          ! x%state index for omega_x
            p%dxIdx_map2_xStateIdx(idx+5) = m%RodStateIs1(l)+4          ! x%state index for omega_y
            p%dxIdx_map2_xStateIdx(idx+6) = m%RodStateIs1(l)+5          ! x%state index for omega_z
            idx = idx + 6
         end if
      END DO      

      ! Free Connnections
      DO l = 1,p%nFreeCons                   ! Point m%ConnectList(m%FreeConIs(l))
         ! corresponds to state indices: (m%ConStateIs1(l):m%ConStateIs1(l)+2)
         p%dx(idx+1:idx+3) = 0.1             ! point translational velocity [m/s]
         InitOut%LinNames_x(idx+1) = 'Point '//trim(num2lstr(m%FreeConIs(l)))//' Vx, m/s'
         InitOut%LinNames_x(idx+2) = 'Point '//trim(num2lstr(m%FreeConIs(l)))//' Vy, m/s'
         InitOut%LinNames_x(idx+3) = 'Point '//trim(num2lstr(m%FreeConIs(l)))//' Vz, m/s'
         p%dxIdx_map2_xStateIdx(idx+1) = m%ConStateIs1(l)+0          ! x%state index for Vx
         p%dxIdx_map2_xStateIdx(idx+2) = m%ConStateIs1(l)+1          ! x%state index for Vy
         p%dxIdx_map2_xStateIdx(idx+3) = m%ConStateIs1(l)+2          ! x%state index for Vz
         idx = idx + 3
      END DO

      ! Lines
      DO l = 1,p%nLines                      ! Line m%LineList(l)         
         ! corresponds to state indices: (m%LineStateIs1(l):m%LineStateIs1(l)+3*N-4) -- NOTE: end nodes not included
         N = m%LineList(l)%N                 ! number of segments in the line
         DO i = 0,N-2
            p%dx(idx+1:idx+3) = 0.1          ! line internal node translational velocity [m/s]
            InitOut%LinNames_x(idx+1) = 'Line '//trim(num2lstr(l))//' node '//trim(num2lstr(i+1))//' Vx, m/s'
            InitOut%LinNames_x(idx+2) = 'Line '//trim(num2lstr(l))//' node '//trim(num2lstr(i+1))//' Vy, m/s'
            InitOut%LinNames_x(idx+3) = 'Line '//trim(num2lstr(l))//' node '//trim(num2lstr(i+1))//' Vz, m/s'
            p%dxIdx_map2_xStateIdx(idx+1) = m%LineStateIs1(l)+3*i+0  ! x%state index for Vx
            p%dxIdx_map2_xStateIdx(idx+2) = m%LineStateIs1(l)+3*i+1  ! x%state index for Vy
            p%dxIdx_map2_xStateIdx(idx+3) = m%LineStateIs1(l)+3*i+2  ! x%state index for Vz
            idx = idx + 3
         END DO
      END DO

      ! If a summary file is ever made...
      !  !Formatting may be needed to make it pretty
      !  if(UnSum > 0) then
      !     write(UnSum,*) ' Lin_Jac_x       idx        x%state idx'
      !     do i=1,p%Jac_nx
      !        write(UnSum,*) InitOut%LinNames_x(i),'  ',i,'   ',p%dxIdx_map2_xStateIdx(i)
      !     enddo
      !  endif

      InitOut%RotFrame_x   = .false.
      InitOut%DerivOrder_x = 2      
   END SUBROUTINE Init_Jacobian_x

   SUBROUTINE Init_Jacobian_u()
      INTEGER(IntKi) :: i, j, idx, nu, i_meshField
      character(10)  :: LinStr      ! for noting which line a DeltaL control is attached to
      logical        :: LinCtrl     ! Is the current DeltaL channel associated with a line?
      ! Number of inputs
      i = 0
      if (allocated(u%DeltaL))   i=size(u%DeltaL)
      nu = u%CoupledKinematics(1)%nNodes * 18 &   ! 3 Translation Displacements + 3 orientations + 6 velocities + 6 accelerations at each node <<<<<<<
         + i*2                                 ! a deltaL and rate of change for each active tension control channel
      
      ! --- Info of linearized inputs (Names, RotFrame, IsLoad)
      call AllocAry(InitOut%LinNames_u, nu, 'LinNames_u', ErrStat2, ErrMsg2); if(ErrStat2/=ErrID_None) return
      call AllocAry(InitOut%RotFrame_u, nu, 'RotFrame_u', ErrStat2, ErrMsg2); if(ErrStat2/=ErrID_None) return
      call AllocAry(InitOut%IsLoad_u  , nu, 'IsLoad_u'  , ErrStat2, ErrMsg2); if(ErrStat2/=ErrID_None) return
      
      InitOut%IsLoad_u   = .false. ! None of MoorDyn's inputs are loads
      InitOut%RotFrame_u = .false. ! every input is on a mesh, which stores values in the global (not rotating) frame
      
      idx = 1
      call PackMotionMesh_Names(u%CoupledKinematics(1), 'CoupledKinematics', InitOut%LinNames_u, idx) ! all 6 motion fields
      
      ! --- Jac_u_indx:  matrix to store index to help us figure out what the ith value of the u vector really means
      ! (see perturb_u ... these MUST match )
      ! column 1 indicates module's mesh and field
      ! column 2 indicates the first index (x-y-z component) of the field
      ! column 3 is the node
      call allocAry( p%Jac_u_indx, nu, 3, 'p%Jac_u_indx', ErrStat2, ErrMsg2); if(ErrStat2/=ErrID_None) return
      p%Jac_u_indx = 0  ! initialize to zero
      idx = 1
      !Module/Mesh/Field: u%CoupledKinematics(1)%TranslationDisp  = 1;
      !Module/Mesh/Field: u%CoupledKinematics(1)%Orientation      = 2;
      !Module/Mesh/Field: u%CoupledKinematics(1)%TranslationVel   = 3;
      !Module/Mesh/Field: u%CoupledKinematics(1)%RotationVel      = 4;
      !Module/Mesh/Field: u%CoupledKinematics(1)%TranslationAcc   = 5;
      !Module/Mesh/Field: u%CoupledKinematics(1)%RotationAcc      = 6;
      do i_meshField = 1,6
         do i=1,u%CoupledKinematics(1)%nNodes
            do j=1,3
               p%Jac_u_indx(idx,1) =  i_meshField     ! mesh field type (indicated by 1-6)
               p%Jac_u_indx(idx,2) =  j               ! x, y, or z
               p%Jac_u_indx(idx,3) =  i               ! node
               idx = idx + 1
            end do !j
         end do !i
      end do
      ! now do the active tensioning commands if there are any
      if (allocated(u%DeltaL)) then
         do i=1,size(u%DeltaL)            ! Signals may be passed in without being requested for control
            ! Figure out if this DeltaL control channel is associated with a line or multiple or none and label
            LinCtrl = .FALSE.
            LinStr = '(lines: '
            do J=1,p%NLines
               if (m%LineList(J)%CtrlChan == i) then
                  LinCtrl = .TRUE.
                  LinStr = LinStr//trim(num2lstr(i))//' '
               endif
            enddo
            if (      LinCtrl)   LinStr = LinStr//' )'
            if (.not. LinCtrl)   LinStr = '(lines: none)'

            p%Jac_u_indx(idx,1) =  10              ! 10-11 mean active tension changes (10: deltaL; 11: deltaLdot)
            p%Jac_u_indx(idx,2) =  0               ! not used
            p%Jac_u_indx(idx,3) =  i               ! indicates DeltaL entry number 
            InitOut%LinNames_u(idx) = 'CtrlChan DeltaL '//trim(num2lstr(i))//', m '//trim(LinStr)
            idx = idx + 1
         
            p%Jac_u_indx(idx,1) =  11
            p%Jac_u_indx(idx,2) =  0
            p%Jac_u_indx(idx,3) =  i               
            InitOut%LinNames_u(idx) = 'CtrlChan DeltaLdot '//trim(num2lstr(i))//', m/s'//trim(LinStr)
            idx = idx + 1
         end do
      endif

      ! --- Default perturbations, p%du:
      call allocAry( p%du, 11, 'p%du', ErrStat2, ErrMsg2); if(ErrStat2/=ErrID_None) return
      p%du( 1) = dl_slack_min  ! u%CoupledKinematics(1)%TranslationDisp  = 1;
      p%du( 2) = 0.1_ReKi      ! u%CoupledKinematics(1)%Orientation      = 2;
      p%du( 3) = 0.1_ReKi      ! u%CoupledKinematics(1)%TranslationVel   = 3;
      p%du( 4) = 0.1_ReKi      ! u%CoupledKinematics(1)%RotationVel      = 4;
      p%du( 5) = 0.1_ReKi      ! u%CoupledKinematics(1)%TranslationAcc   = 5;
      p%du( 6) = 0.1_ReKi      ! u%CoupledKinematics(1)%RotationAcc      = 6;
      p%du(10) = dl_slack_min  ! deltaL [m]
      p%du(11) = 0.2_ReKi      ! deltaLdot [m/s] 
   END SUBROUTINE Init_Jacobian_u

END SUBROUTINE MD_Init_Jacobian
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine perturbs the nth element of the u array (and mesh/field it corresponds to)
!! Do not change this without making sure subroutine MD_init_jacobian is consistant with this routine!
SUBROUTINE MD_Perturb_u( p, n, perturb_sign, u, du )
   TYPE(MD_ParameterType)              , INTENT(IN   ) :: p            !< parameters
   INTEGER( IntKi )                    , INTENT(IN   ) :: n            !< number of array element to use
   INTEGER( IntKi )                    , INTENT(IN   ) :: perturb_sign !< +1 or -1 (value to multiply perturbation by; positive or negative difference)
   TYPE(MD_InputType)                  , INTENT(INOUT) :: u            !< perturbed MD inputs
   REAL( R8Ki )                        , INTENT(  OUT) :: du           !< amount that specific input was perturbed
   ! local variables
   INTEGER :: fieldIndx
   INTEGER :: node
   fieldIndx = p%Jac_u_indx(n,2)
   node      = p%Jac_u_indx(n,3)
   du = p%du(  p%Jac_u_indx(n,1) )
   ! determine which mesh we're trying to perturb and perturb the input:
   SELECT CASE( p%Jac_u_indx(n,1) )
   CASE ( 1)
      u%CoupledKinematics(1)%TranslationDisp( fieldIndx,node) = u%CoupledKinematics(1)%TranslationDisp( fieldIndx,node) + du * perturb_sign
   CASE ( 2)
      CALL PerturbOrientationMatrix( u%CoupledKinematics(1)%Orientation(:,:,node), du * perturb_sign, fieldIndx, UseSmlAngle=.true. )
   CASE ( 3)
      u%CoupledKinematics(1)%TranslationVel( fieldIndx,node) = u%CoupledKinematics(1)%TranslationVel( fieldIndx,node) + du * perturb_sign
   CASE ( 4)
      u%CoupledKinematics(1)%RotationVel(fieldIndx,node) = u%CoupledKinematics(1)%RotationVel(fieldIndx,node) + du * perturb_sign
   CASE ( 5)
      u%CoupledKinematics(1)%TranslationAcc( fieldIndx,node) = u%CoupledKinematics(1)%TranslationAcc( fieldIndx,node) + du * perturb_sign
   CASE ( 6)
      u%CoupledKinematics(1)%RotationAcc(fieldIndx,node) = u%CoupledKinematics(1)%RotationAcc(fieldIndx,node) + du * perturb_sign
   CASE (10)
      u%deltaL(node) = u%deltaL(node) + du * perturb_sign
   CASE (11)
      u%deltaLdot(node) = u%deltaLdot(node) + du * perturb_sign
   END SELECT
END SUBROUTINE MD_Perturb_u
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine uses values of two output types to compute an array of differences.
!! Do not change this packing without making sure subroutine MD_init_jacobian is consistant with this routine!
SUBROUTINE MD_Compute_dY(p, y_p, y_m, delta, dY)
   TYPE(MD_ParameterType)            , INTENT(IN   ) :: p     !< parameters
   TYPE(MD_OutputType)               , INTENT(IN   ) :: y_p   !< MD outputs at \f$ u + \Delta_p u \f$ or \f$ z + \Delta_p z \f$ (p=plus)
   TYPE(MD_OutputType)               , INTENT(IN   ) :: y_m   !< MD outputs at \f$ u - \Delta_m u \f$ or \f$ z - \Delta_m z \f$ (m=minus)
   REAL(R8Ki)                        , INTENT(IN   ) :: delta !< difference in inputs or states \f$ delta_p = \Delta_p u \f$ or \f$ delta_p = \Delta_p x \f$
   REAL(R8Ki)                        , INTENT(INOUT) :: dY(:) !< column of dYdu or dYdx: \f$ \frac{\partial Y}{\partial u_i} = \frac{y_p - y_m}{2 \, \Delta u}\f$ or \f$ \frac{\partial Y}{\partial z_i} = \frac{y_p - y_m}{2 \, \Delta x}\f$
   ! local variables:
   INTEGER(IntKi) :: i              ! loop over outputs
   INTEGER(IntKi) :: indx_first     ! index indicating next value of dY to be filled
   indx_first = 1
   call PackLoadMesh_dY(  y_p%CoupledLoads(1), y_m%CoupledLoads(1), dY, indx_first)
   !call PackMotionMesh_dY(y_p%Y2Mesh, y_m%Y2Mesh, dY, indx_first) ! all 6 motion fields
   do i=1,p%NumOuts
      dY(i+indx_first-1) = y_p%WriteOutput(i) - y_m%WriteOutput(i)
   end do
   dY = dY / (2.0_R8Ki*delta)
END SUBROUTINE MD_Compute_dY
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine perturbs the nth element of the x array (and mesh/field it corresponds to)
!! Do not change this without making sure subroutine MD_init_jacobian is consistant with this routine!
SUBROUTINE MD_Perturb_x( p, i, perturb_sign, x, dx )
   TYPE(MD_ParameterType)      , INTENT(IN   ) :: p            !< parameters
   INTEGER( IntKi )            , INTENT(IN   ) :: i            !< state array index number 
   INTEGER( IntKi )            , INTENT(IN   ) :: perturb_sign !< +1 or -1 (value to multiply perturbation by; positive or negative difference)
   TYPE(MD_ContinuousStateType), INTENT(INOUT) :: x            !< perturbed MD states
   REAL( R8Ki )                , INTENT(  OUT) :: dx           !< amount that specific state was perturbed

   dx=p%dx(i)
   x%states(i)    = x%states(i)    + dx * perturb_sign
END SUBROUTINE MD_Perturb_x
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine uses values of two output types to compute an array of differences.
!! Do not change this packing without making sure subroutine MD_init_jacobian is consistant with this routine!
SUBROUTINE MD_Compute_dX(p, x_p, x_m, delta, dX)
   TYPE(MD_ParameterType)      , INTENT(IN   ) :: p            !< parameters
   TYPE(MD_ContinuousStateType), INTENT(IN   ) :: x_p          !< <D continuous states at \f$ u + \Delta_p u \f$ or \f$ x + \Delta_p x \f$ (p=plus)
   TYPE(MD_ContinuousStateType), INTENT(IN   ) :: x_m          !< <D continuous states at \f$ u - \Delta_m u \f$ or \f$ x - \Delta_m x \f$ (m=minus)
   REAL(R8Ki)                  , INTENT(IN   ) :: delta        !< difference in inputs or states \f$ delta_p = \Delta_p u \f$ or \f$ delta_p = \Delta_p x \f$
   REAL(R8Ki)                  , INTENT(INOUT) :: dX(:)        !< column of dXdu or dXdx: \f$ \frac{\partial X}{\partial u_i} = \frac{x_p - x_m}{2 \, \Delta u}\f$ or \f$ \frac{\partial X}{\partial x_i} = \frac{x_p - x_m}{2 \, \Delta x}\f$
   INTEGER(IntKi) :: i ! loop over modes
   do i=1,p%Jac_nx   ! index to dX 
      ! NOTE: order of entries in dX is different than the x%states, so mapping is required
      dX(i) = x_p%states(p%dxIdx_map2_xStateIdx(i)) - x_m%states(p%dxIdx_map2_xStateIdx(i))
   end do
   dX = dX / (2.0_R8Ki*delta)
END SUBROUTINE MD_Compute_dX



END MODULE MoorDyn
