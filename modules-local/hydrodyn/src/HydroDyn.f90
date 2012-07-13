MODULE HydroDyn

   USE                                 NWTC_Library
   USE                                 SharedDataTypes    !bjj: add this to NWTC_Library????
   
   USE                                 FixedBottomSupportStructure
   USE                                 FloatingPlatform 
   USE                                 HD_Output
   USE                                 Waves       
   
   
   
   IMPLICIT NONE
   PRIVATE
 
   !-------------------------------------------------------------------------------------------------
   ! Definitions of public data parameters
   !-------------------------------------------------------------------------------------------------
   TYPE(ProgDesc), PARAMETER, PUBLIC :: HD_Prog = ProgDesc( 'HydroDyn', '(v1.00.01e-bjj, 06-Jul-2012)' )   ! the name/version/date of the hydrodynamics program

   !-------------------------------------------------------------------------------------------------
   ! Definitions of private parameters
   !-------------------------------------------------------------------------------------------------    
   INTEGER, PARAMETER                           :: Unknown_Type    = -1          ! Unknown/unset value
   INTEGER, PARAMETER                           :: FixedBtm_Type   =  1          ! fixed bottom/monopile tower support structure
   INTEGER, PARAMETER                           :: FloatPltfm_Type =  2          ! floating platform support structure

   !-------------------------------------------------------------------------------------------------
   ! Definitions of private types 
   !-------------------------------------------------------------------------------------------------
 
   TYPE Node_DataType
      REAL(ReKi)                                :: Position(3)             ! Center position of the element in the ground coordinate system (X,Y,Z)
      REAL(ReKi)                                :: DNode                   ! Length of the element
      REAL(ReKi)                                :: D                       ! Diameter (m)
      REAL(ReKi)                                :: CA                      ! Normalized hydrodynamic added mass   coefficient in Morison's equation (-)
      REAL(ReKi)                                :: CD                      ! Normalized hydrodynamic viscous drag coefficient in Morison's equation (-)
   END TYPE       

  
   !-------------------------------------------------------------------------------------------------
   ! Definitions of public types and routines
   !-------------------------------------------------------------------------------------------------
   PUBLIC                                       :: HD_Init                    ! Subroutine to initialize module
   PUBLIC                                       :: HD_CalculateLoads          ! Subroutine to calculate hydrodynamic loads
   PUBLIC                                       :: HD_Terminate               ! Subroutine to terminate module
   PUBLIC                                       :: HD_CheckLin                ! Function to determine if a linearization analysis can be performed with HydroDyn's settings
   
   
   PUBLIC                                       :: HD_GetUndisturbedWaveElev  ! Returns undisturbed wave elevation
   
   INTERFACE HD_GetValue
      MODULE PROCEDURE HD_GetValue_CHAR
      MODULE PROCEDURE HD_GetValue_AllOuts
   END INTERFACE HD_GetValue
      
      
   PUBLIC                                       :: HD_GetValue                ! Returns ancillary values that other codes may want to use
!   PUBLIC                                       :: HD_GetValue_CHAR           ! Returns ancillary values that other codes may want to use
!   PUBLIC                                       :: HD_GetValue_AllOuts        ! Returns ancillary values that other codes may want to use


   TYPE, PUBLIC :: HD_DataType                                                ! The type is public
      PRIVATE                                                                 ! but the data is private
      INTEGER                                   :: StrctType   = Unknown_Type ! Type of support structure (valid choices are parameters FloatPltfm_Type and FixedBtm_Type)
      INTEGER                                   :: UnSum       = -1           ! Unit number for the HydroDyn summary file   
      TYPE(Waves_DataType)                      :: Waves_data                 ! waves module data   
      TYPE(FltPtfm_DataType)                    :: FltPtfm_data               ! floating platform module data   
      
      TYPE(HD_OutputDataType)                   :: HDOut_Data                 ! output data
      REAL(DbKi)                                :: LastCalcTime = 0.0         ! The last time an output was calculated

      INTEGER                                   :: NElements    = 0           ! number of elements in the structure
      TYPE(Node_DataType), ALLOCATABLE          :: Node  (:)                  ! the elements/nodes that make up the structural discretization
      
      REAL(ReKi)                                :: MaxDiam  = 0.0             ! Maximum value of input array LRadAnch (stores a maximum value for graphics capabilities).  ! was MaxLRadAnch
      
   END TYPE HD_DataType


   TYPE, PUBLIC :: HD_InitDataType
      REAL(ReKi)                                :: Gravity                    ! Gravitational acceleration in m/s^2
      
      CHARACTER(1024)                           :: FileName                   ! Name of the input file
      CHARACTER(1024)                           :: OutRootName                ! The name of the root file (without extension) including the full path.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to OutRootName
   END TYPE


      ! Data types developed for HydroDyn

   TYPE, PUBLIC :: AllHydroMarkers
      TYPE(Marker), ALLOCATABLE  :: Substructure (:)  
   END TYPE
   
   TYPE, PUBLIC :: AllHydroLoads
      TYPE(Load),   ALLOCATABLE  :: Substructure (:)  
   END TYPE
      
   TYPE, PUBLIC :: HydroConfig
      TYPE(Marker)               :: Substructure
   END TYPE HydroConfig

   
!   REAL(ReKi), PARAMETER                        :: RePrcsn = PRECISION( 0.0_ReKi )  ! precision of ReKi values (used to test for equality)

   
CONTAINS
!====================================================================================================
SUBROUTINE HD_CalculateLoads( CurrentTime,  HD_Markers, HD_Data, HD_Loads, ErrStat )
!     This public subroutine computes hydrodynamic loads.
! 
!----------------------------------------------------------------------------------------------------

! bjj: do we not care about HD_ConfigMarkers? Shouldn't that also be passed for future reference?

      ! Passed variables

   REAL(DbKi),                 INTENT( IN    )     :: CurrentTime                ! Current simulation time (sec)    

   TYPE(AllHydroMarkers),      INTENT( IN    )     :: HD_Markers                 ! the markers of the loads calculated at this time 
   TYPE(HD_DataType),          INTENT( INOUT )     :: HD_Data                    ! the internal hydrodyn data
   TYPE(AllHydroLoads),        INTENT( INOUT )     :: HD_Loads                   ! the loads calculated at this time (it's INOUT so we don't have to reallocate space each call, but really is an OUTPUT)
   INTEGER,                    INTENT(   OUT )     :: ErrStat                    ! Returns zero if no errors were encountered
   
    
      ! local variables
   
   REAL(ReKi)                                      :: X    (6)                   ! Temp array holding 3 components of translational displacement and the 3 components of the rotational displacement
   REAL(ReKi)                                      :: XD   (6)                   ! Temp array holding 3 components of the translational velocity and the 3 components of the rotational (angular) velocity
   REAL(ReKi)                                      :: Ft   (6)                   ! Temp variable holding the 3 forces and 3 moments ; positive forces are in the direction of motion.
   INTEGER                                         :: J                          ! The number of the current node / element (-) 
   INTEGER                                         :: NMarkers                   ! The number of markers sent to this subroutine
   
   !-------------------------------------------------------------------------------------------------
   ! initialize the error status
   !-------------------------------------------------------------------------------------------------
     
   ErrStat = 0
   
   
   !-------------------------------------------------------------------------------------------------
   ! write previous output if we've passed the last calculation time 
   ! (this still has problems with going backwards in time, but will write the "best" solution 
   ! when called multiple times per time step.)
   !-------------------------------------------------------------------------------------------------
      
   IF ( ( CurrentTime > HD_Data%LastCalcTime ) ) THEN
   
      CALL HDOut_WriteOutputs(HD_Data%HDOut_Data, ErrStat )
      IF ( ErrStat > 0 ) RETURN  !less than zero in this case means that there aren't any outputs
      
      ErrStat = 0
      
   END IF     
      
      
   !-------------------------------------------------------------------------------------------------
   ! check that the HD_Loads and HD_Markers are allocated properly
   !-------------------------------------------------------------------------------------------------
   NMarkers = SIZE( HD_Markers%Substructure )
   
   IF ( .NOT. ALLOCATED( HD_Loads%Substructure ) ) THEN
   
      ALLOCATE( HD_Loads%Substructure( NMarkers ), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgAbort( ' Error allocating substructure load array in HydroDyn.', TrapErrors = .TRUE. )
         RETURN
      END IF
      
   ELSEIF ( SIZE( HD_Loads%Substructure ) /= NMarkers ) THEN
   
      CALL ProgAbort( ' Substructure load array is the wrong size in HydroDyn.', TrapErrors = .TRUE. )
      ErrStat = 1
      RETURN
   
   END IF
         
!   IF ( SIZE( HD_Markers%Substructure ) /= HD_Data%NElements ) THEN
!   
!      CALL ProgAbort( ' Substructure marker array is the wrong size in HydroDyn.', TrapErrors = .TRUE. )
!      ErrStat = 1
!      RETURN
!      
!   END IF
      
      
   !-------------------------------------------------------------------------------------------------
   ! get the loads from the appropriate module
   !-------------------------------------------------------------------------------------------------      
   
   SELECT CASE ( HD_Data%StrctType )
   
      CASE ( FloatPltfm_Type )
      
            ! These loads are point loads

!bjj question for jmj: you made a comment about floating platforms not requiring any elements ... how does that come into play here???

         DO J = 1,NMarkers !there should be only one marker at the moment
      
            X( 1:3) = HD_Markers%Substructure(J)%Position                               ! The 3 components of the translational        displacement (in m  )   of the current marker (platform reference or tower node)
            X( 4:6) = GetSmllRotAngs( HD_Markers%Substructure(J)%Orientation, ErrStat ) ! The 3 components of the rotational           displacement (in rad)   of the current marker (platform or tower element) relative to the inertial frame origin
            XD(1:3) = HD_Markers%Substructure(J)%TranslationVel                         ! The 3 components of the translational        velocity     (in m/s)   of the current marker (platform reference or tower node)
            XD(4:6) = HD_Markers%Substructure(J)%RotationVel                            ! The 3 components of the rotational (angular) velocity     (in rad/s) of the current marker (platform or tower element) relative to the inertial frame origin
                       

            IF ( J == NMarkers ) THEN !the last one is the point load
               CALL FltngPtfmLd( X, XD, CurrentTime, HD_Loads%Substructure(J)%AddedMass, Ft, HD_Data%Waves_data, &
                                        HD_Data%FltPtfm_data, ErrStat )   
            ELSE     ! these are loads per unit length
               CALL MorisonTwrLd ( J, HD_Data%Node(J)%D, HD_Data%Node(J)%CA, HD_Data%Node(J)%CD, &
                                 X, XD, CurrentTime, HD_Loads%Substructure(J)%AddedMass, Ft, HD_Data%Waves_data, ErrStat )
            END IF                                 
            
            HD_Loads%Substructure(J)%Force  = Ft(1:3)          ! the 3                                              components of the portion of the platform force                  (in N    ) acting at the platform reference    associated with everything but the added-mass effects
                                                               ! or the surge/xi (1), sway/yi (2), and heave/zi (3)-components of the portion of the tower    force  per unit length (in N/m  )        at the current tower element associated with everything but the added-mass effects; positive forces are in the direction of motion.
            HD_Loads%Substructure(J)%Moment = Ft(4:6)          ! the 3                                              components of the portion of the platform moment                 (in N-m  ) acting at the platform reference    associated with everything but the added-mass effects
                                                               ! or the roll/xi (1), pitch/yi (2), and yaw/zi   (3)-components of the portion of the tower    moment per unit length (in N-m/m) acting at the current tower element associated with everything but the added-mass effects; 
                                                               
                                                               ! Platform added mass matrix is (kg, kg-m, kg-m^2)
            IF ( ErrStat /= 0 ) RETURN
            
         END DO            
         
      CASE ( FixedBtm_Type )
      
            ! These loads are per unit length
            
         DO J = 1,NMarkers
         
            X( 1:3) = HD_Markers%Substructure(J)%Position                                 ! The 3 components of the translational        displacement (in m  )   of the current node
            X( 4:6) = GetSmllRotAngs( HD_Markers%Substructure(J)%Orientation, ErrStat )   ! The 3 components of the rotational           displacement (in rad)   of the current element relative to the inertial frame origin
            XD(1:3) = HD_Markers%Substructure(J)%TranslationVel                           ! The 3 components of the translational        velocity     (in m/s)   of the current node
            XD(4:6) = HD_Markers%Substructure(J)%RotationVel                              ! The 3 components of the rotational (angular) velocity     (in rad/s) of the current element relative to the inertial frame origin
         
            CALL MorisonTwrLd ( J, HD_Data%Node(J)%D, HD_Data%Node(J)%CA, HD_Data%Node(J)%CD, &
                                X, XD, CurrentTime, HD_Loads%Substructure(J)%AddedMass, Ft, HD_Data%Waves_data, ErrStat )
            HD_Loads%Substructure(J)%Force  = Ft(1:3)          ! The surge/xi (1), sway/yi (2), and heave/zi (3)-components of the portion of the tower force per unit length (in N/m)         at the current tower element associated with everything but the added-mass effects; positive forces are in the direction of motion.
            HD_Loads%Substructure(J)%Moment = Ft(4:6)          ! The roll/xi (1), pitch/yi (2), and yaw/zi (3)-components of the portion of the tower moment per unit length (in N-m/m) acting at the current tower element associated with everything but the added-mass effects; 
                                                               ! Fixed-Bottom added mass matrix is per unit length of element (kg/m, kg-m/m, kg-m^2/m)
            IF ( ErrStat /= 0 ) RETURN

         END DO            
      
      CASE DEFAULT
      
         CALL ProgAbort( ' Unknown support structure type in HD_CalculateLoads().', TrapErrors = .TRUE. )
         ErrStat = 1
         RETURN
   
   END SELECT 
   
   
   !-------------------------------------------------------------------------------------------------
   ! calculate output as requested
   !-------------------------------------------------------------------------------------------------
      
   CALL HDOut_CalcOutput ( CurrentTime, HD_Data%HDOut_Data, HD_Data%FltPtfm_data, HD_Data%Waves_data, ErrStat )   
   HD_Data%LastCalcTime = CurrentTime


END SUBROUTINE HD_CalculateLoads
!====================================================================================================
FUNCTION HD_CheckLin( HD_Data, ErrStat )
!   This public function checks that HydroDyn's data is properly set for linearization analysis.
!   The function returns .TRUE. if all linearization checks have been passed and .FALSE. if HydroDyn
!   will not be able to perform a linearization analysis based on the current settings.
!
! bjj: this function is still a bit of a hack, but I'm just tyring to get HydroDyn separate from FAST
! and ADAMS right now.
!----------------------------------------------------------------------------------------------------
      
      ! Passed variables
      
   TYPE(HD_DataType),         INTENT( INOUT )   :: HD_Data              ! the hydrodyn data 
   INTEGER,                   INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs  

   LOGICAL                                      :: HD_CheckLin
   
   
      ! Initialize the error status
   ErrStat = 0
   HD_CheckLin = .TRUE.
   
      !       
   IF ( HD_Data%Waves_data%WaveMod  /= 0  )  THEN     ! bjj: this should be done in the Waves module
      CALL ProgAbort ( ' HydroDyn can''t linearize a model with incident wave kinematics.  Set WaveMod to 0.', &
                         TrapErrors = .TRUE.  )
      ErrStat = -1
      HD_CheckLin = .FALSE.      
   ELSE IF ( HD_Data%StrctType == FloatPltfm_Type ) THEN
      
      IF ( FP_GetValue( "RdtnTMax", HD_Data%FltPtfm_data, ErrStat ) /= 0.0  ) THEN
         CALL ProgAbort ( ' HydroDyn can''t linearize a model with wave radiation damping.  Set RdtnTMax to 0.0.', &
                            TrapErrors = .TRUE. )
         ErrStat = -1
         HD_CheckLin = .FALSE.      
      END IF
   END IF

END FUNCTION HD_CheckLin
!====================================================================================================
SUBROUTINE HD_GetInput( HD_Data, FileName, Current_Data, Waves_InitData, FP_InitData, HDOut_InitData, ErrStat )
!     This private subroutine reads the input required for HydroDyn from the file whose name is an  
!     input parameter.
!----------------------------------------------------------------------------------------------------   

   
      ! Passed variables
   
   TYPE(HD_DataType),         INTENT( INOUT )   :: HD_Data              ! the hydrodyn data 
   CHARACTER(*),              INTENT( IN    )   :: FileName             ! Name of the HydroDyn input file   
   TYPE(Current_DataType),    INTENT(   OUT )   :: Current_Data         ! The current data from the input file
   TYPE(Waves_InitDataType),  INTENT(   OUT )   :: Waves_InitData       ! The waves data from the input file
   TYPE(FltPtfm_InitDataType),INTENT(   OUT )   :: FP_InitData          ! The floating platform data from the input file
   TYPE(HDO_InitDataType)    ,INTENT(   OUT )   :: HDOut_InitData       ! the values needed to initialize the the HydroDyn Output module   

   INTEGER,                   INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs  
   
  
      ! Local variables
      
   REAL(ReKi)                                   :: LAngAnch             ! Azimuth angle   of the current anchor   relative to the positive xi-axis of the inertial frame.
   REAL(ReKi)                                   :: LAngFair             ! Azimuth angle   of the current fairlead relative to the positive xt-axis of the platform.
   REAL(ReKi)                                   :: LDpthAnch            ! Depth           of the current anchor   relative to the origin           of the inertial frame.
   REAL(ReKi)                                   :: LDrftFair            ! Draft           of the current fairlead relative to the platform reference point.
   REAL(ReKi)                                   :: LRadAnch             ! Radial distance of the current anchor   relative to the origin           of the inertial frame.
   REAL(ReKi)                                   :: LRadFair             ! Radial distance of the current fairlead relative to the platform reference point.
         
   INTEGER                                      :: I                    ! generic integer for counting
   INTEGER                                      :: J                    ! generic integer for counting

   INTEGER                                      :: UnIn                 ! Unit number for the input file
      
   CHARACTER(80)                                :: Line                 ! String to temporarially hold value of read line   
   CHARACTER(1024)                              :: TmpPath              ! Temporary storage for relative path name
   CHARACTER(1024)                              :: TmpFmt             ! Temporary storage for format statement
     
   !-------------------------------------------------------------------------------------------------
   ! Open the file
   !-------------------------------------------------------------------------------------------------   
   CALL GetNewUnit( UnIn )   
   CALL OpenFInpFile( UnIn, TRIM(FileName), ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF


   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------
   
   CALL ReadCom( UnIn, FileName, 'HydroDyn input file header line 1', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF


   CALL ReadCom( UnIn, FileName, 'HydroDyn input file header line 2', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Features section
   !-------------------------------------------------------------------------------------------------
   
      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Features header', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF


      ! monopile or  floating platform


   CALL ReadVar( UnIn, FileName, HD_Data%StrctType, 'SupportStrct', 'Switch for structure type', ErrStat )

   IF (ErrStat /= 0) THEN
      CLOSE( UnIn )
      RETURN
   END IF

   IF ( HD_Data%StrctType /= FloatPltfm_Type .AND. HD_Data%StrctType /= FixedBtm_Type ) THEN
      CALL ProgAbort ( ' Error in file "'//TRIM(FileName)//&
                  '": Support structure must be either '//TRIM(Int2LStr(FloatPltfm_Type))// &
                  ' or '//TRIM(Int2LStr(FixedBtm_Type))//'.', .TRUE. )
      ErrStat = 1      
      CLOSE( UnIn )
      RETURN
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Environmental conditions section
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Environmental conditions header', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF


      ! WtrDens - Water density.
      
   CALL ReadVar ( UnIn, FileName, Waves_InitData%WtrDens, 'WtrDens', 'Water density', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF

   IF ( Waves_InitData%WtrDens < 0.0 )  THEN
      CALL ProgAbort ( ' WtrDens must not be negative.', .TRUE. )
      ErrStat = 1
      CLOSE( UnIn )
      RETURN
   END IF

      
      ! WtrDpth - Water depth   
      
   CALL ReadVar ( UnIn, FileName, Waves_InitData%WtrDpth, 'WtrDpth', 'Water depth', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF

   IF ( Waves_InitData%WtrDpth <= 0.0 )  THEN
      CALL ProgAbort ( ' WtrDpth must be greater than zero.', .TRUE. )
      ErrStat = 1
      CLOSE( UnIn )
      RETURN
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Data section for waves
   !-------------------------------------------------------------------------------------------------
      
      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Wave header', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF


!BJJ: verify that these tests can be performed (i.e., do we have the correct data to compare here?)
!   IF ( HD_Data%StrctType == FixedBtm_Type ) THEN
!
!!!JASON: WHAT LOADING DO WE APPLY TO THE FLEXIBLE PORTION OF THE TOWER EXTENDING BELOW THE SEABED?
! bjj: replace this :
!         IF ( ( TwrDraft - TwrRBHt ) < WtrDpth )  THEN   ! Print out a warning when the flexible portion of the support structure does not extend to the seabed.
! with this:
!         IF ( ( HydroConfig%Substructure%Position(3) ) < -WtrDpth )  THEN   ! Print out a warning when the flexible portion of the support structure does not extend to the seabed.
!            CALL ProgWarn( ' Hydrodynamic loading will only be applied to the flexible portion of the support structure.'// &
!                           ' Make sure that ( TwrDraft - TwrRBHt ) >= WtrDpth if you want hydrodynamic loading applied'// &
!                           ' along the entire submerged portion of the support structure. ')
!         ENDIF
!      
!   ELSE IF ( HD_Data%StrctType == FloatPltfm_Type ) THEN
!   
!      IF ( WtrDpth <= PtfmDraft  )  THEN
!         CALL ProgAbort ( ' WtrDpth must be greater than PtfmDraft.', .TRUE. )
!         CLOSE( UnIn )
!         RETURN
!      END IF
!         
!      IF ( FP_InitData%LineMod == 1 )  THEN  ! .TRUE if we have standard quasi-static mooring lines.
!         DO I = 1,FP_InitData%NumLines ! Loop through all mooring lines
!
!            IF ( WtrDpth < -FP_InitData%MooringLine(I)%LAnchzi )  THEN
!               CALL ProgAbort ( ' WtrDpth must not be less than LDpthAnch('//TRIM( Int2LStr( I ) )//').', .TRUE. )
!               ErrStat = 1
!               CLOSE( UnIn )
!               RETURN
!            END IF
!               
!         ENDDO             ! I - All mooring lines
!      END IF
!      
!   END IF


      
      ! WaveMod - Wave kinematics model switch.

   CALL ReadVar ( UnIn, FileName, Line, 'WaveMod', 'Wave kinematics model switch', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   
   Waves_InitData%WavePhase = 0.0
   Waves_InitData%WaveNDAmp = .FALSE.
   
   IF ( LEN_TRIM(Line) > 1 ) THEN
   
      CALL Conv2UC( Line )    ! Convert Line to upper case.
            
      IF ( Line(1:2) == '1P' )  THEN                     ! The user wants to specify the phase in place of a random phase

         READ (Line(3:),*,IOSTAT=ErrStat)  Waves_InitData%WavePhase
         IF ( ErrStat /= 0 ) THEN
            CALL CheckIOS ( ErrStat, FileName, 'WavePhase', NumType, .TRUE. )
            CLOSE( UnIn )
            RETURN
         END IF
         Waves_InitData%WaveMod   = 10                               ! Internally define WaveMod = 10 to mean regular waves with a specified (nonrandom) phase
         Waves_InitData%WavePhase = Waves_InitData%WavePhase*D2R     ! Convert the phase from degrees to radians

      ELSE                                               ! The user must have specified WaveMod incorrectly.

         ErrStat = 1

      ENDIF

   
   ELSE
   
      READ( Line, *, IOSTAT=ErrStat ) Waves_InitData%WaveMod
      
      IF ( ErrStat /= 0 ) THEN
         CALL CheckIOS ( ErrStat, FileName, 'WaveMod', NumType, .TRUE. )
         CLOSE( UnIn )
         RETURN
      END IF

   END IF ! LEN_TRIM(Line)


   IF ( ErrStat /= 0 .OR. Waves_InitData%WaveMod < 0 .OR. Waves_InitData%WaveMod > 3 ) THEN
   
      IF ( HD_Data%StrctType == FloatPltfm_Type ) THEN
      
         CALL ProgAbort ( ' WaveMod must be 0, 1, 1P#, 2, or 3.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      
      ELSE IF ( ErrStat /= 0 .OR. Waves_InitData%WaveMod /= 4 )  THEN
            
         CALL ProgAbort ( ' WaveMod must be 0, 1, 1P#, 2, 3, or 4.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
         
      END IF
      
   END IF          
      
      
      ! WaveStMod - Model switch for stretching incident wave kinematics to instantaneous free surface. 
      
   IF ( HD_Data%StrctType == FixedBtm_Type .AND. Waves_InitData%WaveMod > 0 ) THEN
   
      CALL ReadVar ( UnIn, FileName, Waves_InitData%WaveStMod, 'WaveStMod', &
         'Model switch for stretching incident wave kinematics to instantaneous free surface', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF        

      IF ( ( Waves_InitData%WaveStMod /= 0 ) .AND. ( Waves_InitData%WaveStMod /= 1 ) .AND. &
           ( Waves_InitData%WaveStMod /= 2 ) .AND. ( Waves_InitData%WaveStMod /= 3 ) )  THEN
         CALL ProgAbort ( ' WaveStMod must be 0, 1, 2, or 3.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF

      IF ( ( Waves_InitData%WaveStMod /= 3 ) .AND. ( Waves_InitData%WaveMod == 4 ) )  THEN
         CALL ProgAbort ( ' WaveStMod must be set to 3 when WaveMod is set to 4.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF
      
   ELSE !don't use this one
   
         ! NOTE: Do not read in WaveStMod for floating platforms since it is
         !       inconsistent to use stretching (which is a nonlinear correction) for
         !       the viscous drag term in Morison's equation while not accounting for
         !       stretching in the diffraction and radiation problems (according to
         !       Paul Sclavounos, there are such corrections).  Instead, the viscous
         !       drag term from Morison's equation is computed by integrating up to
         !       the MSL, regardless of the instantaneous free surface elevation.

      Waves_InitData%WaveStMod = 0

      CALL ReadCom ( UnIn, FileName, 'unused WaveStMod', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
         
   END IF

         
      ! WaveTMax - Analysis time for incident wave calculations.  
      
   IF ( Waves_InitData%WaveMod > 0 )  THEN   ! .TRUE if we have incident waves.
   
      CALL ReadVar ( UnIn, FileName, Waves_InitData%WaveTMax, 'WaveTMax', &
                              'Analysis time for incident wave calculations', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
         
!      IF ( WaveTMax < TMax )  THEN
!         CALL ProgAbort ( ' WaveTMax must not be less than TMax.', .TRUE. )    
!         ErrStat = 1
!         CLOSE( UnIn )
!         RETURN
!      END IF
        
   ELSE

      Waves_InitData%WaveTMax = 0.0
      
      CALL ReadCom ( UnIn, FileName, 'unused WaveTMax',  ErrStat )
   
      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

   END IF   
            
      
      ! WaveDT - Time step for incident wave calculations   
      
   IF ( Waves_InitData%WaveMod > 0 )  THEN   ! .TRUE if we have incident waves.
   
      CALL ReadVar ( UnIn, FileName, Waves_InitData%WaveDT, 'WaveDT', &
                         'Time step for incident wave calculations', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

      IF ( Waves_InitData%WaveDT <= 0.0 )  THEN
         CALL ProgAbort ( ' WaveDT must be greater than zero.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF
      
   ELSE
   
      Waves_InitData%WaveDT = 0.0
      
      CALL ReadCom ( UnIn, FileName, 'unused WaveDT', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

   END IF        

      
      ! WaveHs - Significant wave height    
      
   IF ( ( Waves_InitData%WaveMod /= 0 ) .AND. ( Waves_InitData%WaveMod /= 3 ) .AND. ( Waves_InitData%WaveMod /= 4 ) ) THEN   ! .TRUE. (when WaveMod = 1, 2, or 10) if we have plane progressive (regular) or JONSWAP/Pierson-Moskowitz spectrum (irregular) waves, but not GH Bladed wave data.

      CALL ReadVar ( UnIn, FileName, Waves_InitData%WaveHs, 'WaveHs', 'Significant wave height', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

      IF ( Waves_InitData%WaveHs <= 0.0 )  THEN
         CALL ProgAbort ( ' WaveHs must be greater than zero.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF

   ELSE
   
      Waves_InitData%WaveHs = 0.0
      
      CALL ReadCom ( UnIn, FileName, 'unused WaveHs', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
   
   END IF      

      
      ! WaveTp - Peak spectral period.

   IF ( ( Waves_InitData%WaveMod /= 0 ) .AND. ( Waves_InitData%WaveMod /= 3 ) .AND. ( Waves_InitData%WaveMod /= 4 ) ) THEN   ! .TRUE. (when WaveMod = 1, 2, or 10) if we have plane progressive (regular) or JONSWAP/Pierson-Moskowitz spectrum (irregular) waves, but not GH Bladed wave data.

      CALL ReadVar ( UnIn, FileName, Waves_InitData%WaveTp, 'WaveTp', 'Peak spectral period', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

      IF ( Waves_InitData%WaveTp <= 0.0 )  THEN
         CALL ProgAbort ( ' WaveTp must be greater than zero.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF

   ELSE
   
      Waves_InitData%WaveTp = 0.0
      
      CALL ReadCom ( UnIn, FileName, 'unused WaveTp', ErrStat )
   
      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

   END IF      

      
      ! WavePkShp - Peak shape parameter.
      
   IF ( Waves_InitData%WaveMod == 2 ) THEN   ! .TRUE if we have JONSWAP/Pierson-Moskowitz spectrum (irregular) waves, but not GH Bladed wave data.


      CALL ReadVar ( UnIn, FileName, Line, 'WavePkShp', 'Peak shape parameter', ErrStat )
      
      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
               
      CALL Conv2UC( Line )    ! Convert Line to upper case.

      IF ( TRIM(Line) == 'DEFAULT' )  THEN   ! .TRUE. when one wants to use the default value of the peak shape parameter, conditioned on significant wave height and peak spectral period.

         Waves_InitData%WavePkShp = WavePkShpDefault ( Waves_InitData%WaveHs, Waves_InitData%WaveTp )

      ELSE                                   ! The input must have been specified numerically.

         READ (Line,*,IOSTAT=ErrStat)  Waves_InitData%WavePkShp
         
         IF ( ErrStat /= 0 ) THEN
            CALL CheckIOS ( ErrStat, FileName, 'WavePkShp', NumType, .TRUE. )
            CLOSE( UnIn )
            RETURN
         END IF

         IF ( ( Waves_InitData%WavePkShp < 1.0 ) .OR. ( Waves_InitData%WavePkShp > 7.0 ) )  THEN        
            CALL ProgAbort ( ' WavePkShp must be greater than or equal to 1 and less than or equal to 7.', .TRUE. )
            ErrStat = 1
            CLOSE( UnIn )
            RETURN
         END IF

      END IF

   ELSE

      Waves_InitData%WavePkShp = 1.0
      
      CALL ReadCom ( UnIn, FileName, 'unused WavePkShp', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
   
   END IF      

      
      ! WaveDir - Wave heading direction.  
      
   IF ( ( Waves_InitData%WaveMod > 0 ) .AND. ( Waves_InitData%WaveMod /= 4 ) )  THEN   ! .TRUE if we have incident waves, but not GH Bladed wave data.

      CALL ReadVar ( UnIn, FileName, Waves_InitData%WaveDir, 'WaveDir', 'Wave heading direction', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

      IF ( ( Waves_InitData%WaveDir <= -180.0 ) .OR. ( Waves_InitData%WaveDir > 180.0 ) )  THEN
         CALL ProgAbort ( ' WaveDir must be greater than -180 and less than or equal to 180.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF

   ELSE

      Waves_InitData%WaveDir = 0.0
      
      CALL ReadCom ( UnIn, FileName, 'unused WaveDir', ErrStat)
      
      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
   
   END IF      

      
      ! WaveSeed(1), !WaveSeed(2)

   IF ( ( Waves_InitData%WaveMod > 0 ) .AND. ( Waves_InitData%WaveMod /= 4 ) .AND. ( Waves_InitData%WaveMod /= 10 ) )  THEN   !.TRUE. for plane progressive (regular) with random phase or irregular wave 

      DO I = 1,2
      
         WRITE(Line,'(I2)') I
      
         CALL ReadVar ( UnIn, FileName, Waves_InitData%WaveSeed(I), 'WaveSeed('//TRIM(Line)//')', &
                                       'Random seed #'//TRIM(Line), ErrStat )
      
         IF ( ErrStat /= 0 ) THEN
            CLOSE( UnIn )
            RETURN
         END IF
      
      END DO !I


   ELSE

      DO I = 1,2

         Waves_InitData%WaveSeed(I) = 0

         WRITE(Line,'(I2)') I
      
         CALL ReadCom ( UnIn, FileName, 'unused WaveSeed('//TRIM(Line)//')' , ErrStat)

         IF ( ErrStat /= 0 ) THEN
            CLOSE( UnIn )
            RETURN
         END IF
      
      END DO !I
   
   END IF      
      
      
      
      ! GHWvFile   
      
   IF ( Waves_InitData%WaveMod == 4 ) THEN      ! .TRUE if we are to use GH Bladed wave data.

      CALL ReadVar ( UnIn, FileName, Waves_InitData%GHWvFile, 'GHWvFile', &
                                     'Root name of GH Bladed files containing wave data', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

      IF ( LEN_TRIM( Waves_InitData%GHWvFile ) == 0 )  THEN      
         CALL ProgAbort ( ' GHWvFile must not be an empty string.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF

   ELSE !don't use this one

      Waves_InitData%GHWvFile = ""
      
      CALL ReadCom ( UnIn, FileName, 'unused GHWvFile', ErrStat )
      
      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

   END IF


      ! NWaveElev

   CALL ReadVar ( UnIn, FileName, Waves_InitData%NWaveElev, 'NWaveElev', &
                                  'Number of points where the incident wave elevations can be output', ErrStat )
      
   IF ( ErrStat /= 0 ) THEN
   
      CLOSE( UnIn )
      RETURN
      
   END IF
      
   IF ( Waves_InitData%NWaveElev < 0 ) THEN
       
      CALL ProgAbort ( ' NWaveElev must not be negative.', TrapErrors = .TRUE. )
      CLOSE( UnIn )
      RETURN
      
   ELSE
      
         ! allocate space for the output location arrays: 
      
      ALLOCATE ( Waves_InitData%WaveElevxi(Waves_InitData%NWaveElev), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgAbort ( ' Error allocating space for WaveElevxi array.', TrapErrors = .TRUE. )
         CLOSE( UnIn )
         RETURN
      END IF
      
      ALLOCATE ( Waves_InitData%WaveElevyi(Waves_InitData%NWaveElev), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgAbort ( ' Error allocating space for WaveElevyi array.', TrapErrors = .TRUE. )
         CLOSE( UnIn )
         RETURN
      END IF
      
   END IF


      ! WaveElevxi
         
   IF ( Waves_InitData%NWaveElev > 0 ) THEN
            
      CALL ReadAry ( UnIn, FileName, Waves_InitData%WaveElevxi, Waves_InitData%NWaveElev, 'WaveElevxi', &
                           'List of xi-coordinates for points where the incident wave elevations can be output', ErrStat )         
      
   ELSE
         ! this adds a line to the echo file, even if NWaveElev < 1
         
      CALL ReadCom ( UnIn, FileName, 'unused WaveElevxi', ErrStat )
      
   END IF
         
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
      
      
      ! WaveElevyi

   IF ( Waves_InitData%NWaveElev > 0 ) THEN

      CALL ReadAry ( UnIn, FileName, Waves_InitData%WaveElevyi, Waves_InitData%NWaveElev, 'WaveElevyi', &
                           'List of yi-coordinates for points where the incident wave elevations can be output', ErrStat )         

   ELSE
   
      CALL ReadCom ( UnIn, FileName, 'unused WaveElevyi', ErrStat )
      
   END IF
      
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Data section for current
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Current header', ErrStat )
   
   IF (ErrStat /= 0) THEN
      CLOSE( UnIn )
      RETURN
   END IF


      ! CurrMod - Current profile model switch
      
   CALL ReadVar ( UnIn, FileName, Current_Data%CurrMod, 'CurrMod', 'Current profile model switch', ErrStat )

   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF

   IF ( ( Current_Data%CurrMod /= 0 ) .AND. ( Current_Data%CurrMod /= 1 ) .AND. ( Current_Data%CurrMod /= 2 ) )  THEN
      CALL ProgAbort ( ' CurrMod must be 0, 1, or 2.', .TRUE. )
      ErrStat = 1
      CLOSE( UnIn )
      RETURN
   END IF

   IF ( ( Current_Data%CurrMod /= 0 ) .AND. ( Waves_InitData%WaveMod == 4 ) )  THEN
      CALL ProgAbort ( ' CurrMod must be set to 0 when WaveMod is set to 4.', .TRUE. )
      ErrStat = 1
      CLOSE( UnIn )
      RETURN
   END IF


      ! CurrSSV0 - Sub-surface current velocity at still water level
      
   IF ( Current_Data%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.
    
      CALL ReadVar ( UnIn, FileName, Current_Data%CurrSSV0, 'CurrSSV0', 'Sub-surface current velocity at still water level', &
                         ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

      IF ( Current_Data%CurrSSV0 < 0.0 )  THEN
         CALL ProgAbort ( ' CurrSSV0 must not be less than zero.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN         
      END IF  
  
   ELSE

      Current_Data%CurrSSV0 = 0.0
      
      CALL ReadCom ( UnIn, FileName, 'unused CurrSSV0', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
   
   END IF
   
     
      ! CurrSSDir - Sub-surface current heading direction
   
   IF ( Current_Data%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.
  
  
      CALL ReadVar ( UnIn, FileName, Line, 'CurrSSDir', 'Sub-surface current heading direction', ErrStat )
      
      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
            
            
      CALL Conv2UC( Line )    ! Convert Line to upper case.

      IF ( TRIM(Line) == 'DEFAULT' )  THEN   ! .TRUE. when one wants to use the default value of codirectionality between sub-surface current and incident wave propogation heading directions.

         IF ( Waves_InitData%WaveMod == 0 ) THEN
            CALL ProgAbort ( ' CurrSSDir must not be set to ''DEFAULT'' when WaveMod is set to 0.', .TRUE. )
            ErrStat = 1
            CLOSE( UnIn )
            RETURN         
         END IF  

         Current_Data%CurrSSDir = Waves_InitData%WaveDir

      ELSE                                   ! The input must have been specified numerically.

         READ (Line,*,IOSTAT=ErrStat)  Current_Data%CurrSSDir
         CALL CheckIOS ( ErrStat, FileName, 'CurrSSDir', NumType, .TRUE. )

         IF ( ErrStat /= 0 ) THEN
            CLOSE( UnIn )
            RETURN
         END IF

         IF ( ( Current_Data%CurrSSDir <= -180.0 ) .OR. ( Current_Data%CurrSSDir > 180.0 ) )  THEN 
            CALL ProgAbort ( ' CurrSSDir must be greater than -180 and less than or equal to 180.', .TRUE. ) 
            ErrStat = 1
            CLOSE( UnIn )
            RETURN         
         END IF  

      END IF
  
  
   ELSE

      Current_Data%CurrSSDir = 0.0
      
      CALL ReadCom ( UnIn, FileName, 'unused CurrSSDir', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
   
   END IF
   
   
      ! CurrNSRef - Near-surface current reference depth.
   
   IF ( Current_Data%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.
  
      CALL ReadVar ( UnIn, FileName, Current_Data%CurrNSRef, 'CurrNSRef', 'Near-surface current reference depth', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

      IF ( Current_Data%CurrNSRef <= 0.0 ) THEN
         CALL ProgAbort ( ' CurrNSRef must be greater than zero.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF

  
   ELSE
   
      Current_Data%CurrNSRef = 0.0

      CALL ReadCom ( UnIn, FileName, 'unused CurrNSRef', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
   
   END IF
   
   
      ! CurrNSV0 - Near-surface current velocity at still water level.
   
   IF ( Current_Data%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.
  
      CALL ReadVar ( UnIn, FileName, Current_Data%CurrNSV0, 'CurrNSV0', 'Near-surface current velocity at still water level', &
                           ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
  
      IF ( Current_Data%CurrNSV0 < 0.0 ) THEN
         CALL ProgAbort ( ' CurrNSV0 must not be less than zero.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF
    
   ELSE
   
      Current_Data%CurrNSV0 = 0.0

      CALL ReadCom ( UnIn, FileName, 'unused CurrNSV0' , ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
   
   END IF

   
      ! CurrNSDir - Near-surface current heading direction.

   IF ( Current_Data%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.
  
      CALL ReadVar ( UnIn, FileName, Current_Data%CurrNSDir, 'CurrNSDir', 'Near-surface current heading direction', ErrStat )
  
      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
  
      IF ( ( Current_Data%CurrNSDir <= -180.0 ) .OR. ( Current_Data%CurrNSDir > 180.0 ) )  THEN
         CALL ProgAbort ( ' CurrNSDir must be greater than -180 and less than or equal to 180.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF
           
   ELSE
   
      Current_Data%CurrNSDir = 0.0
      
      CALL ReadCom ( UnIn, FileName, 'unused CurrNSDir', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
   
   END IF
   
   
      ! CurrDIV - Depth-independent current velocity.

   IF ( Current_Data%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.
  
      CALL ReadVar ( UnIn, FileName, Current_Data%CurrDIV, 'CurrDIV', 'Depth-independent current velocity', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

      IF ( Current_Data%CurrDIV < 0.0 ) THEN
         CALL ProgAbort ( ' CurrDIV must not be less than zero.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF
         
   ELSE

      Current_Data%CurrDIV = 0.0
      
      CALL ReadCom ( UnIn, FileName, 'unused CurrDIV'  , ErrStat )
      
      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
   
   END IF
   
   
      ! CurrDIDir - Depth-independent current heading direction.
   
   IF ( Current_Data%CurrMod == 1 )  THEN  ! .TRUE if we have standard current.
  
  
      CALL ReadVar ( UnIn, FileName, Current_Data%CurrDIDir, 'CurrDIDir', 'Depth-independent current heading direction', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

      IF ( ( Current_Data%CurrDIDir <= -180.0 ) .OR. ( Current_Data%CurrDIDir > 180.0 ) ) THEN
         CALL ProgAbort ( ' CurrDIDir must be greater than -180 and less than or equal to 180.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF
  
  
   ELSE

      Current_Data%CurrDIDir = 0.0
      
      CALL ReadCom ( UnIn, FileName, 'unused CurrDIDir', ErrStat )
      
      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
   
   END IF


   !-------------------------------------------------------------------------------------------------
   ! Data section to define the discretization 
   !-------------------------------------------------------------------------------------------------

   CALL HD_GetDiscretization(UnIn, FileName, HD_Data, HDOut_InitData%NumMemberNodes, ErrStat)
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF


!   !-------------------------------------------------------------------------------------------------
!   ! Data section for monopile tower
!   !-------------------------------------------------------------------------------------------------
!
!      ! Header
!      
!   CALL ReadCom( UnIn, FileName, 'Monopile header', ErrStat )
!   
!   IF ( ErrStat /= 0 ) THEN
!      CLOSE( UnIn )
!      RETURN
!   END IF
!
!
!
!      ! MPLdMod - Tower loading model switch
!
!   IF ( HD_Data%StrctType == FixedBtm_Type ) THEN
!      
!      CALL ReadVar ( UnIn, FileName, MPLdMod, 'MPLdMod', 'Tower loading model switch', ErrStat )
!      
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!      
!      IF ( ( MPLdMod /= 0 ) .AND. ( MPLdMod /= 1 ) .AND. ( MPLdMod /= 2 ) )  THEN         
!         CALL ProgAbort ( ' MPLdMod must be 0, 1, or 2.', .TRUE. )
!         CLOSE( UnIn )
!         ErrStat = 1
!         RETURN
!      END IF      
!      
!   ELSE
!   
!      MPLdMod = 0
!      CALL ReadCom ( UnIn, FileName, 'unused MPLdMod', ErrStat )
!
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!      
!   END IF
!
!
!      ! MPNodes
!      
!!   IF ( HD_Data%StrctType == FixedBtm_Type ) THEN
!!
!!      CALL ReadVar ( UnIn, FileName, MPNodes, 'MPNodes', 'Number of tower nodes', ErrStat )
!!      IF ( ErrStat /= 0 ) THEN
!!         CLOSE( UnIn )
!!         RETURN
!!      END IF
!!
!!      IF ( MPNodes < 1 ) THEN
!!      
!!         CALL ProgAbort( 'MPNodes must be greater than 0.', TrapErrors = .TRUE. )
!!         ErrStat = 1
!!         CLOSE( UnIn )
!!         RETURN
!!      
!!      END IF
!!
!!   ELSE
!         
!      MPNodes = 0
!      CALL ReadCom ( UnIn, FileName, 'unused MPNodes', ErrStat )
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!      
!!   END IF   
!      
!      
!      ! MPzi(MPNodes)
!
!!   IF ( MPNodes > 0  ) THEN
!!   
!!      ALLOCATE( MPzi(MPNodes), STAT=ErrStat )
!!      
!!      IF ( ErrStat /= 0 ) THEN
!!         CALL ProgAbort( ' Error allocating space for array MPzi.', TrapErrors = .TRUE. )
!!         CLOSE( UnIn )
!!         RETURN
!!      END IF
!!      
!!      
!!      CALL ReadAry ( UnIn, FileName, MPzi, MPNodes, 'MPzi', 'Tower node locations', ErrStat )
!!
!!   ELSE
!
!      CALL ReadCom ( UnIn, FileName, 'unused MPzi', ErrStat )
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!!   END IF
!         
!      ! MPDiam - Tower diameter in Morison's equation
!
!   IF ( MPLdMod == 1 ) THEN      ! we will be using the built-in Morison's equation         
!   
!      CALL ReadVar ( UnIn, FileName, MPDiam, 'MPDiam', 'Tower diameter in Morison''s equation', ErrStat )
!
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!      IF ( MPDiam < 0.0 ) THEN
!         CALL ProgAbort ( ' MPDiam must not be negative.', .TRUE. )
!         ErrStat = 1
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!   ELSE
!   
!      MPDiam = 0.0
!      
!      CALL ReadCom ( UnIn, FileName, 'unused MPDiam', ErrStat )
!
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!   
!   END IF
!
!
!      ! MPCA - Normalized hydrodynamic added mass coefficient in Morison's equation
!
!   IF ( MPLdMod == 1 ) THEN      ! we will be using the built-in Morison's equation         
!
!      CALL ReadVar ( UnIn, FileName, MPCA, 'MPCA', &
!                        'Normalized hydrodynamic added mass coefficient in Morison''s equation', ErrStat )
!
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!      IF ( MPCA < 0.0 )  THEN
!         CALL ProgAbort ( ' MPCA must not be negative.', .TRUE. )
!         ErrStat = 1
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!   ELSE
!
!      MPCA = 0.0
!      CALL ReadCom ( UnIn, FileName, 'unused MPCA',   ErrStat )
!
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!   
!   END IF
!
!
!      ! MPCD - Normalized hydrodynamic viscous drag coefficient in Morison's equation
!
!   IF ( MPLdMod == 1 ) THEN      ! we will be using the built-in Morison's equation         
!
!      CALL ReadVar ( UnIn, FileName, MPCD, 'MPCD', &
!                        'Normalized hydrodynamic viscous drag coefficient in Morison''s equation', ErrStat )
!
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!      IF ( MPCD < 0.0 ) THEN
!         CALL ProgAbort ( ' MPCD must not be negative.', .TRUE. )
!         ErrStat = 1
!         CLOSE( UnIn )
!         RETURN
!      END IF        
!
!
!   ELSE  ! We must not be using the built-in Morison's equation, so skip these inputs.
!     
!      MPCD = 0.0
!      CALL ReadCom ( UnIn, FileName, 'unused MPCD',   ErrStat )
!      
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!   END IF


   !-------------------------------------------------------------------------------------------------
   ! Data section for floating platform
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Floating platform header', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF
   !!!!!!!!!!!!!!!!!!!!!!!!
   
!   IF ( HD_Data%StrctType == FloatPltfm_Type ) THEN ! (should this be done in FAST?)
!
!instead of the following 2 checks, now check that the HydroConfig marker for platform is at 0,0,0
!      IF ( TwrDraft > 0.0 ) THEN
!         CALL ProgAbort ( ' TwrDraft must be less than or equal to zero when PtfmLdMod is set to "'//TRIM(Line)//'".', .TRUE. )  ! Do not allow the combination of tower hydrodynamics using Morison's equation and platform hydrodynamics using the true form of the using the true form of the hydrodynamics equations since the true equations require that the shape of the platform does not change above the MSL (platform reference point)--Consider the linear hydrostatic restoring matrix, for example.
!         ErrStat = 1
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!      IF ( PtfmRef /= 0.0 ) THEN
!         CALL ProgAbort ( ' PtfmRef must be zero when PtfmLdMod is set to "'//TRIM(Line)//'".', .TRUE. )
!         ErrStat = 1
!         CLOSE( UnIn )
!         RETURN
!      END IF
!bjj: like this:
!      IF ( HydroConfig%Substructure%Position  /= ( 0,0,0 ) .OR.
!           HydroConfig%Substructure%Orientation /= eye(3) ) THEN
!      
!         CALL ProgAbort ( ' HydroConfig%Substructure%Position must be zero and Orientation must be the identity matrix when PtfmLdMod is set to "'//TRIM(Line)//'".', .TRUE. )
!           
!      END IF
!
!   END IF


      ! WAMITFile - Root name of WAMIT output files

   IF ( HD_Data%StrctType == FloatPltfm_Type ) THEN
   
      CALL ReadVar ( UnIn, FileName, FP_InitData%WAMITFile, 'WAMITFile', 'Root name of WAMIT output files', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

      IF ( LEN_TRIM( FP_InitData%WAMITFile ) == 0 ) THEN
         CALL ProgAbort ( ' WAMITFile must not be an empty string.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF
      
      
         ! if this is a relative path, let's make it relative to the location of the main input file
         
      IF ( PathIsRelative( FP_InitData%WAMITFile ) ) THEN      
         CALL GetPath( FileName, TmpPath ) 
         FP_InitData%WAMITFile = TRIM(TmpPath)//TRIM(FP_InitData%WAMITFile)
      END IF
         

   ELSE         
   
      FP_InitData%WAMITFile = ""
   
      CALL ReadCom( UnIn, FileName, 'unused WAMITFile', ErrStat )
   
      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
      
   END IF


      ! WAMITFile - Root name of WAMIT output files

   IF ( HD_Data%StrctType == FloatPltfm_Type ) THEN
   
      CALL ReadVar ( UnIn, FileName, FP_InitData%WAMITULEN, 'WAMITULEN', 'WAMIT characteristic body length scale', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF                     

   ELSE         
   
      FP_InitData%WAMITULEN = 1.0
   
      CALL ReadCom( UnIn, FileName, 'unused WAMITULEN', ErrStat )
   
      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
      
   END IF



      ! PtfmVol0 - Displaced volume of water when the platform is in its undisplaced position

   IF ( HD_Data%StrctType == FloatPltfm_Type ) THEN
   
      CALL ReadVar ( UnIn, FileName, FP_InitData%PtfmVol0, 'PtfmVol0', &
         'Displaced volume of water when the platform is in its undisplaced position', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF


      IF ( FP_InitData%PtfmVol0 < 0.0 ) THEN
         CALL ProgAbort ( ' PtfmVol0 must not be negative.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF
      
   ELSE         
   
      FP_InitData%PtfmVol0 = 0.0
      
      CALL ReadCom( UnIn, FileName, 'unused PtfmVol0', ErrStat )
   
      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
      
   END IF


!      ! PtfmNodes - Number of platform nodes used in calculation of viscous drag term from Morison's equation
!
!   IF ( HD_Data%StrctType == FloatPltfm_Type ) THEN
!   
!      CALL ReadVar ( UnIn, FileName, PtfmNodes, 'PtfmNodes', &
!         'Number of platform nodes used in calculation of viscous drag term from Morison''s equation', ErrStat )
!
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!      IF ( PtfmNodes < 0 ) THEN
!         CALL ProgAbort ( ' PtfmNodes must not be less than 0.', .TRUE. )
!         ErrStat = 1
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!   ELSE         
!   
!      PtfmNodes = 0
!      
!      CALL ReadCom( UnIn, FileName, 'unused PtfmNodes', ErrStat )
!   
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!      
!   END IF
!
!
!      ! PtfmDraft - Effective platform draft in calculation of viscous drag term from Morison's equation
!
!   IF ( HD_Data%StrctType == FloatPltfm_Type ) THEN
!   
!      CALL ReadVar ( UnIn, FileName, PtfmDraft, 'PtfmDraft', &
!         'Effective platform draft in calculation of viscous drag term from Morison''s equation', ErrStat )
!
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!      IF ( PtfmDraft < 0.0 ) THEN
!         CALL ProgAbort ( ' PtfmDraft must not be negative.', .TRUE. )
!         ErrStat = 1
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!   ELSE         
!   
!      PtfmDraft = 0.0
!      
!      CALL ReadCom( UnIn, FileName, 'unused PtfmDraft', ErrStat )
!   
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!      
!   END IF
!
!
!      ! PtfmDiam - Effective platform diameter in calculation of viscous drag term from Morison's equation
!
!   IF ( HD_Data%StrctType == FloatPltfm_Type ) THEN
!   
!      CALL ReadVar ( UnIn, FileName, FP_InitData%PtfmDiam, 'PtfmDiam', &
!         'Effective platform diameter in calculation of viscous drag term from Morison''s equation', ErrStat )
!
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!      IF ( FP_InitData%PtfmDiam < 0.0 ) THEN
!         CALL ProgAbort ( ' PtfmDiam must not be negative.', .TRUE. )
!         ErrStat = 1
!         CLOSE( UnIn )
!         RETURN
!      END IF
!   
!   ELSE         
!   
!      FP_InitData%PtfmDiam = 0.0
!            
!      CALL ReadCom( UnIn, FileName, 'unused PtfmDiam', ErrStat )
!   
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!      
!   END IF
!
!
!      ! PtfmCD - Effective platform normalized hydrodynamic viscous drag coefficient in calculation of viscous drag term from Morison's equation
!
!   IF ( HD_Data%StrctType == FloatPltfm_Type ) THEN
!   
!      CALL ReadVar ( UnIn, FileName, FP_InitData%PtfmCD, 'PtfmCD', &
!         'Effective platform normalized hydrodynamic viscous drag coefficient in Morison''s equation', ErrStat )
!
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!      IF ( FP_InitData%PtfmCD < 0.0 ) THEN
!         CALL ProgAbort ( ' PtfmCD must not be negative.', .TRUE. )
!         ErrStat = 1
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!   ELSE         
!   
!      FP_InitData%PtfmCD = 0.0
!      
!      CALL ReadCom( UnIn, FileName, 'unused PtfmCD', ErrStat )
!   
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!      
!   END IF


      ! RdtnTMax - Analysis time for wave radiation kernel calculations
      ! NOTE: Use RdtnTMax = 0.0 to eliminate wave radiation damping

   IF ( HD_Data%StrctType == FloatPltfm_Type ) THEN
   
      CALL ReadVar ( UnIn, FileName, FP_InitData%RdtnTMax, 'RdtnTMax', &
                                    'Analysis time for wave radiation kernel calculations', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

      IF ( FP_InitData%RdtnTMax < 0.0 ) THEN
         CALL ProgAbort ( ' RdtnTMax must not be negative.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF

   ELSE         
   
      FP_InitData%RdtnTMax = 0.0
   
      CALL ReadCom( UnIn, FileName, 'unused RdtnTMax', ErrStat )
   
      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
      
   END IF


      ! RdtnDT - Time step for wave radiation kernel calculations

   IF ( HD_Data%StrctType == FloatPltfm_Type ) THEN
   
      CALL ReadVar ( UnIn, FileName, FP_InitData%RdtnDT, 'RdtnDT', 'Time step for wave radiation kernel calculations', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

      IF ( FP_InitData%RdtnDT <= 0.0 ) THEN
         CALL ProgAbort ( ' RdtnDT must be greater than zero.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN
      END IF

   ELSE         
   
      FP_InitData%RdtnDT = 0.0
   
      CALL ReadCom( UnIn, FileName, 'unused RdtnDT', ErrStat )
   
      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF
      
   END IF

!bjj: should we add this?
!test for numerical stability
!      IF ( FP_InitData%RdtnDT <= FP_InitData%RdtnTMax*EPSILON(FP_InitData%RdtnDT) )  THEN  ! Test RdtnDT and RdtnTMax to ensure numerical stability -- HINT: see the use of OnePlusEps." 
!         CALL ProgAbort ( ' RdtnDT must be greater than '//TRIM ( Num2LStr( RdtnTMax*EPSILON(RdtnDT) ) )//' seconds.', .TRUE. ) 
!         ErrStat = 1
!         CLOSE( UnIn )
!         RETURN   
!      END IF
   
   
   !-------------------------------------------------------------------------------------------------
   ! Data section for mooring lines (ONLY valid for StrctType == FloatPltfm_Type)
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Mooring lines header', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF


         ! LineMod - Mooring line model switch.
         
   IF ( HD_Data%StrctType == FloatPltfm_Type )  THEN

      CALL ReadVar ( UnIn, FileName, FP_InitData%LineMod, 'LineMod', 'Mooring line model switch', ErrStat )
      
      IF ( ErrStat /= 0 )  THEN
         CLOSE( UnIn )
         RETURN
      END IF
      

      IF ( ( FP_InitData%LineMod /= 1 ) .AND. ( FP_InitData%LineMod /= 2 ) ) THEN
         CALL ProgAbort ( ' LineMod must be 1 or 2.', .TRUE. )
         CLOSE( UnIn )
         ErrStat = 1
         RETURN
      END IF

   ELSE 
   
      FP_InitData%LineMod = 0 ! Set LineMod to zero for "none".

      CALL ReadCom ( UnIn, FileName, 'unused LineMod', ErrStat )

      IF ( ErrStat /= 0 )  THEN
         CLOSE( UnIn )
         RETURN
      END IF
      

   END IF


      ! NumLines - Number of mooring lines
      
   IF ( FP_InitData%LineMod == 1 ) THEN

      CALL ReadVar ( UnIn, FileName, FP_InitData%NumLines, 'NumLines', 'Number of mooring lines', ErrStat )

      IF ( ErrStat /= 0 ) THEN
         CLOSE( UnIn )
         RETURN
      END IF

      IF ( FP_InitData%NumLines < 0 ) THEN   
         CALL ProgAbort ( ' NumLines must not be less than zero.', .TRUE. )
         ErrStat = 1
         CLOSE( UnIn )
         RETURN   
      END IF
      
         ! Allocate the input mooring line array.
         ! NOTE: We must ALLOCATE these arrays even when LineMod <> 1 because the
         !       arrays are passed into the InitFltngPtfm() routine.
     
      ALLOCATE ( FP_InitData%MooringLine ( FP_InitData%NumLines ) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort ( ' Error allocating memory for the MooringLine array.', .TRUE. )
         CLOSE( UnIn )
         RETURN
      END IF
         
   
   ELSE
   
      FP_InitData%NumLines = 0
      
      CALL ReadCom ( UnIn, FileName, 'unused NumLines', ErrStat )

      IF ( ErrStat /= 0 )  THEN
         CLOSE( UnIn )
         RETURN
      END IF
      
   END IF
    
   
      ! Skip the header/comment lines.

   CALL ReadCom ( UnIn, FileName, 'mooring line parameter names', ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CLOSE( UnIn )
      RETURN
   END IF

   
   CALL ReadCom ( UnIn, FileName, 'mooring line parameter units', ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CLOSE( UnIn )
      RETURN
   END IF
   

   IF ( FP_InitData%LineMod == 1 )  THEN  ! .TRUE if we have standard quasi-static mooring lines.

!!bjj: why are we printing this twice? (once in the comments and now again here) probably because it looks nicer to have the headers printed exactly on top of the columns   
!
      IF ( Echo )  THEN

         WRITE (UnEc, "( '   LRadAnch     LAngAnch     LDpthAnch    LRadFair     LAngFair     LDrftFair"// &
                      "    LUnstrLen    LDiam        LMassDen     LEAStff      LSeabedCD    LTenTol     LineNodes   LSNodes:' )" )
         WRITE (UnEc, "( '   --------     --------     ---------    --------     --------     ---------"// &
                      "    ---------    -----        --------     -------      ---------    -------     ---------   --------' )" )
!         Frmt = '( I5, 1X, 12( 2X, '//TRIM( OutFmt )//'), 4X, I5 )'

      ENDIF


         ! Read in the mooring line data.
         ! NOTE: Store the x-, y-, and z-coordinates of each anchor and fairlead,
         !       instead of the radius, angle, and depth/draft.

      
      DO I = 1,FP_InitData%NumLines ! Loop through all mooring lines

         READ (UnIn,*,IOSTAT=ErrStat) LRadAnch                            , LAngAnch                            , &
                                      LDpthAnch                           , LRadFair                            , &
                                      LAngFair                            , LDrftFair                           , &
                                      FP_InitData%MooringLine(I)%LUnstrLen, FP_InitData%MooringLine(I)%LDiam    , &
                                      FP_InitData%MooringLine(I)%LMassDen , FP_InitData%MooringLine(I)%LEAStff  , &
                                      FP_InitData%MooringLine(I)%LSeabedCD, FP_InitData%MooringLine(I)%LTenTol  , &
                                      FP_InitData%MooringLine(I)%LineNodes 
   !bjj: we're now going to assume LSNodes are equally distributed instead of allowing the user to input them like this:                                      
                                                                          !,                                       &
!                                       ( TmpAry(J), J=1, MIN(MaxLineNodes,FP_InitData%MooringLine(I)%LineNodes) )


         CALL CheckIOS ( ErrStat, FileName, 'mooring line '//TRIM( Int2LStr( I ) )//' properties', NumType, .TRUE. )
         IF ( ErrStat /= 0 ) THEN
            CLOSE( UnIn )
            RETURN
         END IF

         IF (     FP_InitData%MooringLine(I)%LineNodes < 0           ) THEN
            CALL ProgAbort( ' LineNodes('//TRIM( Int2LStr( I ) )//') must not be a negative number.', TrapErrors = .TRUE. )
            ErrStat = 1
            CLOSE( UnIn )
            RETURN       
!         ELSEIF ( FP_InitData%MooringLine(I)%LineNodes > MaxLineNodes ) THEN        ! error if we didn't read the whole array
!            CALL ProgAbort( ' LineNodes('//TRIM( Int2LStr( I ) )//') must not be larger than ' &
!                                        //TRIM(Int2LStr(MaxLineNodes))//'.', TrapErrors = .TRUE. )
!            ErrStat = 1
!            CLOSE( UnIn )
!            RETURN       
         END IF

        
         ALLOCATE ( FP_InitData%MooringLine(I)%LSNodes ( FP_InitData%MooringLine(I)%LineNodes ) , STAT=ErrStat )
         IF ( ErrStat /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the LSNodes array.', .TRUE. )
            CLOSE( UnIn )
            RETURN
         END IF
            
         DO J = 1,FP_InitData%MooringLine(I)%LineNodes
               
               ! We'll assume LineNodes equally spaced nodes from 0 to LUnstrLen:
               
            FP_InitData%MooringLine(I)%LSNodes(J) = (J - 0.5) * ( FP_InitData%MooringLine(I)%LUnstrLen / &
                                                                  FP_InitData%MooringLine(I)%LineNodes )
            
!               ! Check the user-input values:
!            
!            FP_InitData%MooringLine(I)%LSNodes(J) = TmpAry(J)                                    
!
!            IF ( ( FP_InitData%MooringLine(I)%LSNodes(J) < 0                                    ) .OR.  &
!                 ( FP_InitData%MooringLine(I)%LSNodes(J) > FP_InitData%MooringLine(I)%LUnstrLen )  ) THEN                 
!               CALL ProgAbort ( ' Error: Check that 0 <= LSNodes('//TRIM( Int2LStr( J ) )//') <= LUnstrLen on mooring line ' &
!                                                                  //TRIM( Int2LStr( I ) )//'.', .TRUE. )
!               CLOSE( UnIn )
!               ErrStat = 1
!               RETURN
!            END IF                             
         END DO !J            

         IF ( Echo )  THEN
            WRITE( TmpFmt, '( 12( 2X, '//TRIM( HD_Data%HDOut_Data%OutFmt )// '), 2X, I5  ) ' )                  &
                                    LRadAnch                            , LAngAnch                            , &
                                    LDpthAnch                           , LRadFair                            , &
                                    LAngFair                            , LDrftFair                           , &
                                    FP_InitData%MooringLine(I)%LUnstrLen, FP_InitData%MooringLine(I)%LDiam    , &
                                    FP_InitData%MooringLine(I)%LMassDen , FP_InitData%MooringLine(I)%LEAStff  , &
                                    FP_InitData%MooringLine(I)%LSeabedCD, FP_InitData%MooringLine(I)%LTenTol  , &
                                    FP_InitData%MooringLine(I)%LineNodes
            
            
            IF ( FP_InitData%MooringLine(I)%LineNodes  > 0 ) THEN       ! break this up for gfortran
            
               WRITE (UnEc,'( A, 4X, ' //TRIM( Int2LStr(FP_InitData%MooringLine(I)%LineNodes) )//'( 2X, '   &
                                       //TRIM( HD_Data%HDOut_Data%OutFmt )//' ) )'                        ) &
                                         TRIM( TmpFmt ),  FP_InitData%MooringLine(I)%LSNodes
                                      
            ELSE
               WRITE (UnEc, '( A )' ) TRIM( TmpFmt)
            END IF                                      
                              
         ENDIF

         IF ( LRadAnch                             <  0.0       ) THEN
            CALL ProgAbort ( ' LRadAnch('//TRIM( Int2LStr( I ) )//') must not be less than zero.', .TRUE. )
            ErrStat = 1
            CLOSE( UnIn )
            RETURN
         END IF

         IF ( LRadFair                             <  0.0       ) THEN
            CALL ProgAbort ( ' LRadFair('//TRIM( Int2LStr( I ) )//') must not be less than zero.', .TRUE. )
            ErrStat = 1
            CLOSE( UnIn )
            RETURN
         END IF

         IF ( LDrftFair                            <  0.0       ) THEN
            CALL ProgAbort ( ' LDrftFair('//TRIM( Int2LStr( I ) )//') must not be less than zero.', .TRUE. )
            ErrStat = 1
            CLOSE( UnIn )
            RETURN
         END IF

         IF ( LDpthAnch                            <  LDrftFair ) THEN
            CALL ProgAbort ( ' LDpthAnch('//TRIM( Int2LStr( I ) )//') must not be less than'// &
                           ' LDrftFair('//TRIM( Int2LStr( I ) )//').', .TRUE.                 )
            ErrStat = 1
            CLOSE( UnIn )
            RETURN
         END IF
         
         IF ( Waves_InitData%WtrDpth < LDpthAnch               )  THEN
            CALL ProgAbort ( ' WtrDpth must not be less than LDpthAnch('//TRIM( Int2LStr( I ) )//').', .TRUE. )
            ErrStat = 1
            CLOSE( UnIn )
            RETURN
         END IF                  

         IF ( FP_InitData%MooringLine(I)%LUnstrLen <= 0.0       ) THEN
            CALL ProgAbort ( ' LUnstrLen('//TRIM( Int2LStr( I ) )//') must be greater than zero.', .TRUE. )
            ErrStat = 1
            CLOSE( UnIn )
            RETURN
         END IF

         IF ( FP_InitData%MooringLine(I)%LDiam     <  0.0       ) THEN
            CALL ProgAbort ( ' LDiam('//TRIM( Int2LStr( I ) )//') must not be less than zero.', .TRUE. )
            ErrStat = 1
            CLOSE( UnIn )
            RETURN
         END IF

         IF ( FP_InitData%MooringLine(I)%LMassDen  <  0.0       ) THEN
            CALL ProgAbort ( ' LMassDen('//TRIM( Int2LStr( I ) )//') must not be less than zero.', .TRUE. )
            ErrStat = 1
            CLOSE( UnIn )
            RETURN
         END IF

         IF ( FP_InitData%MooringLine(I)%LEAStff   <= 0.0       ) THEN
            CALL ProgAbort ( ' LEAStff('//TRIM( Int2LStr( I ) )//') must be greater than zero.', .TRUE. )
            ErrStat = 1
            CLOSE( UnIn )
            RETURN
         END IF

               ! NOTE: Values of LSeabedCD less than zero indicate no seabed interaction (i.e., the mooring line may fall below the anchor)

         IF ( FP_InitData%MooringLine(I)%LTenTol   <= 0.0       ) THEN
            CALL ProgAbort ( ' LTenTol('//TRIM( Int2LStr( I ) )//') must be greater than zero.', .TRUE. )
            ErrStat = 1
            CLOSE( UnIn )
            RETURN
         END IF


         HD_Data%MaxDiam                    = MAX( HD_Data%MaxDiam, LRadAnch )         ! Find the maximum value of the input array LRadAnch

         LAngAnch                           =  LAngAnch*D2R                            ! Convert the azimuth angle of the current
         LAngFair                           =  LAngFair*D2R                            ! anchor and fairlead from degrees to radians.

         FP_InitData%MooringLine(I)%LAnchxi =  LRadAnch *COS(LAngAnch)                 !
         FP_InitData%MooringLine(I)%LAnchyi =  LRadAnch *SIN(LAngAnch)                 ! Convert the radius, azimuth angle, and depth of the current anchor   to x-, y-, and z-coordinates into the inertial frame       coordinate system
         FP_InitData%MooringLine(I)%LAnchzi = -LDpthAnch                               !
         FP_InitData%MooringLine(I)%LFairxt =  LRadFair *COS(LAngFair)                 !
         FP_InitData%MooringLine(I)%LFairyt =  LRadFair *SIN(LAngFair)                 ! Convert the radius, azimuth angle, and draft of the current fairlead to x-, y-, and z-coordinates into the tower base / platform coordinate system
         FP_InitData%MooringLine(I)%LFairzt = -LDrftFair                               !


      END DO             ! I - All mooring lines

!   ELSE
!   
!      DO I = 1,FP_InitData%NumLines ! Loop through all mooring lines
!      
!         CALL ReadCom( UnIn, FileName, 'unused mooring line '//TRIM( Int2LStr( I ) )//' properties', ErrStat )
!
!         IF ( ErrStat /= 0 ) THEN
!            CLOSE( UnIn )
!            RETURN
!         END IF         
!      
!      END DO

   END IF


   !-------------------------------------------------------------------------------------------------
   ! Data section for output
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Output header', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF


!      ! NWaveKin - Number of points where the incident wave kinematics can be output
!
!   CALL ReadVar ( UnIn, FileName, HD_Data%HDOut_Data%NWaveKin, 'NWaveKin', &
!      'Number of points where the incident wave kinematics can be output', ErrStat )
!
!   IF ( ErrStat /= 0 ) THEN
!      CLOSE( UnIn )
!      RETURN
!   END IF
!
!
!   IF ( ( HD_Data%HDOut_Data%NWaveKin < 0 ) .OR. ( HD_Data%HDOut_Data%NWaveKin > 9 ) ) THEN
!      CALL ProgAbort ( ' NWaveKin must be between 0 and 9 (inclusive).', .TRUE. )
!      ErrStat = 1
!      CLOSE( UnIn )
!      RETURN
!   END IF


!      ! WaveKinNd - List of tower nodes that have wave kinematics sensors.
!
!   IF ( HD_Data%HDOut_Data%NWaveKin > 0 ) THEN
!      CALL ReadAry ( UnIn, FileName, HD_Data%HDOut_Data%WaveKinNd(1:HD_Data%HDOut_Data%NWaveKin), HD_Data%HDOut_Data%NWaveKin, &
!                         'WaveKinNd', 'List of tower nodes that have wave kinematics sensors', ErrStat )
!
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!      
!         ! Check to see if all WaveKinNd(:) analysis points are existing analysis points:
!
!      DO I=1,HD_Data%HDOut_Data%NWaveKin
!              
!         IF ( ( HD_Data%HDOut_Data%WaveKinNd(I) < 1 ) .OR. ( HD_Data%HDOut_Data%WaveKinNd(I) > HD_Data%NElements ) )  THEN 
!            CALL ProgAbort ( ' All WaveKinNd values must be between 1 and '//TRIM( Int2LStr( HD_Data%NElements ) )//' (inclusive).', .TRUE. )
!            ErrStat = 1
!            CLOSE( UnIn )
!            RETURN
!         END IF
!                           
!      END DO ! I      
!
!   ELSE
!   
!      CALL ReadCom ( UnIn, FileName, 'unused WaveKinNd', ErrStat )
!
!      IF ( ErrStat /= 0 ) THEN
!         CLOSE( UnIn )
!         RETURN
!      END IF
!
!   END IF


      ! OutList - list of requested parameters to output to a file

   CALL ReadOutputList ( UnIn, FileName, HDOut_InitData%OutList, HD_Data%HDOut_Data%NumOuts, &
                                              'OutList', 'List of outputs requested', ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CLOSE( UnIn )
      RETURN
   END IF


   !-------------------------------------------------------------------------------------------------
   ! This is the end of the input file
   !-------------------------------------------------------------------------------------------------
   CLOSE ( UnIn )
   RETURN    
   
   
END SUBROUTINE HD_GetInput
!====================================================================================================
SUBROUTINE HD_GetDiscretization(UnIn, FileName, HD_Data, NumMemberNodes, ErrStat)
! This private subroutine is used to read the joint, morison member property, and morison member
! sections from the HydroDyn input file.  With this data, it creates the discretization that
! HydroDyn will use, storing in it in the HD_Data%Nodes data structure.
!----------------------------------------------------------------------------------------------------

      ! Passed variables
   
   TYPE(HD_DataType),         INTENT( INOUT )   :: HD_Data              ! the hydrodyn data 
   INTEGER,                   INTENT( IN    )   :: UnIn                 ! Unit number of open input file
   CHARACTER(*),              INTENT( IN    )   :: FileName             ! Name of the HydroDyn input file   
   INTEGER,                   INTENT(   OUT )   :: NumMemberNodes(9)    ! the number of nodes for each of the first [9] members

   INTEGER,                   INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs  
   
   
      ! Local variable types
      
   TYPE JointType
      INTEGER                                   :: ID
      REAL(ReKi)                                :: Position(3)
   END TYPE JointType
   
   TYPE SetType
      INTEGER                                   :: ID
      REAL(ReKi)                                :: D  (2)
      REAL(ReKi)                                :: CA (2)
      REAL(ReKi)                                :: CD (2)
   END TYPE SetType
         
   TYPE MemberType
      INTEGER                                   :: NumElements
      INTEGER                                   :: SetIndx
      INTEGER                                   :: JointIndx (2)
   END TYPE MemberType

      ! Local variables
      
   REAL(ReKi)                                   :: DNode (3)
   REAL(ReKi)                                   :: DSlope
   REAL(ReKi)                                   :: CASlope
   REAL(ReKi)                                   :: CDSlope
   
   TYPE(JointType),   ALLOCATABLE               :: Joint    (:)                  
   TYPE(SetType),     ALLOCATABLE               :: Set      (:)                 
   TYPE(MemberType),  ALLOCATABLE               :: Member   (:)                 

   INTEGER                                      :: I                    ! Generic loop counter      
   INTEGER                                      :: J                    ! Generic loop counter      
   INTEGER                                      :: K                    ! Generic loop counter      
   INTEGER                                      :: JointID(2)
   INTEGER                                      :: NJoints              ! Number of joints defined
   INTEGER                                      :: NSets                ! Number of member property sets
   INTEGER                                      :: NMembers             ! Number of members defined  
   INTEGER                                      :: SetID

   LOGICAL                                      :: FoundIt              ! lets us know if we found a matching ID


   !-------------------------------------------------------------------------------------------------
   ! initialize values
   !-------------------------------------------------------------------------------------------------
   ErrStat   = 0
   NumMemberNodes = 0
   
      !..............................................................................................
      ! Joints
      !..............................................................................................

   CALL ReadCom( UnIn, FileName, 'Joints header', ErrStat )   
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine()
      RETURN
   END IF

      ! NJoints
      
   CALL ReadVar( UnIn, FileName, NJoints, 'NJoints', 'Number of joints in the structure', ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine()
      RETURN
   END IF

   
   IF ( NJoints < 2  .AND. NJoints /= 0 ) THEN
      CALL ExitThisRoutine( ' There must be either no joints or at least two joints defined in the HydroDyn dataset.' )
      ErrStat = 1
      RETURN
   ELSE
         ! Allocate Joint array
         
      ALLOCATE( Joint( NJoints ), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN   
         CALL ExitThisRoutine( ' Error allocating space for HydroDyn joints.' )
         ErrStat = 1
         RETURN
      END IF
         
   END IF
   
      ! Skip the 2 header/comment lines.

   CALL ReadCom ( UnIn, FileName, 'joint parameter names', ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine()
      RETURN
   END IF
   
   CALL ReadCom ( UnIn, FileName, 'joint parameter units', ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine()
      RETURN
   END IF
   
   
   ! Read in joint data: JointID, Jointxi, Jointyi, Jointzi
   
   DO I = 1,NJoints ! Loop through all joints

      READ (UnIn,*,IOSTAT=ErrStat) Joint(I)%ID, ( Joint(I)%Position(J), J=1,3 )           

      CALL CheckIOS ( ErrStat, FileName, 'Joint property (line '//TRIM( Int2LStr( I ) )//')', NumType, .TRUE. )
      IF ( ErrStat /= 0 )  THEN
         CALL ExitThisRoutine()
         RETURN
      END IF
      
      IF ( Echo )  THEN
         WRITE (UnEc,'( I5, 3( 2X, ' //TRIM( HD_Data%HDOut_Data%OutFmt )// ') )' ) & 
                                               Joint(I)%ID, ( Joint(I)%Position(J), J=1,3 )
      ENDIF


         ! check that the joint numbers are unique
         
      DO J = 1,I-1
         IF ( Joint(J)%ID == Joint(I)%ID ) THEN
            CALL ExitThisRoutine( ' Duplicate joint ID# '//TRIM( Int2LStr(Joint(J)%ID) )//' found.')
            ErrStat = 1
            RETURN      
         END IF
      END DO  !J
      
   END DO !I
   
   
      !..............................................................................................
      ! Morison member property sets
      !..............................................................................................


   CALL ReadCom( UnIn, FileName, 'Member properties header', ErrStat )   
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine()
      RETURN
   END IF
   

      ! NSets
      
   CALL ReadVar( UnIn, FileName, NSets, 'NSets', 'Number of member property sets', ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine()
      RETURN
   END IF

   
!   IF ( NSets < 1 ) THEN
!   
!      CALL ExitThisRoutine( ' There must be at least one member property set defined in the HydroDyn dataset.' )
   IF ( NSets < 0 ) THEN
   
      CALL ExitThisRoutine( ' The number of Morison member property sets defined in the HydroDyn dataset must not be negative.' )
      ErrStat = 1
      RETURN
      
   ELSE
         ! Allocate Set array
         
      ALLOCATE( Set( NSets ), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN   
         CALL ExitThisRoutine( ' Error allocating space for HydroDyn property sets.' )
         ErrStat = 1
         RETURN
      END IF
         
   END IF
   
   
      ! Skip the header/comment lines.

   CALL ReadCom ( UnIn, FileName, 'member property set parameter names', ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine()
      RETURN
   END IF
   
   CALL ReadCom ( UnIn, FileName, 'member property set parameter units', ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine()
      RETURN
   END IF
   
   
   ! read member property data
   
   DO I = 1,NSets ! Loop through all joints

      READ (UnIn,*,IOSTAT=ErrStat) Set(I)%ID, Set(I)%D(1), Set(I)%D(2), Set(I)%CA(1), Set(I)%CA(2), Set(I)%CD(1), Set(I)%CD(2)

      CALL CheckIOS ( ErrStat, FileName, 'Member property set (line '//TRIM( Int2LStr( I ) )//')', NumType, .TRUE. )
      IF ( ErrStat /= 0 )  THEN
         CALL ExitThisRoutine()
         RETURN
      END IF
      
      IF ( Echo )  THEN
         WRITE (UnEc,'( I5, 6( 2X, ' //TRIM( HD_Data%HDOut_Data%OutFmt )// ') )' ) & 
                                   Set(I)%ID, Set(I)%D(1), Set(I)%D(2), Set(I)%CA(1), Set(I)%CA(2), Set(I)%CD(1), Set(I)%CD(2)
      ENDIF

       
         ! check that the set id numbers are unique
         
      DO J = 1,I-1
         IF ( Set(J)%ID == Set(I)%ID ) THEN
            CALL ExitThisRoutine( ' Duplicate set ID# '//TRIM( Int2LStr(Set(J)%ID) )//' found.')
            ErrStat = 1
            RETURN      
         END IF
      END DO  ! J
      
         ! check that the properties are valid

      IF ( ( Set(I)%D( 1) < 0 ) .OR. ( Set(I)%D( 2) < 0 ) .OR. &
           ( Set(I)%CA(1) < 0 ) .OR. ( Set(I)%CA(2) < 0 ) .OR. &
           ( Set(I)%CD(1) < 0 ) .OR. ( Set(I)%CD(2) < 0 ) ) THEN  
            CALL ExitThisRoutine( " Values for D1, D2, CA1, CA2, CD1, and CD2 must not be negative in HydroDyn's Set ID# " & 
                                              //TRIM( Int2LStr(Set(I)%ID) )//".")
            ErrStat = 1
            RETURN      
      END IF
      

   END DO !I


      !..............................................................................................
      ! Morison members (with nodes/elements)
      !..............................................................................................

   CALL ReadCom( UnIn, FileName, 'Members header', ErrStat )   
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine()
      RETURN
   END IF
   

      ! NMembers
      
   CALL ReadVar( UnIn, FileName, NMembers, 'NMembers', 'Number of members in the structure', ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine()
      RETURN
   END IF
   
   
!   IF ( NMembers < 1 ) THEN
!      CALL ExitThisRoutine( ' There must be at least one member defined in the HydroDyn dataset.' )
! In the general case, there can be 0 members, but we don't have equations/interface to handle this yet. jmj says ptfmnodes could be 0 in FAST.

   IF ( NMembers < 0 ) THEN
      CALL ExitThisRoutine( ' The number of Morison members defined in the HydroDyn dataset must not be negative.' )
      ErrStat = 1
      RETURN
   ELSE
         ! Allocate Member array
         
      ALLOCATE( Member( NMembers ), STAT = ErrStat )
      IF ( ErrStat /= 0 ) THEN   
         CALL ExitThisRoutine( ' Error allocating space for HydroDyn members.' )
         ErrStat = 1
         RETURN
      END IF
         
   END IF
   
      ! Skip the header/comment lines.

   CALL ReadCom ( UnIn, FileName, 'member parameter names', ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ExitThisRoutine()
      RETURN
   END IF
   
!   CALL ReadCom ( UnIn, FileName, 'member parameter units', ErrStat )
!   IF ( ErrStat /= 0 )  THEN
!      CALL ExitThisRoutine()
!      RETURN
!   END IF
   
   
   ! read member data
   
   HD_Data%NElements = 0           ! initialize total number of elements
   
   DO I = 1,NMembers ! Loop through all members

      READ (UnIn,*,IOSTAT=ErrStat) JointID(1), JointID(2), SetID, Member(I)%NumElements

      CALL CheckIOS ( ErrStat, FileName, 'Member (line '//TRIM( Int2LStr( I ) )//')', NumType, .TRUE. )
      IF ( ErrStat /= 0 )  THEN
         CALL ExitThisRoutine()
         RETURN
      END IF
      
      IF ( Echo )  THEN
!         WRITE (UnEc,'( I5, 3( 2X, ' //TRIM( HD_Data%HDOut_Data%OutFmt )// ') )' ) & 
         WRITE (UnEc,'( I7, 3( 2X, I7) )' ) & 
                                   JointID(1), JointID(2), SetID, Member(I)%NumElements
      ENDIF
      
      
         ! check that NumElements is valid
      
      IF ( Member(I)%NumElements < 1 ) THEN
         CALL ExitThisRoutine( ' Each member in the HydroDyn dataset must contain at least one element.' )
         ErrStat = 1
         RETURN
      END IF

      HD_Data%NElements = HD_Data%NElements + Member(I)%NumElements
      
      
         ! check that the property set number is valid
         
      FoundIt = .FALSE.
      DO J = 1,NSets
         IF ( Set(J)%ID == SetID ) THEN
            Member(I)%SetIndx = J
            FoundIt           = .TRUE.            
            EXIT  ! exit this DO loop
         END IF
      END DO
         
      IF ( .NOT. FoundIt ) THEN  ! we didn't find it
         CALL ExitThisRoutine( ' Member (line '//TRIM( Int2LStr( I ) )//') contains an invalid property set number.' )
         ErrStat = 1
         RETURN      
      END IF
     
     
         ! check that the joint numbers are valid
         
      DO K = 1,2    
           
         FoundIt = .FALSE.
         
         DO J = 1,NJoints
            IF ( Joint(J)%ID == JointID(K) ) THEN
               Member(I)%JointIndx(K) = J
               FoundIt                = .TRUE.            
               EXIT  ! exit this DO loop
            END IF
         END DO
            
         IF ( .NOT. FoundIt ) THEN  ! we didn't find it
            CALL ExitThisRoutine( ' Member (line '//TRIM( Int2LStr( I ) )//') contains an invalid Joint' &
                                                  //TRIM( Int2LStr( K ) )//' number.' )
            ErrStat = 1
            RETURN      
         END IF
         
      END DO !K
      
   END DO !I - members


   !-------------------------------------------------------------------------------------------------
   ! Set return value for HD_Output (this array was already initialized to 0)
   !-------------------------------------------------------------------------------------------------
   DO I = 1, MIN( NMembers, SIZE(NumMemberNodes) )
      NumMemberNodes(I) = Member(I)%NumElements
   END DO      


   !-------------------------------------------------------------------------------------------------
   ! We will check that the discretization follows the criteria that are required for the current
   ! version of the FloatingPlatform and FixedBottomSupportStructure modules.
   ! When these modules are more robust, this check can be removed.
   !-------------------------------------------------------------------------------------------------

   SELECT CASE ( HD_Data%StrctType )
   
      CASE ( FloatPltfm_Type )
      
         ! for now, we check that there is one member, that d1=d2, ca(:)=0, cd1=cd2, and jointxi=jointyi=0
         ! also, ( jointzi(1) = 0 and jointzi(2) < 0 ) OR ( jointzi(2) = 0 and jointzi(1) < 0 )
      
         IF      ( NMembers /= 1                                                                               ) THEN
            CALL ExitThisRoutine( ' Exactly one Morison member must be defined for floating platform substructures.' )
            RETURN
         ELSE IF ( Set( Member(1)%SetIndx )%D(1) /= Set( Member(1)%SetIndx )%D(2)                              ) THEN      ! note comparison of floating points!
            CALL ExitThisRoutine( ' Morison member properties D1 and D2 must be equal for floating platform substructures.' )
            RETURN         
         ELSE IF ( Set( Member(1)%SetIndx )%CA(1) /= 0.0_ReKi .OR. Set( Member(1)%SetIndx )%CA(2) /= 0.0_ReKi  ) THEN      ! note comparison of floating points!
            CALL ExitThisRoutine( ' Morison member properties CA1 and CA2 must be zero for floating platform substructures.' )
            RETURN
         ELSE IF ( Set( Member(1)%SetIndx )%CD(1) /= Set( Member(1)%SetIndx )%CD(2)                            ) THEN      ! note comparison of floating points!
            CALL ExitThisRoutine( ' Morison member properties CD1 and CD2 must be equal for floating platform substructures.' )
            RETURN         
         ELSE IF ( Joint( Member(1)%JointIndx(1) )%Position(1) /= 0.0_ReKi .OR. &
                   Joint( Member(1)%JointIndx(1) )%Position(2) /= 0.0_ReKi .OR. &
                   Joint( Member(1)%JointIndx(2) )%Position(1) /= 0.0_ReKi .OR. &
                   Joint( Member(1)%JointIndx(2) )%Position(2) /= 0.0_ReKi                                     ) THEN      ! note comparison of floating points!
            CALL ExitThisRoutine( ' Jointxi and Jointyi must be zero for floating platform substructures.' )
            RETURN         
         ELSE IF ( Joint( Member(1)%JointIndx(1) )%Position(3) /= 0.0_ReKi                                     ) THEN      ! note comparison of floating points!                                          
            IF ( Joint( Member(1)%JointIndx(1) )%Position(3) >= 0.0_ReKi .OR. &
                 Joint( Member(1)%JointIndx(2) )%Position(3) /= 0.0_ReKi ) THEN
               CALL ExitThisRoutine( &
                         ' The two values of Jointzi must be zero and a negative number for floating platform substructures.' )
               RETURN         
            END IF
         ELSE IF  ( Joint( Member(1)%JointIndx(2) )%Position(3) /= 0.0_ReKi                                    ) THEN      ! note comparison of floating points!
            IF ( Joint( Member(1)%JointIndx(2) )%Position(3) >= 0.0_ReKi .OR. &
                 Joint( Member(1)%JointIndx(1) )%Position(3) /= 0.0_ReKi ) THEN
               CALL ExitThisRoutine( &
                         ' The two values of Jointzi must be zero and a negative number for floating platform substructures.' )
               RETURN         
            END IF           
         END IF            
      
      CASE ( FixedBtm_Type )
      
         ! for now, we check that there is one member, that d1=d2, ca1=ca2, cd1=cd2, and jointxi=jointyi=0
!bjj is this correct????         ! also, jointzi(1) = 0 and jointzi(2)<0

         IF      ( NMembers /= 1                                                                               ) THEN
            CALL ExitThisRoutine( ' Exactly one Morison member must be defined for fixed bottom substructures.' )
            RETURN
         ELSE IF ( Set( Member(1)%SetIndx )%D(1)  /= Set( Member(1)%SetIndx )%D(2)                             ) THEN      ! note comparison of floating points!
            CALL ExitThisRoutine( ' Morison member properties D1 and D2 must be equal for fixed bottom substructures.' )
            RETURN         
         ELSE IF ( Set( Member(1)%SetIndx )%CA(1) /= Set( Member(1)%SetIndx )%CA(2)                            ) THEN      ! note comparison of floating points!
            CALL ExitThisRoutine( ' Morison member properties CA1 and CA2 must be equal for fixed bottom substructures.' )
            RETURN
         ELSE IF ( Set( Member(1)%SetIndx )%CD(1) /= Set( Member(1)%SetIndx )%CD(2)                            ) THEN      ! note comparison of floating points!
            CALL ExitThisRoutine( ' Morison member properties CD1 and CD2 must be equal for fixed bottom substructures.' )
            RETURN         
         ELSE IF ( Joint( Member(1)%JointIndx(1) )%Position(1) /= 0.0_ReKi .OR. &
                   Joint( Member(1)%JointIndx(1) )%Position(2) /= 0.0_ReKi .OR. &
                   Joint( Member(1)%JointIndx(2) )%Position(1) /= 0.0_ReKi .OR. &
                   Joint( Member(1)%JointIndx(2) )%Position(2) /= 0.0_ReKi                                     ) THEN      ! note comparison of floating points!
            CALL ExitThisRoutine( ' Jointxi and Jointyi must be zero for fixed bottom substructures.' )
            RETURN         
!         ELSE IF ( Joint( Member(1)%JointIndx(1) )%Position(3) /= 0.0_ReKi                                     ) THEN      ! note comparison of floating points!                                          
!            IF ( Joint( Member(1)%JointIndx(1) )%Position(3) >= 0.0_ReKi .OR. &
!                 Joint( Member(1)%JointIndx(2) )%Position(3) /= 0.0_ReKi ) THEN
!               CALL ExitThisRoutine( ' The two values of Jointzi must be zero and a negative number for fixed bottom substructures.' )
!               RETURN         
!            END IF
!         ELSE IF  ( Joint( Member(1)%JointIndx(2) )%Position(3) /= 0.0_ReKi                                    ) THEN      ! note comparison of floating points!
!            IF ( Joint( Member(1)%JointIndx(2) )%Position(3) >= 0.0_ReKi .OR. &
!                 Joint( Member(1)%JointIndx(1) )%Position(3) /= 0.0_ReKi ) THEN
!               CALL ExitThisRoutine( ' The two values of Jointzi must be zero and a negative number for fixed bottom substructures.' )
!               RETURN         
!            END IF           
         END IF            
      
      
!      CASE DEFAULT
   
   END SELECT


   !-------------------------------------------------------------------------------------------------
   ! Create the nodes/elements for the discretization using above inputs
   !-------------------------------------------------------------------------------------------------
   
      ! Allocate Node array
   ALLOCATE ( HD_Data%Node( HD_Data%NElements ), STAT = ErrStat )   
   IF ( ErrStat /= 0 ) THEN
      CALL ExitThisRoutine( ' Error allocating space for HydroDyn nodes array.' )
      RETURN
   END IF
   
      ! FOR NOW, WE'LL CREATE EVENLY SPACED ELEMENTS ON EACH MEMBER
      ! this is something we can make a bit more sophisticated in the future
      
   J = 1                      ! counter for the current element
         
   DO I = 1,NMembers
   
!      HD_Data%MaxDiam = MAX( HD_Data%MaxDiam, SQRT( DOT_PRODUCT(Joint( Member(I)%JointIndx(1) )%Position(1:2),    &
!                                                                Joint( Member(I)%JointIndx(1) )%Position(1:2))) , &
!                                              SQRT( DOT_PRODUCT(Joint( Member(I)%JointIndx(2) )%Position(1:2)     &
!                                                                Joint( Member(I)%JointIndx(2) )%Position(1:2)))   )  ! let's make sure the graphics also cover the JOINTs we've defined
   
      DNode   = ( Joint( Member(I)%JointIndx(2) )%Position - Joint( Member(I)%JointIndx(1) )%Position ) / Member(I)%NumElements
      DSlope  = ( Set  ( Member(I)%SetIndx      )%D( 2)    - Set  ( Member(I)%SetIndx      )%D( 1)    ) / Member(I)%NumElements
      CASlope = ( Set  ( Member(I)%SetIndx      )%CA(2)    - Set  ( Member(I)%SetIndx      )%CA(1)    ) / Member(I)%NumElements
      CDSlope = ( Set  ( Member(I)%SetIndx      )%CD(2)    - Set  ( Member(I)%SetIndx      )%CD(1)    ) / Member(I)%NumElements
                  
      DO K = 1,Member(I)%NumElements
            
         HD_Data%Node(J)%DNode     = SQRT( DOT_PRODUCT( DNode, DNode ) )
         HD_Data%Node(J)%Position  = Joint( Member(I)%JointIndx(1) )%Position + ( K - 0.5 )*DNode    ! the center of evenly spaced elements
         HD_Data%Node(J)%D         = Set  ( Member(I)%SetIndx      )%D( 1)    + ( K - 0.5 )*DSlope
         HD_Data%Node(J)%CA        = Set  ( Member(I)%SetIndx      )%CA(1)    + ( K - 0.5 )*CASlope
         HD_Data%Node(J)%CD        = Set  ( Member(I)%SetIndx      )%CD(1)    + ( K - 0.5 )*CDSlope
         
!print *, J, HD_Data%Node(J)         
!bjj check that this interpolation is correct for all of the values
         
         J = J + 1                              ! move on to the next element
      END DO               
   
   END DO !I      
   
   
   !-------------------------------------------------------------------------------------------------
   ! Clean up the allocated data
   !-------------------------------------------------------------------------------------------------  
   ErrStat = 0   
   CALL ExitThisRoutine ()   
   
   RETURN
   
CONTAINS

   !=================================================================================================
   SUBROUTINE ExitThisRoutine( ErrMsg )   
   ! This subroutine cleans up after the enveloping routine
   !-------------------------------------------------------------------------------------------------
         ! local variables
         
      INTEGER                            :: AllocStat
      CHARACTER(*), INTENT(IN), OPTIONAL :: ErrMsg
         
      !----------------------------------------------------------------------------------------------
      ! If this subroutine ends on an error, print the message here:
      !----------------------------------------------------------------------------------------------   
      IF ( PRESENT( ErrMsg ) ) CALL ProgAbort( ErrMsg, TrapErrors = .TRUE. )

      !----------------------------------------------------------------------------------------------
      ! Deallocate space from local variables (in case it's compiled with /QSAVE)
      !----------------------------------------------------------------------------------------------
      
      IF ( ALLOCATED( Joint  ) ) DEALLOCATE ( Joint , STAT = AllocStat )
      IF ( ALLOCATED( Set    ) ) DEALLOCATE ( Set   , STAT = AllocStat )
      IF ( ALLOCATED( Member ) ) DEALLOCATE ( Member, STAT = AllocStat )
      
        
   END SUBROUTINE ExitThisRoutine
   !=================================================================================================
   
END SUBROUTINE HD_GetDiscretization
!====================================================================================================
FUNCTION HD_GetValue_CHAR(VarName, HD_Data, ErrStat)
!  This function returns a real scalar value whose name is listed in the VarName input argument.
!  If the name is not recognized, an error is returned in ErrStat.
!----------------------------------------------------------------------------------------------------

   CHARACTER(*),        INTENT( IN    )   :: VarName
   TYPE( HD_DataType ), INTENT( INOUT )   :: HD_Data
   INTEGER,             INTENT(   OUT )   :: ErrStat
   REAL(ReKi)                             :: HD_GetValue_CHAR        ! This function


   INTEGER                                :: Indx

   CHARACTER(20)                          :: VarNameUC


   !-------------------------------------------------------------------------------------------------
   ! Set the initial error status and return value
   !-------------------------------------------------------------------------------------------------
   
   ErrStat          = 0
   HD_GetValue_CHAR = 0
   
   !-------------------------------------------------------------------------------------------------
   ! Warn if HydroDyn hasn't been properly initialized
   !-------------------------------------------------------------------------------------------------
   IF ( HD_Data%StrctType == Unknown_Type ) THEN
      CALL ProgWarn( ' HydroDyn has not been initialized before calling HD_GetValue().')
      ErrStat = -1
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Return the requested values.
   !-------------------------------------------------------------------------------------------------

   VarNameUC = VarName
   CALL Conv2UC( VarNameUC )

   SELECT CASE ( TRIM(VarNameUC) )

      CASE ( 'PTFMDIAM', 'PLATFORMDIAMETER' )
         IF ( HD_DATA%StrctType == FloatPltfm_Type .AND. ALLOCATED(HD_Data%Node) ) THEN
            IF ( SIZE(HD_Data%Node) > 0 ) THEN
               HD_GetValue_CHAR = HD_Data%Node(1)%D
            ELSE
               HD_GetValue_CHAR = 0.0
               ErrStat = 1
            END IF               
         ELSE
            HD_GetValue_CHAR = 0.0
            ErrStat = 1
         END IF
         
      CASE ( 'WTRDENS', 'WATERDENSITY' )
         HD_GetValue_CHAR = HD_Data%Waves_Data%WtrDens

      CASE ( 'WTRDPTH', 'WATERDEPTH' )
         HD_GetValue_CHAR = HD_Data%Waves_Data%WtrDpth

      CASE ( 'MAXLRADANCH', 'MAXDIAM' )
         HD_GetValue_CHAR = HD_Data%MaxDiam                     

      CASE ( 'NWAVEELEV' )
         HD_GetValue_CHAR = REAL( HD_Data%Waves_Data%NWaveElev, ReKi )

      CASE ( 'WAVEDIR', 'WAVEDIRECTION' )
         HD_GetValue_CHAR = HD_Data%Waves_Data%WaveDir                     

      CASE DEFAULT

!bjj start of proposed change v1.00.00c-bjj
!         IF ( VarNameUC(1:8) == 'WAVEELEV' ) THEN
!         
!         ELSEIF ( VarNameUS(1:7) == 'LINEPOS' ) THEN
!         
!         END IF
!bjj end of proposed change v1.00.00c-bjj
         

      
         CALL WrScr( ' Invalid variable name "'//TRIM(VarName)//'" in HD_GetValue().' )
         ErrStat = 1
         HD_GetValue_CHAR = 0.0

   END SELECT
   
END FUNCTION HD_GetValue_CHAR
!====================================================================================================
FUNCTION HD_GetValue_AllOuts(OutputID, HD_Data, ErrStat)
!  This function returns a real scalar value whose parameter ID is listed in the OutputID input argument.
!  If the ID is not recognized, an error is returned in ErrStat.
!  To use this function, please USE the parameters listed in MODULE HD_Output as the OutputID (these 
!  are the same parameter names as are stored in the file "HydroDynOutListParameters.xlsx"
!----------------------------------------------------------------------------------------------------

   INTEGER,             INTENT( IN    )   :: OutputID
   TYPE( HD_DataType ), INTENT( INOUT )   :: HD_Data
   INTEGER,             INTENT(   OUT )   :: ErrStat
   REAL(ReKi)                             :: HD_GetValue_AllOuts

   !-------------------------------------------------------------------------------------------------
   ! Set the initial error status and return value
   !-------------------------------------------------------------------------------------------------   
   ErrStat             = 0
   HD_GetValue_AllOuts = 0
   
   !-------------------------------------------------------------------------------------------------
   ! Warn if HydroDyn hasn't been properly initialized
   !-------------------------------------------------------------------------------------------------
   IF ( HD_Data%StrctType == Unknown_Type ) THEN
      CALL ProgWarn( ' HydroDyn has not been initialized before calling HD_GetValue().')
      ErrStat = -1
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Check that OutputID is valid
   !-------------------------------------------------------------------------------------------------
   IF ( OutputID < 1 .OR. OutputID > MaxOutPts ) THEN
      CALL WrScr( ' Invalid output variable ID "'//TRIM(Int2LStr(OutputID))//'" in HD_GetValue().' )
      ErrStat = 1
      RETURN
   END IF
   

   !-------------------------------------------------------------------------------------------------
   ! Return the requested values.
   !-------------------------------------------------------------------------------------------------  

   HD_GetValue_AllOuts = HD_Data%HDOut_Data%AllOuts(OutputID)
   

END FUNCTION HD_GetValue_AllOuts
!====================================================================================================
FUNCTION HD_GetUndisturbedWaveElev(Position, Time, HD_Data, ErrStat)

      ! Passed variables   
   REAL(ReKi),          INTENT( IN    )   :: Position(2)                ! the x and y coordinates where the wave elevation is requested, in meters
   REAL(DbKi),          INTENT( IN    )   :: Time                       ! the time at which the wave elevation is desired (need not be CurrentTime), in seconds
   
   TYPE( HD_DataType ), INTENT( INOUT )   :: HD_Data                    ! the module data (INOUT b/c Waves Index is changed in the WaveElevation function)
   INTEGER,             INTENT(   OUT )   :: ErrStat                    ! an error code
   
   REAL(ReKi)                             :: HD_GetUndisturbedWaveElev  ! the wave elevation (z) at Position(:), in meters


   !-------------------------------------------------------------------------------------------------
   ! Set the initial error status and return value
   !-------------------------------------------------------------------------------------------------   
   ErrStat                   = 0
   HD_GetUndisturbedWaveElev = 0
   
   !-------------------------------------------------------------------------------------------------
   ! Warn if HydroDyn hasn't been properly initialized
   !-------------------------------------------------------------------------------------------------
   IF ( HD_Data%StrctType == Unknown_Type ) THEN
      CALL ProgAbort( ' HydroDyn has not been initialized before calling HD_GetUndisturbedWaveElev().', TrapErrors = .TRUE.)
      ErrStat = -1
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! interpolate into wave
   !-------------------------------------------------------------------------------------------------

   HD_GetUndisturbedWaveElev = Waves_GetUndisturbedElev( Position, Time, HD_Data%Waves_Data, ErrStat )
   
   IF ( ErrStat /= 0 ) THEN
      CALL ProgAbort( ' Unable to get wave elevation in HD_GetUndisturbedWaveElev().', TrapErrors = .TRUE.)
   END IF      

   RETURN

END FUNCTION
!====================================================================================================
SUBROUTINE HD_Init( HydroDyn_InitData, HD_ConfigMarkers, HD_AllMarkers, HD_Data, ErrStat )
!     This public subroutine initializes the HydroDyn module.
! 
!----------------------------------------------------------------------------------------------------

!must return (1) discretization and (2) the dimension [0=point load, 1=load per unit length, 2=load per unit area, 3=load per unit volume]

      ! Passed variables

   TYPE(HD_InitDataType),      INTENT( IN  )    :: HydroDyn_InitData    
   TYPE(HydroConfig),          INTENT( IN  )    :: HD_ConfigMarkers
   TYPE(AllHydroMarkers),      INTENT( OUT )    :: HD_AllMarkers
   TYPE(HD_DataType),          INTENT( OUT )    :: HD_Data              ! a return value (the data initialized here)
   INTEGER,                    INTENT( OUT )    :: ErrStat              ! Returns zero if no errors were encountered
                
                
      ! Internal variables                
                
   TYPE(Current_DataType)                       :: Current_Data         ! the values needed to initialize the current in the waves module
   TYPE(FltPtfm_InitDataType)                   :: FltPtfm_InitData     ! the values needed to initialize the the floating platform module
   TYPE(HDO_InitDataType)                       :: HDO_InitData         ! the values needed to initialize the the HydroDyn Output module   
   TYPE(Waves_InitDataType)                     :: Waves_InitData       ! the values needed to initialize the the waves module
   
   
   INTEGER                                      :: I                    ! Internal loop counter
   INTEGER                                      :: NMarkerElements      ! number of markers made up of the HydroDyn Node elements
   INTEGER                                      :: NMarkers             ! number of markers where loads will be calculated
   
   CHARACTER(1024)                              :: SummaryName          ! name of the HydroDyn summary file
                
      ! Initialize NWTC_Library
      
   CALL NWTC_Init( )
   
   !-------------------------------------------------------------------------------------------------
   ! Display the HydroDyn version on the screen
   !-------------------------------------------------------------------------------------------------

   CALL WrScr1 ( ' Using '//TRIM( HD_Prog%Name )//' '//TRIM( HD_Prog%Ver )//'.' )  

   SummaryName = TRIM(HydroDyn_InitData%OutRootName)//'_HydroDyn.hds'
   CALL HDOut_OpenSum( HD_Data%UnSum, SummaryName, HD_Prog, ErrStat )    !this must be called before the Waves_Init() routine so that the appropriate wave data can be written to the summary file
   IF ( ErrStat /= 0 ) RETURN

   !-------------------------------------------------------------------------------------------------
   ! Read data from HydroDyn input file 
   !-------------------------------------------------------------------------------------------------
      
   CALL HD_GetInput( HD_Data, HydroDyn_InitData%FileName, Current_Data, Waves_InitData, FltPtfm_InitData, HDO_InitData, ErrStat )   
   IF ( ErrStat /= 0 ) RETURN
      
   
      ! Check HD input with values in input file
   IF ( HD_Data%StrctType == FloatPltfm_Type ) THEN

      IF ( HydroDyn_InitData%Gravity <= 0.0 ) THEN
         CALL ProgAbort ( ' Gravity must be greater than zero when using a floating platform.', .TRUE. )
         ErrStat = 1
         RETURN
      END IF
      
   END IF         
         
   !-------------------------------------------------------------------------------------------------
   ! Initialize wave module 
   ! (this must happen before initializing other support structure module and output module!)
   !-------------------------------------------------------------------------------------------------
      
   Waves_InitData%Gravity = HydroDyn_InitData%Gravity
   Waves_InitData%DirRoot = HydroDyn_InitData%OutRootName
   

   Waves_InitData%NWaveKin0 = HD_Data%NElements
   
   ALLOCATE ( Waves_InitData%DZNodes(   Waves_InitData%NWaveKin0 ), STAT = ErrStat )
   IF ( ErrStat /= 0 ) THEN
      CALL ProgAbort ( ' Error allocating space for DZNodes array in HD_Init().', TrapErrors = .TRUE. )
      RETURN
   END IF   

   ALLOCATE ( Waves_InitData%WaveKinzi0( Waves_InitData%NWaveKin0 ), STAT = ErrStat )
   IF ( ErrStat /= 0 ) THEN
      CALL ProgAbort ( ' Error allocating space for WaveKinzi0 array in HD_Init().', TrapErrors = .TRUE. )
      RETURN !bjj memory leaks with /QSAVE !!!!
   END IF      


   DO I = 1,HD_Data%NElements 
      Waves_InitData%DZNodes   (I) = HD_Data%Node(I)%DNode
      Waves_InitData%WaveKinzi0(I) = HD_Data%Node(I)%Position(3) ! bjj: setting WaveKinzi0 this way assumes Position(1:2) /= 0; if we allow more compilcated structures, this needs to be changed!!!!
   END DO ! I      
  
      
   CALL InitWaves ( Waves_InitData, Current_Data, HD_Data%Waves_data, HD_Data%UnSum, ErrStat   )      
   IF ( ErrStat /= 0 ) RETURN !bjj memory leaks with /QSAVE !!!!
      
      
      ! clean up the variables we've allocated (though they should get cleared automatically at subroutine exit, unless /QSAVE)
      
   IF ( ALLOCATED( Waves_InitData%DZNodes    ) ) DEALLOCATE ( Waves_InitData%DZNodes    ) 
   IF ( ALLOCATED( Waves_InitData%WaveElevxi ) ) DEALLOCATE ( Waves_InitData%WaveElevxi ) 
   IF ( ALLOCATED( Waves_InitData%WaveElevyi ) ) DEALLOCATE ( Waves_InitData%WaveElevyi ) 
   IF ( ALLOCATED( Waves_InitData%WaveKinzi0 ) ) DEALLOCATE ( Waves_InitData%WaveKinzi0 ) 
      
      
   !-------------------------------------------------------------------------------------------------
   ! Initialize appropriate support structure module
   ! (this must happen after initializing the waves module and before the output module!)
   !-------------------------------------------------------------------------------------------------
      
   SELECT CASE ( HD_Data%StrctType )
   
      CASE ( FloatPltfm_Type )
      
         !...........................................................................................
         ! check that the configuration marker is where it should be for the current implementation
         ! of the floating platform module
         !...........................................................................................

         IF ( ANY( HD_ConfigMarkers%Substructure%Position /= 0.0_ReKi ) ) THEN  !note comparisons of Real numbers here
            CALL ProgAbort( ' The floating platform model in HydroDyn requires that the substructure position be (0, 0, 0).', &
                              TrapErrors = .TRUE. )
            ErrStat = 1
            RETURN
         END IF         
         
         !...........................................................................................
         ! initialize the rest of the data that needs to be sent to the floating platform module
         !...........................................................................................
      
         FltPtfm_InitData%DirRoot  = HydroDyn_InitData%OutRootName
         
            ! bjj: we have already checked that one node exists, but if that restriction changes,
            ! this will need to be changed, too (to prevent ungraceful crashing and other errors!):
                     
         FltPtfm_InitData%PtfmDiam = HD_Data%Node(1)%D          
         FltPtfm_InitData%PtfmCD   = HD_Data%Node(1)%CD                 
         
         !...........................................................................................
         ! call the initialization routine
         !...........................................................................................
         
         CALL InitFltngPtfmLd( FltPtfm_InitData, HD_Data%Waves_data, HD_Data%FltPtfm_data, ErrStat )
         IF ( ErrStat /= 0 ) RETURN
         
                  
         !...........................................................................................
         ! clean up after ourselves here:
         !...........................................................................................
         
         IF ( ALLOCATED( FltPtfm_InitData%MooringLine ) ) THEN
         
            DO I = 1,FltPtfm_InitData%NumLines
               IF ( ALLOCATED(  FltPtfm_InitData%MooringLine(I)%LSNodes ) ) &
                    DEALLOCATE( FltPtfm_InitData%MooringLine(I)%LSNodes )               
            END DO ! I
         
            DEALLOCATE ( FltPtfm_InitData%MooringLine ) 
         END IF
      
      
         !...........................................................................................
         !  Set output markers
         !...........................................................................................
         
         NMarkers        = 1
         NMarkerElements = 0
         
            ! 1 per element
         
         
         
            ! 1 for the fix-point loads calculated at the platform reference
         
      
      CASE ( FixedBtm_Type )
      
         ! bjj: we have already checked that one node exists

         ! nothing to do (initialize) here


         ! set number of markers here
         
         NMarkers        = HD_Data%NElements
         NMarkerElements = HD_Data%NElements
         
      
      CASE DEFAULT
      
         CALL ProgAbort( ' Unknown support structure type in HD_CalculateLoads().', TrapErrors = .TRUE. )
         ErrStat = 1
         RETURN
   
   END SELECT
   
      
   !-------------------------------------------------------------------------------------------------
   ! Initizlize HydroDyn Output module
   ! (this must happen after initializing other support structure modules and wave module!)
   !-------------------------------------------------------------------------------------------------     
      
   HDO_InitData%ProgInfo    = TRIM( HD_Prog%Name )//' '//TRIM( HD_Prog%Ver )
   HDO_InitData%OutFileName = TRIM(HydroDyn_InitData%OutRootName)//'_HydroDyn.out'     ! name of the file to write (bjj: what about overlapping instances of HydroDyn? (perhaps OutRootName should take care of this!)
   HDO_InitData%NWaveElev   = Waves_InitData%NWaveElev

   CALL HDOut_Init( HDO_InitData, HD_Data%HDOut_Data, HD_Data%FltPtfm_data, ErrStat )   


   !-------------------------------------------------------------------------------------------------
   ! Fill the variables that return the discretization 
   !-------------------------------------------------------------------------------------------------

   ALLOCATE( HD_AllMarkers%Substructure( NMarkers ), STAT=ErrStat )
   IF ( ErrStat /= 0 ) THEN
      CALL ProgAbort( ' Error allocating substructure markers in HydroDyn.', TrapErrors = .TRUE. )
      RETURN
   END IF
   
      ! These are the nodes that make up the discretization (values are per unit length)
   DO I = 1,NMarkerElements
   !bjj: ask jmj if this is the same!!!
      HD_AllMarkers%Substructure(I)%Position = HD_Data%Node(I)%Position - HD_ConfigMarkers%Substructure%Position       ! position relative to the substructure marker 
         
   END DO
   
   
      ! These are the fixed-point markers
   DO I = NMarkerElements+1,NMarkers
   
      HD_AllMarkers%Substructure(I)%Position = HD_ConfigMarkers%Substructure%Position
      
   END DO
   

   
END SUBROUTINE HD_Init
!====================================================================================================
SUBROUTINE HD_Terminate ( HD_Data, ErrStat )
! This public subroutine is called at program termination.  It deallocates variables and closes files.
!----------------------------------------------------------------------------------------------------  
   
      ! Passed variables   
   
   TYPE( HD_DataType ), INTENT( INOUT )   :: HD_Data
   INTEGER,             INTENT(   OUT )   :: ErrStat


      ! Internal variables

   LOGICAL                                :: Err


      ! Initialize error status code

   ErrStat = 0
   Err     = .FALSE.

   
   !-------------------------------------------------------------------------------------------------
   ! Write the last line to our output file
   ! bjj: this may produce strange results?
   !-------------------------------------------------------------------------------------------------
   CALL HDOut_WriteOutputs( HD_Data%HDOut_Data, ErrStat )   
   !non-zero error just indicates there is no Output data to write
   ErrStat = 0
   
   !-------------------------------------------------------------------------------------------------
   ! Clean up the modules we have used
   !-------------------------------------------------------------------------------------------------

   CALL Waves_Terminate(HD_Data%Waves_data,   ErrStat)
   IF ( ErrStat /= 0 ) Err = .TRUE.  
   
   CALL FP_Terminate(  HD_Data%FltPtfm_data, ErrStat)
   IF ( ErrStat /= 0 ) Err = .TRUE.
      
   CALL FB_Terminate()
      
   CALL HDOut_Terminate( HD_Data%HDOut_Data, ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.

      
   !-------------------------------------------------------------------------------------------------
   ! Deallocate HD memory               
   !-------------------------------------------------------------------------------------------------
      
   IF ( ALLOCATED( HD_Data%Node ) ) DEALLOCATE( HD_Data%Node, STAT = ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.
      
   !-------------------------------------------------------------------------------------------------
   ! Close any files
   !-------------------------------------------------------------------------------------------------
   CALL HDOut_CloseSum( HD_Data%UnSum, ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.


   !-------------------------------------------------------------------------------------------------
   ! Reset any initialization values
   !-------------------------------------------------------------------------------------------------
   HD_Data%LastCalcTime = 0.0
   HD_Data%NElements    = 0.0
   HD_Data%MaxDiam      = 0.0
   HD_Data%StrctType    = Unknown_Type    
   HD_Data%UnSum        = -1

   
   !-------------------------------------------------------------------------------------------------
   ! Return ErrStat = 1 if any error occurred while cleaning up
   !-------------------------------------------------------------------------------------------------
   IF ( Err ) ErrStat = 1
   
   
END SUBROUTINE HD_Terminate
!====================================================================================================
END MODULE HydroDyn
