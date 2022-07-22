!==================================================================================================================================   
MODULE TurbSim_Types

use NWTC_Library

   TYPE(ProgDesc)            :: TurbSim_Ver = ProgDesc( 'TurbSim', '', '' )

   LOGICAL,        PARAMETER :: MVK         = .FALSE.                       ! This parameter has been added to replace the NON-STANDARD compiler directive previously used
   LOGICAL,        PARAMETER :: PeriodicY   = .FALSE. !.TRUE.


   INTEGER(IntKi), PARAMETER :: MaxMsgLen = 1024 ! Maximum length of error messages

      ! Valid spectral models (i.e., turbulence models; values of TurbModel_ID)
   INTEGER(IntKi), PARAMETER :: SpecModel_NONE    =  0  ! No turbulence
   INTEGER(IntKi), PARAMETER :: SpecModel_IECKAI  =  1  ! IEC Kaimal
   INTEGER(IntKi), PARAMETER :: SpecModel_IECVKM  =  2  ! IEC von Karman 
   INTEGER(IntKi), PARAMETER :: SpecModel_GP_LLJ  =  3  ! Great Plains Low-Level Jet
   INTEGER(IntKi), PARAMETER :: SpecModel_NWTCUP  =  4  ! NWTC (upwind)
   INTEGER(IntKi), PARAMETER :: SpecModel_SMOOTH  =  5  ! Risoe Smooth-Terrain   
   INTEGER(IntKi), PARAMETER :: SpecModel_WF_UPW  =  6  ! Wind Farm Upwind
   INTEGER(IntKi), PARAMETER :: SpecModel_WF_07D  =  7  ! Wind Farm  7 rotor diameters downwind
   INTEGER(IntKi), PARAMETER :: SpecModel_WF_14D  =  8  ! Wind Farm 14 rotor diameters downwind
   INTEGER(IntKi), PARAMETER :: SpecModel_TIDAL   =  9  ! Tidal (Hydro)
   INTEGER(IntKi), PARAMETER :: SpecModel_RIVER   = 10  ! River (Hydro)
   INTEGER(IntKi), PARAMETER :: SpecModel_API     = 11  ! API
   INTEGER(IntKi), PARAMETER :: SpecModel_MODVKM  = 12  ! user-specified scaling in von Karman model
   INTEGER(IntKi), PARAMETER :: SpecModel_USRVKM  = 13  ! user-specified scaling in von Karman model
   INTEGER(IntKi), PARAMETER :: SpecModel_USER    = 14  ! User-defined spectra from file
   INTEGER(IntKi), PARAMETER :: SpecModel_TimeSer = 15  ! time series input from file
   
      ! Spatial Coherence Models (SCMod)
   INTEGER(IntKi), PARAMETER :: CohMod_NONE       = 0   ! no additional spatial coherence
   INTEGER(IntKi), PARAMETER :: CohMod_GENERAL    = 1   ! General spatial coherence model using parameters input from file
   INTEGER(IntKi), PARAMETER :: CohMod_IEC        = 2   ! Spatial coherence specified by IEC standard
   INTEGER(IntKi), PARAMETER :: CohMod_API        = 3   ! Spatial coherence specified by API standard
   
   
      ! IEC turbulence types (IEC_WindType) 
      ! bjj: note that EWM models *MUST* directly follow ETM, and EWM models must be at the end
   INTEGER(IntKi), PARAMETER :: IEC_NTM           = 1   ! Number to indicate the IEC Normal Turbulence Model
   INTEGER(IntKi), PARAMETER :: IEC_ETM           = 2   ! Number to indicate the IEC Extreme Turbulence Model
   INTEGER(IntKi), PARAMETER :: IEC_EWM1          = 3   ! Number to indicate the IEC Extreme Wind speed Model (  1-year)
   INTEGER(IntKi), PARAMETER :: IEC_EWM50         = 4   ! Number to indicate the IEC Extreme Wind speed Model ( 50-year)
   INTEGER(IntKi), PARAMETER :: IEC_EWM100        = 5   ! Number to indicate the IEC Extreme Wind speed Model (100-year)
   
   
      ! distinct output file formats (WrFile()) (listed by extension)
   INTEGER(IntKi), PARAMETER :: FileExt_BTS       =  1  ! .bts file         : AeroDyn FF data (binary) [WrADFF]
   INTEGER(IntKi), PARAMETER :: FileExt_WND       =  2  ! .wnd file         : BLADED FF data (binary)  [WrBLFF]
   INTEGER(IntKi), PARAMETER :: FileExt_HH        =  3  ! .hh file          : AeroDyn HH data (formatted) [WrADHH]
   INTEGER(IntKi), PARAMETER :: FileExt_BIN       =  4  ! .bin file         : binary HH turbulence parameters [WrBHHTP]
   INTEGER(IntKi), PARAMETER :: FileExt_DAT       =  5  ! .dat file         : formatted HH turbulence parameters [WrFHHTP]
   INTEGER(IntKi), PARAMETER :: FileExt_UVW       =  6  ! .u, .v, .w files  : formatted FF data (Traditional SNLWIND-3D format) [WrFMTFF]
   INTEGER(IntKi), PARAMETER :: FileExt_CTS       =  7  ! .cts file         : coherent turbulence
   INTEGER(IntKi), PARAMETER :: FileExt_TWR       =  8  ! .twr file         : AeroDyn tower data (binary)
   INTEGER(IntKi), PARAMETER :: FileExt_HAWC      =  9  ! -u.bin, -v.bin, -w.bin .hawc files: binary HAWC FF data [WrHAWCFF]
   INTEGER(IntKi), PARAMETER :: NumFileFmt        =  9  ! TOTAL number of output file formats (used to dimension array)
   
      ! other parameters:
   REAL(ReKi),     PARAMETER :: ZJetMax_UB        = 490.0_ReKi    ! upper bound on height where jet maximum occurs
   REAL(ReKi),     PARAMETER :: ZJetMax_LB        =  70.0_ReKi    ! lower bound on height where jet maximum occurs
   REAL(ReKi),     PARAMETER :: profileZmax       = 140.          ! Upper height limit for extrapolating GP_LLJ profiles of ustar and zl
   REAL(ReKi),     PARAMETER :: profileZmin       =  50.          ! Lower height limit for extrapolating GP_LLJ profiles of ustar and zl
   REAL(ReKi),     PARAMETER :: Omega             = 7.292116E-05  ! Angular speed of rotation of the earth (rad/s)
   REAL(ReKi),     PARAMETER :: Tolerance         = 0.0001        ! The largest difference between two numbers that are assumed to be equal
                   
   CHARACTER(1),   PARAMETER :: Comp (3) = (/ 'u', 'v', 'w' /)   ! The names of the wind components
   
   
   type :: RandNum_ParameterType
   
      integer(IntKi)                  :: pRNG
      INTEGER(IntKi)                  :: RandSeed   (3)                           ! The array that holds the initial random seeds for the 3 components.
      INTEGER(IntKi),    ALLOCATABLE  :: RandSeedAry(:)                           ! The array that holds the random seeds.
      CHARACTER(  6)                  :: RNG_type                                 ! Type of Random Number Generator to use
            
   end type RandNum_ParameterType

   type :: RandNum_OtherStateType
      INTEGER(IntKi), ALLOCATABLE     :: NextSeed   (:)                           ! The array that holds the next random seed for the 3 components.            
   end type RandNum_OtherStateType
                   
   
   TYPE     :: CohStr_ParameterType   

      REAL(ReKi)                   :: CTLy                                     ! Fractional location of tower centerline from right (looking downwind) to left side of the dataset.
      REAL(ReKi)                   :: CTLz                                     ! Fractional location of hub height from the bottom of the dataset.
      REAL(ReKi)                   :: CTStartTime                              ! Minimum time to add coherent structures
      REAL(ReKi)                   :: DistScl                                  ! Disturbance scale for AeroDyn coherent turbulence events
   
                  
      CHARACTER(200)               :: CTEventPath                              ! String used to store the name of the coherent event definition file
      CHARACTER(200)               :: CTEventFile                              ! String used to store the name of the coherent event definition file
      CHARACTER(  3)               :: CTExt                                    ! String used to determine the type of coherent structures ("dns" or "les")

   END TYPE CohStr_ParameterType
   

   
   type :: Grid_ParameterType
      
      REAL(ReKi)                   :: GridHeight                               ! Grid height
      REAL(ReKi)                   :: GridRes_Z                                ! Distance between two consecutive vertical points on the grid (Vertical resolution)
      INTEGER(IntKi)               :: NumGrid_Z                                ! Grid dimension. (in vertical direction)

      REAL(ReKi)                   :: GridWidth                                ! Grid width.
      REAL(ReKi)                   :: GridRes_Y                                ! Distance between two consecutive horizontal points on the grid (Horizontal resolution)
      INTEGER(IntKi)               :: NumGrid_Y                                ! Grid dimension. (in horizontal direction)

      INTEGER(IntKi)               :: NPoints                                  ! Number of points being simulated.                        
      INTEGER(IntKi)               :: NPacked                                  ! Number of entries stored in the packed version of the symmetric matrix of size NPoints by NPoints
      
      REAL(ReKi)                   :: Zbottom                                  ! The height of the lowest point on the grid (before tower points are added), equal to Z(1)
      REAL(ReKi)                   :: RotorDiameter                            ! The assumed diameter of the rotor
      
      INTEGER(IntKi)               :: HubIndx                                  ! Index that tells where the hub point is in the V matrix
      
      REAL(ReKi)                   :: HubHt                                    ! Hub height.
      LOGICAL                      :: HubOnGrid                                ! Flag to indicate if the hub is on the regular grid (true) or if an extra point must be added (false)
      LOGICAL                      :: ExtraTwrPT                               ! Flag to indicate if the tower is on the regular grid or if an extra point must be added
            
      
      REAL(ReKi),    ALLOCATABLE   :: Y          (:)                           ! The lateral locations of the points (NPoints).
      REAL(ReKi),    ALLOCATABLE   :: Z          (:)                           ! The vertical locations of the points (NPoints).
      
      INTEGER(IntKi),ALLOCATABLE   :: GridPtIndx  (:)                          ! size is (NumGrid_Y * NumGrid_Z): The indices into the velocity array, indicating the points in the cartesian grid for output. Previously was assumed to be the first NumGrid_Y * NumGrid_Z points; now necessary because user-defined time series data may specifiy points on the grid.
      INTEGER(IntKi),ALLOCATABLE   :: TwrPtIndx   (:)                          ! size is number of tower points: The indices into the velocity array, indicating the points to be put in the tower file for output. Previously was assumed to be the last NumTower points; now necessary because user-defined time series data may specifiy points on the tower.

                        
      REAL(ReKi)                   :: AnalysisTime                             ! Analysis Time. (amount of time for analysis, allows user to perform analysis using one time length, but output UsableTime            
      REAL(ReKi)                   :: UsableTime                               ! Usable time.  Program adds GridWidth/MeanHHWS if not specified as "ALL"  AnalysisTime in input file.
      REAL(ReKi)                   :: TimeStep                                 ! Time step.
      REAL(ReKi),    ALLOCATABLE   :: Freq       (:)                           ! The array of frequencies (NumFreq).
      INTEGER(IntKi)               :: NumFreq                                  ! Number of frequencies (=NumSteps/2).
      INTEGER(IntKi)               :: NumSteps                                 ! Number of time steps for the FFT.
      INTEGER(IntKi)               :: NumOutSteps                              ! Number of output time steps.
                  
      LOGICAL                      :: Periodic                                 ! Flag to indicate that output files must contain exactly one full (time) period
      
   end type Grid_ParameterType
   
   
   type IEC_ParameterType
      INTEGER(IntKi)               :: IECedition                               ! The edition number of the IEC 61400-1 standard that is being used (determines the scaling)
      INTEGER(IntKi)               :: IECstandard                              ! The standard number (x) of the IEC 61400-x that is being used
      INTEGER(IntKi)               :: IEC_WindType                             ! Number to indicate the IEC wind type
      INTEGER(IntKi)               :: ScaleIEC                                 ! Switch to indicate if turbulence should be scaled to target value; 0 = NO scaling; 1 = scale based on hub; 2 = scale each point individually
      
      REAL(ReKi)                   :: Lambda         (3)                       ! IEC turbulence scale parameter: defined as wavelength where the non-dimensional power spectral density fS(f)/sigma^2 == 0.05 [m]
      REAL(ReKi)                   :: SigmaIEC       (3)                       ! IEC target standard deviation.
      REAL(ReKi)                   :: IntegralScale  (3)                       ! IEC integral scales (s)
      REAL(ReKi)                   :: LC                                       ! IEC coherency scale parameter
      REAL(ReKi)                   :: SigmaSlope                               ! Slope used with IEC models to determine target sigma and turbulent intensity
      REAL(ReKi)                   :: TurbInt                                  ! IEC target Turbulence Intensity 
      REAL(ReKi)                   :: TurbInt15                                ! Turbulence Intensity at hub height with a mean wind speed of 15 m/s
      REAL(ReKi)                   :: ETMc                                     ! The c parameter in IEC ETM, 61400-1, Ed 3. Section 6.3.2.3, Eq. 19.  Variable per last sentence in section 7.4.1
      REAL(ReKi)                   :: Vave                                     ! The IEC Vave for ETM
      REAL(ReKi)                   :: Vref                                     ! The IEC Vref for ETM
      REAL(ReKi)                   :: PerTurbInt                               ! Percent Turbulence Intensity
                  
      LOGICAL                      :: NumTurbInp                               ! Flag to indicate if turbulence is user-specified (as opposed to IEC standard A, B, or C)
         
      CHARACTER(  1)               :: IECTurbC                                 ! IEC turbulence characteristic
      CHARACTER(  1)               :: IECTurbE                                 ! IEC Extreme turbulence class
      CHARACTER( 35)               :: IEC_WindDesc                             ! The description of the IEC wind type      
      CHARACTER( 25)               :: IECeditionStr                            ! description of the IEC standard being used 
   end type IEC_ParameterType
   
   
   type Meteorology_ParameterType
   
      INTEGER(IntKi)               :: TurbModel_ID                             ! Integer value of spectral model (see SpecModel enum)      
      LOGICAL                      :: KHtest                                   ! Flag to indicate that turbulence should be extreme, to demonstrate effect of KH billows
      CHARACTER(  3)               :: WindProfileType                          ! The wind profile type
      CHARACTER(  6)               :: TurbModel                                ! Turbulence model
      CHARACTER( 50)               :: TMName                                   ! Turbulence model name.
      
      REAL(ReKi)                   :: Fc                                       ! Coriolis parameter in units (1/sec)
     !REAL(ReKi)                   :: h                                        ! Boundary layer depth
      REAL(ReKi)                   :: RICH_NO                                  ! Gradient Richardson number
      REAL(ReKi)                   :: Z0                                       ! Surface roughness length, meters
      REAL(ReKi)                   :: ZI                                       ! Mixing layer depth
      REAL(ReKi)                   :: Latitude                                 ! The site latitude in radians
      REAL(ReKi)                   :: L                                        ! M-O length
      REAL(ReKi)                   :: ZL                                       ! A measure of stability
      REAL(ReKi)                   :: PLExp                                    ! Rotor disk power law exponent
      REAL(ReKi)                   :: Ustar                                    ! Shear or friction velocity (m/s) -- rotor-disk average
      REAL(ReKi)                   :: UstarDiab                                ! The diabatic ustar value
      REAL(ReKi)                   :: UstarOffset                              ! A scaling/offset value used with the Ustar_profile to ensure that the mean hub u'w' and ustar inputs agree with the profile values
      REAL(ReKi)                   :: UstarSlope                               ! A scaling/slope value used with the Ustar_profile to ensure that the mean hub u'w' and ustar inputs agree with the profile values
      REAL(ReKi)                   :: ZLoffset                                 ! An offset to align the zl profile with the mean zl input parameter
   
      REAL(ReKi)                   :: RefHt                                    ! Height for reference wind speed.
      REAL(ReKi)                   :: URef                                     ! The input wind speed at the reference height.  (Added by M. Buhl for API profiles)
      
      REAL(ReKi), ALLOCATABLE      :: ZL_profile(:)                            ! A profile of z/l (measure of stability with height)
      REAL(ReKi), ALLOCATABLE      :: Ustar_profile(:)                         ! A profile of ustar (measure of friction velocity with height)
      !REAL(ReKi)                   :: TurbIntH20                               ! Turbulence intensity used for HYDRO module.

      
      REAL(ReKi)                   :: HH_HFlowAng                              ! Horizontal flow angle at the hub (may be different than HFlowAng if using direction profile).
      REAL(ReKi)                   :: HH_VFlowAng                              ! Vertical flow angle at the hub (may be different than VFlowAng if using vertical angle profile (i.e., user-defined time-series)).
      REAL(ReKi)                   :: HFlowAng                                 ! Horizontal flow angle
      REAL(ReKi)                   :: VFlowAng                                 ! Vertical flow angle
                  
      
         ! coefficients for velocity and direction profiles (currently used with jet profiles only)
      REAL(ReKi)                   :: ChebyCoef_WS(11)                         ! The Chebyshev coefficients for wind speed
      REAL(ReKi)                   :: ChebyCoef_WD(11)                         ! The Chebyshev coefficients for wind direction

      
      REAL(ReKi)                   :: ZJetMax                                  ! The height of the jet maximum (m)
      REAL(ReKi)                   :: UJetMax                                  ! The (horizontal) wind speed at the height of the jet maximum (m/s)
      
      
         ! Coherence
      REAL(ReKi)                   :: COHEXP                                   ! Coherence exponent for general spatial coherence model
      REAL(ReKi)                   :: InCDec     (3)                           ! Contains the coherence decrements for general spatial coherence model
      REAL(ReKi)                   :: InCohB     (3)                           ! Contains the coherence b/L (offset) parameters for general spatial coherence model
      
      INTEGER(IntKi)               :: SCMod      (3)                           ! SCMod_u, SCMod_v, and SCMod_w: switches determining which coherence model to use
      
      LOGICAL                      :: IsIECModel                               ! Flag to determine if we're using IEC scaling (coherence, etc)
      
      
         ! Scaling
      REAL(ReKi)                   :: PC_UW                                    ! u'w' cross-correlation coefficient
      REAL(ReKi)                   :: PC_UV                                    ! u'v' cross-correlation coefficient
      REAL(ReKi)                   :: PC_VW                                    ! v'w' cross-correlation coefficient
         
      LOGICAL                      :: UVskip                                   ! Flag to determine if UV cross-feed term should be skipped or used
      LOGICAL                      :: UWskip                                   ! Flag to determine if UW cross-feed term should be skipped or used
      LOGICAL                      :: VWskip                                   ! Flag to determine if VW cross-feed term should be skipped or used
      
      
         ! user-defined profiles (also used with UsrVKM model):      
      INTEGER(IntKi)               :: NumUSRz                                  ! Number of heights defined in the user-defined profiles.
      REAL(ReKi), ALLOCATABLE      :: USR_Z        (:)                         ! Heights of user-specified variables
      REAL(ReKi), ALLOCATABLE      :: USR_U        (:)                         ! User-specified total wind speed, varying with height
      REAL(ReKi), ALLOCATABLE      :: USR_WindDir  (:)                         ! User-specified wind direction profile, varying with height
      REAL(ReKi), ALLOCATABLE      :: USR_Sigma    (:)                         ! User-specified standard deviation of the wind speed components (isotropic), varying with height
      REAL(ReKi), ALLOCATABLE      :: USR_L        (:)                         ! User-specified von Karman length scale, varying with height
      REAL(ReKi)                   :: USR_StdScale (3)                         ! Scaling for the user-specified standard deviation
            
   end type Meteorology_ParameterType

   
   TYPE UserTSSpec_ParameterType
                     
      integer(intKi)               :: nComp                                   ! number of velocity components in the file (1=u; 2=u&v; 3=u,v,w)
      integer(intKi)               :: nFreq                                   ! number of frequencies in the calculated spectra
      real(dbKi)                   :: DelF                                    ! delta frequenc of the calculated spectra (same as f(1) in double precision)
      integer(intKi)               :: nPoints                                 ! number of points in the time series input
      integer(intKi)               :: RefPtID                                 ! Index of the reference point (1-nPoints)
      integer(intKi)               :: nTimes                                  ! number of rows in the time series input
      real(reki),     allocatable  :: pointyi (:)                             ! y position where each time series was input; size: nPoints
      real(reki),     allocatable  :: pointzi (:)                             ! z position (height) where each time series was input; size: nPoints
      real(reki),     allocatable  :: t(:)
      real(reki),     allocatable  :: v(:,:,:)                                ! velocity time series; size: nTimes, nPoints, { 2 if .not. containsW | 3 otherwise }
      
      real(reKi),     allocatable  :: meanU(:,:)                              ! mean velocity; size: nPoints, nComp [m/s]
      real(reKi),     allocatable  :: meanDir(:)                              ! mean horizontal direction; size: nPoints [degrees]
      real(reKi),     allocatable  :: meanVAng(:)                             ! mean vertical angle; size: nPoints [degrees]      
      real(reKi),     allocatable  :: S(:,:,:)                                ! spectra;   size: nFreq, nPoints, nComp
      real(reKi),     allocatable  :: f(:)                                    ! frequency; size: nFreq [Hz]
      real(reKi),     allocatable  :: phaseAngles(:,:,:)
      
      INTEGER(IntKi)               :: TurbModel_ID   = SpecModel_NONE          ! Integer value of spectral model (see SpecModel enum) for filling in high-frequency content; not all SpecModel values are valid here
      
   END TYPE UserTSSpec_ParameterType
   
   
   type TurbSim_ParameterType
      LOGICAL                          :: WrFile(NumFileFmt)          ! Flag to determine which output files should be generated 
      INTEGER                          :: US  = -1                    ! I/O unit for summary file.
              
      
      CHARACTER(200)                   :: DescStr                     ! String used to describe the run (and the first line of the summary file)
      CHARACTER(197)                   :: RootName                    ! Root name of the I/O files.
      TYPE(RandNum_ParameterType)      :: RNG                         ! parameters for random numbers p_RandNum
      TYPE(Grid_ParameterType)         :: grid                        ! parameters for TurbSim (specify grid/frequency size)
      TYPE(Meteorology_ParameterType)  :: met                         ! parameters for TurbSim 
      TYPE(IEC_ParameterType)          :: IEC                         ! parameters for IEC models
      TYPE(UserTSSpec_ParameterType)   :: usr                         ! parameters for user spectra or time-series input 
      TYPE(CohStr_ParameterType)       :: CohStr                      ! parameters for coherent structures 
      
   !bjj: there probably won't be a need for this later...      
      REAL(ReKi)                       :: UHub                        ! Hub-height (total) wind speed (m/s)
      
      
   end type
   
   
   
   
   
END MODULE TurbSim_Types
!==================================================================================================================================   
