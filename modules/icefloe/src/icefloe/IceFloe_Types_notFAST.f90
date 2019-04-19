!*********************************************************************************************************************************
MODULE IceFloe_Types
! this is a simplified version for use with IceFloe with codes other than FAST
!---------------------------------------------------------------------------------------------------------------------------------
use precision
IMPLICIT NONE
! =========  IceFloe_ParameterType  =======
  TYPE, PUBLIC :: IceFloe_ParameterType
    REAL(ReKi) , DIMENSION(:,:), ALLOCATABLE  :: loadSeries      ! - [precalculated time series of ice loads for each leg]
    REAL(ReKi)  :: iceVel      ! ice floe velocity [m/s]
    REAL(ReKi)  :: iceDirection      ! ice floe direction [degrees]
    REAL(ReKi)  :: minStrength      ! minimum dynamic ice strength [Pa]
    REAL(ReKi)  :: minStrengthNegVel      ! minimum dynamic ice strength for negative velocity [Pa]
    REAL(ReKi)  :: defaultArea       ! width to use for coupled crushing
    REAL(ReKi)  :: crushArea      ! cross sectional area of ice against tower [m^2]
    REAL(ReKi)  :: coeffStressRate      ! coefficient to calc stress rate from relative vellocity [Pa/m]
    REAL(ReKi)  :: C(4)      ! coefficient of cubic transition curve for negative stress rates [-]
    REAL(ReKi)  :: dt        ! time step [sec]
    REAL(ReKi)  :: rampTime  ! ramp time for forces at beginning of simulation [sec]
    REAL(ReKi) , DIMENSION(:), ALLOCATABLE  :: legX      ! - [x position of each leg relative to structure center]
    REAL(ReKi) , DIMENSION(:), ALLOCATABLE  :: legY      ! - [y position of each leg relative to structure center]
    REAL(ReKi) , DIMENSION(:), ALLOCATABLE  :: ks      ! - [shelter factor due to upstream leg]
    INTEGER(IntKi)  :: numLegs      ! Number of tower legs (=1 for monopile) [-]
    INTEGER(IntKi)  :: iceType      ! Type of ice Floe: flex, crush, etc. [-]
    INTEGER(IntKi)  :: logUnitNum      ! Unit number for log file [-]
    LOGICAL  :: singleLoad      ! Flag for load application at single point vs multiple legs [-]
    LOGICAL  :: initFlag      ! Flag for successful initialization [-]
  END TYPE IceFloe_ParameterType
! =======================

! Some constants that would otherwise come from NREL code base
   real(ReKi), parameter :: Pi  = 3.141592654_ReKi
   real(ReKi), parameter :: D2R = Pi/180.0_ReKi
   real(ReKi), parameter :: R2D = 180.0_ReKi/Pi

END MODULE IceFloe_Types
