MODULE read_wind_farm_parameter_data
    IMPLICIT NONE
    REAL     ::  HubHt
    REAL     ::  RotorR
    INTEGER  ::  NumWT       ! The total number of wind turbines
    REAL     ::  Uambient
    REAL     ::  TI
    INTEGER  ::  ppR
    REAL     ::  Domain_R
    REAL     ::  Domain_X
    INTEGER  ::  Mstl        ! Meandering_simulation_time_length
    INTEGER  ::  Mmt         ! Meandering_Moving_time
    INTEGER  ::  smooth_flag ! smoothed out upstream wake profile flag
    REAL     ::  WFLowerBd   ! The lower bound height of the wind file (m)
    REAL     ::  Winddir     ! The ambient wind direction
    INTEGER  ::  Tinfluencer ! The max number of upstream turbines that affects a downstream turbine (-)
    CHARACTER(1024):: DWM_exe_name

    REAL(8), ALLOCATABLE  ::  Xcoordinate (:)       ! wind turbine x location
    REAL(8), ALLOCATABLE  ::  Ycoordinate (:)       ! wind turbine y location
        
END MODULE read_wind_farm_parameter_data
    
!------------------------------------------------------------
MODULE wind_farm_geometry_data
    IMPLICIT NONE
    INTEGER,ALLOCATABLE   ::    turbine_sort(:)            ! the array that stores the order of the turbines from upstream to downstream
    INTEGER,ALLOCATABLE   ::    TurbineInfluenceData(:,:)
    REAL,ALLOCATABLE      ::    wake_sector_angle_array(:)
    REAL(8)               ::    xwind
    REAL(8)               ::    ywind
    INTEGER               ::    scale_factor
    REAL                  ::    Pi
    REAL,ALLOCATABLE      ::    length(:)        ! projected length on the wind direction vector

END MODULE wind_farm_geometry_data
    
!------------------------------------------------------------
MODULE DWM_init_data
    CHARACTER(1024)              :: OutFileRoot
    CHARACTER(1024)              :: InputFile
END MODULE DWM_init_data
