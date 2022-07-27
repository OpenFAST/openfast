MODULE DWM_driver_wind_farm_sub
    USE NWTC_Library
    USE VersionInfo
    IMPLICIT NONE

!PUBLIC SUBROUTINES
    PUBLIC  :: read_wind_farm_parameter
    PUBLIC  :: write_parameter_to_file
    PUBLIC  :: wind_farm_geometry
    PUBLIC  :: projected_length
    PUBLIC  :: AngleBetweenVectors
    PUBLIC  :: turbine_position_LorR
    PUBLIC  :: Driver_init
    PUBLIC  :: delete_temp_files
    PUBLIC  :: Check_DWM_parameter
    
    
CONTAINS
!-------------------------------------------------------------------
SUBROUTINE read_wind_farm_parameter(PriFile)
!...................................................................
! This subroutine is to read the wind farm parameter files
! Including the number of rows of the wind farm, the number of turbine in each row
!   and the turbine spacing
!...................................................................
    USE read_wind_farm_parameter_data
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: PriFile
    INTEGER  ::  UnIn = 0
    INTEGER  ::  UnEc = -1
    INTEGER  ::  I
    CHARACTER(1024) :: DWM_Title,comment
    INTEGER  ::  ErrStat = 0
    CHARACTER(ErrMsgLen) :: ErrMsg
    INTEGER(4) :: IOS

    CALL OpenFInpFile ( UnIn, PriFile, ErrStat, ErrMsg )    

    READ (UnIn,'(//,A,/)',IOSTAT=IOS)  DWM_Title                      ! read the words (title)
    CALL CheckIOS( IOS, PriFile, 'file title', StrType )

    READ (UnIn,'(A)',IOSTAT=IOS)  Comment                             ! read the words (comment)
    CALL CheckIOS( IOS, PriFile, 'simulation control parameters comment', StrType )
    
      ! Read in the hub height
    CALL ReadVar ( UnIn, PriFile, HubHt, 'HubHt', 'The hub height (m)', ErrStat, ErrMsg, UnEc )
    IF ( ErrStat /= 0 ) RETURN
    
      ! Read in the rotor radius
    CALL ReadVar ( UnIn, PriFile, RotorR, 'RotorR', 'The Rotor radius (m)', ErrStat, ErrMsg, UnEc )
    IF ( ErrStat /= 0 ) RETURN

      ! Read in the total number of wind turbines
    CALL ReadVar ( UnIn, PriFile, NumWT, 'NumWT', 'The total number of wind turbines (-)', ErrStat, ErrMsg, UnEc )
    IF ( ErrStat /= 0 ) RETURN
       
      ! Read in the ambient wind velocity
    CALL ReadVar ( UnIn, PriFile, Uambient, 'Uambient', 'The ambient wind velocity (m/s)', ErrStat, ErrMsg, UnEc )
    IF ( ErrStat /= 0 ) RETURN
    
      ! Read in the ambient TI
    CALL ReadVar ( UnIn, PriFile, TI, 'TI', 'TI for first turbine (%)', ErrStat, ErrMsg, UnEc )
    IF ( ErrStat /= 0 ) RETURN
    
      ! Read in the radial domain size
    CALL ReadVar ( UnIn, PriFile, ppR, 'ppR', 'Point per R resolution (-)', ErrStat, ErrMsg, UnEc )
    IF ( ErrStat /= 0 ) RETURN
    
      ! Read in the point per R resolution
    CALL ReadVar ( UnIn, PriFile, Domain_R, 'Domain_R', 'Radial domain size (R)', ErrStat, ErrMsg, UnEc )
    IF ( ErrStat /= 0 ) RETURN
    
      ! Read in the point per R resolution
    CALL ReadVar ( UnIn, PriFile, Domain_X, 'Domain_X', 'Longitudinal domain size (R)', ErrStat, ErrMsg, UnEc )
    IF ( ErrStat /= 0 ) RETURN
    
      ! Read in the simulation time length of the meandering wake model
    CALL ReadVar ( UnIn, PriFile, Mstl, 'Meandering_simulation_time_length', 'The length of the simulation time in the meandering wake model (-)', ErrStat, ErrMsg, UnEc )
    IF ( ErrStat /= 0 ) RETURN
    
      ! Read in the moving time length of the meandering wake model
    CALL ReadVar ( UnIn, PriFile, Mmt, 'Meandering_Moving_time', 'The length of the moving time in the meandering wake model (-)', ErrStat, ErrMsg, UnEc )
    IF ( ErrStat /= 0 ) RETURN
    
      ! Read in the lower bound height of the wind file
    CALL ReadVar ( UnIn, PriFile, WFLowerBd, 'WFLowerBd', 'The lower bound height of the wind file (m)', ErrStat, ErrMsg, UnEc )
    IF ( ErrStat /= 0 ) RETURN
   
      ! Read in the ambient wind direction
    CALL ReadVar ( UnIn, PriFile, Winddir, 'Winddir', 'The ambient wind direction (degree)', ErrStat, ErrMsg, UnEc )
    IF ( ErrStat /= 0 ) RETURN
    
      ! Read in the DWM_FAST.exe file rootname
    CALL ReadVar ( UnIn, PriFile, DWM_exe_name, 'DWM_exe_name', 'The file rootname of the DWM-FAST program', ErrStat, ErrMsg, UnEc )

    ALLOCATE (Xcoordinate(NumWT))
    ALLOCATE (Ycoordinate(NumWT))
    
    CALL ReadCom( UnIn, PriFile, 'Coordinate table headers', ErrStat, ErrMsg)
    IF ( ErrStat /= 0 ) RETURN
       
    DO i = 1, NumWT
        IF (ErrStat == 0) THEN
           READ(UnIn,*,IOSTAT=ErrStat) Xcoordinate(i), Ycoordinate(i)
        END IF
    END DO
    
    Tinfluencer = 1
    !Tinfluencer = 10
    
END SUBROUTINE read_wind_farm_parameter

!-------------------------------------------------------------------
SUBROUTINE Check_DWM_parameter()
!...................................................................
! This subroutine is to check if the input parameters are valid
!...................................................................
    USE read_wind_farm_parameter_data

    IF (HubHt < 0.0) THEN
        PRINT*, 'WARNING: HubHt should be positive'
        CALL EXIT
    END IF
    
    IF (RotorR < 0.0) THEN
        PRINT*, 'WARNING: RotorR should be positive'
        CALL EXIT
    END IF
    IF (RotorR > HubHt) THEN
        PRINT*, 'WARNING: RotorR should be smaller than HubHt'
        CALL EXIT
    END IF
    
    IF (NumWT < 1.0) THEN
        PRINT*, 'WARNING: NumWT should be positive'
        CALL EXIT
    END IF
    
    IF (Uambient < 0.0) THEN
        PRINT*, 'WARNING: Uambient should be positive'
        CALL EXIT
    END IF
        
    IF (TI < 0.0) THEN
        PRINT*, 'WARNING: TI should be positive'
        CALL EXIT
    END IF
    
    IF (ppR < 1.0) THEN
        PRINT*, 'WARNING: ppR should be positive'
        CALL EXIT
    END IF
    
    IF (Domain_R < 0.0) THEN
        PRINT*, 'WARNING: Domain_R should be positive'
        CALL EXIT
    END IF
    
    IF (Domain_X < 0.0) THEN
        PRINT*, 'WARNING: Domain_X should be positive'
        CALL EXIT
    END IF

END SUBROUTINE Check_DWM_parameter

!-------------------------------------------------------------------
SUBROUTINE write_parameter_to_file()
!...................................................................
! This subroutine is to write the turbine spacing from the input file
! to a bin file, which enables the DWM to read the turbine spacing at every
! instances of running
!...................................................................

    USE read_wind_farm_parameter_data
    IMPLICIT NONE
    integer   :: Un
    
    !CALL GetNewUnit(Un)
    Un = 10
    
    OPEN(unit = Un, status='replace',file='DWM_parameter.bin',form='unformatted')
    WRITE(Un)   HubHt,RotorR,NumWT,Uambient,TI,Domain_R,Domain_X,ppR,Mstl,Mmt,WFLowerBd,Winddir,Tinfluencer
    CLOSE(Un)
    
    OPEN(unit = Un, status='replace',file='wind_farm_coordinate.bin',form='unformatted')
    WRITE(Un)   Xcoordinate,Ycoordinate
    CLOSE(Un)
    
END SUBROUTINE write_parameter_to_file

!-------------------------------------------------------------------
SUBROUTINE wind_farm_geometry()
!...................................................................
! This subroutine is to calculate the wind direction and sort the wind 
! turbines as the order from upstream to downstream
! output the angles between turbines
!...................................................................
    USE   read_wind_farm_parameter_data,   ONLY: Xcoordinate, Ycoordinate, NumWT, Winddir, Tinfluencer,Mmt,ppR, Domain_X
    USE   wind_farm_geometry_data,         ONLY: turbine_sort, xwind, ywind,length
    
    INTEGER     ::     I,J,K
    INTEGER     ::     temp, Un
    REAL        ::     origin_turbine_x
    REAL        ::     origin_turbine_y
    REAL        ::     downwind_turbine_x
    REAL        ::     downwind_turbine_y
    REAL        ::     max_xturbine
    REAL        ::     max_yturbine
    REAL        ::     Circle_Radius
    REAL        ::     Pi
    REAL        ::     vector1x
    REAL        ::     vector1y
    REAL        ::     vector2x
    REAL        ::     vector2y
    REAL        ::     Max_spacing
    REAL        ::     spacing_temp
    
    INTEGER,ALLOCATABLE :: position_sign(:,:)
    
    REAL,ALLOCATABLE    ::     length_sort(:)
    REAL,ALLOCATABLE    ::     turbine_angle(:,:)

         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! determine the SIZE of the wind farm
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    max_xturbine = 0
    max_yturbine = 0
    
    
    DO I = 1,NumWT
        IF ( ABS(Xcoordinate(I)) > max_xturbine ) THEN
            max_xturbine = ABS(Xcoordinate(I))
        END IF
        
        IF ( ABS(Ycoordinate(I)) > max_yturbine ) THEN
            max_yturbine = ABS(Ycoordinate(I))
        END IF
    END DO
    
    Circle_Radius = 2*SQRT(max_xturbine*max_xturbine + max_yturbine*max_yturbine)
    
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! determine the inflow wind origin
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    Pi    = ACOS( -1.0 )
    xwind = SIN(Winddir*Pi/180) * circle_radius
    ywind = COS(Winddir*Pi/180) * circle_radius
    
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! sort the wind turbines as the order from upstream to downstream
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    ALLOCATE ( turbine_sort(NumWT) )
    ALLOCATE ( length(NumWT) )
    ALLOCATE ( length_sort(NumWT) )
         
    DO I = 1,NumWT
        length(I) = projected_length(xwind,ywind,Xcoordinate(I),Ycoordinate(I))
    END DO
    
    length_sort = length
       
    DO I=1,NumWT
        turbine_sort(I) = I
    END DO
    
    
    DO I = 1,NumWt-1
        DO J = 1,NumWt-1
            IF ( length_sort(J) > length_sort(J+1) ) THEN
               temp = turbine_sort(J)
               turbine_sort(J) = turbine_sort(J+1)
               turbine_sort(J+1) = temp
               
               temp = length_sort(J)
               length_sort(J) = length_sort(J+1)
               length_sort(J+1) = temp
            END IF
        END DO
    END DO
    
    !CALL GetNewUnit(Un)    
    Un = 10
    OPEN(unit = Un, status='replace',file='wind_farm_turbine_sort.bin',form='unformatted')    
    WRITE(Un)   turbine_sort(:)                                                                                                                                                                                                    
    CLOSE(Un)
    
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! calculate the angles between the line connecting two turbines and the wind direction
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    ALLOCATE (turbine_angle(NumWT,NumWT))
    ALLOCATE (position_sign(NumWT,NumWT))
    turbine_angle = 0.0
    position_sign = 0
         
    DO I = 1,NumWT
        DO J = 1,NumWT
            origin_turbine_x   = Xcoordinate(I)
            origin_turbine_y   = Ycoordinate(I)
            downwind_turbine_x = Xcoordinate(J)
            downwind_turbine_y = Ycoordinate(J)
            
            vector1x = 0 - xwind                                         ! the wind direction
            vector1y = 0 - ywind                       
            vector2x = downwind_turbine_x - origin_turbine_x              ! the vector between two turbines 
            vector2y = downwind_turbine_y - origin_turbine_y 
            
            IF (I /= J) THEN
               turbine_angle(I,J) = AngleBetweenVectors(vector1x,vector1y,vector2x,vector2y)
               position_sign(I,J) = turbine_position_LorR(xwind,ywind,origin_turbine_x,origin_turbine_y,downwind_turbine_x,downwind_turbine_y)
               turbine_angle(I,J) = turbine_angle(I,J)*position_sign(I,J) ! to specify the downstream turbine is on left side (+) or right side (-)
            END IF
        END DO
    END DO
    
    !CALL GetNewUnit(Un)    
    Un = 10
    OPEN(unit = Un, status='replace',file='turbine_angles.bin',form='unformatted')    
    WRITE(Un)   turbine_angle(:,:)                                                                                                                                                                                                    
    CLOSE(Un)
    
    OPEN(unit = Un, status='replace',file='turbine_distance.bin',form='unformatted')    
    WRITE(Un)   length(:)                                                                                                                                                                                                    
    CLOSE(Un)
    
    ! check if the Domain_X is larger than the maximum spacing
    Max_spacing = 0.0
    DO I = 1,NumWT
        IF (I+Tinfluencer < NumWT+1) THEN
            spacing_temp = length(I+Tinfluencer) - length(I)
            IF (spacing_temp > Max_spacing) THEN
                Max_spacing = spacing_temp
            END IF
        END IF
    END DO
    
    
    ! LOGIC WRONG !          2.14.2014
    !IF ( Domain_X/2 < Max_spacing ) THEN
    !    PRINT*, 'WARNING: Domain_x should be larger than the maximum turbine spacing'
    !    CALL EXIT
    !END IF
    
    ! test
    
    OPEN (unit=25,file="turbine_spacing.txt")
    WRITE (25,*) length(:)
    CLOSE(25)
    
END SUBROUTINE wind_farm_geometry

!-------------------------------------------------------------------
SUBROUTINE cal_wake_sector_angle()
!...................................................................
! This subroutine is to calculate the wake sector angle 
! with respect to the change of the downstream distance
!...................................................................
    USE wind_farm_geometry_data,        ONLY: wake_sector_angle_array,scale_factor,Pi,TurbineInfluenceData,length,xwind,ywind,turbine_sort
    USE read_wind_farm_parameter_data,  ONLY: NumWT, Winddir, Tinfluencer,Mmt,ppR,Xcoordinate, Ycoordinate,Mstl,Domain_X
    
    REAL,ALLOCATABLE  ::  distance_array(:)
    INTEGER           ::  I,J
    
    REAL              ::  DownStart 
    REAL              ::  DownEnd              
    REAL              ::  delta     
    INTEGER           ::  node
    
    INTEGER, ALLOCATABLE :: wake_width(:)
    REAL, ALLOCATABLE    :: wake_center_position(:,:,:)
    
    INTEGER     ::     Turbine0_index
    INTEGER     ::     Turbine1_index
    INTEGER     ::     Counter
    REAL        ::     Turbine_Spacing
    REAL        ::     vector1x
    REAL        ::     vector1y
    REAL        ::     vector2x
    REAL        ::     vector2y
    REAL        ::     Angle
    REAL        ::     Spacing_rounding
    REAL        ::     Sector_angle_boundary
    
    INTEGER     ::     T7_wake
    INTEGER     ::     T8_wake
    

         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! read the wake file and wake width
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''  
    ALLOCATE (wake_width          (NINT(ppR * Domain_X/2 )))
    !ALLOCATE (wake_width          (ppR * 36/2 ))
    ALLOCATE (wake_center_position(Mstl  ,Mmt+1,3   ))
    scale_factor = 10
    Pi = ACOS( -1.0 )
    
    ! read the wake file and wake width
    OPEN(unit = 10, status='old',file='Wake_width_Turbine_0.bin',form='unformatted')  ! open an existing file
    READ(10)  wake_width(:)
    CLOSE(10)
    
    OPEN(unit = 10, status='old',file='WC_Turbine_0.bin',form='unformatted')  ! open an existing file
    READ(10)  wake_center_position(:,:,:)
    CLOSE(10)
    
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! calculate the wake sector angle array
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''  
    DownStart = 1
    DownEnd   = 12             ! from 1D to 34D
    delta     = 0.1
    node      = (DownEnd-Downstart)/delta + 1
    
    ALLOCATE (distance_array(node))
    ALLOCATE (wake_sector_angle_array(node))
    distance_array = ((DownEnd-DownStart)/(node-1))*[(I,I=1,node)]+(DownStart-((DownEnd-DownStart)/(node-1)))
    
    DO I = 1,node
       wake_sector_angle_array(I) = Sector_angle(wake_width,wake_center_position,distance_array(I))
    END DO
    
    OPEN(unit = 10, status='replace',file='wake_sector_angle.bin',form='unformatted')    
    WRITE(10)   wake_sector_angle_array(:)                                                                                                                                                                                                    
    CLOSE(10)

         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! for a specfic downstream turbine, determine the wakes from which upstream turbines will make an effect 
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    !WindSectorAngle = 7
    ALLOCATE ( TurbineInfluenceData     (NumWT,Tinfluencer ) )
    
    TurbineInfluenceData = 0

    DO I = 1,NumWt
        Turbine0_index = turbine_sort(NumWt-I+1)               ! Start from the most downwind turbine
        Counter = 0
        DO J = 1,NumWt-1
           IF (Counter <Tinfluencer) THEN
                IF (NumWt-I+1-J>0) THEN
                    Turbine1_index = turbine_sort(NumWt-I+1-J)
                    
                    Turbine_Spacing = length(Turbine0_index) - length(Turbine1_index)     ! Absolute spacing
                    
                    Spacing_rounding = NINT(Turbine_Spacing*10.0)/10.0                    ! the rounding spacing x.x to calculate the sector angle
                    IF (Spacing_rounding>(Mmt/(ppR/scale_factor))) THEN
                        Spacing_rounding = Mmt/(ppR/scale_factor)
                    END IF
                    IF (Spacing_rounding<0.1) THEN
                        Spacing_rounding = 0.1
                    END IF 
                                   
                    vector1x = 0 - xwind                                                  ! the wind direction
                    vector1y = 0 - ywind                       
                    vector2x = Xcoordinate(Turbine0_index) - Xcoordinate(Turbine1_index)  ! the vector between two turbines 
                    vector2y = Ycoordinate(Turbine0_index) - Ycoordinate(Turbine1_index)
                
                    !IF (Turbine_Spacing >= 3) THEN
                        Angle = AngleBetweenVectors(vector1x,vector1y,vector2x,vector2y)
                        Sector_angle_boundary = Sector_angle(wake_width,wake_center_position,Spacing_rounding)
                        IF ( Angle < Sector_angle_boundary ) THEN
                            Counter = Counter + 1
                            TurbineInfluenceData(Turbine0_index,Counter) = Turbine1_index
                        END IF
                    !END IF  
                END IF
                           
           END IF           
        END DO
    END DO
    
    OPEN(unit = 10, status='replace',file='turbine_interaction.bin',form='unformatted')    
    WRITE(10)   TurbineInfluenceData(:,:)                                                                                                                                                                                                    
    CLOSE(10)

END SUBROUTINE cal_wake_sector_angle

!-------------------------------------------------------------------
FUNCTION projected_length(ax,ay,bx,by)
!...................................................................
! This function is to calculate the downstream distance from the specific
! turbine to the origin of the inflow wind
! The distance is used to sort the turbines from ipstream to downstream
!...................................................................

    REAL(8) ::   ax                   ! the x coordinate of the point a (reference point)
    REAL(8) ::   ay                   ! the y coordinate of the point a
    REAL(8) ::   bx                   ! the x coordinate of the point b
    REAL(8) ::   by                   ! the y coordinate of the point b
    REAL    ::   projected_length
    
    REAL  ::   side_a               ! triangle side a length
    REAL  ::   side_b               ! triangle side b length
    REAL  ::   side_c               ! triangle side c length
    REAL  ::   semi                 ! triangle semiperimeter
    REAL  ::   area                 ! triangle area
    REAL  ::   height               ! triangle height
    
    
         !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! calculate the area of the triangle whose corners are the (0,0), turbine location and the wind origin using Heron's method
         !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    side_a = SQRT( (bx-ax)*(bx-ax) + (by-ay)*(by-ay) )
    side_b = SQRT( bx*bx + by*by )
    side_c = SQRT( ax*ax + ay*ay )
    semi = (side_a + side_b + side_c)/2
    
    area = SQRT(ABS( semi * ( semi - side_a ) * ( semi - side_b ) * ( semi - side_c ) ))
    
         !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! calculate the projected length of the turbine on the wind direction vector
         !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    height = 2 * area / side_c
    
    projected_length = SQRT( side_a*side_a - height*height )
    
    
END FUNCTION projected_length

!-------------------------------------------------------------------
FUNCTION AngleBetweenVectors(x1,y1,x2,y2)
!...................................................................
! This function is to calculate the angle between two vectors
! x1,y1,x2,y2 are the vector coordinates
!...................................................................

    REAL       ::        x1
    REAL       ::        y1
    REAL       ::        x2
    REAL       ::        y2
    REAL       ::        Pi = ACOS( -1.0 )
    REAL       ::        AngleBetweenVectors
    REAL       ::        CosAlpha
    REAL       ::        Rad
    
    Pi = ACOS( -1.0 )
    
    CosAlpha = (x1*x2 + y1*y2) / ( SQRT(x1*x1+y1*y1) * SQRT(x2*x2+y2*y2) )
    
    !Rad = ACOS(CosAlpha-0.00001)                         ! ACOS(1) = NaN??!!!!
    
    IF (CosAlpha>1.0) THEN
        CosAlpha = 1.0
    END IF
    
    IF (CosAlpha<-1.0) THEN
        CosAlpha = -1.0
    END IF
    
    Rad = ACOS(CosAlpha)
    
    AngleBetweenVectors = Rad * 180/Pi
    
END FUNCTION AngleBetweenVectors
!-------------------------------------------------------------------
FUNCTION Sector_angle(wake_width,wake_center,downstream_distance)
!...................................................................
! This function is to calculate the maximum sector angle at a certain downstream
! distance behind a wind turbine
!................................................................... 
    USE read_wind_farm_parameter_data,     ONLY: ppR,Mstl,Mmt,RotorR
    USE wind_farm_geometry_data,           ONLY: Pi,scale_factor

    REAL                  ::       Sector_angle
    INTEGER               ::       wake_width(:)
    REAL                  ::       wake_center(:,:,:)
    REAL                  ::       downstream_distance
    
    INTEGER           :: I,J,K
    REAL              :: max_temp
    INTEGER           :: index_I
    REAL              :: Wakewidth
    REAL              :: angle_temp
    REAL              :: y1_temp,y2_temp,y_temp,x_temp
    
         !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! Calculate the max arctan(y/x) as the 1/2 sector angle
         !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    max_temp = 0
    
         ! Calculate the index of the corresponding downstream distance
    IF (MOD(downstream_distance,(1.0/(ppR/scale_factor))) == 0) THEN
       J = FLOOR(downstream_distance/(1.0/(ppR/scale_factor)))+1      ! resolution = 1/(ppR/scale_factor) D
    ELSE
       J = FLOOR(downstream_distance/(1.0/(ppR/scale_factor)))+2       
    END IF
    
    IF (J> Mmt+1) THEN
        J = Mmt+1
    END IF
    
    DO I = 1,Mstl
        !DO J = downstream_distance*ppR/scale_factor+1,Mmt+1
            Wakewidth = wake_width( (J-1)*scale_factor )
            
            y1_temp = ABS(wake_center(I,J,2))
            y2_temp = Wakewidth/ppR*RotorR
            y_temp  = y1_temp + y2_temp
            x_temp  = wake_center(I,J,1)
            
            angle_temp  = ATAN( y_temp / x_temp )
            
            IF (angle_temp > max_temp) THEN
                max_temp = angle_temp
                index_I = I
            END IF
        !END DO
    END DO
    
    Sector_angle = max_temp * 180/Pi

END FUNCTION Sector_angle

!-------------------------------------------------------------------
FUNCTION turbine_position_LorR(xwind,ywind,x0,y0,x1,y1)
!...................................................................
! This function is to see if the downwind turbine is on the left side or right side of the line
! connecting the wind origin and the upwind turbine
!...................................................................
    REAL(8) ::   xwind                ! wind origin
    REAL(8) ::   ywind
    REAL    ::   x0                   ! investigated origin turbine
    REAL    ::   y0
    REAL    ::   x1                   ! turbine other than the investigated turbine
    REAL    ::   y1
    INTEGER ::   turbine_position_LorR

    REAL    ::   flag                 
    
    flag = -xwind*(y1-y0) + ywind*(x1-x0)
    
    IF (flag >= 0) THEN
        turbine_position_LorR = 1         ! the downwind turbine is on the left side 
    ELSE
        turbine_position_LorR = -1        ! the downwind turbine is on the right side
    END IF
        
END FUNCTION turbine_position_LorR

!-------------------------------------------------------------------
SUBROUTINE rename_files(turbine_index)
!...................................................................
! This subroutine is to called to rename the output files of the DWM
! according to the turbine index
!...................................................................
    
    INTEGER           :: turbine_index
    CHARACTER(LEN=30) :: filename_u,filename_wakecenter,filename_meanU,filename_TI,filename_A,filename_FastOutput,filename_FastElm
    CHARACTER(LEN=30) :: filename_u_bin,filename_wakecenter_bin,filename_meanU_bin
    CHARACTER(LEN=3)  :: turbine_index_character
    CHARACTER(LEN=2)  :: Uprefix_txt        = 'U_'          ! wake_velocity
    CHARACTER(LEN=3)  :: WCprefix_txt       = 'WC_'         ! wake center
    CHARACTER(LEN=7)  :: MeanUprefix_txt    = 'Mean_U_'     ! average velocity
    CHARACTER(LEN=3)  :: TIprefix_txt       = 'TI_'         ! TI of each turbine downstream
    CHARACTER(LEN=6)  :: Inducprefix_txt    = 'Induc_'      ! Induction factor of each turbine
    CHARACTER(LEN=11) :: Fastprefix_out     = 'FastOutput_' ! Fast output file
    CHARACTER(LEN=8)  :: FastElmprefix_elm  = 'FastElm_'    ! Fast Elm output file
    !CHARACTER(LEN=11) :: Smoothprefix_bin   = 'Smoothwake_' ! Smoothed out wake for downwind turbines
    CHARACTER(LEN=2)  :: Uprefix_bin        = 'U_'
    CHARACTER(LEN=3)  :: WCprefix_bin       = 'WC_'
    CHARACTER(LEN=7)  :: MeanUprefix_bin    = 'Mean_U_'     
    
    
    CHARACTER(LEN=8)  :: Turbineprefix      = 'Turbine_'
    
    IF (turbine_index <= 9) THEN
        write(turbine_index_character,'(i1)') turbine_index
    ELSEIF (turbine_index <= 99) THEN
        write(turbine_index_character,'(i2)') turbine_index
    ELSE
        write(turbine_index_character,'(i3)') turbine_index
    END IF
    
    filename_u          = trim(Uprefix_txt)//trim(Turbineprefix)//trim(turbine_index_character)//".txt"
    filename_wakecenter = trim(WCprefix_txt)//trim(Turbineprefix)//trim(turbine_index_character)//".txt"
    filename_meanU      = trim(MeanUprefix_txt)//trim(Turbineprefix)//trim(turbine_index_character)//".txt"
    filename_TI         = trim(TIprefix_txt)//trim(Turbineprefix)//trim(turbine_index_character)//".txt"
    filename_A          = trim(Inducprefix_txt)//trim(Turbineprefix)//trim(turbine_index_character)//".txt"
    
    filename_FastOutput = trim(Fastprefix_out)//trim(Turbineprefix)//trim(turbine_index_character)//".out"
    filename_FastElm    = trim(FastElmprefix_elm)//trim(Turbineprefix)//trim(turbine_index_character)//".elm"
    
    !filename_u_bin          = trim(Uprefix_bin)//trim(Turbineprefix)//trim(turbine_index_character)//".bin"
    !filename_wakecenter_bin = trim(WCprefix_bin)//trim(Turbineprefix)//trim(turbine_index_character)//".bin"
    !filename_meanU_bin      = trim(MeanUprefix_bin)//trim(Turbineprefix)//trim(turbine_index_character)//".bin"
    

END SUBROUTINE rename_files

!-------------------------------------------------------------------
SUBROUTINE delete_temp_files()
!...................................................................
! This routine is called to delete the temporary DWM binary files
!...................................................................
    USE read_wind_farm_parameter_data, ONLY : NumWT 
    USE wind_farm_geometry_data,       ONLY : TurbineInfluenceData


    CHARACTER(LEN=3)  :: invetigated_turbine_index_character
    CHARACTER(LEN=3)  :: downwind_turbine_index_character
    CHARACTER(LEN=80) :: filename_u_bin,filename_wakecenter_bin,filename_meanU_bin,filename_TI_bin,filename_smoothWake_bin,filename_wakewidth_bin
    CHARACTER(LEN=80) :: filename_TI_txt,filename_meanU_txt,filename_induction_txt,filename_wake_txt,filename_wakecenter_txt
    CHARACTER(LEN=2)  :: Uprefix_bin        = 'U_'
    CHARACTER(LEN=3)  :: WCprefix_bin       = 'WC_'
    CHARACTER(LEN=7)  :: MeanUprefix_bin    = 'Mean_U_'
    CHARACTER(LEN=3)  :: Tiprefix_bin       = 'TI_'
    CHARACTER(LEN=11) :: SmoothWprefix_bin  = 'Smoothwake_'
    CHARACTER(LEN=10) :: InductionPrefix    = 'Induction_'
    CHARACTER(LEN=6)  :: Wakeprefix         = 'WakeU_'
    CHARACTER(LEN=11) :: WWprefix_bin       = 'Wake_width_'
    CHARACTER(LEN=22) :: Prefix             = 'DWM-results/'
    CHARACTER(LEN=4)  :: connectionprefix   = '_to_'
    CHARACTER(LEN=8)  :: Turbineprefix      = 'Turbine_'
    INTEGER           :: I,J,K
    INTEGER           :: downwindturbine_number
    INTEGER,ALLOCATABLE :: downwind_turbine_index(:)
    
          
    DO I = 1,NumWT
       !bjj: you could call Num2LStr(I) instead of using this if statement:
        IF (I <= 9) THEN
            write(invetigated_turbine_index_character,'(i1)') I
        ELSEIF (I <= 99) THEN
            write(invetigated_turbine_index_character,'(i2)') I
        ELSE
            write(invetigated_turbine_index_character,'(i3)') I
        END IF
    
        filename_meanU_bin = trim(Prefix)//trim(MeanUprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".bin"    
        OPEN (29, file=filename_meanU_bin)
        CLOSE (29, status='delete')
        
        filename_wakecenter_bin = trim(Prefix)//trim(WCprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".bin"
        OPEN (29, file=filename_wakecenter_bin)
        CLOSE (29, status='delete')
    END DO  
    
    
    
    DO K = 1,NumWT
        IF (K <= 9) THEN
            write(invetigated_turbine_index_character,'(i1)') K
        ELSEIF (K <= 99) THEN
            write(invetigated_turbine_index_character,'(i2)') K
        ELSE
            write(invetigated_turbine_index_character,'(i3)') K
        END IF
        
        
        downwindturbine_number = 0
        IF (.NOT. ALLOCATED(downwind_turbine_index) ) THEN 
            ALLOCATE (downwind_turbine_index(NumWT-1))
        END IF
        
        DO I = 1,1       !Tinfluencer
            DO J = 1,NumWT
                IF (TurbineInfluenceData(J,I) == K) THEN
                downwindturbine_number = downwindturbine_number + 1
                downwind_turbine_index(downwindturbine_number) = J
                END IF
            END DO
        END DO
        
        IF (downwindturbine_number /= 0 ) THEN
            DO I = 1,downwindturbine_number
                IF (downwind_turbine_index(I) <= 9) THEN
                   write(downwind_turbine_index_character,'(i1)') downwind_turbine_index(I)
                ELSEIF (downwind_turbine_index(I) <= 99) THEN
                   write(downwind_turbine_index_character,'(i2)') downwind_turbine_index(I)
                ELSE
                   write(downwind_turbine_index_character,'(i3)') downwind_turbine_index(I)
                END IF
            
                filename_u_bin          = trim(Prefix)//trim(Uprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                          //trim(connectionprefix)//trim(downwind_turbine_index_character)//".bin"
                filename_TI_bin         = trim(Prefix)//trim(TIprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                          //trim(connectionprefix)//trim(downwind_turbine_index_character)//".bin"
                filename_smoothWake_bin = trim(Prefix)//trim(SmoothWprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                          //trim(connectionprefix)//trim(downwind_turbine_index_character)//".bin"
    
                OPEN  (29, file=filename_u_bin)
                CLOSE (29, status='delete')
                
                OPEN  (29, file=filename_TI_bin)
                CLOSE (29, status='delete')
                
                OPEN  (29, file=filename_smoothWake_bin)
                CLOSE (29, status='delete')
            END DO
        END IF
    END DO
    
    
    OPEN (29, file='turbine_angles.bin')
    CLOSE (29, status='delete')

    OPEN (29, file='turbine_distance.bin')
    CLOSE (29, status='delete')
    
    OPEN (29, file='turbine_interaction.bin')
    CLOSE (29, status='delete')
    
    OPEN (29, file='wake_sector_angle.bin')
    CLOSE (29, status='delete')
    
    OPEN (29, file='Wake_width_Turbine_0.bin')
    CLOSE (29, status='delete')
    
    OPEN (29, file='WC_Turbine_0.bin')
    CLOSE (29, status='delete')
    
    OPEN (29, file='wind_farm_turbine_sort.bin')
    CLOSE (29, status='delete')
    
    
END SUBROUTINE delete_temp_files

!-------------------------------------------------------------------
SUBROUTINE Driver_init()
!...................................................................
! This routine called to initiate the Driver program
!...................................................................
    USE DWM_init_data, ONLY:OutFileRoot, InputFile
    
    INTEGER                      :: Stat
    CHARACTER(1024)              :: DirName

    CALL CheckArgs( InputFile, Stat ) 
    
    CALL GetRoot( InputFile, OutFileRoot )
    
    CALL Get_CWD  ( DirName, Stat )

END SUBROUTINE Driver_init

!----------------------------------------------------------------------------------
SUBROUTINE rename_FAST_output(SimulationOrder_index)
!............................................................................
! This routine is called to rename the fast output
!............................................................................
#ifdef __INTEL_COMPILER
    USE IFPORT
#endif
    USE wind_farm_geometry_data,         ONLY: turbine_sort
    USE DWM_init_data,                   ONLY: OutFileRoot
    
    CHARACTER(LEN=80) :: filename_FastOutput,filename_FastElm
    CHARACTER(LEN=11) :: Fastprefix         = 'FastOutput_' ! Fast output file
    CHARACTER(LEN=8)  :: FastElmprefix      = 'FastElm_'    ! Fast Elm output file
    INTEGER           :: RESULT
    CHARACTER(LEN=22) :: Prefix             = 'DWM-results/'
    CHARACTER(LEN=8)  :: Turbineprefix      = 'Turbine_'
    CHARACTER(LEN=3)  :: invetigated_turbine_index_character
    INTEGER           :: WT_index
    INTEGER           :: SimulationOrder_index
    
    IF ( SimulationOrder_index > 0 ) THEN            ! exclude the first turbine
        
        WT_index = turbine_sort(SimulationOrder_index)
        
        IF (WT_index <= 9) THEN
            write(invetigated_turbine_index_character,'(i1)') WT_index
        ELSEIF (WT_index <= 99) THEN
            write(invetigated_turbine_index_character,'(i2)') WT_index
        ELSE
            write(invetigated_turbine_index_character,'(i3)') WT_index
        END IF

            ! Rename the FAST output wrt the turbine index
        filename_FastOutput = trim(Prefix)//trim(Fastprefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".out"
        filename_FastElm    = trim(Prefix)//trim(FastElmprefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".AD.out"

        RESULT =rename( TRIM(OutFileRoot)//'.out',   filename_FastOutput)          !bjj: what if I'm using .outb in FAST instead of .out?
        RESULT =rename( TRIM(OutFileRoot)//'.AD.out',filename_FastElm   )          !bjj: *.elm has been renamed *.AD.out in FAST v8. Also, .AD.out files are not always generated.
    END IF

END SUBROUTINE rename_FAST_output

END MODULE DWM_driver_wind_farm_sub