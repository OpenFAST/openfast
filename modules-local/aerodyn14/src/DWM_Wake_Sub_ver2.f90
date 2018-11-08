MODULE DWM_Wake_Sub

   USE DWM_Types
   USE NWTC_Library
   !USE InflowWind
   
   IMPLICIT NONE
   
         ! ..... Public Subroutines ............

   PUBLIC :: turbine_average_velocity
   PUBLIC :: pass_velocity
   PUBLIC :: filter_average_induction_factor
   PUBLIC :: calculate_mean_u
   PUBLIC :: calculate_element_area
   PUBLIC :: calculate_induction_factor
   PUBLIC :: get_initial_condition
   PUBLIC :: calculate_wake
   PUBLIC :: create_F1_filter
   PUBLIC :: create_F2_filter
   PUBLIC :: Gauss
   PUBLIC :: shear_correction
   PUBLIC :: filter_velocity
   PUBLIC :: Get_wake_center
   PUBLIC :: smooth_out_wake
   PUBLIC :: smooth_wake_shifted_velocity
   PUBLIC :: shifted_velocity
   PUBLIC :: TI_downstream_total
   PUBLIC :: smallscale_TI
   PUBLIC :: read_parameter_file
   PUBLIC :: read_turbine_position
   PUBLIC :: read_upwind_result_file
   PUBLIC :: write_result_file
!   PUBLIC :: rename_FAST_output
   PUBLIC :: min_of_array
   PUBLIC :: max_of_array

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE CalVelScale(u,v,y,z)
!..................................................................................................................................
! This routine is to calculat the atmospheric length scale before introducing the TI term (which will be used later)
!.................................................................................................................................. 
 !  IMPLICIT NONE
   
  ! TYPE(DWM_OutputType),     INTENT(INOUT)   :: y
   !TYPE(DWM_ConstraintStateType),     INTENT(INOUT)   :: z

      !! Internal variables
   !REAL(ReKi)            :: u   ! atmospheric U velocity
   !REAL(ReKi)            :: v   ! atmospheric V velocity
   
   !z%CalVelScale_data%counter     = z%CalVelScale_data%counter + 1
   
   !z%CalVelScale_data%Denominator = (z%CalVelScale_data%Denominator * (z%CalVelScale_data%counter-1) + u*v)/z%CalVelScale_data%counter
   !z%CalVelScale_data%Numerator   = (z%CalVelScale_data%Numerator   * (z%CalVelScale_data%counter-1) + u*v)/z%CalVelScale_data%counter
   
   !y%AtmUscale                    = z%CalVelScale_data%Numerator / z%CalVelScale_data%Denominator
   
!END SUBROUTINE CalVelScale 


!!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE turbine_average_velocity( single_velocity, blade_num, element, y,X,z)
!..................................................................................................................................
! This routine is called at every time step of the Aerodyn simuilation.
! To calculate the average of the wind speed of a specific blade ring, the outpout is the average_velocity_array
!..................................................................................................................................
 !  IMPLICIT NONE
   
  ! TYPE(DWM_ConstraintStateType),  INTENT(INOUT)   :: z
   !TYPE(DWM_OutputType),           INTENT(INOUT)   :: y
   
    !  ! Internal variables
   !REAL(ReKi)                ::   single_velocity
   !REAL(ReKi),ALLOCATABLE    ::   y%Mean_FFWS_array(:)
   !INTEGER(IntKi)            ::   element
   !INTEGER(IntKi)            ::   blade_num
   !INTEGER(IntKi)            ::   I
   
   !z%turbine_average_velocity_data%time_step_velocity = z%turbine_average_velocity_data%time_step_velocity + 1
   
   !IF (z%turbine_average_velocity_data%time_step_velocity == 0) THEN
    !  ALLOCATE (y%Mean_FFWS_array                            (X%ElOut%NumElOut))
     ! ALLOCATE (z%turbine_average_velocity_data%time_step_velocity_array(X%ElOut%NumElOut))
      !y%Mean_FFWS_array                                      = 0
      !z%turbine_average_velocity_data%time_step_velocity_array          = 0
      !y%Mean_FFWS_array(element)                             = single_velocity
      !z%turbine_average_velocity_data%time_step_velocity_array(element) = 1
   
   !ELSE IF (z%turbine_average_velocity_data%time_step_velocity > 0) THEN
    !  DO I = 1,X%ElOut%NumElOut
     !    IF ( element == I ) THEN
      !      z%turbine_average_velocity_data%time_step_velocity_array(element) = z%turbine_average_velocity_data%time_step_velocity_array(element) + 1
       !     y%Mean_FFWS_array(element) = ( y%Mean_FFWS_array(element)*( (z%turbine_average_velocity_data%time_step_velocity_array(element)-1) )&
        !                                            +single_velocity)/z%turbine_average_velocity_data%time_step_velocity_array(element)
         !   IF ( I == X%ElOut%NumElOut .AND. blade_num = X%Blade%NB ) THEN
          !     CALL pass_velocity(y%Mean_FFWS_array,Q,X)
           !    z%turbine_average_velocity_data%time_step_velocity = -1
            !   IF ( ALLOCATED( y%Mean_FFWS_array ))                              DEALLOCATE ( y%Mean_FFWS_array )
             !  IF ( ALLOCATED( z%turbine_average_velocity_data%time_step_velocity_array ))  DEALLOCATE ( z%turbine_average_velocity_data%time_step_velocity_array )
            !END IF
         !END IF
      !END DO
   !END IF
   
!END SUBROUTINE turbine_average_velocity

!----------------------------------------------------------------------------------
SUBROUTINE turbine_average_velocity( p, m, single_velocity, blade_num, element, average_velocity_array_local )
!..................................................................................
! This routine is called at every time step of the Aerodyn simuilation.
! To calculate the average of the wind speed of a specific blade ring
! the outpout is the average_velocity_array
!..................................................................................
    !USE TAVD,           ONLY: m%TAVD%time_step_velocity_array, m%TAVD%time_step_velocity
    
    TYPE(DWM_MiscVarType),         INTENT(INOUT)  :: m
    TYPE(DWM_ParameterType),       INTENT(IN   )  :: p
    
    REAL(ReKi),             INTENT(IN)    ::   single_velocity
    REAL(ReKi),ALLOCATABLE, INTENT(INOUT) ::   average_velocity_array_local(:)
    INTEGER, intent(in)             ::   element
    INTEGER, intent(in)             ::   blade_num
    INTEGER             ::   i
    !INTEGER             ::   m%TAVD%time_step_velocity = -1
    !INTEGER,ALLOCATABLE,SAVE   :: m%TAVD%time_step_velocity_array(:) ! counter of each section of the blade
    
    m%TAVD%time_step_velocity = m%TAVD%time_step_velocity +1
    
    IF ( m%TAVD%time_step_velocity==0) THEN
       ALLOCATE ( average_velocity_array_local(p%ElementNum) )
       ALLOCATE ( m%TAVD%time_step_velocity_array(p%ElementNum) )
       average_velocity_array_local(:)       = 0
       m%TAVD%time_step_velocity_array(:)           = 0
       average_velocity_array_local(element) = single_velocity
       m%TAVD%time_step_velocity_array(element)     = 1
       
    ELSE IF (m%TAVD%time_step_velocity > 0) THEN
       DO i=1,p%ElementNum
          IF ( element == i) THEN
             m%TAVD%time_step_velocity_array(element)     = m%TAVD%time_step_velocity_array(element)+1
             average_velocity_array_local(element) = (average_velocity_array_local(element)*( (m%TAVD%time_step_velocity_array(element)-1) ) + single_velocity) &
                                               / (m%TAVD%time_step_velocity_array(element))
             IF ( element == p%ElementNum ) THEN
                IF ( blade_num == p%Bnum) THEN
                   CALL pass_velocity(p, m, average_velocity_array_local)
                   m%TAVD%time_step_velocity = -1
                   IF (ALLOCATED( average_velocity_array_local ))              DEALLOCATE ( average_velocity_array_local )
                   IF (ALLOCATED( m%TAVD%time_step_velocity_array ))          DEALLOCATE ( m%TAVD%time_step_velocity_array )
                END IF            
             END IF                
          END IF
       END DO
    END IF
    
    
END SUBROUTINE turbine_average_velocity

!!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE pass_velocity(array_velocity,z,X)
!..................................................................................................................................
!..................................................................................................................................    
  ! IMPLICIT NONE
   !TYPE(DWM_ConstraintStateType),       INTENT(INOUT)   :: z
     
    !  ! Internal variables 
   !REAL(ReKi),ALLOCATEBLE           ::  array_velocity(:)

    
   !z%turbine_average_velocity_data%time_step_pass_velocity = z%turbine_average_velocity_data%time_step_pass_velocity+1
    
   !IF( z%turbine_average_velocity_data%time_step_pass_velocity==0 ) THEN
    !  ALLOCATE (z%turbine_average_velocity_data%average_velocity_array_temp(X%ElOut%NumElOut))
     ! z%turbine_average_velocity_data%average_velocity_array_temp(:) = array_velocity(:)
   !ELSE IF( z%turbine_average_velocity_data%time_step_pass_velocity>0 ) THEN
    !  z%turbine_average_velocity_data%average_velocity_array_temp(:) = array_velocity(:)
   !END IF
   
!END SUBROUTINE pass_velocity

!----------------------------------------------------------------------------------
SUBROUTINE pass_velocity(p, m, array_velocity)
!..................................................................................
! 
! 
!  
!..................................................................................
    !USE TAVD,            ONLY: m%TAVD%average_velocity_array_temp, m%TAVD%time_step_pass_velocity
    
    TYPE(DWM_MiscVarType),         INTENT(INOUT)  :: m
    TYPE(DWM_ParameterType),       INTENT(IN   )  :: p           ! Parameters
    
    
    REAL(ReKi)                  ::  array_velocity(p%ElementNum)
    !INTEGER,SAVE          ::  m%TAVD%time_step_pass_velocity = -1
    
    m%TAVD%time_step_pass_velocity = m%TAVD%time_step_pass_velocity+1
    
    IF (m%TAVD%time_step_pass_velocity==0) THEN
       ALLOCATE (m%TAVD%average_velocity_array_temp(p%ElementNum))
       m%TAVD%average_velocity_array_temp(:) = array_velocity(:)
    ELSE IF(m%TAVD%time_step_pass_velocity>0) THEN
       m%TAVD%average_velocity_array_temp(:) = array_velocity(:)
    END IF
        
END SUBROUTINE pass_velocity
!!----------------------------------------------------------------------------------------------------------------------------------
!SUBROUTINE filter_average_induction_factor( X, z, y )
!..................................................................................................................................
! This routine is called at every time step of the Aerodyn simuilation.
! The output of the subroutine is induction factor at this time step 
! and the average induction factor through all the time steps which have been simulated
!.................................................................................................................................. 
  ! IMPLICIT NONE
  ! TYPE(DWM_OutputType),       INTENT(INOUT)   :: y
  ! TYPE(DWM_ConstraintStateType),       INTENT(INOUT)   :: z
    
     ! ! Internal variables
   !REAL(ReKi)                     ::  thrust_coefficient          ( X%Element%NELM )
   !REAL(ReKi)                     ::  average_induction_factor    ( X%Element%NELM )
   !REAL(ReKi)                     ::  induction_factor_local_temp ( X%Element%NELM )
    
    
   !z%turbine_average_velocity_data%time_step_force = z%turbine_average_velocity_data%time_step_force +1
    
    
   !IF ( z%turbine_average_velocity_data%time_step_force==0) THEN
    !  ALLOCATE (y%induction_factor ( X%Element%NELM ))
     ! ALLOCATE (z%turbine_average_velocity_data%average_velocity_array ( X%Element%NELM ))
      !ALLOCATE (y%turbine_thrust_force (X%Element%NELM ))
      !ALLOCATE (z%turbine_average_velocity_data%swept_area (X%Element%NELM ))
      !y%turbine_thrust_force (:) = X%Blade%NB * X%Element%DFNSAV(:)
      !CALL calculate_element_area ( X%Blade%R, X%Element%NElm, X%Element%RELM(:), z%turbine_average_velocity_data%swept_area )
      !CALL calculate_induction_factor ( y%turbine_thrust_force , z%turbine_average_velocity_data%swept_area , X%Element%NELM, &
       !                                 z%turbine_average_velocity_data%average_velocity_array_temp, induction_factor_local_temp )
      !y%induction_factor = induction_factor_local_temp
      !z%turbine_average_velocity_data%average_velocity_array = z%turbine_average_velocity_data%average_velocity_array_temp
   !ELSE IF ( z%turbine_average_velocity_data%time_step_force>0) THEN
    !  y%turbine_thrust_force (:) = X%Blade%NB * X%Element%DFNSAV(:)
     ! CALL calculate_induction_factor ( y%turbine_thrust_force , z%turbine_average_velocity_data%swept_area , X%Element%NELM, &
      !                                  z%turbine_average_velocity_data%average_velocity_array_temp, induction_factor_local_temp )
      !y%induction_factor = ( y%induction_factor(:) * z%turbine_average_velocity_data%time_step_force + induction_factor_local_temp(:) ) &
       !                     / ( z%turbine_average_velocity_data%time_step_force+1 )
      !z%turbine_average_velocity_data%average_velocity_array = ( z%turbine_average_velocity_data%average_velocity_array(:) * &
       !                   z%turbine_average_velocity_data%time_step_force + z%turbine_average_velocity_data%average_velocity_array_temp(:) ) &
        !                    / ( z%turbine_average_velocity_data%time_step_force+1 )
   !END IF

!END SUBROUTINE filter_average_induction_factor

!----------------------------------------------------------------------------------
SUBROUTINE filter_average_induction_factor( m, p, y, thrust_force, num_of_element, dr_blade)
!..................................................................................
! This routine is called at every time step of the Aerodyn simuilation.
! The output of the subroutine is induction factor at this time step 
! and the average induction factor through all the time steps which have been simulated
!.................................................................................. 
    !USE TAVD,           ONLY: m%TAVD%time_step_force, m%TAVD%swept_area, m%TAVD%average_velocity_array, m%TAVD%average_velocity_array_temp
    !USE DWN_OutputType, ONLY: m%induction_factor, m%turbine_thrust_force
    !USE Blade,          ONLY: R
    
    TYPE(DWM_MiscVarType),         INTENT(INOUT)  :: m
    TYPE(DWM_ParameterType),       INTENT(IN   )  :: p           ! Parameters
    TYPE(DWM_OutputType),          INTENT(INOUT)  :: y
    
    
    INTEGER, INTENT(IN)          ::  num_of_element                   ! The number of the nodes in the blade
    REAL(ReKi), INTENT(IN)       ::  thrust_force ( num_of_element,p%BNum )  ! Thrust force at each node
    !INTEGER, SAVE                ::  m%TAVD%time_step_force = -1                   ! The time step (save attribute) of the FAST simulation
    REAL(ReKi)                   ::  thrust_coefficient ( num_of_element )
    REAL(ReKi)                   ::  average_induction_factor   ( num_of_element )
    REAL(ReKi)                   ::  induction_factor_local_temp ( num_of_element )
    INTEGER                      ::  I,J
    REAL(ReKi), INTENT(IN)       ::  dr_blade ( num_of_element )
    
    
    m%TAVD%time_step_force = m%TAVD%time_step_force +1
    
    
    IF ( m%TAVD%time_step_force==0) THEN
       ALLOCATE (y%induction_factor ( num_of_element ))
       ALLOCATE (m%TAVD%average_velocity_array ( num_of_element ))
       ALLOCATE (y%turbine_thrust_force (num_of_element ))
       ALLOCATE (m%TAVD%swept_area (num_of_element ))
       
       y%turbine_thrust_force = 0
       DO I = 1,num_of_element
           DO J = 1,p%BNum
           y%turbine_thrust_force (I) = y%turbine_thrust_force (I) + thrust_force(I,J)
           END DO
       END DO
       
       DO I = 1,num_of_element
           y%turbine_thrust_force (I) = y%turbine_thrust_force(I)   !* dr_blade(I)           ! integrate dFn through blade
       END DO
       
       CALL calculate_element_area ( p%RotorR, p%ElementNum, p%ElementRad(:), m%TAVD%swept_area )
       CALL calculate_induction_factor ( p, y%turbine_thrust_force , m%TAVD%swept_area , num_of_element, m%TAVD%average_velocity_array_temp, induction_factor_local_temp )
       y%induction_factor = induction_factor_local_temp
       m%TAVD%average_velocity_array = m%TAVD%average_velocity_array_temp
    
    ELSE IF ( m%TAVD%time_step_force>0) THEN

       y%turbine_thrust_force = 0
       DO J = 1,p%BNum
           DO I = 1,num_of_element 
           y%turbine_thrust_force (I) = y%turbine_thrust_force (I) + thrust_force(I,J)
           END DO
       END DO
       
       DO I = 1,num_of_element
           y%turbine_thrust_force (I) = y%turbine_thrust_force(I)    !* dr_blade(I)           ! integrate dFn through blade
       END DO
       
       CALL calculate_induction_factor ( p, y%turbine_thrust_force , m%TAVD%swept_area , num_of_element, m%TAVD%average_velocity_array_temp, induction_factor_local_temp )
       y%induction_factor = ( y%induction_factor(:) * m%TAVD%time_step_force + induction_factor_local_temp(:) ) / ( m%TAVD%time_step_force+1 )
       m%TAVD%average_velocity_array = ( m%TAVD%average_velocity_array(:) * m%TAVD%time_step_force + m%TAVD%average_velocity_array_temp(:) ) / ( m%TAVD%time_step_force+1 )
    END IF
    
    !print*, y%induction_factor(40), induction_factor_local_temp(40)
    
END SUBROUTINE filter_average_induction_factor

!----------------------------------------------------------------------------------
SUBROUTINE calculate_mean_u( m, p, u, num_element,r_t,turbine_mean_velocity,TI_normalization, FAST_Time )
!..................................................................................
! This routine is called to calculate the mean velocity and the TI of the turbine
! Using weighting method according to the blade ring area 
!.................................................................................. 
    !USE TAVD,                          ONLY : m%TAVD%average_velocity_array
    !USE read_turbine_position_data   , ONLY : m%RTPD%SimulationOrder_index,m%RTPD%upwindturbine_number
    !USE weighting_method             , ONLY : m%weighting_method%sweptarea,m%weighting_method%weighting_denominator
    !USE read_upwind_result_file_data , ONLY : m%Upwind_result%upwind_TI
    !USE Blade                        , ONLY : R
    
    TYPE(DWM_MiscVarType),         INTENT(INOUT)  :: m
    TYPE(DWM_ParameterType),       INTENT(IN   )  :: p           ! Parameters
    TYPE(DWM_InputType),           INTENT(INOUT)  :: u
    
    INTEGER           ::    num_element
    INTEGER           ::    i
    REAL(ReKi)              ::    r_t ( num_element )              ! The distance from the node to the hub
    REAL(ReKi)              ::    turbine_mean_velocity            ! turbine mean velocity
    REAL(ReKi)              ::    node_radius    ( num_element )
    REAL(ReKi)              ::    element_length ( num_element )
    REAL(ReKi)              ::    TI_normalization
    REAL(ReKi)              ::    FAST_Time
    
    ! check if the meandering simulation time is valid
    IF (p%WakePosition_1 < FAST_Time/ ( (20*p%RotorR/p%p_p_r)/(0.32*p%Uambient) ) + 1 ) THEN
         ! bjj: this at least is standard fortran, but calling ProgAbort is not allowed in a module in the FAST framework. Please trap your errors and return an error code.
        CALL ProgAbort('WARNING: Meandering_simulation_time is not valid, please refer to the DWM manual') 
    END IF

    ALLOCATE (m%weighting_method%sweptarea(num_element))
    
    m%weighting_method%weighting_denominator = 0
    turbine_mean_velocity = 0
    
    element_length (num_element) = 2.0*( p%RotorR - r_t(num_element) )
    DO i=num_element-1,1,(-1)
       element_length(i)= 2.0*( r_t(i+1)-r_t(i) ) - element_length (i+1)
    END DO

    node_radius ( num_element ) = p%RotorR - element_length (num_element)
    DO i=num_element-1,1,(-1)
       node_radius (i) = node_radius (i+1) - element_length (i)
    END DO

    DO i=1, num_element-1,1
       m%weighting_method%sweptarea (i) = Pi * (node_radius (i+1) **2 - node_radius (i) **2)
    END DO
    m%weighting_method%sweptarea (num_element) = Pi * (p%RotorR**2-node_radius (num_element)**2) ! ring area
    
    DO i=1,num_element
       m%weighting_method%weighting_denominator = m%weighting_method%weighting_denominator + m%weighting_method%sweptarea (i)   ! denominator
    END DO
    
    ! calculate the mean velocity of the turbine using weighting method
    DO i=1,num_element 
       turbine_mean_velocity = turbine_mean_velocity + m%weighting_method%sweptarea (i) / m%weighting_method%weighting_denominator * m%TAVD%average_velocity_array(i)       
    END DO
    
    IF (p%RTPD%SimulationOrder_index == 1 .OR. p%RTPD%SimulationOrder_index == 0) THEN
       TI_normalization = p%TI_amb
    ELSE
        IF(p%RTPD%upwindturbine_number /= 0) THEN           ! superimpose the TI from upstream wakes
            TI_normalization  =  0
            DO I = 1,p%RTPD%upwindturbine_number               
               TI_normalization = (TI_normalization**2 + u%Upwind_result%upwind_TI(I)**2)**0.5
            END DO
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           TI_normalization = u%Upwind_result%upwind_TI(1)          ! only take the TI effect from the closest upstream turbine
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
           TI_normalization = TI_normalization/(turbine_mean_velocity/p%Uambient)
        ELSE
           TI_normalization = p%TI_amb
        END IF
    END IF 
     
END SUBROUTINE calculate_mean_u

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE calculate_element_area (blade_radius, num_element, r_t, sweptarea)
!..................................................................................................................................
! This routine is called when the Aerodyn simuilation finishes.
! This routine is called to calculate the swept area of each blade section.
! The output of the subroutine is swept_area (:), which is the swept area of of each blade element.
!..................................................................................................................................
   IMPLICIT NONE
    
      ! Internal variables
   INTEGER(IntKi)        :: num_element
   INTEGER(IntKi)        :: I
   REAL(ReKi)            :: blade_radius
   REAL(ReKi)            :: r_t ( num_element )              ! The distance from the node to the hub
   REAL(ReKi)            :: node_radius    ( num_element )
   REAL(ReKi)            :: element_length ( num_element )
   REAL(ReKi)            :: sweptarea(num_element)

   element_length (num_element) = 2.0*( blade_radius - r_t(num_element) )
   DO I=num_element-1,1,(-1)
      element_length(I)= 2.0*( r_t(I+1)-r_t(I) ) - element_length (I+1)
   END DO

   node_radius ( num_element ) = blade_radius - element_length (num_element)
   DO I=num_element-1,1,(-1)
      node_radius (I) = node_radius (I+1) - element_length (I)
   END DO

   DO I=1, num_element-1,1
      sweptarea (I) = Pi * (node_radius (I+1) **2 - node_radius (I) **2)
   END DO
   sweptarea (num_element) = Pi * (blade_radius**2-node_radius (num_element)**2)
    
END SUBROUTINE calculate_element_area

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE calculate_induction_factor ( p, normal_force , element_swept_area , num_element, FFWS_array, induction_factor_local )
!..................................................................................................................................
! This routine is called to calculate thrust coefficient then the induction factor using local thrust force.
! The output of the subroutine is induction_factor (:), which is the induction factor of each blade node.
!..................................................................................................................................
   
   TYPE(DWM_ParameterType),       INTENT(IN   )         :: p
   
    
      ! Internal variables   
   INTEGER(IntKi)  ::   num_element
   INTEGER(IntKi)  ::   I
   REAL(ReKi)      ::   normal_force ( num_element )
   REAL(ReKi)      ::   element_swept_area ( num_element )
   REAL(ReKi)      ::   thrust_coefficient ( num_element )
   REAL(ReKi)      ::   induction_factor_local   ( num_element )
   REAL(ReKi)      ::   FFWS_array( num_element)
   REAL(ReKi)      ::   Ct_1
   REAL(ReKi)      ::   a_t
   REAL(ReKi)      ::   Ct_critical


   !Calculate thrust coefficient
   DO I=1,num_element
     thrust_coefficient (I) = normal_force(I)/(0.5* p%air_density * element_swept_area(I)* FFWS_array(I)**2)
   END DO
    
    
    ! Then calculate the induction factor by solving 4a(1-a)= Ct
    ! Applying the Glauert empirical Ct modification (10.7.2013)
      
    Ct_1       = 1.816
    a_t        = 1 - 0.5*SQRT(Ct_1)
    Ct_critical = 4*a_t*(1-a_t)
    
    DO I=1,num_element
         IF (thrust_coefficient(I)<=Ct_critical) THEN
            induction_factor_local (I) = (-4 + (16-16*thrust_coefficient(I))**(0.5))/(2*(-4))
         ELSE
            induction_factor_local (I) = 1 - (thrust_coefficient(I)-Ct_1) / (-4*(SQRT(Ct_1)-1))
         END IF
    END DO

END SUBROUTINE calculate_induction_factor

!----------------------------------------------------------------------------------
SUBROUTINE get_initial_condition( m, p, u, y, induc_array, r_t, element_num, r_w, U_w )
!..................................................................................
! This routine is called at the end of the subroutine calculate_initial_condition.
! This routine is called to calculate the initial condition of the DWM model.
! The output of the subroutine is r_wake (:) and U_wake (:).
! Which are the the scaled rotor radius and the scaled velocity at the rotor.
!..................................................................................
    !USE     read_turbine_position_data,   ONLY: m%RTPD%SimulationOrder_index,m%RTPD%upwindturbine_number
    !USE     DWM_ParameterType,            ONLY: p%smoothed_wake,p%smooth_flag,p%Uambient
    !USE     DWM_OutputType,               ONLY: m%Mean_FFWS 
    !USE     read_upwind_result_file_data, ONLY: m%Upwind_result%upwind_smoothWake
    !USE     BLADE,                        ONLY: R
    
    TYPE(DWM_MiscVarType),         INTENT(INOUT)  :: m
    TYPE(DWM_ParameterType),       INTENT(IN   )  :: p           ! Parameters
    TYPE(DWM_InputType),           INTENT(INOUT)  :: u
    TYPE(DWM_OutputType),          INTENT(INOUT)  :: y
    
    
    INTEGER           :: element_num
    INTEGER           :: i,J
    REAL(ReKi)              :: induc_array(element_num)
    REAL(ReKi)              :: r_t(element_num)
    REAL(ReKi)              :: dA (element_num-1)
    REAL(ReKi)              :: a_cellC (element_num-1)
    REAL(ReKi)              :: mean_a
    REAL(ReKi)              :: f_w
    REAL(ReKi)              :: fU !fU factor (realised induction for wake depth) {0-1}
    REAL(ReKi)              :: fR !fR factor (realised expansion for wake width) {0-1} 
    REAL(ReKi), ALLOCATABLE :: r_w(:)
    REAL(ReKi), ALLOCATABLE :: U_w(:)
    

    ALLOCATE       (r_w(element_num))
    ALLOCATE       (U_w(element_num))
    fU             = 1.10           !1.10 is not working when induction factor is reaching 0.5: (1-a*(fu+1))<1     !!!! 0.92
    fR             = 0.98
    

   !-------------------------------------------------------------------------------------------------
   ! apply the smoothed wake profile as the input wind profile for downstream turbine
   !-------------------------------------------------------------------------------------------------

   !-------------------------------------------------------------------------------------------------
   ! calculate the boundary condition
   !-------------------------------------------------------------------------------------------------
    !===Initial Condition===
    DO i=1,element_num-1
       dA (i) =( (r_t(i+1)/p%RotorR)**2 - (r_t(i)/p%RotorR)**2 )*Pi
    END DO
    !== simple approximation of cell center value
    DO i=1,element_num-1
       a_cellC (i) = ( induc_array(i+1) +induc_array(i) ) /2
    END DO

    !== Boundary Condition r_w
    mean_a = 0.0
    DO i=1,element_num-1
      mean_a = ( a_cellC (i) * dA (i) ) /Pi + mean_a
    END DO

    !== Uniform expansion
    f_w = ( (1-mean_a) / (1- ((1+fR) * mean_a)) )**0.5
    r_w = r_t /p%RotorR * f_w

    !== Boundary velocity
    !U_w = 1 - (induc_array * (1+fU))
    
    ! superimpose the smoothed wake from upwind turbines
    IF (p%RTPD%SimulationOrder_index == 1 .OR. p%RTPD%SimulationOrder_index == 0) THEN
        U_w = 1 - (induc_array * (1+fU))
    ELSEIF(p%RTPD%SimulationOrder_index > 1) THEN
        IF (p%RTPD%upwindturbine_number == 0) THEN
            U_w = 1 - (induc_array * (1+fU))
        ELSEIF (p%RTPD%upwindturbine_number > 0) THEN
            !!ALLOCATE (p%smoothed_wake(element_num))
            !!p%smoothed_wake = 1
            !!DO I = 1,p%RTPD%upwindturbine_number
                !!DO J = 1,element_num
                    !!p%smoothed_wake(J) = 1- ( (1-p%smoothed_wake(J))**2 + (1-u%Upwind_result%upwind_smoothWake(I,J))**2 )**0.5
                !!END DO
            !!END DO
            
            DO I = 1,element_num
                !U_w(I) = (p%Uambient/y%Mean_FFWS) * u%Upwind_result%upwind_smoothWake(1,I)*(1 - (induc_array(I) * (1+fU)))
                U_w(I) = (y%Mean_FFWS/p%Uambient) *(1 - (induc_array(I) * (1+fU)))
            END DO
        END IF
    END IF
    
    DO i=1,element_num,1           ! modification for low wind speed, high thrust situation
        IF (U_w(i) < 0.0) THEN
            U_w(i) = 0.01
        END IF
    END DO
    
    ! calculate the average induction factor of the rotor plane
    !avg_induction_factor = 0
    
    !DO i=1,element_num
       !avg_induction_factor = avg_induction_factor + m%weighting_method%sweptarea (i) / m%weighting_method%weighting_denominator * induc_array(i)       
    !END DO
    
    !m%skew_angle = 0.60*avg_induction_factor*NacYaw *(-1)    ! minus sign means different direction
    
    
    !--- calculate the average thrust coefficient and the ct_tilde ---
    
    y%avg_ct = 0
    
    !-------test--------
    !u%NacYaw = 0.00
    !-------------------
    
    DO i=1,element_num
       y%avg_ct = y%avg_ct + m%weighting_method%sweptarea (i) / m%weighting_method%weighting_denominator * ( 4*induc_array(i)*(1-induc_array(i)) )       
    END DO
    
    !m%ct_tilde  = 0.5*COS(m%NacYaw)**2*SIN(m%NacYaw)*y%avg_ct 
    m%ct_tilde = y%avg_ct
           
END SUBROUTINE get_initial_condition

!-------------------------------------------------------------------------
SUBROUTINE calculate_wake(m, p, y, r_w, U_w, element_num, U, b)
!..................................................................................
! This routine is the main routine to calculate the wake
! This routine is called after receiving the scaled rotor radius and the scaled velocity at the rotor
! The output of this routine is the wake velocity which is "U"
!   and the wake width which is the "b"
!..................................................................................
    !USE   DWM_Wake_Deficit_Data
    !USE   DWM_ParameterType,     ONLY: p%p_p_r, p%r_domain, p%Uambient, p%x_domain, p%hub_height, p%TI_amb
    
    TYPE(DWM_MiscVarType),         INTENT(INOUT)  :: m
    TYPE(DWM_ParameterType),       INTENT(IN   )  :: p           ! Parameters
    TYPE(DWM_OutputType),          INTENT(INOUT)  :: y
    
    INTEGER              ::    element_num
    REAL(ReKi)                 ::    r_w (element_num)        ! scaled rotor radius r_w
    REAL(ReKi)                 ::    U_w (element_num)        ! scaled velocity  U_w
    REAL(ReKi),ALLOCATABLE     ::    U(:,:)
    INTEGER,ALLOCATABLE     ::    b(:)
    
    ! local variables    
    REAL(ReKi)   ::  mtemp
    REAL(ReKi)   ::  ntemp
    REAL(ReKi)   ::  xtemp
    REAL(ReKi)   ::  ytemp
    
    REAL(ReKi), DIMENSION(2)      ::   filter1
    REAL(ReKi), DIMENSION(2)      ::   filter2
    REAL(ReKi), ALLOCATABLE       ::   F1_vector (:)
    REAL(ReKi), ALLOCATABLE       ::   F2_vector (:)
    !REAL(ReKi)                    ::   m%DWDD%ppR                           ! Point_per_R_resoulution
    REAL(ReKi)                    ::   Domain_R                      ! Domain_size_in_radial_direction
    REAL(ReKi)                    ::   Domain_X                      ! Domain_size_in_flow_direction
    !REAL(ReKi)                    ::   TI_original                ! Turbulence_intensity normalized back to ambient wind speed
    REAL(ReKi)                    ::   k1                            ! Amb turb. coeff.
    REAL(ReKi)                    ::   k2                            ! Shear layer coeff.

    !INTEGER                 ::   m%DWDD%n_x_vector
    !INTEGER                 ::   m%DWDD%n_r_vector
    
    !%%%%% Rolf modification
    INTEGER                 ::   length_F1_vector
    REAL(ReKi)                    ::   L_ABL_vector(3)
    REAL(ReKi)                    ::   UW_UU_vector(3)
    REAL(ReKi)                    ::   L_DEF_vector(3)
    REAL(ReKi)                    ::   UU_DEF_UU_ABL_vector(3)
    REAL(ReKi)                    ::   UW_DEF_UU_DEF_vector(3)
    REAL(ReKi)                    ::   x_ary(3)
    REAL(ReKi)                    ::   L_ABL
    REAL(ReKi)                    ::   UW_UU
    REAL(ReKi)                    ::   L_DEF
    REAL(ReKi)                    ::   UU_DEF_UU_ABL
    REAL(ReKi)                    ::   UW_DEF_UU_DEF
    REAL(ReKi)                    ::   Rotor_fixed_R
    REAL(ReKi)                    ::   l_star_ABL
    REAL(ReKi)                    ::   l_star_DEF
    REAL(ReKi)                    ::   UU_DEF_UU_ABL_fac
    REAL(ReKi)                    ::   u_star_ABL
    REAL(ReKi)                    ::   u_star_DEF
    REAL(ReKi)                    ::   Shear_add_du_dz
    REAL(ReKi),ALLOCATABLE        ::   visc_wake(:,:)
    REAL(ReKi),ALLOCATABLE        ::   visc_wake1(:,:) 
    REAL(ReKi),ALLOCATABLE        ::   visc_wake2(:,:) 
    REAL(ReKi)                    ::   visc_norm_factor
    REAL(ReKi),ALLOCATABLE        ::   alfa_1(:)
    REAL(ReKi),ALLOCATABLE        ::   alfa_2(:)
    REAL(ReKi),ALLOCATABLE        ::   du_dr_tot(:,:)
    INTEGER,ALLOCATABLE     ::   shear_flag(:)
    REAL(ReKi),ALLOCATABLE        ::   One_div_du_dr_DWM(:,:)
    REAL(ReKi),ALLOCATABLE        ::   visc_fac(:)

    REAL(ReKi)                    ::    R_WTG                        ! normalized radius
    REAL(ReKi)                    ::    U0                           ! normalized wind speed
    REAL(ReKi)                    ::    D_WTG                        ! normalized diameter
    REAL(ReKi)                    ::    R_length                     ! normalized length in radial direction
   !REAL(ReKi)                    ::    m%DWDD%X_length                     ! normalized length in axial direction
    INTEGER                 ::    np_r                         ! point per radial distance
    !INTEGER                 ::    m%DWDD%np_x                         ! point per axial distance
    REAL(ReKi)                    ::    delrad                       ! delta r
    REAL(ReKi)                    ::    delaxi                       ! delta x

    REAL(ReKi), ALLOCATABLE       ::    x_vector(:)
    REAL(ReKi), ALLOCATABLE       ::    r_vector(:)

    REAL(ReKi), ALLOCATABLE       ::    V(:,:)
    REAL(ReKi), ALLOCATABLE       ::    visc(:,:)
    REAL(ReKi), ALLOCATABLE       ::    visc_DWM(:,:)
    REAL(ReKi), ALLOCATABLE       ::    du_dr_DWM(:,:)
    REAL(ReKi), ALLOCATABLE       ::    du_dr_total(:,:)
   !REAL(ReKi), ALLOCATABLE       ::    m%DWDD%Turb_Stress_DWM(:,:)
   !REAL(ReKi), ALLOCATABLE       ::    TI_DWM(:,:)
   !REAL(ReKi), ALLOCATABLE       ::    U_face(:,:)
   !REAL(ReKi), ALLOCATABLE       ::    VOL_x_jhigh(:,:)
   !REAL(ReKi), ALLOCATABLE       ::    VOL_x_jlow (:,:)
   !REAL(ReKi), ALLOCATABLE       ::    VOL_r_ihigh(:,:)
   !REAL(ReKi), ALLOCATABLE       ::    VOL_r_ilow (:,:)
    REAL(ReKi), ALLOCATABLE       ::    r_vec_DWM (:)
    REAL(ReKi), ALLOCATABLE       ::    dA_DWM (:)

    INTEGER                       ::    n_r_vec_DWM
    INTEGER                       ::    b_loop
    INTEGER                       ::    b_counter
    REAL(ReKi)                    ::    dr_DWM
    REAL(ReKi)                    ::    Def_DWM
    REAL(ReKi)                    ::    Def_DWM_mixL
    REAL(ReKi)                    ::    A_total
    REAL(ReKi)                    ::    k_wiener

    INTEGER, ALLOCATABLE          ::    counter(:)
    INTEGER                       ::    i
    INTEGER                       ::    j
    INTEGER                       ::    k
    INTEGER                       ::    ILo
    INTEGER                       ::    NumEqu
    INTEGER                       ::    n_xi
    INTEGER                       ::    n_U_tmp_2

    REAL(ReKi), ALLOCATABLE       ::    bin_filter(:)
    REAL(ReKi), ALLOCATABLE       ::    xi(:)
    REAL(ReKi), ALLOCATABLE       ::    U_tmp_1(:)
    REAL(ReKi), ALLOCATABLE       ::    U_tmp_2(:)
    REAL(ReKi), ALLOCATABLE       ::    U_tmp(:)
    REAL(ReKi), ALLOCATABLE       ::    mat(:,:)
    REAL(ReKi), ALLOCATABLE       ::    RHS(:)
    REAL(ReKi), ALLOCATABLE       ::    Soln(:)
    REAL(ReKi), ALLOCATABLE       ::    AugMat(:,:)

    REAL(ReKi)                    ::    LHS1
    REAL(ReKi)                    ::    LHS2
    REAL(ReKi)                    ::    LHS3
    REAL(ReKi)                    ::    LHS11
    REAL(ReKi)                    ::    LHS12
    REAL(ReKi)                    ::    LHS13
    REAL(ReKi)                    ::    LHS21
    REAL(ReKi)                    ::    LHS22
    REAL(ReKi)                    ::    LHS23
    REAL(ReKi)                    ::    LHS31
    REAL(ReKi)                    ::    LHS41
    REAL(ReKi)                    ::    LHS32
    REAL(ReKi)                    ::    LHS33
    REAL(ReKi)                    ::    LHS43
    
    REAL(ReKi),ALLOCATABLE        ::    main_diagonal(:)
    REAL(ReKi),ALLOCATABLE        ::    sub_diagonal(:)
    REAL(ReKi),ALLOCATABLE        ::    sup_diagonal(:)
      
    m%DWDD%ppR       = p%p_p_r
    Domain_R  = p%r_domain       !10.0  domain size in R [R]
    Domain_X  = p%x_domain      !42.0  domain size in X [R]
    filter1   = (/0.0,   4.0 /)
    filter2   = (/0.035, 0.35/)
    k1        = 0.0919  
    k2        = 0.0178
    R_WTG     = 1.0
    U0        = 1.0
    D_WTG     = 2.0
    R_length  = Domain_R
    m%DWDD%X_length  = Domain_X
    
    m%TI_original = y%TI*(y%Mean_FFWS/p%Uambient)    ! calculate the TI if under ambient wind speed

    np_r       =  m%DWDD%ppR          ! per R   ie. R resolution is 50
    m%DWDD%np_x       =  m%DWDD%ppR          ! per D   ie. X resolution is 50
    delrad     =  R_WTG/np_r   ! dr
    delaxi     =  D_WTG/m%DWDD%np_x   ! dx
    m%DWDD%n_x_vector =  floor((m%DWDD%X_length)/D_WTG*m%DWDD%np_x)      ! number of point in equally spaced array x_vector
    m%DWDD%n_r_vector =  floor((R_length)/R_WTG*np_r)      ! number of point in equally spaced array r_vector

    ! create coordinate vectors
    ALLOCATE    (x_vector(m%DWDD%n_x_vector))                             
    ALLOCATE    (r_vector(m%DWDD%n_r_vector))
    ! similar to linspace function
    x_vector   = ( (m%DWDD%X_length-delaxi)/(m%DWDD%n_x_vector-1 ) )*[(i,i=1,m%DWDD%n_x_vector)]+(0-( (m%DWDD%X_length-delaxi)/(m%DWDD%n_x_vector-1) ))        
    r_vector   = ( (R_length-delrad)/(m%DWDD%n_r_vector-1 ) )*[(i,i=1,m%DWDD%n_r_vector)]+(0-( (R_length-delrad)/(m%DWDD%n_r_vector-1) )) 
    
    ! Create the F1 filter
    
    CALL create_F1_filter (F1_vector, filter1, length_F1_vector,m%DWDD%np_x,m%DWDD%X_length)
         !OPEN (unit=25,file="DWM\results\F1_filter.txt")
         !WRITE (25,'(f13.6)'), F1_vector(:)
         !CLOSE(25)
         
    CALL create_F2_filter (F2_vector, filter2, m%DWDD%np_x, length_F1_vector)
    

    !CALL create_filter_vector ( filter1, F1_vector )
    
    !CALL create_filter_vector ( filter2, F2_vector )
    
         !OPEN (unit=25,file="DWM\results\F2_filter.txt")
         !WRITE (25,'(f13.6)'), F2_vector(:)
         !CLOSE(25)
  
     ! Initiate the U, V, visc, TI_add, Turb_stress matrices
    ALLOCATE (V               (m%DWDD%n_x_vector,m%DWDD%n_r_vector))
    ALLOCATE (U               (m%DWDD%n_x_vector,m%DWDD%n_r_vector))    !axial velocity matrix
    ALLOCATE (visc            (m%DWDD%n_x_vector,m%DWDD%n_r_vector))
    ALLOCATE (visc_DWM        (m%DWDD%n_x_vector,m%DWDD%n_r_vector))
    ALLOCATE (du_dr_DWM       (m%DWDD%n_x_vector,m%DWDD%n_r_vector))
    ALLOCATE (du_dr_total     (m%DWDD%n_x_vector,m%DWDD%n_r_vector))
    ALLOCATE (m%DWDD%Turb_Stress_DWM (m%DWDD%n_x_vector,m%DWDD%n_r_vector))
    !ALLOCATE (TI_DWM          (m%DWDD%n_x_vector,m%DWDD%n_r_vector))
    !ALLOCATE (U_face          (m%DWDD%n_x_vector,m%DWDD%n_r_vector-1))
    !ALLOCATE (VOL_x_jhigh     (m%DWDD%n_x_vector,m%DWDD%n_r_vector-1))
    !ALLOCATE (VOL_x_jlow      (m%DWDD%n_x_vector,m%DWDD%n_r_vector-1))
    !ALLOCATE (VOL_r_ihigh     (m%DWDD%n_x_vector,m%DWDD%n_r_vector-1))
    !ALLOCATE (VOL_r_ilow      (m%DWDD%n_x_vector,m%DWDD%n_r_vector-1))
    ALLOCATE (visc_wake         (m%DWDD%n_x_vector,m%DWDD%n_r_vector))
    ALLOCATE (visc_wake1        (m%DWDD%n_x_vector,m%DWDD%n_r_vector))
    ALLOCATE (visc_wake2        (m%DWDD%n_x_vector,m%DWDD%n_r_vector))
    ALLOCATE (alfa_1            (m%DWDD%n_r_vector           ))
    ALLOCATE (alfa_2            (m%DWDD%n_r_vector           ))
    ALLOCATE (du_dr_tot         (m%DWDD%n_x_vector,m%DWDD%n_r_vector))
    ALLOCATE (shear_flag        (m%DWDD%n_r_vector           )) 
    ALLOCATE (One_div_du_dr_DWM (m%DWDD%n_x_vector,m%DWDD%n_r_vector))
    ALLOCATE (visc_fac          (m%DWDD%n_r_vector           ))
    
    ALLOCATE (main_diagonal     (m%DWDD%n_r_vector           ))
    ALLOCATE (sub_diagonal      (m%DWDD%n_r_vector           ))
    ALLOCATE (sup_diagonal      (m%DWDD%n_r_vector           ))

    V                  = 0
    U                  = 0
    visc               = 0
    du_dr_DWM          = 0
    m%DWDD%Turb_Stress_DWM    = 0
    !TI_DWM             = 0
    !U_face             = 0
    !VOL_x_jhigh        = 0
    !VOL_x_jlow         = 0
    !VOL_r_ihigh        = 0
    !VOL_r_ilow         = 0

    !%%%% BOUNDARY CONDITIONS
    ! ROTOR PLANE
    ALLOCATE (bin_filter(m%DWDD%n_r_vector))
    DO i=1,m%DWDD%n_r_vector
      IF (MAXVAL(r_w)>r_vector(i)) THEN
      bin_filter(i) = 1
      ELSE
      bin_filter(i) = 0
      END IF
    END DO


    n_xi=floor(sum(bin_filter))
    ALLOCATE (xi(n_xi))
    xi=r_vector(1:n_xi)*bin_filter(1:n_xi)

    ALLOCATE (U_tmp_1(n_xi))
    ILo = 1
    DO i=1,n_xi
       U_tmp_1(i) = InterpBin( xi(i), r_w, U_w, ILo, size(r_w))
    END DO

    n_U_tmp_2 = size(r_vector)-n_xi 
    ALLOCATE (U_tmp_2(n_U_tmp_2))
    U_tmp_2 = U0

    ALLOCATE (U_tmp(n_xi +n_U_tmp_2))    !  =m%DWDD%n_r_vector
    U_tmp = (/U_tmp_1, U_tmp_2/)
    U (1,:) = U_tmp(1:m%DWDD%n_r_vector)

    ! Centerline
    V (1,:) = 0
    

       
    !%%%% SOLVING FLOW FIELD
    ALLOCATE  (b(m%DWDD%n_x_vector))
    ALLOCATE  (counter(m%DWDD%n_x_vector))
              counter=1
    ALLOCATE (AugMat (m%DWDD%n_r_vector,m%DWDD%n_r_vector+1) )
    ALLOCATE (Soln   (m%DWDD%n_r_vector)              )

    ALLOCATE (mat(m%DWDD%n_r_vector,m%DWDD%n_r_vector))
    ALLOCATE (RHS(m%DWDD%n_r_vector))

    n_r_vec_DWM = floor(m%DWDD%ppR*Domain_R)
    ALLOCATE ( r_vec_DWM (n_r_vec_DWM) )
    ALLOCATE ( dA_DWM (n_r_vec_DWM-1)  )
    
    !%%%%%%%%%%%%%%%  Atmospheric stability effects  %%%%%%%%%%%%%%%%%%
    L_ABL_vector         = (/26.5352,      34.026,      40.7458/)
    UW_UU_vector         = (/-0.27359,    -0.27887,    -0.27935/)
    L_DEF_vector         = (/11.065,       12.9746,     14.4395/)
    UU_DEF_UU_ABL_vector = (/0.63044,      0.57982,      0.5287/)
    UW_DEF_UU_DEF_vector = (/-0.27341,    -0.25684,    -0.24217/)
    
    !%%%%% interpolation wrt to the hub height
    x_ary = (/40,100,160/)
    L_ABL         = InterpBin( p%hub_height, x_ary, L_ABL_vector, ILo, size(x_ary) ) !int_Lww(k3) i.e (Integral length scale (in vertical directions), from ww(k3)) 
    UW_UU         = InterpBin( p%hub_height, x_ary, UW_UU_vector, ILo, size(x_ary) ) ! ratio of UW and UU stresses for whole spectra   
    L_DEF         = InterpBin( p%hub_height, x_ary, L_DEF_vector, ILo, size(x_ary) ) ! Integral length scale (in vertical directions), Meandering length scale subtracted, from ww(k3)   
    UU_DEF_UU_ABL = InterpBin( p%hub_height, x_ary, UU_DEF_UU_ABL_vector, ILo, size(x_ary) ) ! Part of normal stress in the deficit module   
    UW_DEF_UU_DEF = InterpBin( p%hub_height, x_ary, UW_DEF_UU_DEF_vector, ILo, size(x_ary) ) ! ratio of UW and UU stresses for spectra in deficit scales
    
    !%%%%% normalized by the fixed rotor R
    Rotor_fixed_R  = 40 !ATMOSTAB ANALYSIS IS CARRIED OUT OVER R = 40m, which should be used to normalize the length scales
    l_star_ABL     = L_ABL / Rotor_fixed_R;
    l_star_DEF     = L_DEF / Rotor_fixed_R;
    
    !%%%%% Normalize UU_160m to neutral condition
    UU_DEF_UU_ABL_fac = InterpBin( p%hub_height, x_ary, (/0.63044_ReKi, 0.57982_ReKi, 0.5287_ReKi/), ILo, size(x_ary) )
    UU_DEF_UU_ABL     = UU_DEF_UU_ABL / UU_DEF_UU_ABL_fac
    
    !%%%%% CALCULATE u* according to:
    ! 1. u* ~= (mean(u'w')^2 )^0.25 
    ! 2. {mean(u'w') = mean(u'u')*Cuw_uu} 
    ! 3. {u' ~= TI (in normalized form)}
    ! => u* ~= ((TI^2 * Cuw_uu )^2)^0.25 
    u_star_ABL     = ( (  (p%TI_amb/100)**2                      * ABS(UW_UU)         )**2 )**0.25
    u_star_DEF     = ( (  (m%TI_original/100)**2 * UU_DEF_UU_ABL * abs(UW_DEF_UU_DEF) )**2 )**0.25;
    
    Shear_add_du_dz = u_star_ABL / l_star_ABL
    
    !k_wiener = du_dz_ABL + delrad**2

    DO j=2, m%DWDD%n_x_vector, 1  ! start from the plane next to the rotor plane

    !==== Calculating wake width "b" where 95% of the deficit is captured
      dr_DWM = 1.0/m%DWDD%ppR
  
       DO i=1,n_r_vec_DWM    ! build r_vec_DWM
          r_vec_DWM (i) = dr_DWM/2 + (i-1)*dr_DWM
       END DO
   
       DO i=1,n_r_vec_DWM-1    ! build dA_DWM
          dA_DWM (i) = Pi * ( r_vec_DWM (i+1)**2 - r_vec_DWM (i)**2 )
       END DO
   
       Def_DWM = 0                             ! Calculate Def_DWM 
       DO i=1, n_r_vec_DWM-1
          Def_DWM = (1-U (j-1,i+1)) * dA_DWM (i) + Def_DWM
       END DO
   
       Def_DWM_mixL = 0                        ! Calculate the wake width "b"
       DO i = 2,NINT( m%DWDD%ppR )
          Def_DWM_mixL = ( 1- U(j-1,i) ) * dA_DWM(i) + Def_DWM_mixL
       END DO
       DO b_counter = NINT(m%DWDD%ppR)+1, n_r_vec_DWM-1          
          Def_DWM_mixL = ( 1- U(j-1,b_counter) ) * dA_DWM(b_counter) + Def_DWM_mixL
          IF ( Def_DWM_mixL > Def_DWM * 0.99 ) THEN
             EXIT
          ELSE IF (b_counter == n_r_vec_DWM-1) THEN
             EXIT
          END IF
       END DO
       b(j-1) = b_counter
 
       ! %%%%% Calculate eddy viscosity
       ! Include blend between original Prandtl model and Ainslie to avoid issues when wake turbulence goes to 0.
       ! The largest eddy viscosity at each point is applied.
       
       ! Calculate mean flow gradient - du/dr is created with CDS (apart from 1st and last point)
       du_dr_DWM(j-1,1)                  = (U(j-1,2) - U(j-1,1))/delrad
       du_dr_DWM(j-1,2:NINT(R_length*np_r)-1)      = ( U(j-1,3:NINT(R_length*np_r-1)+1) - U(j-1,1:NINT(R_length*np_r-1)-1))/(2*delrad)
       du_dr_DWM(j-1,  NINT(R_length*np_r)  )      = ( U(j-1,  NINT(R_length*np_r  )  ) - U(j-1,  NINT(R_length*np_r-1)  ))/delrad
       
       ! %%% Blend of mixL and Ainslie eddy visc
       DO I = 1,m%DWDD%n_r_vector
           visc_wake1(j-1,I)     = F2_vector(j-1)* k2 *( r_vector(b(j-1))/R_WTG )**2 * ABS(du_dr_DWM(j-1,I));
           visc_wake2(j-1,I)     = F2_vector(j-1)* k2 *( r_vector(b(j-1))/R_WTG )    * ( 1 - min_of_array( U(j-1,:),SIZE(U(j-1,:)) ) );
           visc_wake (j-1,I)     = max( visc_wake1(j-1,I),visc_wake2(j-1,I) );
       END DO
       
       ! %%% Atmospheric eddy visc as u*l*, yields total eddy viscosity
       visc_norm_factor = 6.3918
       DO I = 1,m%DWDD%n_r_vector
           visc(j-1,I)           = F1_vector(j-1) * k1 * visc_norm_factor * u_star_DEF * l_star_DEF + visc_wake(j-1,I);
       END DO
       
       ! %%%%% Include contribution from atmospheric boundary layer on DWM
       ! % 1. Calculate the azimuthally averaged local gradient (du/dr tot) acting of the eddy viscosity as a combination of du/dr in the DWM model and du/dz from ABL
       ! % 2. The du/dr contribution is constant in azimuthal direction. The du/dz part is assumed linear, which gives a sinus curve in a du/dr system
       
       !% Calculate total mean flow gradient - adds shear contribution via
       !% sin function. This gets the stresses right, but sign is wrong in
       !% regions where du/dr_DWM - sign of du/dz_ABL is negative 
       
       DO I = 1,m%DWDD%n_r_vector
          !alfa_1(I)      = ASIN(ABS(du_dr_DWM(j-1,I)) / Shear_add_du_dz)
          !alfa_2(I)      = Pi - alfa_1(I)
          
          ! % condition for added shear gradient (if du/dr_DWM >= du/dz_ABL there are no contribution)
          IF ( ABS(du_dr_DWM(j-1,I)) < Shear_add_du_dz ) THEN
              shear_flag(I) = 1
              alfa_1(I)      = ASIN(ABS(du_dr_DWM(j-1,I)) / Shear_add_du_dz)
              alfa_2(I)      = Pi - alfa_1(I)
          ELSE
              shear_flag(I) = 0
              alfa_1(I)      = 0
              alfa_2(I)      = 0
          END IF
          
          
          du_dr_tot(j-1,I) = (  ABS(du_dr_DWM(j-1,I))*2*Pi + shear_flag(I)*2*&
                                ( Shear_add_du_dz*2*COS(alfa_1(I))-ABS(du_dr_DWM(j-1,I))*(alfa_2(I) - alfa_1(I))) )/(2*Pi)
       END DO
       
       ! %%% Use "wiener filter" for numerical stability:  1/f(x) ~= f(x) / (f(x)^2 + k)
       k_wiener                 = 2*Shear_add_du_dz * delrad**2;
       DO I = 1,m%DWDD%n_r_vector
           One_div_du_dr_DWM(j-1,I) = du_dr_DWM(j-1,I) / (du_dr_DWM(j-1,I)**2 + k_wiener)
           visc_fac(I)              = max(1.0_ReKi, (du_dr_tot(j-1,I) * ABS(One_div_du_dr_DWM(j-1,I))))
           visc(j-1,I)              = visc(j-1,I) * visc_fac(I)
       END DO
       
          
       
       
       
       
       
       
       
       !!!DO i = 1,m%DWDD%n_r_vector
        !!!   IF ( ABS(du_dr_DWM(j-1,i)) >= du_dz_ABL ) THEN
         !! !     A_total = 2*Pi*du_dr_DWM(j-1,i)
          !!! ELSE
           !!!    ytemp = du_dr_DWM(j-1,i)
            !!!   xtemp = 2 * shear_correction( du_dr_DWM(j-1,i),du_dz_ABL )
             !!!  A_total = 2*Pi*du_dr_DWM(j-1,i) + 2 * shear_correction( du_dr_DWM(j-1,i),du_dz_ABL )
           !!!END IF
           !!!du_dr_total(j-1,i) = A_total / (2*Pi)
       !!!END DO
   
       !!!visc_DWM(j-1,:) = F1_vector(j-1)*k1*(TI_original/100) + F2_vector(j-1)* k2 *( r_vector(b(j-1))/R_WTG )**2 * ABS( du_dr_DWM(j-1,:) )
       !!!visc(j-1,:) = visc_DWM(j-1,:)
       
       
       !DO i = 1,m%DWDD%n_r_vector
           !mtemp = du_dr_total(j-1,i)
           !ntemp = du_dr_DWM(j-1,i)
           !IF (ABS( du_dr_DWM(j-1,i) < 0.0001 ) ) THEN              
            !  visc(j-1,i) = visc_DWM(j-1,i) * du_dr_total(j-1,i) * ABS( du_dr_DWM(j-1,i)/(du_dr_DWM(j-1,i)**2 + k_wiener) )
           !ELSE
            !  visc(j-1,i) = visc_DWM(j-1,i) * du_dr_total(j-1,i) / ABS( du_dr_DWM(j-1,i) )
           !END IF           
           !visc(j-1,i) = visc_DWM(j-1,i) * du_dr_total(j-1,i) * ABS( 1/du_dr_DWM(j-1,i))
       !END DO

       mat=0
       ! ====SHORT INSTRUCIONS TO SOLVE RUTINE:
       ! The terms LHS and RHS (left/right hand side) refers to the terms of
       ! the coefficient matrix developed to solve the then shear layer
       ! approximation of NS. The numbers indicate the position in the equation,
       ! ex LHS21 is the 2nd part of the 1st term on the left side in eq.2.8,
       ! see document "Numerical implementation of DWM deficit module" for details.
   
       ! Input BC for wake center
   
       LHS2       = U(j-1,1)/delaxi     + (2*visc(j-1,1)/(delrad**2))
       LHS3       = -(2*visc(j-1,1) /(delrad**2))
       RHS(1)     = (U(j-1,1)**2    / delaxi)
       mat(1,1)   = LHS2
       mat(2,1)   = LHS3
   
       ! Calculation of U for the wake body
       DO i=2,(m%DWDD%n_r_vector-1),1      ! starts from the point next to the hub center
          LHS11             = -V(j-1,i)      / (2*delrad)
          LHS21             = visc(j-1,i)    / (2*r_vector(i)*delrad)
          LHS31             = -visc(j-1,i)   / (delrad**2)
          LHS41             = (visc(j-1,i+1) - visc(j-1,i-1))  / (2*delrad)**2;  ! new term due to d(nu_t)/dr dependence
          LHS12             = U(j-1,i)       / (delaxi)
          LHS22             = 2*visc(j-1,i)  / (delrad**2)
          LHS13             = V(j-1,i)       / (2*delrad)
          LHS23             = -visc(j-1,i)   / (2*r_vector(i)*delrad)
          LHS33             = -visc(j-1,i)   / (delrad**2)
          LHS43             = -(visc(j-1,i+1) - visc(j-1,i-1))  / (2*delrad)**2; ! new term due to d(nu_t)/dr dependence
          LHS1              = LHS11 + LHS21 + LHS31 + LHS41
          LHS2              = LHS12 + LHS22
          LHS3              = LHS13 + LHS23 + LHS33 + LHS43
          RHS(i)            = (U(j-1,i)**2  / delaxi)
          ! Build the matrix for X =A/B
          mat(i-1,i) = LHS1
          mat(i  ,i) = LHS2
          mat(i+1,i) = LHS3
       END DO
   
       ! Input BC for wake edge
       LHS1                     = 0
       LHS2                     = 1/delaxi
       RHS(NINT(R_length*np_r))       = (U(j-1,NINT(R_length*np_r))/ delaxi)
       mat(NINT(R_length*np_r)-1,    NINT(R_length*np_r) )     = LHS1
       mat(NINT(R_length*np_r)  ,    NINT(R_length*np_r) )     = LHS2
   
       ! Solve for the U
       ! Use Gauss-Jordan elimination
       AugMat (1:m%DWDD%n_r_vector, 1:m%DWDD%n_r_vector) = TRANSPOSE(mat)
       AugMat (:           , m%DWDD%n_r_vector+1) = RHS
       NumEqu                              = m%DWDD%n_r_vector
       !CALL Gauss(AugMat, NumEqu, Soln)
       !U(j,:)=Soln
       
       ! === USE Thomas Algorithm to solve the matrix ====  6.30.2014
       main_diagonal (1) = AugMat(1,1)
       sub_diagonal  (1) = 0               ! means it is the diagonal below the main diagonal
       sup_diagonal  (1) = AugMat(1,2)     ! means it is the diagonal above the main diagonal
      
       DO I = 2,m%DWDD%n_r_vector-1
           main_diagonal (I) = AugMat(I,I)
           sub_diagonal  (I) = AugMat(I,I-1)
           sup_diagonal  (I) = AugMat(I,I+1)
       END DO
       
       main_diagonal (m%DWDD%n_r_vector) = AugMat(m%DWDD%n_r_vector, m%DWDD%n_r_vector)
       sub_diagonal  (m%DWDD%n_r_vector) = AugMat(m%DWDD%n_r_vector, m%DWDD%n_r_vector-1)
       sup_diagonal  (m%DWDD%n_r_vector) = 0
       
       CALL Thomas_diagonal (sub_diagonal, main_diagonal, sup_diagonal, RHS, Soln, NumEqu)
       U(j,:)=Soln
            
       ! === Solve for V  
       DO i = 1,NINT(R_length)*np_r-1,1
         V(j,i+1) = (r_vector(i) / r_vector(i+1)) * V(j,i) -(delrad/(2*delaxi))*( (U(j,i+1) - U(j-1,i+1)) + &
                (r_vector(i) / r_vector(i+1)) * ((U(j,i) - U(j-1,i))) )
       END DO
   
       ! POST PROCESSING SIGNAL: Turbulent stress
       DO i=1,m%DWDD%n_r_vector,1
         !m%DWDD%Turb_Stress_DWM(j-1,i) = visc_DWM(j-1,i) * du_dr_total(j-1,i)
         m%DWDD%Turb_Stress_DWM(j-1,i) = visc(j-1,i) * du_dr_DWM(j-1,i)
       END DO
   
       ! Control calculatoins of mass flux over cells
  
       !!DO i=1,m%DWDD%n_r_vector-1,1
          !! VOL_x_jhigh(j-1,i)  = (Pi/3) *((U(j,i) *(r_vector(i+1)**3 - (3 *r_vector(i)**2 *r_vector(i+1)) + 2 *r_vector(i)**3)) +&
          !!                    (U(j,i+1) *(r_vector(i)**3 - (3*r_vector(i+1)**2 *r_vector(i)) + 2 *r_vector(i+1)**3)))/ delrad
          !! VOL_x_jlow(j-1,i)   = (Pi/3) *((U(j-1,i) *(r_vector(i+1)**3 - (3*r_vector(i)**2 *r_vector(i+1)) + 2 *r_vector(i)**3)) +&
            !!                   (U(j-1,i+1) *(r_vector(i)**3   - (3*r_vector(i+1)**2 *r_vector(i)) + 2 *r_vector(i+1)**3)))/ delrad
          !! VOL_r_ilow(j-1,i)   =  Pi * r_vector(i) * (V(j-1,i)+V(j,i)) * delaxi
      
          !! V(j,i+1)            = ((VOL_x_jlow(j-1,i) - VOL_x_jhigh(j-1,i) + VOL_r_ilow(j-1,i)) / (Pi*(r_vector(i+1)) *delaxi)) - V(j-1,i+1)     !! changed to version 2 2012/4/11
          !! VOL_r_ihigh(j-1,i)  = Pi * r_vector(i+1) * (V(j-1,i+1)+V(j,i+1)) * delaxi
          !! specificly U_face for mass flow and momentum calculations
          !! U_face(j-1,i) = VOL_x_jlow(j-1,i)  / (Pi*((r_vector(i+1)**2)-(r_vector(i))**2))
       !!END DO
   
    END DO

    b(m%DWDD%n_x_vector) = b(m%DWDD%n_x_vector-1)
    
    IF (ALLOCATED( V ))                  DEALLOCATE ( V )
    IF (ALLOCATED( visc ))               DEALLOCATE ( visc )
    IF (ALLOCATED( du_dr_DWM ))          DEALLOCATE ( du_dr_DWM )
    !IF (ALLOCATED( m%DWDD%Turb_Stress_DWM ))    DEALLOCATE ( m%DWDD%Turb_Stress_DWM )
    !IF (ALLOCATED( TI_DWM ))             DEALLOCATE ( TI_DWM )
    !IF (ALLOCATED( U_face ))             DEALLOCATE ( U_face )
    !IF (ALLOCATED( VOL_x_jhigh ))        DEALLOCATE ( VOL_x_jhigh )
    !IF (ALLOCATED( VOL_x_jlow ))         DEALLOCATE ( VOL_x_jlow )
    !IF (ALLOCATED( VOL_r_ihigh ))        DEALLOCATE ( VOL_r_ihigh )
    !IF (ALLOCATED( VOL_r_ilow ))         DEALLOCATE ( VOL_r_ilow )
    IF (ALLOCATED( r_vec_DWM ))          DEALLOCATE ( r_vec_DWM )
    IF (ALLOCATED( dA_DWM ))             DEALLOCATE ( dA_DWM )
    IF (ALLOCATED( bin_filter ))         DEALLOCATE ( bin_filter )
    IF (ALLOCATED( xi ))                 DEALLOCATE ( xi )
    IF (ALLOCATED( U_tmp_1 ))            DEALLOCATE ( U_tmp_1 )
    IF (ALLOCATED( U_tmp_2 ))            DEALLOCATE ( U_tmp_2 )
    IF (ALLOCATED( U_tmp ))              DEALLOCATE ( U_tmp )
    IF (ALLOCATED( mat ))                DEALLOCATE ( mat )
    IF (ALLOCATED( U_tmp ))              DEALLOCATE ( U_tmp )
    IF (ALLOCATED( RHS ))                DEALLOCATE ( RHS )
    IF (ALLOCATED( Soln ))               DEALLOCATE ( Soln )
    IF (ALLOCATED( AugMat ))             DEALLOCATE ( AugMat )
    
    !OPEN(unit = 10, status='replace',file='sizeof_Uvelocity_2nd.bin',form='unformatted')   ! create sizeof_Uvelocity_2nd.bin, or overwrite an existing on
    !WRITE(10)   m%DWDD%ppR,Domain_R                                                               ! write the length of the velocity vector                                                                                                                                                                      
    !CLOSE(10)
    
    !OPEN(unit = 10, status='replace',file='DWM\results\Uvelocity.bin',form='unformatted')          
    !WRITE(10)   U(floor(spacing_turbine * m%DWDD%ppR)+1,:)  ! write the wind data of the plane where the downstream turbine locates                                                                                                                                                                                                                                                                                                    
    !CLOSE(10)
    
    !OPEN (unit=25,file="DWM\results\wake_width.txt")
    !WRITE (25,'(I5)'), b(:)
    !close(25)

    
END SUBROUTINE calculate_wake

!---------------------------------------------------------------------------------------------
SUBROUTINE Thomas_diagonal (lowerDia, mainDia, upperDia, RightHS, SolnVec, NumEq)
!.............................................................................................
! This function returns the F1 filter function
!.............................................................................................

    INTEGER  ::   NumEq
    
    REAL(ReKi)     ::   lowerDia(NumEq)
    REAL(ReKi)     ::   mainDia(NumEq)
    REAL(ReKi)     ::   upperDia(NumEq)
    REAL(ReKi)     ::   RightHS(NumEq)
    REAL(ReKi)     ::   SolnVec(NumEq)
    
    REAL(ReKi)     ::   cp_vec(NumEq)
    REAL(ReKi)     ::   dp_vec(NumEq)
    REAL(ReKi)     ::   temp
    INTEGER  ::   I
    
    ! initialize c-prime and d-prime
    cp_vec(1) = upperDia(1) / mainDia(1)
    dp_vec(1) = RightHS(1)  / mainDia(1)
    
    ! solve for vectors c-prime and d-prime
    DO I = 2,NumEq
        temp = mainDia(i) - cp_vec(i-1)*lowerDia(i)
        cp_vec(i) = upperDia(i)/temp
        dp_vec(i) = (RightHS(i)-dp_vec(i-1)*lowerDia(i))/temp
    END DO
    
    ! initialized SolnVec
    SolnVec(NumEq) = dp_vec(NumEq)
    
    ! solve for x from the vectors c-prime and d-prime
    DO I = NumEq-1, 1,-1
        SolnVec(i) = dp_vec(i) - cp_vec(i)*SolnVec(i+1)
    END DO
    
END SUBROUTINE Thomas_diagonal

!---------------------------------------------------------------------------------------------
SUBROUTINE create_F1_filter (F1_vector, filter1, length_F1_vector,np_x,X_length)
!.............................................................................................
! This function returns the F1 filter function
!.............................................................................................
    REAL(ReKi),ALLOCATABLE   ::     F1_vector(:)
    REAL(ReKi)               ::     filter1(2)
    INTEGER                  ::     length_F1_vector
    INTEGER                  ::     np_x
    REAL(ReKi)               ::     X_length
    
    INTEGER                  ::     length_F1_vector_1
    INTEGER                  ::     length_F1_vector_2
    REAL(ReKi),ALLOCATABLE   ::     F1_vector_1(:)
    REAL(ReKi),ALLOCATABLE   ::     F1_vector_2(:)
    INTEGER                  ::     I
    
    length_F1_vector_1 = floor(filter1(2)*np_x/2)
    length_F1_vector_2 = floor(X_length*np_x/2)
    length_F1_vector   = length_F1_vector_1 + length_F1_vector_2
    ALLOCATE    (F1_vector_1(length_F1_vector_1))
    ALLOCATE    (F1_vector_2(length_F1_vector_2))
    ALLOCATE    (F1_vector  (length_F1_vector  ))
    
    F1_vector_1 = ( (1-filter1(1))   /(length_F1_vector_1-1 ) )*[(i,i=1,length_F1_vector_1)]+(0-( (1-filter1(1) )/(length_F1_vector_1-1 ) ))
    !r_vector    = ( (R_length-delrad)/(n_r_vector-1 ) )*[(i,i=1,n_r_vector)]+(0-( (R_length-delrad)/(n_r_vector-1) ))
    F1_vector_2 = 1
    
    F1_vector   = (/F1_vector_1,F1_vector_2/)
  
END SUBROUTINE create_F1_filter

!---------------------------------------------------------------------------------------------
SUBROUTINE create_F2_filter (F2_vector, filter2, np_x, length_F1_vector)
!.............................................................................................
! This function returns the F2 filter function
!.............................................................................................
    REAL(ReKi),ALLOCATABLE   ::    F2_vector(:)
    REAL(ReKi)               ::    filter2(2)
    INTEGER            ::    np_x
    INTEGER            ::    length_F1_vector
    
    REAL(ReKi),ALLOCATABLE   ::    F2_vector_x(:)
    REAL(ReKi),ALLOCATABLE   ::    F2_vector_1(:)
    REAL(ReKi),ALLOCATABLE   ::    F2_vector_2(:)
    INTEGER            ::    length_F2_vector_x
    INTEGER            ::    length_F2_vector_1
    INTEGER            ::    length_F2_vector_2
    INTEGER            ::    length_F2_vector
    INTEGER            ::    I
    
    length_F2_vector_x = floor(( REAL(length_F1_vector,ReKi) * (1/REAL(np_x,ReKi)) - (2+1/REAL(np_x,ReKi)) ) / (1/REAL(np_x,ReKi)) + 1)
    length_F2_vector_1 = 2*np_x
    length_F2_vector_2 = length_F2_vector_x
    length_F2_vector   = length_F2_vector_1 + length_F2_vector_2
    
    ALLOCATE ( F2_vector_x(length_F2_vector_x) )
    ALLOCATE ( F2_vector_1(length_F2_vector_1) )
    ALLOCATE ( F2_vector_2(length_F2_vector_2) )
    ALLOCATE ( F2_vector  (length_F2_vector  ) )
    
    F2_vector_x = ( (length_F1_vector * (1/REAL(np_x,ReKi)) - (2+1/REAL(np_x,ReKi)))   /(length_F2_vector_x-1 ) )*[(i,i=1,length_F2_vector_x)]+(2+1/REAL(np_x,ReKi)-( (length_F1_vector * (1/REAL(np_x,ReKi)) - (2+1/REAL(np_x,ReKi)))   /(length_F2_vector_x-1 ) ))
    F2_vector_1 = filter2(1)
    
    DO I = 1,length_F2_vector_2
        F2_vector_2(I) = 1-(1-filter2(1))*EXP(-filter2(2)*(F2_vector_x(I)-2))
    END DO
    
    F2_vector   = (/F2_vector_1,F2_vector_2/)
  
END SUBROUTINE create_F2_filter

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Gauss( AugMatIn, NumEq, SolnVec )
!..................................................................................................................................
! This routine uses the Gauss-Jordan elimination method for the solution of
!   a given set of simultaneous linear equations.
! NOTE: this routine works if no pivot points are zero and you don't want 
!   the eschelon or reduced eschelon form of the augmented matrix. The form of
!   the original augmented matrix IS preserved in this call.
!..................................................................................................................................
   IMPLICIT NONE


      ! Passed variables:

   INTEGER(4), INTENT(IN )      :: NumEq                                           ! Number of equations in augmented matrix.

   REAL(ReKi), INTENT(IN )      :: AugMatIn (NumEq,( NumEq + 1 ))                  ! Augmented matrix passed into this subroutine.
   REAL(ReKi), INTENT(OUT)      :: SolnVec  (NumEq)                                ! Solution vector.


      ! Local variables:

   REAL(ReKi)                   :: AugMat   (NumEq,( NumEq + 1 ))                  ! A copy of the augmented matrix.

   INTEGER(4)                   :: I                                               ! Steps through columns
   INTEGER(4)                   :: J                                               ! Steps through rows
   INTEGER(4)                   :: L                                               ! Steps through rows
   INTEGER(4)                   :: NAug                                            ! Column dimension of augmented matrix



      ! Transfer the data from AugMatIn to AugMat:

   AugMat = AugMatIn


      ! Find the column dimension of the augmented matrix:

   NAug = NumEq + 1


      ! Perform Gauss-Jordan elimination and store the solution vector
      !   in the last column of the augmented matrix:

   DO L = 1,NumEq             ! Loop through all rows
      DO I = ( L + 1 ), NAug  ! Loop through all columns above current row number
         AugMat(L,I) = AugMat(L,I) / AugMat(L,L)
         DO J = 1,NumEq       ! Loop through all rows except L
            IF ( J /= L )  AugMat(J,I) = AugMat(J,I) - ( AugMat(J,L)*AugMat(L,I) )
         ENDDO                ! J - All rows except L
      ENDDO                   ! I - All columns above current row number
   ENDDO                      ! L - All rows


      ! Transfer the solution vector from AugMat() to SolnVec():

   SolnVec = AugMat(:,NAug)



   RETURN
   
END SUBROUTINE Gauss  

!----------------------------------------------------------------------------------------------------------------------------------
FUNCTION shear_correction (du_dr_dwm,du_dz)
!..................................................................................................................................
! This function returns the shear correction factor A1 and A2
!..................................................................................................................................
   IMPLICIT NONE
    
   REAL(ReKi)                    ::  shear_correction
   REAL(ReKi)                    ::  du_dr_dwm
   REAL(ReKi)                    ::  du_dz             ! du_dz_abl
 
       ! Internal variables
   REAL(ReKi)                    ::  alpha_1                       
   REAL(ReKi)                    ::  alpha_2
   REAL(ReKi)                    ::  temp_integration
   REAL(ReKi)                    ::  correction_factor
   REAL(ReKi),ALLOCATABLE        ::  alpha_array(:)
   REAL(ReKi)                    ::  delta_alpha
   INTEGER(IntKi)                ::  I
   INTEGER(IntKi)                ::  temp_n
    
   alpha_1 = ASIN(du_dr_dwm/du_dz)
   alpha_2 = Pi/2 - alpha_1
   temp_integration = 0
    
   temp_n = 100
   ALLOCATE ( alpha_array (temp_n) )
    
   alpha_array = ((alpha_2-alpha_1)/(temp_n-1))*[(i,i=1,temp_n)]+(alpha_1-((alpha_2-alpha_1)/(temp_n-1)))
   delta_alpha = (alpha_2-alpha_1)/(temp_n-1)
    
   DO I = 1,temp_n
       temp_integration = du_dz * SIN(alpha_array(I)) * delta_alpha + temp_integration
   END DO
    
   shear_correction = temp_integration - (alpha_2 - alpha_1) * du_dr_dwm
 
   IF (ALLOCATED( alpha_array ))                DEALLOCATE ( alpha_array )
    
END FUNCTION shear_correction

!-------------------------------------------------------------------------------
FUNCTION filter_velocity (OS,m,p,u,x,xd,z,y,timestep,y_0,z_0,wake_radius)
!...............................................................................
! This function is called to calculate the filtered wake velocity
! The filter is a low pass filter
! The output is the filtered wake velocity at a certain wake center
!...............................................................................
    !USE DWM_Wake_Deficit_Data, ONLY: m%DWDD%ppR
    !USE filter_velocity_data
    USE InflowWind
    
    TYPE(DWM_OtherStateType),      INTENT(IN   )  :: OS
    TYPE(DWM_MiscVarType),         INTENT(INOUT)  :: m
    TYPE(DWM_ParameterType),       INTENT(IN   )  :: p
    TYPE(DWM_OutputType),          INTENT(INOUT)  :: y
    TYPE(DWM_InputType),           INTENT(INOUT)  :: u           ! An initial guess for the input; input mesh must be defined
    TYPE(DWM_ContinuousStateType), INTENT(INOUT)  :: x           ! Initial continuous states
    TYPE(DWM_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Initial discrete states
    TYPE(DWM_ConstraintStateType), INTENT(INOUT)  :: z           ! Initial guess of the constraint states   
    

    REAL(DbKi)   ::  timestep                     ! upper limit = usable time + grid width /  mean wind speed              %%% will change wrt wind speed  
                                                  ! = second / 0.05    timestep >= 1
    REAL(ReKi)         ::  y_0                          ! wake center point
    REAL(ReKi)         ::  z_0
    INTEGER      ::  wake_radius                  ! b(:) in cal_mixl
    REAL(ReKi)        ::  filter_velocity (3)          ! only v,w components
    
    INTEGER      ::  number_counter               ! counter : how many points are in the circle
    INTEGER      ::  radius_length                ! wake radius (meters)
    INTEGER      ::  y_axis                       
    INTEGER      ::  z_axis
    REAL(ReKi)         ::  temp_filter_velocity (3)     ! interpolation function, has u,v,w three components
    REAL(ReKi)         ::  temp_wind_velocity (3)
    
    INTEGER( IntKi )                   :: ErrStat           ! Error status of the operation
    CHARACTER                          :: ErrMsg            ! Error message if ErrStat /= ErrID_None


    !Print*, y_0
    
    temp_filter_velocity = 0.0
    number_counter       = 0
    !radius_length = NINT( wake_radius/m%DWDD%ppR*p%RotorR )             ! R(m): turbine radius
    radius_length = NINT(2*p%RotorR )
    
    IF (.NOT. ALLOCATED(u%IfW%PositionXYZ) ) THEN
       CALL AllocAry( u%IfW%PositionXYZ, 3, 1, "Position array to send to IfW_CalcOutput", ErrStat, ErrMsg )
       IF (ErrStat >= AbortErrLev)  RETURN
    END IF

    DO y_axis = NINT(y_0-radius_length),NINT(y_0+radius_length),1
       !IF (y_axis > p%WFLowerBd) THEN                                 ! add 9/25/2014
        DO z_axis = NINT(z_0-radius_length),NINT(z_0+radius_length),1
          IF ( z_axis > p%WFLowerBd )  THEN              !(make sure the circle does not exceed wind field)
            IF ( ((y_axis-y_0)**2+(z_axis-z_0)**2)**0.5 <= radius_length )  THEN
            
          
                 
            u%IfW%PositionXYZ(1,1) = 0.0
            u%IfW%PositionXYZ(2,1) = REAL(y_axis,ReKi)
            u%IfW%PositionXYZ(3,1) = REAL(z_axis,ReKi)                                         
            CALL InflowWind_CalcOutput( timestep, u%IfW, p%IfW, x%IfW, xd%IfW, z%IfW, OS%IfW, y%IfW, m%IfW, ErrStat, ErrMsg )
            temp_wind_velocity (:) = y%IfW%VelocityUVW(:,1)                          
          
            !temp_filter_velocity(:) = temp_filter_velocity(:) + AD_WindVelocityWithDisturbance(  REAL(timestep,ReKi), A_u, A_p, A_x, A_xd, A_z, A_O, A_y, ErrStat, ErrMsg,&
                                                                                       !(/0.0,REAL(y_axis,ReKi),REAL(z_axis,ReKi)/) )
            temp_filter_velocity(:) = temp_filter_velocity(:) + temp_wind_velocity(:)
                                                                                       
                                                        !+  AD_GetUndisturbedWind ( (REAL(timestep,ReKi)), (/0.0,&             
                                                                                       !REAL(y_axis,ReKi),REAL(z_axis,ReKi)/), ErrStat)
                                                                
                                                        !   AD_GetUndisturbedWind ( (REAL(timestep,ReKi)/20.0)-315.0, (/0.0,&             
                                                                                   !REAL(y_axis,ReKi),REAL(z_axis,ReKi)/), ErrStat)
            number_counter = number_counter + 1
            END IF
          END IF
        END DO
      !END IF
    END DO

    filter_velocity (1) = temp_filter_velocity(1) / number_counter
       
    filter_velocity (2) = temp_filter_velocity(2) / number_counter                 ! Filtered V velocity in the certain radius circle 
    
    filter_velocity (3) = temp_filter_velocity(3) / number_counter                 ! Filtered W velocity in the certain radius circle

END FUNCTION filter_velocity

!---------------------------------------------------------------------------------
SUBROUTINE Get_wake_center ( OS, m, p, y, u, x, xd, z, wakewidth, wake_center )
!................................................................................
! This routine is called to calculate the wake center of a specific release time 
!   and flying time wind plane.
! The wake center is passed to the filter to calculate the averaged wind velocity for
!   the downstream turbine.
!.................................................................................
    !USE DWM_Wake_Deficit_Data, ONLY : m%DWDD%n_x_vector,m%DWDD%ppR 
    !USE MeanderData 
    !USE DWM_ParameterType,     ONLY : p%WakePosition_1, p%WakePosition_2, p%hub_height, p%TurbRefHt, p%Uambient
    !USE BLADE,                 ONLY : R
    !USE AeroDyn_Types
    USE InflowWind                   
    
    TYPE(DWM_OtherStateType),      INTENT(IN   )  :: OS
    TYPE(DWM_MiscVarType),         INTENT(INOUT)  :: m
    TYPE(DWM_ParameterType),       INTENT(IN   )  :: p           ! Parameters
    TYPE(DWM_OutputType),          INTENT(INOUT)  :: y
    TYPE(DWM_InputType),           INTENT(INOUT)  :: u           ! An initial guess for the input; input mesh must be defined
    TYPE(DWM_ContinuousStateType), INTENT(INOUT)  :: x           ! continuous states
    TYPE(DWM_DiscreteStateType),   INTENT(INOUT)  :: xd          ! discrete states
    TYPE(DWM_ConstraintStateType), INTENT(INOUT)  :: z           ! constraint states
    
    REAL(ReKi), ALLOCATABLE,             INTENT(INOUT)  ::   wake_center (:,:,:)  !bjj: this is actually y%wake_position
    INTEGER, ALLOCATABLE,          INTENT(INOUT)  ::   wakewidth(:)
                                                                         
       ! local variables

    REAL(ReKi)                  ::   Modified_U
    INTEGER               ::   release_time
    INTEGER               ::   flying_time
    INTEGER               ::   simulation_time_length
    REAL(DbKi)            ::   DWM_time_step
    REAL(ReKi)                  ::   temp_center_wake (3)
    REAL(ReKi)                  ::   temp_velocity (3)
    REAL(ReKi)                  ::   U_Scale_Factor
    REAL(ReKi)                  ::   U_factor
    REAL(ReKi)                  ::   x_step
    
    INTEGER( IntKi )                   :: ErrStat           ! Error status of the operation
    CHARACTER                          :: ErrMsg            ! Error message if ErrStat /= ErrID_None
    
    real(ReKi)  :: test_1, test_2
    
    !-------------------------------------------------------------
    !!m%DWDD%n_x_vector = 1700
    !!ALLOCATE  (wakewidth(m%DWDD%n_x_vector)) 
    !!wakewidth = 50
    !!m%DWDD%ppR = 50
    !---------------------------------------------------------------
    
    U_factor  =  1.00 

    Modified_U = y%Mean_FFWS * U_factor
    
    
    !------------------------------ TEST ---------------------------
    !m%DWDD%ppR = 50
    !allocate (wakewidth(1750))
    !wakewidth(:) = 60
    !---------------------------------------------------------------
    
    
    DWM_time_step = (2*p%RotorR/m%DWDD%ppR)/Modified_U          ! resolution (126m/50) / wind speed (8m/s) => make sure there is always a wake width at every time step
                                                  ! D/(DWM_time_step*Mean_FFWS)= 50 which is the X resolution

    U_Scale_Factor =  Modified_U / (p%Uambient*U_factor)       ! modify the wake displacement error caused by the change of Mean_FFWS                                                                                        
    
    U_Scale_Factor = 1                       ! 7.15.2015
    
    simulation_time_length = p%WakePosition_1    !80   in reality, 80*scale_factor*DWM_time_step
                                 ! from 1 to 800 (scale_factor : 800/80=10)         to 16D (16*50)
    m%meandering_data%moving_time       = p%WakePosition_2         !50   from 0 to 49  0: wind turbine plane
                                 !               ppR/scale_factor = 5 presents 1D
                                  
                
    release_time      = simulation_time_length       
    flying_time       = m%meandering_data%moving_time       
    m%meandering_data%scale_factor      = 10       ! to decrease the calculation time
    ALLOCATE (wake_center (release_time,flying_time+1,3) )
                                 ! ex. @8D: (1~release_time,8*[ppR/scale_factor]+1,:)

    DO release_time = 1,simulation_time_length,1               ! wake center position at turbine plane
       wake_center (release_time,1,1) = 0
       wake_center (release_time,1,2) = 0
       wake_center (release_time,1,3) = REAL(p%hub_height,ReKi)
    END DO
    
    x_step = Modified_U * (DWM_time_step*m%meandering_data%scale_factor)
    
    IF (.NOT. ALLOCATED(u%IfW%PositionXYZ) ) THEN
       CALL AllocAry( u%IfW%PositionXYZ, 3, 1, "Position array to send to IfW_CalcOutput", ErrStat, ErrMsg )
       IF (ErrStat >= AbortErrLev)  RETURN
    END IF

    
    ! get the initial wake center position of each cross scetion  (from the velocity at the turbine plane * dt)
    DO release_time=1,simulation_time_length,1                        
       wake_center (release_time,2,1) = Modified_U * (DWM_time_step*m%meandering_data%scale_factor) +0             


       !temp_center_wake (:) = AD_WindVelocityWithDisturbance(  (REAL(((release_time-1)+1)*DWM_time_step*m%meandering_data%scale_factor,ReKi)), &
                                             !A_u, A_p, A_x, A_xd, A_z, A_O, A_y, ErrStat, ErrMsg, (/0.0,REAL(0,ReKi),REAL(p%TurbRefHt,ReKi)/) )
                           !AD_GetUndisturbedWind ( (REAL(((release_time-1)+1)*DWM_time_step*m%meandering_data%scale_factor,ReKi)), (/0.0,&                  
                                                !REAL(0,ReKi),REAL(p%TurbRefHt,ReKi)/), ErrStat)                                ! get the velocity at the turbine plane
                                                
       u%IfW%PositionXYZ(1,1) = (0.0_ReKi)
       u%IfW%PositionXYZ(2,1) = (0.0_ReKi)
       u%IfW%PositionXYZ(3,1) = (p%hub_height)

       
       CALL InflowWind_CalcOutput( ( ((release_time-1)+1)*DWM_time_step*m%meandering_data%scale_factor), u%IfW, p%IfW, &
                            x%IfW, xd%IfW, z%IfW, OS%IfW, y%IfW, m%IfW, ErrStat, ErrMsg )

       temp_center_wake (:) = y%IfW%VelocityUVW(:,1)
       !temp_center_wake (3) = y%IfW%Velocity(3,1)
                                                             
       wake_center (release_time,2,2) = temp_center_wake (2) * (DWM_time_step*m%meandering_data%scale_factor) * U_Scale_Factor + wake_center (release_time,1,2)+ &
                                        local_skew_angle(m%NacYaw, m%ct_tilde, wake_center (release_time,2,1), NINT(m%DWDD%ppR), m%DWDD%ppR) * x_step !+ &
                                        !rotation_lateral_offset( wake_center (release_time,2,1) )
                                        
       wake_center (release_time,2,3) = temp_center_wake (3) * (DWM_time_step*m%meandering_data%scale_factor) * U_Scale_Factor + wake_center (release_time,1,3)
    END DO


    DO flying_time = 2,m%meandering_data%moving_time,1
       DO release_time = 1,simulation_time_length,1
       wake_center (release_time,flying_time+1,1) = wake_center (release_time,flying_time+1-1,1) + Modified_U * (DWM_time_step*m%meandering_data%scale_factor)
   
       temp_velocity(:) = filter_velocity (OS,m,p,u,x,xd,z,y,((release_time-1)+1)*DWM_time_step*m%meandering_data%scale_factor, wake_center (release_time,flying_time+1-1,2), &
                                          wake_center (release_time,flying_time+1-1,3), wakewidth((flying_time-1)*m%meandering_data%scale_factor) )

       !!!--------- temp data------
       test_1 = temp_velocity (2) * (DWM_time_step*m%meandering_data%scale_factor) * U_Scale_Factor  + wake_center (release_time,flying_time,2)+ &
                 local_skew_angle(m%NacYaw, m%ct_tilde, wake_center (release_time,flying_time,1), wakewidth((flying_time-1)*m%meandering_data%scale_factor), m%DWDD%ppR) * x_step
       test_2 = temp_velocity (3) * (DWM_time_step*m%meandering_data%scale_factor) * U_Scale_Factor  + wake_center (release_time,flying_time,3)
       !!!-------------------------
       
       wake_center (release_time,flying_time+1,2) = temp_velocity (2) * (DWM_time_step*m%meandering_data%scale_factor) * U_Scale_Factor  + wake_center (release_time,flying_time,2)+ &
                 local_skew_angle(m%NacYaw, m%ct_tilde, wake_center (release_time,flying_time,1), wakewidth((flying_time-1)*m%meandering_data%scale_factor), m%DWDD%ppR) * x_step  !+ &
                 !rotation_lateral_offset( wake_center (release_time,flying_time+1,1) )                            - &
                 !rotation_lateral_offset( wake_center (release_time,flying_time,  1) )
                                                    
       wake_center (release_time,flying_time+1,3) = temp_velocity (3) * (DWM_time_step*m%meandering_data%scale_factor) * U_Scale_Factor  + wake_center (release_time,flying_time,3)
   
   
       END DO
    END DO
    

    
END SUBROUTINE Get_wake_center

!------------------------------------------------------------------------------------------------ 
FUNCTION shifted_velocity( ZTime, p, m, y, z, upwind_mean_u, Uwake, WakeCenter,spacing,angle)
!............................................................................
! This routine is called to get the DWM wake velocity at a certain point in the downstream turbine plane
! Consideirng the meandered wake center
! Uwake(:) is the axial velocity of the wake at the downstream turbine plane
! WakeCenter(:,:,:) is the wake center (y,z) at the downstream turbine plane 
!............................................................................
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !USE    SimCont,                       ONLY: ZTime
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !USE    DWM_Wake_Deficit_Data,         ONLY: m%DWDD%ppR
    !USE    DWM_ParameterType,             ONLY: p%Wind_file_Mean_u,p%hub_height,p%TurbRefHt
    
    TYPE(DWM_MiscVarType),         INTENT(IN   )  :: m
    TYPE(DWM_ParameterType),       INTENT(IN   )  :: p           ! Parameters
    
    REAL(DbKi),                    INTENT(IN   )  :: ZTime
    REAL(ReKi), intent(in)       ::   y,z                           ! point location on the y,z axis
    REAL(ReKi), intent(in)       ::   Uwake(:)                      ! axial velocity of the wake at the downstream turbine plane
    REAL(ReKi), intent(in)       ::   upwind_mean_u                 ! the mean velocity of the turbine UPstream
    REAL(ReKi), intent(in)       ::   WakeCenter(:,:,:)             ! wake_center
    REAL(ReKi), intent(in)       ::   spacing                       ! the distance from the downstream turbine to the upstream turbine
    REAL(ReKi), intent(in)       ::   angle                         ! the angle between the investigated turbine and the line connecting the upwind turbine and wind origin
    
    REAL(ReKi)       ::   shifted_velocity                   ! the output
                                                  !   the velocity at a certain point
    
    REAL(ReKi)       ::   distance                      ! the distance from the point to the meandered wake center
    REAL(ReKi)       ::   y0                            ! wake center position on y axis
    REAL(ReKi)       ::   z0                            ! wake center position on z axis
    REAL(ReKi)       ::   unit                          ! single unit length  R/ppR
    REAL(ReKi)       ::   scale_factor
    INTEGER    ::   p1
    INTEGER    ::   p2
    INTEGER    ::   time_position                 ! to define which plane's wake center is used
    REAL(ReKi)       ::   Yshifted
    REAL(ReKi)       ::   Zshifted
    
    
    !ALLOCATE  (Uwake(NINT( m%DWDD%ppR*Rdomain ))) ! the axis symmetrical velocity
    !ALLOCATE  (WakeCenter( size_of_WakeCenter1,size_of_WakeCenter2,3 ))
    
    scale_factor = 10
    
    time_position = floor(ZTime/( (2*p%RotorR/p%p_p_r/upwind_mean_u/1.00)*scale_factor ))+1  ! ZTime/(DWM_time_step*scale_factor)
    
    
    
    y0 = WakeCenter(time_position,FLOOR(spacing*p%p_p_r/scale_factor)+1,2) !+ 2*P%RotorR*spacing*TAN(m%skew_angle)
    z0 = WakeCenter(time_position,FLOOR(spacing*p%p_p_r/scale_factor)+1,3) !! - REAL(p%TurbRefHt-p%hub_height)
    
    Yshifted = y + 2*p%RotorR*spacing*TAN(angle*Pi/180)
    Zshifted = z
    
    distance = ( (Yshifted-y0)**2 + (Zshifted-z0)**2 )**(0.5)
    unit=p%RotorR / p%p_p_r
    
    p1=FLOOR(distance/unit)
    p2=p1+1
    IF (p1>0) THEN
       shifted_velocity = Uwake(p1)+( Uwake(p2)-Uwake(p1) )*( (distance/unit)-p1 )    ! Weighting method
    ELSE
       shifted_velocity = Uwake(p2)
    END IF

END FUNCTION shifted_velocity

!----------------------------------------------------------------------------------
SUBROUTINE smooth_out_wake(m, p, Uvelocity,Uwake_center,wake_array,spacing,angle,velocity_matrix)
!..................................................................................
! This routine is called to fillter out the smoothed out upstream wake profile
! The output is the wake_array
! Which is the axisymmetrical wake velocity profile
!..................................................................................
    !USE DWM_Wake_Deficit_Data,    ONLY: m%DWDD%n_x_vector, m%DWDD%n_r_vector, m%DWDD%ppR
    !USE DWM_ParameterType,        ONLY: p%hub_height, p%WakePosition_1
    !USE smooth_out_wake_data,     ONLY: m%SmoothOut%length_velocity_array
    
    TYPE(DWM_ParameterType),         INTENT(IN   )   :: p           ! Parameters
    TYPE(DWM_MiscVarType),           INTENT(INOUT)   :: m
    
    REAL(ReKi),ALLOCATABLE   ::   Uvelocity(:,:)
    REAL(ReKi),ALLOCATABLE   ::   Uwake_center(:,:,:)
    REAL(ReKi)               ::   wake_array(:)
    REAL(ReKi)               ::   spacing   ! the spacing between two turbines
    REAL(ReKi)               ::   angle     ! the angle between the downwind turbine and the line conneting the upwind(investigated) turbine and the wind origin
    REAL(ReKi)               ::   velocity_matrix(:,:)      ! the velocity matrix that store the velocity of the downstream turbine plane
    
    INTEGER,ALLOCATABLE  ::   counter_array(:)
    REAL(ReKi),ALLOCATABLE     ::   velocity_array(:)
    !INTEGER              ::   m%SmoothOut%length_velocity_array     ! the length of velocity_array
    INTEGER              ::   low
    INTEGER              ::   high
    INTEGER              ::   i,j,k,n
    INTEGER              ::   counter
    REAL(ReKi)                 ::   y                         ! y coordinate
    REAL(ReKi)                 ::   z                         ! z coordinate
    
    !m%SmoothOut%length_velocity_array = NINT(1.2*R)
    
    !ALLOCATE   ( Uvelocity(m%DWDD%n_x_vector,m%DWDD%n_r_vector) )
    !ALLOCATE   ( Uwake_center(release_time,flying_time+1,3) )
    IF (ALLOCATED( velocity_array ))                DEALLOCATE ( velocity_array )
    IF (ALLOCATED( counter_array ))                 DEALLOCATE ( counter_array )
    ALLOCATE   ( velocity_array (m%SmoothOut%length_velocity_array) )
    !ALLOCATE   ( velocity_matrix (2*m%SmoothOut%length_velocity_array,2*m%SmoothOut%length_velocity_array) )  ! twice the size of the velocity array
    ALLOCATE   ( counter_array (m%SmoothOut%length_velocity_array) )
    !ALLOCATE   ( wake_array (p%ElementNum) )
    
    velocity_array=0
    velocity_matrix = 0
    counter=0
    counter_array=0
    
   !-------------------------------------------------------------------------------------------------
   ! get the time averaged velocity matrix
   !-------------------------------------------------------------------------------------------------
    DO i=1,p%WakePosition_1,1
       DO j=1,2*m%SmoothOut%length_velocity_array,1                        ! y axis
          DO k=1,2*m%SmoothOut%length_velocity_array,1                     ! z axis
             y = (0-m%SmoothOut%length_velocity_array) + (j-1)             ! y coordinate
             z = (p%hub_height-m%SmoothOut%length_velocity_array) + (k-1)    ! z coordinate
             velocity_matrix(j,k) = velocity_matrix(j,k) + smooth_wake_shifted_velocity( m, p, y, z, Uvelocity(floor(spacing * m%DWDD%ppR)+1,:), Uwake_center(:,:,:),spacing,i,angle)
                       !velocity_matrix(j,k) = velocity_matrix(j,k) + smooth_wake_shifted_velocity( y, z, Uvelocity(floor(4.4 * m%DWDD%ppR)+1,:), Uwake_center(:,:,:),4.4,i)
          END DO
       END DO
       counter = counter+1
    END DO
    
    velocity_matrix = velocity_matrix / counter                

   !-------------------------------------------------------------------------------------------------
   ! get the time averaged axisymmetrical velocity array
   !-------------------------------------------------------------------------------------------------
    DO i=1,m%SmoothOut%length_velocity_array,1                          ! velocity array
       DO j=1,2*m%SmoothOut%length_velocity_array,1                       ! velocity_matrix (j,:)
          DO k=1,2*m%SmoothOut%length_velocity_array,1                    ! velocity_matrix (:,k)    
             y = (0-m%SmoothOut%length_velocity_array) + (j-1)            ! y coordinate
             z = (p%hub_height-m%SmoothOut%length_velocity_array) + (k-1)   ! z coordinate
             IF ( ((y-0)**2+(z-p%hub_height)**2)**0.5>(i-1) .and. ((y-0)**2+(z-p%hub_height)**2)**0.5<=i) THEN        
                velocity_array(i) = velocity_array(i)+velocity_matrix(j,k)
                counter_array(i) = counter_array(i)+1
             END IF
          END DO
       END DO
    END DO
    
    DO i=1,m%SmoothOut%length_velocity_array,1                           
       velocity_array(i) = velocity_array(i) / counter_array(i)
    END DO
    
   !-------------------------------------------------------------------------------------------------
   ! get the wake array at the RELM node point
   !-------------------------------------------------------------------------------------------------
    
    DO i=1,p%ElementNum,1
       low  = FLOOR(p%ElementRad(i))
       high = low+1
       wake_array(i) = velocity_array(low) + ( velocity_array(high)-velocity_array(low) )*(p%ElementRad(i)-low)  ! Weighting method
    END DO

END SUBROUTINE smooth_out_wake

!------------------------------------------------------------------------------------------------ 
FUNCTION smooth_wake_shifted_velocity( m, p, y_coor, z_coor, Uwake, WakeCenter,spacing,time_position,angle)
!............................................................................
! This routine is called to get the DWM wake velocity at a certain point in the downstream turbine plane
! Consideirng the meandered wake center
! Uwake(:) is the axial velocity of the wake at the downstream turbine plane
! WakeCenter(:,:,:) is the wake center 
! Used to calculate the smoothed out wake profile
! (y,z) at the downstream turbine plane 
!............................................................................
 
    !USE    DWM_ParameterType,         ONLY: p%p_p_r, p%hub_height, p%TurbRefHt
    !USE    MeanderData,               ONLY: m%meandering_data%scale_factor
    
    TYPE(DWM_ParameterType),         INTENT(IN   )   :: p           ! Parameters
    TYPE(DWM_MiscVarType),           INTENT(INOUT)   :: m

    REAL(ReKi)       ::   y_coor,z_coor                           ! point location on the y,z axis
    REAL(ReKi)       ::   Uwake(:)                      ! axial velocity of the wake at the downstream turbine plane
    REAL(ReKi)       ::   WakeCenter(:,:,:)             ! wake_center
    REAL(ReKi)       ::   spacing                       ! the distance from the downstream turbine to the upstream turbine
    INTEGER    ::   time_position                 ! to define which plane's wake center is used
    REAL(ReKi)       ::   angle                         ! the angle between the downwind turbine and the line conneting the upwind(investigated) turbine and the wind origin
    
    REAL(ReKi)       ::   smooth_wake_shifted_velocity  ! the output
                                                  !   the velocity at a certain point
    REAL(ReKi)       ::   Yshifted
    REAL(ReKi)       ::   Zshifted
    
    INTEGER    ::   p1
    INTEGER    ::   p2
    REAL(ReKi)       ::   distance                      ! the distance from the point to the meandered wake center
    REAL(ReKi)       ::   y0                            ! wake center position on y axis
    REAL(ReKi)       ::   z0                            ! wake center position on z axis
    REAL(ReKi)       ::   unit                          ! single unit length  R/ppR
    
    y0 = WakeCenter(time_position,FLOOR(spacing*p%p_p_r/m%meandering_data%scale_factor)+1,2) !+ 2*P%RotorR*spacing*TAN(m%skew_angle) 
    z0 = WakeCenter(time_position,FLOOR(spacing*p%p_p_r/m%meandering_data%scale_factor)+1,3) !!- REAL(p%TurbRefHt-p%hub_height)
    
    Yshifted = y_coor + 2*p%RotorR*spacing*TAN(angle*Pi/180)
    Zshifted = z_coor
    
    distance = ( (Yshifted-y0)**2 + (Zshifted-z0)**2 )**(0.5)
    unit=p%RotorR/p%p_p_r
    
    p1=FLOOR(distance/unit)
    p2=p1+1
    IF (p1>0) THEN
       smooth_wake_shifted_velocity = Uwake(p1)+( Uwake(p2)-Uwake(p1) )*( (distance/unit)-p1 )   ! Weighting method
    ELSE
       smooth_wake_shifted_velocity = Uwake(p2)
    END IF
    
    
END FUNCTION smooth_wake_shifted_velocity
!----------------------------------------------------------------------------------
FUNCTION TI_downstream_total (m, p, y, spacing,angle,velocity_matrix)   ! name should be calculate_TI_downstream
!..................................................................................
!  This subroutine is called to calculate the TI of the wake deficit
!  The method is by the paper of Rolf-Erik
!  The output is TI_downstream_matrix which is the TI for each computating node in the DWM domain
!..................................................................................
    !USE  DWM_Wake_Deficit_Data,           ONLY: m%DWDD%Turb_Stress_DWM, m%DWDD%n_x_vector, m%DWDD%n_r_vector, m%DWDD%ppR
    !USE  DWM_OutputType,                  ONLY: m%wake_position,m%wake_u
    !USE  MeanderData,                     ONLY: m%meandering_data%moving_time, m%meandering_data%scale_factor
    !USE  DWM_ParameterType,               ONLY: p%TI_amb, p%hub_height,p%TurbRefHt,p%Uambient,p%WakePosition_1
    !USE  smooth_out_wake_data,            ONLY: m%SmoothOut%length_velocity_array
    !USE  BLADE,                           ONLY: R
    
    TYPE(DWM_ParameterType),         INTENT(IN   )   :: p           ! Parameters
    TYPE(DWM_MiscVarType),           INTENT(INOUT)   :: m
    TYPE(DWM_OutputType),            INTENT(INOUT)   :: y
    
    REAL(ReKi)             ::       TI_downstream_total     ! TI of a downstream turbine
    REAL(ReKi)             ::       spacing                 ! the spacing between the downwind turbine and this upwind turbine
    REAL(ReKi)             ::       angle                   ! the angle between the downwind turbine and the line connecting this upwind turbine and the wind direction
    REAL(ReKi)             ::       velocity_matrix(:,:)    ! the velocity matrix at the certain downswind turbine
    
    ! local variables
    REAL(ReKi), ALLOCATABLE  ::  TI_downstream_matrix(:,:)
    INTEGER            ::  i,j,k
    INTEGER            ::  cross_plane_position_ds   ! the cross plane position which to be investigated in term of the flying time
    INTEGER            ::  cross_plane_position_TI   ! the cross plane position which to be investigated in term of the m%DWDD%n_x_vector
    INTEGER            ::  distance_index            ! the index of the distance in the TI axisymmetric array
    INTEGER            ::  counter1
    INTEGER            ::  counter2
    INTEGER            ::  initial_timestep
    REAL(ReKi)               ::  y_axis_turbine
    REAL(ReKi)               ::  z_axis_turbine
    REAL(ReKi)               ::  distance                  ! the distance between one point to the meandered wake center
    REAL(ReKi)               ::  TI_downstream_node        ! the TI at a specfic point in the inbestigated cross plane
    REAL(ReKi)               ::  TI_node_temp
    REAL(ReKi)               ::  TI_node
    REAL(ReKi)               ::  TI_accumulation
    REAL(ReKi)               ::  TI_apprant_accumulation
    REAl(ReKi)               ::  TI_average                ! THE AVERAGE TI OF THE CROSS PLANE
    REAL(ReKi)               ::  TI_apprant                ! The TI due to the meadering
    REAL(ReKi)               ::  HubHt 
    REAL(ReKi)               ::  wake_center_y
    REAL(ReKi)               ::  wake_center_z
    REAL(ReKi)               ::  Rscale
    REAL(ReKi)               ::  y_coor
    REAL(ReKi)               ::  z_coor  
    REAL(ReKi)               ::  zero_spacing
    REAL(ReKi)               ::  temp1,temp2,temp3
    REAL(ReKi)               ::       c_uw
    
   !-------------------------------------------------------------------------------------------------
   ! calculate the TI at each node at the downstream turbine plane from the wake deficit calculation
   !------------------------------------------------------------------------------------------------- 
    IF (ALLOCATED( TI_downstream_matrix ))                DEALLOCATE ( TI_downstream_matrix )
    ALLOCATE (TI_downstream_matrix(m%DWDD%n_x_vector,m%DWDD%n_r_vector))
    c_uw = (0.7550 - m%TI_original/100 *1.75) / 2
    
    DO i=1,m%DWDD%n_x_vector
       DO j=1,m%DWDD%n_r_vector     
          TI_downstream_matrix(i,j) = abs( ( 1/c_uw *m%DWDD%Turb_Stress_DWM(i,j)   ) )**0.5
       END DO
    END DO
    
   !-------------------------------------------------------------------------------------------------
   ! calculate the TI of the downstream turbine grid considering the meandering effect
   !------------------------------------------------------------------------------------------------- 
    cross_plane_position_TI = ANINT( m%DWDD%ppR*spacing+1 )
                     !cross_plane_position_TI = ANINT( m%DWDD%ppR*4.4+1 )
    cross_plane_position_ds = ANINT( (m%DWDD%ppR/m%meandering_data%scale_factor)*spacing+1 )  ! the moving time index of the cross plane in the wake_position(:,:,:)
                     !cross_plane_position_ds = ANINT( (m%DWDD%ppR/m%meandering_data%scale_factor)*4.4+1 )
    
    Rscale = 2  !1.3
    HubHt = p%hub_height
    counter1 = 0
    counter2 = 0
    TI_accumulation = 0
    TI_apprant_accumulation = 0
    
    DO i=1,p%WakePosition_1,1
       DO j=NINT(HubHt-Rscale*p%RotorR),NINT(HubHt+Rscale*p%RotorR),1      ! Z direction
          DO k=1,NINT(2*Rscale*p%RotorR)+1,1                         ! Y direction
             y_axis_turbine = k-(p%RotorR+1) + 2*p%RotorR*spacing*TAN(angle*Pi/180)          ! shift effect
             z_axis_turbine = j
             
             wake_center_y=y%wake_position(i,cross_plane_position_ds,2) !+ 2*P%RotorR*spacing*TAN(m%skew_angle)
             wake_center_z=y%wake_position(i,cross_plane_position_ds,3) !!- REAL(p%TurbRefHt-p%hub_height,ReKi)
             
             distance = ( (y_axis_turbine-wake_center_y)**2 + (z_axis_turbine-wake_center_z)**2)**0.5
             
             distance_index = FLOOR(distance/(p%RotorR/m%DWDD%ppR)) + 1
             
             TI_node_temp = TI_downstream_matrix( cross_plane_position_TI,distance_index )
             
             IF ( TI_node_temp > (y%TI/100*(y%Mean_FFWS/p%Uambient)) ) THEN
                 TI_node = TI_node_temp
             ELSE
                 TI_node = y%TI/100*(y%Mean_FFWS/p%Uambient)
             END IF
             
             TI_accumulation = TI_accumulation + TI_node
             counter1 = counter1+1
             
          END DO
       END DO
    END DO
    
    TI_average = TI_accumulation / REAL(counter1,ReKi) 
   
   !-------------------------------------------------------------------------------------------------
   ! calculate the apprant TI due to the meadering
   !------------------------------------------------------------------------------------------------- 
    zero_spacing = 0
    initial_timestep = 1
    !ALLOCATE (velocity_matrix(2*m%SmoothOut%length_velocity_array,2*m%SmoothOut%length_velocity_array))
    
    DO i=1,2*m%SmoothOut%length_velocity_array,1                                                        ! velocity_matrix (i,:)
      DO j=1,2*m%SmoothOut%length_velocity_array,1                                                        ! velocity_matrix (i,:)  
          y_coor = (0-m%SmoothOut%length_velocity_array) + (i-1)                                              ! y coordinate
          z_coor = (p%hub_height-m%SmoothOut%length_velocity_array) + (j-1)                                     ! z coordinate 
          TI_apprant_accumulation = TI_apprant_accumulation + &
                                    !(smooth_wake_shifted_velocity(y,z,m%wake_u(floor(spacing_turbine * m%DWDD%ppR)+1,:), m%wake_position(:,:,:),spacing_turbine,k) - &
                                    !velocity_matrix(i,j))**2
                                    (velocity_matrix(i,j) - &
                                    smooth_wake_shifted_velocity(m, p, y_coor, z_coor, y%wake_u(floor(spacing * m%DWDD%ppR)+1,:), y%wake_position(:,:,:),zero_spacing,initial_timestep,angle))**2
                                           !smooth_wake_shifted_velocity(y,z,m%wake_u(floor(4.4 * m%DWDD%ppR)+1,:), m%wake_position(:,:,:),zero_spacing,initial_timestep))**2
          counter2=counter2+1
      END DO
    END DO
    
    TI_apprant = ((TI_apprant_accumulation / REAL(counter2,ReKi))**0.5)
  
    
   !-------------------------------------------------------------------------------------------------
   ! calculate the total TI
   !------------------------------------------------------------------------------------------------- 
    !TI_downstream_total = (TI_average**2 + TI_apprant**2)**0.5*100
    
    TI_downstream_total = TI_average * 100
    
    !OPEN(unit = 10, status='replace',file='DWM\results\Downstream_TI_b4normalization.bin',form='unformatted')          
    !WRITE(10)   TI_total                                                                                                                                                                                                                                                                                                         
    !CLOSE(10)
    
    !open (unit=25,file="D:\5MW_simulation\after_release\results\TI.txt")
    !write (25,'(f13.7)'), TI_total
    
    !print*, TI_downstream_matrix(300,120)
    
END FUNCTION TI_downstream_total

!----------------------------------------------------------------------------------
FUNCTION smallscale_TI (m, p, y, spacing,angle,velocity_matrix)   
!..................................................................................
!  This subroutine is called to calculate the smalle scale TI of the wake deficit
!  
!  
!..................................................................................
    !USE  DWM_Wake_Deficit_Data,           ONLY: m%DWDD%Turb_Stress_DWM, m%DWDD%n_x_vector, m%DWDD%n_r_vector, m%DWDD%ppR
    !USE  DWM_OutputType,                  ONLY: m%wake_position,m%wake_u
    !USE  MeanderData,                     ONLY: m%meandering_data%moving_time, m%meandering_data%scale_factor
    !USE  DWM_ParameterType,               ONLY: p%TI_amb, p%hub_height,p%TurbRefHt,p%Uambient,p%WakePosition_1
    !USE  smooth_out_wake_data,            ONLY: m%SmoothOut%length_velocity_array
    !USE  BLADE,                           ONLY: R
    
    TYPE(DWM_ParameterType),         INTENT(IN   )   :: p           ! Parameters
    TYPE(DWM_MiscVarType),           INTENT(INOUT)   :: m
    TYPE(DWM_OutputType),            INTENT(INOUT)   :: y
    
    REAL(ReKi)             ::       spacing                 ! the spacing between the downwind turbine and this upwind turbine
    REAL(ReKi)             ::       angle                   ! the angle between the downwind turbine and the line connecting this upwind turbine and the wind direction
    REAL(ReKi)             ::       velocity_matrix(:,:)    ! the velocity matrix at the certain downswind turbine
    REAL(ReKi)             ::       smallscale_TI
    
    ! local variables
    REAL(ReKi), ALLOCATABLE  ::  TI_downstream_matrix(:,:)
    INTEGER            ::  i,j,k
    INTEGER            ::  cross_plane_position_ds   ! the cross plane position which to be investigated in term of the flying time
    INTEGER            ::  cross_plane_position_TI   ! the cross plane position which to be investigated in term of the m%DWDD%n_x_vector
    INTEGER            ::  distance_index            ! the index of the distance in the TI axisymmetric array
    INTEGER            ::  counter1
    INTEGER            ::  counter2
    INTEGER            ::  initial_timestep
    REAL(ReKi)               ::  y_axis_turbine
    REAL(ReKi)               ::  z_axis_turbine
    REAL(ReKi)               ::  distance                  ! the distance between one point to the meandered wake center
    REAL(ReKi)               ::  TI_downstream_node        ! the TI at a specfic point in the inbestigated cross plane
    REAL(ReKi)               ::  TI_node_temp
    REAL(ReKi)               ::  TI_node
    REAL(ReKi)               ::  TI_accumulation
    REAL(ReKi)               ::  TI_apprant_accumulation
    REAl(ReKi)               ::  TI_average                ! THE AVERAGE TI OF THE CROSS PLANE
    REAL(ReKi)               ::  TI_apprant                ! The TI due to the meadering
    REAL(ReKi)               ::  HubHt 
    REAL(ReKi)               ::  wake_center_y
    REAL(ReKi)               ::  wake_center_z
    REAL(ReKi)               ::  Rscale
    REAL(ReKi)               ::  y_coor
    REAL(ReKi)               ::  z_coor  
    REAL(ReKi)               ::  zero_spacing
    REAL(ReKi)               ::  temp1,temp2,temp3
    
   !-------------------------------------------------------------------------------------------------
   ! calculate the TI at each node at the downstream turbine plane from the wake deficit calculation
   !------------------------------------------------------------------------------------------------- 
    IF (ALLOCATED( TI_downstream_matrix ))                DEALLOCATE ( TI_downstream_matrix )
    ALLOCATE (TI_downstream_matrix(m%DWDD%n_x_vector,m%DWDD%n_r_vector))
    
    DO i=1,m%DWDD%n_x_vector
       DO j=1,m%DWDD%n_r_vector    
          TI_downstream_matrix(i,j) = abs( ( 10/3*m%DWDD%Turb_Stress_DWM(i,j)   ) )**0.5
       END DO
    END DO
    
   !-------------------------------------------------------------------------------------------------
   ! calculate the TI of the downstream turbine grid considering the meandering effect
   !------------------------------------------------------------------------------------------------- 
    cross_plane_position_TI = ANINT( m%DWDD%ppR*spacing+1 )
                     !cross_plane_position_TI = ANINT( m%DWDD%ppR*4.4+1 )
    cross_plane_position_ds = ANINT( (m%DWDD%ppR/m%meandering_data%scale_factor)*spacing+1 )  ! the moving time index of the cross plane in the wake_position(:,:,:)
                     !cross_plane_position_ds = ANINT( (m%DWDD%ppR/m%meandering_data%scale_factor)*4.4+1 )
    
    Rscale = 2 !1.3
    HubHt = p%hub_height
    counter1 = 0
    counter2 = 0
    TI_accumulation = 0
    TI_apprant_accumulation = 0
    
    DO i=1,p%WakePosition_1,1
       DO j=NINT(HubHt-Rscale*p%RotorR),NINT(HubHt+Rscale*p%RotorR),1      ! Z direction
          DO k=1,NINT(2*Rscale*p%RotorR)+1,1                         ! Y direction
             y_axis_turbine = k-(p%RotorR+1) + 2*p%RotorR*spacing*TAN(angle*Pi/180)          ! shift effect
             z_axis_turbine = j
             
             wake_center_y=y%wake_position(i,cross_plane_position_ds,2) !+ 2*P%RotorR*spacing*TAN(m%skew_angle)
             wake_center_z=y%wake_position(i,cross_plane_position_ds,3) !!- REAL(p%TurbRefHt-p%hub_height)
             
             distance = ( (y_axis_turbine-wake_center_y)**2 + (z_axis_turbine-wake_center_z)**2)**0.5
             
             distance_index = FLOOR(distance/(p%RotorR/m%DWDD%ppR)) + 1
             
             TI_node_temp = TI_downstream_matrix( cross_plane_position_TI,distance_index )
             
             IF ( TI_node_temp > (y%TI/100*(y%Mean_FFWS/p%Uambient)) ) THEN
                 TI_node = TI_node_temp
             ELSE
                 TI_node = y%TI/100*(y%Mean_FFWS/p%Uambient)
             END IF
             
             TI_accumulation = TI_accumulation + TI_node
             counter1 = counter1+1
             
          END DO
       END DO
    END DO
    
    TI_average = TI_accumulation / REAL(counter1, ReKi)
    
    smallscale_TI = TI_average * 100
    
END FUNCTION smallscale_TI

!------------------------------------------------------------------------------------------------ 
SUBROUTINE read_parameter_file( p )
!............................................................................
! This routine is called to read the parameter file from the DWM simulation of upstream turbine
! read wake velocity @ the downstream turbine from the upstream wake
! read the meandered wake center
! read the mean velocity of the upstream turbine
!............................................................................
    TYPE(DWM_ParameterType),         INTENT(INOUT)   :: p
    INTEGER            ::       OPENNUM
    INTEGER ErrStat
    
    CALL GetNewUnit(OPENNUM)

    OPEN(unit = OPENNUM, status='old',file='DWM-driver'//trim(PathSep)//'DWM_parameter.bin',form='unformatted',IOSTAT=ErrStat)  ! open an existing file
    IF (ErrStat /= 0) CALL ProgAbort('Error opening existing file, "'//'DWM-driver'//trim(PathSep)//'DWM_parameter.bin"' )
    
    READ(OPENNUM) p%hub_height, p%RotorR, p%NumWT, p%Uambient, &
                  p%TI_amb, p%r_domain, p%x_domain, p%p_p_r, &
                  p%WakePosition_1, p%WakePosition_2, p%WFLowerBd, &
                  p%Winddir 
    CLOSE(OPENNUM) ! close the file
    
    p%Tinfluencer = 1
    
           
    
END SUBROUTINE read_parameter_file

!----------------------------------------------------------------------------------
SUBROUTINE read_turbine_position( m, p, u )
!..................................................................................
! This routine is called at the first of the FAST.
! To decide the position of the turbine in a row
!   if it is the first turbine or not
!..................................................................................
    
    TYPE(DWM_ParameterType),         INTENT(INOUT)   :: p
    TYPE(DWM_MiscVarType),           INTENT(INOUT)   :: m
    TYPE(DWM_InputType),             INTENT(INOUT)   :: u
    
    
    INTEGER(IntKi) :: N_ARGU
    CHARACTER(20)  :: SimulationOrder_index_char    
    !USE read_turbine_position_data
    !USE DWM_ParameterType,           ONLY: p%NumWT,p%Tinfluencer
    

    
    INTEGER   ::    I,J
    INTEGER   ::    MyUn
    
    N_ARGU = COMMAND_ARGUMENT_COUNT()
    IF (N_ARGU < 2) THEN
       CALL ProgAbort('Incorrect number of command arguments in DWM. Arg_1=InputFileName, Arg_2=SimulationOrder_Index, Arg_n="DWM"')
    END IF
    
    
    CALL GET_COMMAND_ARGUMENT(2, SimulationOrder_index_char)
    !print*,SimulationOrder_index_char
    READ(SimulationOrder_index_char,*) p%RTPD%SimulationOrder_index
    !print*,p%RTPD%SimulationOrder_index
    !CALL WrScr(' simulation order index = '//TRIM(Num2LStr(p%RTPD%SimulationOrder_index) )
    
    CALL GetNewUnit(MyUn)
    
    IF ( p%RTPD%SimulationOrder_index /= 0 ) THEN              ! exclude the base turbine
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! obtain the wind turbine index
         !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' 
         
        ALLOCATE (p%RTPD%Turbine_sort_order(p%NumWT))
        OPEN(unit = MyUn, status='old',file='DWM-results'//trim(PathSep)//'wind_farm_turbine_sort.bin',form='unformatted')  ! open an existing file
        READ(MyUn)  p%RTPD%Turbine_sort_order 
        CLOSE(MyUn) ! close the file
    

        p%RTPD%WT_index = p%RTPD%Turbine_sort_order(p%RTPD%SimulationOrder_index)

        
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         !''''''''''''''''''''''''''''''''''''UPWIND  DIRECTION''''''''''''''''''''''''''''''''''''''
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! obtain the index of upwind turbines that affecting this turbine, and the distance/angle
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
        IF (p%RTPD%SimulationOrder_index > 0) THEN
            ALLOCATE (p%RTPD%TurbineInfluenceData(p%NumWT,p%Tinfluencer ) )
        
            OPEN(unit = MyUn, status='old',file='DWM-results'//trim(PathSep)//'turbine_interaction.bin',form='unformatted')  
            READ(MyUn)  p%RTPD%TurbineInfluenceData
            CLOSE(MyUn) 
        
            ! obtain the upwind turbine index
            ALLOCATE (p%RTPD%upwind_turbine_index(p%Tinfluencer))
            DO I = 1,p%Tinfluencer
                p%RTPD%upwind_turbine_index(I) = p%RTPD%TurbineInfluenceData(p%RTPD%WT_index,I)
            END DO
        
            ! calculate the number of upwind turbines affecting the downwind turbine
            p%RTPD%upwindturbine_number = 0
            DO I = 1,p%Tinfluencer
                IF (p%RTPD%upwind_turbine_index(I) /=0) THEN
                p%RTPD%upwindturbine_number = p%RTPD%upwindturbine_number + 1
                END IF
            END DO
        
            ! obtain the upwind turbine coordinates
            ALLOCATE (p%RTPD%wind_farm_Xcoor (p%NumWT))
            ALLOCATE (p%RTPD%wind_farm_Ycoor (p%NumWT))
        
            OPEN(unit = MyUn, status='old',file='DWM-driver'//trim(PathSep)//'wind_farm_coordinate.bin',form='unformatted')  
            READ(MyUn)  p%RTPD%wind_farm_Xcoor,p%RTPD%wind_farm_Ycoor
            CLOSE(MyUn)
        
        
            IF (p%RTPD%upwindturbine_number /= 0) THEN
                ALLOCATE(p%RTPD%upwind_turbine_Xcoor(p%RTPD%upwindturbine_number))
                ALLOCATE(p%RTPD%upwind_turbine_Ycoor(p%RTPD%upwindturbine_number))
                DO I = 1,p%RTPD%upwindturbine_number
                    p%RTPD%upwind_turbine_Xcoor(I) = p%RTPD%wind_farm_Xcoor(p%RTPD%upwind_turbine_index(I))
                    p%RTPD%upwind_turbine_Ycoor(I) = p%RTPD%wind_farm_Ycoor(p%RTPD%upwind_turbine_index(I))
                END DO
            END IF
        
            ! obtain the distance beween the upwind turbine and this turbine
            ALLOCATE (p%RTPD%turbine_windorigin_length (p%NumWT      ))
        
            OPEN(unit = MyUn, status='old',file='DWM-results'//trim(PathSep)//'turbine_distance.bin',form='unformatted')  
            READ(MyUn)  p%RTPD%turbine_windorigin_length
            CLOSE(MyUn)
        
        
            IF (p%RTPD%upwindturbine_number /= 0) THEN
                ALLOCATE (p%RTPD%upwind_turbine_projected_distance(p%RTPD%upwindturbine_number))
                DO I = 1,p%RTPD%upwindturbine_number
                    p%RTPD%upwind_turbine_projected_distance(I) = p%RTPD%turbine_windorigin_length(p%RTPD%WT_index) - p%RTPD%turbine_windorigin_length(p%RTPD%upwind_turbine_index(I))
                END DO
            END IF
        
        
            ! obtain the angle beween the line connecting the upwind turbine and this turbine and the wind direction vector
            ALLOCATE (p%RTPD%turbine_angle(p%NumWT,p%NumWT))
        
            OPEN(unit = MyUn, status='old',file='DWM-results'//trim(PathSep)//'turbine_angles.bin',form='unformatted')  
            READ(MyUn)  p%RTPD%turbine_angle
            CLOSE(MyUn)
        
            IF (p%RTPD%upwindturbine_number /= 0) THEN
                ALLOCATE (p%RTPD%upwind_align_angle(p%RTPD%upwindturbine_number))
                DO I = 1,p%RTPD%upwindturbine_number
                       p%RTPD%upwind_align_angle(I) = p%RTPD%turbine_angle(p%RTPD%upwind_turbine_index(I),p%RTPD%WT_index)
                END DO
            END IF 
        
        END IF
    
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         !''''''''''''''''''''''''''''''''''DOWNWIND  DIRECTION''''''''''''''''''''''''''''''''''''''
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! obtain the index of downwind turbines that being affected by this turbine, and the distance/angle
         !''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''' 
    
        IF (p%RTPD%SimulationOrder_index/=p%NumWT) THEN    ! WHEN index = 0, there will not such files (which are generated by the 0 turbine)
        
           ! obtain the downwind turbine index and calculate the downwind turbine numbers
            p%RTPD%downwindturbine_number = 0
            ALLOCATE (p%RTPD%downwind_turbine_index(p%NumWT-1))
            DO I = 1,p%Tinfluencer
                DO J = 1,p%NumWT
                    IF (p%RTPD%TurbineInfluenceData(J,I) == p%RTPD%WT_index) THEN
                        p%RTPD%downwindturbine_number = p%RTPD%downwindturbine_number + 1
                        p%RTPD%downwind_turbine_index(p%RTPD%downwindturbine_number) = J
                    END IF
                END DO
            END DO
        
            ! obtain the downwind turbine coordinates
            IF (p%RTPD%downwindturbine_number /= 0) THEN
                ALLOCATE(p%RTPD%downwind_turbine_Xcoor(p%RTPD%downwindturbine_number))
                ALLOCATE(p%RTPD%downwind_turbine_Ycoor(p%RTPD%downwindturbine_number))
                DO I = 1,p%RTPD%downwindturbine_number
                    p%RTPD%downwind_turbine_Xcoor(I) = p%RTPD%wind_farm_Xcoor(p%RTPD%downwind_turbine_index(I))
                    p%RTPD%downwind_turbine_Ycoor(I) = p%RTPD%wind_farm_Ycoor(p%RTPD%downwind_turbine_index(I))
                END DO
            END IF
        
            ! obtain the distance beween the upwind turbine and this turbine
            IF (p%RTPD%downwindturbine_number/=0) THEN
                ALLOCATE (p%RTPD%downwind_turbine_projected_distance(p%RTPD%downwindturbine_number))
                DO I = 1,p%RTPD%downwindturbine_number
                   p%RTPD%downwind_turbine_projected_distance(I) = p%RTPD%turbine_windorigin_length(p%RTPD%downwind_turbine_index(I)) - p%RTPD%turbine_windorigin_length(p%RTPD%WT_index)
                END DO
            END IF
        
            ! obtain the angle beween the line connecting the downwind turbine and this turbine and the wind direction vector
            IF (p%RTPD%downwindturbine_number/=0) THEN
                ALLOCATE (p%RTPD%downwind_align_angle(p%RTPD%downwindturbine_number))
                DO I = 1,p%RTPD%downwindturbine_number
                    p%RTPD%downwind_align_angle(I) = p%RTPD%turbine_angle(p%RTPD%WT_index,p%RTPD%downwind_turbine_index(I))
                END DO
            END IF
        
        END IF
    END IF
    
    ! check if the Meandering_Moving_time is valid
    DO I = 1,p%RTPD%downwindturbine_number
        IF (p%WakePosition_2 < p%RTPD%downwind_turbine_projected_distance(I) * p%p_p_r/10 + 1) THEN
         ! bjj: ProgAbort at least is standard fortran (as opposed to "CALL EXIT"), but calling ProgAbort is not allowed in a module in the FAST framework. Please trap your errors and return an error code.
            CALL ProgAbort('WARNING: Meandering_Moving_time should be larger than the maximum turbine spacing')
        END IF
    END DO
              
END SUBROUTINE read_turbine_position

!------------------------------------------------------------------------------------------------ 
SUBROUTINE read_upwind_result_file( m, p, u )
!............................................................................
! This routine is called to read the results from the DWM simulation of upwind turbines
! and to generate the output variables
!............................................................................
    TYPE(DWM_ParameterType),         INTENT(INOUT)   :: p
    TYPE(DWM_MiscVarType),           INTENT(INOUT)   :: m
    TYPE(DWM_InputType),             INTENT(INOUT)   :: u
    
    !USE read_turbine_position_data, ONLY:m%RTPD%upwindturbine_number,m%RTPD%upwind_turbine_index,m%RTPD%WT_index,m%RTPD%downwindturbine_number,m%RTPD%SimulationOrder_index,m%RTPD%Turbine_sort_order
    !USE DWM_ParameterType,          ONLY:p%p_p_r,p%r_domain,p%WakePosition_1,p%WakePosition_2,p%Wind_file_Mean_u
    !USE read_upwind_result
    !USE smooth_out_wake_data,       ONLY:m%SmoothOut%length_velocity_array
    
    CHARACTER(LEN=3)  :: invetigated_turbine_index_character
    CHARACTER(LEN=3)  :: upwind_turbine_index_character
    CHARACTER(LEN=80) :: filename_u_bin,filename_wakecenter_bin,filename_meanU_bin,filename_TI_bin,filename_smoothWake_bin,filename_smallTI_bin
    CHARACTER(LEN=80) :: filename_meanU_txt
    CHARACTER(LEN=2)  :: Uprefix_bin        = 'U_'
    CHARACTER(LEN=3)  :: WCprefix_bin       = 'WC_'
    CHARACTER(LEN=7)  :: MeanUprefix_bin    = 'Mean_U_'
    CHARACTER(LEN=3)  :: Tiprefix_bin       = 'TI_'
    CHARACTER(LEN=8)  :: smallTIprefix_bin  = 'SmallTI_'
    CHARACTER(LEN=11) :: SmoothWprefix_bin  = 'Smoothwake_'
    CHARACTER(LEN=22) :: Prefix             = 'DWM-results'//trim(PathSep)
    CHARACTER(LEN=4)  :: connectionprefix   = '_to_'
    CHARACTER(LEN=8)  :: Turbineprefix      = 'Turbine_'
    INTEGER           :: I
    INTEGER           :: MyUn
    CHARACTER(LEN=3)  :: turbine_sort_order_char
    
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! read the wind file mean velocity at the turbine plane from the very first turbine
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    IF (p%RTPD%SimulationOrder_index > 1) THEN
        IF (p%RTPD%Turbine_sort_order(1) <= 9) THEN
            write(turbine_sort_order_char,'(i1)') p%RTPD%Turbine_sort_order(1)
        ELSEIF (p%RTPD%Turbine_sort_order(1) <= 99) THEN
            write(turbine_sort_order_char,'(i2)') p%RTPD%Turbine_sort_order(1)
        ELSE
            write(turbine_sort_order_char,'(i3)') p%RTPD%Turbine_sort_order(1)
        END IF
        
        filename_meanU_txt = trim(Prefix)//trim(MeanUprefix_bin)//trim(Turbineprefix)//trim(turbine_sort_order_char)//".txt"
        CALL GetNewUnit(MyUn)
        OPEN(unit = MyUn, status='old',file=filename_meanU_txt,form='formatted')
        READ(MyUn,'(f13.7)') p%Wind_file_Mean_u
        CLOSE(MyUn)
    END IF
    
    
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! read the upwind results if have any
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    IF (p%RTPD%upwindturbine_number > 0) THEN
        ALLOCATE (u%Upwind_result%upwind_U          (p%RTPD%upwindturbine_number,floor(p%p_p_r*p%r_domain)))                   ! declare the input
        ALLOCATE (u%Upwind_result%upwind_wakecenter (p%RTPD%upwindturbine_number,p%WakePosition_1,p%WakePosition_2,3))
        ALLOCATE (u%Upwind_result%upwind_meanU      (p%RTPD%upwindturbine_number))
        ALLOCATE (u%Upwind_result%upwind_TI         (p%RTPD%upwindturbine_number))
        ALLOCATE (u%Upwind_result%upwind_small_TI   (p%RTPD%upwindturbine_number))
        ALLOCATE (u%Upwind_result%upwind_smoothWake (p%RTPD%upwindturbine_number,p%ElementNum))
        ALLOCATE (u%Upwind_result%velocity_aerodyn  (p%RTPD%upwindturbine_number))                       ! the temp velocity used by the aerodyn
        
        ! transfer the turbine index from integer to character
        IF (p%RTPD%WT_index <= 9) THEN
            write(invetigated_turbine_index_character,'(i1)') p%RTPD%WT_index
        ELSEIF (p%RTPD%WT_index <= 99) THEN
            write(invetigated_turbine_index_character,'(i2)') p%RTPD%WT_index
        ELSE
            write(invetigated_turbine_index_character,'(i3)') p%RTPD%WT_index
        END IF
        
        DO I = 1,p%RTPD%upwindturbine_number
            
           IF (p%RTPD%upwind_turbine_index(I) <= 9) THEN
              write(upwind_turbine_index_character,'(i1)') p%RTPD%upwind_turbine_index(I)
           ELSEIF (p%RTPD%upwind_turbine_index(I) <= 99) THEN
              write(upwind_turbine_index_character,'(i2)') p%RTPD%upwind_turbine_index(I)
           ELSE
              write(upwind_turbine_index_character,'(i3)') p%RTPD%upwind_turbine_index(I)
           END IF
           
           ! obtain the coresponded profile name
            
           filename_u_bin          = trim(Prefix)//trim(Uprefix_bin)//trim(Turbineprefix)//trim(upwind_turbine_index_character)&
                                     //trim(connectionprefix)//trim(invetigated_turbine_index_character)//".bin"         ! the file name needs to be read

           filename_TI_bin         = trim(Prefix)//trim(TIprefix_bin)//trim(Turbineprefix)//trim(upwind_turbine_index_character)&
                                     //trim(connectionprefix)//trim(invetigated_turbine_index_character)//".bin"
           
           filename_smallTI_bin    = trim(Prefix)//trim(smallTIprefix_bin)//trim(Turbineprefix)//trim(upwind_turbine_index_character)&
                                     //trim(connectionprefix)//trim(invetigated_turbine_index_character)//".bin"
           
           filename_smoothWake_bin = trim(Prefix)//trim(SmoothWprefix_bin)//trim(Turbineprefix)//trim(upwind_turbine_index_character)&
                                     //trim(connectionprefix)//trim(invetigated_turbine_index_character)//".bin"
           
           filename_wakecenter_bin = trim(Prefix)//trim(WCprefix_bin)//trim(Turbineprefix)//trim(upwind_turbine_index_character)//".bin"
           
           filename_meanU_bin      = trim(Prefix)//trim(MeanUprefix_bin)//trim(Turbineprefix)//trim(upwind_turbine_index_character)//".bin"
           
           ! open the file and read
           OPEN(unit = MyUn, status='old',file=filename_u_bin,form='unformatted')  
           READ(MyUn) u%Upwind_result%upwind_U(I,:)          
           CLOSE(MyUn)
           
           OPEN(unit = MyUn, status='old',file=filename_TI_bin,form='unformatted')  
           READ(MyUn) u%Upwind_result%upwind_TI(I) 
           CLOSE(MyUn)
           
           OPEN(unit = MyUn, status='old',file=filename_smallTI_bin,form='unformatted')  
           READ(MyUn) u%Upwind_result%upwind_small_TI(I) 
           CLOSE(MyUn)
           
           OPEN(unit = MyUn, status='old',file=filename_smoothWake_bin,form='unformatted')  
           READ(MyUn) u%Upwind_result%upwind_smoothWake(I,:) 
           CLOSE(MyUn)
           
           OPEN(unit = MyUn, status='old',file=filename_wakecenter_bin,form='unformatted')  
           READ(MyUn) u%Upwind_result%upwind_wakecenter(I,:,:,:) 
           CLOSE(MyUn)
           
           OPEN(unit = MyUn, status='old',file=filename_meanU_bin,form='unformatted')  
           READ(MyUn) u%Upwind_result%upwind_meanU(I) 
           CLOSE(MyUn)
        
        END DO       
    END IF
   
    
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
         ! declare the downwind output variables
         !'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    
    IF (p%RTPD%SimulationOrder_index > 0) THEN                ! not for the base 0 turbine
        IF (p%RTPD%downwindturbine_number /= 0 ) THEN         ! not for the turbines that don't have a downwind turbine
            m%SmoothOut%length_velocity_array = NINT(1.2*p%RotorR)
            ALLOCATE ( u%Upwind_result%TI_downstream             (p%RTPD%downwindturbine_number                                                ) )
            ALLOCATE ( u%Upwind_result%small_scale_TI_downstream (p%RTPD%downwindturbine_number                                                ) )
            ALLOCATE ( u%Upwind_result%smoothed_velocity_array   (p%RTPD%downwindturbine_number,p%ElementNum                                           ) )
            ALLOCATE ( u%Upwind_result%vel_matrix                (p%RTPD%downwindturbine_number,2*m%SmoothOut%length_velocity_array,2*m%SmoothOut%length_velocity_array) ) 
        END IF
    END IF

END SUBROUTINE read_upwind_result_file

!------------------------------------------------------------------------------------------------ 
SUBROUTINE write_result_file(m,p,y,u)
!............................................................................
! This routine is called to write the results from the DWM simulation of this turbine
!............................................................................

    TYPE(DWM_ParameterType),         INTENT(INOUT)   :: p
    TYPE(DWM_MiscVarType),           INTENT(INOUT)   :: m
    TYPE(DWM_OutputType),            INTENT(INOUT)   :: y
    TYPE(DWM_InputType),             INTENT(INOUT)   :: u
    
    CHARACTER(LEN=3)  :: invetigated_turbine_index_character
    CHARACTER(LEN=3)  :: downwind_turbine_index_character
    CHARACTER(LEN=80) :: filename_u_bin,filename_wakecenter_bin,filename_meanU_bin,filename_TI_bin,filename_smallTI_bin,filename_smoothWake_bin,filename_wakewidth_bin,filename_wake_bin
    CHARACTER(LEN=80) :: filename_TI_txt,filename_meanU_txt,filename_induction_txt,filename_wake_txt,filename_wakecenter_txt,filename_SDpower_txt,filename_Ct_txt
    CHARACTER(LEN=2)  :: Uprefix_bin        = 'U_'
    CHARACTER(LEN=3)  :: WCprefix_bin       = 'WC_'
    CHARACTER(LEN=7)  :: MeanUprefix_bin    = 'Mean_U_'
    CHARACTER(LEN=3)  :: Tiprefix_bin       = 'TI_'
    CHARACTER(LEN=8)  :: smallTIprefix_bin  = 'SmallTI_'
    CHARACTER(LEN=11) :: SmoothWprefix_bin  = 'Smoothwake_'
    CHARACTER(LEN=10) :: InductionPrefix    = 'Induction_'
    CHARACTER(LEN=6)  :: Wakeprefix         = 'WakeU_'
    CHARACTER(LEN=11) :: WWprefix_bin       = 'Wake_width_'
    CHARACTER(LEN=22) :: Prefix             = 'DWM-results'//trim(PathSep)
    CHARACTER(LEN=4)  :: connectionprefix   = '_to_'
    CHARACTER(LEN=8)  :: Turbineprefix      = 'Turbine_'
    CHARACTER(LEN=8)  :: Powerprefix        = 'SDpower_'
    CHARACTER(LEN=7)  :: MeanCtPrefix       = 'MeanCt_'
    INTEGER           :: I,RESULT
    
    INTEGER           :: write_unit
    
    CHARACTER(LEN=80) :: filename_TI_to_txt
    CHARACTER(LEN=80) :: filename_SmallTI_to_txt
    
!bjj: use GetNewUnit() and parameters instead of "25" and "10". "10" is especially bad, because some other file is bound to be using it!    
    
    IF ( p%RTPD%SimulationOrder_index > 0 ) THEN            ! exclude the first turbine
        
        IF (p%RTPD%WT_index <= 9) THEN
            write(invetigated_turbine_index_character,'(i1)') p%RTPD%WT_index
        ELSEIF (p%RTPD%WT_index <= 99) THEN
            write(invetigated_turbine_index_character,'(i2)') p%RTPD%WT_index
        ELSE
            write(invetigated_turbine_index_character,'(i3)') p%RTPD%WT_index
        END IF
    
        ! Write the TI of this turbine
        filename_TI_txt    = trim(Prefix)//trim(Tiprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"
          
          CALL  GetNewUnit (write_unit)
          OPEN  (unit = write_unit,file=filename_TI_txt)         
          WRITE (write_unit,'(f13.7)') y%TI
          CLOSE (write_unit) 
    
        ! Write the mean velocity of this turbine
        filename_meanU_bin = trim(Prefix)//trim(MeanUprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".bin"
        filename_meanU_txt = trim(Prefix)//trim(MeanUprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"
          
          CALL  GetNewUnit (write_unit)
          OPEN  (unit = write_unit,file=filename_meanU_txt)         
          WRITE (write_unit,'(f13.7)') y%Mean_FFWS
          CLOSE (write_unit)
          
          CALL  GetNewUnit (write_unit)
          OPEN  (unit = write_unit, status='replace',file=filename_meanU_bin,form='unformatted')    
          WRITE (write_unit)   y%Mean_FFWS                                                                                                                                                                                                                                 
          CLOSE (write_unit)
    
        ! Write the induction factor of this turbine
        filename_induction_txt = trim(Prefix)//trim(InductionPrefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"

          CALL  GetNewUnit (write_unit)    
          OPEN  (unit = write_unit,file=filename_induction_txt)
          WRITE (write_unit,'(f13.7)') y%induction_factor(:)
          CLOSE (write_unit)
          
        ! Write the mean Ct of this turbine
        filename_Ct_txt = trim(Prefix)//trim(MeanCtPrefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"
 
          CALL  GetNewUnit (write_unit)         
          OPEN  (unit = write_unit,file=filename_Ct_txt)
          WRITE (write_unit,'(f13.7)') y%avg_ct
          CLOSE (write_unit)
          
        ! Write the averaged SD power of this turbine
        filename_SDpower_txt = trim(Prefix)//trim(Powerprefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"
        
          CALL  GetNewUnit (write_unit)
          OPEN  (unit = write_unit,file=filename_SDpower_txt,POSITION = 'APPEND')
          WRITE (write_unit,'(f13.7)') y%mean_SDgenpwr
          CLOSE (write_unit)
        
    
        ! Write the wake deficit profile of this turbine
        filename_wake_txt = trim(Prefix)//trim(Wakeprefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"
        filename_wake_bin = trim(Prefix)//trim(Wakeprefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".bin"
    
          CALL  GetNewUnit (write_unit)
          OPEN  (unit = write_unit,file=filename_wake_txt)
          WRITE (write_unit,'(f10.7)') y%wake_u(:,:)
          CLOSE (write_unit)
    
          CALL   GetNewUnit (write_unit)      
          OPEN  (unit = write_unit, status='replace',file=filename_wake_bin,form='unformatted')    
          WRITE (write_unit)   y%wake_u(:,:)                                                                                                                                                                                              
          CLOSE (write_unit)
           
        ! Write the meandered wake center result of this turbine
        filename_wakecenter_bin = trim(Prefix)//trim(WCprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".bin"
        filename_wakecenter_txt = trim(Prefix)//trim(WCprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".txt"
   
          CALL  GetNewUnit (write_unit) 
          OPEN  (unit = write_unit,file=filename_wakecenter_txt)
          WRITE (write_unit,'(f13.7)') y%wake_position(:,:,:)
          CLOSE (write_unit)

          CALL  GetNewUnit (write_unit)    
          OPEN  (unit = write_unit, status='replace',file=filename_wakecenter_bin,form='unformatted')    
          WRITE (write_unit)   y%wake_position(:,:,:)                                                                                                                                                                                              
          CLOSE (write_unit)
    
        ! Write the downstream turbine customized output files
        IF (p%RTPD%downwindturbine_number /= 0) THEN
        
            DO I = 1,p%RTPD%downwindturbine_number
               IF (p%RTPD%downwind_turbine_index(I) <= 9) THEN
                  write(downwind_turbine_index_character,'(i1)') p%RTPD%downwind_turbine_index(I)
               ELSEIF (p%RTPD%downwind_turbine_index(I) <= 99) THEN
                  write(downwind_turbine_index_character,'(i2)') p%RTPD%downwind_turbine_index(I)
               ELSE
                  write(downwind_turbine_index_character,'(i3)') p%RTPD%downwind_turbine_index(I)
               END IF

               filename_u_bin          = trim(Prefix)//trim(Uprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                         //trim(connectionprefix)//trim(downwind_turbine_index_character)//".bin"
               filename_TI_bin         = trim(Prefix)//trim(TIprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                         //trim(connectionprefix)//trim(downwind_turbine_index_character)//".bin"
               filename_smallTI_bin    = trim(Prefix)//trim(smallTIprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                         //trim(connectionprefix)//trim(downwind_turbine_index_character)//".bin" 
               filename_smoothWake_bin = trim(Prefix)//trim(SmoothWprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                         //trim(connectionprefix)//trim(downwind_turbine_index_character)//".bin"
               
               
               filename_TI_to_txt         = trim(Prefix)//trim(TIprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                         //trim(connectionprefix)//trim(downwind_turbine_index_character)//".txt"
               filename_SmallTI_to_txt    = trim(Prefix)//trim(smallTIprefix_bin)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)&
                                         //trim(connectionprefix)//trim(downwind_turbine_index_character)//".txt"
           
               ! Write the wake velocity at the certain downstream turbine plane
                CALL  GetNewUnit (write_unit)
                OPEN  (unit = write_unit, status='replace',file=filename_u_bin,form='unformatted')          
                WRITE (write_unit)   y%wake_u(floor(p%RTPD%downwind_turbine_projected_distance(I) * p%p_p_r)+1,:)               
                CLOSE (write_unit)
            
               ! Write the TI at the certain downstream turbine plane
                CALL  GetNewUnit (write_unit)
                OPEN  (unit = write_unit, status='replace',file=filename_TI_bin,form='unformatted')          
                WRITE (write_unit)   u%Upwind_result%TI_downstream (I)                                                                                                                                                                                                                                                                                                        
                CLOSE (write_unit)
                
                CALL  GetNewUnit (write_unit)
                OPEN  (unit = write_unit, status='replace',file=filename_smallTI_bin,form='unformatted')          
                WRITE (write_unit)   u%Upwind_result%small_scale_TI_downstream (I)                                                                                                                                                                                                                                                                                                        
                CLOSE (write_unit)
        
                CALL  GetNewUnit (write_unit)        
                OPEN  (unit = write_unit,file=filename_TI_to_txt)
                WRITE (write_unit,'(f14.7)') u%Upwind_result%TI_downstream (I)
                CLOSE (write_unit)
                
                CALL  GetNewUnit (write_unit)  
                OPEN  (unit = write_unit,file=filename_smallTI_to_txt)
                WRITE (write_unit,'(f14.7)') u%Upwind_result%small_scale_TI_downstream (I)
                CLOSE (write_unit)
            
               ! Write the smoothed wake profile at the certain downstream turbine plane
                CALL  GetNewUnit (write_unit)
                OPEN  (unit = write_unit, status='replace',file=filename_smoothWake_bin,form='unformatted')          
                WRITE (write_unit)   u%Upwind_result%smoothed_velocity_array(I,:)          ! write the wind data of the plane where the downstream turbine locates,                                                                                                                                                                                                                                                                
                CLOSE (write_unit)
            END DO
        END IF
    END IF
    
        ! Write the meandered wake center and the wake width result from the 0 base wind turbine
    IF (p%RTPD%SimulationOrder_index == 0) THEN
        
            filename_wakecenter_bin = trim(Prefix)//trim(WCprefix_bin)//trim(Turbineprefix)//"0"//".bin"
            filename_wakewidth_bin  = trim(Prefix)//trim(WWprefix_bin)//trim(Turbineprefix)//"0"//".bin"
            
            CALL  GetNewUnit (write_unit)            
            OPEN  (unit = write_unit, status='replace',file=filename_wakecenter_bin,form='unformatted')    
            WRITE (write_unit)   y%wake_position(:,:,:)                                                                                                                                                                                              
            CLOSE (write_unit)
            
            CALL  GetNewUnit (write_unit)
            OPEN  (unit = write_unit, status='replace',file=filename_wakewidth_bin,form='unformatted')               
            WRITE (write_unit)   m%WMC%wake_width(:)                                                                                                                                                                                              
            CLOSE (write_unit)
            
    END IF
       
END SUBROUTINE write_result_file

!!----------------------------------------------------------------------------------
!!BJJ: THIS routine uses non-standard Fortran. IFPORT and DFLIB are incompatible with gfortran. Since this routine isn't used, I'm commenting it out.
!SUBROUTINE rename_FAST_output(m, u, p)
!!............................................................................
!! This routine is called to rename the fast output
!!............................................................................
!    USE IFPORT
!    USE DFLIB
!    
!    TYPE(DWM_InputType),           INTENT(INOUT)  :: u           ! Inputs at Time
!    TYPE(DWM_MiscVarType),         INTENT(INOUT)  :: m
!    TYPE(DWM_ParameterType),       INTENT(INOUT)  :: p
!    
!    CHARACTER(LEN=80) :: filename_FastOutput,filename_FastElm
!    CHARACTER(LEN=11) :: Fastprefix         = 'FastOutput_' ! Fast output file
!    CHARACTER(LEN=8)  :: FastElmprefix      = 'FastElm_'    ! Fast Elm output file
!    INTEGER           :: RESULT
!    CHARACTER(LEN=22) :: Prefix             = 'DWM-results'//trim(PathSep)
!    CHARACTER(LEN=8)  :: Turbineprefix      = 'Turbine_'
!    CHARACTER(LEN=3)  :: invetigated_turbine_index_character
!    
!    IF ( p%RTPD%SimulationOrder_index > 0 ) THEN            ! exclude the first turbine
!        
!        IF (p%RTPD%WT_index <= 9) THEN
!            write(invetigated_turbine_index_character,'(i1)') p%RTPD%WT_index
!        ELSEIF (p%RTPD%WT_index <= 99) THEN
!            write(invetigated_turbine_index_character,'(i2)') p%RTPD%WT_index
!        ELSE
!            write(invetigated_turbine_index_character,'(i3)') p%RTPD%WT_index
!        END IF
!    
!            ! Rename the FAST output wrt the turbine index
!        filename_FastOutput = trim(Prefix)//trim(Fastprefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".out"
!        filename_FastElm    = trim(Prefix)//trim(FastElmprefix)//trim(Turbineprefix)//trim(invetigated_turbine_index_character)//".AD.out"
!    
!        !RESULT =rename('V80_2MW.out',filename_FastOutput)          
!        !RESULT =rename('V80_2MW.elm',filename_FastElm   )
!        
!        RESULT =rename('FAST_V80.out',filename_FastOutput)          
!        RESULT =rename('FAST_V80.elm',filename_FastElm   )
!    END IF
!
!END SUBROUTINE rename_FAST_output

!------------------------------------------------------------------------------------------------ 
FUNCTION min_of_array(ary, ary_length)
!............................................................................
! This routine is called to return the minmum value in an array
!............................................................................
    
    INTEGER ::    ary_length
    REAL(ReKi)    ::    ary(ary_length)
    REAL(ReKi)    ::    min_of_array
    INTEGER ::    I
    

    min_of_array = ary(1) 
    
    DO I = 2,ary_length
        IF (ary(I) < min_of_array) THEN
            min_of_array = ary(I)
        END IF
    END DO

END FUNCTION min_of_array

!------------------------------------------------------------------------------------------------ 
FUNCTION max_of_array(ary, ary_length)
!............................................................................
! This routine is called to return the maximum value in an array
!............................................................................
    
    INTEGER ::    ary_length
    REAL    ::    ary(ary_length)
    REAL    ::    max_of_array
    INTEGER ::    I
    
    max_of_array = ary(1)
    
    DO I = 2,ary_length
        IF (ary(I) > max_of_array) THEN
            max_of_array = ary(I)
        END IF
    END DO

END FUNCTION max_of_array

!------------------------------------------------------------------------------------------------ 
FUNCTION rotation_lateral_offset(x_spacing)
!............................................................................
! This routine is called to return the wake lateral offset due to the turbine rotation
! (assume the rotor spins clockwise) --> always shift to right (negative)
!............................................................................

    REAL      ::       x_spacing
    REAL      ::       rotation_lateral_offset
    
    ! parameters
    REAL      ::       ad = -4.5
    REAL      ::       bd = -0.01
    
    rotation_lateral_offset = ad + bd*x_spacing
 
END FUNCTION rotation_lateral_offset

!------------------------------------------------------------------------------------------------ 
FUNCTION local_skew_angle(yaw_angle, tilde_ct, x_spacing, wake_width, ppr)
!............................................................................
! This routine is called to return the local skew angle at a certain downstream location
!............................................................................
    
    REAL(ReKi)     ::     yaw_angle
    REAL(ReKi)     ::     tilde_ct
    REAL(ReKi)     ::     x_spacing
    INTEGER  ::     wake_width
    REAL(ReKi)     ::     ppr
    REAL(ReKi)     ::     local_skew_angle
    
    IF ( ABS(yaw_angle) > 0.000001 ) THEN
        local_skew_angle = (ppr/wake_width)**2 *COS(yaw_angle)**2 *SIN(yaw_angle) *tilde_ct/2
        local_skew_angle = -local_skew_angle        ! the direction (positive or negative) is opposite to the turbine yaw angle
    ELSE
        local_skew_angle = 0.0
    END IF
    
    local_skew_angle = TAN(local_skew_angle)
 
END FUNCTION local_skew_angle

!------------------------------------------------------------------------------------------------
SUBROUTINE calculate_SD_averagePower( m,y )
!------------------------------------------------------------------------------------------------
! This routine is used to calculate the average time step power from SD subroutine
!------------------------------------------------------------------------------------------------

    TYPE(DWM_MiscVarType),           INTENT(IN   )   :: m
    TYPE(DWM_OutputType),            INTENT(INOUT)   :: y
    
    
    y%mean_SDgenpwr = y%total_SDgenpwr / m%SDtimestep
    
END SUBROUTINE calculate_SD_averagePower
    
END MODULE DWM_Wake_Sub
   
   
   
   
