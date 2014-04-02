!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    This file is part of module IceLoad.
!
!    IceLoad is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with Module2.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!    IceLoad is a module describing ice load on offshore wind turbine supporting structure.  
! 
!    The module is given the module name ModuleName = IceLoad and the abbreviated name ModName = Ice. The mathematical
!    formulation of this module is a subset of the most general form permitted by the FAST modularization framework in tight
!    coupling, thus, the module is developed to support both loose and tight coupling (tight coupling for both time marching and
!    linearization).
!
!
!    References:
!
!    Ice Module Manual, by Bingbin Yu, Dale Karr.
!
!**********************************************************************************************************************************
MODULE IceDyn

   USE IceDyn_Types
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER  :: ID_Ver = ProgDesc( 'IceDyn', 'v1.00.01', '27-Oct-2013' )

   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: ID_Init                           ! Initialization routine
   PUBLIC :: ID_End                            ! Ending routine (includes clean up)

   PUBLIC :: ID_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                 !   continuous states, and updating discrete states
   PUBLIC :: ID_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: ID_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: ID_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: ID_UpdateDiscState                ! Tight coupling routine for updating discrete states

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!..................................................................................................................................

   TYPE(ID_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(ID_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
   TYPE(ID_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
   TYPE(ID_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
   TYPE(ID_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
   TYPE(ID_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
   TYPE(ID_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states
   TYPE(ID_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                               !   only the output mesh is initialized)
   REAL(DbKi),                   INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                               !   (1) ID_UpdateStates() is called in loose coupling &
                                                               !   (2) ID_UpdateDiscState() is called in tight coupling.
                                                               !   Input is the suggested time from the glue code;
                                                               !   Output is the actual coupling interval that will be used
                                                               !   by the glue code.
   TYPE(ID_InitOutputType),      INTENT(  OUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! Local variables

   TYPE(ID_InputFile)                           :: InputFileData           ! Data stored in the module's input file
   INTEGER(IntKi)                               :: ErrStat2                ! temporary Error status of the operation
   CHARACTER(LEN(ErrMsg))                       :: ErrMsg2                 ! temporary Error message if ErrStat /= ErrID_None


      ! Initialize variables for this routine

   ErrStat = ErrID_None
   ErrMsg  = ""


      ! Initialize the NWTC Subroutine Library

   CALL NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information

   CALL DispNVD( ID_Ver )

      !............................................................................................
      ! Read the input file and validate the data
      !............................................................................................
   p%RootName = TRIM(InitInp%RootName)//'_'//TRIM(ID_Ver%Name) ! all of the output file names from this module will end with '_ElastoDyn'

   CALL ID_ReadInput( InitInp, InputFileData, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF ( ErrStat >= AbortErrLev ) RETURN

   ! CALL ID_ValidateInput( InputFileData, ErrStat2, ErrMsg2 )
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF (ErrStat >= AbortErrLev) RETURN

      !............................................................................................
      ! Define parameters here:
      !............................................................................................
   CALL ID_SetParameters( InputFileData, p, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
      
      p%verif = 1
      p%method = 1
      p%dt = Interval
      
   !p%DT  = Interval !bjj: this is set in ID_SetParameters


      !............................................................................................
      ! Define initial system states here:
      !............................................................................................
   
   z%DummyConstrState         = 0                                             ! we don't have constraint states


      ! initialize the continuous states:
   ! CALL Init_ContStates( x, p, InputFileData, OtherState, ErrStat2, ErrMsg2 )     ! initialize the continuous states
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF (ErrStat >= AbortErrLev) RETURN
   
    IF (p%ModNo /= 6) THEN
       
       x%q              = 0                                             ! we don't have continuous states for ice model 1-5 
       x%dqdt           = 0                                             ! we don't have continuous states for ice model 1-5
       
    ELSE
    
       x%q    = p%InitLoc    ! Initial ice floe location 
       x%dqdt = p%v          ! Initial ice velocity
    
    ENDIF
      
      ! initialize the discrete states:
   CALL ID_Init_DiscrtStates( xd, p, InputFileData, ErrStat2, ErrMsg2 )     ! initialize the continuous states
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
   
      ! Initialize other states:
   ! CALL Init_OtherStates( OtherState, p, x, InputFileData, ErrStat2, ErrMsg2 )    ! initialize the other states
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF (ErrStat >= AbortErrLev) RETURN
      if ( p%method .eq. 2) then       

         Allocate( OtherState%xdot(4), STAT=ErrStat )
         IF (ErrStat /= 0) THEN
            ErrStat = ErrID_Fatal
            ErrMsg = ' Error in Module2: could not allocate OtherStat%xdot.'
            RETURN
         END IF

      elseif ( p%method .eq. 3) then       

         Allocate( OtherState%xdot(4), STAT=ErrStat )
         IF (ErrStat /= 0) THEN
            ErrStat = ErrID_Fatal
            ErrMsg = ' Error in Module2: could not allocate OtherStat%xdot.'
            RETURN
         END IF

      endif

      !............................................................................................
      ! Define initial guess for the system inputs here:
      !............................................................................................

         ! allocate all the arrays that store data in the input type:
   ! CALL ID_AllocInput( u, p, ErrStat2, ErrMsg2 )      
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF (ErrStat >= AbortErrLev) RETURN                  
   
    ! Define initial guess for the system inputs here:

      u%q    = 0.
      u%dqdt = 0.

      !............................................................................................
      ! Define system output initializations (set up meshes) here:
      !............................................................................................

   ! CALL ID_AllocOutput(u, y, p,ErrStat2,ErrMsg2) !u is sent so we can create sibling meshes
   !   CALL CheckError( ErrStat2, ErrMsg2 )
   !   IF (ErrStat >= AbortErrLev) RETURN
   
   y%fice = 0.
    
   
      !............................................................................................
      ! Define initialization-routine output here:
      !............................................................................................
   CALL AllocAry( InitOut%WriteOutputHdr, p%NumOuts, 'WriteOutputHdr', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN
   CALL AllocAry( InitOut%WriteOutputUnt, p%NumOuts, 'WriteOutputUnt', ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

   InitOut%WriteOutputHdr = p%OutName
   InitOut%WriteOutputUnt = p%OutUnit

      !............................................................................................
      ! If you want to choose your own rate instead of using what the glue code suggests, tell the glue code the rate at which
      !   this module must be called here:
      !............................................................................................

  ! Interval = p%DT


  !     ! Print the summary file if requested:
  ! IF (InputFileData%SumPrint) THEN
  !    CALL ID_PrintSum( p, OtherState, GetAdamsVals, ErrStat2, ErrMsg2 )
  !       CALL CheckError( ErrStat2, ErrMsg2 )
  !       IF (ErrStat >= AbortErrLev) RETURN
  ! END IF
       
       ! Destroy the InputFileData structure (deallocate arrays)

   CALL ID_DestroyInputFile(InputFileData, ErrStat2, ErrMsg2 )
      CALL CheckError( ErrStat2, ErrMsg2 )
      IF (ErrStat >= AbortErrLev) RETURN

CONTAINS
   !...............................................................................................................................
   SUBROUTINE CheckError(ErrID,Msg)
   ! This subroutine sets the error message and level and cleans up if the error is >= AbortErrLev
   !...............................................................................................................................

         ! Passed arguments
      INTEGER(IntKi), INTENT(IN) :: ErrID       ! The error identifier (ErrStat)
      CHARACTER(*),   INTENT(IN) :: Msg         ! The error message (ErrMsg)

      INTEGER(IntKi)             :: ErrStat3    ! The error identifier (ErrStat)
      CHARACTER(1024)            :: ErrMsg3     ! The error message (ErrMsg)

      !............................................................................................................................
      ! Set error status/message;
      !............................................................................................................................

      IF ( ErrID /= ErrID_None ) THEN

         IF ( LEN_TRIM(ErrMsg) > 0 ) ErrMsg = TRIM(ErrMsg)//NewLine
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(Msg)
         ErrStat = MAX(ErrStat, ErrID)

         !.........................................................................................................................
         ! Clean up if we're going to return on error: close files, deallocate local arrays
         !.........................................................................................................................
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL ID_DestroyInputFile(InputFileData, ErrStat3, ErrMsg3 )
         END IF

      END IF


   END SUBROUTINE CheckError

END SUBROUTINE ID_Init
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!
! This routine is called at the end of the simulation.
!..................................................................................................................................

      TYPE(ID_InputType),           INTENT(INOUT)  :: u           ! System inputs
      TYPE(ID_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
      TYPE(ID_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
      TYPE(ID_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
      TYPE(ID_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
      TYPE(ID_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(ID_OutputType),          INTENT(INOUT)  :: y           ! System outputs
      INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Place any last minute operations or calculations here:

      ! Close files here:

      ! Destroy the input data:

      CALL ID_DestroyInput( u, ErrStat, ErrMsg )

      ! Destroy the parameter data:

      CALL ID_DestroyParam( p, ErrStat, ErrMsg )

      ! Destroy the state data:

      CALL ID_DestroyContState(   x,           ErrStat, ErrMsg )
      CALL ID_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      CALL ID_DestroyConstrState( z,           ErrStat, ErrMsg )
      CALL ID_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )

      ! Destroy the output data:

      CALL ID_DestroyOutput( y, ErrStat, ErrMsg )


END SUBROUTINE ID_End
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input t; Continuous and discrete states are updated for t + p%dt
! (stepsize dt assumed to be in ModName parameter)
!..................................................................................................................................

      REAL(DbKi),                         INTENT(IN   ) :: t          ! Current simulation time in seconds
      INTEGER(IntKi),                     INTENT(IN   ) :: n          ! Current simulation time step n = 0,1,...
      TYPE(ID_InputType),                 INTENT(IN   ) :: u(:)       ! Inputs at utimes
      REAL(DbKi),                         INTENT(IN   ) :: utimes(:)  ! Times associated with u(:), in seconds
      TYPE(ID_ParameterType),             INTENT(IN   ) :: p          ! Parameters
      TYPE(ID_ContinuousStateType),       INTENT(INOUT) :: x          ! Input: Continuous states at t;
                                                                      !   Output: Continuous states at t + Interval
      TYPE(ID_DiscreteStateType),         INTENT(INOUT) :: xd         ! Input: Discrete states at t;
                                                                      !   Output: Discrete states at t  + Interval
      TYPE(ID_ConstraintStateType),       INTENT(INOUT) :: z          ! Input: Initial guess of constraint states at t+dt;
                                                                      !   Output: Constraint states at t+dt
      TYPE(ID_OtherStateType),            INTENT(INOUT) :: OtherState ! Other/optimization states
      INTEGER(IntKi),                     INTENT(  OUT) :: ErrStat    ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(ID_InputType)                :: u_interp                     ! input interpolated from given u at utimes
      TYPE(ID_ContinuousStateType)      :: xdot                         ! continuous state time derivative
      INTEGER(IntKi)                    :: I                            ! Loop count
      REAL(ReKi)                        :: Del2                         ! Deflection of the current ice tooth, for model 2,3
      REAL(ReKi)                        :: Del(p%Zn)                    ! Deflection of ice tooth in each zone, for model 4
      REAL(ReKi)                        :: Tolerence = 1e-6
      REAL(ReKi)                        :: StrRt                        ! Strain rate (s^-1) 
      REAL(ReKi)                        :: SigCrp                       ! Creep stress (Pa)
   
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 
      
      IF (p%ModNo == 2) THEN
      
        Del2 = p%InitLoc + p%v * t - p%Pitch * (xd%IceTthNo2 -1) - u(1)%q
            
            IF ( Del2 >= (p%Delmax + Tolerence) ) THEN 
                
                xd%IceTthNo2 = xd%IceTthNo2 + 1
                !CALL WrScr( 'IceTthNo=' // Num2LStr(xd%IceTthNo2) )
                
            ENDIF
            
      ENDIF
      
      IF  ( p%ModNo == 3 ) THEN
       
           IF (p%SubModNo == 1) THEN
           
               IF (t >= xd%ten) THEN !At the end of the current event, generate random parameters for the next event
                   
                   xd%Nc        = xd%Nc + 1
                   xd%t0n       = xd%ten    ! The beginning time of the next event is the ending time the previous one
                   
                   CALL ID_Generate_RandomStat ( xd, p, ErrStat, ErrMsg)
                   
                   StrRt        = xd%vn / 4 / p%StrWd
                   SigCrp       = p%Cstr * StrRt **(1.0/3.0) 
                   xd%Fmaxn     = p%Ikm * p%StrWd * xd%hn * SigCrp
                   xd%tmn       = xd%t0n + p%Ikm * SigCrp / p%EiPa / StrRt
                   
               ENDIF
               
           ELSEIF (p%SubModNo == 2) THEN
               
               IF (t >= xd%ten) THEN !At the end of the current event, generate random parameters for the next event
                   
                   xd%Nc        = xd%Nc + 1
                   xd%t0n       = xd%ten    ! The beginning time of the next event is the ending time the previous one
                   
                   CALL ID_Generate_RandomStat ( xd, p, ErrStat, ErrMsg)
                   
                   xd%Fmaxn    = xd%hn * p%StrWd * xd%sign
                   xd%tmn      = xd%t0n + xd%sign / p%EiPa / StrRt
                   
                   !CALL WrScr ('t0' // Num2Lstr(xd%t0n) )
                   !CALL WrScr ('tm' // Num2Lstr(xd%tmn) )
                   !CALL WrScr ('te' // Num2Lstr(xd%ten) )
                   
               ENDIF
            
           
           ELSEIF (p%SubModNo == 3) THEN
           
               Del2 = p%InitLoc + p%v * t - xd%Psum - u(1)%q ! Deflection of the current ice tooth
               
               IF ( Del2 >= (p%Delmax + Tolerence) ) THEN 
                
                    xd%Nc       = xd%Nc + 1          ! The current ice tooth breaks
                    xd%Psum     = xd%Psum + xd%Pchn  ! Add the pitch of the current tooth to the sum
                    xd%Dmaxn    = xd%Dmaxnext        ! Copy the properties of the current ice tooth from last "next" ice tooth
                    xd%Pchn     = xd%Pchnext
                    xd%Kn       = xd%Knext
                
                    CALL ID_Generate_RandomStat ( xd, p, ErrStat, ErrMsg) ! Generate random properties of the next ice tooth
                    xd%Knext    = p%h * p%StrWd * xd%signext / xd%Dmaxnext
                
               ENDIF
           
           ENDIF
       
      ENDIF
      
      IF (p%ModNo == 4) THEN
      
        DO I = 1, p%Zn
            
            Del (I) = p%Y0 (I) + p%v * t - p%ZonePitch * (xd%IceTthNo (I)-1) - u(1)%q
            
            IF ( Del(I) >= (p%Delmax + Tolerence) ) THEN 
                
                xd%IceTthNo (I) = xd%IceTthNo (I)+1
                
            ENDIF
            
        END DO
      
      END IF
      
      IF (p%ModNo == 5) THEN
      
      		xd%Beta   = SolveBeta( p%alphaR, p%v, t - xd%Tinit, p%Lbr)
      		
      		IF (xd%Beta >= p%alphaR) THEN
		    	xd%Tinit = t
            END IF
      
      ENDIF
      
      IF (p%ModNo == 6) THEN
      
	        if (p%method .eq. 1) then
	 
	           	CALL ID_RK4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
	
	      	elseif (p%method .eq. 2) then
	
	        	CALL ID_AB4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
	
	      	elseif (p%method .eq. 3) then
	
	         	CALL ID_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
	
	      	else
	
	         	ErrStat = ErrID_Fatal
	         	ErrMsg  = ' Error in Mod2_UpdateStates: p%method must be 1 (RK4), 2 (AB4), or 3 (ABM4)'
	         	RETURN
	         
	      	endif
	      	
	      	! IF ((x%q - u(1)%q) >= xd%dxc) THEN
      		!     xd%dxc = x%q - u(1)%q
      		! ENDIF
      
		      ! Determine whether the splitting failure happens
		      IF (xd%Splitf == 0) THEN
		          CALL Ice_Split (t, u(1), p, x, xd, ErrStat, ErrMsg)
		      ENDIF
	      
      END IF

      IF ( ErrStat >= AbortErrLev ) RETURN
      
      CONTAINS
      
	      FUNCTION SolveBeta (alpha, vice, t, l) Result (beta1)
			
				!SOLVEBETA Solve for ice wedge uprising angle
			
				IMPLICIT NONE
				
				! Input values
				REAL(ReKi) :: alpha 		! Cone angle (rad)
				REAL(ReKi) :: vice 	    ! Ice velocity (m/s)
				REAL(DbKi) :: t 			! Ice thickness
				REAL(ReKi) :: l 			! Ice breaking length
				REAL(ReKi) :: beta = 0 	! Initial Ice wedge uprising angle
	            REAL(ReKi) :: beta1       ! Ice wedge uprising angle
				
				REAL(ReKi) :: Equ		
				REAL(ReKi) :: Derv
				
				DO i = 1,100
				
					Equ = sin(beta) - tan(alpha) * cos(beta) + tan(alpha) * (1-vice*t/l);
			        Derv = cos(beta) + tan(alpha) * sin(beta); 
			        
			        IF ( abs(Equ) <= 1e-6) EXIT	
			        
			        beta = beta - Equ / Derv
			        	
	            END DO
			
	            beta1 = beta
	            
			END FUNCTION SolveBeta
      
END SUBROUTINE ID_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
!
! Routine for computing outputs, used in both loose and tight coupling.
!..................................................................................................................................

      REAL(DbKi),                    INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(ID_InputType),            INTENT(IN   )  :: u           ! Inputs at t
      TYPE(ID_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(ID_ContinuousStateType),  INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(ID_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ID_ConstraintStateType),  INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(ID_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(ID_OutputType),           INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                    !   nectivity information does not have to be recalculated)
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
      INTEGER(IntKi)                    :: I         ! Loop count
      REAL(ReKi)                        :: Del2      ! Deflection of the current ice tooth, for model 2
      REAL(ReKi)                        :: Del(p%Zn) ! Deflection of each ice tooth
      REAL(ReKi)                        :: ZoneF(p%Zn)   ! Ice force of each ice tooth
      REAL(ReKi)                        :: Tolerence = 1e-6
      REAL(ReKi)						:: R
      REAL(ReKi)						:: Pn1
      REAL(ReKi)						:: Pn2
      REAL(ReKi)						:: Xb
      REAL(ReKi)						:: Fb
      REAL(ReKi)						:: pbeta
      REAL(ReKi)						:: qbeta
      REAL(ReKi)						:: Rh
      REAL(ReKi)						:: Rv
      REAL(ReKi)						:: Wr
      REAL(ReKi)						:: gamma
      
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Compute outputs here:
      
      ! Ice Model 1 --------------------------------
      IF (p%ModNo == 1) THEN
          
          IF (p%SubModNo == 1) THEN
              
              IF (t <= p%t0) THEN
                  
                  y%fice = 0
                  
              ELSEIF (t<= p%t0 + p%tm1a) THEN
                  
                  y%fice = (t-p%t0) / p%tm1a * p%Fmax1a
                  
              ELSE
                  
                  y%fice = p%Fmax1a
                  
              ENDIF
              
          ELSEIF (p%SubModNo == 2) THEN
              
              IF (t <= p%t0) THEN
                  
                  y%fice = 0
                  
              ELSEIF (t<= p%t0 + p%tm1b) THEN
                  
                  y%fice = (t-p%t0) / p%tm1b * p%Fmax1b
                  
              ELSE
                  
                  y%fice = 0
                  
              ENDIF
              
          ELSEIF (p%SubModNo == 3) THEN
              
              IF (t <= p%t0) THEN
                  
                  y%fice = 0
                  
              ELSEIF (t<= p%t0 + p%tm1c) THEN
                  
                  y%fice = (t-p%t0) / p%tm1c * p%Fmax1c
                  
              ELSE
                  
                  y%fice = p%Fmax1c
                  
              ENDIF
              
          ENDIF
              
      ENDIF 
      
      ! Ice Model 2 --------------------------------
      
      IF (p%ModNo == 2) THEN
         
          Del2  = p%InitLoc + p%v * t - p%Pitch * (xd%IceTthNo2 -1) - u%q
            
            IF ( Del2 <= 0) THEN 
                
                y%fice = 0.0
                
            ELSEIF (Del2 < p%Delmax2 + Tolerence) THEN
                
                y%fice = Del2 * p%Kice2
                
            ELSE 
                
                ErrStat = ErrID_Fatal
                ErrMsg  = ' Error in IceDyn Model 2: IceToothDel > Delmax'
                
            ENDIF          
       
      END IF
      
      ! Ice Model 3 --------------------------------
      
      IF (p%ModNo == 3) THEN
          
          IF (p%SubModNo == 1) THEN
              
              IF (t <= xd%t0n) THEN
                  
                  y%fice = 0
                   
              ELSEIF (t <= xd%tmn) THEN
                  
                  y%fice = (t-xd%t0n) / (xd%tmn-xd%t0n) * xd%Fmaxn
                  
              ELSE
                  
                  y%fice = xd%Fmaxn
                  
              ENDIF
              
          ELSEIF (p%SubModNo == 2) THEN
              
              IF (t <= xd%t0n) THEN
                  
                  y%fice = 0
                   
              ELSEIF (t <= xd%tmn) THEN
                  
                  y%fice = (t-xd%t0n) / (xd%tmn-xd%t0n) * xd%Fmaxn
                  
              ELSE
                  
                  y%fice = 0
                  
              ENDIF
              
          ELSEIF (p%SubModNo == 3) THEN
              
		        Del2  = p%InitLoc + p%v * t - xd%Psum - u%q     ! Determine the contact state between ice sheet and the tower
              
                IF (Del2 >= xd%Dmaxn) THEN
                    ErrStat = ErrID_Fatal
                    ErrMsg  = ' Error in IceDyn Model 3c: two ice teeth break at once'
                ENDIF
                
		        IF (Del2 <= 0) THEN
                
			       y%fice = 0
               
                ELSE IF (Del2 > 0 .AND. Del2 <= xd%Pchn) THEN
		       
                   y%fice = xd%Kn * Del2
                
                ELSE
                
                   y%fice = xd%Kn * Del2 + xd%Knext * (Del2-xd%Pchn) ! Two teeth in contact
               
                ENDIF	
              
          ENDIF
          
          
          
      ENDIF
      
      ! Ice Model 4 --------------------------------
      
      IF (p%ModNo == 4) THEN
      
        DO I = 1, p%Zn
            
            Del (I) = p%Y0 (I) + p%v * t - p%ZonePitch * (xd%IceTthNo (I)-1) - u%q
            
            IF ( Del (I) <= 0) THEN 
                
                ZoneF (I) = 0.0
                
            ELSEIF (Del (I) < p%Delmax + Tolerence) THEN
                
                ZoneF (I) = Del (I) * p%Kice
                
            ELSE 
                
                ErrStat = ErrID_Fatal
                ErrMsg  = ' Error in IceDyn Model 4: ZoneDel > Delmax'
                
            ENDIF          
            
        END DO
        
        y%fice = sum(ZoneF)
      
      END IF
      
      ! Ice Model 5 --------------------------------
      
      IF (p%ModNo == 5) THEN
      
      	  IF (t <= xd%Tinit) THEN
      	  
      	  	y%fice = 0
      	  
      	  ELSE IF (t == xd%Tinit) THEN  ! Ice load at breakage
      	  
      	  	y%fice = p%RHbr
      	  
      	  ELSE ! Ice load after breakage
      	  
      	  	Wr = p%Wri * ( p%Zr - p%Lbr * sin(xd%Beta) ) / p%Zr * cos(p%alphaR)
            Pn1 = Wr * cos(p%alphaR)
            gamma = p%rhoi / p%rhow
            
            IF (xd%Beta < gamma * p%h / p%Lbr) THEN 
        
                Xb = p%Lbr /3.0 *( ( 3.0*gamma*p%h - p%Lbr*tan(xd%Beta) ) / ( 2.0 *gamma*p%h - p%Lbr*tan(xd%Beta) ) )
                Fb = p%rhow * 9.81 * p%Dwl * (0.5 * p%Lbr**2 * tan(xd%Beta) + p%Lbr*(gamma*p%h - p%Lbr*tan(xd%Beta)) )
        
            ELSE
        
                Xb = 1.0/3.0 * gamma * p%h / sin(xd%Beta)
                Fb = p%rhow * 9.81 * p%Dwl * (0.5 * (gamma*p%h)**2 / tan(xd%Beta) )
        
            END IF
    
            pbeta = sin(xd%Beta) * ( sin(p%alphaR) + p%mu*cos(p%alphaR) ) + cos(xd%Beta) * ( cos(p%alphaR) - p%mu*sin(p%alphaR) )
            qbeta = ( sin(p%alphaR) + p%mu*cos(p%alphaR) ) * sin( p%alphaR - xd%Beta )
    
            Pn2 = ( p%WL * p%Lbr /2.0 * cos(xd%Beta) + Wr * p%Lbr * qbeta - Fb*Xb) / pbeta / p%Lbr
    
            Rh = ( Pn1 + Pn2 ) * ( sin(p%alphaR) + p%mu*cos(p%alphaR) )
            Rv = ( Pn1 + Pn2 ) * ( cos(p%alphaR) - p%mu*sin(p%alphaR) )
            
      	  	y%fice = Rh
      	  
      	  ENDIF
      
      ENDIF
      
      ! Ice Model 6 --------------------------------
      IF (p%ModNo == 6) THEN
      
      	  R = p%StrWd/2
      
	      IF ( xd%Splitf == 0 ) THEN
	          
	        IF ((x%q - u%q) >= xd%dxc .AND. (x%q - u%q) < xd%dxc+R ) THEN
	      
	             y%fice =  p%Cpa * ( 2 * p%h * ( R**2 - (R - x%q + u%q)**2 )**0.5 )**( p%dpa + 1 ) * 1.0e6 
	         
	        ELSE IF (  (x%q - u%q) >= xd%dxc+R ) THEN
	          
	            y%fice = p%Cpa * ( 2 * R *  p%h )**( p%dpa + 1 ) * 1.0e6 
	          
	        ELSE
	          
	            y%fice = 0 
	           
	        ENDIF
	        
	      ELSE 
	          
	           y%fice = 0 
	        
	      ENDIF
      
      ENDIF

END SUBROUTINE ID_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_CalcContStateDeriv( t, u, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )
!
! Routine for computing derivatives of continuous states.
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(ID_InputType),             INTENT(IN   )  :: u           ! Inputs at t
      TYPE(ID_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(ID_ContinuousStateType),   INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(ID_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ID_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(ID_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(ID_ContinuousStateType),   INTENT(  OUT)  :: xdot        ! Continuous state derivatives at t
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
      ! REAL(ReKi) :: mass		! Mass of ice feature (kg)
      REAL(ReKi) :: force		! Ice force (N)
	  REAL(ReKi) :: R			! Structure radius
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 
      
      IF (p%ModNo == 6) THEN

      ! Compute the first time derivatives of the continuous states here:

	  ! When the ice and the structure is in contact, there is ice force.     

	     R = p%StrWd/2 
	     
	     IF ( xd%Splitf == 0 ) THEN
	
	        IF ((x%q - u%q) >= xd%dxc .AND. (x%q - u%q) < xd%dxc+R ) THEN
	
	             force = -p%Cpa * ( 2 * p%h * ( R**2 - (R - x%q + u%q)**2 )**0.5 )**( p%dpa + 1 ) * 1.0e6 + p%FdrN
	      
	        ELSE IF (  (x%q - u%q) >= xd%dxc+R ) THEN
	          
	             force = -p%Cpa * ( 2 * R *  p%h )**( p%dpa + 1 ) * 1.0e6  + p%FdrN
	         
	        ELSE
	
	            force = 0. + p%FdrN
	
	        ENDIF
	        
	      ELSE
	          
	          force = 0
	          
	      ENDIF
	      
	      xdot%q = x%dqdt
	
	      xdot%dqdt = force / p%Mice
      
      ENDIF

END SUBROUTINE ID_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for updating discrete states
!..................................................................................................................................

      REAL(DbKi),                    INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                INTENT(IN   )  :: n           ! Current step of the simulation: t = n*Interval
      TYPE(ID_InputType),            INTENT(IN   )  :: u           ! Inputs at t
      TYPE(ID_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(ID_ContinuousStateType),  INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(ID_DiscreteStateType),    INTENT(INOUT)  :: xd          ! Input: Discrete states at t;
                                                                    !   Output: Discrete states at t + Interval
      TYPE(ID_ConstraintStateType),  INTENT(IN   )  :: z           ! Constraint states at t
      TYPE(ID_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                  INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      ! TYPE(Ice_InputType)             :: u_interp                     ! input interpolated from given u at utimes
      ! TYPE(Ice_ContinuousStateType)   :: xdot                         ! continuous state time derivative
      INTEGER(IntKi)                    :: I                            ! Loop count
      REAL(ReKi)                        :: Del2                         ! Deflection of the current ice tooth, for model 2,3
      REAL(ReKi)                        :: Del(p%Zn)                    ! Deflection of ice tooth in each zone, for model 4
      REAL(ReKi)                        :: Tolerence = 1e-6
      REAL(ReKi)                        :: StrRt                        ! Strain rate (s^-1) 
      REAL(ReKi)                        :: SigCrp                       ! Creep stress (Pa)
   
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 
      
      ! Update discrete states here:
      
      IF (p%ModNo == 2) THEN
      
        Del2 = p%InitLoc + p%v * t - p%Pitch * (xd%IceTthNo2 -1) - u%q
            
            IF ( Del2 >= (p%Delmax + Tolerence) ) THEN 
                
                xd%IceTthNo2 = xd%IceTthNo2 + 1
                !CALL WrScr( 'IceTthNo=' // Num2LStr(xd%IceTthNo2) )
                
            ENDIF
            
      ENDIF
      
      IF  ( p%ModNo == 3 ) THEN
       
           IF (p%SubModNo == 1) THEN
           
               IF (t >= xd%ten) THEN !At the end of the current event, generate random parameters for the next event
                   
                   xd%Nc        = xd%Nc + 1
                   xd%t0n       = xd%ten    ! The beginning time of the next event is the ending time the previous one
                   
                   CALL ID_Generate_RandomStat ( xd, p, ErrStat, ErrMsg)
                   
                   StrRt        = xd%vn / 4 / p%StrWd
                   SigCrp       = p%Cstr * StrRt **(1.0/3.0) 
                   xd%Fmaxn     = p%Ikm * p%StrWd * xd%hn * SigCrp
                   xd%tmn       = xd%t0n + p%Ikm * SigCrp / p%EiPa / StrRt
                   
               ENDIF
               
           ELSEIF (p%SubModNo == 2) THEN
               
               IF (t >= xd%ten) THEN !At the end of the current event, generate random parameters for the next event
                   
                   xd%Nc        = xd%Nc + 1
                   xd%t0n       = xd%ten    ! The beginning time of the next event is the ending time the previous one
                   
                   CALL ID_Generate_RandomStat ( xd, p, ErrStat, ErrMsg)
                   
                   xd%Fmaxn    = xd%hn * p%StrWd * xd%sign
                   xd%tmn      = xd%t0n + xd%sign / p%EiPa / StrRt
                   
                   CALL WrScr ('t0' // Num2Lstr(xd%t0n) )
                   CALL WrScr ('tm' // Num2Lstr(xd%tmn) )
                   CALL WrScr ('te' // Num2Lstr(xd%ten) )
                   
               ENDIF
            
           
           ELSEIF (p%SubModNo == 3) THEN
           
               Del2 = p%InitLoc + p%v * t - xd%Psum - u%q ! Deflection of the current ice tooth
               
               IF ( Del2 >= (p%Delmax + Tolerence) ) THEN 
                
                    xd%Nc       = xd%Nc + 1          ! The current ice tooth breaks
                    xd%Psum     = xd%Psum + xd%Pchn  ! Add the pitch of the current tooth to the sum
                    xd%Dmaxn    = xd%Dmaxnext        ! Copy the properties of the current ice tooth from last "next" ice tooth
                    xd%Pchn     = xd%Pchnext
                    xd%Kn       = xd%Knext
                
                    CALL ID_Generate_RandomStat ( xd, p, ErrStat, ErrMsg) ! Generate random properties of the next ice tooth
                    xd%Knext    = p%h * p%StrWd * xd%signext / xd%Dmaxnext
                
               ENDIF
           
           ENDIF
       
      ENDIF
      
      IF (p%ModNo == 4) THEN
      
        DO I = 1, p%Zn
            
            Del (I) = p%Y0 (I) + p%v * t - p%ZonePitch * (xd%IceTthNo (I)-1) - u%q
            
            IF ( Del(I) >= (p%Delmax + Tolerence) ) THEN 
                
                xd%IceTthNo (I) = xd%IceTthNo (I)+1
                
            ENDIF
            
        END DO
      
      END IF
	
	  IF (p%ModNo == 5) THEN
      
      		xd%Beta   = SolveBeta( p%alphaR, p%v, t - xd%Tinit, p%Lbr)
      		IF (xd%Beta >= p%alphaR) THEN
		    	xd%Tinit = t
            END IF
      
      ENDIF
	
      IF (p%ModNo == 6) THEN
      
      	! IF ((x%q - u(1)%q) >= xd%dxc) THEN
      	!     xd%dxc = x%q - u(1)%q
      	! ENDIF
      
	      ! Determine whether the splitting failure happens
	      IF (xd%Splitf == 0) THEN
	          CALL Ice_Split (t, u, p, x, xd, ErrStat, ErrMsg)
	      ENDIF
      
      ENDIF
      
      IF ( ErrStat >= AbortErrLev ) RETURN
      
      
      CONTAINS
      
	      FUNCTION SolveBeta (alpha, vice, t, l) Result (beta1)
			
				!SOLVEBETA Solve for ice wedge uprising angle
			
				IMPLICIT NONE
				
				! Input values
				REAL(ReKi) :: alpha 		! Cone angle (rad)
				REAL(ReKi) :: vice 	    ! Ice velocity (m/s)
				REAL(DbKi) :: t 			! Ice thickness
				REAL(ReKi) :: l 			! Ice breaking length
				REAL(ReKi) :: beta = 0 	! Initial Ice wedge uprising angle
	            REAL(ReKi) :: beta1       ! Ice wedge uprising angle
				
				REAL(ReKi) :: Equ		
				REAL(ReKi) :: Derv
				
				DO i = 1,100
				
					Equ = sin(beta) - tan(alpha) * cos(beta) + tan(alpha) * (1-vice*t/l);
			        Derv = cos(beta) + tan(alpha) * sin(beta); 
			        
			        IF ( abs(Equ) <= 1e-6) EXIT	
			        
			        beta = beta - Equ / Derv
			        	
	            END DO
			
	            beta1 = beta
	            
		   END FUNCTION SolveBeta
      

!      xd%DummyDiscState = 0.0

END SUBROUTINE ID_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )
!
! Routine for solving for the residual of the constraint state functions
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(ID_InputType),            INTENT(IN   )  :: u           ! Inputs at t
      TYPE(ID_ParameterType),        INTENT(IN   )  :: p           ! Parameters
      TYPE(ID_ContinuousStateType),  INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(ID_DiscreteStateType),    INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ID_ConstraintStateType),  INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(ID_OtherStateType),       INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(ID_ConstraintStateType),  INTENT(  OUT)  :: Z_residual  ! Residual of the constraint state functions using
                                                                    !     the input values described above
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 


      ! Solve for the constraint states here:

      Z_residual%DummyConstrState = 0

END SUBROUTINE ID_CalcConstrStateResidual
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_ReadInput( InitInp, InputFileData, ErrStat, ErrMsg )
!     This public subroutine reads the input required for IceDyn from the file whose name is an  
!     input parameter.
!----------------------------------------------------------------------------------------------------   

   
      ! Passed variables
   
   TYPE(ID_InitInputType),        INTENT( IN    )   :: InitInp              ! the IceDyn initial input data 
   TYPE(ID_InputFile),            INTENT(   OUT )   :: InputFileData        ! Data stored in the IceDyn's input file
   INTEGER,                       INTENT(   OUT )   :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                  INTENT(   OUT )   :: ErrMsg               ! Error message if ErrStat /= ErrID_None

  
      ! Local variables  
         
   INTEGER                                          :: UnIn                 ! Unit number for the input file
   CHARACTER(1024)                                  :: FileName             ! Name of HydroDyn input file  
   
   
   
      ! Initialize ErrStat
         
   ErrStat = ErrID_None         
   ErrMsg  = ""   
   
   
   !-------------------------------------------------------------------------------------------------
   ! Open the file
   !-------------------------------------------------------------------------------------------------   
   FileName = TRIM(InitInp%InputFile)
   
   CALL GetNewUnit( UnIn )   
   CALL OpenFInpFile( UnIn, FileName, ErrStat )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to open IceDyn input file: '//FileName
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   
   !CALL WrScr( 'Opening HydroDyn input file:  '//FileName )
   
   
   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------
   
   CALL ReadCom( UnIn, FileName, 'IceDyn input file header line 1', ErrStat )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read HydroDyn input file header line 1.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF


   CALL ReadCom( UnIn, FileName, 'IceDyn input file header line 2', ErrStat )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read HydroDyn input file header line 2.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   CALL ReadCom( UnIn, FileName, 'IceDyn input file header line 3', ErrStat )
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read HydroDyn input file header line 3.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Ice Model Number
   !-------------------------------------------------------------------------------------------------

      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Ice Models header', ErrStat, ErrMsg)
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Comment line - Ice Models'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF


      ! IceModel - Number that represents different ice models.
      
   CALL ReadVar ( UnIn, FileName, InputFileData%IceModel, 'IceModel', 'Number that represents different ice models', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read IceModel parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

    
    ! IceSubModel - Number that represents different ice sub models.
      
   CALL ReadVar ( UnIn, FileName, InputFileData%IceSubModel, 'IceSubModel', 'Number that represents different ice sub-models', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read IceSubModel parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   !-------------------------------------------------------------------------------------------------
   ! Ice properties - General
   !-------------------------------------------------------------------------------------------------
      
      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Ice properties - General header', ErrStat, ErrMsg)
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice properties - General header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

     ! IceVel - Velocity of ice sheet movement (m/s)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%v, 'IceVel', 'Velocity of ice sheet movement', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read IceVel parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! IceThks  - Thickness of the ice sheet (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%h, 'IceThks', 'Thickness of the ice sheet (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read IceThks parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! StWidth  - Width of the structure in contact with the ice (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%StrWd, 'StWidth', 'Width of the structure in contact with the ice (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read StWidth parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! WtDen     - Mass density of water (kg/m3)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%rhow, 'WtDen', 'Mass density of water (kg/m3)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WtDen parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

      ! IceDen  - Mass density of ice (kg/m3)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%rhoi, 'IceDen', 'Mass density of ice (kg/m3)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read IceDen parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

          ! InitLoc - Ice sheet initial location (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%InitLoc, 'InitLoc', 'Ice sheet initial location (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read InitLoc parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

          ! InitTm - Ice load starting time (s)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%t0, 'InitTm', 'Ice load starting time (s)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read InitTm parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Ice properties - Ice Model 1
   !-------------------------------------------------------------------------------------------------
   
        ! Header
      
   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 1 SubModel 1', ErrStat, ErrMsg)
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice properties - Ice Model 1 SubModel 1 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
    ! Ikm - Indentation factor
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Ikm, 'Ikm', 'Indentation factor', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ikm parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
    ! Ag - Constant depends only on ice crystal type, used in calculating uniaxial stress (MPa-3s-1)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Ag, 'Ag', 'Constant depends only on ice crystal type, used in calculating uniaxial stress (MPa-3s-1)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ag parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
    ! Qg - Activation Energy (kJ•mol^-1)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Qg, 'Qg', 'Activation Energy (kJ•mol^-1)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Qg parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
    ! Rg - Universal gas constant (J•mol-1K-1)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Rg, 'Rg', 'Universal gas constant (J•mol-1K-1)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Rg parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
    ! Tice - Ice temperature (K)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Tice, 'Tice', 'Ice temperature (K)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Tice parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   
    ! Header
      
   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 1 SubModel 2', ErrStat, ErrMsg)
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice properties - Ice Model 1 SubModel 2 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
    ! Poisson  - Poisson's ratio of ice
      
   CALL ReadVar ( UnIn, FileName, InputFileData%nu, 'Poisson', 'Poisson ratio of ice', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Poisson parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
    ! WgAngle - Wedge Angel, degree. Default 45 Degrees.
      
   CALL ReadVar ( UnIn, FileName, InputFileData%phi, 'WgAngle', ' Wedge Angel, degree. Default 45 Degrees', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read WgAngle parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
    ! EIce - Young's modulus of ice (GPa)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Eice, 'EIce', 'Youngs modulus of ice (GPa)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read EIce parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
    ! Header
      
   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 1 SubModel 3', ErrStat, ErrMsg)
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice properties - Ice Model 1 SubModel 3 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
    ! SigN - Nominal ice stress (MPa)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%SigN, 'SigN', 'Nominal ice stress (MPa)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read SigN parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Ice properties - Ice Model 2
   !-------------------------------------------------------------------------------------------------
   
   ! Header
      
   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 2 SubModel 1,2', ErrStat, ErrMsg)
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice properties - Ice Model 2 SubModel 1,2 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
    ! Pitch - Distance between sequential ice teeth (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Pitch, 'Pitch', 'Distance between sequential ice teeth (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Pitch parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
    ! IceStr2 - Ice failure stress (MPa)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%IceStr2, 'IceStr2', 'Ice failure stress (MPa)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read IceStr2 parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
    ! Delmax2 - Ice tooth maximum elastic deformation (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Delmax2, 'Delmax2', 'Ice tooth maximum elastic deformation (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Delmax2 parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Ice properties - Ice Model 3
   !-------------------------------------------------------------------------------------------------
   
     ! Header
      
   CALL ReadCom( UnIn, FileName, 'Ice PROPERTIES -Ice Model 3, SubModel 1,2', ErrStat, ErrMsg)
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice PROPERTIES -Ice Model 3, SubModel 1,2 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! ThkMean - Mean value of ice thickness (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%miuh, 'ThkMean', 'Mean value of ice thickness (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read ThkMean parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! ThkVar - Variance of ice thickness (m^2)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%varh, 'ThkVar', 'Variance of ice thickness (m^2)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read ThkVar parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! VelMean - Mean value of ice velocity (m/s)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%miuv, 'VelMean', 'Mean value of ice velocity (m/s)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read VelMean parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! VelVar - Variance of ice velocity (m^2/s^2)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%varv, 'VelVar', 'Variance of ice velocity (m^2/s^2)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read VelVar parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! TeMean - Mean value of ice loading event duration (s)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%miut, 'TeMean', 'Mean value of ice loading event duration (s)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read TeMean parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
    ! Header
      
   CALL ReadCom( UnIn, FileName, 'Ice PROPERTIES -Ice Model 3, SubModel 2,3', ErrStat, ErrMsg)
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice PROPERTIES -Ice Model 3, SubModel 2,3 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! StrMean - Mean value of ice thickness (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%miubr, 'StrMean', 'Mean value of ice strength (MPa)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read StrMean parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! StrVar - Variance of ice thickness (m^2)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%varbr, 'StrVar', 'Variance of ice strength (MPa)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read StrVar parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! Header
      
   CALL ReadCom( UnIn, FileName, 'Ice PROPERTIES -Ice Model 3, SubModel 3', ErrStat, ErrMsg)
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice PROPERTIES -Ice Model 3, SubModel 3 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! DelMean - Mean value of maximum ice tooth tip displacement (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%miuDelm, 'DelMean', 'Mean value of maximum ice tooth tip displacement (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read DelMean parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! DelVar - Variance of maximum ice tooth tip displacement (m^2)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%varDelm, 'DelVar', 'Variance of maximum ice tooth tip displacement (m^2)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read DelVar parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! PMean - Mean value of the distance between sequential ice teeth (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%miuP, 'PMean', 'Mean value of the distance between sequential ice teeth (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PMean parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! PVar - Variance of the distance between sequential ice teeth (m^2)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%varP, 'PVar', 'Variance of the distance between sequential ice teeth (m^2)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PVar parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Ice properties - Ice Model 4
   !-------------------------------------------------------------------------------------------------
      
      ! Header
      
   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 4', ErrStat, ErrMsg)
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice properties - Ice Model 4 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

     ! PrflMean - Mean value of ice contact face position (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%PrflMean, 'PrflMean', 'Mean value of ice contact face position (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PrflMean parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! PrflSig - Standard deviation of ice contact face position (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%PrflSig, 'PrflSig', 'Standard deviation of ice contact face position (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read PrflSig parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! ZoneNo1 - Number of failure zones along contact width 
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Zn1, 'ZoneNo1', 'Number of failure zones along contact width', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read ZoneNo1 parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
      
     ! ZoneNo2 - Number of failure zones along contact height/thickness
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Zn2, 'ZoneNo2', 'Number of failure zones along contact height/thickness', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read ZoneNo2 parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
      ! ZonePitch - Distance between sequential ice teeth (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%ZonePitch, 'ZonePitch', 'Distance between sequential ice teeth (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read ZonePitch parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
      ! IceStr - Ice failure stress within each failure region (MPa)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%IceStr, 'IceStr', 'Ice failure stress within each failure region (MPa)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read IceStr parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF

   ! Delmax - Ice teeth maximum elastic deformation (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Delmax, 'Delmax', 'Ice teeth maximum elastic deformation (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Delmax parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Ice properties - Ice Model 5
   !-------------------------------------------------------------------------------------------------
   
   ! Header
      
   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 5, Submodel 1,2', ErrStat, ErrMsg)
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice properties - Ice Model 5 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
     ! ConeAgl - Slope angle of the cone (degree)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%alpha, 'ConeAgl', 'Slope angle of the cone (degree)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read ConeAgl parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
     ! ConeDwl - Cone waterline diameter (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Dwl, 'ConeDwl', 'Cone waterline diameter (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read ConeDwl parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
     ! ConeDtp - Cone top diameter (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Dtp, 'ConeDtp', 'Cone top diameter (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read ConeDtp parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
     ! RdupThk - Ride-up ice thickness (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%hr, 'RdupThk', 'Ride-up ice thickness (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read RdupThk parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
     ! mu - Friction coefficient between structure and ice (-)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%mu, 'mu', 'Friction coefficient between structure and ice (-)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read mu parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   
     ! FlxStr - Flexural strength of ice (MPa)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%sigf, 'FlxStr', 'Flexural strength of ice (MPa)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read FlxStr parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! StrLim - Limit strain for ice fracture failure (-)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%StrLim, 'StrLim', 'Limit strain for ice fracture failure (-)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read StrLim parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   ! StrRtLim - Limit strain rate for ice brittle behavior (s^-1)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%StrRtLim, 'StrRtLim', 'Limit strain rate for ice brittle behavior (s^-1)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read StrRtLim parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
   !-------------------------------------------------------------------------------------------------
   ! Ice properties - Ice Model 6
   !-------------------------------------------------------------------------------------------------
   
   ! Header
      
   CALL ReadCom( UnIn, FileName, 'Ice properties - Ice Model 6', ErrStat, ErrMsg)
   
   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Ice properties - Ice Model 6 header comment line.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
   
     ! FloeLth - Ice floe length (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Ll, 'FloeLth', 'Ice floe length (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read FloeLth parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
 
   ! FloeWth - Ice floe width (m)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Lw, 'FloeWth', 'Ice floe width (m)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read FloeWth parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
 
   ! CPrAr - Ice crushing strength pressure-area relation constant
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Cpa, 'CPrAr', 'Ice crushing strength pressure-area relation constant', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read CPrAr parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
 
   ! dPrAr - Ice crushing strength pressure-area relation order
      
   CALL ReadVar ( UnIn, FileName, InputFileData%dpa, 'dPrAr', 'Ice crushing strength pressure-area relation order', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read dPrAr parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
 
   ! Fdr - Constant external driving force (MN)
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Fdr, 'Fdr', 'Constant external driving force (MN)', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Fdr parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
 
   ! Kic - Fracture toughness of ice (kNm^(-3/2))
      
   CALL ReadVar ( UnIn, FileName, InputFileData%Kic, 'Kic', 'Fracture toughness of ice (kNm^(-3/2))', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read Kic parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
 
	! FspN - Non-dimensional splitting load
	      
   CALL ReadVar ( UnIn, FileName, InputFileData%FspN, 'FspN', 'Non-dimensional splitting load', ErrStat, ErrMsg )

   IF ( ErrStat /= ErrID_None ) THEN
      ErrMsg  = ' Failed to read FspN parameter.'
      ErrStat = ErrID_Fatal
      CLOSE( UnIn )
      RETURN
   END IF
 
   
   !-------------------------------------------------------------------------------------------------
   ! This is the end of the input file
   !-------------------------------------------------------------------------------------------------
   CLOSE ( UnIn )
   
   
   RETURN    
   
   
END SUBROUTINE ID_ReadInput
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_SetParameters( InputFileData, p, ErrStat, ErrMsg  )
! This takes the primary input file data and sets the corresponding parameters.
!..................................................................................................................................

   IMPLICIT                        NONE


      ! Passed variables

   TYPE(ID_ParameterType),   INTENT(INOUT)  :: p                            ! Parameters of the IceDyn module
   TYPE(ID_InputFile),       INTENT(IN)     :: InputFileData                ! Data stored in the module's input file
   INTEGER(IntKi),           INTENT(OUT)    :: ErrStat                      ! Error status
   CHARACTER(*),             INTENT(OUT)    :: ErrMsg                       ! Error message

     ! Local variables
   
     ! Ice Model 1
     
   REAL(ReKi)                               :: StrRt                        ! Strain rate (s^-1) 
   REAL(ReKi)                               :: SigCrp                       ! Creep stress (Pa) 
   REAL(ReKi)                               :: Bf                           ! Flexural rigidity of ice plate (Nm)
   REAL(ReKi)                               :: kappa                        ! Constants in Ice Model 1b
   REAL(ReKi)                               :: PhiR                         ! Phi in radius
   
     ! Ice Model 4
   REAL(ReKi)                               :: ZoneWd                       ! Width of a single failure zone
   REAL(ReKi)                               :: ZoneHt                       ! Height of a single failuer zone
   INTEGER(IntKi)                           :: I
   
     ! Ice Model 5
   REAL(ReKi)								:: flexStrPa					! Ice flexural strength (Pa)  
   REAL(ReKi)								:: A(6)							! Coefficients when calculating ice breaking force using sub-model 1
   REAL(ReKi)								:: Pn1
   REAL(ReKi)								:: Pn2
   REAL(ReKi)								:: F1
   REAL(ReKi)								:: Lxlim1
   REAL(ReKi)								:: Lxlim2
   REAL(ReKi)                               :: Pbr
   
   
!bjj: ERROR CHECKING!!!

      ! Initialize error data
   ErrStat = ErrID_None
   ErrMsg  = ''

   !...............................................................................................................................
   ! Direct copy of variables:
   !...............................................................................................................................
   p%ModNo     = InputFileData%IceModel
   p%SubModNo  = InputFileData%IceSubModel
   p%v         = InputFileData%v
   p%h         = InputFileData%h 
   p%InitLoc   = InputFileData%InitLoc
   p%t0        = InputFileData%t0
   p%StrWd     = InputFileData%StrWd

   ! Ice Model 1
   p%Ikm       = InputFileData%Ikm
   
   ! Ice Model 2
   p%Delmax2   = InputFileData%Delmax2
   p%Pitch     = InputFileData%Pitch
   
   ! Ice Model 3
   p%miuh      = InputFileData%miuh
   p%varh      = InputFileData%varh
   p%miuv      = InputFileData%miuv
   p%varv      = InputFileData%varv
   p%miut      = InputFileData%miut
   p%miubr     = InputFileData%miubr
   p%varbr     = InputFileData%varbr
   p%miuDelm   = InputFileData%miuDelm
   p%varDelm   = InputFileData%varDelm
   p%miuP      = InputFileData%miuP
   p%varP      = InputFileData%varP
   
   ! Ice Model 4
   p%Delmax    = InputFileData%Delmax
   p%ZonePitch = InputFileData%ZonePitch
   
   ! Ice Model 5
   p%rhoi	   = InputFileData%rhoi
   p%rhow	   = InputFileData%rhow
   p%Dwl	   = InputFileData%Dwl
   p%mu		   = InputFileData%mu
      
   ! Ice Model 6
   p%Cpa	   = InputFileData%Cpa
   p%dpa	   = InputFileData%dpa
   
   !...............................................................................................................................
   ! Calculate some indirect inputs:
   !...............................................................................................................................
   
   StrRt       = p%v / 4 / p%StrWd
   p%EiPa      = InputFileData%EIce * 1.0e9
   
   ! Ice Model 1a
   
   SigCrp      = ( 1/InputFileData%Ag * exp( InputFileData%Qg / InputFileData%Rg / InputFileData%Tice ) * StrRt )**(1.0/3.0) * 1e6
   p%Cstr      = ( 1/InputFileData%Ag * exp( InputFileData%Qg / InputFileData%Rg / InputFileData%Tice ) )**(1.0/3.0) * 1e6
   p%Fmax1a    = InputFileData%Ikm * p%StrWd * p%h * SigCrp
   p%tm1a      = InputFileData%Ikm * SigCrp / p%EiPa / StrRt
   
   ! Ice Model 1b
   
   Bf          = p%EiPa * p%h**3 / 12.0 / (1.0-InputFileData%nu**2)
   kappa       = ( InputFileData%rhow * 9.81 / 4.0 / Bf ) ** 0.25
   PhiR        = InputFileData%phi / 180.0 * 3.1415927
   p%Fmax1b    = 5.3 * Bf * kappa * ( kappa * p%StrWd + 2 * tan(PhiR/2.0) )
   p%tm1b      = p%Fmax1b / p%StrWd / p%h / p%EiPa / StrRt
   
   ! Ice Model 1c
   
   p%Fmax1c     = p%StrWd * p%h * InputFileData%SigN *1e6
   p%tm1c       = InputFileData%SigN / p%EiPa / StrRt
   
   ! Ice Model 2
   
   p%Kice2      = InputFileData%IceStr2 *1e6 * p%StrWd * p%h / p%Delmax2
   
   ! Ice Model 4
   ZoneWd      = InputFileData%StrWd / REAL(InputFileData%Zn1)
   ZoneHt      = InputFileData%h     / REAL(InputFileData%Zn2)
   p%Kice      = InputFileData%IceStr *1e6 * ZoneWd * ZoneHt / InputFileData%Delmax
   
   p%Zn        = InputFileData%Zn1   * InputFileData%Zn2                           ! Total number of failure zones
   
   CALL AllocAry( p%ContPrfl, p%Zn, 'ContPrfl', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   
    CALL AllocAry( p%Y0, p%Zn, 'ContPrfl', ErrStat, ErrMsg )
   IF ( ErrStat >= AbortErrLev ) RETURN
   
   p%Y0        = InputFileData%InitLoc + p%ContPrfl - MAXVAL(p%ContPrfl)
   
   ! Ice Model 5
   flexStrPa   = InputFileData%sigf * 1e6
   
   p%alphaR	   = InputFileData%alpha / 180.0 * 3.1415927
   p%Zr		   = (InputFileData%Dwl - InputFileData%Dtp) / 2 * tan(p%alphaR)
   p%Wri 	   = p%rhoi * 9.81 * p%Dwl * p%h * p%Zr / sin(p%alphaR)
   
   IF (p%SubModNo == 1) THEN
   		
   		p%LovR = SolveLambda ( p%rhoi, p%h, p%Dwl, flexStrPa )
   		p%Lbr  = p%LovR * p%Dwl / 2
   		A	   = BrkLdPar (p%alphaR, p%LovR, InputFileData%mu)
   		
   		p%RHbr = ( A(1) * flexStrPa * p%h**2 + A(2) * p%rhoi * 9.81 * p%h * p%Dwl**2 + A(3) * p%rhoi * 9.81 * p%h * (p%Dwl**2 - InputFileData%Dtp**2) ) * A(4)
		p%RVbr = A(5) * p%RHbr + A(6) * p%rhoi * 9.81 * p%h * (p%Dwl**2 - InputFileData%Dtp**2)
   
        CALL WrScr(Num2LStr(A(1)))
        CALL WrScr(Num2LStr(A(2)))
        CALL WrScr(Num2LStr(A(3)))
        CALL WrScr(Num2LStr(A(4)))
        CALL WrScr(Num2LStr(A(5)))
        CALL WrScr(Num2LStr(A(6)))
        !CALL WrScr(Num2LStr(p%RVbr))
        
   ELSEIF (p%SubModNo == 2) THEN
   
   		Pbr	   = 8.0 * sqrt(2.0) * ( ( flexStrPa * p%h**2) / 4 )
   		Pn1	   = p%Wri * cos(p%alphaR)
   		F1     = p%Wri * ( sin(p%alphaR) + p%mu * cos(p%alphaR) );
        Pn2    = ( Pbr + F1*sin(p%alphaR) ) / ( cos(p%alphaR) - p%mu * sin(p%alphaR) )
   		
   		p%RHbr = (Pn1 + Pn2) * ( sin(p%alphaR) + p%mu * cos(p%alphaR) )
   		p%RVbr = (Pn1 + Pn2) * ( cos(p%alphaR) - p%mu * sin(p%alphaR) )
   		
   		Lxlim1 = ( 3.0 * sqrt(6.0) ) / 8.0 * ( p%v * tan(p%alphaR) ) / InputFileData%StrRtLim   !Limit strain rate criteria
   		Lxlim2 = sqrt(6.0) *  ( ( flexStrPa * p%h**2) / 4 / (p%rhow * 9.81) / InputFileData%StrLim )**(1.0/3.0)
   		
   		IF (Lxlim1 <= Lxlim2) THEN
   		
   			p%Lbr = ( 3.0 * sqrt(2.0) ) / 8.0 * ( p%v * tan(p%alphaR) ) / InputFileData%StrRtLim
   		
   		ELSE
   		
   			p%Lbr = 2.0 *  ( ( flexStrPa * p%h**2) / 4 / (p%rhow * 9.81) / InputFileData%StrLim )**(1.0/3.0)
   		
   		ENDIF
   
   ELSE
   		
   		ErrMsg	= 'Sub-model number for model 5 should be 1 or 2'
   		ErrStat = ErrID_Fatal
   
   ENDIF 
   
   p%WL		   = p%rhoi * 9.81 * p%Dwl * p%h * p%Lbr
   
   ! Ice Model 6
   p%FdrN	   = InputFileData%Fdr * 1e6
   p%Mice	   = InputFileData%rhoi * p%h * InputFileData%Lw * InputFileData%Ll
   p%Fsp 	   = InputFileData%FspN * p%h * InputFileData%Kic * 1e3 * sqrt(InputFileData%Ll) 
   
   !...............................................................................................................................
   ! Calculate Output variables:
   !...............................................................................................................................

   p%NumOuts    = 3
      
      CALL AllocAry( p%OutName, p%NumOuts, 'OutName', ErrStat, ErrMsg )
      IF (ErrStat >= AbortErrLev) RETURN
      
      p%OutName (1) = 'Time'
      p%OutName (2) = 'IceDisp'
      p%OutName (3) = 'IceForce'
      
      CALL AllocAry( p%OutUnit, p%NumOuts, 'OutUnit', ErrStat, ErrMsg )
      IF (ErrStat >= AbortErrLev) RETURN
      
      p%OutUnit (1) = '(s)'
      p%OutUnit (2) = '(m)'
      p%OutUnit (3) = '(kN)'
      
   CONTAINS 

!Functions that generate random number with respect to certain distributions

        FUNCTION random_normal() RESULT(fn_val)
		
			! Adapted from the following Fortran 77 code
			!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
			!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
			!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
			
			!  The function random_normal() returns a normally distributed pseudo-random
			!  number with zero mean and unit variance.
			
			!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
			!  and J.F. Monahan augmented with quadratic bounding curves.
			
			REAL(ReKi) :: fn_val
			
			!     Local variables
			REAL(ReKi)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
			            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
			
			!     Generate P = (u,v) uniform in rectangle enclosing acceptance region
			
			DO
			  CALL RANDOM_NUMBER(u)
			  CALL RANDOM_NUMBER(v)
			  v = 1.7156 * (v - 0.5)
			
			!     Evaluate the quadratic form
			  x = u - s
			  y = ABS(v) - t
			  q = x**2 + y*(a*y - b*x)
			
			!     Accept P if inside inner ellipse
			  IF (q < r1) EXIT
			!     Reject P if outside outer ellipse
			  IF (q > r2) CYCLE
			!     Reject P if outside acceptance region
			  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
			END DO
			
			!     Return ratio of P's coordinates as the normal deviate
			fn_val = v/u
			RETURN
			
        END FUNCTION random_normal
        
        	 
        FUNCTION  SolveLambda(rhoi, t, D, sigf) Result (rho)
		
		!	SOLVERHO Solve for rho according to Ralston model (Ralston 1978)
		!   Rho = A/R, A is the first circumferential crack radius, R is the cone
		!   structure waterline radius. According to first equation on Ralston 
		!   paper (P301), calculate rho.
			
			IMPLICIT NONE
		
			! Input values
			REAL(ReKi) :: rhoi		! Mass density of ice, (kg/m^3)
			REAL(ReKi) :: t 			! Ice thickness (m)
			REAL(ReKi) :: D 			! Cone waterline diameter (m)
			REAL(ReKi) :: sigf 		! Ice flextural strength (Pa)
			
			REAL(ReKi) :: rho			! Rho = A/R
			REAL(ReKi) :: x = 1.01    ! Initial value of rho
			REAL(ReKi) :: Equ		
			REAL(ReKi) :: Derv
			
			DO i = 1,100
			
				Equ = x - log(x) + 0.0922 * rhoi * 9.81 * t * D**2 / sigf / t**2 * (2*x+1) * (x-1)**2 - 1.369
				Derv = 1 - 1/x + 0.0922 * rhoi * 9.81 * t * D**2 / sigf / t**2 * (2*(x-1)**2 + (2*x+1)*2*(x-1))
				
				IF ( abs(Equ) <= 1e-6) THEN
				
					rho = x
					EXIT
				
				END IF 
				
				x = x - Equ / Derv
			
			END DO
		          
		END FUNCTION SolveLambda
		
		
		FUNCTION BrkLdPar (alpha, lambda, mu) Result (A)
		
		!BRKLDPAR Calculates Ralston's horizontal force paramters A1, A2, A3, A4 and B1, B2.
		!   Detailed explanation in Ralston's paper: Ice Force Desgin Consideration 
		!	for Conical Offshore Structure and Ice Module Manual
		
			IMPLICIT NONE
		
			! Input values
			REAL(ReKi) :: alpha			! Cone angle, (rad)
			REAL(ReKi) :: lambda 		! Ratio of breaking length over cone waterline radius
			REAL(ReKi) :: mu 			! Friction coefficient between structure and ice
			
			REAL(ReKi) :: A(6)   		! Coefficients when calculating ice breaking force
			
			! Local variables
			REAL(ReKi) :: f
			REAL(ReKi) :: g
			REAL(ReKi) :: h
			REAL(ReKi) :: pi = 3.1415927
			
			A(1) = 1.0/3.0 * ( lambda/(lambda-1) + (1-lambda+lambda*log(lambda))/(lambda-1) + 2.422* (lambda*log(lambda))/(lambda-1) )

			A(2) = ( lambda**2 + lambda -2.0 )/12.0

			f = pi/2.0 + pi/8.0 * (sin(alpha))**2 / (1-(sin(alpha))**2) - pi/16.0 * (sin(alpha))**4 / (1-(sin(alpha))**4)
			g = ( 1.0/2.0 + alpha/sin(2*alpha) ) / ( pi/4.0*sin(alpha) + mu*alpha*cos(alpha)/sin(alpha) )
			A(3) = 1.0/4.0 * ( 1/cos(alpha) + mu*Esina(alpha,5)/sin(alpha) - mu*f*g/tan(alpha) )

			A(4) = tan(alpha) / ( 1 - mu * g)

			h = cos(alpha) - mu/sin(alpha) * ( Esina(alpha,5) - cos(alpha)**2 * Fsina(alpha) )

			A(5) = h / ( pi/4.0 * sin(alpha) + mu * alpha / tan(alpha) )
			A(6) = 1.0/4.0 * (pi/2.0*cos(alpha) - mu*alpha - f*h/ ( pi/4.0 * sin(alpha) + mu * alpha / tan(alpha) ))
            
            !CALL WrScr ('fac='// Num2LStr(factorial(5)))
            !CALL WrScr(Num2LStr(Esina(alpha,10)))
            !CALL WrScr(Num2LStr(Fsina(alpha)))
					
		END FUNCTION BrkLdPar
		
		
		FUNCTION Esina (alpha, n) Result (Esin)
		!ESINA calculates E(sin(alpha)). Detailed explanation in Ice Module Manual, Model 5, Submodel 1
		
			IMPLICIT NONE
			
			!Input variable
			REAL(ReKi)     :: alpha			! Cone angle, (rad)
			INTEGER(IntKi) :: n		
			
			!Output
			REAL(ReKi)     :: Esin 
            
			
			!Local variable 
			INTEGER(IntKi) :: i	
			REAL(ReKi)     :: E = 0
            REAL(ReKi)     :: pi = 3.1415927
			
			DO i = 1,n
			
				E = E + pi/2.0*( factorial(2*(i-1)) / 2**(2*(i-1)) / (factorial(i-1))**2 )**2 * (sin(alpha))**(2*(i-1)) / (1-2*(i-1))
			
            END DO
            
            Esin = E
		    
		END FUNCTION Esina
		
		
		FUNCTION Fsina (alpha) Result (F)
		!ESINA calculates F(sin(alpha)). Detailed explanation in Ice Module Manual, Model 5, Submodel 1
		
			IMPLICIT NONE
			
			!Input variable
			REAL(ReKi)     :: alpha			! Cone angle, (rad)
			!Output
			REAL(ReKi)     :: F
			!Local variable 
			REAL(ReKi)     :: pi = 3.1415927
			
			F = pi/2.0 + pi/8.0 * sin(alpha)**2 / (1-sin(alpha)**2) - pi/16.0 * sin(alpha)**4 / (1-sin(alpha)**4)
		
		END FUNCTION Fsina
		

        FUNCTION factorial (n) Result (fac)
        ! FACTORIAL calculates the factorial of n 
        
        	IMPLICIT NONE
        	
        	!Input variable
        	INTEGER(IntKi),INTENT(IN) :: n
        	
        	!Output
        	REAL(ReKi)	   :: fac
        	
        	!Local variables
        	INTEGER(IntKi) :: i	
        	REAL(ReKi) :: M 
            
            M = 1
        	
        	DO i = 1,n
        	
        		M = M * i  
        	
            ENDDO
        	
            !CALL WrScr ('new')
            !CALL WrScr ('n='// Num2LStr(n))
            !CALL WrScr ('M='// Num2LStr(M))
        	fac = REAL(M)
            !CALL WrScr ('fac='// Num2LStr(fac))
        
        END FUNCTION factorial
        
        
   END SUBROUTINE ID_SetParameters
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_Init_DiscrtStates( xd, p, InputFileData, ErrStat, ErrMsg  )
! This routine initializes the continuous states of the module.
! It assumes the parameters are set and that InputFileData contains initial conditions for the continuous states.
!..................................................................................................................................
   IMPLICIT                        NONE
   
   TYPE(ID_DiscreteStateType),   INTENT(OUT)    :: xd                ! Initial continuous states
   TYPE(ID_ParameterType),       INTENT(IN)     :: p                 ! Parameters of the structural dynamics module
   TYPE(ID_InputFile),           INTENT(IN)     :: InputFileData     ! Data stored in the module's input file
   INTEGER(IntKi),               INTENT(OUT)    :: ErrStat           ! Error status
   CHARACTER(*),                 INTENT(OUT)    :: ErrMsg            ! Error message

      ! local variables
   INTEGER(IntKi)                               :: I                 ! loop counter
   REAL(ReKi)                                   :: StrRt             ! Strain rate (s^-1) 
   REAL(ReKi)                                   :: SigCrp            ! Creep stress (Pa) 
      ! Initialize error data
   ErrStat = ErrID_None
   ErrMsg  = ''

   IF  ( p%ModNo == 2 ) THEN
       
       xd%IceTthNo2 = 1 ! Initialize first ice tooth number
       
   ENDIF 
   
   IF  ( p%ModNo == 3 ) THEN
       
       xd%Nc = 1 ! Initialize first loading event/ice tooth number
       xd%t0n = p%t0
       CALL ID_Generate_RandomStat ( xd, p, ErrStat, ErrMsg)
       
       StrRt       = xd%vn / 4 / p%StrWd
       
       IF (p%SubModNo == 1) THEN
           
           SigCrp      = p%Cstr * StrRt **(1.0/3.0) 
           xd%Fmaxn    = p%Ikm * p%StrWd * xd%hn * SigCrp
           xd%tmn      = xd%t0n + p%Ikm * SigCrp / p%EiPa / StrRt
           
       ELSEIF (p%SubModNo == 2) THEN
           
           xd%Fmaxn    = xd%hn * p%StrWd * xd%sign
           xd%tmn      = xd%t0n + xd%sign / p%EiPa / StrRt
           
       ELSEIF (p%SubModNo == 3) THEN
           
           xd%Kn       = p%h * p%StrWd * xd%sign / xd%Dmaxn
           xd%Knext    = p%h * p%StrWd * xd%signext / xd%Dmaxnext
           xd%Psum     = 0
           
       ENDIF
       
   ENDIF
   
   IF ( p%ModNo == 4 ) THEN
   
      ! First allocate the arrays stored here:

       CALL AllocAry( xd%IceTthNo, p%Zn,   'IceTthNo',   ErrStat, ErrMsg )
       IF ( ErrStat /= ErrID_None ) RETURN ! Initialize first ice tooth number for each zone
       DO I = 1,p%Zn
          xd%IceTthNo (I) = 1
       END DO
   
   END IF
   
   IF ( p%ModNo == 5 ) THEN
   
      xd%beta = 0.    ! ice crushed depth
      xd%Tinit = p%t0
            
   END IF
   
   IF ( p%ModNo == 6 ) THEN
   
      xd%dxc = 0.    ! ice crushed depth
      xd%Splitf = 0. ! flag to indicate if the ice floe has splitted (0 not splitted, 1 splitted)
      
   END IF

END SUBROUTINE ID_Init_DiscrtStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_Generate_RandomStat ( xd, p, ErrStat, ErrMsg)
! This routine initializes the continuous states of the module.
! It assumes the parameters are set and that InputFileData contains initial conditions for the continuous states.
!..................................................................................................................................
   IMPLICIT                        NONE   

   TYPE(ID_DiscreteStateType),   INTENT(INOUT)  :: xd                ! Initial continuous states
   TYPE(ID_ParameterType),       INTENT(IN)     :: p                 ! Parameters of the structural dynamics module
   INTEGER(IntKi),               INTENT(OUT)    :: ErrStat           ! Error status
   CHARACTER(*),                 INTENT(OUT)    :: ErrMsg            ! Error message

      ! local variables
   INTEGER(IntKi)                               :: I                 ! loop counter
   REAL(ReKi)                                   :: SigLogh           ! sigma_log(h), standard deviation of log(h)
   REAL(ReKi)                                   :: MiuLogh           ! miu_log(h), mean value of log(h)
   REAL(ReKi)                                   :: VelSig            ! parameter for a Rayleigh distribution
   REAL(ReKi)                                   :: TeLamb            ! parameter for a exponential distribution
   
      ! Initialize error data
   ErrStat = ErrID_None
   ErrMsg  = ''
   
      !Ice thickness has a lognormal distribution
   SigLogh = SQRT( LOG ( p%varh / p%miuh + 1) )
   MiuLogh = LOG ( p%miuh ) - 0.5 * SigLogh **2 
   xd%hn   = EXP( random_normal() * SigLogh + MiuLogh )
   
      !Ice velocity has a Rayleigh distribution
   VelSig = p%miuv / SQRT(Pi/2)
   xd%vn  = random_rayleigh (VelSig)
   
      !Iceloading event time has a exponential distribution
   TeLamb = 1 / p%miut
   xd%ten = xd%t0n + random_exponential(TeLamb)
   
      !Ice strength has a Weibull distribution
   xd%sign = random_weibull (p%miubr, p%varbr) * 1e6
   
      !Ice teeth Delmax and pitch have normal distributions
    xd%signext  = random_weibull (p%miubr, p%varbr) * 1e6
    xd%Dmaxnext = p%miuDelm + p%varDelm ** 0.5 * random_normal() 
    xd%Pchnext  = p%miuP + p%varP ** 0.5 * random_normal() 
    
    IF( xd%Nc == 1 ) THEN
        xd%Dmaxn = p%miuDelm + p%varDelm ** 0.5 * random_normal() 
        xd%Pchn  = p%miuP + p%varP ** 0.5 * random_normal()
    ENDIF
   
    CONTAINS 

!Functions that generate random number with respect to certain distributions

			FUNCTION random_normal() RESULT(fn_val)
		
			! Adapted from the following Fortran 77 code
			!      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
			!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
			!      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.
			
			!  The function random_normal() returns a normally distributed pseudo-random
			!  number with zero mean and unit variance.
			
			!  The algorithm uses the ratio of uniforms method of A.J. Kinderman
			!  and J.F. Monahan augmented with quadratic bounding curves.
			
			REAL :: fn_val
			
			!     Local variables
			REAL     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
			            r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
			
			!     Generate P = (u,v) uniform in rectangle enclosing acceptance region
			
			DO
			  CALL RANDOM_NUMBER(u)
			  CALL RANDOM_NUMBER(v)
			  v = 1.7156 * (v - 0.5)
			
			!     Evaluate the quadratic form
			  x = u - s
			  y = ABS(v) - t
			  q = x**2 + y*(a*y - b*x)
			
			!     Accept P if inside inner ellipse
			  IF (q < r1) EXIT
			!     Reject P if outside outer ellipse
			  IF (q > r2) CYCLE
			!     Reject P if outside acceptance region
			  IF (v**2 < -4.0*LOG(u)*u**2) EXIT
			END DO
			
			!     Return ratio of P's coordinates as the normal deviate
			fn_val = v/u
			RETURN
			
		END FUNCTION random_normal
		
		
		FUNCTION random_exponential(Lambda) RESULT(fn_val)

			! Adapted from Fortran 77 code from the book:
			!     Dagpunar, J. 'Principles of random variate generation'
			!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
			
			! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
			! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
			! TO EXP(-random_exponential), USING INVERSION.
			
			REAL  :: fn_val
			REAL  :: Lambda
			
			!     Local variable
			REAL  :: r
			
			DO
			  CALL RANDOM_NUMBER(r)
			  IF (r > 0.0) EXIT
			END DO
			
			fn_val = -LOG(1-r)/lambda
			RETURN
		
		END FUNCTION random_exponential
		
		
		FUNCTION random_rayleigh(Sigma) RESULT(fn_val)

			! Adapted from Fortran 77 code from the book:
			!     Dagpunar, J. 'Principles of random variate generation'
			!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9
			
			! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
			! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
			! TO EXP(-random_exponential), USING INVERSION.
			
			REAL  :: fn_val
			REAL  :: Sigma
			
			!     Local variable
			REAL  :: r
			
			DO
			  CALL RANDOM_NUMBER(r)
			  IF (r > 0.0) EXIT
			END DO
			
			fn_val = SQRT(-LOG(1-r)*2*sigma**2)
			RETURN
		
		END FUNCTION random_rayleigh
		
		
		FUNCTION random_weibull (mean,var) RESULT(fn_val)
		
			!Function generates a random variate in (0,infinity) from Weibull 
			!distribution with input mean value and variance
			
			IMPLICIT NONE
			
			REAL :: mean 
			REAL :: var 
			REAL :: fn_val
			
			!Local variables
			REAL :: k
			REAL :: Lambda
			REAL :: u
			
			k = wbpar (mean, var)
			lambda = mean / gamma(1+1/k)
			
			CALL RANDOM_NUMBER(u)
			fn_val = lambda * (-log(1-u)) ** (1/k)

		END FUNCTION random_weibull
			
		
		FUNCTION wbpar (mean, var) Result (k1)
		
			!Calculate Weibull distribution parameters due to mean value and variance of the data
			IMPLICIT NONE

			REAL :: mean 
			REAL :: var 
			REAL :: k = 10
            REAL :: k1, F1, dFdk
			REAL :: error = 1e-6
	
			INTEGER :: I
	
			DO i = 1,10000
	
				F1 = (gamma(1+1/k))**2 / gamma(1+2/k) - mean**2/(mean**2+var);
		
				IF (abs(F1) < error) EXIT
		
				dFdk = 2* (gamma(1+1/k))**2 * (-1/k**2) / gamma(1+2/k) * (digamma(1+1/k) -digamma(1+2/k));
        		k = k - F1/dFdk;
	
            END DO
            
            !IF (abs(F1) >= error)THEN
                
            !    WrScr('Weibull parameters never found')
                
            !ENDIF
            
 
		    k1 = k
            
        END FUNCTION wbpar
        
        FUNCTION digamma(z) RESULT(phy)
		
			!Calculate the value of digamma function of z
			REAL, INTENT(IN) :: z
			REAL             :: phy
			
			phy = log(z) - 1/2/z - 1/12/z**2 + 1/120/z**4 - 1/252/z**6 + 1/240/z**8 - 5/660/z**10;
		
		END FUNCTION digamma
		
END SUBROUTINE ID_Generate_RandomStat
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Ice_Split (t, u, p, x, xd, ErrStat, ErrMsg)

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      TYPE(ID_InputType),             INTENT(IN   )  :: u           ! Inputs at t
      TYPE(ID_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(ID_ContinuousStateType),   INTENT(IN   )  :: x           ! Continuous states at t
      TYPE(ID_DiscreteStateType),     INTENT(INOUT)  :: xd          ! Discrete states at t
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
      
      ! Local variable
      REAL(ReKi)                                     :: IceForce 
	  REAL(ReKi)									 :: R
      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Compute outputs here:
      
      R = p%StrWd/2
      
      IF ( xd%Splitf == 0 ) THEN
          
        IF ((x%q - u%q) >= xd%dxc .AND. (x%q - u%q) < xd%dxc+R ) THEN
      
             IceForce =  p%Cpa * ( 2 * p%h * ( R**2 - (R - x%q + u%q)**2 )**0.5 )**( p%dpa + 1 ) * 1.0e6 
         
        ELSE IF (  (x%q - u%q) >= xd%dxc+R ) THEN
          
             IceForce = p%Cpa * ( 2 * R *  p%h )**( p%dpa + 1 ) * 1.0e6 
          
        ELSE
          
             IceForce = 0 
           
        ENDIF
        
      ELSE 
          
           IceForce = 0 
        
      ENDIF
      
      IF ( IceForce >= p%Fsp ) THEN
          xd%Splitf = 1
      ENDIF 

END SUBROUTINE Ice_Split
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_RK4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!   Define constants k1, k2, k3, and k4 as 
!        k1 = dt * f(t        , x_t        )
!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!        k4 = dt * f(t + dt   , x_t + k3   ).
!   Then the continuous states at t = t + dt are
!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!
! For details, see:
! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for 
!   Runge-Kutta." §16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
!   Cambridge University Press, pp. 704-716, 1992.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(ID_InputType),             INTENT(IN   )  :: u(:)        ! Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(ID_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(ID_ContinuousStateType),   INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(ID_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ID_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(ID_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables
         
      TYPE(ID_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states      
      TYPE(ID_ContinuousStateType)                 :: k1          ! RK4 constant; see above
      TYPE(ID_ContinuousStateType)                 :: k2          ! RK4 constant; see above 
      TYPE(ID_ContinuousStateType)                 :: k3          ! RK4 constant; see above 
      TYPE(ID_ContinuousStateType)                 :: k4          ! RK4 constant; see above 
      TYPE(ID_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
      TYPE(ID_InputType)                           :: u_interp    ! interpolated value of inputs 

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 


      ! interpolate u to find u_interp = u(t)
      CALL ID_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat, ErrMsg )

      ! find xdot at t
      CALL ID_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k1%q    = p%dt * xdot%q
      k1%dqdt = p%dt * xdot%dqdt
  
      x_tmp%q    = x%q    + 0.5 * k1%q
      x_tmp%dqdt = x%dqdt + 0.5 * k1%dqdt

      ! interpolate u to find u_interp = u(t + dt/2)
      CALL ID_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%dt, ErrStat, ErrMsg)

      ! find xdot at t + dt/2
      CALL ID_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k2%q    = p%dt * xdot%q
      k2%dqdt = p%dt * xdot%dqdt

      x_tmp%q    = x%q    + 0.5 * k2%q
      x_tmp%dqdt = x%dqdt + 0.5 * k2%dqdt

      ! find xdot at t + dt/2
      CALL ID_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )
     
      k3%q    = p%dt * xdot%q
      k3%dqdt = p%dt * xdot%dqdt

      x_tmp%q    = x%q    + k3%q
      x_tmp%dqdt = x%dqdt + k3%dqdt

      ! interpolate u to find u_interp = u(t + dt)
      CALL ID_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat, ErrMsg)

      ! find xdot at t + dt
      CALL ID_CalcContStateDeriv( t + p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      k4%q    = p%dt * xdot%q
      k4%dqdt = p%dt * xdot%dqdt

      x%q    = x%q    +  ( k1%q    + 2. * k2%q    + 2. * k3%q    + k4%q    ) / 6.      
      x%dqdt = x%dqdt +  ( k1%dqdt + 2. * k2%dqdt + 2. * k3%dqdt + k4%dqdt ) / 6.      

END SUBROUTINE ID_RK4
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_AB4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Adams-Bashforth Method (RK4) for numerically integrating ordinary differential 
! equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!
!   x(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!
!  See, e.g.,
!  http://en.wikipedia.org/wiki/Linear_multistep_method
!
!  or
!
!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(ID_InputType),             INTENT(IN   )  :: u(:)        ! Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(ID_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(ID_ContinuousStateType),   INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(ID_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ID_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(ID_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


      ! local variables
      TYPE(ID_ContinuousStateType) :: xdot       ! Continuous state derivs at t
      TYPE(ID_InputType)           :: u_interp
         

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! need xdot at t
      CALL ID_Input_ExtrapInterp(u, utimes, u_interp, t, ErrStat, ErrMsg)
      CALL ID_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )

      if (n .le. 2) then

         OtherState%n = n

         OtherState%xdot ( 3 - n ) = xdot

         CALL ID_RK4(t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

      else

         if (OtherState%n .lt. n) then

            OtherState%n = n
            OtherState%xdot(4)    = OtherState%xdot(3)
            OtherState%xdot(3)    = OtherState%xdot(2)
            OtherState%xdot(2)    = OtherState%xdot(1)

         elseif (OtherState%n .gt. n) then
 
            ErrStat = ErrID_Fatal
            ErrMsg = ' Backing up in time is not supported with a multistep method '
            RETURN

         endif

         OtherState%xdot ( 1 )     = xdot  ! make sure this is most up to date

         x%q    = x%q    + (p%dt / 24.) * ( 55.*OtherState%xdot(1)%q - 59.*OtherState%xdot(2)%q    + 37.*OtherState%xdot(3)%q  &
                                       - 9. * OtherState%xdot(4)%q )

         x%dqdt = x%dqdt + (p%dt / 24.) * ( 55.*OtherState%xdot(1)%dqdt - 59.*OtherState%xdot(2)%dqdt  &
                                          + 37.*OtherState%xdot(3)%dqdt  - 9.*OtherState%xdot(4)%dqdt )

      endif

END SUBROUTINE ID_AB4
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ID_ABM4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! This subroutine implements the fourth-order Adams-Bashforth-Moulton Method (RK4) for numerically integrating ordinary 
! differential equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!
!   Adams-Bashforth Predictor:
!   x^p(t+dt) = x(t)  + (dt / 24.) * ( 55.*f(t,x) - 59.*f(t-dt,x) + 37.*f(t-2.*dt,x) - 9.*f(t-3.*dt,x) )
!
!   Adams-Moulton Corrector:
!   x(t+dt) = x(t)  + (dt / 24.) * ( 9.*f(t+dt,x^p) + 19.*f(t,x) - 5.*f(t-dt,x) + 1.*f(t-2.*dt,x) )
!
!  See, e.g.,
!  http://en.wikipedia.org/wiki/Linear_multistep_method
!
!  or
!
!  K. E. Atkinson, "An Introduction to Numerical Analysis", 1989, John Wiley & Sons, Inc, Second Edition.
!
!..................................................................................................................................

      REAL(DbKi),                     INTENT(IN   )  :: t           ! Current simulation time in seconds
      INTEGER(IntKi),                 INTENT(IN   )  :: n           ! time step number
      TYPE(ID_InputType),             INTENT(IN   )  :: u(:)        ! Inputs at t
      REAL(DbKi),                     INTENT(IN   )  :: utimes(:)   ! times of input
      TYPE(ID_ParameterType),         INTENT(IN   )  :: p           ! Parameters
      TYPE(ID_ContinuousStateType),   INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
      TYPE(ID_DiscreteStateType),     INTENT(IN   )  :: xd          ! Discrete states at t
      TYPE(ID_ConstraintStateType),   INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
      TYPE(ID_OtherStateType),        INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! local variables

      TYPE(ID_InputType)            :: u_interp        ! Continuous states at t
      TYPE(ID_ContinuousStateType)  :: x_pred          ! Continuous states at t
      TYPE(ID_ContinuousStateType)  :: xdot_pred       ! Continuous states at t

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      CALL ID_CopyContState(x, x_pred, 0, ErrStat, ErrMsg)

      CALL ID_AB4( t, n, u, utimes, p, x_pred, xd, z, OtherState, ErrStat, ErrMsg )

      if (n .gt. 2) then

         CALL ID_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat, ErrMsg)

         CALL ID_CalcContStateDeriv(t + p%dt, u_interp, p, x_pred, xd, z, OtherState, xdot_pred, ErrStat, ErrMsg )

         x%q    = x%q    + (p%dt / 24.) * ( 9. * xdot_pred%q +  19. * OtherState%xdot(1)%q - 5. * OtherState%xdot(2)%q &
                                          + 1. * OtherState%xdot(3)%q )

         x%dqdt = x%dqdt + (p%dt / 24.) * ( 9. * xdot_pred%dqdt + 19. * OtherState%xdot(1)%dqdt - 5. * OtherState%xdot(2)%dqdt &
                                          + 1. * OtherState%xdot(3)%dqdt )
     
      else

         x%q    = x_pred%q
         x%dqdt = x_pred%dqdt

       endif

END SUBROUTINE ID_ABM4
!----------------------------------------------------------------------------------------------------------------------------------
!..................................................................................................................................
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! WE ARE NOT YET IMPLEMENTING THE JACOBIANS...
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
END MODULE IceDyn
!**********************************************************************************************************************************
