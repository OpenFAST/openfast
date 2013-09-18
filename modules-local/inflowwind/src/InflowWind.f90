!**********************************************************************************************************************************
! $Id$
!
! This module is used to read and process the (undisturbed) inflow winds.  It must be initialized
! using InflowWind_Init() with the name of the file, the file type, and possibly reference height and
! width (depending on the type of wind file being used).  This module calls appropriate routines
! in the wind modules so that the type of wind becomes seamless to the user.  IfW_End()
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
!
!..................................................................................................................................
! LICENSING
! Copyright (C) 2012  National Renewable Energy Laboratory
!
!    This file is part of InflowWind.
!
!    InflowWind is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with InflowWind.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
MODULE InflowWind


   USE                              InflowWind_Types   !FIXME: this file will replace SharedInflowDefs when I can get it to work with the framework registry generator.
   USE                              NWTC_Library

      !-------------------------------------------------------------------------------------------------
      ! The included wind modules
      !-------------------------------------------------------------------------------------------------

   USE                              IfW_HHWind_Types           ! Types for IfW_HHWind
   USE                              IfW_HHWind                 ! hub-height text wind files
   USE                              IfW_FFWind_Types           ! Types for IfW_FFWind
   USE                              IfW_FFWind                 ! full-field binary wind files
!   USE                              HAWCWind                   ! full-field binary wind files in HAWC format
!   USE                              FDWind                     ! 4-D binary wind files
!   USE                              CTWind                     ! coherent turbulence from KH billow - binary file superimposed on another wind type
!   USE                              UserWind                   ! user-defined wind module


      !-------------------------------------------------------------------------------------------------
      ! The subroutines
      !-------------------------------------------------------------------------------------------------

   USE                              InflowWind_Subs             ! all the subroutines live here now.




   IMPLICIT NONE
   PRIVATE

   TYPE(ProgDesc), PARAMETER            :: IfW_Ver = ProgDesc( 'InflowWind', 'v1.00.00', '17-Sept-2013' )



      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: IfW_Init                                          ! Initialization routine
   PUBLIC :: IfW_CalcOutput                                    ! Calculate the wind velocities
   PUBLIC :: IfW_End                                           ! Ending routine (includes clean up)

   PUBLIC :: WindInf_ADhack_diskVel

      ! These routines satisfy the framework, but do nothing at present.
   PUBLIC :: IfW_UpdateStates                                  ! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
   PUBLIC :: IfW_CalcConstrStateResidual                       ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: IfW_CalcContStateDeriv                            ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: IfW_UpdateDiscState                               ! Tight coupling routine for updating discrete states


      ! Not coded

!   PUBLIC :: InflowWind_JacobianPInput                 ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations all with respect to the inputs (u)
!   PUBLIC :: InflowWind_JacobianPContState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations all with respect to the continuous states (x)
!   PUBLIC :: InflowWind_JacobianPDiscState             ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations all with respect to the discrete states (xd)
!   PUBLIC :: InflowWind_JacobianPConstrState           ! Routine to compute the Jacobians of the output (Y), continuous- (X), discrete- (Xd), and constraint-state (Z) equations all with respect to the constraint states (z)



CONTAINS
!====================================================================================================
SUBROUTINE IfW_Init( InitData,   InputGuess,    ParamData,                          &
                     ContStates, DiscStates,    ConstrStateGuess,    OtherStates,   &
                     OutData,    TimeInterval,  InitOutData,                        &
                     ErrStat,    ErrMsg )
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
! Since this module acts as an interface to other modules, on some things are set before initiating
! calls to the lower modules.
!----------------------------------------------------------------------------------------------------


!      USE CTWind

         ! Initialization data and guesses

      TYPE( IfW_InitInputType ),          INTENT(IN   )  :: InitData          ! Input data for initialization
      TYPE( IfW_InputType ),              INTENT(  OUT)  :: InputGuess        ! An initial guess for the input; the input mesh must be defined
      TYPE( Ifw_ParameterType ),          INTENT(  OUT)  :: ParamData         ! Parameters
      TYPE( IfW_ContinuousStateType ),    INTENT(  OUT)  :: ContStates        ! Initial continuous states
      TYPE( IfW_DiscreteStateType ),      INTENT(  OUT)  :: DiscStates        ! Initial discrete states
      TYPE( IfW_ConstraintStateType ),    INTENT(  OUT)  :: ConstrStateGuess  ! Initial guess of the constraint states
      TYPE( IfW_OtherStateType ),         INTENT(  OUT)  :: OtherStates       ! Initial other/optimization states
      TYPE( IfW_OutputType ),             INTENT(  OUT)  :: OutData           ! Initial output (outputs are not calculated; only the output mesh is initialized)
      REAL(DbKi),                         INTENT(IN   )  :: TimeInterval      ! Coupling time interval in seconds: InflowWind does not change this.
      TYPE( IfW_InitOutputType ),         INTENT(  OUT)  :: InitOutData       ! Initial output data -- Names, units, and version info.


         ! Error Handling

      INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat           ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg            ! Error message if ErrStat /= ErrID_None


         ! Local variables

      TYPE(IfW_HHWind_InitInputType)                     :: HH_InitData       ! initialization info   
      TYPE(IfW_HHWind_InputType)                         :: HH_InitGuess      ! input positions.      
      TYPE(IfW_HHWind_ContinuousStateType)               :: HH_ContStates     ! Unused                
      TYPE(IfW_HHWind_DiscreteStateType)                 :: HH_DiscStates     ! Unused                
      TYPE(IfW_HHWind_ConstraintStateType)               :: HH_ConstrStates   ! Unused                
      TYPE(IfW_HHWind_OutputType)                        :: HH_OutData        ! output velocities     

      TYPE(IfW_FFWind_InitInputType)                     :: FF_InitData       ! initialization info   
      TYPE(IfW_FFWind_InputType)                         :: FF_InitGuess      ! input positions.      
      TYPE(IfW_FFWind_ContinuousStateType)               :: FF_ContStates     ! Unused                
      TYPE(IfW_FFWind_DiscreteStateType)                 :: FF_DiscStates     ! Unused                
      TYPE(IfW_FFWind_ConstraintStateType)               :: FF_ConstrStates   ! Unused                
      TYPE(IfW_FFWind_OutputType)                        :: FF_OutData        ! output velocities     


!     TYPE(CT_Backgr)                                    :: BackGrndValues


!NOTE: It isn't entirely clear what the purpose of Height is. Does it sometimes occur that Height  /= ParamData%ReferenceHeight???
      REAL(ReKi)                                         :: Height            ! Retrieved from FF
      REAL(ReKi)                                         :: HalfWidth         ! Retrieved from FF

         ! Temporary variables for error handling
      INTEGER(IntKi)                                     :: TmpErrStat
      CHARACTER(1024)                                    :: TmpErrMsg

!NOTE: I may need to revamp how data is passed to the lower modules. Might need to do that before going any further.



         !----------------------------------------------------------------------------------------------
         ! Initialize variables and check to see if this module has been initialized before.
         !----------------------------------------------------------------------------------------------

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! check to see if we are already initialized. Return if it has.
         ! If for some reason a different type of windfile should be used, then call InflowWind_End first, then reinitialize.

      IF ( ParamData%Initialized ) THEN
         ErrMsg   = TRIM(ErrMsg)//NewLine//' Wind inflow has already been initialized.'
         ErrStat  = ErrID_Warn
         RETURN
      ENDIF


         !----------------------------------------------------------------------------------------------
         ! Define the parameters
         !----------------------------------------------------------------------------------------------

      ParamData%DT            = TimeInterval             ! InflowWind does not require a specific time interval, so this is never changed.
      ParamData%WindFileType  = InitData%WindFileType    
      ParamData%WindFileName  = InitData%WindFileName    

      CALL NWTC_Init()                                   ! This might not be needed
      CALL DispNVD( IfW_Ver )                            ! This might be changed later



         !----------------------------------------------------------------------------------------------
         ! State definitions -- only need to define OtherStates, but that can be handled elsewhere.
         !----------------------------------------------------------------------------------------------

         ! At this point in a standard module, we would define the other states (ContStates, DiscStates, etc). Those aren't used here, so we don't.
         ! We would also define the initial guess for the Input_Type, and any meshtypes needed, but we don't need one here.
         ! This is also where we would initialize the output data, but for this module, it will occur within the CalcOutput routine instead.



         !----------------------------------------------------------------------------------------------
         ! Get default wind type, based on file name, if requested. Otherwise store what we are given for the type
         !----------------------------------------------------------------------------------------------

      IF ( InitData%WindFileType == DEFAULT_WindNumber ) THEN
         CALL GetWindType( ParamData, TmpErrStat, TmpErrMsg )
         ErrStat  = MAX( ErrStat, TmpErrStat )
         IF (TmpErrStat /= 0 )         ErrMsg   = TRIM(ErrMsg)//NewLine//TRIM(ErrMsg)
         IF (ErrStat >= AbortErrLev)   RETURN
      ELSE
         ParamData%WindFileType = InitData%WindFileType
      END IF


         !----------------------------------------------------------------------------------------------
         ! Check for coherent turbulence file (KH superimposed on a background wind file)
         ! Initialize the CTWind module and initialize the module of the other wind type.
         !----------------------------------------------------------------------------------------------

      IF ( ParamData%WindFileType == CTP_WindNumber ) THEN

!FIXME: remove this error message when we add CTP_Wind in
         ErrMsg   = TRIM(ErrMsg)//NewLine//'InflowWind cannot currently handle the CTP_Wind type.'
         ErrStat  = ErrID_Fatal

!         CALL CT_Init(UnWind, ParamData%WindFileName, BackGrndValues, ErrStat, ErrMsg)
!         IF (ErrStat /= 0) THEN
!   !         CALL IfW_End( ParamData, ErrStat )
!   !FIXME: cannot call IfW_End here -- requires InitData to be INOUT. Not allowed by framework.
!   !         CALL IfW_End( InitData, ParamData, ContStates, DiscStates, ConstrStateGuess, OtherStates, &
!   !                       OutData, ErrStat, ErrMsg )
!            ParamData%WindFileType = Undef_Wind
!            ErrStat  = 1
!            RETURN
!         END IF
!
!   !FIXME: check this
!         ParamData%WindFileName = BackGrndValues%WindFile
!         ParamData%WindFileType = BackGrndValues%WindFileType
!   !      CT_Flag  = BackGrndValues%CoherentStr
!         ParamData%CT_Flag  = BackGrndValues%CoherentStr    ! This might be wrong

      ELSE

         ParamData%CT_Flag  = .FALSE.

      END IF

         !----------------------------------------------------------------------------------------------
         ! Initialize based on the wind type
         !----------------------------------------------------------------------------------------------

      SELECT CASE ( ParamData%WindFileType )

         CASE (HH_WindNumber)

            HH_InitData%ReferenceHeight = InitData%ReferenceHeight
            HH_InitData%Width           = InitData%Width
            HH_InitData%WindFileName    = InitData%WindFileName

            CALL IfW_HHWind_Init(HH_InitData,   HH_InitGuess,  ParamData%HHWind,                         &
                                 HH_ContStates, HH_DiscStates, HH_ConstrStates,     OtherStates%HHWind,  &
                                 HH_OutData,    TimeInterval,  InitOutData%HHWind,  TmpErrStat,          TmpErrMsg)

            ErrStat = MAX( ErrStat, TmpErrStat )
            IF ( TmpErrStat /= ErrID_None )  ErrMsg = TRIM(ErrMsg)//NewLine//TRIM(TmpErrMsg)

              ! Copy Relevant info over to InitOutData
            IF ( TmpErrStat < AbortErrLev ) THEN

                  ! Allocate and copy over the WriteOutputHdr info
               CALL AllocAry( InitOutData%WriteOutputHdr, SIZE(InitOutData%HHWind%WriteOutputHdr,1), &
                              'Empty array for names of outputable information.', TmpErrStat, TmpErrMsg )
               ErrStat  = MAX(TmpErrStat, ErrStat)
               ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
               IF ( ErrStat >= AbortErrLev ) RETURN
               InitOutData%WriteOutputHdr    =  InitOutData%HHWind%WriteOutputHdr

                  ! Allocate and copy over the WriteOutputUnt info
               CALL AllocAry( InitOutData%WriteOutputUnt, SIZE(InitOutData%HHWind%WriteOutputUnt,1), &
                              'Empty array for units of outputable information.', TmpErrStat, TmpErrMsg )
               ErrStat  = MAX(TmpErrStat, ErrStat)
               ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
               IF ( ErrStat >= AbortErrLev ) RETURN

               InitOutData%WriteOutputUnt    =  InitOutData%HHWind%WriteOutputUnt

                  ! Copy the hub height info over
               InitOutData%HubHeight         =  InitOutData%HHWind%HubHeight
            ENDIF

!           IF (CT_Flag) CALL CT_SetRefVal(FileInfo%ReferenceHeight, 0.5*FileInfo%Width, ErrStat)  !FIXME: check if this was originally used
!           IF (ErrStat == ErrID_None .AND. ParamData%CT_Flag) &
!              CALL CT_SetRefVal(InitData%ReferenceHeight, REAL(0.0, ReKi), ErrStat, ErrMsg)      !FIXME: will need to put this routine in the Init of CT


         CASE (FF_WindNumber)

            FF_InitData%ReferenceHeight = InitData%ReferenceHeight
            FF_InitData%Width           = InitData%Width
            FF_InitData%WindFileName    = InitData%WindFileName

            CALL IfW_FFWind_Init(FF_InitData,   FF_InitGuess,  ParamData%FFWind,                         &
                                 FF_ContStates, FF_DiscStates, FF_ConstrStates,     OtherStates%FFWind,  &
                                 FF_OutData,    TimeInterval,  InitOutData%FFWind,  TmpErrStat,          TmpErrMsg)

            ErrStat = MAX( ErrStat, TmpErrStat )
            IF ( TmpErrStat /= ErrID_None )  ErrMsg = TRIM(ErrMsg)//NewLine//TRIM(TmpErrMsg)

              ! Copy Relevant info over to InitOutData
            IF ( TmpErrStat < AbortErrLev ) THEN

                  ! Allocate and copy over the WriteOutputHdr info
               CALL AllocAry( InitOutData%WriteOutputHdr, SIZE(InitOutData%FFWind%WriteOutputHdr,1), &
                              'Empty array for names of outputable information.', TmpErrStat, TmpErrMsg )
               ErrStat  = MAX(TmpErrStat, ErrStat)
               ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
               IF ( ErrStat >= AbortErrLev ) RETURN
               InitOutData%WriteOutputHdr    =  InitOutData%FFWind%WriteOutputHdr

                  ! Allocate and copy over the WriteOutputUnt info
               CALL AllocAry( InitOutData%WriteOutputUnt, SIZE(InitOutData%FFWind%WriteOutputUnt,1), &
                              'Empty array for units of outputable information.', TmpErrStat, TmpErrMsg )
               ErrStat  = MAX(TmpErrStat, ErrStat)
               ErrMsg   = TRIM(ErrMsg)//TRIM(TmpErrMsg)//NewLine
               IF ( ErrStat >= AbortErrLev ) RETURN

               InitOutData%WriteOutputUnt    =  InitOutData%FFWind%WriteOutputUnt

                  ! Copy the hub height info over
               InitOutData%HubHeight         =  InitOutData%FFWind%HubHeight
            ENDIF

            !FIXME: Fix this when CT_Wind is available
!               ! Set CT parameters
!            IF ( ErrStat == ErrID_None .AND. ParamData%CT_Flag ) THEN
!               Height     = FF_GetValue('HubHeight', ErrStat, ErrMsg)
!               IF ( ErrStat /= 0 ) Height = InitData%ReferenceHeight
!
!               HalfWidth  = 0.5*FF_GetValue('GridWidth', ErrStat, ErrMsg)
!               IF ( ErrStat /= 0 ) HalfWidth = 0
!
!               CALL CT_SetRefVal(Height, HalfWidth, ErrStat, ErrMsg)
!            END IF


         CASE (UD_WindNumber)

               !FIXME: remove this error message when we add UD_Wind in
            ErrMsg   = TRIM(ErrMsg)//NewLine//'InflowWind cannot currently handle the UD_Wind type.'
            ErrStat  = ErrID_Fatal

!            CALL UsrWnd_Init(ErrStat)


         CASE (FD_WindNumber)

               !FIXME: remove this error message when we add FD_Wind in
            ErrMsg   = TRIM(ErrMsg)//NewLine//'InflowWind cannot currently handle the FD_Wind type.'
            ErrStat  = ErrID_Fatal

!            CALL IfW_FDWind_Init(UnWind, ParamData%WindFileName, InitData%ReferenceHeight, ErrStat)


         CASE (HAWC_WindNumber)

               !FIXME: remove this error message when we add HAWC_Wind in
            ErrMsg   = TRIM(ErrMsg)//NewLine//'InflowWind cannot currently handle the HAWC_Wind type.'
            ErrStat  = ErrID_Fatal

!            CALL HW_Init( UnWind, ParamData%WindFileName, ErrStat )


         CASE DEFAULT

            ErrMsg   = TRIM(ErrMsg)//NewLine//' Error: Undefined wind type in WindInflow_Init()'
            ErrStat  = ErrID_Fatal

      END SELECT


         ! check error status. If the error is recoverable (ErrID_Warn or less), set flag to indicate we are initialized.

      IF ( ErrStat >= ErrID_Severe ) THEN
         ParamData%Initialized   = .FALSE.
         ParamData%WindFileType  = Undef_WindNumber

            ! Just in case something was allocated
         CALL IfW_End(  InputGuess,    ParamData,                                       &
                        ContStates,    DiscStates,    ConstrStateGuess,    OtherStates, &
                        OutData,       TmpErrStat,    TmpErrMsg )
         ErrStat = MAX( ErrStat, TmpErrStat )
         IF ( TmpErrStat /= ErrID_None )  ErrMsg = TRIM(ErrMsg)//NewLine//TRIM(TmpErrMsg)

      ELSE
         ParamData%Initialized = .TRUE.
      END IF


         ! Set the version information in InitOutData
      InitOutData%Ver   = IfW_Ver



      RETURN

 END SUBROUTINE IfW_Init
 !====================================================================================================
 SUBROUTINE IfW_CalcOutput( Time, InputData, ParamData, &
                              ContStates, DiscStates, ConstrStates, OtherStates, &   ! States -- none in this case
                              OutputData, ErrStat, ErrMsg )
   ! This routine takes an input dataset of type InputType which contains a position array of dimensions 3*n. It then calculates
   ! and returns the output dataset of type OutputType which contains a corresponding velocity array of dimensions 3*n. The input

   ! array contains XYZ triplets for each position of interest (first index is X/Y/Z for values 1/2/3, second index is the point
   ! number to evaluate). The returned values in the OutputData are similar with U/V/W for the first index of 1/2/3.
   !----------------------------------------------------------------------------------------------------

         ! Inputs / Outputs

      REAL( DbKi ),                       INTENT(IN   )  :: Time              ! Current simulation time in seconds
      TYPE( IfW_InputType ),              INTENT(IN   )  :: InputData         ! Inputs at Time
      TYPE( Ifw_ParameterType ),          INTENT(IN   )  :: ParamData         ! Parameters
      TYPE( IfW_ContinuousStateType ),    INTENT(IN   )  :: ContStates        ! Continuous states at Time
      TYPE( IfW_DiscreteStateType ),      INTENT(IN   )  :: DiscStates        ! Discrete states at Time
      TYPE( IfW_ConstraintStateType ),    INTENT(IN   )  :: ConstrStates      ! Constraint states at Time
      TYPE( IfW_OtherStateType ),         INTENT(INOUT)  :: OtherStates       ! Other/optimization states at Time
      TYPE( IfW_OutputType ),             INTENT(  OUT)  :: OutputData        ! Outputs computed at Time (IN for mesh reasons -- not used here)

      INTEGER( IntKi ),                   INTENT(  OUT)  :: ErrStat           ! Error status of the operation
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg            ! Error message if ErrStat /= ErrID_None


         ! Local variables

      TYPE(IfW_HHWind_InitInputType)                     :: HH_InitData       ! initialization info  
      TYPE(IfW_HHWind_InputType)                         :: HH_InData         ! input positions.      
      TYPE(IfW_HHWind_ContinuousStateType)               :: HH_ContStates     ! Unused               
      TYPE(IfW_HHWind_DiscreteStateType)                 :: HH_DiscStates     ! Unused               
      TYPE(IfW_HHWind_ConstraintStateType)               :: HH_ConstrStates   ! Unused               
      TYPE(IfW_HHWind_OutputType)                        :: HH_OutData        ! output velocities   

      TYPE(IfW_FFWind_InitInputType)                     :: FF_InitData       ! initialization info 
      TYPE(IfW_FFWind_InputType)                         :: FF_InData         ! input positions.     
      TYPE(IfW_FFWind_ContinuousStateType)               :: FF_ContStates     ! Unused             
      TYPE(IfW_FFWind_DiscreteStateType)                 :: FF_DiscStates     ! Unused               
      TYPE(IfW_FFWind_ConstraintStateType)               :: FF_ConstrStates   ! Unused               
      TYPE(IfW_FFWind_OutputType)                        :: FF_OutData        ! output velocities    




!NOTE: It isn't entirely clear what the purpose of Height is. Does it sometimes occur that Height  /= ParamData%ReferenceHeight???
      REAL(ReKi)                                         :: Height      ! Retrieved from FF
      REAL(ReKi)                                         :: HalfWidth   ! Retrieved from FF


         ! Sub modules use the InflIntrpOut derived type to store the wind information
!     TYPE(CT_Backgr)                                    :: BackGrndValues
!      TYPE(InflIntrpOut)                                 :: CTWindSpeed       ! U, V, W velocities to superimpose on background wind
!      TYPE(InflIntrpOut)                                 :: TempWindSpeed     ! U, V, W velocities returned
!      REAL(ReKi)                                         :: CTWindSpeed(3)     ! U, V, W velocities to superimpose on background wind
!      REAL(ReKi)                                         :: TempWindSpeed(3)   ! Temporary U, V, W velocities



         ! Temporary variables for error handling
      INTEGER(IntKi)                                     :: TmpErrStat
      CHARACTER(1024)                                    :: TmpErrMsg



         ! Initialize ErrStat
      ErrStat  = ErrID_None
      ErrMsg   = ""


         ! Allocate the velocity array to get out
      CALL AllocAry( OutputData%Velocity, 3, SIZE(InputData%Position,2), &
                     "Velocity array returned from IfW_CalcOutput", TmpErrStat, TmpErrMsg )
      IF (TmpErrStat /= ErrID_None) THEN
         ErrMsg   = TRIM(TmpErrMsg)
         ErrStat  = MAX(ErrStat,TmpErrStat)
      ENDIF
      IF (ErrStat >= AbortErrLev) RETURN

         ! Compute the wind velocities by stepping through all the data points and calling the appropriate GetWindSpeed routine
      SELECT CASE ( ParamData%WindFileType )
         CASE (HH_WindNumber)

               ! Allocate the position array to pass in
            CALL AllocAry( HH_InData%Position, 3, SIZE(InputData%Position,2), &
                           "Position grid for passing to IfW_HHWind_CalcOutput", TmpErrStat, TmpErrMsg )
            IF (TmpErrStat /= ErrID_None) THEN
               ErrMsg   = TRIM(ErrMsg)//NewLine//TRIM(TmpErrMsg)
               ErrStat  = MAX(ErrStat,TmpErrStat)
            ENDIF
            IF (ErrStat >= AbortErrLev) RETURN
               ! Copy positions over
            HH_InData%Position   = InputData%Position

            CALL  IfW_HHWind_CalcOutput(  Time,          HH_InData,     ParamData%HHWind,                         &
                                          HH_ContStates, HH_DiscStates, HH_ConstrStates,     OtherStates%HHWind,  &
                                          HH_OutData,    TmpErrStat,    TmpErrMsg)

               ! Copy the velocities over
            OutputData%Velocity  = HH_OutData%Velocity


         CASE (FF_WindNumber)

               ! Allocate the position array to pass in
            CALL AllocAry( FF_InData%Position, 3, SIZE(InputData%Position,2), &
                           "Position grid for passing to IfW_FFWind_CalcOutput", TmpErrStat, TmpErrMsg )
            IF (TmpErrStat /= ErrID_None) THEN
               ErrMsg   = TRIM(ErrMsg)//NewLine//TRIM(TmpErrMsg)
               ErrStat  = MAX(ErrStat,TmpErrStat)
            ENDIF
            IF (ErrStat >= AbortErrLev) RETURN
               ! Copy positions over
            FF_InData%Position   = InputData%Position

            CALL  IfW_FFWind_CalcOutput(  Time,          FF_InData,     ParamData%FFWind,                         &
                                          FF_ContStates, FF_DiscStates, FF_ConstrStates,     OtherStates%FFWind,  &
                                          FF_OutData,    TmpErrStat,    TmpErrMsg)

               ! Copy the velocities over
            OutputData%Velocity  = FF_OutData%Velocity

!               OutputData%Velocity(:,PointCounter) = FF_GetWindSpeed(     Time, InputData%Position(:,PointCounter), ErrStat, ErrMsg)


!         CASE (UD_WindNumber)

!               OutputData%Velocity(:,PointCounter) = UsrWnd_GetWindSpeed( Time, InputData%Position(:,PointCounter), ErrStat )!, ErrMsg)


!         CASE (FD_WindNumber)

!               OutputData%Velocity(:,PointCounter) = FD_GetWindSpeed(     Time, InputData%Position(:,PointCounter), ErrStat )



!         CASE (HAWC_WindNumber)

!               OutputData%Velocity(:,PointCounter) = HW_GetWindSpeed(     Time, InputData%Position(:,PointCounter), ErrStat )



            ! If it isn't one of the above cases, we have a problem and won't be able to continue

         CASE DEFAULT

            ErrMsg = ' Error: Undefined wind type in IfW_CalcOutput. ' &
                      //'Call WindInflow_Init() before calling this function.'
            ErrStat = ErrID_Fatal               ! No data returned
            OutputData%Velocity(:,:) = 0.0
            RETURN

      END SELECT


         ! If we had a severe or fatal error, we need to make sure we zero out the result and return.

      IF (ErrStat >= ErrID_Severe) THEN

         OutputData%Velocity(:,:) = 0.0
         RETURN

      ELSE

            ! Add coherent turbulence to background wind

!         IF (ParamData%CT_Flag) THEN
!
!            DO PointCounter = 1, SIZE(InputData%Position, 2)
!
!               TempWindSpeed = CT_GetWindSpeed(     Time, InputData%Position(:,PointCounter), ErrStat, ErrMsg )
!
!                  ! Error Handling -- move ErrMsg inside CT_GetWindSPeed and simplify
!               IF (ErrStat >= ErrID_Severe) THEN
!                  ErrMsg   = 'IfW_CalcOutput: Error in CT_GetWindSpeed for point number '//TRIM(Num2LStr(PointCounter))
!                  EXIT        ! Exit the loop
!               ENDIF
!
!               OutputData%Velocity(:,PointCounter) = OutputData%Velocity(:,PointCounter) + TempWindSpeed
!
!            ENDDO
!
!               ! If something went badly wrong, Return
!            IF (ErrStat >= ErrID_Severe ) RETURN
!
!         ENDIF
!
      ENDIF



END SUBROUTINE IfW_CalcOutput
!====================================================================================================
SUBROUTINE IfW_End( InitData, ParamData, ContStates, DiscStates, ConstrStateGuess, OtherStates, &
                       OutData, ErrStat, ErrMsg )
   ! Clean up the allocated variables and close all open files.  Reset the initialization flag so
   ! that we have to reinitialize before calling the routines again.
   !----------------------------------------------------------------------------------------------------

         ! Initialization data and guesses

      TYPE( IfW_InputType ),              INTENT(INOUT)  :: InitData          ! Input data for initialization
      TYPE( Ifw_ParameterType ),          INTENT(INOUT)  :: ParamData         ! Parameters
      TYPE( IfW_ContinuousStateType ),    INTENT(INOUT)  :: ContStates        ! Continuous states
      TYPE( IfW_DiscreteStateType ),      INTENT(INOUT)  :: DiscStates        ! Discrete states
      TYPE( IfW_ConstraintStateType ),    INTENT(INOUT)  :: ConstrStateGuess  ! Guess of the constraint states
      TYPE( IfW_OtherStateType ),         INTENT(INOUT)  :: OtherStates       ! Other/optimization states
      TYPE( IfW_OutputType ),             INTENT(INOUT)  :: OutData           ! Output data


         ! Error Handling

      INTEGER( IntKi ),                   INTENT(  OUT)  :: ErrStat
      CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg

         ! Local variables

      TYPE(IfW_HHWind_InputType)                         :: HH_InitData       ! input positions.     
      TYPE(IfW_HHWind_ContinuousStateType)               :: HH_ContStates     ! Unused              
      TYPE(IfW_HHWind_DiscreteStateType)                 :: HH_DiscStates     ! Unused             
      TYPE(IfW_HHWind_ConstraintStateType)               :: HH_ConstrStates   ! Unused              
      TYPE(IfW_HHWind_OutputType)                        :: HH_OutData        ! output velocities    

      TYPE(IfW_FFWind_InputType)                         :: FF_InitData       ! input positions.   
      TYPE(IfW_FFWind_ContinuousStateType)               :: FF_ContStates     ! Unused              
      TYPE(IfW_FFWind_DiscreteStateType)                 :: FF_DiscStates     ! Unused             
      TYPE(IfW_FFWind_ConstraintStateType)               :: FF_ConstrStates   ! Unused               
      TYPE(IfW_FFWind_OutputType)                        :: FF_OutData        ! output velocities    


!     TYPE(CT_Backgr)                                    :: BackGrndValues


!NOTE: It isn't entirely clear what the purpose of Height is. Does it sometimes occur that Height  /= ParamData%ReferenceHeight???
      REAL(ReKi)                                         :: Height      ! Retrieved from FF
      REAL(ReKi)                                         :: HalfWidth   ! Retrieved from FF



         ! End the sub-modules (deallocates their arrays and closes their files):

      SELECT CASE ( ParamData%WindFileType )

         CASE (HH_WindNumber)
            CALL IfW_HHWind_End( HH_InitData,   ParamData%HHWind,                                        &
                                 HH_ContStates, HH_DiscStates,    HH_ConstrStates,  OtherStates%HHWind,  &
                                 HH_OutData,    ErrStat,          ErrMsg )

         CASE (FF_WindNumber)
            CALL IfW_FFWind_End( FF_InitData,   ParamData%FFWind,                                        &
                                 FF_ContStates, FF_DiscStates,    FF_ConstrStates,  OtherStates%FFWind,  &
                                 FF_OutData,    ErrStat,          ErrMsg )

!         CASE (UD_WindNumber)
!            CALL UsrWnd_Terminate( ErrStat )

!         CASE (FD_WindNumber)
!            CALL FD_Terminate(     ErrStat )

!         CASE (HAWC_WindNumber)
!            CALL HW_Terminate(     ErrStat )

         CASE ( Undef_WindNumber )
            ! Do nothing

         CASE DEFAULT  ! keep this check to make sure that all new wind types have a terminate function
            ErrMsg   = TRIM(ErrMsg)//NewLine//' InflowWind: Undefined wind type in IfW_End().'
            ErrStat  = ErrID_Severe

      END SELECT

!  !   IF (CT_Flag) CALL CT_Terminate( ErrStat ) !FIXME: should it be this line or the next?
!         CALL CT_Terminate( ErrStat, ErrMsg )


         ! Reset the wind type so that the initialization routine must be called
      ParamData%WindFileType = Undef_WindNumber
      ParamData%Initialized = .FALSE.
      ParamData%CT_Flag  = .FALSE.


END SUBROUTINE IfW_End
!====================================================================================================
! The following routines were added to satisfy the framework, but do nothing useful.
!====================================================================================================
SUBROUTINE IfW_UpdateStates( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input Time; Continuous and discrete states are updated for Time + Interval
!..................................................................................................................................

      REAL(DbKi),                            INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(IfW_InputType),                   INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(IfW_ParameterType),               INTENT(IN   )  :: p           ! Parameters
      TYPE(IfW_ContinuousStateType),         INTENT(INOUT)  :: x           ! Input: Continuous states at Time;
                                                                           ! Output: Continuous states at Time + Interval
      TYPE(IfW_DiscreteStateType),           INTENT(INOUT)  :: xd          ! Input: Discrete states at Time;
                                                                           ! Output: Discrete states at Time  + Interval
      TYPE(IfW_ConstraintStateType),         INTENT(INOUT)  :: z           ! Input: Initial guess of constraint states at Time;
                                                                           ! Output: Constraint states at Time
      TYPE(IfW_OtherStateType),              INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

         ! Local variables

      TYPE(IfW_ContinuousStateType)                         :: dxdt        ! Continuous state derivatives at Time
      TYPE(IfW_ConstraintStateType)                         :: z_Residual  ! Residual of the constraint state equations (Z)

      INTEGER(IntKi)                                        :: ErrStat2    ! Error status of the operation (occurs after initial error)
      CHARACTER(LEN(ErrMsg))                                :: ErrMsg2     ! Error message if ErrStat2 /= ErrID_None

         ! Initialize ErrStat

      ErrStat = ErrID_Warn
      ErrMsg  = "IfW_UpdateStates was called.  That routine does nothing useful."



         ! Solve for the constraint states (z) here:

         ! Check if the z guess is correct and update z with a new guess.
         ! Iterate until the value is within a given tolerance.

      CALL IfW_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_Residual, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL IfW_DestroyConstrState( z_Residual, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN
      ENDIF

      ! DO WHILE ( z_Residual% > tolerance )
      !
      !  z =
      !
      !  CALL IfW_FFWind_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_Residual, ErrStat, ErrMsg )
      !  IF ( ErrStat >= AbortErrLev ) THEN
      !     CALL IfW_FFWind_DestroyConstrState( z_Residual, ErrStat2, ErrMsg2)
      !     ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
      !     RETURN
      !  ENDIF
      !
      ! END DO


         ! Destroy z_Residual because it is not necessary for the rest of the subroutine:

      CALL IfW_DestroyConstrState( z_Residual, ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN



         ! Get first time derivatives of continuous states (dxdt):

      CALL IfW_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL IfW_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN
      ENDIF


         ! Update discrete states:
         !   Note that xd [discrete state] is changed in IfW_FFWind_UpdateDiscState(), so IfW_FFWind_CalcOutput(),
         !   IfW_FFWind_CalcContStateDeriv(), and IfW_FFWind_CalcConstrStates() must be called first (see above).

      CALL IfW_UpdateDiscState(Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
      IF ( ErrStat >= AbortErrLev ) THEN
         CALL IfW_DestroyContState( dxdt, ErrStat2, ErrMsg2)
         ErrMsg = TRIM(ErrMsg)//' '//TRIM(ErrMsg2)
         RETURN
      ENDIF


         ! Integrate (update) continuous states (x) here:

      !x = function of dxdt and x


         ! Destroy dxdt because it is not necessary for the rest of the subroutine

      CALL IfW_DestroyContState( dxdt, ErrStat, ErrMsg)
      IF ( ErrStat >= AbortErrLev ) RETURN



END SUBROUTINE IfW_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IfW_CalcContStateDeriv( Time, u, p, x, xd, z, OtherState, dxdt, ErrStat, ErrMsg )
! Tight coupling routine for computing derivatives of continuous states
!..................................................................................................................................

      REAL(DbKi),                            INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(IfW_InputType),                   INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(IfW_ParameterType),               INTENT(IN   )  :: p           ! Parameters
      TYPE(IfW_ContinuousStateType),         INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(IfW_DiscreteStateType),           INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(IfW_ConstraintStateType),         INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(IfW_OtherStateType),              INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(IfW_ContinuousStateType),         INTENT(  OUT)  :: dxdt        ! Continuous state derivatives at Time
      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_Warn
      ErrMsg  = "IfW_CalcContStateDrive was called.  That routine does nothing."


         ! Compute the first time derivatives of the continuous states here:

      dxdt%DummyContState = 0


END SUBROUTINE IfW_CalcContStateDeriv
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IfW_UpdateDiscState( Time, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
! Tight coupling routine for updating discrete states
!..................................................................................................................................

      REAL(DbKi),                            INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(IfW_InputType),                   INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(IfW_ParameterType),               INTENT(IN   )  :: p           ! Parameters
      TYPE(IfW_ContinuousStateType),         INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(IfW_DiscreteStateType),           INTENT(INOUT)  :: xd          ! Input: Discrete states at Time;
                                                                           !   Output: Discrete states at Time + Interval
      TYPE(IfW_ConstraintStateType),         INTENT(IN   )  :: z           ! Constraint states at Time
      TYPE(IfW_OtherStateType),              INTENT(INOUT)  :: OtherState  ! Other/optimization states
      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_Warn
      ErrMsg  = "IfW_UpdateDiscState was called.  That routine does nothing."


         ! Update discrete states here:

      ! StateData%DiscState =

END SUBROUTINE IfW_UpdateDiscState
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE IfW_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, z_residual, ErrStat, ErrMsg )
! Tight coupling routine for solving for the residual of the constraint state equations
!..................................................................................................................................

      REAL(DbKi),                            INTENT(IN   )  :: Time        ! Current simulation time in seconds
      TYPE(IfW_InputType),                   INTENT(IN   )  :: u           ! Inputs at Time
      TYPE(IfW_ParameterType),               INTENT(IN   )  :: p           ! Parameters
      TYPE(IfW_ContinuousStateType),         INTENT(IN   )  :: x           ! Continuous states at Time
      TYPE(IfW_DiscreteStateType),           INTENT(IN   )  :: xd          ! Discrete states at Time
      TYPE(IfW_ConstraintStateType),         INTENT(IN   )  :: z           ! Constraint states at Time (possibly a guess)
      TYPE(IfW_OtherStateType),              INTENT(INOUT)  :: OtherState  ! Other/optimization states
      TYPE(IfW_ConstraintStateType),         INTENT(  OUT)  :: z_residual  ! Residual of the constraint state equations using
                                                                           ! the input values described above
      INTEGER(IntKi),                        INTENT(  OUT)  :: ErrStat     ! Error status of the operation
      CHARACTER(*),                          INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "IfW_CalcConstrStateResidual was called.  That routine does nothing."


         ! Solve for the constraint states here:

      z_residual%DummyConstrState = 0

END SUBROUTINE IfW_CalcConstrStateResidual
!====================================================================================================
!====================================================================================================
FUNCTION WindInf_ADhack_diskVel( Time,ParamData, OtherStates,ErrStat, ErrMsg )
! This function should be deleted ASAP.  It's purpose is to reproduce results of AeroDyn 12.57;
! when a consensus on the definition of "average velocity" is determined, this function will be
! removed.  InpPosition(2) should be the rotor radius; InpPosition(3) should be hub height
!----------------------------------------------------------------------------------------------------

      ! Passed variables

   REAL(DbKi),                  INTENT(IN)     :: Time
   TYPE( Ifw_ParameterType ),   INTENT(IN)     :: ParamData         ! Parameters
   TYPE( IfW_OtherStateType ),  INTENT(INOUT)  :: OtherStates       ! Other/optimization states

   INTEGER(intKi), INTENT(OUT)       :: ErrStat
   CHARACTER(*),INTENT(OUT)   :: ErrMsg

      ! Function definition
   REAL(ReKi)                 :: WindInf_ADhack_diskVel(3)

      ! Local variables
   REAL(ReKi)                    :: Delta_tmp            ! interpolated Delta   at input TIME
   REAL(ReKi)                    :: P                    ! temporary storage for slope (in time) used in linear interpolation
   REAL(ReKi)                    :: V_tmp                ! interpolated V       at input TIME
   REAL(ReKi)                    :: VZ_tmp               ! interpolated VZ      at input TIME

   
   
   ErrStat = ErrID_None

   SELECT CASE ( ParamData%WindFileType )
      CASE (HH_WindNumber)
            
         !-------------------------------------------------------------------------------------------------
         ! Linearly interpolate in time (or use nearest-neighbor to extrapolate) 
         ! (compare with NWTC_Num.f90\InterpStpReal)
         !-------------------------------------------------------------------------------------------------

      
            ! Let's check the limits.
         IF ( Time <= OtherStates%HHWind%Tdata(1) .OR. OtherStates%HHWind%NumDataLines == 1 )  THEN

            OtherStates%HHWind%TimeIndex      = 1
            V_tmp         = OtherStates%HHWind%V      (1)
            Delta_tmp     = OtherStates%HHWind%Delta  (1)
            VZ_tmp        = OtherStates%HHWind%VZ     (1)

         ELSE IF ( Time >= OtherStates%HHWind%Tdata(OtherStates%HHWind%NumDataLines) )  THEN

            OtherStates%HHWind%TimeIndex = OtherStates%HHWind%NumDataLines - 1
            V_tmp                 = OtherStates%HHWind%V      (OtherStates%HHWind%NumDataLines)
            Delta_tmp             = OtherStates%HHWind%Delta  (OtherStates%HHWind%NumDataLines)
            VZ_tmp                = OtherStates%HHWind%VZ     (OtherStates%HHWind%NumDataLines)
               
         ELSE
   
              ! Let's interpolate!

            OtherStates%HHWind%TimeIndex = MAX( MIN( OtherStates%HHWind%TimeIndex, OtherStates%HHWind%NumDataLines-1 ), 1 )

            DO

               IF ( Time < OtherStates%HHWind%Tdata(OtherStates%HHWind%TimeIndex) )  THEN

                  OtherStates%HHWind%TimeIndex = OtherStates%HHWind%TimeIndex - 1

               ELSE IF ( Time >= OtherStates%HHWind%Tdata(OtherStates%HHWind%TimeIndex+1) )  THEN

                  OtherStates%HHWind%TimeIndex = OtherStates%HHWind%TimeIndex + 1

               ELSE
                  P           = ( Time - OtherStates%HHWind%Tdata(OtherStates%HHWind%TimeIndex) )/( OtherStates%HHWind%Tdata(OtherStates%HHWind%TimeIndex+1) &
                                 - OtherStates%HHWind%Tdata(OtherStates%HHWind%TimeIndex) )
                  V_tmp       = ( OtherStates%HHWind%V(      OtherStates%HHWind%TimeIndex+1) - OtherStates%HHWind%V(      OtherStates%HHWind%TimeIndex) )*P  &
                                + OtherStates%HHWind%V(      OtherStates%HHWind%TimeIndex)
                  Delta_tmp   = ( OtherStates%HHWind%Delta(  OtherStates%HHWind%TimeIndex+1) - OtherStates%HHWind%Delta(  OtherStates%HHWind%TimeIndex) )*P  &
                                + OtherStates%HHWind%Delta(  OtherStates%HHWind%TimeIndex)
                  VZ_tmp      = ( OtherStates%HHWind%VZ(     OtherStates%HHWind%TimeIndex+1) - OtherStates%HHWind%VZ(     OtherStates%HHWind%TimeIndex) )*P  &
                                + OtherStates%HHWind%VZ(     OtherStates%HHWind%TimeIndex)
                  EXIT

               END IF

            END DO

         END IF         
                                      
      !-------------------------------------------------------------------------------------------------
      ! calculate the wind speed at this time
      !-------------------------------------------------------------------------------------------------
      
         WindInf_ADhack_diskVel(1) =  V_tmp * COS( Delta_tmp )
         WindInf_ADhack_diskVel(2) = -V_tmp * SIN( Delta_tmp )
         WindInf_ADhack_diskVel(3) =  VZ_tmp      
                                  


      CASE (FF_WindNumber)

         WindInf_ADhack_diskVel(1)   = OtherStates%FFWind%MeanFFWS
         WindInf_ADhack_diskVel(2:3) = 0.0                 

      CASE DEFAULT
         ErrStat = ErrID_Fatal
         ErrMsg = ' WindInf_ADhack_diskVel: Undefined wind type.'

   END SELECT

   RETURN

END FUNCTION WindInf_ADhack_diskVel                       
                       



!====================================================================================================
END MODULE InflowWind

!!----Removed during conversion to new framework: may put back in as part of OtherStates
!!       PUBLIC                         :: InflowWind_GetMean        ! function to get the mean wind speed at a point in space
!!       PUBLIC                         :: InflowWind_GetStdDev      ! function to calculate standard deviation at a point in space
!!       PUBLIC                         :: InflowWind_GetTI          ! function to get TI at a point in space
