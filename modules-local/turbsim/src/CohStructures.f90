!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2014, 2016  National Renewable Energy Laboratory
!
!    This file is part of TurbSim.
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
MODULE TS_CohStructures

   USE TurbSim_Types
   
   use TS_Profiles
   use TS_RandNum

   IMPLICIT NONE

   
   REAL(ReKi), PARAMETER   ::  KHT_LES_dT = 0.036335           ! The average time step in the LES test file, used here for the KH test
   REAL(ReKi), PARAMETER   ::  KHT_LES_Zm = 6.35475            ! The non-dimensional z dimension defined in LES test file, used here for the KH test
   
   TYPE                         :: Event                                    ! Coherent turbulent event to add to the background wind
      INTEGER                   :: EventNum                                 ! The event number (index into EventID() array)
      REAL(ReKi)                :: TStart                                   ! The time at which to add this event
      REAL(ReKi)                :: delt                                     ! The delta time before the event begins (for interpolation in AeroDyn)
      LOGICAL                   :: Connect2Prev = .FALSE.                   ! Whether this event is connected to the next, otherwise there is space between them
      TYPE(Event), POINTER      :: Next         => NULL()                   ! The next event to add
   END TYPE
      
   TYPE :: CohStr_OutputType
      REAL(ReKi)                   :: CTKE                                     ! Maximum predicted Coherent Turbulent Kenetic Energy at the center of the billow
      REAL(ReKi)                   :: lambda                                   ! The expected value of interarrival times for the Poisson process
      INTEGER(IntKi)               :: NumCTEvents                              ! Number of events to be inserted into the .cts file
      INTEGER(IntKi)               :: NumCTEvents_separate                     ! Number of separate events inserted into the .cts file (# events with .Connect2Prev = .false.) 
      REAL(ReKi)                   :: ExpectedTime                             ! Amount of time the coherent structures should take
      REAL(ReKi)                   :: EventTimeSum                             ! Amount of time the coherent structure takes
      REAL(ReKi)                   :: EventTimeStep                            ! The average length of timesteps in output events   
      
      REAL(ReKi)                   :: Zbottom                                  ! The height of the lowest point on the grid (before tower points are added), equal to Z(1)
      REAL(ReKi)                   :: ScaleWid                                 ! Scaling width for LE coherent turbulence (RotDiam in AeroDyn FD_Wind)
      REAL(ReKi)                   :: ScaleVel                                 ! Scaling velocity for LE coherent turbulence, U0.  2*U0 is the difference in wind speed between the top and bottom of the wave.   
            
      REAL(ReKi)                   :: Uwave                                    ! Wind speed at center of the k-h billow (wave)
      REAL(ReKi)                   :: Wsig                                     ! Standard deviation of the w-component wind speed
      
      
      INTEGER(IntKi)               :: NumCTt                                   ! Number of data points to be printed in the output coherent event timestep file
            
      TYPE (Event), POINTER        :: PtrHead      => NULL()                   ! Pointer to the first event
      TYPE (Event), POINTER        :: PtrTail      => NULL()                   ! Pointer to the last event   
   
   END TYPE CohStr_OutputType   
   
        
      ! local type, used only in two subroutines here
   TYPE :: CohStr_EventType
      
      REAL(ReKi)                   :: Ym_max                                   ! The nondimensional lateral width of the coherent turbulence dataset
      REAL(ReKi)                   :: Zm_max                                   ! The nondimensional vertical height of the coherent turbulence dataset      
      
      REAL(ReKi),     ALLOCATABLE  :: pkCTKE     (:)                           ! Array containing the peak CTKE of each coherent event
      REAL(ReKi),     ALLOCATABLE  :: EventLen   (:)                           ! The length of each event stored in EventStart() (non-dimensional time)
      INTEGER(IntKi), ALLOCATABLE  :: EventID    (:)                           ! The timestep where the event starts, which determines the name of the event file
      INTEGER(IntKi), ALLOCATABLE  :: EventTS    (:)                           ! The length of each event stored in EventStart() (number of timesteps)
      INTEGER(IntKi)               :: NumEvents                                ! Number of events in the event data file (length of the Event* arrays)      
      
   END TYPE CohStr_EventType

   
   
CONTAINS

!=======================================================================
SUBROUTINE CohStr_ReadEventFile( p_CohStr, y_CohStr, e_CohStr, TSclFact, ErrStat, ErrMsg )

      ! This subroutine reads the events definitions from the event data file


   IMPLICIT                NONE


      ! Passed Variables
      
TYPE(CohStr_ParameterType), INTENT(IN   ) :: p_CohStr
TYPE(CohStr_OutputType),    INTENT(INOUT) :: y_CohStr
TYPE(CohStr_EventType),     INTENT(  OUT) :: e_CohStr

!REAL(ReKi),                      INTENT(IN)    :: CTKE           ! Predicted maximum CTKE
!REAL(ReKi),                      INTENT(INOUT) :: ScaleVel       ! The shear we're scaling for
!REAL(ReKi),                      INTENT(IN)    :: ScaleWid       ! The height of the wave we're scaling with
REAL(ReKi),                      INTENT(  OUT) :: TsclFact       ! Scale factor for time (h/U0) in coherent turbulence events
INTEGER(IntKi),                  intent(  out) :: ErrStat        ! Error level
CHARACTER(*),                    intent(  out) :: ErrMsg         ! Message describing error

      ! Local variables
REAL(ReKi)              :: MaxEvtCTKE        ! The maximum CTKE in the dataset of events

INTEGER                 :: I              ! DO loop counter
INTEGER                 :: IOS            ! I/O Status
INTEGER                 :: Un             ! I/O Unit

INTEGER(IntKi)                                 :: ErrStat2                         ! Error level (local)
CHARACTER(MaxMsgLen)                           :: ErrMsg2                          ! Message describing error (local)


   ErrStat = ErrID_None
   ErrMsg  = ""
   
   MaxEvtCTKE = 0.0  ! initialize the MAX variable

   CALL GetNewUnit( Un, ErrStat2, ErrMsg2 )
   
   CALL OpenFInpFile ( Un,  p_CohStr%CTEventFile, ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CohStr_ReadEventFile')
      IF (ErrStat >= AbortErrLev) RETURN
   

         ! Read the nondimensional lateral width of the dataset, Ym_max

   CALL ReadVar( Un, p_CohStr%CTEventFile, e_CohStr%Ym_max, "Ym_max", "Nondimensional lateral dataset width", ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CohStr_ReadEventFile')

         ! Read the nondimensional vertical height of the dataset, Zm_max

   CALL ReadVar( Un, p_CohStr%CTEventFile, e_CohStr%Zm_max, "Zm_max", "Nondimensional vertical dataset height", ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CohStr_ReadEventFile')


         ! Read the rest of the header

   CALL ReadVar( Un, p_CohStr%CTEventFile, e_CohStr%NumEvents, "NumEvents", "the number of coherent structures.", ErrStat2, ErrMsg2)
      CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CohStr_ReadEventFile')
      
      
   IF (ErrStat >= AbortErrLev) THEN
      CLOSE(Un)
      RETURN
   END IF
      

   IF ( e_CohStr%NumEvents > 0 ) THEN


      CALL AllocAry( e_CohStr%EventID,  e_CohStr%NumEvents , 'EventID',  ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'CohStr_ReadEventFile')
      CALL AllocAry( e_CohStr%EventTS,  e_CohStr%NumEvents , 'EventTS',  ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'CohStr_ReadEventFile')
      CALL AllocAry( e_CohStr%EventLen, e_CohStr%NumEvents , 'EventLen', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'CohStr_ReadEventFile')
      CALL AllocAry( e_CohStr%pkCTKE,   e_CohStr%NumEvents , 'pkCTKE',   ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'CohStr_ReadEventFile')

      IF (ErrStat >= AbortErrLev) THEN
         CLOSE(Un)
         RETURN
      END IF
      
            ! Read the last header lines

      CALL ReadCom( Un, p_CohStr%CTEventFile, 'the fourth header line', ErrStat2, ErrMsg2) ! A blank line
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CohStr_ReadEventFile')  
      CALL ReadCom( Un, p_CohStr%CTEventFile, 'the fifth header line', ErrStat2, ErrMsg2) ! The column heading lines
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CohStr_ReadEventFile')   


            ! Read the event definitions and scale times by TScale

      DO I=1,e_CohStr%NumEvents

         READ ( Un, *, IOSTAT=IOS )  e_CohStr%EventID(I),  e_CohStr%EventTS(I), e_CohStr%EventLen(I), e_CohStr%pkCTKE(I)

         IF ( IOS /= 0 )  THEN
            CALL SetErrStat(ErrID_Fatal, 'Error reading event '//TRIM( Int2LStr( I ) )//' from the coherent event data file.', ErrStat, ErrMsg, 'CohStr_ReadEventFile')
            CLOSE(UN)
            RETURN
         ENDIF
         MaxEvtCTKE = MAX( MaxEvtCTKE, e_CohStr%pkCTKE(I) )

      ENDDO

      IF ( MaxEvtCTKE > 0.0 ) THEN
            y_CohStr%ScaleVel = MAX( y_CohStr%ScaleVel, SQRT( y_CohStr%CTKE / MaxEvtCTKE ) )
            ! Calculate the Velocity Scale Factor, based on the requested maximum CTKE
      ENDIF

         ! Calculate the TimeScaleFactor, based on the Zm_max in the Events file.

      TSclFact = y_CohStr%ScaleWid / (y_CohStr%ScaleVel * e_CohStr%Zm_max)

         ! Scale the time based on TSclFact

      DO I=1,e_CohStr%NumEvents
         e_CohStr%EventLen(I) = e_CohStr%EventLen(I)*TSclFact
      ENDDO

   ELSE

      TSclFact = y_CohStr%ScaleWid / (y_CohStr%ScaleVel * e_CohStr%Zm_max)

   ENDIF  ! FileNum > 0

   CLOSE ( Un )  

END SUBROUTINE CohStr_ReadEventFile
!=======================================================================
SUBROUTINE CohStr_CalcEvents( p, e_CohStr, Height, OtherSt_RandNum, y_cohStr, ErrStat, ErrMsg )

      ! This subroutine calculates what events to use and when to use them.
      ! It computes the number of timesteps in the file, NumCTt.

   IMPLICIT                    NONE

      ! passed variables
   TYPE(TurbSim_ParameterType),     INTENT(IN)     :: P
   TYPE(CohStr_EventType)         , INTENT(IN)     :: e_CohStr            ! event parameters for coherent structures 
   REAL(ReKi),                      INTENT(IN)     :: Height              ! Height for expected length PDF equation
   TYPE(RandNum_OtherStateType),    INTENT(INOUT)  :: OtherSt_RandNum     ! other states for random numbers (next seed, etc)
   TYPE(CohStr_OutputType),         INTENT(INOUT)  :: y_CohStr
   INTEGER(IntKi),                  INTENT(  OUT)  :: ErrStat ! Error level
   CHARACTER(*),                    INTENT(  OUT)  :: ErrMsg  ! Message describing error

      ! local variables
REAL(ReKi)                  :: iA                  ! Variable used to calculate IAT
REAL(ReKi)                  :: iB                  ! Variable used to calculate IAT
REAL(ReKi)                  :: iC                  ! Variable used to calculate IAT
REAL(ReKi)                  :: rn                  ! random number
REAL(ReKi)                  :: TEnd                ! End time for the current event
REAL(ReKi)                  :: TStartNext   = 0.0  ! temporary start time for next event
REAL(ReKi)                  :: MaxCTKE             ! Maximum CTKE of events we've picked

INTEGER                     :: ErrStat2            ! temp error status
INTEGER                     :: NewEvent            ! event number of the new event
INTEGER                     :: NumCompared         ! Number of events we've compared

LOGICAL(1)                  :: Inserted            ! Whether an event was inserted here

TYPE(Event), POINTER        :: PtrCurr  => NULL()  ! Pointer to the current event in the list
TYPE(Event), POINTER        :: PtrNew   => NULL()  ! A new event to be inserted into the list


   ErrStat = ErrID_None
   ErrMsg = ""

      ! Compute the mean interarrival time and the expected length of events

   SELECT CASE ( p%met%TurbModel_ID )

      CASE ( SpecModel_NWTCUP, SpecModel_NONE, SpecModel_USRVKM )
         y_CohStr%lambda = -0.000904*p%met%Rich_No + 0.000562*y_CohStr%Uwave + 0.001389
         y_CohStr%lambda = 1.0 / y_CohStr%lambda

         IF ( p%met%TurbModel_ID == SpecModel_NONE ) THEN
            y_CohStr%ExpectedTime = 600.0
         ELSE
            CALL RndModLogNorm( p%RNG, OtherSt_RandNum, y_CohStr%ExpectedTime, Height )
         ENDIF

      CASE ( SpecModel_GP_LLJ, SpecModel_SMOOTH, SpecModel_TIDAL, SpecModel_RIVER) ! HYDRO: added 'TIDAL' and 'RIVER' to the spectral models that get handled this way.
         iA     =        0.001797800 + (7.17399E-10)*Height**3.021144723
         iB     =  EXP(-10.590340100 - (4.92440E-05)*Height**2.5)
         iC     = SQRT(  3.655013599 + (8.91203E-06)*Height**3  )
         y_CohStr%lambda = iA + iB*MIN( (y_CohStr%Uwave**iC), HUGE(iC) )  ! lambda = iA + iB*(WindSpeed**iC)
         y_CohStr%lambda = 1.0 / y_CohStr%lambda

         CALL RndTcohLLJ( p%RNG, OtherSt_RandNum, y_CohStr%ExpectedTime, Height )

      CASE ( SpecModel_WF_UPW )
        y_CohStr%lambda = 0.000529*y_CohStr%Uwave + 0.000365*p%met%Rich_No - 0.000596
        y_CohStr%lambda = 1.0 / y_CohStr%lambda

         CALL RndTcoh_WF( p%RNG, OtherSt_RandNum, y_CohStr%ExpectedTime, SpecModel_WF_UPW )

      CASE ( SpecModel_WF_07D )
         y_CohStr%lambda = 0.000813*y_CohStr%Uwave - 0.002642*p%met%Rich_No + 0.002676
         y_CohStr%lambda = 1.0 / y_CohStr%lambda

         CALL RndTcoh_WF( p%RNG, OtherSt_RandNum, y_CohStr%ExpectedTime, SpecModel_WF_07D )

      CASE ( SpecModel_WF_14D )
         y_CohStr%lambda = 0.001003*y_CohStr%Uwave - 0.00254*p%met%Rich_No - 0.000984
         y_CohStr%lambda = 1.0 / y_CohStr%lambda

         CALL RndTcoh_WF( p%RNG, OtherSt_RandNum, y_CohStr%ExpectedTime, SpecModel_WF_14D )

      CASE DEFAULT
         !This should not happen

   END SELECT

   y_CohStr%ExpectedTime = y_CohStr%ExpectedTime * ( p%grid%UsableTime - p%CohStr%CTStartTime ) / 600.0_ReKi  ! Scale for use with the amount of time we've been given


!BONNIE: PERHAPS WE SHOULD JUST PUT IN A CHECK THAT TURNS OFF THE COHERENT TIME STEP FILE IF THE
!        CTSTARTTIME IS LESS THAN THE USABLETIME... MAYBE WHEN WE'RE READING THE INPUT FILE...
y_CohStr%ExpectedTime = MAX( y_CohStr%ExpectedTime, 0.0_ReKi )  ! This occurs if CTStartTime = 0

      ! We start by adding events at random times

   y_CohStr%NumCTEvents = 0                                    ! Number of events = length of our linked list
   y_CohStr%NumCTt      = 0                                    ! Total number of time steps in the events we've picked
   MaxCTKE              = 0.0                                  ! Find the maximum CTKE for the events that we've selected

   y_CohStr%EventTimeSum = 0.0 

   CALL RndExp(p%RNG, OtherSt_RandNum, rn, y_CohStr%lambda)                            ! Assume the last event ended at time zero

   TStartNext = rn / 2.0

   IF ( p%met%KHtest ) THEN
      y_CohStr%ExpectedTime = p%grid%UsableTime     / 2                 ! When testing, add coherent events for half of the record
      TStartNext            = y_CohStr%ExpectedTime / 2                 ! When testing, start about a quarter of the way into the record
   ENDIF

   IF ( TStartNext < p%CohStr%CTStartTime ) THEN
      TStartNext = TStartNext + p%CohStr%CTStartTime           ! Make sure the events start after time specified by CTStartTime
   ENDIF

   IF ( TStartNext > 0 ) y_CohStr%NumCTt = y_CohStr%NumCTt + 1          ! Add a point before the first event

   DO WHILE ( TStartNext < p%grid%UsableTime .AND. y_CohStr%EventTimeSum < y_CohStr%ExpectedTime )

      CALL RndUnif( p%RNG, OtherSt_RandNum, rn )

      NewEvent = INT( rn*( e_CohStr%NumEvents - 1 ) ) + 1
      NewEvent = MAX( 1, MIN( NewEvent, e_CohStr%NumEvents ) ) ! take care of possible rounding issues....


      IF ( .NOT. ASSOCIATED ( y_CohStr%PtrHead ) ) THEN

         ALLOCATE ( y_CohStr%PtrHead, STAT=ErrStat2 )             ! The pointer %Next is nullified in allocation

         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for new event.' , ErrStat, ErrMsg, 'CohStr_CalcEvents')
            RETURN
         ENDIF

         y_CohStr%PtrTail => y_CohStr%PtrHead

      ELSE

         ALLOCATE ( y_CohStr%PtrTail%Next, STAT=ErrStat2 )     ! The pointer PtrTail%Next%Next is nullified in allocation

         IF ( ErrStat2 /= 0 ) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for new event.' , ErrStat, ErrMsg, 'CohStr_CalcEvents')
            RETURN
         ENDIF

         y_CohStr%PtrTail => y_CohStr%PtrTail%Next                   ! Move the pointer to point to the last record in the list

      ENDIF

      y_CohStr%PtrTail%EventNum     = NewEvent
      y_CohStr%PtrTail%TStart       = TStartNext
      y_CohStr%PtrTail%delt         = e_CohStr%EventLen( NewEvent ) / e_CohStr%EventTS( NewEvent )          ! the average delta time in the event
      y_CohStr%PtrTail%Connect2Prev = .FALSE.

      MaxCTKE                       = MAX( MaxCTKE, e_CohStr%pkCTKE( NewEvent ) )
      y_CohStr%NumCTEvents          = y_CohStr%NumCTEvents + 1

      TEnd = TStartNext + e_CohStr%EventLen( NewEvent )


      IF ( p%met%KHtest ) THEN
         TStartNext   = p%grid%UsableTime + TStartNext !TEnd + PtrTail%delt ! Add the events right after each other
      ELSE

         DO WHILE ( TStartNext <= TEnd )

            CALL RndExp(p%RNG, OtherSt_RandNum, rn, y_CohStr%lambda)    ! compute the interarrival time
            TStartNext        = TStartNext + rn !+ EventLen( NewEvent )

         ENDDO

      ENDIF


      IF ( (TStartNext - TEnd) > y_CohStr%PtrTail%delt ) THEN
         y_CohStr%NumCTt = y_CohStr%NumCTt + e_CohStr%EventTS( NewEvent ) + 2                                  ! add a zero-line (essentially a break between events)
      ELSE
         y_CohStr%NumCTt = y_CohStr%NumCTt + e_CohStr%EventTS( NewEvent ) + 1
      ENDIF

      y_CohStr%EventTimeSum     = y_CohStr%EventTimeSum + e_CohStr%EventLen( NewEvent )

   ENDDO

   y_CohStr%NumCTEvents_separate = y_CohStr%NumCTEvents

      ! Next, we start concatenating events until there is no space or we exceed the expected time

   IF ( p%met%TurbModel_ID /= SpecModel_NONE ) THEN

      NumCompared = 0

      DO WHILE ( y_CohStr%EventTimeSum < y_CohStr%ExpectedTime .AND. NumCompared < y_CohStr%NumCTEvents )

         CALL RndUnif( p%RNG, OtherSt_RandNum, rn )

         NewEvent = INT( rn*( e_CohStr%NumEvents - 1.0 ) ) + 1
         NewEvent = MAX( 1, MIN( NewEvent, e_CohStr%NumEvents ) )    ! take care of possible rounding issues....

         NumCompared = 0
         Inserted    = .FALSE.

         DO WHILE ( NumCompared < y_CohStr%NumCTEvents .AND. .NOT. Inserted )

            IF ( .NOT. ASSOCIATED ( PtrCurr ) ) THEN        ! Wrap around to the beginning of the list
               PtrCurr => y_CohStr%PtrHead
            ENDIF


               ! See if the NewEvent fits between the end of event pointed to by PtrCurr and the
               ! beginning of the event pointed to by PtrCurr%Next

            IF ( ASSOCIATED( PtrCurr%Next ) ) THEN
               TStartNext = PtrCurr%Next%TStart
            ELSE !We're starting after the last event in the record
               TStartNext = p%grid%UsableTime + 0.5 * e_CohStr%EventLen( NewEvent )  ! We can go a little beyond the end...
            ENDIF

            IF ( TStartNext - (PtrCurr%TStart + e_CohStr%EventLen( PtrCurr%EventNum ) + PtrCurr%delt) > e_CohStr%EventLen( NewEvent ) ) THEN

               Inserted = .TRUE.

               ALLOCATE ( PtrNew, STAT=ErrStat2 )           ! The pointer %Next is nullified in allocation

               IF ( ErrStat2 /= 0 ) THEN
                  CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for new event.' , ErrStat, ErrMsg, 'CohStr_CalcEvents')
               ENDIF

               PtrNew%EventNum      = NewEvent
               PtrNew%TStart        = PtrCurr%TStart + e_CohStr%EventLen( PtrCurr%EventNum )
               PtrNew%delt          = e_CohStr%EventLen( NewEvent ) / e_CohStr%EventTS( NewEvent )          ! the average delta time in the event
               PtrNew%Connect2Prev  = .TRUE.

               PtrNew%Next  => PtrCurr%Next
               PtrCurr%Next => PtrNew
               PtrCurr      => PtrCurr%Next    ! Let's try to add the next event after the other events

               MaxCTKE                       = MAX( MaxCTKE, e_CohStr%pkCTKE( NewEvent ) )
               y_CohStr%NumCTEvents          = y_CohStr%NumCTEvents + 1
               y_CohStr%NumCTt               = y_CohStr%NumCTt + e_CohStr%EventTS( NewEvent )  ! there is no break between events
                                   !(we may have one too many NumCTt here, so we'll deal with it when we write the file later)
               y_CohStr%EventTimeSum         = y_CohStr%EventTimeSum + e_CohStr%EventLen( NewEvent )


            ELSE

               NumCompared = NumCompared + 1

            ENDIF

            PtrCurr => PtrCurr%Next

         ENDDO ! WHILE (NumCompared < NumCTEvents .AND. .NOT. Inserted)

      ENDDO ! WHILE (EventTimeSum < ExpectedTime .AND. NumCompared < NumCTEvents)

   ENDIF ! SpecModel /= SpecModel_NONE

   IF ( y_CohStr%NumCTt > 0 ) THEN
      y_CohStr%EventTimeStep = y_CohStr%EventTimeSum / y_CohStr%NumCTt                                          ! Average timestep of coherent event data
   ELSE
      y_CohStr%EventTimeStep = 0.0
   ENDIF


END SUBROUTINE CohStr_CalcEvents
!=======================================================================
!> This subroutine writes the coherent events CTS file
SUBROUTINE CohStr_WriteCTS(p, WSig, OtherSt_RandNum, ErrStat, ErrMsg)

   TYPE(TurbSim_ParameterType),     INTENT(IN   )  :: p                   ! parameters for TurbSim (out only b/c it doesn't generate file for certain cases...)
   TYPE(RandNum_OtherStateType),    INTENT(INOUT)  :: OtherSt_RandNum     ! other states for random numbers (next seed, etc)
   REAL(ReKi),                      intent(in   )  :: WSig                ! Standard deviation of the vertical component

   INTEGER(IntKi),                  intent(  out)  :: ErrStat             ! Error level
   CHARACTER(*),                    intent(  out)  :: ErrMsg              ! Message describing error


      ! local variables
   TYPE(CohStr_OutputType)                         :: y_CohStr
   type(CohStr_EventType)                          :: e_CohStr            ! coherent structure events
   REAL(ReKi)                                      :: TmpVel              ! A temporary variable holding a velocity
   REAL(ReKi)                                      :: TmpRndNum           ! A temporary variable holding a random variate
   REAL(ReKi)                                      :: TsclFact            ! Scale factor for time (h/U0) in coherent turbulence events
   
   INTEGER(IntKi)                                  :: ErrStat2
   CHARACTER(MaxMsgLen)                            :: ErrMsg2
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   y_CohStr%WSig=WSig
   

   CALL WrScr ( ' Generating coherent turbulent time step file "'//TRIM( p%RootName )//'.cts"' )
   
   
   y_CohStr%ScaleWid = p%grid%RotorDiameter * p%CohStr%DistScl           !  This is the scaled height of the coherent event data set
   y_CohStr%Zbottom  = p%grid%HubHt - p%CohStr%CTLz*y_CohStr%ScaleWid    !  This is the height of the bottom of the wave in the scaled/shifted coherent event data set
   
   CALL getVelocity(p, p%UHub,p%grid%HubHt,y_CohStr%Zbottom + 0.5_ReKi*y_CohStr%ScaleWid, y_CohStr%Uwave, ErrStat2, ErrMsg2)                 ! y_CohStr%Uwave =WindSpeed at center of wave
      CALL SetErrStat(  ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CohStr_WriteCTS')
   
   !-------------------------
   ! compute ScaleVel:
   !-------------------------

   IF ( p%met%KHtest ) THEN      
         ! for LES test case....
      y_CohStr%ScaleVel = y_CohStr%ScaleWid * KHT_LES_dT /  KHT_LES_Zm    
      y_CohStr%ScaleVel = 50 * y_CohStr%ScaleVel                  ! We want 25 hz bandwidth so multiply by 50
   ELSE
      
      CALL getVelocity(p, p%UHub,p%grid%HubHt,y_CohStr%Zbottom,                   TmpVel,            ErrStat2, ErrMsg2)    ! Velocity at bottom of billow
         CALL SetErrStat(  ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CohStr_WriteCTS') 
      CALL getVelocity(p, p%UHub,p%grid%HubHt,y_CohStr%Zbottom+y_CohStr%ScaleWid, y_CohStr%ScaleVel, ErrStat2, ErrMsg2)    ! Velocity at the top of the billow
         CALL SetErrStat(  ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CohStr_WriteCTS')
            
      y_CohStr%ScaleVel =  y_CohStr%ScaleVel - TmpVel   ! Shear across the wave
      y_CohStr%ScaleVel = 0.5 * y_CohStr%ScaleVel       ! U0 is half the difference between the top and bottom of the billow
      
      
         ! If the coherent structures do not cover the whole disk, increase the shear

      IF ( p%CohStr%DistScl < 1.0 ) THEN ! Increase the shear by up to two when the wave is half the size of the disk...
         CALL RndUnif( p%RNG, OtherSt_RandNum, TmpRndNum ) !returns TmpRndNum, a random variate
         y_CohStr%ScaleVel = y_CohStr%ScaleVel * ( 1.0 + TmpRndNum * (1 - p%CohStr%DistScl) / p%CohStr%DistScl )
      ENDIF

         !Apply a scaling factor to account for short inter-arrival times getting wiped out due to long events

      y_CohStr%ScaleVel =  y_CohStr%ScaleVel*( 1.0 + 323.1429 * EXP( -MAX(y_CohStr%Uwave,10.0_ReKi) / 2.16617 ) )
            
   ENDIF

   IF (y_CohStr%ScaleVel < 0. ) THEN
      CALL SetErrStat( ErrID_Warn, 'A coherent turbulence time step file cannot be generated with negative shear.', ErrStat, ErrMsg, 'CohStr_WriteCTS')
      !p%WrFile(FileExt_CTS) = .FALSE.
      RETURN
   ENDIF
      

   !-------------------------
   ! compute maximum predicted CTKE:
   !-------------------------
         
   SELECT CASE ( p%met%TurbModel_ID )
         
      CASE ( SpecModel_NWTCUP,  SpecModel_NONE, SpecModel_USRVKM )
            
         IF (p%met%KHtest) THEN
            y_CohStr%CTKE = 30.0 !Scale for large coherence
            CALL RndNWTCpkCTKE( p%RNG, OtherSt_RandNum, y_CohStr%CTKE )
         ELSE    
               
               ! Increase the Scaling Velocity for computing U,V,W in AeroDyn
               ! These numbers are based on LIST/ART data (58m-level sonic anemometer)

            y_CohStr%CTKE =  0.616055*p%met%Rich_No - 0.242143*y_CohStr%Uwave + 23.921801*y_CohStr%WSig - 11.082978
            
               ! Add up to +/- 10% or +/- 6 m^2/s^2 (uniform distribution)
            CALL RndUnif( p%RNG, OtherSt_RandNum, TmpRndNum )
            y_CohStr%CTKE = MAX( y_CohStr%CTKE + (2.0_ReKi * TmpRndNum - 1.0_ReKi) * 6.0_ReKi, 0.0_ReKi )

            IF ( y_CohStr%CTKE > 0.0 ) THEN
               IF ( y_CohStr%CTKE > 20.0)  THEN    ! Correct with residual
                  y_CohStr%CTKE = y_CohStr%CTKE + ( 0.11749127 * (y_CohStr%CTKE**1.369025) - 7.5976449 )
               ENDIF

               IF ( y_CohStr%CTKE >= 30.0 .AND. p%met%Rich_No >= 0.0 .AND. p%met%Rich_No <= 0.05 ) THEN
                  CALL RndNWTCpkCTKE( p%RNG, OtherSt_RandNum, y_CohStr%CTKE )
               ENDIF
            ENDIF
                  
         ENDIF !p%met%KHtest
               
      CASE ( SpecModel_GP_LLJ, SpecModel_SMOOTH, SpecModel_TIDAL, SpecModel_RIVER )          

         y_CohStr%CTKE = pkCTKE_LLJ( p, OtherSt_RandNum, y_CohStr%Zbottom+0.5_ReKi*y_CohStr%ScaleWid, p%met%ZL, p%met%UStar )
               
      CASE ( SpecModel_WF_UPW )
         y_CohStr%CTKE = -2.964523*p%met%Rich_No - 0.207382*y_CohStr%Uwave + 25.640037*y_CohStr%WSig - 10.832925
               
      CASE ( SpecModel_WF_07D )
         y_CohStr%CTKE = 9.276618*p%met%Rich_No + 6.557176*p%met%Ustar + 3.779539*y_CohStr%WSig - 0.106633

         IF ( (p%met%Rich_No > -0.025) .AND. (p%met%Rich_No < 0.05) .AND. (p%met%Ustar > 1.0) .AND. (p%met%Ustar < 1.56) ) THEN
            CALL RndpkCTKE_WFTA( p%RNG, OtherSt_RandNum, TmpRndNum )  ! Add a random residual
            y_CohStr%CTKE = y_CohStr%CTKE + TmpRndNum
         ENDIF
               
               
      CASE ( SpecModel_WF_14D )
         y_CohStr%CTKE = 1.667367*p%met%Rich_No - 0.003063*y_CohStr%Uwave + 19.653682*y_CohStr%WSig - 11.808237
               
      CASE DEFAULT   ! This case should not happen   
         CALL SetErrStat( ErrID_Fatal, 'Invalid turbulence model in coherent structure analysis.', ErrStat, ErrMsg, 'CohStr_WriteCTS')
         CALL Cleanup()
         RETURN               
   END SELECT                    

   y_CohStr%CTKE   = MAX( y_CohStr%CTKE, 1.0_ReKi )     ! make sure CTKE is not negative and, so that we don't divide by zero in ReadEventFile, set it to some arbitrary low number   
   
   !-------------------------
   ! Read and allocate coherent event start times and lengths, calculate TSclFact:
   !-------------------------   
   CALL CohStr_ReadEventFile( p%CohStr, y_CohStr, e_CohStr, TSclFact, ErrStat2, ErrMsg2 ) !y_CohStr%%ScaleWid, y_CohStr%ScaleVel
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CohStr_WriteCTS')
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF
   
   CALL CohStr_CalcEvents( p, e_CohStr, y_CohStr%Zbottom+0.5_ReKi*y_CohStr%ScaleWid, OtherSt_RandNum, y_cohStr, ErrStat2, ErrMsg2) 
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CohStr_WriteCTS')
      IF (ErrStat >= AbortErrLev) THEN
         CALL Cleanup()
         RETURN
      END IF

   !-------------------------
   ! Write the file:
   !-------------------------
   
   CALL CohStr_WriteEvents ( p%RootName, p%CohStr, e_CohStr, y_CohStr, TSclFact, p%UHub, ErrStat2, ErrMsg2) 
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CohStr_WriteCTS')

   
   !-------------------------
   ! Write some summary information:
   !-------------------------
   IF (p%met%KHtest) THEN
      WRITE ( p%US,'(/)' )
   ELSE
      WRITE ( p%US,'(//A,F8.3," seconds")' ) 'Average expected time between events = ',y_CohStr%lambda
   ENDIF

   WRITE ( p%US, '(A,I8)'   )            'Number of coherent events            = ', y_CohStr%NumCTEvents_separate
   WRITE ( p%US, '(A,F8.3," seconds")')  'Predicted length of coherent events  = ', y_CohStr%ExpectedTime
   WRITE ( p%US, '(A,F8.3," seconds")')  'Length of coherent events            = ', y_CohStr%EventTimeSum
   WRITE ( p%US, '(A,F8.3," (m/s)^2")')  'Maximum predicted event CTKE         = ', y_CohStr%CTKE
IF ( y_CohStr%EventTimeStep > 0.0_ReKi ) THEN
   WRITE ( p%US, '(A,F8.3," Hz")'     )  'Nyquist frequency of coherent events = ',  0.5_ReKi / y_CohStr%EventTimeStep
ENDIF
      
      
   !-------------------------
   ! Deallocate the coherent event arrays.
   !-------------------------
   CALL Cleanup()
   RETURN
   
CONTAINS
   SUBROUTINE Cleanup()
      IF ( ALLOCATED( e_CohStr%EventID   ) )  DEALLOCATE( e_CohStr%EventID   )
      IF ( ALLOCATED( e_CohStr%EventTS   ) )  DEALLOCATE( e_CohStr%EventTS   )
      IF ( ALLOCATED( e_CohStr%EventLen  ) )  DEALLOCATE( e_CohStr%EventLen  )
      IF ( ALLOCATED( e_CohStr%pkCTKE    ) )  DEALLOCATE( e_CohStr%pkCTKE    )
   END SUBROUTINE  Cleanup  
         
END SUBROUTINE CohStr_WriteCTS
!=======================================================================
FUNCTION pkCTKE_LLJ(p, OtherSt_RandNum, Ht, ZL, UStar)
   
   IMPLICIT                 NONE

   TYPE(TurbSim_ParameterType),     INTENT(IN   ):: p                   ! parameters
   TYPE(RandNum_OtherStateType),    INTENT(INOUT):: OtherSt_RandNum     ! other states for random numbers (next seed, etc)
   REAL(ReKi),                      INTENT(IN)   :: Ht                  ! The height at the billow center
   REAL(ReKi),                      INTENT(IN)   :: ZL                  ! The height at the billow center
   REAL(ReKi),                      INTENT(IN)   :: Ustar               ! The height at the billow center
   
   ! local variables

   REAL(ReKi)               :: A                                        ! A constant/offset term in the pkCTKE calculation
   REAL(ReKi)               :: A_uSt                                    ! The scaling term for Ustar
   REAL(ReKi)               :: A_zL                                     ! The scaling term for z/L
   REAL(ReKi)               :: pkCTKE_LLJ                               ! The max CTKE expected for LLJ coh structures
   REAL(ReKi)               :: rndCTKE                                  ! The random residual

   REAL(ReKi), PARAMETER    :: RndParms(5) = (/0.252510525, -0.67391279, 2.374794977, 1.920555797, -0.93417558/) ! parameters for the Pearson IV random residual
   REAL(ReKi), PARAMETER    :: z_Ary(4)    = (/54., 67., 85., 116./)    ! Aneomoeter heights

   INTEGER                  :: Zindx_mn (1)

   

   Zindx_mn = MINLOC( ABS(z_Ary-Ht) )

   SELECT CASE ( Zindx_mn(1) )
      CASE ( 1 )  ! 54 m
         A     = -0.051
         A_zL  = -0.0384
         A_uSt =  9.9710

      CASE ( 2 )  ! 67 m
         A     = -0.054
         A_zL  = -0.1330
         A_uSt = 10.2460

      CASE ( 3 )  ! 85 m
         A     = -0.062
         A_zL  = -0.1320
         A_uSt = 10.1660

      CASE ( 4 )  !116 m
         A     = -0.092
         A_zL  = -0.3330
         A_uSt = 10.7640

      !CASE DEFAULT !This should not occur
      !   ErrStat = ErrID_Fatal
      !   ErrMsg  = 'Error in pkCTKE_LLJ():: Height index is invalid.' 
   END SELECT

   CALL RndPearsonIV( p%RNG, OtherSt_RandNum, rndCTKE, RndParms, (/ -10.0_ReKi, 17.5_ReKi /) )

   pkCTKE_LLJ = MAX(0.0_ReKi, A + A_uSt*UStar + A_zL*ZL + rndCTKE)

END FUNCTION pkCTKE_LLJ
!=======================================================================
SUBROUTINE CohStr_WriteEvents( RootName, p_CohStr, e_CohStr, y_CohStr, TScale, UHub, ErrStat, ErrMsg )

    ! This subroutine writes the events as calculated in CalcEvents.

   IMPLICIT                NONE

         ! Passed Variables
   TYPE(CohStr_ParameterType)     , INTENT(IN   )  :: p_CohStr            ! parameters for coherent structures
   TYPE(CohStr_EventType)         , INTENT(IN   )  :: e_CohStr            ! parameters for coherent structure events
   TYPE(CohStr_OutputType),         INTENT(INOUT)  :: y_CohStr

   REAL(ReKi),                      INTENT(IN)     :: TScale              ! Time scaling factor
   REAL(ReKi),                      INTENT(IN)     :: UHub                ! Mean wind speed at hub height (advection speed)
   CHARACTER(*),                    INTENT(IN)     :: RootName

   INTEGER(IntKi),                  intent(  out)  :: ErrStat             ! Error level
   CHARACTER(*),                    intent(  out)  :: ErrMsg              ! Message describing error

 
      ! Local Variables

   REAL(ReKi)              :: CurrentTime = 0.0      ! the current time (in seconds)
   REAL(ReKi)              :: CTTime                 ! Time from beginning of event file
   REAL(ReKi)              :: deltaTime = 0.0        ! difference between two time steps in the event files

   INTEGER                 :: FileNum                ! File Number in the event file
   INTEGER                 :: IE                     ! Loop counter for event number
   INTEGER                 :: IT                     ! Loop counter for time step

   INTEGER(IntKi)          :: UnIn                   ! I/O Unit for input file
   INTEGER(IntKi)          :: UnOut                  ! I/O Unit for output file

   INTEGER(IntKi)          :: ErrStat2                        ! Error level (local)
  CHARACTER(MaxMsgLen)     :: ErrMsg2                         ! Message describing error (local)


   CHARACTER(200)          :: InpFile                ! Name of the input file
   TYPE (Event), POINTER   :: PtrCurr  => NULL()     ! Pointer to the current event
   TYPE (Event), POINTER   :: PtrPrev  => NULL()     ! Pointer to the previous event (for deallocation purposes)


   ErrStat = ErrID_None
   ErrMsg  = ""
   
   UnOut = -1
   CALL GetNewUnit( UnOut, ErrStat2, ErrMsg2 )      
   CALL OpenFOutFile ( UnOut, TRIM( RootName )//'.cts', ErrStat2, ErrMsg2 ) 
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CohStr_WriteEvents' )

   UnIn = -1
   CALL GetNewUnit( UnIn, ErrStat2, ErrMsg2 )
      

      ! Write event data to the time step output file (opened at the beginnig)

   WRITE (UnOut, "( A14,   ' = FileType')")     p_CohStr%CTExt
   WRITE (UnOut, "( G14.7, ' = ScaleVel')")     y_CohStr%ScaleVel
   WRITE (UnOut, "( G14.7, ' = MHHWindSpeed')") UHub
   WRITE (UnOut, "( G14.7, ' = Ymax')")         y_CohStr%ScaleWid*e_CohStr%Ym_max/e_CohStr%Zm_max
   WRITE (UnOut, "( G14.7, ' = Zmax')")         y_CohStr%ScaleWid
   WRITE (UnOut, "( G14.7, ' = DistScl')")      p_CohStr%DistScl
   WRITE (UnOut, "( G14.7, ' = CTLy')")         p_CohStr%CTLy
   WRITE (UnOut, "( G14.7, ' = CTLz')")         p_CohStr%CTLz
   WRITE (UnOut, "( G14.7, ' = NumCTt')")       y_CohStr%NumCTt


   PtrCurr => y_CohStr%PtrHead


   DO IE = 1,y_CohStr%NumCTEvents

      IF ( .NOT. ASSOCIATED ( PtrCurr ) ) EXIT     ! This shouldn't be necessary, given the way we created the list


      IF ( .NOT. PtrCurr%Connect2Prev ) THEN

         IF ( CurrentTime < PtrCurr%TStart ) THEN

            WRITE ( UnOut, '(G14.7,1x,I5.5)') CurrentTime, 0     ! Print end of previous event

            y_CohStr%NumCTt = y_CohStr%NumCTt - 1  ! Let's make sure the right number of points have been written to the file.

            IF ( CurrentTime < PtrCurr%TStart - PtrCurr%delt ) THEN  !This assumes a ramp of 1 delta t for each structure....

               WRITE ( UnOut, '(G14.7,1x,I5.5)') MAX(PtrCurr%TStart - PtrCurr%delt, REAL(0.0, ReKi) ), 0
               y_CohStr%NumCTt = y_CohStr%NumCTt - 1

            ENDIF

         ENDIF

      ENDIF  ! NOT Connect2Prev


      WRITE ( InpFile, '(I5.5)' ) e_CohStr%EventID( PtrCurr%EventNum )
      InpFile = TRIM( p_CohStr%CTEventPath )//PathSep//'Event'//TRIM( InpFile)//'.dat'

      CALL OpenFInpFile( UnIn, InpFile, ErrStat2, ErrMsg2 ) 
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'CohStr_WriteEvents' )
         IF ( ErrStat >= AbortErrLev ) THEN
            CALL Cleanup()
            RETURN
         END IF
         

      DO IT = 1,e_CohStr%EventTS( PtrCurr%EventNum )

         READ  ( UnIn, *, IOSTAT=ErrStat2 ) FileNum, CTTime, deltaTime

         IF (ErrStat2 /= 0) THEN
            CALL SetErrStat( ErrID_Fatal, 'Error reading event file'//TRIM( InpFile ), ErrStat, ErrMsg, 'CohStr_WriteEvents')
            CALL Cleanup()
            RETURN
         ENDIF

         CurrentTime = PtrCurr%TStart + CTTime*TScale

         WRITE ( UnOut, '(G14.7,1x,I5.5)') CurrentTime, FileNum
         y_CohStr%NumCTt = y_CohStr%NumCTt - 1

      ENDDO    ! IT: Event timestep


      CLOSE ( UnIn )


         ! Add one (delta time) space between events

      CurrentTime = CurrentTime + deltaTime*TScale

      PtrPrev => PtrCurr
      PtrCurr => PtrCurr%Next
      
      DEALLOCATE ( PtrPrev, STAT=ErrStat2 )

   ENDDO !IE: number of events

   WRITE ( UnOut, '(G14.7,1x,I5.5)') CurrentTime, 0     !Add the last line
   y_CohStr%NumCTt = y_CohStr%NumCTt - 1

      ! Let's append zero lines at the end of the file if we haven't output NumCTt lines, yet.
      ! We've subtracted from NumCTt every time we wrote a line so now the number in NumCTt is
      ! how many lines short we are.

   IF ( deltaTime > 0 ) THEN
      deltaTime = deltaTime*TScale
   ELSE
      deltaTime = 0.5
   ENDIF

   DO IE = 1, y_CohStr%NumCTt  ! Write zeros at the end if we happened to insert an event that overwrote one of our zero lines
      CurrentTime = CurrentTime + deltaTime
      WRITE ( UnOut, '(G14.7,1x,I5.5)') CurrentTime, 0
   ENDDO

   CLOSE ( UnOut )

CONTAINS
!............................
SUBROUTINE Cleanup()

   IF (UnIn  > 0) CLOSE( UnIn  )
   IF (UnOut > 0) CLOSE( UnOut )

END SUBROUTINE Cleanup
!............................
END SUBROUTINE CohStr_WriteEvents
!=======================================================================
END MODULE TS_CohStructures
