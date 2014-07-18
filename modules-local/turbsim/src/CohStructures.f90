!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2014  National Renewable Energy Laboratory
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

   USE                     NWTC_Library   
use ts_errors
use TS_Profiles
use TS_RandNum


   IMPLICIT                NONE

   REAL(ReKi), PARAMETER   ::  KHT_LES_dT = 0.036335           ! The average time step in the LES test file, used here for the KH test
   REAL(ReKi), PARAMETER   ::  KHT_LES_Zm = 6.35475            ! The non-dimensional z dimension defined in LES test file, used here for the KH test
   

CONTAINS

SUBROUTINE CohStr_Open()
use TSMods

   

   IF ( WrACT ) THEN
      p_CohStr%ScaleWid = RotorDiameter * DistScl           !  This is the scaled height of the coherent event data set
      p_grid%Zbottom  = HubHt - CTLz*p_CohStr%ScaleWid             !  This is the height of the bottom of the wave in the scaled/shifted coherent event data set

      IF ( KHtest ) THEN      
            ! for LES test case....
         p_CohStr%ScaleVel = p_CohStr%ScaleWid * KHT_LES_dT /  KHT_LES_Zm    
         p_CohStr%ScaleVel = 50 * p_CohStr%ScaleVel                  ! We want 25 hz bandwidth so multiply by 50
      ELSE
      !   TmpPLExp = PLExp 
      !   PLExp    = MIN( 2.0, 1.35*PLExp )        ! Increase the shear of the background (?)

         p_CohStr%ScaleVel =                     getWindSpeed(UHub,HubHt,p_grid%Zbottom+p_CohStr%ScaleWid,RotorDiameter,PROFILE=WindProfileType)   ! Velocity at the top of the wave
         p_CohStr%ScaleVel = p_CohStr%ScaleVel - getWindSpeed(UHub,HubHt,p_grid%Zbottom,                  RotorDiameter,PROFILE=WindProfileType)   ! Shear across the wave
         p_CohStr%ScaleVel = 0.5 * p_CohStr%ScaleVel                                                                               ! U0 is half the difference between the top and bottom of the billow
      
      !   PLExp = TmpPLExp
      ENDIF

      p_CohStr%Uwave = getWindSpeed(UHub,HubHt,p_grid%Zbottom+0.5*p_CohStr%ScaleWid,RotorDiameter,PROFILE=WindProfileType)                 ! WindSpeed at center of wave

   !BONNIE: MAYBE WE SHOULDN'T OPEN THIS FILE UNTIL WE NEED TO WRITE TO IT
      IF (p_CohStr%ScaleVel < 0. ) THEN
         CALL TS_Warn( ' A coherent turbulence time step file cannot be generated with negative shear.', .TRUE. )
         WrACT = .FALSE.
      ENDIF
   ENDIF

END SUBROUTINE CohStr_Open
!=======================================================================
SUBROUTINE CohStr_ReadEventFile( Un, ScaleWid, ScaleVel, CTKE )

      ! This subroutine reads the events definitions from the event data file

   USE                     TSMods


   IMPLICIT                NONE


      ! Passed Variables

INTEGER,    INTENT(IN)  :: Un             ! I/O Unit
REAL(ReKi),INTENT(IN)   :: CTKE           ! Predicted maximum CTKE
REAL(ReKi),INTENT(INOUT):: ScaleVel       ! The shear we're scaling for
REAL(ReKi),INTENT(IN)   :: ScaleWid       ! The height of the wave we're scaling with

      ! Local variables
REAL(ReKi)              :: MaxEvtCTKE        ! The maximum CTKE in the dataset of events

INTEGER                 :: AllocStat      ! Array allocation status
INTEGER                 :: I              ! DO loop counter
INTEGER                 :: IOS            ! I/O Status


   MaxEvtCTKE = 0.0  ! initialize the MAX variable

   CALL OpenFInpFile ( Un,  p_CohStr%CTEventFile  )
   

         ! Read the nondimensional lateral width of the dataset, Ym_max

   CALL ReadVar( Un, p_CohStr%CTEventFile, Ym_max, "the nondimensional lateral width of the dataset", &
                                                   "Nondimensional lateral dataset width")

         ! Read the nondimensional vertical height of the dataset, Zm_max

   CALL ReadVar( Un, p_CohStr%CTEventFile, Zm_max, "the nondimensional vertical height of the dataset", &
                                          "Nondimensional vertical dataset height")


         ! Read the rest of the header

   CALL ReadVar( Un, p_CohStr%CTEventFile, NumEvents, "NumEvents", "the number of coherent structures.")


   IF ( NumEvents > 0 ) THEN

            ! Allocate memory for coherent event start times and lengths

      ALLOCATE ( EventName(NumEvents) , STAT=AllocStat )

      IF ( AllocStat /= 0 )  THEN
         CALL TS_Abort ( 'Error allocating memory for the coherent event name.' )
      ENDIF

      ALLOCATE ( EventTS(NumEvents) , STAT=AllocStat )

      IF ( AllocStat /= 0 )  THEN
         CALL TS_Abort ( 'Error allocating memory for the coherent event timestep lengths.' )
      ENDIF

      ALLOCATE ( EventLen(NumEvents) , STAT=AllocStat )

      IF ( AllocStat /= 0 )  THEN
         CALL TS_Abort ( 'Error allocating memory for the coherent event lengths.' )
      ENDIF

      ALLOCATE ( pkCTKE(NumEvents) , STAT=AllocStat )

      IF ( AllocStat /= 0 )  THEN
         CALL TS_Abort ( 'Error allocating memory for the coherent event peak CTKE.' )
      ENDIF


            ! Read the last header lines

      CALL ReadCom( Un, p_CohStr%CTEventFile, 'the fourth header line')  ! A blank line
      CALL ReadCom( Un, p_CohStr%CTEventFile, 'the fifth header line')   ! The column heading lines


            ! Read the event definitions and scale times by TScale

      DO I=1,NumEvents

         READ ( Un, *, IOSTAT=IOS )  EventName(I),  EventTS(I), EventLen(I), pkCTKE(I)

         IF ( IOS /= 0 )  THEN
            CALL TS_Abort ( 'Error reading event '//TRIM( Int2LStr( I ) )//' from the coherent event data file.' )
         ENDIF
         MaxEvtCTKE = MAX( MaxEvtCTKE, pkCTKE(I) )

      ENDDO

      IF ( MaxEvtCTKE > 0.0 ) THEN
            ScaleVel = MAX( ScaleVel, SQRT( CTKE / MaxEvtCTKE ) )
            ! Calculate the Velocity Scale Factor, based on the requested maximum CTKE
      ENDIF

         ! Calculate the TimeScaleFactor, based on the Zm_max in the Events file.

      TSclFact = ScaleWid / (ScaleVel * Zm_max)

         ! Scale the time based on TSclFact

      DO I=1,NumEvents
         EventLen(I) = EventLen(I)*TSclFact
      ENDDO

   ELSE

      TSclFact = ScaleWid / (ScaleVel * Zm_max)

   ENDIF  ! FileNum > 0

   CLOSE ( Un )  

END SUBROUTINE CohStr_ReadEventFile
!=======================================================================
SUBROUTINE CohStr_CalcEvents( WindSpeed, MaxCTKE, Height )

      ! This subroutine calculates what events to use and when to use them.
      ! It computes the number of timesteps in the file, NumCTt.

   USE                         TSMods

   IMPLICIT                    NONE

      ! passed variables
REAL(ReKi), INTENT(IN)      :: WindSpeed           ! Hub height wind speed
REAL(ReKi), INTENT(OUT)     :: MaxCTKE             ! Maximum CTKE of events we've picked
REAL(ReKi), INTENT(IN)      :: Height              ! Height for expected length PDF equation

      ! local variables
REAL(ReKi)                  :: iA                  ! Variable used to calculate IAT
REAL(ReKi)                  :: iB                  ! Variable used to calculate IAT
REAL(ReKi)                  :: iC                  ! Variable used to calculate IAT
REAL(ReKi)                  :: rn                  ! random number
REAL(ReKi)                  :: TEnd                ! End time for the current event
REAL(ReKi)                  :: TStartNext   = 0.0  ! temporary start time for next event

INTEGER                     :: IStat               ! Status of memory allocation
INTEGER                     :: NewEvent            ! event number of the new event
INTEGER                     :: NumCompared         ! Number of events we've compared

LOGICAL(1)                  :: Inserted            ! Whether an event was inserted here

TYPE(Event), POINTER        :: PtrCurr  => NULL()  ! Pointer to the current event in the list
TYPE(Event), POINTER        :: PtrNew   => NULL()  ! A new event to be inserted into the list


      ! Compute the mean interarrival time and the expected length of events

   SELECT CASE ( SpecModel )

      CASE ( SpecModel_NWTCUP, SpecModel_NONE, SpecModel_USRVKM )
         y_CohStr%lambda = -0.000904*Rich_No + 0.000562*WindSpeed + 0.001389
         y_CohStr%lambda = 1.0 / y_CohStr%lambda

         IF ( SpecModel == SpecModel_NONE ) THEN
            y_CohStr%ExpectedTime = 600.0
         ELSE
            CALL RndModLogNorm( p_RandNum, OtherSt_RandNum, y_CohStr%ExpectedTime, Height )
         ENDIF

      CASE ( SpecModel_GP_LLJ, SpecModel_SMOOTH, SpecModel_TIDAL, SpecModel_RIVER) ! HYDRO: added 'TIDAL' and 'RIVER' to the spectral models that get handled this way.
         iA     =        0.001797800 + (7.17399E-10)*Height**3.021144723
         iB     =  EXP(-10.590340100 - (4.92440E-05)*Height**2.5)
         iC     = SQRT(  3.655013599 + (8.91203E-06)*Height**3  )
         y_CohStr%lambda = iA + iB*MIN( (WindSpeed**iC), HUGE(iC) )  ! lambda = iA + iB*(WindSpeed**iC)
         y_CohStr%lambda = 1.0 / y_CohStr%lambda

         CALL RndTcohLLJ( p_RandNum, OtherSt_RandNum, y_CohStr%ExpectedTime, Height )

      CASE ( SpecModel_WF_UPW )
        y_CohStr%lambda = 0.000529*WindSpeed + 0.000365*Rich_No - 0.000596
        y_CohStr%lambda = 1.0 / y_CohStr%lambda

         CALL RndTcoh_WF( p_RandNum, OtherSt_RandNum, y_CohStr%ExpectedTime, SpecModel_WF_UPW )

      CASE ( SpecModel_WF_07D )
         y_CohStr%lambda = 0.000813*WindSpeed - 0.002642*Rich_No + 0.002676
         y_CohStr%lambda = 1.0 / y_CohStr%lambda

         CALL RndTcoh_WF( p_RandNum, OtherSt_RandNum, y_CohStr%ExpectedTime, SpecModel_WF_07D )

      CASE ( SpecModel_WF_14D )
         y_CohStr%lambda = 0.001003*WindSpeed - 0.00254*Rich_No - 0.000984
         y_CohStr%lambda = 1.0 / y_CohStr%lambda

         CALL RndTcoh_WF( p_RandNum, OtherSt_RandNum, y_CohStr%ExpectedTime, SpecModel_WF_14D )

      CASE DEFAULT
         !This should not happen

   END SELECT

   y_CohStr%ExpectedTime = y_CohStr%ExpectedTime * ( UsableTime - CTStartTime ) / 600.0  ! Scale for use with the amount of time we've been given


!BONNIE: PERHAPS WE SHOULD JUST PUT IN A CHECK THAT TURNS OFF THE COHERENT TIME STEP FILE IF THE
!        CTSTARTTIME IS LESS THAN THE USABLETIME... MAYBE WHEN WE'RE READING THE INPUT FILE...
y_CohStr%ExpectedTime = MAX( y_CohStr%ExpectedTime, REAL(0.0,ReKi) )  ! This occurs if CTStartTime = 0

      ! We start by adding events at random times

   y_CohStr%NumCTEvents = 0                                    ! Number of events = length of our linked list
   NumCTt      = 0                                    ! Total number of time steps in the events we've picked
   MaxCTKE     = 0.0                                  ! Find the maximum CTKE for the events that we've selected

   y_CohStr%EventTimeSum = 0.0 

   CALL RndExp(p_RandNum, OtherSt_RandNum, rn, y_CohStr%lambda)                            ! Assume the last event ended at time zero

   TStartNext = rn / 2.0

   IF ( KHtest ) THEN
      y_CohStr%ExpectedTime = UsableTime   / 2                 ! When testing, add coherent events for half of the record
      TStartNext   = y_CohStr%ExpectedTime / 2                 ! When testing, start about a quarter of the way into the record
   ENDIF

   IF ( TStartNext < CTStartTime ) THEN
      TStartNext = TStartNext + CTStartTime           ! Make sure the events start after time specified by CTStartTime
   ENDIF

   IF ( TStartNext > 0 ) NumCTt = NumCTt + 1          ! Add a point before the first event

   DO WHILE ( TStartNext < UsableTime .AND. y_CohStr%EventTimeSum < y_CohStr%ExpectedTime )

      CALL RndUnif( p_RandNum, OtherSt_RandNum, rn )

      NewEvent = INT( rn*( NumEvents - 1.0 ) ) + 1
      NewEvent = MAX( 1, MIN( NewEvent, NumEvents ) ) ! take care of possible rounding issues....


      IF ( .NOT. ASSOCIATED ( PtrHead ) ) THEN

         ALLOCATE ( PtrHead, STAT=IStat )             ! The pointer %Next is nullified in allocation

         IF ( IStat /= 0 ) THEN
            CALL TS_Abort ( 'Error allocating memory for new event.' )
         ENDIF

         PtrTail => PtrHead

      ELSE

         ALLOCATE ( PtrTail%Next, STAT=IStat )     ! The pointer PtrTail%Next%Next is nullified in allocation

         IF ( IStat /= 0 ) THEN
            CALL TS_Abort ( 'Error allocating memory for new event.' )
         ENDIF

         PtrTail => PtrTail%Next                   ! Move the pointer to point to the last record in the list

      ENDIF

      PtrTail%EventNum     = NewEvent
      PtrTail%TStart       = TStartNext
      PtrTail%delt         = EventLen( NewEvent ) / EventTS( NewEvent )          ! the average delta time in the event
      PtrTail%Connect2Prev = .FALSE.

      MaxCTKE              = MAX( MaxCTKE, pkCTKE( NewEvent ) )
      y_CohStr%NumCTEvents          = y_CohStr%NumCTEvents + 1

      TEnd = TStartNext + EventLen( NewEvent )


      IF ( KHtest ) THEN
         TStartNext   = UsableTime + TStartNext !TEnd + PtrTail%delt ! Add the events right after each other
      ELSE

         DO WHILE ( TStartNext <= TEnd )

            CALL RndExp(p_RandNum, OtherSt_RandNum, rn, y_CohStr%lambda)    ! compute the interarrival time
            TStartNext        = TStartNext + rn !+ EventLen( NewEvent )

         ENDDO

      ENDIF


      IF ( (TStartNext - TEnd) > PtrTail%delt ) THEN
         NumCTt = NumCTt + EventTS( NewEvent ) + 2                                  ! add a zero-line (essentially a break between events)
      ELSE
         NumCTt = NumCTt + EventTS( NewEvent ) + 1
      ENDIF

      y_CohStr%EventTimeSum     = y_CohStr%EventTimeSum + EventLen( NewEvent )

   ENDDO


      ! Next, we start concatenating events until there is no space or we exceed the expected time

   IF ( SpecModel /= SpecModel_NONE ) THEN

      NumCompared = 0

      DO WHILE ( y_CohStr%EventTimeSum < y_CohStr%ExpectedTime .AND. NumCompared < y_CohStr%NumCTEvents )

         CALL RndUnif( p_RandNum, OtherSt_RandNum, rn )

         NewEvent = INT( rn*( NumEvents - 1.0 ) ) + 1
         NewEvent = MAX( 1, MIN( NewEvent, NumEvents ) )    ! take care of possible rounding issues....

         NumCompared = 0
         Inserted    = .FALSE.

         DO WHILE ( NumCompared < y_CohStr%NumCTEvents .AND. .NOT. Inserted )

            IF ( .NOT. ASSOCIATED ( PtrCurr ) ) THEN        ! Wrap around to the beginning of the list
               PtrCurr => PtrHead
            ENDIF


               ! See if the NewEvent fits between the end of event pointed to by PtrCurr and the
               ! beginning of the event pointed to by PtrCurr%Next

            IF ( ASSOCIATED( PtrCurr%Next ) ) THEN
               TStartNext = PtrCurr%Next%TStart
            ELSE !We're starting after the last event in the record
               TStartNext = UsableTime + 0.5 * EventLen( NewEvent )  ! We can go a little beyond the end...
            ENDIF

            IF ( TStartNext - (PtrCurr%TStart + EventLen( PtrCurr%EventNum ) + PtrCurr%delt) > EventLen( NewEvent ) ) THEN

               Inserted = .TRUE.

               ALLOCATE ( PtrNew, STAT=IStat )           ! The pointer %Next is nullified in allocation

               IF ( IStat /= 0 ) THEN
                  CALL TS_Abort ( 'Error allocating memory for new event.' )
               ENDIF

               PtrNew%EventNum      = NewEvent
               PtrNew%TStart        = PtrCurr%TStart + EventLen( PtrCurr%EventNum )
               PtrNew%delt          = EventLen( NewEvent ) / EventTS( NewEvent )          ! the average delta time in the event
               PtrNew%Connect2Prev  = .TRUE.

               PtrNew%Next  => PtrCurr%Next
               PtrCurr%Next => PtrNew
               PtrCurr      => PtrCurr%Next    ! Let's try to add the next event after the other events

               MaxCTKE              = MAX( MaxCTKE, pkCTKE( NewEvent ) )
               y_CohStr%NumCTEvents          = y_CohStr%NumCTEvents + 1
               NumCTt               = NumCTt + EventTS( NewEvent )  ! there is no break between events
                                   !(we may have one too many NumCTt here, so we'll deal with it when we write the file later)
               y_CohStr%EventTimeSum         = y_CohStr%EventTimeSum + EventLen( NewEvent )


            ELSE

               NumCompared = NumCompared + 1

            ENDIF

            PtrCurr => PtrCurr%Next

         ENDDO ! WHILE (NumCompared < NumCTEvents .AND. .NOT. Inserted)

      ENDDO ! WHILE (EventTimeSum < ExpectedTime .AND. NumCompared < NumCTEvents)

   ENDIF ! SpecModel /= SpecModel_NONE

IF ( NumCTt > 0 ) THEN
   EventTimeStep = y_CohStr%EventTimeSum / NumCTt                                          ! Average timestep of coherent event data
ELSE
   EventTimeStep = 0.0
ENDIF



END SUBROUTINE CohStr_CalcEvents
!=======================================================================
END MODULE TS_CohStructures
