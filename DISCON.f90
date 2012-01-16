!=======================================================================
SUBROUTINE DISCON ( avrSWAP, aviFAIL, accINFILE, avcOUTNAME, avcMSG )
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:'DISCON' :: DISCON


   ! This is a dummy routine for holding the place of a user-specified
   ! Bladed-style DLL controller.  Modify this code to create your own
   ! DLL controller.


IMPLICIT                        NONE


   ! Passed Variables:

REAL(4),    INTENT(INOUT)    :: avrSWAP   (*)                                   ! The swap array, used to pass data to, and receive data from, the DLL controller.

INTEGER(4), INTENT(  OUT)    :: aviFAIL                                         ! A flag used to indicate the success of this DLL call set as follows: 0 if the DLL call was successful, >0 if the DLL call was successful but cMessage should be issued as a warning messsage, <0 if the DLL call was unsuccessful or for any other reason the simulation is to be stopped at this point with cMessage as the error message.

INTEGER(1), INTENT(IN   )    :: accINFILE (*)                                   ! The address of the first record of an array of 1-byte CHARACTERs giving the name of the parameter input file, 'DISCON.IN'.
INTEGER(1), INTENT(  OUT)    :: avcMSG    (*)                                   ! The address of the first record of an array of 1-byte CHARACTERS giving the message contained in cMessage, which will be displayed by the calling program if aviFAIL <> 0.
INTEGER(1), INTENT(IN   )    :: avcOUTNAME(*)                                   ! The address of the first record of an array of 1-byte CHARACTERS giving the simulation run name without extension.


   ! Local Variables:

REAL(4)                      :: rMeasuredPitch                                  ! Current value of the pitch angle of blade 1 (rad).
REAL(4)                      :: rMeasuredSpeed                                  ! Current value of the generator (HSS) speed (rad/s).
REAL(4), SAVE                :: rPitchDemand                                    ! Desired collective pitch angles returned by this DLL (rad).
REAL(4)                      :: rTime                                           ! Current simulation time (sec).

INTEGER(4)                   :: I                                               ! Generic index.
INTEGER(4)                   :: iStatus                                         ! A status flag set by the simulation as follows: 0 if this is the first call, 1 for all subsequent time steps, -1 if this is the final call at the end of the simulation.
 
INTEGER(1)                   :: iInFile   ( 256)                                ! CHARACTER string cInFile  stored as a 1-byte array.
INTEGER(1)                   :: iMessage  ( 256)                                ! CHARACTER string cMessage stored as a 1-byte array.
INTEGER(1), SAVE             :: iOutName  (1024)                                ! CHARACTER string cOutName stored as a 1-byte array.

CHARACTER( 256)              :: cInFile                                         ! CHARACTER string giving the name of the parameter input file, 'DISCON.IN'
CHARACTER( 256)              :: cMessage                                        ! CHARACTER string giving a message that will be displayed by the calling program if aviFAIL <> 0.
CHARACTER(1024), SAVE        :: cOutName                                        ! CHARACTER string giving the simulation run name without extension.


   ! Set EQUIVALENCE relationships between INTEGER(1) byte arrays and CHARACTER strings:

EQUIVALENCE (iInFile , cInFile )
EQUIVALENCE (iMessage, cMessage)
EQUIVALENCE (iOutName, cOutName)



   ! Load variables from calling program (See Appendix A of Bladed User's Guide):

iStatus        = NINT( avrSWAP( 1) )

rTime          =       avrSWAP( 2)
rMeasuredPitch =       avrSWAP( 4)
rMeasuredSpeed =       avrSWAP(20)


   ! Initialize aviFAIL to 0:

aviFAIL        = 0


   ! Read any External Controller Parameters specified in the User Interface
   !   and initialize variables:

IF ( iStatus == 0 )  THEN  ! .TRUE. if were on the first call to the DLL


   ! Convert byte arrays to CHARACTER strings, for convenience:

   DO I = 1,MIN(  256, NINT( avrSWAP(50) ) )
      iInFile (I) = accINFILE (I)   ! Sets cInfile  by EQUIVALENCE
   ENDDO
   DO I = 1,MIN( 1024, NINT( avrSWAP(51) ) )
      iOutName(I) = avcOUTNAME(I)   ! Sets cOutName by EQUIVALENCE
   ENDDO


   ! Read any External Controller Parameters specified in the User Interface:

   ! READ IN DATA CONTAINED IN FILE cInFile HERE
   aviFAIL  = 0   ! SET aviFAIL AND cMessage IF ERROR RESULTS
   cMessage = ''  !


   ! Variable initializations:
   
   rPitchDemand = rMeasuredPitch

ENDIF


   ! Set return values using previous demand if a sample delay is required:

avrSWAP(45) = rPitchDemand


   ! Main control calculations:

IF ( ( iStatus >= 0 ) .AND. ( aviFAIL >= 0 ) )  THEN  ! Only compute control calculations if no error has occured and we are not on the last time step

   ! PLACE MAIN CONTROL CALCULATIONS HERE BASED ON VALUES IN avrSWAP
   ! WRITE OUTPUT DATA TO FILE cOutName (APPENDED WITH AN APPROPRIATE EXTENSION) HERE IF DESIRED
   aviFAIL  = 0   ! SET aviFAIL AND cMessage IF ERROR RESULTS
   cMessage = ''  !

ENDIF


   ! Convert CHARACTER string to byte array for the return message:

DO I = 1,MIN(  256, NINT( avrSWAP(49) ) )
   avcMSG(I) = iMessage(I) ! Same as cMessage by EQUIVALENCE
ENDDO



RETURN
END SUBROUTINE DISCON
!=======================================================================
