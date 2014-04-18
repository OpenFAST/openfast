MODULE BeamSubs

   USE NWTC_Library
   USE BeamDyn_Types


   IMPLICIT NONE

   CONTAINS

   SUBROUTINE BD_GetInput(InitInp,p,x,xd,z,O,y,ErrStat,ErrMess)

   ! Passed Variables:
   TYPE(BDyn_InitInputType),       INTENT(INOUT)  :: InitInp
   TYPE(BDyn_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
   TYPE(BDyn_ContinuousStateType), INTENT(INOUT)  :: x           ! Initial continuous states
   TYPE(BDyn_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Initial discrete states
   TYPE(BDyn_ConstraintStateType), INTENT(INOUT)  :: z           ! Initial guess of the constraint states
   TYPE(BDyn_OtherStateType),      INTENT(INOUT)  :: O !therState  ! Initial other/optimization states
   TYPE(BDyn_OutputType),          INTENT(INOUT)  :: y           ! Initial system outputs (outputs are not calculated;
   INTEGER, INTENT(OUT)                   :: ErrStat
   CHARACTER(*), INTENT(OUT)              :: ErrMess


   ! Local variables:
   CHARACTER(1024)    :: FilePath
   !-------------------------------------------------------------------------------------------------
   ! Open the BeamDyn input file
   !-------------------------------------------------------------------------------------------------
   CALL OpenFInpFile(UnIn, TRIM(InitInp%BDFileName), ErrStat)
   IF (ErrStat /= ErrID_None ) RETURN

   CALL GetPath( InitInp%BDFileName, FilePath )

   !-------------------------------------------------------------------------------------------------
   ! If the echo file is open, write the header...
   !-------------------------------------------------------------------------------------------------
   IF ( p%Echo ) THEN
      WRITE( p%UnEc, '(// A /)' ) 'BeamDyn input data from file "'//TRIM( InitInp%bDFileName )//'":'
   ENDIF

   !-------------------------------------------------------------------------------------------------
   ! Read the BeamDyn input file
   !-------------------------------------------------------------------------------------------------

   ! Read in the title line
   CALL ReadStr( UnIn, InitInp%BDFileName, InitInp%Title, VarName='Title', VarDescr='File title', ErrStat=ErrStat)
   IF ( ErrStat /= ErrID_None ) RETURN

   CALL WrScr( ' Heading of the BeamDyn input file: '//NewLine//'   '//TRIM(InitInp%Title) )
   p%TITLE = InitInp%TITLE


   END SUBROUTINE BD_GetInput


END MODULE BeamSubs
