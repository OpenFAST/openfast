   SUBROUTINE BeamDyn_ReadInput(InputFileName,&
!             &BladeFileName,
             &InputFileData,ErrStat,ErrMsg)

   ! Passed Variables:
   CHARACTER(*),                 INTENT(IN   )  :: InputFileName    ! Name of the input file
!   CHARACTER(*),                 INTENT(IN   )  :: BladeFileName    ! File that contains the blade properties information

   TYPE(BD_InputFile),           INTENT(  OUT)  :: InputFileData    ! Data stored in the module's input file
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat          ! The error status code
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg           ! The error message, if an error occurred


   ! Local variables:
!   INTEGER(4)         :: allo_stat
!   CHARACTER(1024)    :: FilePath
   LOGICAL            :: Echo
   INTEGER(IntKi)                               :: UnEcho
   INTEGER(IntKi)                               :: ErrStat2
   CHARACTER(LEN(ErrMsg))                       :: ErrMsg2
   CHARACTER(1024)                        :: BldFile ! File that contains the blade information (specified in the primary input file)

   INTEGER(IntKi)     :: i

   ErrStat = ErrID_None
   ErrMsg = ''

   CALL ReadPrimaryFile(InputFileName,InputFileData,&
                       &BldFile,
                       &UnEcho,ErrStat2,ErrMsg2)

   CALL ReadBladeInput
   
   END SUBROUTINE BeamDyn_ReadInput


