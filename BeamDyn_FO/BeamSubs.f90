MODULE BeamSubs

   USE NWTC_Library
   USE BeamDyn_Types


   IMPLICIT NONE

   CONTAINS

   SUBROUTINE BD_GetInput(InitInp,p,ErrStat,ErrMsg)

   ! Passed Variables:
   TYPE(BD_InitInputType),       INTENT(INOUT)  :: InitInp
   TYPE(BD_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg


   ! Local variables:
   INTEGER(4)         :: allo_stat
   CHARACTER(1024)    :: FilePath
   LOGICAL            :: Echo
   INTEGER(IntKi)     :: UnIn
   INTEGER(IntKi)     :: UnEc

   INTEGER(IntKi)     :: i

   UnEc = -1
   !-------------------------------------------------------------------------------------------------
   ! Open the BeamDyn input file
   !-------------------------------------------------------------------------------------------------
   CALL GetNewUnit(UnIn)
   CALL OpenFInpFile(UnIn, TRIM(InitInp%BDInputFile), ErrStat)
   IF ( ErrStat /= ErrID_None ) THEN
       ErrStat = ErrID_Fatal
       ErrMsg  = 'Could not open BeamDyn input file: '//InitInp%BDInputFile
       CLOSE(UnIn)
       RETURN
   ENDIF

   CALL GetPath(InitInp%BDInputFile, FilePath)

   !-------------------------- HEADER ---------------------------------------------
   ! Skip header lines
   DO i = 1,3
       CALL ReadCom(UnIn,InitInp%BDInputFile, 'BeamDyn input file header line '//TRIM(Int2LStr(i)), ErrStat, ErrMsg)
       IF (ErrStat /= ErrID_None) THEN
           ErrMsg  = 'Could not read BeamDyn input file header line'
           CLOSE(UnIn)
           RETURN
       ENDIF
   ENDDO 

   CALL ReadLVar(UnIn,InitInp%BDInputFile,Echo,'Echo','Echo Input File Logic Variable',ErrStat, ErrMsg)

   IF ( ErrStat /= ErrID_None ) THEN
        ErrMsg  = 'Error reading Echo Input File Logic Variable'
        CLOSE(UnIn)
        RETURN
   ENDIF
   !-------------------------- GEOMETRY PARAMETERS ----------------------
   ! Skip the comment line.
   CALL ReadCom(UnIn,InitInp%BDInputFile,'GEOMETRY PARAMETERS',ErrStat,ErrMsg)
   IF(ErrStat /= ErrID_None) THEN
      ErrMsg  = 'Could not read BeamDyn GEOMETRY PARAMETERS header line'
      CLOSE(UnIn)
      RETURN
   ENDIF

   ! Read Number of Member
   CALL ReadVar(UnIn,InitInp%BDInputFile,p%member_total,'member_total','BeamDyn Number of Member',ErrStat,ErrMsg,UnEc)
   IF(ErrStat /= ErrID_None) THEN
       CLOSE(UnIn)
       IF(Echo) CLOSE(UnEc)
       RETURN
   ENDIF
   ALLOCATE(p%member_number(p%member_total),STAT=allo_stat)
   IF(allo_stat /= 0) THEN
       ErrMsg = 'Error allocating memeory for the p%member_number in BD_GetInput'
       ErrStat = ErrID_Fatal
       CLOSE(UnIn)
       IF(Echo) CLOSE( UnEc )
       RETURN
   ENDIF
   p%member_number = 0
   ALLOCATE(p%member_coord(3,3,p%member_total),STAT=allo_stat)
   IF(allo_stat /= 0) THEN
       ErrMsg = 'Error allocating memeory for the p%member_coord in BD_GetInput'
       ErrStat = ErrID_Fatal
       CLOSE(UnIn)
       IF(Echo) CLOSE( UnEc )
       RETURN
   ENDIF
   p%member_coord = 0.0D0
   ALLOCATE(p%member_twist(p%member_total,2),STAT=allo_stat)
   IF(allo_stat /= 0) THEN
       ErrMsg = 'Error allocating memeory for the p%member_twist in BD_GetInput'
       ErrStat = ErrID_Fatal
       CLOSE(UnIn)
       IF(Echo) CLOSE( UnEc )
       RETURN
   ENDIF
   p%member_twist = 0.0D0

   DO i=1,p%member_total
       CALL ReadVar(UnIn,InitInp%BDInputFile,p%member_number(i),'member number','BeamDyn Member Number',ErrStat,ErrMsg,UnEc)
       IF(p%member_number(i) /= i) THEN
           ErrMsg = 'Member Number must be consecutive'
           CLOSE(UnIn)
           IF(Echo) CLOSE(UnEc)
           RETURN
       ENDIF
       DO j=1,3
           CALL ReadAry(UnIn,InitInp%BDInputFile,p%member_coord(j,1:3,i),3,'member_coord','Key point coordinates',ErrStat,ErrMsg,UnEc)
       ENDDO
       IF(i /= 1) THEN
           DO j=1,3
               IF(p%member_coord(1,j,i) /= p%member_coord(2,j,i-1)) THEN
                   ErrMsg = 'The staring point of current member must be the ending point of previous member'
                   CLOSE(UnIn)
                   IF(Echo) CLOSE(UnEc)
                   RETURN
               ENDIF
           ENDDO
       ENDIF
       CALL ReadAry(UnIn,InitInp%BDInputFile,p%member_twist(i,1:2),2,'member_twist','Initial twist angle at starting and ending key points',ErrStat,ErrMsg,UnEc)
       IF(i /= 1) THEN
           IF(p%member_twist(i,1) /= p%member_twist(i-1,2)) THEN
               ErrMsg = 'The current member initial twist phi_1 must be equal to phi_2 in the previous member'
               CLOSE(UnIn)
               IF(Echo) CLOSE(UnEc)
               RETURN
           ENDIF
       ENDIF
   ENDDO

   
   END SUBROUTINE BD_GetInput


END MODULE BeamSubs
