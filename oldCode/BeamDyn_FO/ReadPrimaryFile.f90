   SUBROUTINE ReadPrimaryFile(InputFile,InputFileData,&
              OutFileRoot,UnEc,ErrStat,ErrMsg)
   !------------------------------------------------------------------------------------
   ! This routine reads in the primary BeamDyn input file and places the values it reads
   ! in the InputFileData structure.
   !   It opens an echo file if requested and returns the (still-open) echo file to the
   !     calling routine.
   !   It also returns the names of the BldFile, FurlFile, and TrwFile for further 
   !     reading of inputs.
   !------------------------------------------------------------------------------------

   ! Passed variables
   INTEGER(IntKi),               INTENT(  OUT) :: UnEc
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat
   CHARACTER(*),                 INTENT(IN   ) :: InputFile
   CHARACTER(*),                 INTENT(IN   ) :: OutFileRoot
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg

   TYPE(BD_InputFile),           INTENT(INOUT) :: InputFileData
  
   ! Local variables:
   INTEGER(IntKi)               :: UnIn                         ! Unit number for reading file
   INTEGER(IntKi)               :: ErrStat2                     ! Temporary Error status
   LOGICAL                      :: Echo                         ! Determines if an echo file should be written
   CHARACTER(LEN(ErrMsg))       :: ErrMsg2                      ! Temporary Error message
   CHARACTER(1024)              :: PriPath                      ! Path name of the primary file
   CHARACTER(1024)              :: FTitle              ! "File Title": the 2nd line of the input file, which contains a description of its contents
   CHARACTER(1024)              :: BldFile

   INTEGER(IntKi)               :: i
   INTEGER(IntKi)               :: j
   INTEGER(IntKi)               :: temp_int 

   Echo = .FALSE.
   UnEc = -1

   CALL GetNewUnit(UnIn,ErrStat,ErrMsg)

   CALL OpenFInpFile(UnIn,InputFile,ErrStat2,ErrMsg2)

   !-------------------------- HEADER ---------------------------------------------
   CALL ReadCom(UnIn,InputFile,'File Header: Module Version (line 1)',ErrStat2,ErrMsg2,UnEc)
   CALL ReadStr(UnIn,InputFile,FTitle,'FTitle','File Header: File Description (line 2)',ErrStat2, ErrMsg2, UnEc)

   !---------------------- SIMULATION CONTROL --------------------------------------
   CALL ReadCom(UnIn,InputFile,'Section Header: Simulation Control',ErrStat2,ErrMsg2,UnEc)
   CALL ReadVar(UnIn,InputFile,Echo,'Echo','Echo switch',ErrStat2,ErrMsg2,UnEc)
   IF(Echo) THEN
       CALL OpenEcho(UnEc,OutFileRoot//'.ech',ErrStat2,ErrMsg2)
   ENDIF
   IF ( UnEc > 0 )  WRITE(UnEc,*)  'test'
   !---------------------- GEOMETRY PARAMETER --------------------------------------
   CALL ReadCom(UnIn,InputFile,'Section Header: Geometry Parameter',ErrStat2,ErrMsg2,UnEc)
   CALL ReadVar(UnIn,InputFile,InputFileData%member_total,"member_total", "Total number of member",ErrStat2,ErrMsg2,UnEc)
   CALL ReadVar(UnIn,InputFile,InputFileData%ini_curv,"ini_curv", "Initial curvature flag",ErrStat2,ErrMsg2,UnEc)
   CALL ReadCom(UnIn,InputFile,'key point x,y,z locations and initial twist angles',ErrStat2,ErrMsg2,UnEc)
   CALL ReadCom(UnIn,InputFile,'key point and initial twist units',ErrStat2,ErrMsg2,UnEc)

   temp_int = 2*InputFileData%member_total+1
   CALL AllocAry(InputFileData%kp_coordinate,temp_int,3,'Key point coordinates input array',ErrStat2,ErrMsg2)
   CALL AllocAry(InputFileData%initial_twist,InputFileData%member_total+1,'Key point initial twist array',ErrStat2,ErrMsg2)
   InputFileData%kp_coordinate(:,:) = 0.0D0
   InputFileData%initial_twist(:)   = 0.0D0

   j = 0
   DO i=1,temp_int
       IF(MOD(i,2) .NE. 0) THEN
           j = j + 1
           READ(UnIn,*) InputFileData%kp_coordinate(i,2),InputFileData%kp_coordinate(i,3),&
                        InputFileData%kp_coordinate(i,1),InputFileData%initial_twist(j)
       ELSE
           READ(UnIn,*) InputFileData%kp_coordinate(i,2),InputFileData%kp_coordinate(i,3),&
                        InputFileData%kp_coordinate(i,1)
       ENDIF
!       IF(MOD(i,2)==0) THEN
!           IF(InputFileData%initial_twist(i)/=270.0) THEN
!               WRITE(*,*) "Incorrect initial twist angle at mid point:",i
!           ENDIF
!       ENDIF
       IF(UnEc>0) THEN
           IF(i==1) WRITE(UnEc,'(/,A,/)') 'Key points coordinates and initial twist angle'
           WRITE(UnEc,'(/,A,/)') InputFileData%kp_coordinate(i,2),InputFileData%kp_coordinate(i,3),&
                                 InputFileData%kp_coordinate(i,1),InputFileData%initial_twist(i)
       ENDIF
   ENDDO

   !---------------------- MESH PARAMETER -----------------------------------------
   CALL ReadCom(UnIn,InputFile,'Section Header: Mesh Parameter',ErrStat2,ErrMsg2,UnEc)
   CALL ReadVar(UnIn,InputFile,InputFileData%order_elem,"order_elem","Order of basis function",&
                ErrStat2,ErrMsg2,UnEc)
   
   !---------------------- BLADE PARAMETER ----------------------------------------
   CALL ReadCom(UnIn,InputFile,'Section Header: Blade Parameter',ErrStat2,ErrMsg2,UnEc)
   CALL ReadVar ( UnIn, InputFile, InputFileData%BldFile, 'BldFile', 'Name of the file containing properties for blade', ErrStat2, ErrMsg2, UnEc )


   END SUBROUTINE ReadPrimaryFile
