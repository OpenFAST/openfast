   SUBROUTINE BD_ReadBladeFile(BldFile,BladeInputFileData,UnEc,ErrStat,ErrMsg)

   ! Passed variables:
   TYPE(BladeInputData), INTENT(  OUT):: BladeInputFileData
   CHARACTER(*),         INTENT(IN   ):: BldFile
   INTEGER(IntKi),       INTENT(IN   ):: UnEc
   INTEGER(IntKi),       INTENT(  OUT):: ErrStat                             ! Error status
   CHARACTER(*),         INTENT(  OUT):: ErrMsg                              ! Error message
   
   ! Local variables:
   INTEGER(IntKi)             :: UnIn                                            ! Unit number for reading file
   INTEGER(IntKi)             :: ErrStat2                                        ! Temporary Error status
   CHARACTER(LEN(ErrMsg))     :: ErrMsg2                                         ! Temporary Err msg   
   REAL(ReKi)                 :: temp_xm2
   REAL(ReKi)                 :: temp_xm3
   REAL(ReKi)                 :: temp_mass(5)
   REAL(ReKi)                 :: temp_bend(6)
   REAL(ReKi)                 :: temp_sher(6)
   REAL(ReKi)                 :: temp_edg
   REAL(ReKi)                 :: temp_flp
   REAL(ReKi)                 :: temp_crs
   INTEGER(IntKi)             :: i
   INTEGER(IntKi)             :: j

   CALL GetNewUnit(UnIn,ErrStat,ErrMsg)

   CALL OpenFInpFile (UnIn,BldFile,ErrStat2,ErrMsg2)

   !  -------------- HEADER -------------------------------------------------------
   ! Skip the header.
   CALL ReadCom(UnIn,BldFile,'unused blade file header line 1',ErrStat2,ErrMsg2,UnEc)
   CALL ReadCom(UnIn,BldFile,'unused blade file header line 2',ErrStat2,ErrMsg2,UnEc)

   !  -------------- BLADE PARAMETER-----------------------------------------------
   CALL ReadCom(UnIn,BldFile,'blade parameters',ErrStat2,ErrMsg2,UnEc)

   CALL ReadVar(UnIn,BldFile,BladeInputFileData%station_total,'station_total','Number of blade input stations',ErrStat2,ErrMsg2,UnEc)

   CALL AllocAry(BladeInputFileData%stiff0,6,6,BladeInputFileData%station_total,'Cross-sectional 6 by 6 stiffness matrix',ErrStat2,ErrMsg2)
   BladeInputFileData%stiff0(:,:,:) = 0.0D0 
   CALL AllocAry(BladeInputFileData%mass0,6,6,BladeInputFileData%station_total,'Cross-sectional 6 by 6 mass matrix',ErrStat2,ErrMsg2)
   BladeInputFileData%mass0(:,:,:) = 0.0D0 
   CALL AllocAry(BladeInputFileData%station_eta,BladeInputFileData%station_total,'Station eta array',ErrStat2,ErrMsg2)
   BladeInputFileData%station_eta(:) = 0.0D0 

   CALL ReadVar(UnIn,BldFile,BladeInputFileData%damp_flag,'damp_flag','Damping flag',ErrStat2,ErrMsg2,UnEc)
   !  -------------- DAMPING PARAMETER-----------------------------------------------
   CALL ReadCom(UnIn,BldFile,'damping parameters',ErrStat2,ErrMsg2,UnEc)
   CALL ReadCom(UnIn,BldFile,'mu1 to mu6',ErrStat2,ErrMsg2,UnEc)
   CALL ReadCom(UnIn,BldFile,'units',ErrStat2,ErrMsg2,UnEc)
   CALL AllocAry(BladeInputFileData%beta,6,'Number of damping coefficient',ErrStat2,ErrMsg2)
   CALL ReadAry(UnIn,BldFile,BladeInputFileData%beta(:),6,'damping coefficient','damping coefficient',ErrStat2,ErrMsg2,UnEc)
!  -------------- DISTRIBUTED PROPERTIES--------------------------------------------
   CALL ReadCom(UnIn,BldFile,'Distributed properties',ErrStat2,ErrMsg2,UnEc)
   DO i=1,BladeInputFileData%station_total
       READ(UnIn,*) BladeInputFileData%station_eta(i)
       DO j=1,6
           CALL ReadAry(UnIn,BldFile,BladeInputFileData%stiff0(j,:,i),6,'siffness_matrix',&
                   'Blade C/S stiffness matrix',ErrStat2,ErrMsg2,UnEc)
       ENDDO
       DO j=1,6
           CALL ReadAry(UnIn,BldFile,BladeInputFileData%mass0(j,:,i),6,'mass_matrix',&
                   'Blade C/S mass matrix',ErrStat2,ErrMsg2,UnEc)
       ENDDO
   ENDDO

   END SUBROUTINE BD_ReadBladeFile
