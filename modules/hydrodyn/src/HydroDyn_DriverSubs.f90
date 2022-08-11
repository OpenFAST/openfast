!**********************************************************************************************************************************
! HydroDyn_DriverSubs
!..................................................................................................................................
! LICENSING
! Copyright (C) 2022 Envision Energy USA LTD
!
!    This file is part of HydroDyn.
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

MODULE HydroDynDriverSubs

   USE NWTC_Library
   use SeaState
   use SeaState_Types
   USE HydroDyn
   USE HydroDyn_Types
   USE HydroDyn_Output
   USE ModMesh_Types
   USE VersionInfo
   
   IMPLICIT NONE
   
   TYPE HD_Drvr_MappingData
      Integer(IntKi)                   :: Ind                                     ! index for interpolation
      type(MeshMapType)                :: HD_Ref_2_WB_P                           ! Mesh mapping between HD Reference pt mesh and WAMIT body(ies) mesh
      type(MeshMapType)                :: HD_Ref_2_M_P                            ! Mesh mapping between HD Reference pt mesh and Morison mesh
      
      ! For 6x6 linearization
      type(MeshType)                   :: EDRPtMesh                               ! 1-node Point mesh located at (0,0,zRef) in global system where ElastoDyn Reference point is
      type(MeshType)                   :: ZZZPtMeshMotion                         ! 1-node Point mesh located at (0,0,0) in global system and never moving
      type(MeshType)                   :: ZZZPtMeshLoads                          ! 1-node Point mesh located at (0,0,0) in global system and never moving
      type(MeshMapType)                :: ED_Ref_2_HD_Ref                         ! Mesh mapping between ED Reference pt mesh and HD PRP mesh
      type(MeshMapType)                :: HD_Ref_2_ED_Ref                         ! Mesh mapping between HD Reference pt mesh and ED ref poing mesh
      type(MeshMapType)                :: HD_RefLoads_2_ED_Ref                    ! Mesh mapping between HDHdroOrigin pt mesh and ED ref point mesh for loads
      type(MeshMapType)                :: HD_RefLoads_2_ZZZLoads                  ! Mesh mapping between HDHdroOrigin pt mesh and ZZZPtMesh
      
   END TYPE HD_Drvr_MappingData
   
   TYPE HD_Drvr_OutputFile
      INTEGER                          :: NumOuts
      INTEGER                          :: NumOutsMods(2)
      CHARACTER(ChanLen), ALLOCATABLE  :: WriteOutputHdr(:)
      CHARACTER(ChanLen), ALLOCATABLE  :: WriteOutputUnt(:)
      REAL(ReKi),         ALLOCATABLE  :: Storage(:,:)
      CHARACTER(500)                   :: FileDescLines(3)
      INTEGER                          :: unOutFile = -1
      CHARACTER(20)                    :: OutFmt
      CHARACTER(20)                    :: OutFmt_t
      INTEGER                          :: n_Out = 0
      REAL(DbKi)                       :: TimeData(2)
   END TYPE HD_Drvr_OutputFile
   
   TYPE HD_Drvr_Data
      LOGICAL                          :: Echo
      REAL(ReKi)                       :: Gravity
      REAL(ReKi)                       :: WtrDens
      REAL(ReKi)                       :: WtrDpth
      REAL(ReKi)                       :: MSL2SWL
      CHARACTER(1024)                  :: HDInputFile
      CHARACTER(1024)                  :: SeaStateInputFile
      CHARACTER(1024)                  :: OutRootName
      LOGICAL                          :: Linearize
      LOGICAL                          :: WrTxtOutFile = .true.
      LOGICAL                          :: WrBinOutFile = .false.
      INTEGER                          :: NSteps
      REAL(DbKi)                       :: TimeInterval
      REAL(DbKi)                       :: TMax
      INTEGER                          :: PRPInputsMod
      CHARACTER(1024)                  :: PRPInputsFile
      REAL(ReKi)                       :: uPRPInSteady(6)
      REAL(ReKi)                       :: uDotPRPInSteady(6)
      REAL(ReKi)                       :: uDotDotPRPInSteady(6)
      REAL(ReKi), ALLOCATABLE          :: PRPin(:,:)           ! Variable for storing time, forces, and body velocities, in m/s or rad/s for PRP
      REAL(ReKi), ALLOCATABLE          :: PRPinTime(:)         ! Variable for storing time, forces, and body velocities, in m/s or rad/s for PRP
      INTEGER(IntKi)                   :: NBody                ! Number of WAMIT bodies to work with if prescribing kinematics on each body (PRPInputsMod<0)
      REAL(ReKi)                       :: PtfmRefzt
      TYPE(HD_Drvr_OutputFile)         :: OutData
      character(500)                   :: FTitle                  ! description from 2nd line of driver file
      
   END TYPE HD_Drvr_Data
   
! -----------------------------------------------------------------------------------   
! NOTE:  this module and the ModMesh.f90 modules must use the Fortran compiler flag:  
!        /fpp                  because of they both have preprocessor statements
! ----------------------------------------------------------------------------------- 
   TYPE(ProgDesc), PARAMETER        :: version   = ProgDesc( 'HydroDyn Driver', '', '' )  ! The version number of this program.
   character(*), parameter          :: Delim = Tab


CONTAINS

SUBROUTINE ReadDriverInputFile( FileName, drvrData, ErrStat, ErrMsg )

   CHARACTER(*),                  INTENT( IN    ) :: FileName
   TYPE(HD_Drvr_Data),            INTENT( INOUT ) :: drvrData
   INTEGER,                       INTENT(   OUT ) :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                  INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
      ! Local variables  

   INTEGER                                          :: UnIn                 ! Unit number for the input file
   INTEGER                                          :: UnEchoLocal          ! The local unit number for this module's echo file
   CHARACTER(1024)                                  :: EchoFile             ! Name of HydroDyn echo file  
   CHARACTER(1024)                                  :: PriPath              ! Temporary storage for relative path name

   integer(IntKi)                                   :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                             :: errMsg2       ! temporary error message 
   character(*), parameter                          :: RoutineName = 'ReadDriverInputFile'
   
   
      ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEchoLocal = -1
   ErrStat = ErrID_None
   ErrMsg = ""
      
   CALL GetNewUnit( UnIn ) 
   CALL OpenFInpFile ( UnIn, FileName, ErrStat2, ErrMsg2 ) 
   if (Failed()) return


   CALL WrScr( 'Opening HydroDyn Driver input file:  '//trim(FileName) )
   call GetPath( TRIM(FileName), PriPath ) ! store path name in case any of the file names are relative to the primary input file

   
   !-------------------------------------------------------------------------------------------------
   ! File header
   !-------------------------------------------------------------------------------------------------
   
   CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 1', ErrStat2, ErrMsg2 )
   if (Failed()) return


   CALL ReadStr( UnIn, FileName, drvrData%FTitle, 'FTitle', 'HydroDyn Driver input file header line 2', ErrStat2, ErrMsg2 )
   if (Failed()) return


     ! Echo Input Files.
   CALL ReadVar ( UnIn, FileName, drvrData%Echo, 'Echo', 'Echo Input', ErrStat2, ErrMsg2 )
   if (Failed()) return
   
   
      ! If we are Echoing the input then we should re-read the first three lines so that we can echo them
      ! using the NWTC_Library routines.  The echoing is done inside those routines via a global variable
      ! which we must store, set, and then replace on error or completion.
      
   IF ( drvrData%Echo ) THEN
      
      EchoFile = TRIM(FileName)//'.ech'
      CALL GetNewUnit( UnEchoLocal )   
      CALL OpenEcho ( UnEchoLocal, EchoFile, ErrStat2, ErrMsg2 )
      if (Failed()) return

      
      REWIND(UnIn)
      
      CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 1', ErrStat2, ErrMsg2, UnEchoLocal )
      if (Failed()) return

      CALL ReadCom( UnIn, FileName, 'HydroDyn Driver input file header line 2', ErrStat2, ErrMsg2, UnEchoLocal )
      if (Failed()) return

         ! Echo Input Files. Note this line is prevented from being echoed by the ReadVar routine.
      CALL ReadVar ( UnIn, FileName, drvrData%Echo, 'Echo', 'Echo the input file data', ErrStat2, ErrMsg2, UnEchoLocal )
      if (Failed()) return

      
   END IF
   !-------------------------------------------------------------------------------------------------
   ! Environmental conditions section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName, 'Environmental conditions header', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return

      ! Gravity - Gravity.
   CALL ReadVar ( UnIn, FileName, drvrData%Gravity, 'Gravity', 'Gravity', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return

      ! WtrDens - Water density.
   CALL ReadVar ( UnIn, FileName, drvrData%WtrDens, 'WtrDens', 'Water density', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return

      ! WtrDpth - Water depth.
   CALL ReadVar ( UnIn, FileName, drvrData%WtrDpth, 'WtrDpth', 'Water depth', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return

      ! MSL2SWL - Offset between still-water level and mean sea level.
   CALL ReadVar ( UnIn, FileName, drvrData%MSL2SWL, 'MSL2SWL', 'Offset between still-water level and mean sea level', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   
   !-------------------------------------------------------------------------------------------------
   ! HYDRODYN section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName, 'HYDRODYN header', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   
      ! HDInputFile
   CALL ReadVar ( UnIn, FileName, drvrData%HDInputFile, 'HDInputFile', 'HydroDyn input filename', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   IF ( PathIsRelative( drvrData%HDInputFile ) ) drvrData%HDInputFile = TRIM(PriPath)//TRIM(drvrData%HDInputFile)

       ! SeaStInputFile
   CALL ReadVar ( UnIn, FileName, drvrData%SeaStateInputFile, 'SeaStateInputFile', 'SeaState input filename', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   IF ( PathIsRelative( drvrData%SeaStateInputFile ) ) drvrData%SeaStateInputFile = TRIM(PriPath)//TRIM(drvrData%SeaStateInputFile)

      ! OutRootName
   CALL ReadVar ( UnIn, FileName, drvrData%OutRootName, 'OutRootName', 'HydroDyn output root filename', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   IF ( PathIsRelative( drvrData%OutRootName ) ) drvrData%OutRootName = TRIM(PriPath)//TRIM(drvrData%OutRootName)

       ! Linearize
   CALL ReadVar ( UnIn, FileName, drvrData%Linearize, 'Linearize', 'Linearize parameter', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
  
      ! NSteps
   CALL ReadVar ( UnIn, FileName, drvrData%NSteps, 'NSteps', 'Number of time steps in the HydroDyn simulation', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   
      ! TimeInterval   
   CALL ReadVar ( UnIn, FileName, drvrData%TimeInterval, 'TimeInterval', 'Time interval for any HydroDyn inputs', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   
   
   !-------------------------------------------------------------------------------------------------
   ! PRP INPUTS section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName, 'PRP INPUTS header', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   
      ! PRPInputsMod
   CALL ReadVar ( UnIn, FileName, drvrData%PRPInputsMod, 'PRPInputsMod', 'Model for the PRP (principal reference point) inputs', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   
       ! PtfmRefzt
   CALL ReadVar ( UnIn, FileName, drvrData%PtfmRefzt, 'PtfmRefzt', 'Vertical distance from the ground level to the platform reference point', ErrStat, ErrMsg, UnEchoLocal )
   if (Failed()) return
   
      ! PRPInputsFile
   CALL ReadVar ( UnIn, FileName, drvrData%PRPInputsFile, 'PRPInputsFile', 'Filename for the PRP HydroDyn inputs', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return
   IF ( PathIsRelative( drvrData%PRPInputsFile ) ) drvrData%PRPInputsFile = TRIM(PriPath)//TRIM(drvrData%PRPInputsFile)
   
   
   !-------------------------------------------------------------------------------------------------
   ! PRP STEADY STATE INPUTS section
   !-------------------------------------------------------------------------------------------------

      ! Header
   CALL ReadCom( UnIn, FileName, 'PRP STEADY STATE INPUTS header', ErrStat2, ErrMsg2, UnEchoLocal )
   if (Failed()) return

      ! uPRPInSteady
   CALL ReadAry ( UnIn, FileName, drvrData%uPRPInSteady, 6, 'uPRPInSteady', 'PRP Steady-state displacements and rotations.', ErrStat2,  ErrMsg2, UnEchoLocal)
   if (Failed()) return
   
      ! uDotPRPInSteady
   CALL ReadAry ( UnIn, FileName, drvrData%uDotPRPInSteady, 6, 'uDotPRPInSteady', 'PRP Steady-state translational and rotational velocities.', ErrStat2,  ErrMsg2, UnEchoLocal)
   if (Failed()) return
      
      ! uDotDotPRPInSteady
   CALL ReadAry ( UnIn, FileName, drvrData%uDotDotPRPInSteady, 6, 'uDotDotPRPInSteady', 'PRP Steady-state translational and rotational accelerations.', ErrStat2,  ErrMsg2, UnEchoLocal)
   if (Failed()) return

      
   IF ( drvrData%PRPInputsMod /= 1 ) THEN
      drvrData%uPRPInSteady       = 0.0
      drvrData%uDotPRPInSteady    = 0.0
      drvrData%uDotDotPRPInSteady = 0.0
   END IF

   drvrData%WrTxtOutFile = .true.
   drvrData%WrBinOutFile = .false.
   

   CALL cleanup()
   
CONTAINS

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
      Failed =  ErrStat >= AbortErrLev
      if (Failed) call Cleanup()
        
   end function Failed

   subroutine Cleanup()
      CLOSE( UnIn )
      IF ( UnEchoLocal > 0 ) CLOSE( UnEchoLocal )
   end subroutine Cleanup
   
END SUBROUTINE ReadDriverInputFile
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE ReadPRPInputsFile( drvrData, ErrStat, ErrMsg )

   TYPE(HD_Drvr_Data),            INTENT( INOUT ) :: drvrData
   INTEGER,                       INTENT(   OUT ) :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                  INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
      ! Local variables  

   INTEGER                                          :: UnIn                 ! Unit number for the input file
   INTEGER                                          :: UnEchoLocal          ! The local unit number for this module's echo file
!   CHARACTER(1024)                                  :: EchoFile             ! Name of HydroDyn echo file  

   integer(IntKi)                                   :: n, sizeAry
   integer(IntKi)                                   :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                             :: errMsg2       ! temporary error message 
   character(*), parameter                          :: RoutineName = 'ReadDriverInputFile'
   real(ReKi), allocatable                          :: TmpAry(:)
   
   
      ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEchoLocal = -1
   UnIn = -1
   
   ErrStat = ErrID_None
   ErrMsg = ""
      
   drvrData%NBody= 0
   
   IF ( drvrData%PRPInputsMod == 2 ) THEN
      sizeAry = 19
   ELSEIF ( drvrData%PRPInputsMod < 0 ) THEN
      ! multi-body kinematics driver option (time, PRP DOFs 1-6, body1 DOFs 1-6, body2 DOFs 1-6...)
      
      drvrData%NBody = -drvrData%PRPInputsMod
      sizeAry = 7 + 6*drvrData%NBody

      call WrScr( 'NBody is '//trim(Num2LStr(drvrData%NBody))//' and planning to read in  '//trim(Num2LStr(sizeAry))//' columns from the input file' )
      
   ELSE
   
      RETURN
      
   END IF
   
   CALL AllocAry(TmpAry, sizeAry, 'TmpAry', ErrStat2, ErrMsg2)
      if (Failed()) return
   CALL AllocAry(drvrData%PRPin, drvrData%NSteps, sizeAry-1, 'PRPin', ErrStat2, ErrMsg2)
      if (Failed()) return
   CALL AllocAry(drvrData%PRPinTime, drvrData%NSteps, 'PRPinTime', ErrStat2, ErrMsg2)
      if (Failed()) return
      
   
      ! Open the (PRP or WAMIT) inputs data file
   CALL GetNewUnit( UnIn ) 
   CALL OpenFInpFile ( UnIn, trim(drvrData%PRPInputsFile), ErrStat2, ErrMsg2 )
      if (Failed()) return
   
      !seems like it would be more efficient to switch the indices on drvrData%PRPin
   DO n = 1,drvrData%NSteps
      CALL ReadAry ( UnIn, drvrData%PRPInputsFile, TmpAry, sizeAry, 'Line', 'drvrData%PRPin', ErrStat2, ErrMsg2, UnEchoLocal )
      drvrData%PRPin(n,:) = TmpAry(2:sizeAry)
      drvrData%PRPinTime(n) = TmpAry(1)
      if (Failed()) return
   END DO
   
   
   
   call Cleanup()
   
CONTAINS

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName) 
      Failed =  ErrStat >= AbortErrLev
      if (Failed) call Cleanup()
        
   end function Failed

   subroutine Cleanup()
      IF ( ALLOCATED(TmpAry) ) DEALLOCATE(TmpAry)
      IF ( UnIn > 0 ) CLOSE( UnIn )
      IF ( UnEchoLocal > 0 ) CLOSE( UnEchoLocal )
   end subroutine Cleanup
   
   

END SUBROUTINE ReadPRPInputsFile
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE InitOutputFile(InitOutData_HD, InitOutData_SeaSt, drvrData, ErrStat, ErrMsg)

   TYPE(HydroDyn_InitOutputType),   INTENT(IN)      :: InitOutData_HD          ! Output data from initialization
   TYPE(SeaSt_InitOutputType),      INTENT(IN)      :: InitOutData_SeaSt       ! Output data from initialization
   TYPE(HD_Drvr_Data),              INTENT( INOUT ) :: drvrData
   INTEGER,                         INTENT(   OUT ) :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                    INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None

   integer(IntKi)                                   :: FmtWidth, TChanLen
   integer(IntKi)                                   :: i, Indx
   integer(IntKi)                                   :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                             :: errMsg2       ! temporary error message 
   character(*), parameter                          :: RoutineName = 'InitOutputFile'
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   drvrData%OutData%n_Out = 0
   drvrData%OutData%OutFmt = "ES15.6E2"
   CALL ChkRealFmtStr( drvrData%OutData%OutFmt, 'OutFmt', FmtWidth, ErrStat2, ErrMsg2 )
   !IF ( drvrData%WrTxtOutFile .and. FmtWidth < MinChanLen ) CALL SetErrStat( ErrID_Warn, 'OutFmt produces a column width of '// &
   !      TRIM(Num2LStr(FmtWidth))//'), which may be too small.', ErrStat, ErrMsg, RoutineName )
   
   if (drvrData%TMax < 1.0_DbKi) then ! log10(0) gives floating point divide-by-zero error
      TChanLen = MinChanLen
   else
      TChanLen = max( MinChanLen, int(log10(drvrData%TMax))+7 )
   end if
   drvrData%OutData%OutFmt_t = 'F'//trim(num2lstr( TChanLen ))//'.4' ! 'F10.4'
   
   

   drvrData%OutData%NumOutsMods = 0
   if (Allocated(InitOutData_SeaSt%WriteOutputHdr)) drvrData%OutData%NumOutsMods(1) = size(InitOutData_SeaSt%WriteOutputHdr)
   if (Allocated(InitOutData_HD%WriteOutputHdr   )) drvrData%OutData%NumOutsMods(2) = size(InitOutData_HD%WriteOutputHdr)
   drvrData%OutData%NumOuts = sum(drvrData%OutData%NumOutsMods) + 1 ! add 1 for time channel
   
   call AllocAry(drvrData%OutData%WriteOutputHdr, drvrData%OutData%NumOuts, ' DriverWriteOutputHdr', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(drvrData%OutData%WriteOutputUnt, drvrData%OutData%NumOuts, ' DriverWriteOutputUnt', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   call AllocAry(drvrData%OutData%Storage, drvrData%OutData%NumOuts-1, drvrData%NSteps, ' DriverWriteOutputStorage', ErrStat2, ErrMsg2); call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   
   IF (ErrStat >= AbortErrLev) RETURN
   
   ! Fill concatenated WriteOuput header and unit arrays:
   drvrData%OutData%WriteOutputHdr(1) = 'Time'
   drvrData%OutData%WriteOutputUnt(1) = '(s)'
   Indx = 1
   do i=1,drvrData%OutData%NumOutsMods(1)
      Indx = Indx + 1
      drvrData%OutData%WriteOutputHdr(Indx) = InitOutData_SeaSt%WriteOutputHdr(i)
      drvrData%OutData%WriteOutputUnt(Indx) = InitOutData_SeaSt%WriteOutputUnt(i)
   end do
   
   do i=1,drvrData%OutData%NumOutsMods(2)
      Indx = Indx + 1
      drvrData%OutData%WriteOutputHdr(Indx) = InitOutData_HD%WriteOutputHdr(i)
      drvrData%OutData%WriteOutputUnt(Indx) = InitOutData_HD%WriteOutputUnt(i)
   end do
   
   ! get lines for output file:
   drvrData%OutData%FileDescLines(1)  = 'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//TRIM(GetVersion(version))
   drvrData%OutData%FileDescLines(2)  = 'linked with ' //' '//TRIM(GetNVD(NWTC_Ver            ))  ! we'll get the rest of the linked modules in the section below
   drvrData%OutData%FileDescLines(3)  = 'Description from the driver input file: '//TRIM(drvrData%FTitle)
   
   
   IF (drvrData%WrTxtOutFile) THEN

      call GetNewUnit(drvrData%OutData%unOutFile, ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            drvrData%OutData%unOutFile = -1
            return
         end if
            
      call OpenFOutFile ( drvrData%OutData%unOutFile, trim(drvrData%OutRootName)//'.out', ErrStat2, ErrMsg2)
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
         if ( ErrStat >= AbortErrLev ) then
            drvrData%OutData%unOutFile = -1
            return
         end if

         ! Add some file information:

      WRITE (drvrData%OutData%unOutFile,'(/,A)')  TRIM( drvrData%OutData%FileDescLines(1) )
      WRITE (drvrData%OutData%unOutFile,'(1X,A)') TRIM( drvrData%OutData%FileDescLines(2) )
      WRITE (drvrData%OutData%unOutFile,'()' )    !print a blank line
      WRITE (drvrData%OutData%unOutFile,'(A)'   ) TRIM( drvrData%OutData%FileDescLines(3) )
      WRITE (drvrData%OutData%unOutFile,'()' )    !print a blank line


         !......................................................
         ! Write the names of the output parameters on one line:
         !......................................................
      CALL WrFileNR ( drvrData%OutData%unOutFile, trim(drvrData%OutData%WriteOutputHdr(1)) )
      DO I=2,drvrData%OutData%NumOuts
         CALL WrFileNR ( drvrData%OutData%unOutFile, Delim//trim(drvrData%OutData%WriteOutputHdr(I)) )
      ENDDO ! I

      WRITE (drvrData%OutData%unOutFile,'()')

         !......................................................
         ! Write the units of the output parameters on one line:
         !......................................................
      CALL WrFileNR ( drvrData%OutData%unOutFile, trim(drvrData%OutData%WriteOutputUnt(1)) )
      DO I=2,drvrData%OutData%NumOuts
         CALL WrFileNR ( drvrData%OutData%unOutFile, Delim//trim(drvrData%OutData%WriteOutputUnt(I)) )
      ENDDO ! I

      WRITE (drvrData%OutData%unOutFile,'()')
         
   END IF
   
   IF (drvrData%WrBinOutFile) THEN
      drvrData%OutData%TimeData(1) = 0.0_DbKi                  ! This is the first output time, which we will set later
      drvrData%OutData%TimeData(2) = drvrData%TimeInterval      ! This is the (constant) time between subsequent writes to the output file
   END IF

   
END SUBROUTINE InitOutputFile
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE FillOutputFile(time, y_SeaSt, y_HD, drvrData, ErrStat, ErrMsg)
   REAL(DbKi),                      INTENT( IN    ) :: time
   TYPE(SeaSt_OutputType),          INTENT( IN    ) :: y_SeaSt                 ! SeaState outputs
   TYPE(HydroDyn_OutputType),       INTENT( IN    ) :: y_HD                    ! HydroDyn outputs
   TYPE(HD_Drvr_Data),              INTENT( INOUT ) :: drvrData
   INTEGER,                         INTENT(   OUT ) :: ErrStat                ! returns a non-zero value when an error occurs  
   CHARACTER(*),                    INTENT(   OUT ) :: ErrMsg                 ! Error message if ErrStat /= ErrID_None
   
   character(60)                                    :: TmpStr
   integer(IntKi)                                   :: i, Indx
   integer(IntKi)                                   :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                             :: errMsg2       ! temporary error message 
   character(*), parameter                          :: RoutineName = 'FillOutputFile'
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   IF ( drvrData%OutData%n_Out < drvrData%NSteps ) THEN
      drvrData%OutData%n_Out = drvrData%OutData%n_Out + 1
   ELSE IF (drvrData%WrBinOutFile) THEN
      ErrStat = ErrID_Warn
      ErrMsg = 'Not all data could be written to the binary output file.'
   END IF

   ! Fill data array with concatenated writeOutput data:
   Indx = 0
   do i=1,drvrData%OutData%NumOutsMods(1)
      Indx = Indx + 1
      drvrData%OutData%Storage(Indx, drvrData%OutData%n_Out) = y_SeaSt%WriteOutput(i)
   end do
   do i=1,drvrData%OutData%NumOutsMods(2)
      Indx = Indx + 1
      drvrData%OutData%Storage(Indx, drvrData%OutData%n_Out) = y_HD%WriteOutput(i)
   end do


   IF (drvrData%WrTxtOutFile) THEN
            ! Write one line of tabular output:

            ! time
      WRITE( TmpStr, '('//trim(drvrData%OutData%OutFmt_t)//')' ) time
      CALL WrFileNR( drvrData%OutData%unOutFile, trim(TmpStr) )

         ! write the individual module output (convert to SiKi if necessary, so that we don't need to print so many digits in the exponent)
      CALL WrNumAryFileNR ( drvrData%OutData%unOutFile, REAL(drvrData%OutData%Storage(:,drvrData%OutData%n_Out),SiKi), '"'//Delim//'"'//drvrData%OutData%OutFmt, ErrStat2, ErrMsg2 )
         CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

         ! write a new line (advance to the next line)
      WRITE (drvrData%OutData%unOutFile,'()')
   END IF


   IF (drvrData%WrBinOutFile) THEN
         ! store time data
      IF ( drvrData%OutData%n_Out == 1_IntKi ) THEN
         drvrData%OutData%TimeData(drvrData%OutData%n_Out) = time   ! First time in the output file
      END IF
   END IF


END SUBROUTINE FillOutputFile
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE WriteOutputFile(drvrData, ErrStat, ErrMsg)
   TYPE(HD_Drvr_Data),              INTENT( IN    ) :: drvrData
   INTEGER,                         INTENT(   OUT ) :: ErrStat              ! returns a non-zero value when an error occurs  
   CHARACTER(*),                    INTENT(   OUT ) :: ErrMsg               ! Error message if ErrStat /= ErrID_None
   
   ErrStat = ErrID_None
   ErrMsg = ""

   IF (drvrData%WrTxtOutFile) THEN
      IF (drvrData%OutData%unOutFile > 0) CLOSE(drvrData%OutData%unOutFile)
   END IF
   
   IF (drvrData%WrBinOutFile .AND. drvrData%OutData%n_Out > 0) THEN

      CALL WrBinFAST(TRIM(drvrData%OutRootName)//'.outb', FileFmtID_ChanLen_In, TRIM(drvrData%OutData%FileDescLines(1))//' '//TRIM(drvrData%OutData%FileDescLines(2))//'; '//TRIM(drvrData%OutData%FileDescLines(3)), &
            drvrData%OutData%WriteOutputHdr, drvrData%OutData%WriteOutputUnt, drvrData%OutData%TimeData, drvrData%OutData%Storage(:,1:drvrData%OutData%n_Out), ErrStat, ErrMsg)

   END IF
   
   
END SUBROUTINE WriteOutputFile
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetHDInputs_Constant(u_HD, mappingData, drvrData, ErrStat, ErrMsg)
   TYPE(HydroDyn_InputType),        INTENT( INOUT ) :: u_HD                    ! HydroDyn inputs
   TYPE(HD_Drvr_MappingData),       INTENT( INOUT ) :: mappingData
   TYPE(HD_Drvr_Data),              INTENT( IN    ) :: drvrData
   
   INTEGER,                         INTENT(   OUT ) :: ErrStat                ! returns a non-zero value when an error occurs  
   CHARACTER(*),                    INTENT(   OUT ) :: ErrMsg                 ! Error message if ErrStat /= ErrID_None
   
   integer(IntKi)                                   :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                             :: errMsg2       ! temporary error message 
   character(*), parameter                          :: RoutineName = 'SetHDInputs_Constant'
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   IF (( drvrData%PRPInputsMod /= 2 ) .AND. ( drvrData%PRPInputsMod >= 0 )) THEN
                
      u_HD%PRPMesh%TranslationDisp(:,1)   = drvrData%uPRPInSteady(1:3) 

         ! Compute direction cosine matrix from the rotation angles
      CALL SmllRotTrans( 'InputRotation', drvrData%uPRPInSteady(4), drvrData%uPRPInSteady(5), drvrData%uPRPInSteady(6), u_HD%PRPMesh%Orientation(:,:,1), 'Junk', ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      u_HD%PRPMesh%TranslationVel(:,1)    = drvrData%uDotPRPInSteady(1:3)
      u_HD%PRPMesh%RotationVel(:,1)       = drvrData%uDotPRPInSteady(4:6)
      u_HD%PRPMesh%TranslationAcc(:,1)    = drvrData%uDotDotPRPInSteady(1:3)
      u_HD%PRPMesh%RotationAcc(:,1)       = drvrData%uDotDotPRPInSteady(4:6)
      
         ! Map PRP kinematics to the WAMIT mesh with 1 to NBody nodes
      IF ( u_HD%WAMITMesh%Initialized ) THEN 
         CALL Transfer_Point_to_Point( u_HD%PRPMesh, u_HD%WAMITMesh, mappingData%HD_Ref_2_WB_P, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      END IF
      
         ! Map PRP kinematics to the Morison mesh
      if ( u_HD%Morison%Mesh%Initialized ) then
         CALL Transfer_Point_to_Point( u_HD%PRPMesh, u_HD%Morison%Mesh, mappingData%HD_Ref_2_M_P, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      end if
      
   END IF
   
END SUBROUTINE SetHDInputs_Constant
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE SetHDInputs(time, n, u_HD, mappingData, drvrData, ErrStat, ErrMsg)
   REAL(DbKi),                      INTENT( IN    ) :: time
   INTEGER(IntKi),                  INTENT( IN    ) :: n
   TYPE(HydroDyn_InputType),        INTENT( INOUT ) :: u_HD                    ! HydroDyn inputs
   TYPE(HD_Drvr_MappingData),       INTENT( INOUT ) :: mappingData
   TYPE(HD_Drvr_Data),              INTENT( IN    ) :: drvrData
   
   INTEGER,                         INTENT(   OUT ) :: ErrStat                ! returns a non-zero value when an error occurs  
   CHARACTER(*),                    INTENT(   OUT ) :: ErrMsg                 ! Error message if ErrStat /= ErrID_None
   
   integer(IntKi)                                   :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                             :: errMsg2       ! temporary error message 
   character(*), parameter                          :: RoutineName = 'SetHDInputs_Constant'
   real(ReKi)                                       :: yInterp(size(drvrData%PRPin,2))
   integer(intKi)                                   :: indxHigh, indxMid, indxLow
   integer(intKi)                                   :: i
   
   ErrStat = ErrID_None
   ErrMsg = ""

   ! PRPInputsMod 2: Reads time series of positions, velocities, and accelerations for the platform reference point
   IF ( drvrData%PRPInputsMod == 2 ) THEN
      call InterpStpMat( real(time,ReKi), drvrData%PRPinTime, drvrData%PRPin, mappingData%Ind, size(drvrData%PRPinTime), yInterp )
      
      u_HD%PRPMesh%TranslationDisp(:,1)   = yInterp(1:3) 

         ! Compute direction cosine matrix from the rotation angles
               
!         maxAngle = max( maxAngle, abs(yInterp(4:6)) )
            
      CALL SmllRotTrans( 'InputRotation', yInterp(4), yInterp(5), yInterp(6), u_HD%PRPMesh%Orientation(:,:,1), 'Junk', ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

      u_HD%PRPMesh%TranslationVel(:,1)    = yInterp( 7: 9)
      u_HD%PRPMesh%RotationVel(:,1)       = yInterp(10:12)
      u_HD%PRPMesh%TranslationAcc(:,1)    = yInterp(13:15)
      u_HD%PRPMesh%RotationAcc(:,1)       = yInterp(16:18)
            
      IF ( u_HD%WAMITMesh%Initialized ) THEN
            ! Map kinematics to the WAMIT mesh with 1 to NBody nodes
         CALL Transfer_Point_to_Point( u_HD%PRPMesh, u_HD%WAMITMesh, mappingData%HD_Ref_2_WB_P, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      END IF
         
      IF ( u_HD%Morison%Mesh%Initialized ) THEN
            ! Map kinematics to the WAMIT mesh with 1 to NBody nodes
         CALL Transfer_Point_to_Point( u_HD%PRPMesh, u_HD%Morison%Mesh, mappingData%HD_Ref_2_M_P, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      END IF
          
   ELSEIF ( drvrData%PRPInputsMod < 0 ) THEN
      
      !@mhall: new kinematics input for moving bodies individually
      ! PRPInputsMod < 0: Reads time series of positions for each body individually, and uses finite differences to also get velocities and accelerations.
      ! The number of bodies is the negative of PRPInputsMod.
      
      i = min(n,drvrData%NSteps)
      if (n <= drvrData%NSteps .and. .not. EqualRealNos( REAL(time,ReKi), drvrData%PRPinTime(i) ) ) then
         call SetErrStat(ErrID_Fatal, 'time does not match PRP input file data', ErrStat, ErrMsg, RoutineName)
         return
      end if
               
      ! platform reference point (PRP), and body 1-NBody displacements
      u_HD%PRPMesh%TranslationDisp(:,1)   = drvrData%PRPin(n,1:3) 
      DO I=1,drvrData%NBody
         u_HD%WAMITMesh%TranslationDisp(:,I)   = drvrData%PRPin(n, 6*I+1:6*I+3) 
      END DO
               
      ! PRP and body 1-NBody orientations (skipping the maxAngle stuff)
      CALL SmllRotTrans(    'InputRotation', drvrData%PRPin(n,    4), drvrData%PRPin(n,    5), drvrData%PRPin(n,    6), u_HD%PRPMesh%Orientation(:,:,1), 'PRP orientation', ErrStat2, ErrMsg2 )
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
      DO I=1, drvrData%NBody
         CALL SmllRotTrans( 'InputRotation', drvrData%PRPin(n,6*I+4), drvrData%PRPin(n,6*I+5), drvrData%PRPin(n,6*I+6), u_HD%WAMITMesh%Orientation(:,:,I), 'body orientation', ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      END DO

      ! use finite differences for velocities and accelerations
      IF (n == 1) THEN   ! use forward differences for first time step
         indxHigh = n+2
         indxMid  = n+1
         indxLow  = n
            
         u_HD%PRPMesh%TranslationVel(:,1) = (drvrData%PRPin(indxMid , 1:3) -   drvrData%PRPin(indxLow , 1:3))/drvrData%TimeInterval
         u_HD%PRPMesh%RotationVel(   :,1) = (drvrData%PRPin(indxMid , 4:6) -   drvrData%PRPin(indxLow , 4:6))/drvrData%TimeInterval
               
         DO I=1,drvrData%NBody
            u_HD%WAMITMesh%TranslationVel(:,I) = (drvrData%PRPin(indxMid,  6*I+1:6*I+3) -   drvrData%PRPin(indxLow, 6*I+1:6*I+3))/drvrData%TimeInterval
            u_HD%WAMITMesh%RotationVel(   :,I) = (drvrData%PRPin(indxMid,  6*I+4:6*I+6) -   drvrData%PRPin(indxLow, 6*I+4:6*I+6))/drvrData%TimeInterval
         END DO

      ELSE IF (n >= drvrData%NSteps) THEN  ! use backward differences for last time step
         indxHigh = n
         indxMid  = n-1
         indxLow  = n-2
            
         u_HD%PRPMesh%TranslationVel(:,1) = (drvrData%PRPin(indxHigh, 1:3) -   drvrData%PRPin(indxMid, 1:3))/drvrData%TimeInterval
         u_HD%PRPMesh%RotationVel(   :,1) = (drvrData%PRPin(indxHigh, 4:6) -   drvrData%PRPin(indxMid, 4:6))/drvrData%TimeInterval
               
         DO I=1,drvrData%NBody
            u_HD%WAMITMesh%TranslationVel(:,I) = (drvrData%PRPin(indxHigh, 6*I+1:6*I+3) -   drvrData%PRPin(indxMid, 6*I+1:6*I+3))/drvrData%TimeInterval
            u_HD%WAMITMesh%RotationVel(   :,I) = (drvrData%PRPin(indxHigh, 6*I+4:6*I+6) -   drvrData%PRPin(indxMid, 6*I+4:6*I+6))/drvrData%TimeInterval
         END DO
            
      ELSE   ! otherwise use central differences for intermediate time steps
         indxHigh = n+1
         indxMid  = n
         indxLow  = n -1
                     
         u_HD%PRPMesh%TranslationVel(:,1) = (drvrData%PRPin(indxHigh, 1:3) - drvrData%PRPin(indxLow, 1:3))*0.5/drvrData%TimeInterval
         u_HD%PRPMesh%RotationVel(   :,1) = (drvrData%PRPin(indxHigh, 4:6) - drvrData%PRPin(indxLow, 4:6))*0.5/drvrData%TimeInterval
               
         DO I=1,drvrData%NBody
            u_HD%WAMITMesh%TranslationVel(:,I) = (drvrData%PRPin(indxHigh, 6*I+1:6*I+3) - drvrData%PRPin(indxLow, 6*I+1:6*I+3))*0.5/drvrData%TimeInterval
            u_HD%WAMITMesh%RotationVel(   :,I) = (drvrData%PRPin(indxHigh, 6*I+4:6*I+6) - drvrData%PRPin(indxLow, 6*I+4:6*I+6))*0.5/drvrData%TimeInterval
         END DO
               
      END IF
            
      ! calculate accelerations based on displacements:
      u_HD%PRPMesh%TranslationAcc(:,1)      = (drvrData%PRPin(indxHigh, 1:3)         - 2*drvrData%PRPin(indxMid, 1:3)         + drvrData%PRPin(indxLow, 1:3))        /(drvrData%TimeInterval**2)
      u_HD%PRPMesh%RotationAcc(   :,1)      = (drvrData%PRPin(indxHigh, 4:6)         - 2*drvrData%PRPin(indxMid, 4:6)         + drvrData%PRPin(indxLow, 4:6))        /(drvrData%TimeInterval**2)

      DO I=1,drvrData%NBody
         u_HD%WAMITMesh%TranslationAcc(:,I) = (drvrData%PRPin(indxHigh, 6*I+1:6*I+3) - 2*drvrData%PRPin(indxMid, 6*I+1:6*I+3) + drvrData%PRPin(indxLow, 6*I+1:6*I+3))/(drvrData%TimeInterval**2)
         u_HD%WAMITMesh%RotationAcc(   :,I) = (drvrData%PRPin(indxHigh, 6*I+4:6*I+6) - 2*drvrData%PRPin(indxMid, 6*I+4:6*I+6) + drvrData%PRPin(indxLow, 6*I+4:6*I+6))/(drvrData%TimeInterval**2)
      END DO
               
            
      IF ( u_HD%Morison%Mesh%Initialized ) THEN
         ! Map kinematics to the WAMIT mesh with 1 to NBody nodes
         CALL Transfer_Point_to_Point( u_HD%PRPMesh, u_HD%Morison%Mesh, mappingData%HD_Ref_2_M_P, ErrStat2, ErrMsg2 )
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      END IF

   ELSE
      ! constant inputs are not recalculated at each time step. Instead this is called at initialization
      ! CALL SetHDInputs_Constant()
   END IF
      
END SUBROUTINE
!----------------------------------------------------------------------------------------------------------------------------------
!>
subroutine CreatePointMesh(mesh, posInit, orientInit, HasMotion, HasLoads, errStat, errMsg)
   type(MeshType),               intent(inout) :: mesh
   real(ReKi),                   intent(in   ) :: PosInit(3)                                             !< Xi,Yi,Zi, coordinates of node
   real(R8Ki),                   intent(in   ) :: orientInit(3,3)                                        !< Orientation (direction cosine matrix) of node; identity by default
   logical,                      intent(in   ) :: HasMotion   !< include displacements in mesh
   logical,                      intent(in   ) :: HasLoads   !< include loads in mesh
   integer(IntKi)              , intent(out)   :: errStat       ! Status of error message
   character(*)                , intent(out)   :: errMsg        ! Error message if ErrStat /= ErrID_None
   integer(IntKi)       :: errStat2      ! local status of error message
   character(ErrMsgLen) :: errMsg2       ! local error message if ErrStat /= ErrID_None
   errStat = ErrID_None
   errMsg  = ''

   call MeshCreate(mesh, COMPONENT_INPUT, 1, errStat2, errMsg2,  &
      Orientation=HasMotion, TranslationDisp=HasMotion, TranslationVel=HasMotion, RotationVel=HasMotion, TranslationAcc=HasMotion, RotationAcc=HasMotion, &
      Force = HasLoads, Moment = HasLoads)
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'CreatePointMesh')
   if (ErrStat >= AbortErrLev) return

   call MeshPositionNode(mesh, 1, posInit, errStat2, errMsg2, orientInit); 
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'CreatePointMesh')

   call MeshConstructElement(mesh, ELEMENT_POINT, errStat2, errMsg2, p1=1); 
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'CreatePointMesh')

   call MeshCommit(mesh, errStat2, errMsg2);
   call SetErrStat(errStat2, errMsg2, errStat, errMsg, 'CreatePointMesh')

end subroutine CreatePointMesh
!----------------------------------------------------------------------------------------------------------------------------------
!> Compute Rigid body loads at the PRP, after a perturbation of the PRP
SUBROUTINE PRP_CalcOutput(t, u, p, x, xd, z, OtherState, y, m, EDRPMesh, Loads, mappingData, ErrStat, ErrMsg)
   TYPE(MeshType)          ,             INTENT(INOUT) :: EDRPMesh !<
   
   REAL(DbKi),                         INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(HydroDyn_InputType),           INTENT(INOUT)  :: u           !< Inputs at Time (note that this is intent out because we're copying the u%WAMITMesh into m%u_wamit%mesh)
   TYPE(HydroDyn_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(HydroDyn_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at Time
   TYPE(HydroDyn_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at Time
   TYPE(HydroDyn_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at Time
   TYPE(HydroDyn_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at Time
   TYPE(HydroDyn_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at Time (Input only so that mesh con-
                                                                     !!   nectivity information does not have to be recalculated)
   TYPE(HydroDyn_MiscVarType),         INTENT(INOUT)  :: m           !< Initial misc/optimization variables           
   INTEGER(IntKi),                     INTENT(  OUT)  :: ErrStat     !! Error status of the operation
   CHARACTER(*),                       INTENT(  OUT)  :: ErrMsg      !! Error message if ErrStat /= ErrID_None
   
   Real(ReKi)               ,            INTENT(OUT)   :: Loads(18) !< Loads at PRP and EDRP
   TYPE(HD_Drvr_MappingData),            INTENT(INOUT) :: mappingData
   
   INTEGER(IntKi)                                     :: ErrStat2     ! Status of error message
   CHARACTER(ErrMsgLen)                               :: ErrMsg2       ! Error message if ErrStat /= ErrID_None
   CHARACTER(*), PARAMETER                            :: RoutineName = 'PRP_CalcOutput'
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   call HydroDyn_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat2, ErrMsg2 )
   call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   ! Integrate all the mesh loads onto the platfrom reference Point (PRP) at (0,0,0)
   Loads(1:6) = m%F_Hydro ! NOTE this is mapped to PRP using m%AllHdroOrigin 

   ! --- Transfer loads from HydroOrigin to EDRPMesh
   call Transfer_Point_to_Point( m%AllHdroOrigin, EDRPMesh, mappingData%HD_RefLoads_2_ED_Ref, ErrStat2, ErrMsg2, u%PRPMesh, EDRPMesh )
   call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   Loads(7:9)   = EDRPMesh%Force(:,1)
   Loads(10:12) = EDRPMesh%Moment(:,1)

   ! --- Transfer loads from HydroOrigin to (0,0,0)
   call Transfer_Point_to_Point( m%AllHdroOrigin, mappingData%ZZZPtMeshLoads, mappingData%HD_RefLoads_2_ZZZLoads, ErrStat2, ErrMsg2, u%PRPMesh, mappingData%ZZZPtMeshMotion )
   call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,RoutineName)
   Loads(13:15) = mappingData%ZZZPtMeshLoads%Force(:,1)
   Loads(16:18) = mappingData%ZZZPtMeshLoads%Moment(:,1)

   !print*,'LoadsPRP',Loads(1:6)
   !print*,'LoadsEDP',Loads(7:12)
   !print*,'Loads000',Loads(13:18)

END SUBROUTINE PRP_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> Pertub the "PRP" inputs and trigger the rigid body motion on the other HydroDyn meshes
SUBROUTINE PRP_Perturb_u( n, perturb_sign, p, u, EDRPMesh, du, Motion_HDRP, mappingData, ErrStat, ErrMsg)
   INTEGER( IntKi )                    , INTENT(IN   ) :: n                      !< number of array element to use 
   INTEGER( IntKi )                    , INTENT(IN   ) :: perturb_sign           !< +1 or -1 (value to multiply perturbation by; positive or negative difference)
   TYPE(HydroDyn_ParameterType),         INTENT(IN   )  :: p           !< Parameters
   TYPE(HydroDyn_InputType), target    , INTENT(INOUT) :: u                      !< perturbed HD inputs
   TYPE(MeshType)          , target    , INTENT(INOUT) :: EDRPMesh !<
   REAL( R8Ki )                        , INTENT(  OUT) :: du                     !< amount that specific input was perturbed
   logical                             , INTENT(IN   ) :: Motion_HDRP   !< If True, perturb the PRP otherwise perturb the EDRP for motion
   TYPE(HD_Drvr_MappingData),            INTENT(INOUT) :: mappingData
   integer(IntKi)                      , intent(  out) :: errStat       ! Status of error message
   character(*)                        , intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   
   type(MeshType), pointer :: pointMesh !Alias

   ! local variables
   integer                                             :: fieldType ! 1=TranslationDisp, 2=Orientation, 3=TranslationVel etc. 6
   integer                                             :: fieldIndx
   integer                                             :: fieldIndx6
   integer , parameter                                 :: node =1
   Real(R8Ki)   perturb_t, perturb
!   REAL(R8Ki) :: dcm (3,3)            ! The resulting transformation matrix from X to x, (-).
!   Real(R8Ki) :: theta(3)
   INTEGER(IntKi)                                      :: ErrStat2     ! Status of error message
   CHARACTER(ErrMsgLen)                                :: ErrMsg2       ! Error message if ErrStat /= ErrID_None
   character(*), parameter                             :: RoutineName = 'PRP_Perturb_u'
   
   ErrStat = ErrID_None
   ErrMsg = ""

   ! From "n" to: field type, axis, variable
   fieldType = int((n-1)/3)+1  ! 1=TranslationDisp, 2=Orientation, 3=TranslationVel etc. 6
   fieldIndx = mod(n-1,3)+1    ! 1=x, 2=y 3=z (axis)
   fieldIndx6= mod(n-1,6)+1    ! 1=x, 2=y 3=z 4=theta_x, 5=theta_y 3=theta_z (variable)

   ! Perturbation amplitude
   perturb_t = 0.02_ReKi*D2R * max(p%WtrDpth,1.0_ReKi) ! translation input scaling  
   perturb   = 2*D2R                 ! rotational input scaling
   !perturb_t = 1.0
   !perturb   = 0.1
   if (fieldIndx6<=3) then
     du = perturb_t    ! TranslationDisp,TranslationVel, TranslationAcc
   elseif (fieldIndx<=6) then
     du = perturb      ! Orientation, TranslationVel
   else
      call SetErrStat(ErrID_Fatal, 'Wrong field index', ErrStat, ErrMsg, RoutineName)
      return
   endif

   if (Motion_HDRP) then
      pointMesh => u%PRPMesh
   else
      pointMesh => EDRPMesh
   endif

   ! --- Perturbing the point mesh
   !print*,''
   !print*,'Perturb',n, perturb_sign
   SELECT CASE(fieldType)      
      CASE ( 1) !Module/Mesh/Field: u%PRPMesh%TranslationDisp = 1     
         pointMesh%TranslationDisp (fieldIndx,node) = pointMesh%TranslationDisp (fieldIndx,node) + du * perturb_sign       
      CASE ( 2) !Module/Mesh/Field: u%PRPMesh%Orientation = 2
         CALL PerturbOrientationMatrix( pointMesh%Orientation(:,:,node), du * perturb_sign, fieldIndx, UseSmlAngle=.true. )
      CASE ( 3) !Module/Mesh/Field: u%PRPMesh%TranslationVel = 3
         pointMesh%TranslationVel( fieldIndx,node) = pointMesh%TranslationVel( fieldIndx,node) + du * perturb_sign         
      CASE ( 4) !Module/Mesh/Field: u%PRPMesh%RotationVel = 4
         pointMesh%RotationVel (fieldIndx,node) = pointMesh%RotationVel (fieldIndx,node) + du * perturb_sign               
      CASE ( 5) !Module/Mesh/Field: u%PRPMesh%TranslationAcc = 5
         pointMesh%TranslationAcc( fieldIndx,node) = pointMesh%TranslationAcc( fieldIndx,node) + du * perturb_sign       
      CASE ( 6) !Module/Mesh/Field: u%PRPMesh%RotationAcc = 6
         pointMesh%RotationAcc(fieldIndx,node) = pointMesh%RotationAcc(fieldIndx,node) + du * perturb_sign               
      CASE default
         call SetErrStat(ErrID_Fatal, 'Wrong fieldType', ErrStat, ErrMsg, RoutineName)
   END SELECT   

   ! --- Trigger ED->PRP or PRP->ED
   if (Motion_HDRP) then
      ! PRP->ED
      call Transfer_Point_to_Point( pointMesh, EDRPMesh, mappingData%HD_Ref_2_ED_Ref, ErrStat2, ErrMsg2 );
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   else
      ! ED->PRP
      call Transfer_Point_to_Point( pointMesh, u%PRPMesh, mappingData%ED_Ref_2_HD_Ref, ErrStat2, ErrMsg2 );
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
      !print*,'--------------------------------------------  EDRP  -------------------------------------'
      !print*,'--------------------------------------------  EDRP  -------------------------------------'
      !print*,'--------------------------------------------  EDRP  -------------------------------------'
      !call MeshPrintInfo (6, EDRPMesh)
      !print*,''
      !print*,''
      !print*,''
      !print*,'--------------------------------------------  PRP  -------------------------------------'
      !print*,'--------------------------------------------  PRP  -------------------------------------'
      !print*,'--------------------------------------------  PRP  -------------------------------------'
      !call MeshPrintInfo (6, u%PRPMesh)
   endif

END SUBROUTINE PRP_Perturb_u
!----------------------------------------------------------------------------------------------------------------------------------
!> Calculate the partial derivative of the output functions (Y) with respect to the inputs (u)
SUBROUTINE PRP_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, dYdu, Motion_HDRP, mappingData, ErrStat, ErrMsg)
   REAL(DbKi),                                 INTENT(IN   ) :: t          !< Time in seconds at operating point
   TYPE(HydroDyn_InputType),                   INTENT(INOUT) :: u           !< Inputs at Time (note that this is intent out because we're copying the u%WAMITMesh into m%u_wamit%mesh)
   TYPE(HydroDyn_ParameterType),               INTENT(IN   ) :: p           !< Parameters
   TYPE(HydroDyn_ContinuousStateType),         INTENT(IN   ) :: x           !< Continuous states at Time
   TYPE(HydroDyn_DiscreteStateType),           INTENT(IN   ) :: xd          !< Discrete states at Time
   TYPE(HydroDyn_ConstraintStateType),         INTENT(IN   ) :: z           !< Constraint states at Time
   TYPE(HydroDyn_OtherStateType),              INTENT(IN   ) :: OtherState  !< Other states at Time
   TYPE(HydroDyn_OutputType),                  INTENT(INOUT) :: y           !< Outputs computed at Time (Input only so that mesh con-
                                                                            !!   nectivity information does not have to be recalculated)
   TYPE(HydroDyn_MiscVarType),                 INTENT(INOUT) :: m           !< Initial misc/optimization variables           
   INTEGER(IntKi),                             INTENT(  OUT) :: ErrStat    !< Error status of the operation
   CHARACTER(*),                               INTENT(  OUT) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   REAL(R8Ki), ALLOCATABLE, OPTIONAL,          INTENT(INOUT) :: dYdu(:,:)  !< Partial derivatives of output functions (Y) with respect
   logical,                                    INTENT(IN   ) :: Motion_HDRP   !< If True, perturb the PRP otherwise perturb the EDRP for motion
   TYPE(HD_Drvr_MappingData),                  INTENT(INOUT) :: mappingData
   
   ! local variables
   TYPE(HydroDyn_OutputType)                                 :: y_tmp
   TYPE(HydroDyn_InputType)                                  :: u_perturb
   TYPE(MeshType)                                            :: EDRPtMesh_perturb
   Real(ReKi)                                                :: Loads_p(18)
   Real(ReKi)                                                :: Loads_m(18)
   REAL(R8Ki)                                                :: delta        ! delta change in input or state
   integer(IntKi)                                            :: i
   INTEGER(IntKi)                                            :: ErrStat2
   CHARACTER(ErrMsgLen)                                      :: ErrMsg2
   CHARACTER(*), PARAMETER                                   :: RoutineName = 'PRP_JacobianPInput'
   
   ErrStat = ErrID_None
   ErrMsg  = ''
   
    ! make a copy of the inputs to perturb
    call HydroDyn_CopyInput( u, u_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2);      call SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg,RoutineName)
    call MeshCopy(mappingData%EDRPtMesh, EDRPtMesh_perturb, MESH_NEWCOPY, ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg,RoutineName)
    
    ! allocate dYdu if necessary
    if (.not. allocated(dYdu)) then
       call AllocAry(dYdu, size(Loads_p), 18, 'dYdu', ErrStat2, ErrMsg2)
       call SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg,RoutineName)
       !if (ErrStat >= AbortErrLev) call cleanup() 
       dYdu=0.0_ReKi
    endif
    ! make a copy of outputs because we will need two for the central difference computations (with orientations)
    call HydroDyn_CopyOutput( y, y_tmp, MESH_NEWCOPY, ErrStat2, ErrMsg2);

    do i=1,size(dYdu,2)
       ! get u_op + delta u
       call HydroDyn_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg,RoutineName)
       call MeshCopy(mappingData%EDRPtMesh, EDRPtMesh_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg,RoutineName)
       call PRP_Perturb_u(i, 1, p, u_perturb, EDRPtMesh_perturb, delta, Motion_HDRP, mappingData, ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg,RoutineName)
       ! compute y at u_op + delta u
       call PRP_CalcOutput( t, u_perturb, p, x, xd, z, OtherState, y_tmp, m, EDRPtMesh_perturb, Loads_p, mappingData, ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg,RoutineName)

       ! get u_op - delta u
       call HydroDyn_CopyInput( u, u_perturb, MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg,RoutineName)
       call MeshCopy(mappingData%EDRPtMesh, EDRPtMesh_perturb,  MESH_UPDATECOPY, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg,RoutineName)
       call PRP_Perturb_u( i, -1, p, u_perturb, EDRPtMesh_perturb, delta , Motion_HDRP, mappingData, ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg,RoutineName)
       ! compute y at u_op - delta u
       call PRP_CalcOutput( t, u_perturb, p, x, xd, z, OtherState, y_tmp, m, EDRPtMesh_perturb, Loads_m, mappingData, ErrStat2, ErrMsg2); call SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg,RoutineName)

       ! get central difference:            
       dYdu(:,i) = (Loads_p-Loads_m) / (2.0_R8Ki*delta)
       !if(i==4) STOP
    end do

    call HydroDyn_DestroyOutput(      y_tmp, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg,RoutineName)
    call HydroDyn_DestroyInput (  u_perturb, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg,RoutineName)
    call MeshDestroy(EDRPtMesh_perturb, ErrStat2, ErrMsg2 ); call SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg,RoutineName)
   
END SUBROUTINE PRP_JacobianPInput
!----------------------------------------------------------------------------------------------------------------------------------
! --- Rigid body Linearization at t=0
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Linearization(t, u, p, x, xd, z, OtherState, y, m, Motion_HDRP, mappingData, ErrStat, ErrMsg)
   REAL(DbKi),                                 INTENT(IN   ) :: t             !< Time in seconds at operating point
   TYPE(HydroDyn_InputType),                   INTENT(INOUT) :: u             !< Inputs at Time (note that this is intent out because we're copying the u%WAMITMesh into m%u_wamit%mesh)
   TYPE(HydroDyn_ParameterType),               INTENT(IN   ) :: p             !< Parameters
   TYPE(HydroDyn_ContinuousStateType),         INTENT(IN   ) :: x             !< Continuous states at Time
   TYPE(HydroDyn_DiscreteStateType),           INTENT(IN   ) :: xd            !< Discrete states at Time
   TYPE(HydroDyn_ConstraintStateType),         INTENT(IN   ) :: z             !< Constraint states at Time
   TYPE(HydroDyn_OtherStateType),              INTENT(IN   ) :: OtherState    !< Other states at Time
   TYPE(HydroDyn_OutputType),                  INTENT(INOUT) :: y             !< Outputs computed at Time (Input only so that mesh con-
                                                                              !!   nectivity information does not have to be recalculated)
   TYPE(HydroDyn_MiscVarType),                 INTENT(INOUT) :: m             !< Initial misc/optimization variables           
   logical   ,                                 INTENT(IN   ) :: Motion_HDRP   !< If True, perturb the PRP otherwise perturb the EDRP for motion
   TYPE(HD_Drvr_MappingData),                  INTENT(INOUT) :: mappingData
   integer(IntKi)              ,               intent(  out) :: errStat       ! Status of error message
   character(*)                ,               intent(  out) :: errMsg        ! Error message if ErrStat /= ErrID_None
   
   real(R8Ki), allocatable, dimension(:,:) :: dYdu
   integer :: i,j
   character(40):: sMotion
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   !print*,'>>>> Linearize', drvrData%PtfmRefzt
   if (Motion_HDRP) then
      sMotion ='motions at PRP'
   else
      sMotion ='motions at EDRP'
   endif
   
   print'(A,F13.6,A)','   Performing rigid-body linearization at t=',t,' with '//trim(sMotion)
   call PRP_JacobianPInput( t, u, p, x, xd, z, OtherState, y, m, dYdu, Motion_HDRP, mappingData, ErrStat, ErrMsg)

   do i=1,size(dYdu,1)
      do j=1,size(dYdu,2)
         if(abs(dYdu(i,j))<1e-5) then
            dYdu(i,j)=0.0_ReKi
         endif
      enddo
   enddo

   call WrMatrix( dYdu( 1: 6,  1: 6), CU, 'F13.6', 'K: (Loads at PRP, '//trim(sMotion)//')' )
   call WrMatrix( dYdu( 7:12,  1: 6), CU, 'F13.6', 'K: (Loads at EDRP, '//trim(sMotion)//')' )
   call WrMatrix( dYdu(13:18,  1: 6), CU, 'F13.6', 'K: (Loads at 0,0,0, fixed, '//trim(sMotion)//')' )
   call WrMatrix( dYdu( 1: 6,  7:12), CU, 'F13.6', 'C:' )
   call WrMatrix( dYdu( 1: 6, 13:18), CU, 'F13.6', 'M:' )

END SUBROUTINE LINEARIZATION
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE HydroDynDriverSubs

