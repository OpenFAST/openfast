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
      type(MeshMapType)                :: HD_Ref_2_WB_P                           ! Mesh mapping between HD Reference pt mesh and WAMIT body(ies) mesh
      type(MeshMapType)                :: HD_Ref_2_M_P                            ! Mesh mapping between HD Reference pt mesh and Morison mesh
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

   integer(IntKi)                                   :: n
   integer(IntKi)                                   :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                             :: errMsg2       ! temporary error message 
   character(*), parameter                          :: RoutineName = 'ReadDriverInputFile'
   
   
      ! Initialize the echo file unit to -1 which is the default to prevent echoing, we will alter this based on user input
   UnEchoLocal = -1
   UnIn = -1
   
   ErrStat = ErrID_None
   ErrMsg = ""
      
   drvrData%NBODY= 0
   
   IF ( drvrData%PRPInputsMod == 2 ) THEN
   
      CALL AllocAry(drvrData%PRPin, drvrData%NSteps, 19, 'PRPin', ErrStat2, ErrMsg2)
         if (Failed()) return
      
   ELSEIF ( drvrData%PRPInputsMod < 0 ) THEN
      ! multi-body kinematics driver option (time, PRP DOFs 1-6, body1 DOFs 1-6, body2 DOFs 1-6...)
      
      drvrData%NBODY = -drvrData%PRPInputsMod
       
      CALL AllocAry(drvrData%PRPin, drvrData%NSteps, 7+6*drvrData%NBODY, 'PRPin', ErrStat2, ErrMsg2)
         if (Failed()) return
      
      call WrScr( 'NBody is '//trim(Num2LStr(drvrData%NBody))//' and planning to read in  '//trim(Num2LStr(SIZE(drvrData%PRPin,2)))//' columns from the input file' )
      
   ELSE
   
      RETURN
      
   END IF
   
      ! Open the (PRP or WAMIT) inputs data file
   CALL GetNewUnit( UnIn ) 
   CALL OpenFInpFile ( UnIn, trim(drvrData%PRPInputsFile), ErrStat2, ErrMsg2 )
      if (Failed()) return
   
      !seems like it would be more efficient to switch the indices on drvrData%PRPin
   DO n = 1,drvrData%NSteps
      CALL ReadAry ( UnIn, drvrData%PRPInputsFile, drvrData%PRPin(n,:), SIZE(drvrData%PRPin,2), 'Line', 'drvrData%PRPin', ErrStat2, ErrMsg2, UnEchoLocal )
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
!----------------------------------------------------------------------------------------------------------------------------------
END MODULE HydroDynDriverSubs

