!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
! Copyright (C) 2016-2017  Envision Energy USA, LTD
!
!    This file is part of the NWTC Subroutine Library.
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
module BeamDyn_driver_subs
  
   USE BeamDyn
   USE BeamDyn_Subs

   IMPLICIT NONE
   
  ! Variables for multi-point loads
   TYPE , PUBLIC :: BD_DriverInternalType
      REAL(ReKi)     , DIMENSION(1:6)               :: DistrLoad        !< Constant distributed load along beam axis, 3 forces and 3 moments [-]
      REAL(ReKi)     , DIMENSION(1:6)               :: TipLoad          !< Constant point load applied at tip, 3 forces and 3 moments [-]
      INTEGER(IntKi)                                :: NumPointLoads    !< Number of constant point loads applied along beam axis, 3 forces and 3 moments, from the driver input file [-]
      REAL(BDKi) ,     DIMENSION(:,:), ALLOCATABLE  :: MultiPointLoad   !< Constant point loads applied along beam axis, size (NumPointLoads,7); (dimension 2: index 1=Relative position along blade span; indices 2-7 = Fx, Fy, Fz, Mx, My, Mz) [-]
            
      TYPE(MeshType)                                :: mplMotion        ! Mesh for blade motion at multipoint loads locations
      TYPE(MeshType)                                :: mplLoads         ! Mesh for multipoint loads
      TYPE(MeshMapType)                             :: Map_BldMotion_to_mplMotion
      TYPE(MeshMapType)                             :: Map_mplLoads_to_PointLoad
      TYPE(MeshType)                                :: y_BldMotion_at_u_point ! Intermediate mesh to transfer motion from output mesh to input mesh
      TYPE(MeshMapType)                             :: Map_y_BldMotion_to_u_point
                                                    
      TYPE(MeshType)                                :: RotationCenter
      TYPE(MeshMapType)                             :: Map_RotationCenter_to_RootMotion

      LOGICAL                                       :: DynamicSolve 
      LOGICAL                                       :: GlbRotBladeT0    ! Initial blade root orientation is also the GlbRot reference frame
      REAL(DbKi)                                    :: RootRelInit(3,3) ! Initial root orientation relative to GlbRot
      REAL(DbKi)                                    :: t_initial
      REAL(DbKi)                                    :: t_final
      REAL(R8Ki)                                    :: w                ! magnitude of rotational velocity vector

      INTEGER(IntKi)                                :: WrVTK            ! VTK visualization data output: (switch) {0=none; 1=init; 2=animation}
      INTEGER(IntKi)                                :: VTK_fps          ! Frame rate for VTK output (frames per second) {will use closest integer multiple of DT} [used only if WrVTK=2]
      INTEGER(IntKi)                                :: n_VTKTime        ! Number of time steps between writing VTK files
      INTEGER(IntKi)                                :: VTK_tWidth       ! number of digits in the time part of file name
      character(1024)                               :: VTK_OutFileRoot  ! rootname for the output file
      
   END TYPE
   
   
   
  contains
  
!------------------------------------------------------------------------------------
!> This routine reads in the primary BeamDyn input file and places the values it reads
!! in the InputFileData structure.
!!   It opens an echo file if requested and returns the (still-open) echo file to the
!!     calling routine.
!!   It also returns the names of the BldFile, FurlFile, and TrwFile for further
!!     reading of inputs.
!------------------------------------------------------------------------------------
  SUBROUTINE BD_ReadDvrFile(DvrInputFile,dt,InitInputData,DvrData,&
                          ErrStat,ErrMsg)

   ! Passed variables
   CHARACTER(*),                 INTENT(IN   ) :: DvrInputFile
   INTEGER(IntKi),               INTENT(  OUT) :: ErrStat
   CHARACTER(*),                 INTENT(  OUT) :: ErrMsg
   TYPE(BD_InitInputType),       INTENT(  OUT) :: InitInputData
   TYPE(BD_DriverInternalType),  INTENT(  OUT) :: DvrData
   REAL(DbKi),                   INTENT(  OUT) :: dt

   ! Local variables:
   REAL(BDKi)                   :: TmpReAry(7)
   INTEGER(IntKi)               :: UnIn                         ! Unit number for reading file
   INTEGER(IntKi)               :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)         :: ErrMsg2                      ! Temporary Error message
   character(*), parameter      :: RoutineName = 'BD_ReadDvrFile'
   character(1024)              :: line
   INTEGER(IntKi)               :: UnEc
   
   CHARACTER(1024)              :: FTitle                       ! "File Title": the 2nd line of the input file, which contains a description of its contents
   CHARACTER(1024)              :: PriPath                      ! Path name of the primary file

   INTEGER(IntKi)               :: i
   INTEGER(IntKi)               :: IOS
!------------------------------------------------------------------------------------

   ! Initialize some variables:
   ErrStat = ErrID_None
   ErrMsg  = ""
   UnEc = -1
   
   CALL GetNewUnit(UnIn,ErrStat2,ErrMsg2);   if (Failed())  return;
   CALL OpenFInpFile(UnIn,DvrInputFile,ErrStat2,ErrMsg2);   if (Failed())  return;
      
      
   CALL GetPath( DvrInputFile, PriPath )     ! Input files will be relative to the path where the primary input file is located.
      
   !-------------------------- HEADER ---------------------------------------------
   CALL ReadCom(UnIn,DvrInputFile,'File Header: Module Version (line 1)',ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadStr(UnIn,DvrInputFile,FTitle,'FTitle','File Header: File Description (line 2)',ErrStat2, ErrMsg2, UnEc);   if (Failed())  return;

   !---------------------- SIMULATION CONTROL --------------------------------------
   CALL ReadCom(UnIn,DvrInputFile,'Section Header: Simulation Control',ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%DynamicSolve,'DynamicSolve','Use Dynamic solve (false for static solve).',ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%t_initial,'t_initial','Starting time of simulation',ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%t_final,"t_final", "Ending time of simulation",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,dt,"dt", "Time increment size",ErrStat2,ErrMsg2,UnEc)
      
   !---------------------- GRAVITY PARAMETER --------------------------------------
   CALL ReadCom(UnIn,DvrInputFile,'Section Header: Gravity Parameter',ErrStat2,ErrMsg2,UnEc)
   InitInputData%gravity(:) = 0.0_ReKi
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%gravity(1),"InitInputData%gravity(1)", "gravity vector X",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%gravity(2),"InitInputData%gravity(2)", "gravity vector Y",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%gravity(3),"InitInputData%gravity(3)", "gravity vector Z",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;

   !---------------------- FRAME PARAMETER --------------------------------------
   CALL ReadCom(UnIn,DvrInputFile,'Section Header: Frame Parameter',ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   InitInputData%GlbPos(:)   = 0.0_ReKi
   InitInputData%GlbRot(:,:) = 0.0_R8Ki
   InitInputData%RootOri(:,:) = 0.0_R8Ki
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%GlbPos(1),"InitInputData%GlbPos(1)", "position vector X",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%GlbPos(2),"InitInputData%GlbPos(2)", "position vector Y",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%GlbPos(3),"InitInputData%GlbPos(3)", "position vector Z",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;

   CALL ReadCom(UnIn,DvrInputFile,'Comments on DCM',ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadCom(UnIn,DvrInputFile,'Comments on DCM',ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   DO i=1,3
       CALL ReadAry(UnIn,DvrInputFile,InitInputData%RootOri(i,:),3,"InitInputData%RootOri",&
               "Initial root orientation",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   ENDDO

   CALL ReadVar(UnIn,DvrInputFile,DvrData%GlbRotBladeT0,"DvrData%GlbRotBladeT0","Is the blade initial orientation also the GlbRot calculation frame",ErrStat2,ErrMSg2,UnEc);   if (Failed())  return;

      ! Use the initial blade root orientation as the GlbRot reference orientation for all calculations?
   if ( DvrData%GlbRotBladeT0 ) then
         ! Set the GlbRot matrix
      InitInputData%GlbRot = InitInputData%RootOri
      CALL eye( DvrData%RootRelInit, ErrStat2, ErrMsg2 );   if (Failed())  return;
   else
         ! Initialize the GlbRot matrix as the identity.  Relative rotation for root to GlbRot
      DvrData%RootRelInit = InitInputData%RootOri
       CALL eye( InitInputData%GlbRot, ErrStat2, ErrMsg2 );   if (Failed())  return;
   end if

   !---------------------- INITIAL VELOCITY PARAMETER --------------------------------
   CALL ReadCom(UnIn,DvrInputFile,'Section Header: Initial Velocity Parameter',ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%RootVel(4),"InitInputData%IniRootVel(1)", "angular velocity vector X",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%RootVel(5),"InitInputData%IniRootVel(2)", "angular velocity vector Y",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,InitInputData%RootVel(6),"InitInputData%IniRootVel(3)", "angular velocity vector Z",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   InitInputData%RootVel(1:3) = cross_product(InitInputData%RootVel(4:6),InitInputData%GlbPos(:))
  
   !---------------------- APPLIED FORCE --------------------------------
   CALL ReadCom(UnIn,DvrInputFile,'Section Header: Applied Force',ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%DistrLoad(1),"InitInputData%DistrLoad(1)", "Distributed load vector X",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%DistrLoad(2),"InitInputData%DistrLoad(2)", "Distributed load vector Y",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%DistrLoad(3),"InitInputData%DistrLoad(3)", "Distributed load vector Z",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%DistrLoad(4),"InitInputData%DistrLoad(4)", "Distributed load vector X",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%DistrLoad(5),"InitInputData%DistrLoad(5)", "Distributed load vector Y",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%DistrLoad(6),"InitInputData%DistrLoad(6)", "Distributed load vector Z",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%TipLoad(1),"InitInputData%TipLoad(1)", "Tip load vector X",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%TipLoad(2),"InitInputData%TipLoad(2)", "Tip load vector Y",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%TipLoad(3),"InitInputData%TipLoad(3)", "Tip load vector Z",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%TipLoad(4),"InitInputData%TipLoad(4)", "Tip load vector X",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%TipLoad(5),"InitInputData%TipLoad(5)", "Tip load vector Y",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%TipLoad(6),"InitInputData%TipLoad(6)", "Tip load vector Z",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;

      !---------------------- MULTI-POINT LOAD INPUTS ----------------------------------------
   !First read into temporary "line" variable so we can check if this is numeric or not (for backward compatibility)
   CALL ReadVar(UnIn,DvrInputFile,line,"DvrData%NumPointLoads", "Number of Point Loads (primary input file section header)",ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
         
   READ( Line, *, IOSTAT=IOS) DvrData%NumPointLoads
   if (IOS == 0) then !this is numeric, so we can go ahead with the multi-point loads
            
      CALL ReadCom(UnIn,DvrInputFile,'Multiple Point Loads Table',ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
      CALL ReadCom(UnIn,DvrInputFile,'Multiple Point Loads Table Units',ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
      CALL AllocAry(DvrData%MultiPointLoad,max(1,DvrData%NumPointLoads),7,'Point loads input array',ErrStat2,ErrMsg2);   if (Failed())  return;
      DvrData%MultiPointLoad = 0.0_ReKi      ! this must have at least one node, and it will be initialized to 0      
         
      DO i = 1,DvrData%NumPointLoads
          CALL ReadAry( UnIn, DvrInputFile, TmpReAry, 7, 'PointLoad', 'Nodal point loads - Node No., DOF No., ', ErrStat2, ErrMsg2, UnEc );   if (Failed())  return;
          DvrData%MultiPointLoad(i,1) =  TmpReAry(1)
          DvrData%MultiPointLoad(i,2) =  TmpReAry(2)
          DvrData%MultiPointLoad(i,3) =  TmpReAry(3)
          DvrData%MultiPointLoad(i,4) =  TmpReAry(4)
          DvrData%MultiPointLoad(i,5) =  TmpReAry(5)
          DvrData%MultiPointLoad(i,6) =  TmpReAry(6)
          DvrData%MultiPointLoad(i,7) =  TmpReAry(7)
      ENDDO  
      
      DvrData%NumPointLoads = max(1,DvrData%NumPointLoads) 
   
      !---------------------- BEAM SECTIONAL PARAMETER ----------------------------------------
      CALL ReadCom(UnIn,DvrInputFile,'Section Header: Primary input file',ErrStat2,ErrMsg2,UnEc);   if (Failed())  return;
         
   else
      DvrData%NumPointLoads = 1
      CALL AllocAry(DvrData%MultiPointLoad,DvrData%NumPointLoads,7,'Point loads input array',ErrStat2,ErrMsg2);   if (Failed())  return;
      DvrData%MultiPointLoad = 0.0_ReKi           
   end if ! we read the header already
      
      !---------------------- BEAM SECTIONAL PARAMETER ----------------------------------------
   CALL ReadVar ( UnIn, DvrInputFile, InitInputData%InputFile, 'InputFile', 'Name of the primary input file', ErrStat2,ErrMsg2, UnEc );   if (Failed())  return;
      IF ( PathIsRelative( InitInputData%InputFile ) ) InitInputData%InputFile = TRIM(PriPath)//TRIM(InitInputData%InputFile)
      
      !---------------------- Outputs ---------------------------------------------------------
   CALL ReadCom(UnIn,DvrInputFile,'Section Header: Outputs',ErrStat2,ErrMsg2,UnEc); if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%WrVTK,'WrVTK','WrVTK',ErrStat2,ErrMsg2,UnEc); if (Failed())  return;
   CALL ReadVar(UnIn,DvrInputFile,DvrData%VTK_fps,'VTK_fps','VTK_fps',ErrStat2,ErrMsg2,UnEc); if (Failed())  return;

      ! FIXME: added error check here, but probably should be done with more comprehensive error checks on the input file
   if (DvrData%WrVTK < 0 .or. DvrData%WrVTK > 2) then
      ErrStat2 = ErrID_Fatal;    ErrMsg2  = "WrVTK must be 0=none; 1=init; 2=animation";
      if (Failed())  return;
   endif



   call cleanup()
   return
      
contains
   subroutine cleanup() 
      close(UnIn)
      return
   end subroutine cleanup
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'BD_ReadDvrFile')
      Failed = ErrStat>=ErrID_Fatal
      if (Failed) call cleanup()
      return
   end function Failed
END SUBROUTINE BD_ReadDvrFile

SUBROUTINE Dvr_InitializeOutputFile(OutUnit,IntOutput,RootName,ErrStat,ErrMsg)


   INTEGER(IntKi),              INTENT(  OUT):: OutUnit
   TYPE(BD_InitOutputType),     INTENT(IN   ):: IntOutput     ! Output for initialization routine
   INTEGER(IntKi),              INTENT(  OUT):: ErrStat     ! Error status of the operation
   CHARACTER(*),                INTENT(  OUT):: ErrMsg      ! Error message if ErrStat /= ErrID_None
   CHARACTER(*),                INTENT(IN   ):: RootName

   integer(IntKi)                            :: i      
   integer(IntKi)                            :: numOuts
   INTEGER(IntKi)                            :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)                      :: ErrMsg2                      ! Temporary Error message
   character(*), parameter                   :: RoutineName = 'Dvr_InitializeOutputFile'

   ErrStat = ErrID_none
   ErrMsg  = ""
   
   CALL GetNewUnit(OutUnit,ErrStat2,ErrMsg2)
   CALL OpenFOutFile ( OutUnit, trim(RootName)//'.out', ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      if (ErrStat >= AbortErrLev) return

   write (OutUnit,'(/,A)')  'Predictions were generated on '//CurDate()//' at '//CurTime()//' using '//trim(GetNVD(IntOutput%Ver))
   write (OutUnit,'()' )    !print a blank line
   
   numOuts = size(IntOutput%WriteOutputHdr)
   !......................................................
   ! Write the names of the output parameters on one line:
   !......................................................

   write (OutUnit,'()')
   write (OutUnit,'()')
   write (OutUnit,'()')

   call WrFileNR ( OutUnit, 'Time' )

   do i=1,NumOuts
      call WrFileNR ( OutUnit, tab//IntOutput%WriteOutputHdr(i) )
   end do ! i

   write (OutUnit,'()')

      !......................................................
      ! Write the units of the output parameters on one line:
      !......................................................

   call WrFileNR ( OutUnit, '(s)' )

   do i=1,NumOuts
      call WrFileNR ( Outunit, tab//trim(IntOutput%WriteOutputUnt(i)) )
   end do ! i

   write (OutUnit,'()')  


END SUBROUTINE Dvr_InitializeOutputFile

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Dvr_WriteOutputLine(t,OutUnit, OutFmt, Output)


   real(DbKi)             ,  intent(in   )   :: t                    ! simulation time (s)
   INTEGER(IntKi)         ,  intent(in   )   :: OutUnit              ! Status of error message
   CHARACTER(*)           ,  intent(in   )   :: OutFmt
!   real(ReKi)             ,  intent(in   )   :: output(:)            ! Rootname for the output file
   TYPE(BD_OutputType),      INTENT(IN   )   :: Output
      
   ! Local variables.

   integer(IntKi)                            :: errStat              ! Status of error message (we're going to ignore errors in writing to the file)
   character(ErrMsgLen)                      :: errMsg               ! Error message if ErrStat /= ErrID_None
   character(200)                            :: frmt                 ! A string to hold a format specifier
   character(15)                             :: tmpStr               ! temporary string to print the time output as text

   frmt = '"'//tab//'"'//trim(OutFmt)      ! format for array elements from individual modules
   
      ! time
   write( tmpStr, '(F15.6)' ) t
   call WrFileNR( OutUnit, tmpStr )
   call WrNumAryFileNR ( OutUnit, Output%WriteOutput,  frmt, errStat, errMsg )
   
     ! write a new line (advance to the next line)
   write (OutUnit,'()')
      
end subroutine Dvr_WriteOutputLine
  
!----------------------------------------------------------------------------------------------------------------------------------
subroutine CreateMultiPointMeshes(DvrData,BD_InitInput,BD_InitOutput,BD_Parameter,y, u, ErrStat,ErrMsg)

   TYPE(BD_DriverInternalType), INTENT(INOUT) :: DvrData
   TYPE(BD_InitInputType)     , INTENT(IN   ) :: BD_InitInput
   TYPE(BD_InitOutputType)    , INTENT(IN   ) :: BD_InitOutput
   TYPE(BD_ParameterType)     , INTENT(IN   ) :: BD_Parameter
   TYPE(BD_OutputType),         INTENT(IN   ) :: y
   TYPE(BD_InputType),          INTENT(INOUT) :: u                ! sets pointLoad with values from BD driver input file
   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat          ! Error status of the operation
   CHARACTER(*),                INTENT(  OUT) :: ErrMsg           ! Error message if ErrStat /= ErrID_None

   integer(intKi)                             :: i
   integer(intKi)                             :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                       :: ErrMsg2           ! temporary Error message
   character(*), parameter                    :: RoutineName = 'CreateMultiPointMeshes'
   
   REAL(BDKi)                                 :: temp_POS(3)
   REAL(BDKi)                                 :: temp_CRV(3)
   REAL(BDKi)                                 :: temp_CRV2(3)
   REAL(R8Ki)                                 :: DCM(3,3)          ! must be same type as mesh orientation fields
   REAL(ReKi)                                 :: Pos(3)            ! must be same type as mesh position fields
   REAL(BDKi)                                 :: TmpDCM(3,3)
   
   
   ErrStat = ErrID_None
   ErrMsg = ""
   
   ! DvrData%NumPointLoads is at least 1
   
   !.......................
   ! Mesh for multi-point loading on blades
   !.......................
   CALL MeshCreate( BlankMesh        = DvrData%mplMotion  &
                   ,IOS              = COMPONENT_INPUT    &
                   ,NNodes           = DvrData%NumPointLoads      &
                   ,TranslationDisp  = .TRUE.             &
                   ,Orientation      = .TRUE.             &
                   ,TranslationVel   = .TRUE.             &
                   ,RotationVel      = .TRUE.             &
                   ,TranslationAcc   = .TRUE.             &
                   ,RotationAcc      = .TRUE.             &
                   ,ErrStat          = ErrStat2           &
                   ,ErrMess          = ErrMsg2             )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) return   
      
      ! these nodes are placed along the key point line (as are the GLL nodes)
   DO i = 1,DvrData%NumPointLoads

       !FIXME Mike Sprague; given DvrData%MultiPointLoad(i,1), find temp_POS and temp_CRV; we should have p%uuN0
       !    Find_IniNode(kp_coordinate,                          p, member_first_kp, member_last_kp,                                 eta, POS,           CRV, ErrStat, ErrMsg)
       !call Find_IniNode(BD_InitOutput%kp_coordinate, BD_Parameter,               1, BD_InitOutput%kp_total, DvrData%MultiPointLoad(i,1), temp_POS, temp_CRV, ErrStat2, ErrMsg2)
       call BD_Interp_Pos_CRV(BD_Parameter, DvrData%MultiPointLoad(i,1), temp_POS, temp_CRV, ErrStat, ErrMsg)
       CALL SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
       if (ErrStat >= AbortErrLev) return
       
       Pos = BD_Parameter%GlbPos + MATMUL(BD_Parameter%GlbRot,temp_POS)
       
       temp_CRV2 = MATMUL(BD_Parameter%GlbRot,temp_CRV)
       CALL BD_CrvCompose(temp_CRV,BD_Parameter%Glb_crv,temp_CRV2,FLAG_R1R2) !temp_CRV = p%Glb_crv composed with temp_CRV2

       CALL BD_CrvMatrixR(temp_CRV,TmpDCM) ! returns TmpDCM (the transpose of the DCM orientation matrix)

       ! possible type conversions here:
       DCM = TRANSPOSE(TmpDCM)

       ! set the reference position and orientation for each node.
       CALL MeshPositionNode ( Mesh    = DvrData%mplMotion     &
                              ,INode   = i             &
                              ,Pos     = Pos           &
                              ,ErrStat = ErrStat2      &
                              ,ErrMess = ErrMsg2       &
                              ,Orient  = DCM           )
       CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ENDDO
   
   DO i = 1,DvrData%NumPointLoads
       CALL MeshConstructElement( Mesh     = DvrData%mplMotion      &
                                 ,Xelement = ELEMENT_POINT    &
                                 ,P1       = i                &
                                 ,ErrStat  = ErrStat2         &
                                 ,ErrMess  = ErrMsg2          )
         CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ENDDO
   
   CALL MeshCommit ( Mesh    = DvrData%mplMotion       &
                    ,ErrStat = ErrStat2        &
                    ,ErrMess = ErrMsg2         )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   !.......................
   ! sibling mesh for motions, which are needed to transfer loads
   !.......................

   CALL MeshCopy ( SrcMesh  = DvrData%mplMotion    &
                 , DestMesh = DvrData%mplLoads     &
                 , CtrlCode = MESH_SIBLING         &
                 , IOS      = COMPONENT_INPUT      &
                 , Force    = .TRUE.               &
                 , Moment   = .TRUE.               &
                 , ErrStat  = ErrStat2             &
                 , ErrMess  = ErrMsg2              )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

      
   !.......................
   ! "sibling" mesh for BD motions, which are needed to transfer loads
   ! BD's input and output meshes may not be at the same location, so we need another copy
   !.......................
      
   CALL MeshCopy ( SrcMesh  = u%PointLoad                     &
                 , DestMesh = DvrData%y_BldMotion_at_u_point       &
                 , CtrlCode = MESH_COUSIN                  &  ! Like a sibling, except using new memory for position/refOrientation and elements
                 , IOS      = COMPONENT_OUTPUT             &
                 , Orientation     = .TRUE.                &
                 , TranslationDisp = .TRUE.                &
                 , ErrStat         = ErrStat2              &
                 , ErrMess         = ErrMsg2               )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   if (ErrStat >= AbortErrLev) return      
      
   !.......................
   ! initialize the mapping between the BD output motions and driver mpl motion mesh:
   !.......................
      
   CALL MeshMapCreate( y%BldMotion, DvrData%mplMotion, DvrData%Map_BldMotion_to_mplMotion, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   !.......................
   ! initialize the mapping between the driver mpl loads and BD input point loads mesh:
   !.......................
   CALL MeshMapCreate( DvrData%mplLoads, u%PointLoad, DvrData%Map_mplLoads_to_PointLoad, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   !.......................
   ! initialize the mapping between the BD output motions and BD motions at the nodes on the y%BldMotion mesh:
   !.......................
      
   CALL MeshMapCreate( y%BldMotion, DvrData%y_BldMotion_at_u_point, DvrData%Map_y_BldMotion_to_u_point, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
      
   DvrData%mplMotion%remapFlag = .false.
   DvrData%mplLoads%remapFlag = .false.
   DvrData%y_BldMotion_at_u_point%remapFlag = .false.
      
end subroutine CreateMultiPointMeshes
!----------------------------------------------------------------------------------------------------------------------------------
subroutine Transfer_MultipointLoads(DvrData, y, u, ErrStat, ErrMsg)

   TYPE(BD_DriverInternalType), INTENT(INOUT) :: DvrData
   TYPE(BD_OutputType),         INTENT(IN   ) :: y
   TYPE(BD_InputType),          INTENT(INOUT) :: u                ! sets pointLoad with values from BD driver input file
   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat          ! Error status of the operation
   CHARACTER(*),                INTENT(  OUT) :: ErrMsg           ! Error message if ErrStat /= ErrID_None

   integer(intKi)                             :: i
   integer(intKi)                             :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                       :: ErrMsg2           ! temporary Error message
   character(*), parameter                    :: RoutineName = 'Transfer_MultipointLoads'

   
   ErrStat = ErrID_None
   ErrMsg = ""
   
      ! the mpl loads are constant over the course of the simulation, but they could be set up to vary
   DO i = 1,DvrData%NumPointLoads
      DvrData%mplLoads%Force(1:3,i)  = DvrData%MultiPointLoad(i,2:4)
      DvrData%mplLoads%Moment(1:3,i) = DvrData%MultiPointLoad(i,5:7)
   ENDDO
      
      ! get the motions from BD to the two meshes needed in this loads transfer:
   CALL Transfer_Line2_to_Point( y%BldMotion, DvrData%y_BldMotion_at_u_point, DvrData%Map_y_BldMotion_to_u_point, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   CALL Transfer_Line2_to_Point( y%BldMotion, DvrData%mplMotion, DvrData%Map_BldMotion_to_mplMotion, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   ! we'll do this each time step in case we want to make this time-series input at some point
      ! transfer the mpl loads to BD's input point mesh at the GLL nodes:
   CALL Transfer_Point_to_Point( DvrData%mplLoads, u%PointLoad, DvrData%Map_mplLoads_to_PointLoad, ErrStat2, ErrMsg2, DvrData%mplMotion, DvrData%y_BldMotion_at_u_point)  
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

end subroutine Transfer_MultipointLoads

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Init_RotationCenterMesh(DvrData, InitInputData, RootMotionMesh, ErrStat, ErrMsg)
   TYPE(BD_InitInputType),      INTENT(IN   ) :: InitInputData
   TYPE(BD_DriverInternalType), INTENT(INOUT) :: DvrData
   TYPE(MeshType),              INTENT(INOUT) :: RootMotionMesh
   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat           ! Error status of the operation
   CHARACTER(*),                INTENT(  OUT) :: ErrMsg            ! Error message if ErrStat /= ErrID_None

   integer(intKi)                             :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                       :: ErrMsg2           ! temporary Error message
   character(*), parameter                    :: RoutineName = 'Init_RotationCenterMesh'
   
   real(r8Ki)                                 :: orientation(3,3)
   real(ReKi)                                 :: position(3)
   real(R8Ki)                                 :: z_hat(3)  ! unit-magnitude rotational velocity vector
   real(R8Ki)                                 :: Z_unit(3) ! unit vector in the Z direction
   real(R8Ki)                                 :: vec(3)    ! temporary vector
   
   ErrStat = ErrID_None
   ErrMsg = ''
   
   
   position = 0.0_ReKi  ! center of rotation

   DvrData%w = TwoNorm( InitInputData%RootVel(4:6) )
   
   if (EqualRealNos(DvrData%w,0.0_R8Ki)) then
      DvrData%w = 0.0_R8Ki
         ! the beam is not rotating, so pick an orientation
      call eye(orientation, ErrStat2, ErrMsg2)
   else
      z_hat = InitInputData%RootVel(4:6) / DvrData%w
      
      if ( EqualRealNos( z_hat(3), 1.0_R8Ki ) ) then
         call eye(orientation, ErrStat2, ErrMsg2)
      elseif ( EqualRealNos( z_hat(3), -1.0_R8Ki ) ) then
         orientation = 0.0_ReKi
         orientation(1,1) = -1.0_R8Ki
         orientation(2,2) =  1.0_R8Ki
         orientation(3,3) = -1.0_R8Ki
      else
         
         Z_unit = (/0.0_R8Ki, 0.0_R8Ki, 1.0_R8Ki/)
         
         vec = Z_unit - z_hat*z_hat(3) ! vec = matmul( eye(3) - outerproduct(z_hat,z_hat), (/ 0,0,1/) )
         vec = vec / TwoNorm(vec)      ! we've already checked that this is not zero
         orientation(1,:) = vec

         vec = cross_product(z_hat,Z_unit)
         vec = vec / TwoNorm(vec)
         orientation(2,:) = vec
                  
         orientation(3,:) = z_hat
         
      end if
      
   end if
   
   !.......................
   ! Mesh for center of rotation
   !.......................
   CALL MeshCreate( BlankMesh        = DvrData%RotationCenter &
                   ,IOS              = COMPONENT_OUTPUT       &
                   ,NNodes           = 1                      &
                   ,TranslationDisp  = .TRUE.                 &
                   ,Orientation      = .TRUE.                 &
                   ,TranslationVel   = .TRUE.                 &
                   ,RotationVel      = .TRUE.                 &
                   ,TranslationAcc   = .TRUE.                 &
                   ,RotationAcc      = .TRUE.                 &
                   ,ErrStat          = ErrStat2               &
                   ,ErrMess          = ErrMsg2                 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) return   
      
       ! set the reference position and orientation.
   CALL MeshPositionNode ( DvrData%RotationCenter, 1, position, ErrStat2, ErrMsg2, orientation )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL MeshConstructElement( DvrData%RotationCenter, ELEMENT_POINT, ErrStat2, ErrMsg2, 1 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

   CALL MeshCommit (DvrData%RotationCenter, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      IF (ErrStat >= AbortErrLev) return   
   
      
   ! note that the following fields do not change during the simulation;  orientation will change in BD_InputSolve.
   DvrData%RotationCenter%TranslationDisp  = 0.0_ReKi
   DvrData%RotationCenter%TranslationVel   = 0.0_ReKi
   DvrData%RotationCenter%RotationVel(:,1) = InitInputData%RootVel(4:6)
   DvrData%RotationCenter%TranslationAcc   = 0.0_ReKi
   DvrData%RotationCenter%RotationAcc      = 0.0_ReKi
   
   
   !.......................
   ! initialize the mapping between the BD center of rotation and BD root motion input mesh:
   !.......................
      
   CALL MeshMapCreate( DvrData%RotationCenter, RootMotionMesh, DvrData%Map_RotationCenter_to_RootMotion, ErrStat2, ErrMsg2 )
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
      
   DvrData%RotationCenter%remapFlag = .false.
   RootMotionMesh%remapFlag = .false.
   
END SUBROUTINE Init_RotationCenterMesh

!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD_InputSolve( t, u, DvrData, ErrStat, ErrMsg)
 
   REAL(DbKi),                  INTENT(IN   ) :: t
   TYPE(BD_InputType),          INTENT(INOUT) :: u
   TYPE(BD_DriverInternalType), INTENT(INOUT) :: DvrData
   INTEGER(IntKi),              INTENT(  OUT) :: ErrStat          ! Error status of the operation
   CHARACTER(*),                INTENT(  OUT) :: ErrMsg           ! Error message if ErrStat /= ErrID_None
                                         
   ! local variables                     
   INTEGER(IntKi)                             :: i                ! do-loop counter
   REAL(R8Ki)                                 :: Orientation(3,3)
   REAL(R8Ki)                                 :: wt               ! time from start start of simulation multiplied by magnitude of rotational velocity
   REAL(R8Ki)                                 :: swt, cwt         ! sine and cosine of w*t
   
   integer(intKi)                             :: ErrStat2         ! temporary Error status
   character(ErrMsgLen)                       :: ErrMsg2          ! temporary Error message
   character(*), parameter                    :: RoutineName = 'BD_InputSolve'
   
   ErrStat = ErrID_None
   ErrMsg  = ''

   !.............................
   ! Set up u%RootMotion: 
   !.............................
   
   ! Compute the orientation at the center of rotation 
   wt = DvrData%w * (t - DvrData%t_initial)
   swt = sin( wt )
   cwt = cos( wt )
   Orientation(1,1) = cwt
   Orientation(2,1) =-swt
   Orientation(3,1) =   0.0_R8Ki
   
   Orientation(1,2) = swt
   Orientation(2,2) = cwt
   Orientation(3,2) = 0.0_R8Ki
   
   Orientation(1,3) = 0.0_R8Ki
   Orientation(2,3) = 0.0_R8Ki
   Orientation(3,3) = 1.0_R8Ki
      
   DvrData%RotationCenter%Orientation(:,:,1) = matmul(Orientation, matmul(DvrData%RotationCenter%RefOrientation(:,:,1),DvrData%RootRelInit))

   CALL Transfer_Point_to_Point( DvrData%RotationCenter, u%RootMotion, DvrData%Map_RotationCenter_to_RootMotion, ErrStat2, ErrMsg2)  
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )         
      
   !.............................
   ! set up the point load input:
   ! @VA: if we want to apply these at different positions, we should call Transfer_MultipointLoads(); 
   ! putting the calculation here so we can have changing loads at some point...
   !.............................
   CALL Transfer_Point_to_Point( DvrData%mplLoads, u%PointLoad, DvrData%Map_mplLoads_to_PointLoad, ErrStat2, ErrMsg2, DvrData%mplMotion, DvrData%y_BldMotion_at_u_point)  
      CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
   
   u%PointLoad%Force(1:3,u%PointLoad%NNodes)  = u%PointLoad%Force(1:3,u%PointLoad%NNodes)  + DvrData%TipLoad(1:3)
   u%PointLoad%Moment(1:3,u%PointLoad%NNodes) = u%PointLoad%Moment(1:3,u%PointLoad%NNodes) + DvrData%TipLoad(4:6)
   
   !.............................
   ! LINE2 mesh: DistrLoad
   !.............................
   DO i=1,u%DistrLoad%NNodes
      u%DistrLoad%Force(:,i) =  DvrData%DistrLoad(1:3)
      u%DistrLoad%Moment(:,i)=  DvrData%DistrLoad(4:6)
   ENDDO

END SUBROUTINE BD_InputSolve


end module BeamDyn_driver_subs
