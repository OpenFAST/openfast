module FVW_IO

  USE FVW_Types
  USE FVW_Subs
  use FVW_VortexTools
  implicit none

contains

! ==============================================================================
!> Reads the input file for FVW
SUBROUTINE FVW_ReadInputFile( FileName, p, Inp, ErrStat, ErrMsg )
   character(len=*),             intent(in)    :: FileName !< Input file name for FVW
   type( FVW_ParameterType ),    intent(inout) :: p        !< Parameters
   type(FVW_InputFile),          intent(out)   :: Inp      !< Data stored in the module's input file
   integer(IntKi),               intent(  out) :: ErrStat  !< Error status of the operation
   character(*),                 intent(  out) :: ErrMsg   !< Error message if ErrStat /= ErrID_None
   ! Local variables
   character(1024)      :: PriPath                         ! the path to the primary input file
   character(1024)      :: VTK_fps_line                    ! string to temporarially hold value of read line for VTK_fps
   integer(IntKi)       :: UnIn
   integer(IntKi)       :: ErrStat2
   character(ErrMsgLen) :: ErrMsg2
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! Open file
   CALL GetNewUnit( UnIn )   
   CALL OpenFInpfile(UnIn, TRIM(FileName), ErrStat2, ErrMsg2)
   if (Check( ErrStat2 /= ErrID_None , 'Could not open input file')) return
   CALL GetPath( FileName, PriPath )    ! Input files will be relative to the path where the primary input file is located.
   !------------------------------------- HEADER ---------------------------------------------------
   CALL ReadCom(UnIn, FileName, 'FVW input file header line 1', ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL ReadCom(UnIn, FileName, 'FVW input file header line 2', ErrStat2, ErrMsg2 ); if(Failed()) return
   !------------------------ GENERAL OPTIONS  -------------------------------------------
   CALL ReadCom        (UnIn,FileName,                         'General option header'                                  , ErrStat2,ErrMsg2); if(Failed()) return
   CALL ReadVarWDefault(UnIn,FileName,Inp%IntMethod           ,'Integration method' ,'', idEuler1                       , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%DTfvw               ,'DTfvw'              ,'',  p%DTaero                      , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%FreeWakeStart       ,'FreeWakeStart'       ,'', 0.0_ReKi                      , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%FullCirculationStart,'FullCirculationStart','', real(20.0_ReKi*Inp%DTfvw,ReKi), ErrStat2,ErrMsg2); if(Failed())return
   !------------------------ CIRCULATION SPECIFICATIONS  -------------------------------------------
   CALL ReadCom(UnIn,FileName,                  'Circulation specification header', ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL ReadVarWDefault(UnIn,FileName,Inp%CirculationMethod ,'CirculationMethod' ,'', idCircPolarData, ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%CircSolvConvCrit  ,'CircSolvConvCrit ' ,'', 0.001          , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%CircSolvRelaxation,'CircSolvRelaxation','', 0.1            , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%CircSolvMaxIter   ,'CircSolvMaxIter'   ,'', 30             , ErrStat2,ErrMsg2); if(Failed())return
   !CALL ReadVar(UnIn,FileName,Inp%CircSolvPolar     ,'CircSolvPolar'   ,'',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%CirculationFile   ,'CirculationFile'   ,'',ErrStat2,ErrMsg2); if(Failed())return
   !------------------------ WAKE OPTIONS -------------------------------------------
   CALL ReadCom        (UnIn,FileName,                  'Wake options header'                         , ErrStat2,ErrMsg2); if(Failed()) return
   CALL ReadVar        (UnIn,FileName,Inp%nNWPanels          ,'nNWPanels'         ,''                 , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar        (UnIn,FileName,Inp%nFWPanels          ,'nFWPanels'         ,''                 , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%nFWPanelsFree      ,'nFWPanelsFree'     ,'', Inp%nFWPanels  , ErrStat2,ErrMsg2); if(Failed())return

   CALL ReadVarWDefault(UnIn,FileName,Inp%FWShedVorticity    ,'FWShedVorticity'   ,'', .False.        , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%DiffusionMethod    ,'DiffusionMethod'   ,'',idDiffusionNone , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%RegDeterMethod     ,'RegDeterMethod'    ,'',idRegDeterManual, ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%RegFunction        ,'RegFunction'       ,'',idRegVatistas   , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%WakeRegMethod      ,'WakeRegMethod'     ,'',idRegConstant   , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar        (UnIn,FileName,Inp%WakeRegParam       ,'WakeRegParam'      ,''                 , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar        (UnIn,FileName,Inp%WingRegParam       ,'WingRegParam'      ,''                 , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%CoreSpreadEddyVisc ,'CoreSpreadEddyVisc','',100.0_ReKi      , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%ShearModel         ,'ShearModel'        ,'',idShearNone     , ErrStat2,ErrMsg2); if(Failed())return
   !CALL ReadVarWDefault(UnIn,FileName,Inp%TwrShadowOnWake    ,'TwrShadowOnWake'   ,'',.false.     , ErrStat2,ErrMsg2); if(Failed())return
   !CALL ReadVarWDefault(UnIn,FileName,Inp%TreeModel          ,'TreeModel'         ,'',idTreeNone  , ErrStat2,ErrMsg2); if(Failed())return
   !CALL ReadVarWDefault(UnIn,FileName,Inp%TreeBranchFactor   ,'TreeBranchFactor'  ,'',3.0_ReKi    , ErrStat2,ErrMsg2); if(Failed())return
   !CALL ReadVarWDefault(UnIn,FileName,Inp%TreeBranchSmall    ,'TreeBranchSmall'   ,'',0.1_ReKi    , ErrStat2,ErrMsg2); if(Failed())return
   Inp%TwrShadowOnWake  = .False.
   Inp%TreeModel        = idTreeNone
   Inp%TreeBranchFactor = 3.0_ReKi
   Inp%TreeBranchSmall  = 0.1_ReKi
   !------------------------ OUTPUT OPTIONS -----------------------------------------
   CALL ReadCom        (UnIn,FileName,                  'Output options header'              ,ErrStat2,ErrMsg2); if(Failed()) return
   CALL ReadVarWDefault(UnIn,FileName,Inp%WrVTK       , 'WrVTK'              ,'',     0      ,ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%VTKBlades   , 'VTKBlades'          ,'',     1      ,ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%VTKCoord    , 'VTKCoord'           ,'',     1      ,ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar        (UnIn,FileName,VTK_fps_line    , 'VTK_fps'            ,''             ,ErrStat2,ErrMsg2); if(Failed())return

   ! --- Validation of inputs
   if (PathIsRelative(Inp%CirculationFile)) Inp%CirculationFile = TRIM(PriPath)//TRIM(Inp%CirculationFile)

   if (Check(.not.(ANY(idCircVALID ==Inp%CirculationMethod)), 'Circulation method (CircSolvingMethod) not implemented')) return
   if (Check(.not.(ANY(idIntMethodVALID==Inp%IntMethod    )) , 'Time integration method (IntMethod) not yet implemented. Use Euler 1st order method for now.')) return
   if (Check(.not.(ANY(idDiffusionVALID==Inp%DiffusionMethod)) , 'Diffusion method (DiffusionMethod) not yet implemented. Use None for now.')) return
   if (Check(.not.(ANY(idRegDeterVALID ==Inp%RegDeterMethod))  , 'Regularization determination method (RegDeterMethod) not yet implemented. Use Manual method for now.')) return
   if (Check(.not.(ANY(idRegVALID      ==Inp%RegFunction  )), 'Regularization function (RegFunction) not implemented')) return
   if (Check(.not.(ANY(idRegMethodVALID==Inp%WakeRegMethod)), 'Wake regularization method (WakeRegMethod) not implemented')) return
   if (Check(.not.(ANY(idShearVALID    ==Inp%ShearModel   )), 'Shear model (`ShearModel`) not valid')) return
   if (Check(.not.(ANY(idTreeVALID     ==Inp%TreeModel    )), 'Shear model (`ShearModel`) not valid')) return

   if (Check( Inp%DTfvw < p%DTaero, 'DTfvw must be >= DTaero from AD15.')) return
   if (abs(Inp%DTfvw-p%DTaero)>epsilon(1.0_ReKi)) then
      ! subcycling
      if (Check(Inp%IntMethod/=idEuler1 , 'Sub-cycling (DTfvw>DTaro) is only possible with Forward Euler `IntMethod`')) return
   endif
   if (Inp%CirculationMethod == idCircPolarData) then
      if (Check( Inp%nNWPanels<1 , 'Number of near wake panels (`nNWPanels`) must be >=1 when using circulation solving with polar data (`CircSolvingMethod=1`)')) return
   endif

   if (Check( Inp%nNWPanels<0     , 'Number of near wake panels must be >=0')) return
   if (Check( Inp%nFWPanels<0     , 'Number of far wake panels must be >=0')) return
   if (Check( Inp%nFWPanelsFree<0 , 'Number of free far wake panels must be >=0')) return
   if (Check( Inp%nFWPanelsFree>Inp%nFWPanels , 'Number of free far wake panels must be <=Number of far wake panels')) return

   if (Check(Inp%WakeRegParam<0             , 'Wake regularization parameter (WakeRegParam) should be positive')) return
   if (Check(Inp%WingRegParam<0             , 'Wing regularization parameter (WakeRegParam) should be positive')) return
   if (Check(Inp%CoreSpreadEddyVisc<0       , 'Core spreading eddy viscosity (CoreSpreadEddyVisc) should be positive')) return

   ! Removing the shed vorticity is a dangerous option if this is done too close to the blades. 
   ! To be safe, we will no matter what ensure that the last segments of NW are 0 if FWShedVorticity is False (see PackPanelsToSegments)
   ! Still we force the user to be responsible.
   if (Check((.not.(Inp%FWShedVorticity)) .and. Inp%nNWPanels<30, '`FWShedVorticity` should be true if `nNWPanels`<30. Alternatively, use a larger number of NWPanels  ')) return

   Inp%DTvtk = Get_DTvtk( VTK_fps_line, p%DTaero, Inp%DTfvw )

   ! At least one NW panel if FW, this shoudln't be a problem since the LL is in NW, but safety for now
   !if (Check( (Inp%nNWPanels<=0).and.(Inp%nFWPanels>0)      , 'At least one near wake panel is required if the number of far wake panel is >0')) return
   call CleanUp()

CONTAINS
   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'FVW_ReadInputFile') 
      Failed =  ErrStat >= AbortErrLev
      if (Failed) call CleanUp()
   end function Failed

   logical function Check(Condition, ErrMsg_in)
      logical, intent(in) :: Condition
      character(len=*), intent(in) :: ErrMsg_in
      Check=Condition
      if (Check) then
         call SetErrStat(ErrID_Fatal, 'Error in file '//TRIM(FileName)//': '//trim(ErrMsg_in), ErrStat, ErrMsg, 'FVW_ReadInputFile');
         call CleanUp()
      endif
   end function Check

   subroutine CleanUp()
      close( UnIn )
   end subroutine

   real(DbKi) function Get_DTvtk( VTK_fps_line, DTaero, DTfvw )
      character(len=*), intent(inout)  :: VTK_fps_line
      real(DbKi),       intent(in   )  :: DTaero
      real(DbKi),       intent(in   )  :: DTfvw
      real(DbKi)                       :: VTK_fps
      integer(IntKi)                   :: IOS
      integer(IntKi)                   :: TmpRate
      real(DbKi)                       :: TmpTime

      call Conv2UC( VTK_fps_line )
      if ( index(VTK_fps_line, "DEFAULT" ) == 1 ) then   ! at DTfvw frequency
         Get_DTvtk = DTfvw
      elseif ( index(VTK_fps_line, "ALL" ) == 1 ) then   ! at DTaero frequency
         Get_DTvtk = DTaero
      else  ! read a number.  Calculate this later. {will use closest integer multiple of DT}
         read( VTK_fps_line, *, IOSTAT=IOS) VTK_fps
            CALL CheckIOS ( IOS, FileName, 'VTK_fps', NumType, ErrStat2, ErrMsg2 ); if (Failed()) return;

         ! convert frames-per-second to seconds per sample:
         if ( EqualRealNos(VTK_fps, 0.0_DbKi) ) then
            Get_DTvtk = HUGE(1.0_DbKi)
         else
            TmpTime = 1.0_DbKi / VTK_fps
            TmpRate = max( NINT( TmpTime / DTaero ),1_IntKi )      ! Can't be smaller that DTaero
            Get_DTvtk = TmpRate * DTaero
            ! warn if DTvtk is not TmpTime
            if (.not. EqualRealNos(Get_DTvtk, TmpTime)) then
               call SetErrStat(ErrID_Info, '1/VTK_fps is not an integer multiple of DT. FVW will output VTK information at '//&
                              trim(num2lstr(1.0_DbKi/(TmpRate*DTaero)))//' fps, the closest rate possible.',ErrStat,ErrMsg,'FVW_ReadInputFile')
            end if
         end if
      end if
   end function Get_DTvtk


END SUBROUTINE FVW_ReadInputFile

!=================================================
!> Export FVW variables to VTK
!! NOTE: when entering this function nNW and nFW has been ncremented by 1
subroutine WrVTK_FVW(p, x, z, m, FileRootName, VTKcount, Twidth)
   use VTK ! for all the vtk_* functions
   type(FVW_ParameterType),        intent(in   ) :: p !< Parameters
   type(FVW_ContinuousStateType),  intent(in   ) :: x !< States
   type(FVW_ConstraintStateType),  intent(in   ) :: z !< Constraints
   type(FVW_MiscVarType),          intent(in   ) :: m !< MiscVars
   character(*),    intent(in)           :: FileRootName    !< Name of the file to write the output in (excluding extension)
   integer(IntKi),  intent(in)           :: VTKcount        !< Indicates number for VTK output file (when 0, the routine will also write reference information)
   integer(IntKi),  intent(in)           :: Twidth          !< Number of digits in the maximum write-out step (used to pad the VTK write-out in the filename with zeros)
   ! local variables
   integer:: iW
   character(1024)                       :: FileName
   character(255)                        :: Label
   character(Twidth)                     :: Tstr          ! string for current VTK write-out step (padded with zeros)
   character(1), dimension(3) :: I2ABC =(/'A','B','C'/)
   integer(IntKi)       :: nSeg, nSegP, nSegNW
   logical              :: bMirror
   integer(IntKi)       :: ErrStat2
   character(ErrMsgLen) :: ErrMsg2
   real(Reki), dimension(:,:,:), allocatable :: dxdt_0 !<

   if (DEV_VERSION) then
      print*,'------------------------------------------------------------------------------'
      print'(A,L1,A,I0,A,I0,A,I0)','VTK Output  -      First call ',m%FirstCall, '                                nNW:',m%nNW,' nFW:',m%nFW,'  i:',VTKCount
   endif
   !
   call set_vtk_binary_format(.false.) ! TODO binary fails

   ! TimeStamp
   write(Tstr, '(i' // trim(Num2LStr(Twidth)) //'.'// trim(Num2LStr(Twidth)) // ')') VTKcount

   ! --------------------------------------------------------------------------------}
   ! --- Blade 
   ! --------------------------------------------------------------------------------{
   ! --- Blade Quarter chord points (AC)
   do iW=1,p%VTKBlades
      write(Label,'(A,A)') 'BldPointCP.Bld', i2ABC(iW)
      Filename = TRIM(FileRootName)//'.'//trim(Label)//'.'//Tstr//'.vtk'
      if ( vtk_new_ascii_file(trim(filename),Label) ) then
         call vtk_dataset_polydata(m%CP_LL(1:3,1:p%nSpan,iW))
         call vtk_point_data_init()
         call vtk_point_data_scalar(m%Gamma_ll(    1:p%nSpan,iW),'Gamma_ll')
         call vtk_point_data_vector(m%Vind_ll (1:3,1:p%nSpan,iW),'Vind_ll')
         call vtk_point_data_vector(m%Vtot_ll (1:3,1:p%nSpan,iW),'Vtot_ll')
         call vtk_point_data_vector(m%Vstr_ll (1:3,1:p%nSpan,iW),'Vstr_ll')
         call vtk_point_data_vector(m%Vwnd_ll (1:3,1:p%nSpan,iW),'Vwnd_ll')
         call vtk_point_data_vector(m%Tang    (1:3,1:p%nSpan,iW),'Tangent')
         call vtk_point_data_vector(m%Norm    (1:3,1:p%nSpan,iW),'Normal')
         call vtk_point_data_vector(m%Orth    (1:3,1:p%nSpan,iW),'Orth')
         call vtk_close_file()
      endif
   enddo
   ! --- Lifting line panels
   do iW=1,p%VTKBlades
      write(Label,'(A,A)') 'LL.Bld', i2ABC(iW)
      Filename = TRIM(FileRootName)//'.'//trim(Label)//'.'//Tstr//'.vtk'
      call WrVTK_Lattice(FileName, m%r_LL(1:3,:,:,iW), m%Gamma_LL(:,iW:iW))
   enddo
   ! --------------------------------------------------------------------------------}
   ! --- Near wake 
   ! --------------------------------------------------------------------------------{
   ! --- Near wake panels
   do iW=1,p%VTKBlades
      write(Label,'(A,A)') 'NW.Bld', i2ABC(iW)
      Filename = TRIM(FileRootName)//'.'//trim(Label)//'.'//Tstr//'.vtk'
      if (m%FirstCall) then ! Small Hack - At t=0, NW not set, but first NW panel is the LL panel
         allocate(dxdt_0(3, size(m%dxdt_NW,2) , m%nNW+1)); dxdt_0=0.0_ReKi
         call WrVTK_Lattice(FileName, m%r_LL(1:3,:,1:2,iW), m%Gamma_LL(:,iW:iW),dxdt_0)
         deallocate(dxdt_0)
      else
         call WrVTK_Lattice(FileName, x%r_NW(1:3,:,1:m%nNW+1,iW), x%Gamma_NW(:,1:m%nNW,iW), m%dxdt_NW(:,:,1:m%nNW+1,iW))
      endif
   enddo
   ! --------------------------------------------------------------------------------}
   ! --- Far wake 
   ! --------------------------------------------------------------------------------{
   ! --- Far wake panels
   do iW=1,p%VTKBlades
      write(Label,'(A,A)') 'FW.Bld', i2ABC(iW)
      Filename = TRIM(FileRootName)//'.'//trim(Label)//'.'//Tstr//'.vtk'
      call WrVTK_Lattice(FileName, x%r_FW(1:3,1:FWnSpan+1,1:m%nFW+1,iW), x%Gamma_FW(1:FWnSpan,1:m%nFW,iW),m%dxdt_FW(:,:,1:m%nFW+1,iW))
   enddo
   ! --------------------------------------------------------------------------------}
   ! --- All Segments
   ! --------------------------------------------------------------------------------{
   ! NOTE: now we rely on the fact that the segments in Misc are well set
   !       These segments are correct after a call to CalcOutput
   !       The alternative is to call PackPanelsToSegments as was done before
   !       This would require to allocate some local SegPoints,SegConnct here.
   ! False below is to avoid writing the mirrored vorticity, this could be an option though
   bMirror= (p%ShearModel==idShearMirror) .and. (p%VTKBlades<0) ! NOTE: temporary hack to output mirrored vorticity
   call CountSegments(p, m%nNW, m%nFW, 1, nSeg, nSegP, nSegNW)
   if (bMirror) then
      nSeg  = 2*nSeg
      nSegP = 2*nSegP
   endif
   Filename = TRIM(FileRootName)//'.AllSeg.'//Tstr//'.vtk'
   CALL WrVTK_Segments(Filename, m%SegPoints(:,1:nSegP), m%SegConnct(:,1:nSeg), m%SegGamma(1:nSeg), m%SegEpsilon(1:nSeg)) 

   if(.false.) print*,z%Gamma_LL(1,1) ! unused var for now
end subroutine WrVTK_FVW


subroutine WrVTK_Segments(filename, SegPoints, SegConnct, SegGamma, SegEpsilon) 
   use VTK
   character(len=*),intent(in)                 :: filename
   real(ReKi), dimension(:,:),      intent(in) :: SegPoints  !< 
   integer(IntKi), dimension(:,:),  intent(in) :: SegConnct  !< 
   real(ReKi),     dimension(:)  ,  intent(in) :: SegGamma   !< 
   real(ReKi),     dimension(:)  ,  intent(in) :: SegEpsilon !< 
   if ( vtk_new_ascii_file(filename,'Sgmt') ) then
      call vtk_dataset_polydata(SegPoints(1:3,:))
      call vtk_lines(SegConnct(1:2,:)-1) ! NOTE: VTK indexing at 0
      call vtk_cell_data_init()
      call vtk_cell_data_scalar(SegGamma  ,'Gamma')
      call vtk_cell_data_scalar(SegEpsilon,'Epsilon')
!       call vtk_cell_data_scalar(real(SegConnct(3,:), ReKi),'Age')
      !call vtk_cell_data_scalar(real(SegConnct(4,:), ReKi),'Span')
      call vtk_close_file()
   endif
end subroutine

subroutine WrVTK_Lattice(filename, LatticePoints, LatticeGamma, LatticeData3d)
   use VTK ! for all the vtk_* functions
   character(len=*), intent(in)                         :: filename
   real(Reki), dimension(:,:,:), intent(in  )           :: LatticePoints !< Array of points 3 x nSpan x nDepth
   real(Reki), dimension(:,:), intent(in  )             :: LatticeGamma  !< Array of            nSpan x nDepth
   real(Reki), dimension(:,:,:), intent(in  ), optional :: LatticeData3d !< Array of n x nSpan x nDepth KEEP ME
   !
   integer(IntKi), dimension(:,:), allocatable :: Connectivity
   real(ReKi), dimension(:,:), allocatable     :: Points

   CALL LatticeToPanlConnectivity(LatticePoints, Connectivity, Points)

   if ( vtk_new_ascii_file(filename,'')) then
      call vtk_dataset_polydata(Points)
      call vtk_quad(Connectivity)
      call vtk_cell_data_init()
      call vtk_cell_data_scalar(LatticeGamma,'Gamma')
      if (present(LatticeData3d)) then
         call vtk_point_data_init()
         call vtk_point_data_vector(LatticeData3d,'Uconv')
      endif
      call vtk_close_file()
   endif

end subroutine WrVTK_Lattice

subroutine LatticeToPanlConnectivity(LatticePoints, Connectivity, Points)
   real(Reki), dimension(:,:,:), intent(in   )  :: LatticePoints  !< Array of points 3 x nSpan x nDepth
   integer(IntKi), dimension(:,:), allocatable :: Connectivity
   real(ReKi), dimension(:,:), allocatable     :: Points
   ! Local
   integer(IntKi) :: nSpan, nDepth
   integer(IntKi) :: iSpan, iDepth, k
   nSpan  = size(LatticePoints,2)
   nDepth = size(LatticePoints,3)

   if (allocated(Connectivity)) deallocate(Connectivity)
   allocate(Connectivity(1:4, 1:(nSpan-1)*(nDepth-1)))
   if (allocated(Points)) deallocate(Points)
   allocate(Points(1:3, 1:nSpan*nDepth))

   k=1
   do iDepth=1,nDepth-1; do iSpan=1,nSpan-1
      Connectivity(1,k)=(iDepth-1)*nSpan+(iSpan-1)
      Connectivity(2,k)=(iDepth-1)*nSpan+(iSpan )
      Connectivity(3,k)=(iDepth  )*nSpan+(iSpan)
      Connectivity(4,k)=(iDepth  )*nSpan+(iSpan-1)
      k=k+1
   enddo; enddo

   k=1
   do iDepth=1,nDepth; do iSpan=1,nSpan
      Points(1:3,k) = LatticePoints(1:3,iSpan,iDepth)
      k=k+1
   enddo; enddo

!     do iWing=1,p%NumBlades
!         if ( vtk_new_ascii_file(trim(filename),Label) ) then
!             ! Buffer for points
!             k=1; do iNW=1,nNW; do iSpan=1,nSpan
!                 Buffer(1:3,k) = Misc%NWake%r_nearj(1:3,iSpan,iNW,iWing)
!                 k=k+1
!             enddo; enddo
!             call vtk_dataset_polydata(Buffer)
!             call vtk_quad(Connectivity)
!             call vtk_cell_data_init()
!             ! Buffer for Gammas m1
!             k=1; do iNW=1,(nNW-1); do iSpan=1,(nSpan-1)
!                 if (iSpan<p%NumBlNds_start) then
!                     Buffer1d(k)=0
!                 else if (iSpan==p%NumBlNds_start) then
!                     Buffer1d(k)=-Misc%NWake%Gamma_nearjm1(iNW,iSpan,iWing)
!                 else
!                     Buffer1d(k)=-Misc%NWake%Gamma_nearjm1(iNW,iSpan,iWing)+Buffer1d(k-1)
!                 endif
!                 k=k+1
!             enddo; enddo
!             call vtk_cell_data_scalar(Buffer1d,'Gamma_NW_p1')
end subroutine



END MODULE FVW_IO
