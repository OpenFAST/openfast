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
   integer :: iLine
   real(ReKi) :: TODO_Re
   character(1024)      :: PriPath                                         ! the path to the primary input file
   character(1024)      :: line                                            ! string to temporarially hold value of read line
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
   CALL ReadCom(UnIn,FileName,                         'General option header', ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL ReadVar(UnIn,FileName,Inp%IntMethod           ,'Integration method' ,'',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%FreeWakeStart       ,'FreeWakeStart'      ,'',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%FullCirculationStart,'FullCirculationStart'  ,'',ErrStat2,ErrMsg2); if(Failed())return
   !------------------------ CIRCULATION SPECIFICATIONS  -------------------------------------------
   CALL ReadCom(UnIn,FileName,                  'Circulation specification header', ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL ReadVar(UnIn,FileName,Inp%CirculationMethod ,'CirculationMethod' ,'',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%CircSolvConvCrit  ,'CircSolvConvCrit ' ,'',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%CircSolvRelaxation,'CircSolvRelaxation','',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%CircSolvMaxIter   ,'CircSolvMaxIter'   ,'',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%CircSolvPolar     ,'CircSolvPolar'   ,'',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%CirculationFile   ,'CirculationFile'   ,'',ErrStat2,ErrMsg2); if(Failed())return
   !------------------------ WAKE OPTIONS -------------------------------------------
   CALL ReadCom(UnIn,FileName,                  'Wake options header', ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL ReadVar(UnIn,FileName,Inp%nNWPanels     ,'nNWPanels'       ,'',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%nFWPanels     ,'nFWPanels'       ,'',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%nFWPanelsFree ,'nFWPanelsFree'   ,'',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%RegFunction   ,'RegFunction'     ,'',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%WakeRegMethod ,'WakeRegMethod'   ,'',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%WakeRegFactor ,'WakeRegFactor'   ,'',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%WingRegFactor ,'WingRegFactor'   ,'',ErrStat2,ErrMsg2); if(Failed())return
   !------------------------ OUTPUT OPTIONS -----------------------------------------
   CALL ReadCom(UnIn,FileName,                  'Output options header', ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL ReadVar(UnIn,FileName,Inp%WrVTK       , 'WrVTK'              ,'',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%VTKBlades   , 'VTKBlades'          ,'',ErrStat2,ErrMsg2); if(Failed())return
   !------------------------ HACK OPTIONS -----------------------------------------
   CALL ReadCom(UnIn,FileName,                  'Hack options header', ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL ReadVar(UnIn,FileName,Inp%Uinf        , 'Uinf'   ,'',ErrStat2,ErrMsg2); if(Failed())return

   ! --- Validation of inputs
   if (PathIsRelative(Inp%CirculationFile)) Inp%CirculationFile = TRIM(PriPath)//TRIM(Inp%CirculationFile)

   if (Check(.not.(ANY((/idCircPrescribed,idCircPolarData/)==Inp%CirculationMethod)), 'Circulation method not implemented')) return

   if (Check( Inp%IntMethod/=idEuler1 , 'Time integration method not implemented')) return

   if (Check( Inp%nNWPanels<0     , 'Number of near wake panels must be >=0')) return
   if (Check( Inp%nFWPanels<0     , 'Number of far wake panels must be >=0')) return
   if (Check( Inp%nFWPanelsFree<0 , 'Number of free far wake panels must be >=0')) return
   if (Check( Inp%nFWPanelsFree>Inp%nFWPanels , 'Number of free far wake panels must be <=Number of far wake panels')) return

   if (Check(.not.(ANY(idRegVALID      ==Inp%RegFunction  )), 'Regularization function not implemented')) return
   if (Check(.not.(ANY(idRegMethodVALID==Inp%WakeRegMethod)), 'Wake regularization method not implemented')) return
   if (Check(Inp%WakeRegFactor<0                            , 'Wake regularization factor should be positive')) return
   if (Check(Inp%WingRegFactor<0                            , 'Wing regularization factor should be positive')) return

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
   integer :: iSeg
   integer :: iSpan, iNW, iFW
   integer :: k
   real(ReKi), dimension(:,:), allocatable :: Buffer2d
   character(1), dimension(3) :: I2ABC =(/'A','B','C'/)
   !
   integer(IntKi),dimension(:,:), allocatable :: SegConnct !< Segment connectivity
   real(ReKi),    dimension(:,:), allocatable :: SegPoints !< Segment Points
   real(ReKi),    dimension(:)  , allocatable :: SegGamma  !< Segment Circulation
   integer(IntKi) :: iHeadC, iHeadP, nC, nP
   !real(ReKi),    dimension(:),   allocatable :: SegSmooth !< 

   print*,'------------------------------------------------------------------------------'
   print'(A,L1,A,I0,A,I0,A,I0)','VTK Output  -      First call ',m%FirstCall, '                                nNW:',m%nNW,' nFW:',m%nFW,'  i:',VTKCount
   !
   call set_vtk_binary_format(.false.)

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
         call WrVTK_Lattice(FileName, m%r_LL(1:3,:,1:2,iW), m%Gamma_LL(:,iW:iW))
      else
         call WrVTK_Lattice(FileName, x%r_NW(1:3,:,1:m%nNW+1,iW), x%Gamma_NW(:,1:m%nNW,iW))
      endif
   enddo
   ! --------------------------------------------------------------------------------}
   ! --- Far wake 
   ! --------------------------------------------------------------------------------{
   ! --- Far wake panels
   do iW=1,p%VTKBlades
      write(Label,'(A,A)') 'FW.Bld', i2ABC(iW)
      Filename = TRIM(FileRootName)//'.'//trim(Label)//'.'//Tstr//'.vtk'
      call WrVTK_Lattice(FileName, x%r_FW(1:3,1:FWnSpan+1,1:m%nFW+1,iW), x%Gamma_FW(1:FWnSpan,1:m%nFW,iW))
   enddo
   ! --------------------------------------------------------------------------------}
   ! --- All Segments
   ! --------------------------------------------------------------------------------{
   nP =      p%nWings * (  (p%nSpan+1)*(m%nNW+1)            )
   nC =      p%nWings * (2*(p%nSpan+1)*(m%nNW+1)-(p%nSpan+1)-(m%nNW+1))
   if (m%nFW>0) then
      nP = nP + p%nWings * (  (FWnSpan+1)*(m%nFW+1) )
      nC = nC + p%nWings * (2*(FWnSpan+1)*(m%nFW+1)-(FWnSpan+1)-(m%nFW+1))  
   endif
   allocate(SegConnct(1:2,1:nC)); SegConnct=-1
   allocate(SegPoints(1:3,1:nP)); SegPoints=-1
   allocate(SegGamma (1:nC)); SegGamma =-1
   iHeadP=1
   iHeadC=1
   do iW=1,p%nWings
      if (m%nNW==1) then ! Small Hack - At t=0, NW not set, but first NW panel is the LL panel
         CALL LatticeToSegments(m%r_LL(1:3,:,1:2,iW), m%Gamma_LL(:,iW:iW), 1, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC )
      else
         CALL LatticeToSegments(x%r_NW(1:3,:,1:m%nNW+1,iW), x%Gamma_NW(:,1:m%nNW,iW), 1, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC )
      endif
   enddo
   if (m%nFW>0) then !TODO TODO TODO
      do iW=1,p%nWings
         CALL LatticeToSegments(x%r_FW(1:3,:,1:m%nFW+1,iW), x%Gamma_FW(:,1:m%nFW,iW), 1, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC )
      enddo
   endif
   Filename = TRIM(FileRootName)//'.AllSeg.'//Tstr//'.vtk'
   CALL WrVTK_Segments(Filename, SegPoints, SegConnct, SegGamma) 

   if ((iHeadP-1)/=nP) then
      print*,'IO: Number of points wrongly estimated',nP, iHeadP-1
      STOP
   endif
   if ((iHeadC-1)/=nC) then
      print*,'IO: Number of segments wrongly estimated',nC, iHeadC-1
      STOP
   endif

end subroutine WrVTK_FVW


subroutine WrVTK_Segments(filename, SegPoints, SegConnct, SegGamma) 
   use VTK
   character(len=*),intent(in)                 :: filename
   real(ReKi), dimension(:,:),      intent(in) :: SegPoints !< 
   integer(IntKi), dimension(:,:),  intent(in) :: SegConnct !< 
   real(ReKi),     dimension(:)  ,  intent(in) :: SegGamma !< 
   if ( vtk_new_ascii_file(filename,'Sgmt') ) then
      call vtk_dataset_polydata(SegPoints(1:3,:))
      call vtk_lines(SegConnct(1:2,:)-1) ! NOTE: VTK indexing at 0
      call vtk_cell_data_init()
      call vtk_cell_data_scalar(SegGamma,'Gamma')
      !call vtk_point_data_init()
      !call vtk_point_data_vector(Sgmt%UconvP(1:3,1:Sgmt%nP_Storage),'Uconv')
      call vtk_close_file()
   endif
end subroutine

subroutine WrVTK_Lattice(filename, LatticePoints, LatticeGamma, LatticeData3d)
   use VTK ! for all the vtk_* functions
   character(len=*), intent(in)                         :: filename
   real(Reki), dimension(:,:,:), intent(in  )           :: LatticePoints !< Array of points 3 x nSpan x nDepth
   real(Reki), dimension(:,:), intent(in  )             :: LatticeGamma  !< Array of            nSpan x nDepth
   real(Reki), dimension(:,:,:), intent(in  ), optional :: LatticeData3d !< Array of n x nSpan x nDepth
   !
   integer(IntKi), dimension(:,:), allocatable :: Connectivity
   real(ReKi), dimension(:,:), allocatable     :: Points

   CALL LatticeToPanlConnectivity(LatticePoints, Connectivity, Points)

   if ( vtk_new_ascii_file(filename,'')) then
      call vtk_dataset_polydata(Points)
      call vtk_quad(Connectivity)
      call vtk_cell_data_init()
      call vtk_cell_data_scalar(LatticeGamma,'Gamma')
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

!     do iWing=1,p%NumBl
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
