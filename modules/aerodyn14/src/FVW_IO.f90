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
   CALL ReadCom(UnIn,FileName,                  'General option header', ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL ReadVar(UnIn,FileName,Inp%IntMethod    ,'Integration method' ,'',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%FreeWake     ,'FreeWake'         ,'',ErrStat2,ErrMsg2); if(Failed())return
   !------------------------ CIRCULATION SPECIFICATIONS  -------------------------------------------
   CALL ReadCom(UnIn,FileName,                  'Circulation specification header', ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL ReadVar(UnIn,FileName,Inp%CirculationMethod,'CirculationMethod','',ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%CirculationFile  ,'CirculationFile'  ,'',ErrStat2,ErrMsg2); if(Failed())return

   ! Post pro and validation of inputs
   if (PathIsRelative(Inp%CirculationFile)) Inp%CirculationFile = TRIM(PriPath)//TRIM(Inp%CirculationFile)

   if (Check(.not.(ANY((/idCircNoFlowThrough,idCircPrescribed/)==Inp%CirculationMethod)), 'Circulation method not implemented')) return

   if (Check( Inp%IntMethod/=idEuler1 , 'Time integration method not implemented')) return

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
   integer :: nSpan, nNW, nWings, nFW
   integer :: k
   real(ReKi), dimension(:,:), allocatable :: Buffer2d
   character(1), dimension(3) :: I2ABC =(/'A','B','C'/)
   !
   integer(IntKi),dimension(:,:), allocatable :: SegConnct !< Segment connectivity
   real(ReKi),    dimension(:,:), allocatable :: SegPoints !< Segment Points
   real(ReKi),    dimension(:)  , allocatable :: SegGamma  !< Segment Circulation
   integer(IntKi) :: iHeadC, iHeadP, nC, nP
   !real(ReKi),    dimension(:),   allocatable :: SegSmooth !< 

   ! TimeStamp
   write(Tstr, '(i' // trim(Num2LStr(Twidth)) //'.'// trim(Num2LStr(Twidth)) // ')') VTKcount

   nSpan  = p%nSpan
   nWings = p%nWings
   nNW    = m%nNW
   nFW    = m%nFW
   ! --------------------------------------------------------------------------------}
   ! --- Blade 
   ! --------------------------------------------------------------------------------{
   ! --- Blade Quarter chord points (AC)
   do iW=1,nWings
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
   allocate(Buffer2d(1,nSpan))
   do iW=1,nWings
      write(Label,'(A,A)') 'LL.Bld', i2ABC(iW)
      Filename = TRIM(FileRootName)//'.'//trim(Label)//'.'//Tstr//'.vtk'
      Buffer2d(1,:)=m%Gamma_LL(:,iW)
      call WrVTK_Lattice(FileName, m%r_LL(1:3,:,:,iW), Buffer2d)
   enddo
   ! --------------------------------------------------------------------------------}
   ! --- Near wake 
   ! --------------------------------------------------------------------------------{
   ! --- Near wake panels
   do iW=1,nWings
      write(Label,'(A,A)') 'NW.Bld', i2ABC(iW)
      Filename = TRIM(FileRootName)//'.'//trim(Label)//'.'//Tstr//'.vtk'
      call WrVTK_Lattice(FileName, x%r_NW(1:3,:,1:m%nNW+1,iW), x%Gamma_NW(:,1:m%nNW,iW))
   enddo
   ! --------------------------------------------------------------------------------}
   ! --- Far wake 
   ! --------------------------------------------------------------------------------{
   ! --- Far wake panels
   !do iW=1,nWings
   !   write(Label,'(A,A)') 'FW.Bld', i2ABC(iW)
   !   Filename = TRIM(FileRootName)//'.'//trim(Label)//'.'//Tstr//'.vtk'
   !   call WrVTK_Lattice(FileName, x%r_FW(1:3,:,:,iW), x%Gamma_FW(:,:,iW))
   !enddo
   ! --------------------------------------------------------------------------------}
   ! --- All Segments
   ! --------------------------------------------------------------------------------{
   nP =      nWings * (  (nSpan+1)*(nNW+1)            )
   nC =      nWings * (2*(nSpan+1)*(nNW+1)-nSpan-nNW-2)
!    nP = nP + nWings * (nSpan+1)*2
!    nC = nC + nWings * (2*(nSpan+1)*(2)-nSpan-1-2)
   allocate(SegConnct(1:2,1:nC)); SegConnct=-1
   allocate(SegPoints(1:3,1:nP)); SegPoints=-1
   allocate(SegGamma (1:nC)); SegGamma =-1
   iHeadP=1
   iHeadC=1
   do iW=1,nWings
      CALL LatticeToSegments(x%r_NW(1:3,:,1:m%nNW+1,iW), x%Gamma_NW(:,1:m%nNW,iW), SegPoints, SegConnct, SegGamma, iHeadP, iHeadC )
   enddo
!    if (allocated(Buffer2d)) deallocate(Buffer2d)
!    allocate(Buffer2d(1,nSpan))
!    do iW=1,nWings
!       Buffer2d(1,:)=m%Gamma_LL(:,iW)
!       CALL LatticeToSegments(m%r_LL(1:3,:,1:2,iW), Buffer2d, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC )
!    enddo
   Filename = TRIM(FileRootName)//'.AllSeg.'//Tstr//'.vtk'
   CALL WrVTK_Segments(Filename, SegPoints, SegConnct, SegGamma) 

   if ((iHeadP-1)/=nP) then
      print*,'IO: Number of points wrongly estimated',nP, iHeadP-1
!       STOP
   endif
   if ((iHeadC-1)/=nC) then
      print*,'IO: Number of segments wrongly estimated',nC, iHeadC-1
!       STOP
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
      call vtk_cell_data_scalar(SegGamma,'SegGamma')
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
      call vtk_cell_data_scalar(LatticeGamma,'Gamma_NW')
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
