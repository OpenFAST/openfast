module FVW_IO

  USE FVW_Types
  USE FVW_Subs
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
   integer:: iWing
   character(1024)                       :: FileName
   character(255)                        :: Label
   character(Twidth)                     :: Tstr          ! string for current VTK write-out step (padded with zeros)
   real(ReKi), dimension(:,:), allocatable :: Buffer
   real(ReKi), dimension(:), allocatable :: Buffer1d
   integer(IntKi), dimension(:,:), allocatable :: Connectivity
   integer :: iSeg
   integer :: iSpan, iNW, iFW
   integer :: nSpan, nNW, nWings, nFW
   integer :: k
   character(1), dimension(3) :: I2ABC =(/'A','B','C'/)

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
   do iWing=1,nWings
      write(Label,'(A,A)') 'BldPointCP.Bld', i2ABC(iWing)
      Filename = TRIM(FileRootName)//'.'//trim(Label)//'.'//Tstr//'.vtk'
      if ( vtk_new_ascii_file(trim(filename),Label) ) then
         call vtk_dataset_polydata(m%CP_LL(1:3,1:p%nSpan,iWing))
         call vtk_point_data_init()
         call vtk_point_data_scalar(m%Gamma_ll(    1:p%nSpan,iWing),'Gamma_ll')
         call vtk_point_data_vector(m%Vind_ll (1:3,1:p%nSpan,iWing),'Vind_ll')
         call vtk_point_data_vector(m%Vtot_ll (1:3,1:p%nSpan,iWing),'Vtot_ll')
         call vtk_point_data_vector(m%Vstr_ll (1:3,1:p%nSpan,iWing),'Vstr_ll')
         call vtk_point_data_vector(m%Vwnd_ll (1:3,1:p%nSpan,iWing),'Vwnd_ll')
         call vtk_point_data_vector(m%Tang    (1:3,1:p%nSpan,iWing),'Tangent')
         call vtk_point_data_vector(m%Norm    (1:3,1:p%nSpan,iWing),'Normal')
         call vtk_point_data_vector(m%Orth    (1:3,1:p%nSpan,iWing),'Orth')
         call vtk_close_file()
      endif
   enddo
end subroutine WrVTK_FVW

END MODULE FVW_IO
