module FVW_IO

  USE FVW_Types
  USE FVW_Subs
  use FVW_VortexTools
  implicit none

contains

! ==============================================================================
!> Reads the input file for FVW
SUBROUTINE FVW_ReadInputFile( FileName, p, m, Inp, ErrStat, ErrMsg )
   character(len=*),             intent(in)    :: FileName !< Input file name for FVW
   type(FVW_ParameterType ),     intent(inout) :: p        !< Parameters
   type(FVW_MiscVarType),        intent(inout) :: m        !< Misc
   type(FVW_InputFile),          intent(out)   :: Inp      !< Data stored in the module's input file
   integer(IntKi),               intent(  out) :: ErrStat  !< Error status of the operation
   character(*),                 intent(  out) :: ErrMsg   !< Error message if ErrStat /= ErrID_None
   ! Local variables
   character(1024)      :: PriPath                         ! the path to the primary input file
   character(1024)      :: sDummy, sLine, Key, Val         ! string to temporarially hold value of read line 
   integer(IntKi)       :: UnIn, i
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
   CALL ReadCom        (UnIn,FileName,                        '--- General option header'                                  , ErrStat2,ErrMsg2); if(Failed()) return
   CALL ReadVarWDefault(UnIn,FileName,Inp%IntMethod           ,'Integration method' ,'', idEuler1                       , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%DTfvw               ,'DTfvw'              ,'',  p%DTaero                      , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%FreeWakeStart       ,'FreeWakeStart'       ,'', 0.0_ReKi                      , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%FullCirculationStart,'FullCirculationStart','', real(20.0_ReKi*Inp%DTfvw,ReKi), ErrStat2,ErrMsg2); if(Failed())return
   !------------------------ CIRCULATION SPECIFICATIONS  -------------------------------------------
   CALL ReadCom(UnIn,FileName,                               '--- Circulation specification header'  , ErrStat2, ErrMsg2 ); if(Failed()) return
   CALL ReadVarWDefault(UnIn,FileName,Inp%CirculationMethod ,'CirculationMethod' ,'', idCircPolarData, ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%CircSolvConvCrit  ,'CircSolvConvCrit ' ,'', 0.001          , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%CircSolvRelaxation,'CircSolvRelaxation','', 0.1            , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%CircSolvMaxIter   ,'CircSolvMaxIter'   ,'', 30             , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar(UnIn,FileName,Inp%CirculationFile   ,'CirculationFile'   ,'',ErrStat2,ErrMsg2); if(Failed())return
   !------------------------ WAKE OPTIONS -------------------------------------------
   CALL ReadCom        (UnIn,FileName,                        '=== Separator'                         , ErrStat2,ErrMsg2); if(Failed()) return
   CALL ReadCom        (UnIn,FileName,                        '--- Wake options header'               , ErrStat2,ErrMsg2); if(Failed()) return
   CALL ReadCom        (UnIn,FileName,                        '--- Wake extent header'                , ErrStat2,ErrMsg2); if(Failed()) return
   CALL ReadVar        (UnIn,FileName,Inp%nNWPanels          ,'nNWPanels'         ,''                 , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar        (UnIn,FileName,Inp%nFWPanels          ,'nFWPanels'         ,''                 , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%nFWPanelsFree      ,'nFWPanelsFree'     ,'', Inp%nFWPanels  , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%FWShedVorticity    ,'FWShedVorticity'   ,'', .False.        , ErrStat2,ErrMsg2); if(Failed())return

   CALL ReadCom        (UnIn,FileName,                        '--- Wake regularization header'        , ErrStat2,ErrMsg2); if(Failed()) return
   CALL ReadVarWDefault(UnIn,FileName,Inp%DiffusionMethod    ,'DiffusionMethod'   ,'',idDiffusionNone , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%RegDeterMethod     ,'RegDeterMethod'    ,'',idRegDeterConstant, ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%RegFunction        ,'RegFunction'       ,'',idRegVatistas   , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%WakeRegMethod      ,'WakeRegMethod'     ,'',idRegConstant   , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar        (UnIn,FileName,Inp%WakeRegParam       ,'WakeRegParam'      ,''                 , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar        (UnIn,FileName,Inp%WingRegParam       ,'WingRegParam'      ,''                 , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%CoreSpreadEddyVisc ,'CoreSpreadEddyVisc','',100.0_ReKi      , ErrStat2,ErrMsg2); if(Failed())return

   CALL ReadCom        (UnIn,FileName,                        '--- Wake treatment header'             , ErrStat2,ErrMsg2); if(Failed()) return
   CALL ReadVarWDefault(UnIn,FileName,Inp%TwrShadowOnWake    ,'TwrShadowOnWake'   ,'',.false.         , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%ShearModel         ,'ShearModel'        ,'',idShearNone     , ErrStat2,ErrMsg2); if(Failed())return

   CALL ReadCom        (UnIn,FileName,                        '--- Speed up header      '             , ErrStat2,ErrMsg2); if(Failed()) return
   CALL ReadVarWDefault(UnIn,FileName,Inp%VelocityMethod     ,'VelocityMethod'    ,'',idVelocityBasic , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%TreeBranchFactor   ,'TreeBranchFactor'  ,'',2.0_ReKi        , ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%PartPerSegment     ,'PartPerSegment'    ,'',  1             , ErrStat2,ErrMsg2); if(Failed())return
!    Inp%TwrShadowOnWake  = .False.
!    Inp%VelocityMethod   = idVelocityBasic
!    Inp%TreeBranchFactor = 3.0_ReKi
!    Inp%PartPerSegment   = 1
   !------------------------ OUTPUT OPTIONS -----------------------------------------
   CALL ReadCom        (UnIn,FileName,                  '=== Separator'                      ,ErrStat2,ErrMsg2); if(Failed()) return
   CALL ReadCom        (UnIn,FileName,                  '--- Output options header'          ,ErrStat2,ErrMsg2); if(Failed()) return
   CALL ReadVarWDefault(UnIn,FileName,Inp%WrVTK       , 'WrVTK'              ,'',     0      ,ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%VTKBlades   , 'VTKBlades'          ,'',     1      ,ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVarWDefault(UnIn,FileName,Inp%VTKCoord    , 'VTKCoord'           ,'',     1      ,ErrStat2,ErrMsg2); if(Failed())return
   CALL ReadVar        (UnIn,FileName,sDummy          , 'VTK_fps'            ,''             ,ErrStat2,ErrMsg2); if(Failed())return
   Inp%DTvtk = Get_DTvtk( sDummy, p%DTaero, Inp%DTfvw )

   CALL ReadVarWDefault(UnIn,FileName,p%nGridOut      , 'nGridOut'           ,'',     0      ,ErrStat2,ErrMsg2);
   if (ErrStat2/=ErrID_None) then
      call WarnSyntax('Grid output missing')
      p%nGridOut = 0 ! Important
   else
      allocate(m%GridOutputs(p%nGridOut), stat=ErrStat2);
      CALL ReadCom (UnIn,FileName,  'GridOutHeaders', ErrStat2,ErrMsg2); if(Failed()) return
      CALL ReadCom (UnIn,FileName,  'GridOutUnits', ErrStat2,ErrMsg2); if(Failed()) return
      do i =1, p%nGridOut 
         ErrMsg2='Error reading OLAF grid outputs line '//trim(num2lstr(i))
         read(UnIn, fmt='(A)', iostat=ErrStat2) sLine  ; if(Failed()) return
         call ReadGridOut(sLine, m%GridOutputs(i)); if(Failed()) return
         if (Check(m%GridOutputs(i)%nx<1, 'Grid output nx needs to be >=1')) return
         if (Check(m%GridOutputs(i)%ny<1, 'Grid output ny needs to be >=1')) return
         if (Check(m%GridOutputs(i)%nz<1, 'Grid output nz needs to be >=1')) return
      enddo
   endif

   ! --- Advanced Options
   ! NOTE: no error handling since this is for debug
   ! Default options are typically "true"
   p%InductionAtCP = .true.  ! Compute the induced velocities at Control Points, otherwise, at nodes
   p%WakeAtTE      = .true.  ! The wake starts at the trailing edge, otherwise, directly at the lifting line
   p%Induction     = .true.  ! Compute induced velocities, otherwise 0 induced velocities on the lifting line!
   p%DStallOnWake  = .false.
   CALL ReadCom(UnIn,FileName,                  '=== Separator'                      ,ErrStat2,ErrMsg2); 
   CALL ReadCom(UnIn,FileName,                  '--- Advanced options header'        ,ErrStat2,ErrMsg2);
   if(ErrStat2==ErrID_None) then
      call WrScr(' - Reading advanced options for OLAF:')
      do while(ErrStat2==ErrID_None)
         read(UnIn, '(A)',iostat=ErrStat2) sDummy
         call Conv2UC(sDummy)  ! to uppercase
         if (index(sDummy, 'INDUCTIONATCP')>1) then
            read(sDummy, '(L1)') p%InductionAtCP
            print*,'   >>> InductionAtCP',p%InductionAtCP
         elseif (index(sDummy, 'WAKEATTE')>1) then
            read(sDummy, '(L1)') p%WakeAtTE
            print*,'   >>> WakeAtTE     ',p%WakeAtTE
         elseif (index(sDummy, 'DSTALLONWAKE')>1) then
            read(sDummy, '(L1)') p%DStallOnWake
            print*,'   >>> DStallOnWake ',p%DStallOnWake
         elseif (index(sDummy, 'INDUCTION')>1) then
            read(sDummy, '(L1)') p%Induction
            print*,'   >>> Induction    ',p%Induction
         else
            print*,'   >>> Line ignored, starting with'//trim(sDummy)
         endif
      enddo
   endif


   ! --- Validation of inputs
   if (PathIsRelative(Inp%CirculationFile)) Inp%CirculationFile = TRIM(PriPath)//TRIM(Inp%CirculationFile)

   if (Check(.not.(ANY(idCircVALID ==Inp%CirculationMethod)), 'Circulation method (CircSolvingMethod) not implemented: '//trim(Num2LStr(Inp%CirculationMethod)))) return
   if (Check(.not.(ANY(idIntMethodVALID==Inp%IntMethod    )) , 'Time integration method (IntMethod) not yet implemented. Use Euler 1st order method for now.')) return
   if (Check(.not.(ANY(idDiffusionVALID==Inp%DiffusionMethod)) , 'Diffusion method (DiffusionMethod) not implemented: '//trim(Num2LStr(Inp%DiffusionMethod)))) return
   if (Check(.not.(ANY(idRegDeterVALID ==Inp%RegDeterMethod))  , 'Regularization determination method (RegDeterMethod) not yet implemented: '//trim(Num2LStr(Inp%RegDeterMethod)))) return
   if (Check(.not.(ANY(idRegVALID      ==Inp%RegFunction  )), 'Regularization function (RegFunction) not implemented: '//trim(Num2LStr(Inp%RegFunction)))) return
   if (Check(.not.(ANY(idRegMethodVALID==Inp%WakeRegMethod)), 'Wake regularization method (WakeRegMethod) not implemented: '//trim(Num2LStr(Inp%WakeRegMethod)))) return
   if (Check(.not.(ANY(idShearVALID    ==Inp%ShearModel   )), 'Shear model (ShearModel) not valid: '//trim(Num2LStr(Inp%ShearModel)))) return
   if (Check(.not.(ANY(idVelocityVALID ==Inp%VelocityMethod    )), 'Velocity method (VelocityMethod) not valid: '//trim(Num2LStr(Inp%VelocityMethod)))) return

   if (Check( Inp%DTfvw < p%DTaero, 'DTfvw must be >= DTaero from AD15.')) return
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

   subroutine WarnSyntax(msg)
      character(len=*), intent(in) :: msg
      call WrScr('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
      call WrScr('OLAF input file is not at its latest format')
      call WrScr('Error: '//trim(msg))
      call WrScr('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
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

   subroutine ReadGridOut(sLine, GridOut)
      character(len=*),  intent(in)  :: sLine  !< full line
      type(GridOutType), intent(out) :: GridOut
      character(255) :: StrArray(14) ! Array of strings extracted from line
      real(ReKi) :: DummyFloat
      ! Convert line to array of strings
      StrArray(:)='';
      CALL ReadCAryFromStr(sLine, StrArray, 14, 'StrArray', 'StrArray', ErrStat2, ErrMsg2)! NOTE:No Error handling!
      ! Default to error
      ErrStat2=ErrID_Fatal
      ErrMsg2='Error reading OLAF grid outputs line: '//trim(sLine)
      ! Name
      GridOut%name =StrArray(1) 
      ! Type
      if (.not. is_int    (StrArray(2), GridOut%type  ) ) then
         ErrMsg2=trim(ErrMsg2)//achar(13)//achar(10)//'GridType needs to be an integer.'
         return
      endif
      ! tStart
      call Conv2UC( StrArray(3) )
      if ( index(StrArray(3), "DEFAULT" ) == 1 ) then
         GridOut%tStart  = 0.0_ReKi
      else
         if (.not. is_numeric(StrArray(3), GridOut%tStart) ) then 
            ErrMsg2=trim(ErrMsg2)//achar(13)//achar(10)//'TStart needs to be numeric or "default".'
            return
         endif
      endif
      ! tEnd
      call Conv2UC( StrArray(4) )
      if ( index(StrArray(4), "DEFAULT" ) == 1 ) then
         GridOut%tEnd  = 99999.0_ReKi ! TODO
      else
         if (.not. is_numeric(StrArray(4), GridOut%tEnd) ) then
            ErrMsg2=trim(ErrMsg2)//achar(13)//achar(10)//'TEnd needs to be numeric or "default".'
            return
         endif
      endif
      ! Dtout
      call Conv2UC( StrArray(5) )
      if ( index(StrArray(5), "DEFAULT" ) == 1 ) then
         GridOut%DTout  = p%DTfvw
      else if ( index(StrArray(5), "ALL" ) == 1 ) then
         GridOut%DTout  = p%DTaero
      else
         if (.not. is_numeric(StrArray(5), GridOut%DTout) ) then
            ErrMsg2=trim(ErrMsg2)//achar(13)//achar(10)//'DTout needs to be numeric, "default" or "all".'
            return
         endif
      endif
      ! x,y,z
      ErrMsg2='Error reading OLAF "x" inputs for grid outputs line: '//trim(sLine)
      if (.not. is_numeric(StrArray( 6), GridOut%xStart) ) return
      if (.not. is_numeric(StrArray( 7), GridOut%xEnd  ) ) return
      if (.not. is_int    (StrArray( 8), GridOut%nx    ) ) return
      ErrMsg2='Error reading OLAF "y" inputs for grid outputs line: '//trim(sLine)
      if (.not. is_numeric(StrArray( 9), GridOut%yStart) ) return
      if (.not. is_numeric(StrArray(10), GridOut%yEnd  ) ) return
      if (.not. is_int    (StrArray(11), GridOut%ny    ) ) return
      ErrMsg2='Error reading OLAF "z" inputs for grid outputs line: '//trim(sLine)
      if (.not. is_numeric(StrArray(12), GridOut%zStart) ) return
      if (.not. is_numeric(StrArray(13), GridOut%zEnd  ) ) return
      if (.not. is_int    (StrArray(14), GridOut%nz    ) ) return
      ! Success
      ErrStat2=ErrID_None
      ErrMsg2=''
   end subroutine ReadGridOut

END SUBROUTINE FVW_ReadInputFile

function is_numeric(string, x)
   implicit none
   character(len=*), intent(in) :: string
   real(reki), intent(out) :: x
   logical :: is_numeric
   integer :: e,n
   character(len=12) :: fmt
   x = 0.0_reki
   n=len_trim(string)
   write(fmt,'("(F",I0,".0)")') n
   read(string,fmt,iostat=e) x
   is_numeric = e == 0
end function is_numeric

function is_int(string, x)
   implicit none
   character(len=*), intent(in) :: string
   integer(IntKi), intent(out) :: x
   logical :: is_int
   integer :: e,n
   character(len=12) :: fmt
   x = 0
   n=len_trim(string)
   write(fmt,'("(I",I0,")")') n
   read(string,fmt,iostat=e) x
   is_int = e == 0
end function is_int


!=================================================
!> Export FVW variables to VTK
!! NOTE: when entering this function nNW and nFW has been incremented by 1
subroutine WrVTK_FVW(p, x, z, m, FileRootName, VTKcount, Twidth, bladeFrame, HubOrientation, HubPosition)
   use FVW_VTK ! for all the vtk_* functions
   type(FVW_ParameterType),        intent(in   ) :: p !< Parameters
   type(FVW_ContinuousStateType),  intent(in   ) :: x !< States
   type(FVW_ConstraintStateType),  intent(in   ) :: z !< Constraints
   type(FVW_MiscVarType),          intent(in   ) :: m !< MiscVars
   character(*),    intent(in)           :: FileRootName    !< Name of the file to write the output in (excluding extension)
   integer(IntKi),  intent(in)           :: VTKcount        !< Indicates number for VTK output file (when 0, the routine will also write reference information)
   integer(IntKi),  intent(in)           :: Twidth          !< Number of digits in the maximum write-out step (used to pad the VTK write-out in the filename with zeros)
   logical,         intent(in   )        :: bladeFrame      !< Output in blade coordinate frame
   real(ReKi),optional,dimension(3,3), intent(in) :: HubOrientation
   real(ReKi),optional,dimension(3)  , intent(in) :: HubPosition
   ! local variables
   integer:: iW
   character(1024)                       :: FileName
   character(255)                        :: Label
   character(Twidth)                     :: Tstr          ! string for current VTK write-out step (padded with zeros)
   character(1), dimension(26) :: I2ABC =(/'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/)
   integer(IntKi)       :: nSeg, nSegP, nSegNW
   logical              :: bMirror
   !integer(IntKi)       :: ErrStat2
   !character(ErrMsgLen) :: ErrMsg2
   real(Reki), dimension(:,:,:), allocatable :: Arr3D !<
   real(Reki), dimension(:,:), allocatable :: Arr2D !<

   type(FVW_VTK_Misc)   :: mvtk

   call vtk_misc_init(mvtk)

   if (bladeFrame) then
      if (present(HubOrientation) .and. present(HubPosition)) then
         call set_vtk_coordinate_transform(HubOrientation,HubPosition,mvtk)
      else
         Call ProgAbort('Programming error in WrVTK_FVW call: Cannot use the WrVTK_FVW with bladeFrame==TRUE without the optional arguments of HubOrientation and HubPosition')
      endif
   endif
 
   if (DEV_VERSION) then
      print*,'------------------------------------------------------------------------------'
      print'(A,L1,A,I0,A,I0,A,I0)','VTK Output  -      First call ',m%FirstCall, '                                nNW:',m%nNW,' nFW:',m%nFW,'  i:',VTKCount
   endif
   !
   call set_vtk_binary_format(.false.,mvtk) ! TODO binary fails

   ! TimeStamp
   write(Tstr, '(i' // trim(Num2LStr(Twidth)) //'.'// trim(Num2LStr(Twidth)) // ')') VTKcount

   ! --------------------------------------------------------------------------------}
   ! --- Blade 
   ! --------------------------------------------------------------------------------{
   ! --- Blade Quarter chord points (AC)
   do iW=1,p%VTKBlades
      write(Label,'(A,A)') 'BldPointCP.Bld', i2ABC(iW)
      Filename = TRIM(FileRootName)//'.'//trim(Label)//'.'//Tstr//'.vtk'
      if ( vtk_new_ascii_file(trim(filename),Label,mvtk) ) then
         call vtk_dataset_polydata(m%W(iW)%CP(1:3,1:p%W(iW)%nSpan),mvtk,bladeFrame)
         call vtk_point_data_init(mvtk)
         call vtk_point_data_scalar(z%W(iW)%Gamma_LL(    1:p%W(iW)%nSpan),'Gamma_LL',mvtk)
         call vtk_point_data_vector(m%W(iW)%Vind_CP (1:3,1:p%W(iW)%nSpan),'Vind_CP',mvtk)
         call vtk_point_data_vector(m%W(iW)%Vtot_CP (1:3,1:p%W(iW)%nSpan),'Vtot_CP',mvtk)
         call vtk_point_data_vector(m%W(iW)%Vstr_CP (1:3,1:p%W(iW)%nSpan),'Vstr_CP',mvtk)
         call vtk_point_data_vector(m%W(iW)%Vwnd_CP (1:3,1:p%W(iW)%nSpan),'Vwnd_CP',mvtk)
         call vtk_point_data_vector(m%W(iW)%Tang    (1:3,1:p%W(iW)%nSpan),'Tangent',mvtk)
         call vtk_point_data_vector(m%W(iW)%Norm    (1:3,1:p%W(iW)%nSpan),'Normal',mvtk)
         call vtk_point_data_vector(m%W(iW)%Orth    (1:3,1:p%W(iW)%nSpan),'Orth',mvtk)
         call vtk_close_file(mvtk)
      endif
   enddo
   ! --- Lifting line panels
   ! TODO
   ! do iW=1,p%VTKBlades
   !    write(Label,'(A,A)') 'LL.Bld', i2ABC(iW)
   !    Filename = TRIM(FileRootName)//'.'//trim(Label)//'.'//Tstr//'.vtk'
   !    call WrVTK_Lattice(FileName, mvtk, m%W(iW)%r_LL(1:3,:,:), m%W(iW)%Gamma_LL(:), bladeFrame=bladeFrame)
   ! enddo
   ! --------------------------------------------------------------------------------}
   ! --- Near wake 
   ! --------------------------------------------------------------------------------{
   ! --- Near wake panels
   do iW=1,p%VTKBlades
      write(Label,'(A,A)') 'NW.Bld', i2ABC(iW)
      Filename = TRIM(FileRootName)//'.'//trim(Label)//'.'//Tstr//'.vtk'
      if (m%FirstCall) then ! Small Hack - At t=0, NW not set, but first NW panel is the LL panel
         allocate(Arr3D(3, size(m%dxdt%W(iW)%r_NW,2) ,2)); Arr3D=0.0_ReKi ! Convection velocity
         allocate(Arr2D(size(z%W(iW)%Gamma_LL), 1) )     ; Arr2D=0.0_ReKi ! Gamma
         Arr2D(:,1)=z%W(iW)%Gamma_LL(:)
         call WrVTK_Lattice(FileName, mvtk, m%W(iW)%r_LL(1:3,:,1:2), Arr2D(:,1:1), Arr3D, bladeFrame=bladeFrame)
         deallocate(Arr3D)
         deallocate(Arr2D)
      else
         call WrVTK_Lattice(FileName, mvtk, x%W(iW)%r_NW(1:3,:,1:m%nNW+1), x%W(iW)%Gamma_NW(:,1:m%nNW), m%dxdt%W(iW)%r_NW(:,:,1:m%nNW+1), bladeFrame=bladeFrame)
      endif
   enddo
   ! --------------------------------------------------------------------------------}
   ! --- Far wake 
   ! --------------------------------------------------------------------------------{
   ! --- Far wake panels
   do iW=1,p%VTKBlades
      write(Label,'(A,A)') 'FW.Bld', i2ABC(iW)
      Filename = TRIM(FileRootName)//'.'//trim(Label)//'.'//Tstr//'.vtk'
      call WrVTK_Lattice(FileName, mvtk, x%W(iW)%r_FW(1:3,1:FWnSpan+1,1:m%nFW+1), x%W(iW)%Gamma_FW(1:FWnSpan,1:m%nFW),m%dxdt%W(iW)%r_FW(:,:,1:m%nFW+1), bladeFrame=bladeFrame)
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
   CALL WrVTK_Segments(Filename, mvtk, m%Sgmt%Points(:,1:nSegP), m%Sgmt%Connct(:,1:nSeg), m%Sgmt%Gamma(1:nSeg), m%Sgmt%Epsilon(1:nSeg), bladeFrame) 

   if(.false.) print*,z%W(1)%Gamma_LL(1) ! unused var for now
end subroutine WrVTK_FVW

!> Export Grid velocity field to VTK
subroutine WrVTK_FVW_Grid(p, x, z, m, iGrid, FileRootName, VTKcount, Twidth, HubOrientation, HubPosition)
   use FVW_VortexTools, only: curl_regular_grid
   use FVW_VTK ! for all the vtk_* functions
   type(FVW_ParameterType),        intent(in   ) :: p !< Parameters
   type(FVW_ContinuousStateType),  intent(in   ) :: x !< States
   type(FVW_ConstraintStateType),  intent(in   ) :: z !< Constraints
   type(FVW_MiscVarType), target,  intent(in   ) :: m !< MiscVars
   integer(IntKi),  intent(in)           :: iGrid           !< Grid out index
   character(*),    intent(in)           :: FileRootName    !< Name of the file to write the output in (excluding extension)
   integer(IntKi),  intent(in)           :: VTKcount        !< Indicates number for VTK output file (when 0, the routine will also write reference information)
   integer(IntKi),  intent(in)           :: Twidth          !< Number of digits in the maximum write-out step (used to pad the VTK write-out in the filename with zeros)
   real(ReKi),optional,dimension(3,3), intent(in) :: HubOrientation
   real(ReKi),optional,dimension(3)  , intent(in) :: HubPosition
   ! local variables
   character(1024)   :: FileName
   character(255)    :: Label
   character(Twidth) :: Tstr     ! string for current VTK write-out step (padded with zeros)
   real(ReKi), dimension(3) :: dx
   type(GridOutType), pointer :: g
   type(FVW_VTK_Misc)   :: mvtk

   call vtk_misc_init(mvtk)
   call set_vtk_binary_format(.false.,mvtk) ! TODO binary fails

   ! TimeStamp
   write(Tstr, '(i' // trim(Num2LStr(Twidth)) //'.'// trim(Num2LStr(Twidth)) // ')') VTKcount

   ! --- Grid
   g => m%GridOutputs(iGrid)
   Label=trim(g%name)
   Filename = TRIM(FileRootName)//'.'//trim(Label)//'.'//Tstr//'.vtk'
   if ( vtk_new_ascii_file(trim(filename),Label,mvtk) ) then
      dx(1) = (g%xEnd- g%xStart)/max(g%nx-1,1)
      dx(2) = (g%yEnd- g%yStart)/max(g%ny-1,1)
      dx(3) = (g%zEnd- g%zStart)/max(g%nz-1,1)
      call vtk_dataset_structured_points((/g%xStart, g%yStart, g%zStart/),dx,(/g%nx,g%ny,g%nz/),mvtk)
      call vtk_point_data_init(mvtk)
      call vtk_point_data_vector(g%uGrid(1:3,:,:,:),'Velocity',mvtk) 
      ! Compute vorticity on the fly
      if (g%type==idGridVelVorticity) then
         call curl_regular_grid(g%uGrid, g%omgrid, 1,1,1, g%nx,g%ny,g%nz, dx(1),dx(2),dx(3))
         call vtk_point_data_vector(g%omGrid(1:3,:,:,:),'Vorticity',mvtk) 
      endif
      !
      call vtk_close_file(mvtk)
   endif

   if(.false.) print*,z%W(1)%Gamma_LL(1) ! unused var for now
end subroutine WrVTK_FVW_Grid



subroutine WrVTK_Segments(filename, mvtk, SegPoints, SegConnct, SegGamma, SegEpsilon, bladeFrame) 
   use FVW_VTK
   character(len=*),intent(in)                 :: filename
   type(FVW_VTK_Misc),           intent(inout) :: mvtk       !< miscvars for VTK output
   real(ReKi), dimension(:,:),      intent(in) :: SegPoints  !< 
   integer(IntKi), dimension(:,:),  intent(in) :: SegConnct  !< 
   real(ReKi),     dimension(:)  ,  intent(in) :: SegGamma   !< 
   real(ReKi),     dimension(:)  ,  intent(in) :: SegEpsilon !< 
   logical,                      intent(in   ) :: bladeFrame !< Output in blade coordinate frame
   if ( vtk_new_ascii_file(filename,'Sgmt',mvtk) ) then
      call vtk_dataset_polydata(SegPoints(1:3,:),mvtk,bladeFrame)
      call vtk_lines(SegConnct(1:2,:)-1,mvtk) ! NOTE: VTK indexing at 0
      call vtk_cell_data_init(mvtk)
      call vtk_cell_data_scalar(SegGamma  ,'Gamma',mvtk)
      call vtk_cell_data_scalar(SegEpsilon,'Epsilon',mvtk)
!       call vtk_cell_data_scalar(real(SegConnct(3,:), ReKi),'Age',mvtk)
      !call vtk_cell_data_scalar(real(SegConnct(4,:), ReKi),'Span',mvtk)
      call vtk_close_file(mvtk)
   endif
end subroutine

subroutine WrVTK_Lattice(filename, mvtk, LatticePoints, LatticeGamma, LatticeData3d, bladeFrame)
   use FVW_VTK ! for all the vtk_* functions
   character(len=*), intent(in)                         :: filename
   type(FVW_VTK_Misc),           intent(inout)          :: mvtk          !< miscvars for VTK output
   real(Reki), dimension(:,:,:), intent(in  )           :: LatticePoints !< Array of points 3 x nSpan x nDepth
   real(Reki), dimension(:,:), intent(in  )             :: LatticeGamma  !< Array of            nSpan x nDepth
   real(Reki), dimension(:,:,:), intent(in  ), optional :: LatticeData3d !< Array of n x nSpan x nDepth KEEP ME
   logical,                      intent(in   )          :: bladeFrame    !< Output in blade coordinate frame
   !
   integer(IntKi), dimension(:,:), allocatable :: Connectivity
   real(ReKi), dimension(:,:), allocatable     :: Points

   CALL LatticeToPanlConnectivity(LatticePoints, Connectivity, Points)

   if ( vtk_new_ascii_file(filename,'',mvtk)) then
      call vtk_dataset_polydata(Points,mvtk,bladeFrame)
      call vtk_quad(Connectivity,mvtk)
      call vtk_cell_data_init(mvtk)
      call vtk_cell_data_scalar(LatticeGamma,'Gamma',mvtk)
      if (present(LatticeData3d)) then
         call vtk_point_data_init(mvtk)
         call vtk_point_data_vector(LatticeData3d,'Uconv',mvtk)
      endif
      call vtk_close_file(mvtk)
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

!     do iW=1,p%NumBlades
!         if ( vtk_new_ascii_file(trim(filename),Label,mvtk) ) then
!             ! Buffer for points
!             k=1; do iNW=1,nNW; do iSpan=1,nSpan
!                 Buffer(1:3,k) = Misc%NWake%r_nearj(1:3,iSpan,iNW)
!                 k=k+1
!             enddo; enddo
!             call vtk_dataset_polydata(Buffer,mvtk)
!             call vtk_quad(Connectivity)
!             call vtk_cell_data_init()
!             ! Buffer for Gammas m1
!             k=1; do iNW=1,(nNW-1); do iSpan=1,(nSpan-1)
!                 if (iSpan<p%NumBlNds_start) then
!                     Buffer1d(k)=0
!                 else if (iSpan==p%NumBlNds_start) then
!                     Buffer1d(k)=-Misc%NWake%W(iW)%Gamma_nearjm1(iNW,iSpan)
!                 else
!                     Buffer1d(k)=-Misc%NWake%W(iW)%Gamma_nearjm1(iNW,iSpan)+Buffer1d(k-1)
!                 endif
!                 k=k+1
!             enddo; enddo
!             call vtk_cell_data_scalar(Buffer1d,'Gamma_NW_p1')
end subroutine



END MODULE FVW_IO
