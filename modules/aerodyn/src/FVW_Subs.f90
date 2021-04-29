module FVW_SUBS

   use NWTC_LIBRARY
   use FVW_TYPES
   use FVW_VortexTools
   use FVW_BiotSavart

   implicit none

   ! --- Module parameters
   ! Circulation solving methods
   integer(IntKi), parameter :: idCircPolarData     = 1
   integer(IntKi), parameter :: idCircNoFlowThrough = 2
   integer(IntKi), parameter :: idCircPrescribed    = 3
   integer(IntKi), parameter, dimension(2) :: idCircVALID = (/idCircPolarData, idCircPrescribed /)
   ! Integration method
   integer(IntKi), parameter :: idRK4      = 1 
   integer(IntKi), parameter :: idAB4      = 2
   integer(IntKi), parameter :: idABM4     = 3
   integer(IntKi), parameter :: idPredictor= 4
   integer(IntKi), parameter :: idEuler1   = 5
   integer(IntKi), parameter, dimension(2) :: idIntMethodVALID      = (/idEuler1, idRK4 /)
   ! Diffusion method
   integer(IntKi), parameter :: idDiffusionNone       = 0
   integer(IntKi), parameter :: idDiffusionCoreSpread = 1
   integer(IntKi), parameter :: idDiffusionPSE        = 2
   integer(IntKi), parameter, dimension(1) :: idDiffusionVALID      = (/idDiffusionNone /)
   ! Regularization Method
   integer(IntKi), parameter :: idRegConstant   = 1
   integer(IntKi), parameter :: idRegStretching = 2
   integer(IntKi), parameter :: idRegAge        = 3
   integer(IntKi), parameter, dimension(2) :: idRegMethodVALID      = (/idRegConstant,idRegAge/)
   ! Regularization determination method
   integer(IntKi), parameter :: idRegDeterConstant  = 0
   integer(IntKi), parameter :: idRegDeterAuto    = 1
   integer(IntKi), parameter :: idRegDeterChord     = 2
   integer(IntKi), parameter :: idRegDeterSpan      = 3
   integer(IntKi), parameter, dimension(4) :: idRegDeterVALID      = (/idRegDeterConstant, idRegDeterAuto, idRegDeterChord, idRegDeterSpan /)
   ! Shear model
   integer(IntKi), parameter :: idShearNone   = 0
   integer(IntKi), parameter :: idShearMirror = 1
   integer(IntKi), parameter, dimension(2) :: idShearVALID         = (/idShearNone, idShearMirror /)
   ! Velocity calculation method
   integer(IntKi), parameter :: idVelocityBasic = 1
   integer(IntKi), parameter :: idVelocityTree  = 2
   integer(IntKi), parameter :: idVelocityPart  = 3
   integer(IntKi), parameter, dimension(3) :: idVelocityVALID      = (/idVelocityBasic, idVelocityTree, idVelocityPart /)

   real(ReKi), parameter :: CoreSpreadAlpha = 1.25643 

   ! Implementation 
   integer(IntKi), parameter :: iNWStart=2 !< Index in r%NW where the near wake start (if >1 then the Wing panels are included in r_NW)
   integer(IntKi), parameter :: FWnSpan=1  !< Number of spanwise far wake panels ! TODO make it an input later
   logical       , parameter :: DEV_VERSION=.False.
contains

!==========================================================================
!> Helper function for 1d interpolation (interp1d)
function interpolation_array( xvals, yvals, xi, nOut, nIn )
   integer nOut, nIn, arindx, ilo
   real(ReKi), dimension( nOut ) :: interpolation_array, xi
   real(ReKi), dimension( nIn ) :: xvals, yvals, tmp2, tmp3
   real(ReKi)                         :: tmp1
   ilo = 1
   DO arindx = 1, nOut
      IF ( xi( arindx ) .LT. xvals( 1 )) THEN
         interpolation_array( arindx ) = yvals( 1 ) + ( xi( arindx ) - xvals( 1 )) / &
            & ( xvals( 2 ) - xvals( 1 )) * ( yvals( 2 ) - yvals( 1 ))
      ELSE IF ( xi( arindx ) .GT. xvals( nIn )) THEN
         interpolation_array( arindx ) = yvals( nIn - 1 ) + ( xi( arindx ) - &
            & xvals( nIn - 1 )) / ( xvals( nIn ) - xvals( nIn - 1 )) * &
            & ( yvals( nIn ) - yvals( nIn - 1 ))
      ELSE
         tmp1 = real( xi( arindx ), ReKi)
         tmp2 = real( xvals , ReKi)
         tmp3 = real( yvals , ReKi)
         interpolation_array( arindx ) = InterpBinReal( tmp1, tmp2, tmp3, ilo, nIn )
      END IF
   END DO
END FUNCTION interpolation_array
!==========================================================================

! =====================================================================================
!> Output blade circulation 
subroutine Output_Gamma(CP, Gamma_LL, iWing, iStep, iLabel, iIter)
   real( ReKi ), dimension( :, : ), intent(in   ) :: CP       !< Control Points
   real( ReKi ), dimension( : ),    intent(in   ) :: Gamma_LL !< Circulation on the lifting line
   integer( IntKi ),                intent(in   ) :: iWing    !< Wing index
   integer( IntKi ),                intent(in   ) :: iStep    !< Call ID
   integer( IntKi ),                intent(in   ) :: iLabel    !< Call ID
   integer( IntKi ),                intent(in   ) :: iIter    !< Call ID
   character(len=255) :: filename
   integer :: i
   integer :: iUnit
   real(ReKi) :: norm
   call GetNewUnit(iUnit)
   ! TODO output folder
   CALL MKDIR('Gamma')
   write(filename,'(A,I0,A,I0,A,I0,A,I0,A)')'Gamma/Gamma_step',int(iStep),'_lab',iLabel,'_it',iIter,'_Wing',int(iWing),'.txt'
   OPEN(unit = iUnit, file = trim(filename), status="unknown", action="write")
   write(iUnit,'(A)') 'norm_[m],x_[m],y_[m],z_[m], Gamma_[m^2/s]'
   do i=1,size(Gamma_LL)
      norm=sqrt(CP(1,i)**2+CP(2,i)**2+CP(3,i)**2)
      write(iUnit,'(E14.7,A,E14.7,A,E14.7,A,E14.7,A,E14.7)') norm,',', CP(1,i),',',CP(2,i),',',CP(3,i),',', Gamma_LL(i)
   enddo
   close(iUnit)
endsubroutine Output_Gamma
! =====================================================================================
!> Read a delimited file  containing a circulation and interpolate it on the requested Control Points
!! The input file is a delimited file with one line of header. 
!! Each following line consists of two columns: r/R_[-] and Gamma_[m^2/s]
subroutine ReadAndInterpGamma(CirculationFileName, s_CP_LL, L, Gamma_CP_LL, ErrStat, ErrMsg)
   character(len=*),           intent(in   ) :: CirculationFileName !< Input file to read
   real(ReKi), dimension(:),   intent(in   ) :: s_CP_LL             !< Spanwise location of the lifting CP [m]
   real(ReKi),                 intent(in   ) :: L                   !< Full span of lifting line
   real(ReKi), dimension(:),   intent(out  ) :: Gamma_CP_LL         !< Interpolated circulation of the LL CP
   integer(IntKi),             intent(  out) :: ErrStat             !< Error status of the operation
   character(*),               intent(  out) :: ErrMsg              !< Error message if ErrStat /= ErrID_None
   ! Local
   integer(IntKi)       :: nLines
   integer(IntKi)       :: i
   integer(IntKi)       :: iStat
   integer(IntKi)       :: iUnit
   character(len=1054)  :: line
   integer(IntKi)       :: ErrStat2                                                           ! temporary Error status
   character(ErrMsgLen) :: ErrMsg2                                                            ! temporary Error message
   real(ReKi), dimension(:), allocatable :: sPrescr, GammaPrescr !< Radius
   real(ReKi), parameter :: ReNaN = huge(1.0_ReKi)
   ErrStat = ErrID_None
   ErrMsg  = ''
   ! --- 
   call GetNewUnit(iUnit)
   call OpenFInpFile(iUnit, CirculationFileName, errStat2, errMsg2); if(Failed()) return
   nLines=line_count(iUnit)-1
   ! Read Header
   read(iUnit,*, iostat=errStat2) line ; if(Failed()) return
   ! Read table:  s/L [-], GammaPresc [m^2/s]
   call AllocAry(sPrescr    , nLines, 'sPrecr'    , errStat2, errMsg2); if(Failed()) return
   call AllocAry(GammaPrescr, nLines, 'GammaPrecr', errStat2, errMsg2); if(Failed()) return
   sPrescr     = ReNaN
   GammaPrescr = ReNaN
   do i=1,nLines
      read(iUnit,*, iostat=istat) sPrescr(i), GammaPrescr(i)
      if (istat/=0) then
         errStat2=ErrID_Fatal
         errMsg2='Error occured while reading line '//num2lstr(i+1)//' of circulation file: '//trim(CirculationFileName)
         if(Failed()) return
      endif
   enddo
   if (any(GammaPrescr>=ReNaN).or.any(sPrescr>=ReNaN)) then
      errStat2=ErrID_Fatal
      errMsg2='Not all values were read properly (check the format) while reading the circulation file: '//trim(CirculationFileName)
      if(Failed()) return
   endif
   sPrescr = sPrescr * L
   ! NOTE: TODO TODO TODO THIS ROUTINE PERFORMS NASTY EXTRAPOLATION, SHOULD BE PLATEAUED
   Gamma_CP_LL =  interpolation_array(sPrescr, GammaPrescr, s_CP_LL, size(s_CP_LL), nLines)

   call CleanUp()
contains
   subroutine CleanUp()
      if(allocated(sPrescr)) deallocate(sPrescr)
      if(allocated(GammaPrescr)) deallocate(GammaPrescr)
      if (iUnit>0) close(iUnit)
   end subroutine

   logical function Failed()
      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ReadAndInterpGamma') 
      Failed =  ErrStat >= AbortErrLev
      if (Failed) call CleanUp()
   end function Failed

   !> Counts number of lines in a file
   integer function line_count(iunit)
      integer(IntKi), intent(in) :: iunit
      character(len=1054) :: line
      ! safety for infinite loop..
      integer(IntKi), parameter :: nline_max=100000000 ! 100 M
      integer(IntKi) :: i
      line_count=0
      do i=1,nline_max 
         line=''
         read(iunit,'(A)',END=100)line
         line_count=line_count+1
      enddo
      if (line_count==nline_max) then
         print*,'Error: maximum number of line exceeded'
      endif
      100 if(len(trim(line))>0) then
         line_count=line_count+1
      endif
      rewind(iunit)
   end function

endsubroutine ReadAndInterpGamma
! =====================================================================================

! --------------------------------------------------------------------------------
! --- Mapping functions 
! --------------------------------------------------------------------------------

!> Make sure the First panel of the NW match the last panel of the Trailing edge
!!  - Same position of points
!!  - Same circulation 
subroutine Map_LL_NW(p, m, z, x, ShedScale, ErrStat, ErrMsg )
   type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
   type(FVW_MiscVarType),           intent(in   )  :: m              !< Initial misc/optimization variables
   type(FVW_ConstraintStateType),   intent(in   )  :: z              !< Constraints states
   type(FVW_ContinuousStateType),   intent(inout)  :: x              !< Continuous states
   real(ReKi),                      intent(in)     :: ShedScale      !< Time scaling of shed vorticity
   integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   real(ReKi) :: Gamma_Prev, Gamma_new
   ! Local
   integer(IntKi) :: iSpan , iW
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! First panel of NW is the last lifting line panel
   do iW = 1,p%nWings
      do iSpan = 1,p%nSpan+1
         x%r_NW(1:3, iSpan, iNWStart-1, iW) = m%r_LL(1:3, iSpan, 1, iW)  ! iAge=1
         x%r_NW(1:3, iSpan, iNWStart  , iW) = m%r_LL(1:3, iSpan, 2, iW)  ! iAge=2
      enddo
   enddo
   ! First panel of NW is the last lifting line panel
   do iW = 1,p%nWings
      do iSpan = 1,p%nSpan
         x%Gamma_NW(iSpan, iNWStart-1, iW) = z%Gamma_LL(iSpan,iW)  ! iAge=1
      enddo
   enddo
   ! Circulations are the same on both side of the TE 
   if (p%nNWMax>iNWStart-1) then
      do iW = 1,p%nWings
         do iSpan = 1,p%nSpan
            x%Gamma_NW(iSpan, iNWStart  , iW) = z%Gamma_LL(iSpan,iW)  ! iAge=2
         enddo
      enddo
   endif
   ! When subcycling, we make sure the new circulation progressively ramps up from the old one
   ! NOTE: subcycling needs improvement. 
   !       Frequencies are introduced, even for prescribed circulation, when wake roll up is included
   !       If the wake is not free, the convection velocity is constant and there is no issue.
   !       As a test case, the elliptical wing with constant circulation can be used, with roll up
   !       The error seems to be bigger near the tip/root for this case.
   if(.false.) then
      if ((ShedScale<1.0_ReKi) .and. (m%nNW>=3)) then
         print*,'Scaling'
         do iW = 1,p%nWings
            do iSpan = 1,p%nSpan
               Gamma_Prev =  x%Gamma_NW(iSpan, iNWStart+1, iW) ! Previous circulation
               Gamma_New  =  x%Gamma_NW(iSpan, iNWStart  , iW)
               x%Gamma_NW(iSpan, iNWStart  , iW)  = Gamma_New*ShedScale + (1.0_ReKi-ShedScale) * Gamma_Prev
            enddo
         enddo
      endif
   endif
end subroutine Map_LL_NW

!>  Map the last NW panel with the first FW panel
subroutine Map_NW_FW(p, m, z, x, ErrStat, ErrMsg)
   type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
   type(FVW_MiscVarType),           intent(inout)  :: m              !< Initial misc/optimization variables
   type(FVW_ConstraintStateType),   intent(in   )  :: z              !< Constraints states
   type(FVW_ContinuousStateType),   intent(inout)  :: x              !< Continuous states
   integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   integer(IntKi)            :: iW, iRoot, iTip, iMax
   real(ReKi), dimension(p%nWings) :: FWGamma
   real(ReKi), dimension(p%nSpan+1) :: Gamma_t
   real(ReKi), dimension(p%nSpan) :: sCoord
   real(ReKi) :: FWEpsTip, FWEpsRoot
   real(ReKi) :: ltip, rTip, Gamma_max
   integer(IntKi), parameter :: iAgeFW=1   !< we update the first FW panel
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! First Panel of Farwake has coordinates of last panel of near wake always
   if (p%nFWMax>0) then
      if (m%nNW==p%nNWMax) then
         ! First circulation of Farwake is taken as the max circulation of last NW column
         FWGamma(:)=0.0_ReKi
         do iW=1,p%nWings
            if (p%FullCirculationStart>0 .and. m%nFW<3) then
               ! we might run into the issue that the circulation is 0
               m%iTip(iW) =-1
               m%iRoot(iW)=-1
            endif
            ! NOTE: on the first pass, m%iTip and m%iRoot are computed, TODO per blade
            call PlaceTipRoot(p%nSpan, x%Gamma_NW(:,m%nNW,iW), x%r_NW(1:3,:,m%nNW,iW), x%Eps_NW(1:3,:,m%nNW,iW),& ! inputs
               m%iRoot(iW), m%iTip(iW), FWGamma(iW), FWEpsTip, FWEpsRoot) ! outputs
            x%Gamma_FW(1:FWnSpan,iAgeFW,iW) = FWGamma(iW)
            x%Eps_FW(3,1:FWnSpan,iAgeFW,iW) = FWEpsTip  ! HACK tip put in third
            x%Eps_FW(2,1:FWnSpan,iAgeFW,iW) = FWEpsRoot ! HACK root put in second
            x%Eps_FW(1,1:FWnSpan,iAgeFW,iW) = FWEpsTip  ! For shed vorticity..
         enddo
      endif
      ! Far wake point always mapped to last near wake
      do iW=1,p%nWings
         if (m%nNW==p%nNWMax) then
            iTip  = m%iTip(iW)
            iRoot = m%iRoot(iW)
         else
            iRoot = 1
            iTip  = p%nSpan+1
         endif
         x%r_FW(1:3,1        ,iAgeFW,iW) =  x%r_NW(1:3,iRoot,p%nNWMax+1,iW) ! Point 1 (root)
         x%r_FW(1:3,FWnSpan+1,iAgeFW,iW) =  x%r_NW(1:3,iTip ,p%nNWMax+1,iW) ! Point FWnSpan (tip)
         !if ((FWnSpan==2)) then
         !   ! in between point
         !   x%r_FW(1:3,2,iAgeFW,iW) =  x%r_NW(1:3,int(p%nSpan+1)/4 ,p%nNWMax+1,iW) ! Point (mid)
         !else if ((FWnSpan>2)) then
         !   ErrMsg='Error: FWnSpan>2 not implemented.'
         !   ErrStat=ErrID_Fatal
         !   return
         !endif
      enddo
   endif
   if (.false.) print*,z%Gamma_LL(1,1) ! Just to avoid unused var warning
endsubroutine Map_NW_FW

!> Propagate the positions and circulation one index forward (loop from end to start) 
subroutine PropagateWake(p, m, z, x, ErrStat, ErrMsg)
   type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
   type(FVW_MiscVarType),           intent(inout)  :: m              !< Initial misc/optimization variables
   type(FVW_ConstraintStateType),   intent(in   )  :: z              !< Constraints states
   type(FVW_ContinuousStateType),   intent(inout)  :: x              !< Continuous states
   integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   integer(IntKi) :: iSpan, iAge, iW
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! -- Propagate far wake
      do iW=1,p%nWings
         do iAge=p%nFWMax+1,2,-1 ! 
            do iSpan=1,FWnSpan+1
               x%r_FW(1:3,iSpan,iAge,iW) = x%r_FW(1:3,iSpan,iAge-1,iW)
            enddo
         enddo
         x%r_FW(1:3,1:FWnSpan+1,1,iW) = -999.9_ReKi ! Nullified
      enddo
   if (p%nFWMax>0) then
      do iW=1,p%nWings
         do iAge=p%nFWMax,2,-1
            do iSpan=1,FWnSpan
               x%Gamma_FW(iSpan,iAge,iW) = x%Gamma_FW(iSpan,iAge-1,iW)
               x%Eps_FW(:,iSpan,iAge,iW) = x%Eps_FW(:,iSpan,iAge-1,iW)
            enddo
         enddo
         x%Gamma_FW(1,1:FWnSpan-1,iW) = -999.9_ReKi ! Nullified
         !x%Gamma_FW(:,1,iW) = -999.9_ReKi ! Nullified  ! TODO TODO TODO FIX BUG
      enddo
   endif
   ! --- Propagate near wake
   do iW=1,p%nWings
      do iAge=p%nNWMax+1,iNWStart+1,-1
         do iSpan=1,p%nSpan+1
            x%r_NW(1:3,iSpan,iAge,iW) = x%r_NW(1:3,iSpan,iAge-1,iW)
         enddo
      enddo
      x%r_NW(1:3,:,1:iNWStart,iW) = -999.9_ReKi ! Nullified
   enddo
   if (p%nNWMax>1) then
      do iW=1,p%nWings
         do iAge=p%nNWMax,iNWStart+1,-1
            do iSpan=1,p%nSpan
               x%Gamma_NW(iSpan,iAge,iW) = x%Gamma_NW(iSpan,iAge-1,iW)
               x%Eps_NW(:,iSpan,iAge,iW) = x%Eps_NW(:,iSpan,iAge-1,iW)
            enddo
         enddo
         x%Gamma_NW(:,1:iNWStart,iW) = -999.9_ReKi ! Nullified
      enddo
   endif

   ! Temporary hack for sub-cycling since straight after wkae computation, the wake size will increase
   ! So we do a "fake" propagation here
   do iW=1,p%nWings
      do iAge=p%nFWMax+1,2,-1 ! 
         do iSpan=1,FWnSpan+1
            m%dxdt%r_FW(1:3,iSpan,iAge,iW) = m%dxdt%r_FW(1:3,iSpan,iAge-1,iW)
         enddo
      enddo
      !m%dxdt_FW(1:3,1:FWnSpan+1,1,iW) = -999999_ReKi ! Important not nullified. The best would be to map the last NW convection velocity for this first row.
   enddo
   do iW=1,p%nWings
      do iAge=p%nNWMax+1,iNWStart+1,-1 
         do iSpan=1,p%nSpan+1
            m%dxdt%r_NW(1:3,iSpan,iAge,iW) = m%dxdt%r_NW(1:3,iSpan,iAge-1,iW)
         enddo
      enddo
      m%dxdt%r_NW(1:3,:,1:iNWStart,iW) = 0.0_ReKi ! Nullified, wing do no convect, handled by LL,NW mapping
   enddo

   if (.false.) print*,m%nNW,z%Gamma_LL(1,1) ! Just to avoid unused var warning
end subroutine PropagateWake


!> Print the states, useful for debugging
subroutine print_x_NW_FW(p, m, x, label)
   type(FVW_ParameterType),         intent(in)  :: p              !< Parameters
   type(FVW_MiscVarType),           intent(in)  :: m              !< Initial misc/optimization variables
   type(FVW_ContinuousStateType),   intent(in)  :: x              !< Continuous states
   character(len=*),intent(in) :: label
   integer(IntKi) :: iAge
   character(len=1):: flag
   print*,'------------------------------------------------------------------'
   print'(A,I0,A,I0)',' NW .....................iNWStart:',iNWStart,' nNW:',m%nNW
   do iAge=1,p%nNWMax+1
      flag='X'
      if ((iAge)<= m%nNW+1) flag='.'
      print'(A,A,I0,A)',flag,'iAge ',iAge,'      Root              Tip'
      print*,trim(label)//'x', x%r_NW(1, 1, iAge,1), x%r_NW(1, p%nSpan+1, iAge,1)
      print*,trim(label)//'y', x%r_NW(2, 1, iAge,1), x%r_NW(2, p%nSpan+1, iAge,1)
      print*,trim(label)//'z', x%r_NW(3, 1, iAge,1), x%r_NW(3, p%nSpan+1, iAge,1)
      if (iAge<p%nNWMax+1) then
         print*,trim(label)//'g', x%Gamma_NW(1, iAge,1), x%Gamma_NW(p%nSpan, iAge,1)
         print*,trim(label)//'e', x%Eps_NW(1,1, iAge,1), x%Eps_NW(1,p%nSpan, iAge,1)
      endif
   enddo
   print'(A,I0)','FW <<<<<<<<<<<<<<<<<<<< nFW:',m%nFW
   do iAge=1,p%nFWMax+1
      flag='X'
      if ((iAge)<= m%nFW+1) flag='.'
      print'(A,A,I0,A)',flag,'iAge ',iAge,'      Root              Tip'
      print*,trim(label)//'x', x%r_FW(1, 1, iAge,1), x%r_FW(1, FWnSpan+1, iAge,1)
      print*,trim(label)//'y', x%r_FW(2, 1, iAge,1), x%r_FW(2, FWnSpan+1, iAge,1)
      print*,trim(label)//'z', x%r_FW(3, 1, iAge,1), x%r_FW(3, FWnSpan+1, iAge,1)
      if (iAge<p%nFWMax+1) then
         print*,trim(label)//'g', x%Gamma_FW(1,iAge,1), x%Gamma_FW(FWnSpan, iAge,1)
         print*,trim(label)//'e', x%Eps_FW(1,1, iAge,1), x%Eps_FW(1,FWnSpan, iAge,1)
      endif
   enddo
endsubroutine

!> Debug function to figure out if data have nan
logical function have_nan(p, m, x, u, label)
   type(FVW_ParameterType),         intent(in) :: p !< Parameters
   type(FVW_MiscVarType),           intent(in) :: m !< Initial misc/optimization variables
   type(FVW_ContinuousStateType),   intent(in) :: x !< Continuous states
   type(FVW_InputType),             intent(in) :: u(:) !< Input states
   character(len=*),                intent(in) :: label !< label for print
   have_nan=.False.
   if (any(isnan(x%r_NW))) then
      print*,trim(label),'NaN in r_NW'
      have_nan=.True.
   endif
   if (any(isnan(x%r_FW))) then
      print*,trim(label),'NaN in r_FW'
      have_nan=.True.
   endif
   if (any(isnan(x%Gamma_NW))) then
      print*,trim(label),'NaN in G_NW'
      have_nan=.True.
   endif
   if (any(isnan(x%Gamma_FW))) then
      print*,trim(label),'NaN in G_FW'
      have_nan=.True.
   endif
   if (any(isnan(x%Eps_NW))) then
      print*,trim(label),'NaN in G_FW'
      have_nan=.True.
   endif
   if (any(isnan(x%Eps_FW))) then
      print*,trim(label),'NaN in G_FW'
      have_nan=.True.
   endif
   if (any(isnan(u(1)%V_wind))) then
      print*,trim(label),'NaN in Vwind1'
      have_nan=.True.
   endif
   if (any(isnan(u(2)%V_wind))) then
      print*,trim(label),'NaN in Vwind2'
      have_nan=.True.
   endif
endfunction


! --------------------------------------------------------------------------------
! --- PACKING/UNPACKING FUNCTIONS
! --------------------------------------------------------------------------------
!> Establish the list of points where we will need the free stream
!! The r_wind array is allocated at initialization to the largest size possible.  This is to
!! ensure that we do not violate requirements in the framework later for changing the size
!! of input and output arrays.
subroutine SetRequestedWindPoints(r_wind, x, p, m)
   real(ReKi), dimension(:,:), allocatable,      intent(inout) :: r_wind  !< Position where wind is requested
   type(FVW_ContinuousStateType),   intent(inout)              :: x       !< States
   type(FVW_ParameterType),         intent(in   )              :: p       !< Parameters
   type(FVW_MiscVarType),           intent(in   ), target      :: m       !< Initial misc/optimization variables
   integer(IntKi) :: iP_start,iP_end   ! Current index of point, start and end of range
   integer(IntKi) :: iGrid,i,j,k
   real(ReKi) :: xP,yP,zP,dx,dy,dz
   type(GridOutType), pointer :: g

   ! Using array reshaping to ensure a given near or far wake point is always at the same location in the array.
   ! NOTE: Maximum number of points are passed, whether they "exist" or not. 
   ! NOTE: InflowWind ignores points at (0,0,0)
   !if (DEV_VERSION) then
   !   ! Removing points that don't exist
   !   !call print_x_NW_FW(p,m,x,'wind befr')
   !   if (m%nNW<=p%nNWMax) then
   !      x%r_NW(1:3, 1:p%nSpan+1, m%nNW+2:p%nNWMax+1, 1:p%nWings) = 0.0_ReKi
   !   endif
   !   if ( ((p%nNWMax<=1) .and. (m%nFW==0)) .or. ((m%nFW>0) .and. (m%nFW<=p%nFWMax))) then
   !      x%r_FW(1:3, 1:FWnSpan+1, m%nFW+2:p%nFWMax+1, 1:p%nWings) = 0.0_ReKi
   !   else 
   !      x%r_FW(1:3, 1:FWnSpan+1, m%nFW+1:p%nFWMax+1, 1:p%nWings) = 0.0_ReKi
   !   endif
   !   !call print_x_NW_FW(p,m,x,'wind after')
   !endif

   ! --- LL CP
   iP_start=1
   iP_end=p%nWings*p%nSpan
   r_wind(1:3,iP_start:iP_end) = reshape( m%CP_LL(1:3,1:p%nSpan,1:p%nWings), (/ 3, p%nSpan*p%nWings /))
   ! --- NW points
   iP_start=iP_end+1
   iP_end=iP_start-1+(p%nSpan+1)*(p%nNWMax+1)*p%nWings
   r_wind(1:3,iP_start:iP_end) = reshape( x%r_NW(1:3,1:p%nSpan+1,1:p%nNWMax+1,1:p%nWings), (/ 3, (p%nSpan+1)*(p%nNWMax+1)*p%nWings /))
   ! --- FW points
   if (p%nFWMax>0) then
      iP_start=iP_end+1
      iP_end=iP_start-1+(FWnSpan+1)*(p%nFWMax+1)*p%nWings
      r_wind(1:3,iP_start:iP_end) = reshape( x%r_FW(1:3,1:FWnSpan+1,1:p%nFWMax+1,1:p%nWings), (/ 3, (FWnSpan+1)*(p%nFWMax+1)*p%nWings /))
   endif
   ! --- VTK points
   ! TODO optimize this, and do it only once
   iP_start=iP_end+1
   do iGrid=1,p%nGridOut
      g => m%GridOutputs(iGrid)
      dx = (g%xEnd- g%xStart)/max(g%nx-1,1)
      dy = (g%yEnd- g%yStart)/max(g%ny-1,1)
      dz = (g%zEnd- g%zStart)/max(g%nz-1,1)
      do k=1,g%nz
         zP = g%zStart  + (k-1)*dz
         do j=1,g%ny
            yP = g%yStart  + (j-1)*dy
            do i=1,g%nx
               xP = g%xStart  + (i-1)*dx
               r_wind(1:3,iP_start) = (/xP,yP,zP/)
               iP_start=iP_start+1
            enddo
         enddo
      enddo ! Loop on z
   enddo ! Loop on grids

   !if (DEV_VERSION) then
   !   ! Additional checks
   !   if (any(r_wind(3,:)<=-99999_ReKi)) then
   !      call print_x_NW_FW(p,m,x,'wind after')
   !      print*,'Error in wind'
   !      STOP
   !   endif
   !   ! Removing points that don't exist
   !   if (m%nNW<=p%nNWMax) then
   !      x%r_NW(1:3, 1:p%nSpan+1, m%nNW+2:p%nNWMax+1, 1:p%nWings) = -999999.0_ReKi
   !   endif
   !   if ( ((p%nNWMax<=1) .and. (m%nFW==0)) .or. ((m%nFW>0) .and. (m%nFW<=p%nFWMax))) then
   !      x%r_FW(1:3, 1:FWnSpan+1, m%nFW+2:p%nFWMax+1, 1:p%nWings) =-999999.0_ReKi
   !   else 
   !      x%r_FW(1:3, 1:FWnSpan+1, m%nFW+1:p%nFWMax+1, 1:p%nWings) =-999999.0_ReKi
   !   endif
   !endif

end subroutine SetRequestedWindPoints


!> Set the requested wind into the correponding misc variables
subroutine DistributeRequestedWind_LL(V_wind, p, Vwnd_LL)
   real(ReKi), dimension(:,:),      intent(in   ) :: V_wind  !< Requested wind, packed
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   real(ReKi), dimension(:,:,:),    intent(inout) :: Vwnd_LL !< Wind on lifting line
   integer(IntKi)          :: iP_start,iP_end   ! Current index of point, start and end of range
   ! Using array reshaping to ensure a given near or far wake point is always at the same location in the array.
   ! NOTE: Maximum number of points are passed, whether they "exist" or not. 
   ! --- LL CP
   iP_start=1
   iP_end=p%nWings*p%nSpan
   Vwnd_LL(1:3,1:p%nSpan,1:p%nWings) = reshape( V_wind(1:3,iP_start:iP_end), (/ 3, p%nSpan, p%nWings /))
end subroutine DistributeRequestedWind_LL

subroutine DistributeRequestedWind_NWFW(V_wind, p, Vwnd_NW, Vwnd_FW)
   real(ReKi), dimension(:,:),      intent(in   ) :: V_wind  !< Requested wind, packed
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   real(ReKi), dimension(:,:,:,:),  intent(inout) :: Vwnd_NW !< Wind on near wake panels
   real(ReKi), dimension(:,:,:,:),  intent(inout) :: Vwnd_FW !< Wind on near wake panels
   integer(IntKi)          :: iP_start,iP_end   ! Current index of point, start and end of range
   ! --- NW points
   iP_start=p%nWings*p%nSpan+1
   iP_end=iP_start-1+(p%nSpan+1)*(p%nNWMax+1)*p%nWings
   Vwnd_NW(1:3,1:p%nSpan+1,1:p%nNWMax+1,1:p%nWings) = reshape( V_wind(1:3,iP_start:iP_end), (/ 3, p%nSpan+1, p%nNWMax+1, p%nWings/))
   ! --- FW points
   if (p%nFWMax>0) then
      iP_start=iP_end+1
      iP_end=iP_start-1+(FWnSpan+1)*(p%nFWMax+1)*p%nWings
      Vwnd_FW(1:3,1:FWnSpan+1,1:p%nFWMax+1,1:p%nWings) = reshape( V_wind(1:3,iP_start:iP_end), (/ 3, FWnSpan+1, p%nFWMax+1, p%nWings /))
   endif
end subroutine DistributeRequestedWind_NWFW

!> Set the requested wind into the correponding misc variables
subroutine DistributeRequestedWind_Grid(V_wind, p, m)
   real(ReKi), dimension(:,:),      intent(in   ) :: V_wind  !< Requested wind, packed
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_MiscVarType), target,   intent(inout) :: m       !< Initial misc/optimization variables
   integer(IntKi)          :: iP_start,iP_end   ! Current index of point, start and end of range
   integer(IntKi) :: iGrid,i,j,k
   type(GridOutType), pointer :: g
   ! --- LL CP
   iP_end  =p%nWings*p%nSpan+1-1+(p%nSpan+1)*(p%nNWMax+1)*p%nWings
   ! --- FW points
   if (p%nFWMax>0) then
      iP_end=iP_end+1-1+(FWnSpan+1)*(p%nFWMax+1)*p%nWings
   endif
   ! --- VTK points
   ! TODO optimize this
   iP_start=iP_end+1
   do iGrid=1,p%nGridOut
      g => m%GridOutputs(iGrid)
      do k=1,g%nz
         do j=1,g%ny
            do i=1,g%nx
               g%uGrid(1:3,i,j,k) = V_wind(1:3,iP_start)
               iP_start=iP_start+1
            enddo
         enddo
      enddo ! Loop on x
   enddo ! Loop on grids
end subroutine DistributeRequestedWind_Grid



!> Count how many segments are needed to represent the Near wake and far wakes, starting at a given depth
subroutine CountSegments(p, nNW, nFW, iDepthStart, nSeg, nSegP, nSegNW)
   type(FVW_ParameterType), intent(in   ) :: p    !< Parameters
   integer(IntKi),          intent(in   ) :: nNW  !< Number of NW panels
   integer(IntKi),          intent(in   ) :: nFW  !< Number of FW panels
   integer(IntKi),          intent(in   ) :: iDepthStart !< Index where we start packing for NW panels
   integer(IntKi),          intent(  out) :: nSeg   !< Total number of segments after packing
   integer(IntKi),          intent(  out) :: nSegP  !< Total number of segments points after packing
   integer(IntKi),          intent(  out) :: nSegNW !< Total number of segments points for the near wake only
   logical        :: LastNWShed
   ! If the FW contains Shed vorticity, we include the last shed vorticity from the NW, otherwise, we don't!
   ! It's important not to include it, otherwise a strong vortex will be present there with no compensating vorticity from the FW
   LastNWShed = (p%FWShedVorticity ) .or. ((.not.p%FWShedVorticity) .and. (nNW<p%nNWMax))
   ! --- Counting total number of segments
   nSegP=0; nSeg=0; nSegNW=0
   ! NW segments
   if ((nNW-iDepthStart)>=0) then
      nSegP  =      p%nWings * (  (p%nSpan+1)*(nNW-iDepthStart+2)            )
      nSegNW =      p%nWings * (2*(p%nSpan+1)*(nNW-iDepthStart+2)-(p%nSpan+1)-(nNW-iDepthStart+1+1))  
      if (.not.LastNWShed) then
         nSegNW =   nSegNW - p%nWings * (p%nSpan) ! Removing last set of shed segments
      endif
   endif
   nSeg=nSegNW
   ! FW segments
   if (nFW>0) then
      nSegP  = nSegP + p%nWings * (  (FWnSpan+1)*(nFW+1) )
      if (p%FWShedVorticity) then
         nSeg = nSeg + p%nWings * (2*(FWnSpan+1)*(nFW+1)-(FWnSpan+1)-(nFW+1))  
      else
         nSeg = nSeg + p%nWings * (  (FWnSpan+1)*(nFW)                    )   ! No Shed vorticity
      endif
   endif
end subroutine CountSegments

!> Count how many control points are convecting (needed to compute the wake convection)
pure integer(IntKi) function CountCPs(p, nNW, nFWEff) result(nCPs)
   type(FVW_ParameterType), intent(in   ) :: p       !< Parameters
   integer(IntKi),          intent(in   ) :: nNW     !< Number of NW panels
   integer(IntKi),          intent(in   ) :: nFWEff  !< Number of effective (ie. convecting) FW panels
   nCPs =  p%nWings * (  (p%nSpan+1)*(nNW+1) )
   if (nFWEff>0)  nCPs = nCPs + p%nWings * ((FWnSpan+1)*(nFWEff+1) )
end function CountCPs


subroutine PackPanelsToSegments(p, x, iDepthStart, bMirror, nNW, nFW, SegConnct, SegPoints, SegGamma, SegEpsilon, nSeg, nSegP)
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_ContinuousStateType),   intent(in   ) :: x       !< States
   integer(IntKi),                  intent(in   ) :: iDepthStart !< Index where we start packing for NW panels
   logical,                         intent(in   ) :: bMirror !< Mirror the vorticity wrt the ground
   integer(IntKi),                  intent(in   ) :: nNW, NFW !< Number of near/far wake panels
   integer(IntKi),dimension(:,:), intent(inout) :: SegConnct !< Segment connectivity
   real(ReKi),    dimension(:,:), intent(inout) :: SegPoints !< Segment Points
   real(ReKi),    dimension(:)  , intent(inout) :: SegGamma  !< Segment Circulation
   real(ReKi),    dimension(:)  , intent(inout) :: SegEpsilon  !< Segment Circulation
   integer(IntKi), intent(out)                :: nSeg      !< Total number of segments after packing
   integer(IntKi), intent(out)                :: nSegP     !< Total number of segments points after packing
   ! Local
   integer(IntKi) :: iHeadC, iHeadP, nC, nCNW, nP, iW, iHeadC_bkp, i, iMirror
   logical        :: LastNWShed

   ! If the FW contains Shed vorticity, we include the last shed vorticity form the NW, orhtwerise, we don't!
   ! It's important not to include it, otherwise a strong vortex will be present there with no compensating vorticity from the FW
   LastNWShed = (p%FWShedVorticity ) .or. ((.not.p%FWShedVorticity) .and. (nNW<p%nNWMax))

   ! Counting total number of segments
   ! Returns nC, nP, nCNW, number of segments (without accounting for mirroring)
   call CountSegments(p, nNW, nFW, iDepthStart, nC, nP, nCNW)

   if (nP>0) then
      ! Nullifying for safety
      SegConnct=-1
      SegPoints=-1
      SegGamma =-1
      !
      iHeadP=1
      iHeadC=1
      if (nCNW>0) then
         do iW=1,p%nWings
            call LatticeToSegments(x%r_NW(1:3,:,1:nNW+1,iW), x%Gamma_NW(:,1:nNW,iW), x%Eps_NW(1:3,:,1:nNW,iW), iDepthStart, SegPoints, SegConnct, SegGamma, SegEpsilon, iHeadP, iHeadC, .True., LastNWShed, .false.)
         enddo
      endif
      if (nFW>0) then
         iHeadC_bkp = iHeadC
         do iW=1,p%nWings
            call LatticeToSegments(x%r_FW(1:3,:,1:nFW+1,iW), x%Gamma_FW(:,1:nFW,iW), x%Eps_FW(1:3,:,1:nFW,iW), 1, SegPoints, SegConnct, SegGamma, SegEpsilon, iHeadP, iHeadC , p%FWShedVorticity, p%FWShedVorticity, .true.)
         enddo
         SegConnct(3,iHeadC_bkp:) = SegConnct(3,iHeadC_bkp:) + nNW ! Increasing iDepth (or age) to account for NW
      endif
      if (DEV_VERSION) then
         ! Safety checks
         if ((iHeadP-1)/=nP) then
            print*,'PackPanelsToSegments: Number of points wrongly estimated',nP, iHeadP-1
            STOP ! Keep me. The check will be removed once the code is well established
         endif
         if ((iHeadC-1)/=nC) then
            print*,'PackPanelsToSegments: Number of segments wrongly estimated',nC, iHeadC-1
            STOP ! Keep me. The check will be removed once the code is well established
         endif
         if (any(SegPoints(3,:)<-99._ReKi)) then
            print*,'PackPanelsToSegments: some segments are NAN'
            STOP ! Keep me. The check will be removed once the code is well established
         endif
      endif
      nSeg  = iHeadC-1
      nSegP = iHeadP-1

      if (bMirror) then
         ! Mirroring the segments directly
         ! NOTE: an alternative is to handle this in the Biot-Savart law directly...
         do i=1,nSeg
            iMirror = i + nSeg
            SegConnct(1:2, iMirror) =  SegConnct(1:2, i) + nSegP ! Increased point indices
            SegConnct(3:4, iMirror) =  SegConnct(3:4, i) ! Span and age is copied
            SegGamma(iMirror)       = -SegGamma(i)       ! Vorticity needs mirroring
         enddo
         do i=1,nSegP
            iMirror = i + nSegP
            SegPoints(1:2, iMirror) =   SegPoints(1:2, i) ! Same x and y
            SegPoints(3  , iMirror) = - SegPoints(3  , i) ! Mirror with respect to z=0
         enddo
         ! We now have double the amount of segments and points
         nSeg  = nSeg*2
         nSegP = nSegP*2
      endif
   else
      nSeg  = 0
      nSegP = 0
   endif
end subroutine PackPanelsToSegments

!> Set up regularization parameter based on diffusion method and regularization method
!! NOTE: - reg param is now stored at panel level
!!       - continuous variables are used, only the LL and NW panel needs to be set at t=0
subroutine FVW_InitRegularization(x, p, m, ErrStat, ErrMsg)
   type(FVW_ContinuousStateType),   intent(inout) :: x       !< States
   type(FVW_ParameterType),         intent(inout) :: p       !< Parameters
   type(FVW_MiscVarType),           intent(inout) :: m       !< Initial misc/optimization variables
   integer(IntKi),                  intent(  out) :: ErrStat    !< Error status of the operation
   character(*),                    intent(  out) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   ! Local variables
   real(ReKi) :: ds_min, ds_max, ds_mean, ds !< min,max and mean of spanwise sections
   real(ReKi) :: c_min, c_max, c_mean !< min,max and mean of chord
   real(ReKi) :: d_min, d_max, d_mean !< min,max and mean of panel diagonal
   real(ReKi) :: RegParam
   real(ReKi) :: Span !< "Blade span"
   integer :: iW, iSpan
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! --- Compute min max and mean spanwise section lengths
   iW =1
   ds_min  = minval(p%s_LL(2:p%nSpan+1,iW)-p%s_LL(1:p%nSpan,iW))
   ds_max  = maxval(p%s_LL(2:p%nSpan+1,iW)-p%s_LL(1:p%nSpan,iW))
   ds_mean = sum(p%s_LL(2:p%nSpan+1,iW)-p%s_LL(1:p%nSpan,iW))/(p%nSpan+1)
   c_min  = minval(p%chord_LL(:,iW))
   c_max  = maxval(p%chord_LL(:,iW))
   c_mean = sum   (p%chord_LL(:,iW))/(p%nSpan+1)
   d_min  = minval(m%diag_LL(:,iW))
   d_max  = maxval(m%diag_LL(:,iW))
   d_mean = sum   (m%diag_LL(:,iW))/(p%nSpan+1)
   Span    = p%s_LL(p%nSpan+1,iW)-p%s_LL(1,iW)
   RegParam = ds_mean*2

   ! Default init of reg param
   x%Eps_NW(1:3,:,:,:) = 0.001_ReKi
   x%Eps_FW(1:3,:,:,:) = 0.001_ReKi
   if (DEV_VERSION) then
      write(*,'(A)')'-----------------------------------------------------------------------------------------'
      write(*,'(A)')'Regularization Info'
      write(*,'(A,1F8.4,A)') 'Span                   : ',Span
      write(*,'(A,3F8.4,A)') 'Chord                  : ',c_min,c_mean,c_max,' (min, mean, max)'
      write(*,'(A,3F8.4,A)') 'Spanwise discretization: ',ds_min,ds_mean,ds_max,' (min, mean, max)'
      write(*,'(A,3F8.4,A)') 'Diagonal discretization: ',d_min,d_mean,d_max,' (min, mean, max)'
      write(*,'(A,1F8.4)')   'RegParam (Recommended) : ',RegParam
      write(*,'(A,1F8.4)')   'RegParam (Input      ) : ',p%WakeRegParam
   endif

   if (p%RegDeterMethod==idRegDeterConstant) then
      ! Constant reg param throughout the wake
      if (p%WakeRegMethod==idRegAge) then ! NOTE: age method implies a division by rc
         p%WingRegParam=max(0.01, p%WingRegParam)
         p%WakeRegParam=max(0.01, p%WakeRegParam)
      endif

      ! Set reg param on wing and first NW
      ! NOTE: setting the same in all three directions for now, TODO!
      x%Eps_NW(1:3,:,1,:) = p%WingRegParam ! First age is always WingRegParam (LL)
      if (p%nNWMax>1) then
         x%Eps_NW(1:3,:,2,:) = p%WakeRegParam ! Second age is always WakeRegParam
      endif

   else if (p%RegDeterMethod==idRegDeterAuto) then
      ! TODO this is beta
      print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print*,'!!! NOTE: using optimized wake regularization parameters is still a beta feature!'
      print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      p%WakeRegMethod      = idRegAge
      p%RegFunction        = idRegVatistas
      p%WakeRegParam       = RegParam
      p%WingRegParam       = RegParam
      p%CoreSpreadEddyVisc = 1000
      m%Sgmt%RegFunction    = p%RegFunction
      write(*,'(A)'   )   'The following regularization parameters will be used:'
      write(*,'(A,I0)'   )   'WakeRegMethod     : ', p%WakeRegMethod
      write(*,'(A,I0)'   )   'RegFunction       : ', p%RegFunction
      write(*,'(A,1F8.4)')   'WakeRegParam      : ', p%WakeRegParam
      write(*,'(A,1F8.4)')   'BladeRegParam     : ', p%WingRegParam
      write(*,'(A,1F9.4)')   'CoreSpreadEddyVisc: ', p%CoreSpreadEddyVisc
   ! Set reg param on wing and first NW
   ! NOTE: setting the same in all three directions for now, TODO!
   x%Eps_NW(1:3,:,1,:) = p%WingRegParam ! First age is always WingRegParam (LL)
   if (p%nNWMax>1) then
      x%Eps_NW(1:3,:,2,:) = p%WakeRegParam ! Second age is always WakeRegParam
   endif

   else if (p%RegDeterMethod==idRegDeterChord) then
      ! Using chord to scale the reg param
      do iW=1,p%nWings
         do iSpan=1,p%nSpan
            x%Eps_NW(1:3, iSpan, 1, iW) = p%WingRegParam * p%chord_CP_LL(iSpan, iW)
            if (p%nNWMax>1) then
               x%Eps_NW(1:3, iSpan, 2, iW) = p%WakeRegParam * p%chord_CP_LL(iSpan, iW)
            endif
         enddo
      enddo

   else if (p%RegDeterMethod==idRegDeterSpan) then
      ! Using dr to scale the reg param
      do iW=1,p%nWings
         do iSpan=1,p%nSpan
            ds = p%s_LL(iSpan+1,iW)-p%s_LL(iSpan,iW)
            x%Eps_NW(1:3, iSpan, 1, iW) = p%WingRegParam * ds
            if (p%nNWMax>1) then
               x%Eps_NW(1:3, iSpan, 2, iW) = p%WakeRegParam * ds
            endif
         enddo
      enddo
   else ! Should never happen (caught earlier)
      ErrStat = ErrID_Fatal
      ErrMsg ='Regularization determination method not implemented' 
   endif

   call WrScr(' - Regularization parameters:')
   write(*,'(A,2F8.4)') '    BladeReg (min/max): ', minval(x%Eps_NW(:, :, 1, :)), maxval(x%Eps_NW(:, :, 1, :))
   if (p%nNWMax>1) then
      write(*,'(A,2F8.4)')    '    WakeReg (min/max) : ', minval(x%Eps_NW(:,:, 2, :)), maxval(x%Eps_NW(:,:, 2, :))
   endif
   write(*,'(A,2F8.4)') '    k = alpha delta nu: ', CoreSpreadAlpha * p%CoreSpreadEddyVisc * p%KinVisc

end subroutine FVW_InitRegularization


!> Compute induced velocities from all vortex elements onto nPoints
!! In : x, x%r_NW, x%r_FW, x%Gamma_NW, x%Gamma_FW
!! Out: Vind
subroutine InducedVelocitiesAll_OnGrid(g, p, x, m, ErrStat, ErrMsg)
   type(GridOutType),               intent(inout) :: g       !< Grid on whcih to compute the velocity
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_ContinuousStateType),   intent(in   ) :: x       !< States
   type(FVW_MiscVarType),           intent(inout) :: m       !< Initial misc/optimization variables
   integer(IntKi),                  intent(  out) :: ErrStat !< Error status of the operation
   character(*),                    intent(  out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi) :: nCPs, iHeadP
   integer(IntKi) :: i,j,k
   real(ReKi) :: xP,yP,zP,dx,dy,dz
   ! TODO new options
   type(T_Tree)   :: Tree
   type(T_Part)   :: Part
   real(ReKi), dimension(:,:), allocatable :: CPs  ! TODO get rid of me with dedicated functions
   real(ReKi), dimension(:,:), allocatable :: Uind ! TODO get rid of me with dedicated functions
   ErrStat= ErrID_None
   ErrMsg =''

   ! --- Packing control points
   nCPs = g%nx * g%ny * g%nz
   allocate(CPs(3, nCPs))
   iHeadP=1
   dx = (g%xEnd- g%xStart)/max(g%nx-1,1)
   dy = (g%yEnd- g%yStart)/max(g%ny-1,1)
   dz = (g%zEnd- g%zStart)/max(g%nz-1,1)
   do k=1,g%nz
      zP = g%zStart  + (k-1)*dz
      do j=1,g%ny
         yP = g%yStart  + (j-1)*dy
         do i=1,g%nx
            xP = g%xStart  + (i-1)*dx
            CPs(1:3,iHeadP) = (/xP,yP,zP/)
            iHeadP=iHeadP+1
         enddo
      enddo
   enddo ! Loop on z

   ! --- Packing Uind points
   allocate(Uind(3, nCPs)); Uind=0.0_ReKi
   iHeadP=1
   call FlattenValues(g%uGrid, Uind, iHeadP); ! NOTE: Uind contains uGrid now (Uwnd)

   ! --- Compute induced velocity
   ! Convert Panels to segments, segments to particles, particles to tree
   call InducedVelocitiesAll_Init(p, x, m, m%Sgmt, Part, Tree, ErrStat, ErrMsg)
   call InducedVelocitiesAll_Calc(CPs, nCPs, Uind, p, m%Sgmt, Part, Tree, ErrStat, ErrMsg)
   call InducedVelocitiesAll_End(p, m, Tree, Part, ErrStat, ErrMsg)

   ! --- Unpacking induced velocity points
   iHeadP=1
   call DeflateValues(Uind, g%uGrid, iHeadP)

   deallocate(CPs)
   deallocate(Uind)

end subroutine InducedVelocitiesAll_OnGrid



!> Perform initialization steps before requesting induced velocities from All vortex elements
!! In : x%r_NW, x%r_FW, x%Gamma_NW, x%Gamma_FW
!! Out: Tree, Part, m
subroutine InducedVelocitiesAll_Init(p, x, m, Sgmt, Part, Tree,  ErrStat, ErrMsg)
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_ContinuousStateType),   intent(in   ) :: x       !< States
   type(FVW_MiscVarType),           intent(in   ) :: m       !< Misc
   type(T_Sgmt),                    intent(inout) :: Sgmt    !< Segments
   type(T_Part),                    intent(out)   :: Part    !< Particle storage if needed
   type(T_Tree),                    intent(out)   :: Tree    !< Tree of particles if needed
   integer(IntKi),                  intent(  out) :: ErrStat !< Error status of the operation
   character(*),                    intent(  out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi) :: iHeadP, nSeg, nSegP
   logical        :: bMirror ! True if we mirror the vorticity wrt ground
   integer(IntKi) :: nPart
   ErrStat= ErrID_None
   ErrMsg =''

   bMirror = p%ShearModel==idShearMirror ! Whether or not we mirror the vorticity wrt ground

   ! --- Packing all vortex elements into a list of segments
   call PackPanelsToSegments(p, x, 1, bMirror, m%nNW, m%nFW, Sgmt%Connct, Sgmt%Points, Sgmt%Gamma, Sgmt%Epsilon, nSeg, nSegP)
   Sgmt%RegFunction=p%RegFunction
   Sgmt%nAct  = nSeg
   Sgmt%nActP = nSegP

   ! --- Converting to particles
   if ((p%VelocityMethod==idVelocityTree) .or. (p%VelocityMethod==idVelocityPart)) then
      iHeadP=1
      nPart = p%PartPerSegment * nSeg 
      allocate(Part%P(3,nPart), Part%Alpha(3,nPart), Part%RegParam(nPart))
      Part%Alpha(:,:)  = -99999.99_ReKi
      Part%P(:,:)      = -99999.99_ReKi
      Part%RegParam(:) = -99999.99_ReKi
      call SegmentsToPart(Sgmt%Points, Sgmt%Connct, Sgmt%Gamma, Sgmt%Epsilon, 1, nSeg, p%PartPerSegment, Part%P, Part%Alpha, Part%RegParam, iHeadP)
      if (p%RegFunction/=idRegNone) then
         Part%RegFunction = idRegExp ! TODO need to find a good equivalence and potentially adapt Epsilon in SegmentsToPart
      endif
      if (DEV_VERSION) then
         if (any(Part%RegParam(:)<-9999.99_ReKi)) then
            print*,'Error in Segment to part conversion'
            STOP
         endif
      endif
   endif

   ! Grow tree if needed
   if (p%VelocityMethod==idVelocityTree) then
      Tree%DistanceDirect = 2*sum(Part%RegParam)/size(Part%RegParam) ! 2*mean(eps), below that distance eps has a strong effect
      call grow_tree(Tree, Part%P, Part%Alpha, Part%RegFunction, Part%RegParam, 0)
   endif

end subroutine InducedVelocitiesAll_Init

!> Compute induced velocity on flat CPs 
subroutine InducedVelocitiesAll_Calc(CPs, nCPs, Uind, p, Sgmt, Part, Tree, ErrStat, ErrMsg)
   real(ReKi), dimension(:,:),      intent(in)    :: CPs     !< Control points (3 x nCPs++)
   integer(IntKi)                 , intent(in)    :: nCPs    !< Number of control points on which to compute (nCPs <= size(CPs,2))
   real(ReKi), dimension(:,: )    , intent(inout) :: Uind    !< Induced velocity vector - Side effects!!! (3 x nCPs++)
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(T_Sgmt),                    intent(in   ) :: Sgmt    !< Tree of particles if needed
   type(T_Part),                    intent(in   ) :: Part    !< Particle storage if needed
   type(T_Tree),                    intent(inout) :: Tree    !< Tree of particles if needed
   integer(IntKi),                  intent(  out) :: ErrStat !< Error status of the operation
   character(*),                    intent(  out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
   ! Local variables
   ErrStat= ErrID_None
   ErrMsg =''

   if (p%VelocityMethod==idVelocityBasic) then
      call ui_seg( 1, nCPs, CPs, 1, Sgmt%nAct, Sgmt%nAct, Sgmt%nActP, Sgmt%Points, Sgmt%Connct, Sgmt%Gamma, Sgmt%RegFunction, Sgmt%Epsilon, Uind)

   elseif (p%VelocityMethod==idVelocityTree) then
      ! Tree has already been grown with InducedVelocitiesAll_Init
      !call print_tree(Tree)
      call ui_tree(Tree, CPs, 0, 1, nCPs, p%TreeBranchFactor, Tree%DistanceDirect, Uind, ErrStat, ErrMsg)

   elseif (p%VelocityMethod==idVelocityPart) then
      call ui_part_nograd(CPs ,Part%P, Part%Alpha, Part%RegFunction, Part%RegParam, Uind, nCPs, size(Part%P,2))
   endif
end subroutine InducedVelocitiesAll_Calc


!> Perform termination steps after velocity was requested from all vortex elements
!! InOut: Tree, Part, m
subroutine InducedVelocitiesAll_End(p, m, Tree, Part, ErrStat, ErrMsg)
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_MiscVarType),           intent(inout) :: m       !< Initial misc/optimization variables
   type(T_Tree),                    intent(inout) :: Tree    !< Tree of particles if needed
   type(T_Part),                    intent(inout) :: Part    !< Particle storage if needed
   integer(IntKi),                  intent(  out) :: ErrStat !< Error status of the operation
   character(*),                    intent(  out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
   ! Local variables
   ErrStat= ErrID_None
   ErrMsg =''

   if (p%VelocityMethod==idVelocityBasic) then
      ! Nothing

   elseif (p%VelocityMethod==idVelocityTree) then
      call cut_tree(Tree)
      deallocate(Part%P, Part%Alpha, Part%RegParam)

   elseif (p%VelocityMethod==idVelocityPart) then
      deallocate(Part%P, Part%Alpha, Part%RegParam)
   endif

end subroutine InducedVelocitiesAll_End




!> Compute induced velocities from all vortex elements onto all the vortex elements
!! In : x%r_NW, x%r_FW, x%Gamma_NW, x%Gamma_FW
!! Out: m%Vind_NW, m%Vind_FW
subroutine WakeInducedVelocities(p, x, m, ErrStat, ErrMsg)
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_ContinuousStateType),   intent(in   ) :: x       !< States
   type(FVW_MiscVarType),           intent(inout) :: m       !< Initial misc/optimization variables
   integer(IntKi),                  intent(  out) :: ErrStat !< Error status of the operation
   character(*),                    intent(  out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi) :: iW, nCPs, iHeadP
   integer(IntKi) :: nFWEff  ! Number of farwake panels that are free at current tmie step
   type(T_Tree)   :: Tree
   type(T_Part)   :: Part
   ErrStat= ErrID_None
   ErrMsg =''

   nFWEff = min(m%nFW, p%nFWFree)

   ! --- Pack control points
   call PackConvectingPoints() ! m%CPs

   ! --- Compute induced velocity
   ! Convert Panels to segments, segments to particles, particles to tree
   m%Uind=0.0_ReKi ! very important due to side effects of ui_* methods
   call InducedVelocitiesAll_Init(p, x, m, m%Sgmt, Part, Tree, ErrStat, ErrMsg)
   call InducedVelocitiesAll_Calc(m%CPs, nCPs, m%Uind, p, m%Sgmt, Part, Tree, ErrStat, ErrMsg)
   call InducedVelocitiesAll_End(p, m, Tree, Part, ErrStat, ErrMsg)
   call UnPackInducedVelocity()

   if (DEV_VERSION) then
      print'(A,I0,A,I0,A,I0)','Convection - nSeg:',m%Sgmt%nAct,' - nSegP:',m%Sgmt%nActP, ' - nCPs:',nCPs
   endif
contains
   !> Pack all the points that convect 
   subroutine PackConvectingPoints()
      ! Counting total number of control points that convects
      nCPs = CountCPs(p, m%nNW, nFWEff)
      m%CPs=-999.9_ReKi
      ! Packing
      iHeadP=1
      do iW=1,p%nWings
         CALL LatticeToPoints(x%r_NW(1:3,:,1:m%nNW+1,iW), 1, m%CPs, iHeadP)
      enddo
      if (nFWEff>0) then
         do iW=1,p%nWings
            CALL LatticeToPoints(x%r_FW(1:3,:,1:nFWEff+1,iW), 1, m%CPs, iHeadP)
         enddo
      endif
      if (DEV_VERSION) then
         ! Additional checks
         if (any(m%CPs(1,1:nCPs)<=-99)) then
            call print_x_NW_FW(p,m,x,'pack')
            ErrMsg='PackConvectingPoints: Problem in Control points'; ErrStat=ErrID_Fatal; return
         endif
         if ((iHeadP-1)/=nCPs) then
            print*,'PackConvectingPoints: Number of points wrongly estimated',nCPs, iHeadP-1
            STOP ! Keep me. The check will be removed once the code is well established
            ErrMsg='PackConvectingPoints: Number of points wrongly estimated '; ErrStat=ErrID_Fatal; return
         endif
      endif
   end subroutine
   !> Distribute the induced velocity to the proper location 
   subroutine UnPackInducedVelocity()
      m%Vind_NW = -9999._ReKi !< Safety
      m%Vind_FW = -9999._ReKi !< Safety
      iHeadP=1
      do iW=1,p%nWings
         CALL VecToLattice(m%Uind, 1, m%Vind_NW(:,:,1:m%nNW+1,iW), iHeadP)
      enddo
      if (nFWEff>0) then 
         do iW=1,p%nWings
            CALL VecToLattice(m%Uind, 1, m%Vind_FW(1:3,1:FWnSpan+1,1:nFWEff+1,iW), iHeadP)
         enddo
         if (DEV_VERSION) then
            if (any(m%Vind_FW(1:3,1:FWnSpan+1,1:nFWEff+1,:)<-99)) then
               ErrMsg='UnPackInducedVelocity: Problem in FW induced velocity on FW points'; ErrStat=ErrID_Fatal; return
            endif
         endif
      endif
      if (DEV_VERSION) then
         if ((iHeadP-1)/=nCPs) then
            print*,'UnPackInducedVelocity: Number of points wrongly estimated',nCPs, iHeadP-1
            STOP ! Keep me. The check will be removed once the code is well established
            ErrMsg='UnPackInducedVelocity: Number of points wrongly estimated'; ErrStat=ErrID_Fatal; return
         endif
      endif
   end subroutine

end subroutine WakeInducedVelocities

!> Compute induced velocities from all vortex elements onto the lifting line control points
!! In : x%r_NW, x%r_FW, x%Gamma_NW, x%Gamma_FW
!! Out: m%Vind_LL
subroutine LiftingLineInducedVelocities(CP_LL, p, x, iDepthStart, m, Vind_LL, ErrStat, ErrMsg)
   real(ReKi), dimension(:,:,:),    intent(in   ) :: CP_LL   !< Control points where velocity is to be evaluated
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_ContinuousStateType),   intent(in   ) :: x       !< States
   integer(IntKi),                  intent(in   ) :: iDepthStart !< Index where we start packing for NW panels
   type(FVW_MiscVarType),           intent(inout) :: m       !< Initial misc/optimization variables
   real(ReKi), dimension(:,:,:),    intent(  out) :: Vind_LL !< Control points where velocity is to be evaluated
   ! Local variables
   integer(IntKi) :: iW, nSeg, nSegP, nCPs, iHeadP
   real(ReKi),    dimension(:,:), allocatable :: CPs   !< ControlPoints
   real(ReKi),    dimension(:,:), allocatable :: Uind  !< Induced velocity
   integer(IntKi),              intent(  out) :: ErrStat    !< Error status of the operation
   character(*),                intent(  out) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   logical ::  bMirror 
   ErrStat = ErrID_None
   ErrMsg  = ""
   Vind_LL = -9999._ReKi !< Safety
   bMirror = p%ShearModel==idShearMirror ! Whether or not we mirror the vorticity wrt ground

   ! --- Packing all vortex elements into a list of segments
   call PackPanelsToSegments(p, x, iDepthStart, bMirror, m%nNW, m%nFW, m%Sgmt%Connct, m%Sgmt%Points, m%Sgmt%Gamma, m%Sgmt%Epsilon, nSeg, nSegP)

   ! --- Computing induced velocity
   if (nSegP==0) then
      nCPs=0
      Vind_LL = 0.0_ReKi
      if (DEV_VERSION) then
         print'(A,I0,A,I0,A,I0,A)','Induction -  nSeg:',nSeg,' - nSegP:',nSegP, ' - nCPs:',nCPs, ' -> No induction'
      endif
   else
      nCPs=p%nWings * p%nSpan
      allocate(CPs (1:3,1:nCPs)) ! NOTE: here we do allocate CPs and Uind insteadof using Misc 
      allocate(Uind(1:3,1:nCPs)) !       The size is reasonably small, and m%Uind then stay filled with "rollup velocities" (for export)
      Uind=0.0_ReKi !< important due to side effects of ui_seg
      ! ---
      call PackLiftingLinePoints()
      if (DEV_VERSION) then
         print'(A,I0,A,I0,A,I0)','Induction -  nSeg:',nSeg,' - nSegP:',nSegP, ' - nCPs:',nCPs
      endif
      call ui_seg( 1, nCPs, CPs, 1, nSeg, nSeg, nSegP, m%Sgmt%Points, m%Sgmt%Connct, m%Sgmt%Gamma, m%Sgmt%RegFunction, m%Sgmt%Epsilon, Uind)
      call UnPackLiftingLineVelocities()

      deallocate(Uind)
      deallocate(CPs)
   endif
contains
   !> Pack all the control points
   subroutine PackLiftingLinePoints()
      iHeadP=1
      do iW=1,p%nWings
         CALL LatticeToPoints(CP_LL(1:3,:,iW:iW), 1, CPs, iHeadP)
      enddo
      if (DEV_VERSION) then
         if ((iHeadP-1)/=size(CPs,2)) then
            print*,'PackLLPoints: Number of points wrongly estimated',size(CPs,2), iHeadP-1
            STOP ! Keep me. The check will be removed once the code is well established
         endif
      endif
      nCPs=iHeadP-1
   end subroutine

   !> Distribute the induced velocity to the proper location 
   subroutine UnPackLiftingLineVelocities()
      iHeadP=1
      do iW=1,p%nWings
         CALL VecToLattice(Uind, 1, Vind_LL(1:3,:,iW:iW), iHeadP)
      enddo
      if (DEV_VERSION) then
         if ((iHeadP-1)/=size(Uind,2)) then
            print*,'UnPackLiftingLineVelocities: Number of points wrongly estimated',size(Uind,2), iHeadP-1
            STOP ! Keep me. The check will be removed once the code is well established
         endif
      endif
   end subroutine
end subroutine

!> Fake ground effect handling to prevents vortices to enter the ground
!! For now a crude bounding is done, engineering models may follow
!! True account of the ground effect (using mirroring or panels) should be done elsewhere
!! This assumes that the ground is at z=0, in harmony with inflow wind
subroutine FakeGroundEffect(p, x, m, ErrStat, ErrMsg)
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_ContinuousStateType),   intent(inout) :: x       !< States
   type(FVW_MiscVarType),           intent(in   ) :: m       !< Initial misc/optimization variables
   integer(IntKi),                  intent(  out) :: ErrStat !< Error status of the operation
   character(*),                    intent(  out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
   integer(IntKi) :: iAge, iWing, iSpan
   integer(IntKi) :: nBelow
   real(ReKi), parameter:: GROUND         = 1.e-4_ReKi
   real(ReKi), parameter:: ABOVE_GROUND   = 0.1_ReKi
   ErrStat = ErrID_None
   ErrMsg  = ""

   nBelow=0
   do iWing = 1,p%nWings
      do iAge = 1,m%nNW+1
         do iSpan = 1,p%nSpan+1
            if (x%r_NW(3, iSpan, iAge, iWing) < GROUND) then
               x%r_NW(3, iSpan, iAge, iWing) = ABOVE_GROUND ! could use m%dxdt
               nBelow=nBelow+1
            endif
         enddo
      enddo
   enddo
   if (m%nFW>0) then
      do iWing = 1,p%nWings
         do iAge = 1,m%nFW+1
            do iSpan = 1,FWnSpan
               if (x%r_FW(3, iSpan, iAge, iWing) < GROUND) then
                  x%r_FW(3, iSpan, iAge, iWing) = ABOVE_GROUND ! could use m%dxdt
                  nBelow=nBelow+1
               endif
            enddo
         enddo
      enddo
   endif
   if (nBelow>0) then
      print*,'[WARN] Check the simulation, some vortices were found below the ground: ',nBelow
   endif
end subroutine FakeGroundEffect

!> Compute typical aerodynamic outputs based on:
!! - the lifting line velocities in global coordinates
!! - some transformation matrices
!!      - M_ag : from global to airfoil (this is well defined, also called "n-t" system in AeroDyn)
!!      - M_sg : from global to section (this is ill-defined), this coordinate is used to define the "axial" and "tangential" inductions
subroutine FVW_AeroOuts( M_sg, M_ag, PitchAndTwist, Vstr_g,  Vind_g, Vwnd_g, KinVisc, Chord, &
                         AxInd, TanInd, Vrel_norm, phi, alpha, Re, Urel_s, ErrStat, ErrMsg )
   real(ReKi),             intent(in   )  :: M_sg(3,3)               ! m%WithoutSweepPitchTwist                               global  coord to "section" coord
   real(R8Ki),             intent(in   )  :: M_ag(3,3)               ! u%BladeMotion(k)%Orientation(1:3,1:3,j)                global  coord to airfoil coord
   real(ReKi),             intent(in   )  :: PitchAndTwist           ! Pitch and twist of section
   real(ReKi),             intent(in   )  :: Vstr_g(3)               ! Structural velocity                                    global  coord
   real(ReKi),             intent(in   )  :: Vind_g(3)               ! Induced wind velocity                                  global  coord
   real(ReKi),             intent(in   )  :: Vwnd_g(3)               ! Disturbed inflow                                       global  coord
   real(ReKi),             intent(in   )  :: KinVisc                 ! Viscosity
   real(ReKi),             intent(in   )  :: Chord                   ! chord length
   real(ReKi),             intent(  out)  :: AxInd                   ! axial induction
   real(ReKi),             intent(  out)  :: TanInd                  ! Tangential induction
   real(ReKi),             intent(  out)  :: Vrel_norm               ! Relative velocity norm
   real(Reki),             intent(  out)  :: phi                     ! Flow angle
   real(Reki),             intent(  out)  :: alpha                   ! angle of attack
   real(ReKi),             intent(  out)  :: Re                      ! Reynolds number
   real(ReKi),             intent(  out)  :: Urel_s(3)               ! Relative wind of the airfoil (Vwnd - Vstr)             section coord
   integer(IntKi),         intent(  out)  :: ErrStat
   character(ErrMsgLen),   intent(  out)  :: ErrMsg

   ! Local vars
   real(ReKi)                             :: Vstr_s(3)               ! Struct Velocity,                                       section coord
   real(ReKi)                             :: Vind_s(3)               ! Induced Velocity,                                      section coord
   real(ReKi)                             :: Vwnd_s(3)               ! Disturbed wind velocity,                               section coord
   real(ReKi)                             :: Vtot_g(3)               ! Vector of total relative velocity                      section coord
   real(ReKi)                             :: Vtot_a(3)               ! Vector of total relative velocity                      global  coord
   real(ReKi)                             :: Vtot_s(3)               ! Vector of total relative velocity                      global  coord
   ErrStat = ErrID_None
   ErrMsg  = ""
   !real(DbKi), dimension(3,3)             :: M_sa                    !< Transformation matrix from airfoil to section  coord
   !real(DbKi), dimension(3,3)             :: M_sg2                   !< Transformation matrix from global  to section  coord
   ! --- Transformation from airfoil to section (KEEP ME)
   !M_sa(1,1:3) = (/  cos(PitchAndTwist*1._DbKi), sin(PitchAndTwist*1._DbKi), 0.0_DbKi /)
   !M_sa(2,1:3) = (/ -sin(PitchAndTwist*1._DbKi), cos(PitchAndTwist*1._DbKi), 0.0_DbKi /)
   !M_sa(3,1:3) = (/                   0.0_DbKi,                  0.0_DbKi, 1.0_DbKi /)
   !M_sg= matmul(M_sa, M_ag ) 

   ! --- Airfoil coordinates: used to define alpha, and Vrel, also called "n-t" system
   Vtot_g    = Vwnd_g - Vstr_g + Vind_g
   Vtot_a    = matmul(M_ag, Vtot_g)
   alpha     = atan2( Vtot_a(1), Vtot_a(2) )
   Vrel_norm = sqrt(Vtot_a(1)**2 + Vtot_a(2)**2) ! NOTE: z component shoudn't be used
   Re        = Chord * Vrel_norm / KinVisc       ! Reynolds number (not in million)

   ! Section coordinates: used to define axial induction andflow angle
   Vstr_s = matmul(M_sg, Vstr_g)
   Vind_s = matmul(M_sg, Vind_g)
   Vwnd_s = matmul(M_sg, Vwnd_g)
   Urel_s = Vwnd_s - Vstr_s          ! relative wind 
   Vtot_s = Vwnd_s - Vstr_s + Vind_s
   AxInd  = -Vind_s(1)/Urel_s(1) 
   TanInd =  Vind_s(2)/Urel_s(2)
   phi    = atan2( Vtot_s(1), Vtot_s(2) )        ! flow angle

   if(.false.) print*,PitchAndTwist ! just to avoid unused var for now
end subroutine FVW_AeroOuts

!> Generic function to compute alpha, Vrel and Re based on global data
subroutine AlphaVrel_Generic(M_ag, Vstr_g,  Vind_g, Vwnd_g, KinVisc, Chord, Vrel_norm, alpha, Re)
   real(R8Ki),             intent(in   )  :: M_ag(3,3)               ! u%BladeMotion(k)%Orientation(1:3,1:3,j)                global  coord to airfoil coord
   real(ReKi),             intent(in   )  :: Vstr_g(3)               ! Structural velocity                                    global  coord
   real(ReKi),             intent(in   )  :: Vind_g(3)               ! Induced wind velocity                                  global  coord
   real(ReKi),             intent(in   )  :: Vwnd_g(3)               ! Disturbed inflow                                       global  coord
   real(ReKi),             intent(in   )  :: KinVisc                 ! Viscosity
   real(ReKi),             intent(in   )  :: Chord                   ! chord length
   real(ReKi),             intent(  out)  :: Vrel_norm               ! Relative velocity norm
   real(Reki),             intent(  out)  :: alpha                   ! angle of attack
   real(ReKi),             intent(  out)  :: Re                      ! Reynolds number
   ! Local vars
   real(ReKi)                             :: Vtot_g(3)               ! Vector of total relative velocity                      section coord
   real(ReKi)                             :: Vtot_a(3)               ! Vector of total relative velocity                      global  coord
   ! --- Airfoil coordinates: used to define alpha, and Vrel, also called "n-t" system
   Vtot_g    = Vwnd_g - Vstr_g + Vind_g
   Vtot_a    = matmul(M_ag, Vtot_g)
   alpha     = atan2( Vtot_a(1), Vtot_a(2) )
   Vrel_norm = sqrt(Vtot_a(1)**2 + Vtot_a(2)**2) ! NOTE: z component shoudn't be used
   Re        = Chord * Vrel_norm / KinVisc       ! Reynolds number NOTE: not in million
end subroutine AlphaVrel_Generic


end module FVW_Subs
