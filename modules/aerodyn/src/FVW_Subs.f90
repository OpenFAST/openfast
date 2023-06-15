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
   integer(IntKi), parameter :: idRegDeterAuto      = 1
   integer(IntKi), parameter :: idRegDeterChord     = 2
   integer(IntKi), parameter :: idRegDeterSpan      = 3
   integer(IntKi), parameter, dimension(4) :: idRegDeterVALID      = (/idRegDeterConstant, idRegDeterAuto, idRegDeterChord, idRegDeterSpan /)
   ! Shear model
   integer(IntKi), parameter :: idShearNone   = 0
   integer(IntKi), parameter :: idShearMirror = 1
   integer(IntKi), parameter, dimension(2) :: idShearVALID         = (/idShearNone, idShearMirror /)
   ! Velocity calculation method
   integer(IntKi), parameter :: idVelocityBasic    = 1
   integer(IntKi), parameter :: idVelocityTreePart = 2
   integer(IntKi), parameter :: idVelocityPart     = 3
   integer(IntKi), parameter :: idVelocityTreeSeg  = 4
   integer(IntKi), parameter, dimension(4) :: idVelocityVALID      = (/idVelocityBasic, idVelocityTreePart, idVelocityPart,&
                                                                       idVelocityTreeSeg/)

   real(ReKi), parameter :: CoreSpreadAlpha = 1.25643

   ! Implementation
   integer(IntKi), parameter :: FWnSpan=1  !< Number of spanwise far wake panels ! TODO make it an input later
   logical       , parameter :: DEV_VERSION=.False.
   logical       , parameter :: OLAF_PROFILING=.False.
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
subroutine Output_Gamma(CP, Gamma_LL, iW, iStep, iLabel, iIter)
   real( ReKi ), dimension( :, : ), intent(in   ) :: CP       !< Control Points
   real( ReKi ), dimension( : ),    intent(in   ) :: Gamma_LL !< Circulation on the lifting line
   integer( IntKi ),                intent(in   ) :: iW    !< Wing index
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
   write(filename,'(A,I0,A,I0,A,I0,A,I0,A)')'Gamma/Gamma_step',int(iStep),'_lab',iLabel,'_it',iIter,'_Wing',int(iW),'.txt'
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

   if (p%WakeAtTE) then
      ! First panel of NW is the last lifting line panel
      do iW = 1,p%nWings
         do iSpan = 1,p%W(iW)%nSpan+1
            x%W(iW)%r_NW(1:3, iSpan, p%iNWStart-1) = m%W(iW)%r_LL(1:3, iSpan, 1)  ! iAge=1 (LL)
            x%W(iW)%r_NW(1:3, iSpan, p%iNWStart  ) = m%W(iW)%r_LL(1:3, iSpan, 2)  ! iAge=2 (TE)
         enddo
      enddo
      ! First panel of NW is the last lifting line panel
      do iW = 1,p%nWings
         do iSpan = 1,p%W(iW)%nSpan
            x%W(iW)%Gamma_NW(iSpan, p%iNWStart-1) = z%W(iW)%Gamma_LL(iSpan)  ! iAge=1
         enddo
      enddo
   else
      ! First panel of NW is the last lifting line panel
      do iW = 1,p%nWings
         do iSpan = 1,p%W(iW)%nSpan+1
            x%W(iW)%r_NW(1:3, iSpan, p%iNWStart  ) = m%W(iW)%r_LL(1:3, iSpan, 1)  ! iAge=1 (LL)
         enddo
      enddo
   endif

   ! Circulations are the same on both side of the TE
   if (p%nNWMax>p%iNWStart-1) then
      do iW = 1,p%nWings
         do iSpan = 1,p%W(iW)%nSpan
            x%W(iW)%Gamma_NW(iSpan, p%iNWStart  ) = z%W(iW)%Gamma_LL(iSpan)  ! iAge=2
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
            do iSpan = 1,p%W(iW)%nSpan
               Gamma_Prev =  x%W(iW)%Gamma_NW(iSpan, p%iNWStart+1) ! Previous circulation
               Gamma_New  =  x%W(iW)%Gamma_NW(iSpan, p%iNWStart  )
               x%W(iW)%Gamma_NW(iSpan, p%iNWStart  )  = Gamma_New*ShedScale + (1.0_ReKi-ShedScale) * Gamma_Prev
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
   integer(IntKi)            :: iW, iRoot, iTip
   real(ReKi), dimension(p%nWings) :: FWGamma
   real(ReKi), dimension(:),allocatable :: Gamma_t
   real(ReKi), dimension(:),allocatable :: sCoord
!    real(ReKi), dimension(p%W(iW)%nSpan+1) :: Gamma_t
!    real(ReKi), dimension(p%W(iW)%nSpan) :: sCoord
   real(ReKi) :: FWEpsTip, FWEpsRoot

   integer(IntKi), parameter :: iAgeFW=1   !< we update the first FW panel
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! First Panel of Farwake has coordinates of last panel of near wake always
   if (p%nFWMax>0) then
      if (m%nNW==p%nNWMax) then
         ! First circulation of Farwake is taken as the max circulation of last NW column
         FWGamma(:)=0.0_ReKi
         do iW=1,p%nWings
            allocate(Gamma_t(p%W(iW)%nSpan+1)) ! TODO TODO TODO, store as misc
            allocate(sCoord(p%W(iW)%nSpan))
            if (p%FullCircStart>0 .and. m%nFW<3) then
               ! we might run into the issue that the circulation is 0
               m%W(iW)%iTip =-1
               m%W(iW)%iRoot=-1
            endif
            ! NOTE: on the first pass, m%iTip and m%iRoot are computed
            call PlaceTipRoot(p%W(iW)%nSpan, x%W(iW)%Gamma_NW(:,m%nNW), x%W(iW)%r_NW(1:3,:,m%nNW), x%W(iW)%Eps_NW(1:3,:,m%nNW),& ! inputs
               m%W(iW)%iRoot, m%W(iW)%iTip, FWGamma(iW), FWEpsTip, FWEpsRoot) ! outputs
            x%W(iW)%Gamma_FW(1:FWnSpan,iAgeFW) = FWGamma(iW)
            x%W(iW)%Eps_FW(3,1:FWnSpan,iAgeFW) = FWEpsTip  ! HACK tip put in third
            x%W(iW)%Eps_FW(2,1:FWnSpan,iAgeFW) = FWEpsRoot ! HACK root put in second
            x%W(iW)%Eps_FW(1,1:FWnSpan,iAgeFW) = FWEpsTip  ! For shed vorticity..
            deallocate(Gamma_t)
            deallocate(sCoord)
         enddo
      endif
      ! Far wake point always mapped to last near wake
      do iW=1,p%nWings
         if (m%nNW==p%nNWMax) then
            iTip  = m%W(iW)%iTip
            iRoot = m%W(iW)%iRoot
         else
            iRoot = 1
            iTip  = p%W(iW)%nSpan+1
         endif
         x%W(iW)%r_FW(1:3,1        ,iAgeFW) =  x%W(iW)%r_NW(1:3,iRoot,p%nNWMax+1) ! Point 1 (root)
         x%W(iW)%r_FW(1:3,FWnSpan+1,iAgeFW) =  x%W(iW)%r_NW(1:3,iTip ,p%nNWMax+1) ! Point FWnSpan (tip)
         !if ((FWnSpan==2)) then
         !   ! in between point
         !   x%W(iW)%r_FW(1:3,2,iAgeFW) =  x%W(iW)%r_NW(1:3,int(p%W(iW)%nSpan+1)/4 ,p%nNWMax+1) ! Point (mid)
         !else if ((FWnSpan>2)) then
         !   ErrMsg='Error: FWnSpan>2 not implemented.'
         !   ErrStat=ErrID_Fatal
         !   return
         !endif
      enddo
   endif
   if (.false.) print*,z%W(iW)%Gamma_LL(1) ! Just to avoid unused var warning
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
               x%W(iW)%r_FW(1:3,iSpan,iAge) = x%W(iW)%r_FW(1:3,iSpan,iAge-1)
            enddo
         enddo
         x%W(iW)%r_FW(1:3,1:FWnSpan+1,1) = -999.9_ReKi ! Nullified
      enddo
   if (p%nFWMax>0) then
      do iW=1,p%nWings
         do iAge=p%nFWMax,2,-1
            do iSpan=1,FWnSpan
               x%W(iW)%Gamma_FW(iSpan,iAge) = x%W(iW)%Gamma_FW(iSpan,iAge-1)
               x%W(iW)%Eps_FW(:,iSpan,iAge) = x%W(iW)%Eps_FW(:,iSpan,iAge-1)
            enddo
         enddo
         x%W(iW)%Gamma_FW(1,1:FWnSpan-1) = -999.9_ReKi ! Nullified
         !x%W(iW)%Gamma_FW(:,1) = -999.9_ReKi ! Nullified  ! TODO TODO TODO FIX BUG
      enddo
   endif
   ! --- Propagate near wake
   do iW=1,p%nWings
      do iAge=p%nNWMax+1,p%iNWStart+1,-1
         do iSpan=1,p%W(iW)%nSpan+1
            x%W(iW)%r_NW(1:3,iSpan,iAge) = x%W(iW)%r_NW(1:3,iSpan,iAge-1)
         enddo
      enddo
      x%W(iW)%r_NW(1:3,:,1:p%iNWStart) = -999.9_ReKi ! Nullified
   enddo
   if (p%nNWMax>1) then
      do iW=1,p%nWings
         do iAge=p%nNWMax,p%iNWStart+1,-1
            do iSpan=1,p%W(iW)%nSpan
               x%W(iW)%Gamma_NW(iSpan,iAge) = x%W(iW)%Gamma_NW(iSpan,iAge-1)
               x%W(iW)%Eps_NW(:,iSpan,iAge) = x%W(iW)%Eps_NW(:,iSpan,iAge-1)
            enddo
         enddo
         x%W(iW)%Gamma_NW(:,1:p%iNWStart) = -999.9_ReKi ! Nullified
      enddo
   endif

   ! Temporary hack for sub-cycling since straight after wkae computation, the wake size will increase
   ! So we do a "fake" propagation here
   do iW=1,p%nWings
      do iAge=p%nFWMax+1,2,-1 !
         do iSpan=1,FWnSpan+1
            m%dxdt%W(iW)%r_FW(1:3,iSpan,iAge) = m%dxdt%W(iW)%r_FW(1:3,iSpan,iAge-1)
         enddo
      enddo
      !m%dxdt_FW(1:3,1:FWnSpan+1,1) = -999999_ReKi ! Important not nullified. The best would be to map the last NW convection velocity for this first row.
   enddo
   do iW=1,p%nWings
      do iAge=p%nNWMax+1,p%iNWStart+1,-1
         do iSpan=1,p%W(iW)%nSpan+1
            m%dxdt%W(iW)%r_NW(1:3,iSpan,iAge) = m%dxdt%W(iW)%r_NW(1:3,iSpan,iAge-1)
         enddo
      enddo
      m%dxdt%W(iW)%r_NW(1:3,:,1:p%iNWStart) = 0.0_ReKi ! Nullified, wing do no convect, handled by LL,NW mapping
   enddo

   if (.false.) print*,m%nNW,z%W(iW)%Gamma_LL(1) ! Just to avoid unused var warning
end subroutine PropagateWake


!> Print the states, useful for debugging
subroutine print_x_NW_FW(p, m, x, label)
   type(FVW_ParameterType),         intent(in)  :: p              !< Parameters
   type(FVW_MiscVarType),           intent(in)  :: m              !< Initial misc/optimization variables
   type(FVW_ContinuousStateType),   intent(in)  :: x              !< Continuous states
   character(len=*),intent(in) :: label
   integer(IntKi) :: iAge, iW
   character(len=1):: flag
   print*,'------------------------------------------------------------------'
   print'(A,I0,A,I0)',' NW .....................iNWStart:',p%iNWStart,' nNW:',m%nNW
   iW=1
   do iAge=1,p%nNWMax+1
      flag='X'
      if ((iAge)<= m%nNW+1) flag='.'
      print'(A,A,I0,A)',flag,'iAge ',iAge,'      Root              Tip'
      print*,trim(label)//'x', x%W(iW)%r_NW(1, 1, iAge), x%W(iW)%r_NW(1, p%W(iW)%nSpan+1, iAge)
      print*,trim(label)//'y', x%W(iW)%r_NW(2, 1, iAge), x%W(iW)%r_NW(2, p%W(iW)%nSpan+1, iAge)
      print*,trim(label)//'z', x%W(iW)%r_NW(3, 1, iAge), x%W(iW)%r_NW(3, p%W(iW)%nSpan+1, iAge)
      if (iAge<p%nNWMax+1) then
         print*,trim(label)//'g', x%W(iW)%Gamma_NW(1, iAge), x%W(iW)%Gamma_NW(p%W(iW)%nSpan, iAge)
         print*,trim(label)//'e', x%W(iW)%Eps_NW(1,1, iAge), x%W(iW)%Eps_NW(1,p%W(iW)%nSpan, iAge)
      endif
   enddo
   print'(A,I0)','FW <<<<<<<<<<<<<<<<<<<< nFW:',m%nFW
   do iAge=1,p%nFWMax+1
      flag='X'
      if ((iAge)<= m%nFW+1) flag='.'
      print'(A,A,I0,A)',flag,'iAge ',iAge,'      Root              Tip'
      print*,trim(label)//'x', x%W(iW)%r_FW(1, 1, iAge), x%W(iW)%r_FW(1, FWnSpan+1, iAge)
      print*,trim(label)//'y', x%W(iW)%r_FW(2, 1, iAge), x%W(iW)%r_FW(2, FWnSpan+1, iAge)
      print*,trim(label)//'z', x%W(iW)%r_FW(3, 1, iAge), x%W(iW)%r_FW(3, FWnSpan+1, iAge)
      if (iAge<p%nFWMax+1) then
         print*,trim(label)//'g', x%W(iW)%Gamma_FW(1,iAge), x%W(iW)%Gamma_FW(FWnSpan, iAge)
         print*,trim(label)//'e', x%W(iW)%Eps_FW(1,1, iAge), x%W(iW)%Eps_FW(1,FWnSpan, iAge)
      endif
   enddo
endsubroutine

!> Debug function to figure out if data have nan
logical function have_nan(p, m, x, z, u, label)
   type(FVW_ParameterType),         intent(in) :: p !< Parameters
   type(FVW_MiscVarType),           intent(in) :: m !< Initial misc/optimization variables
   type(FVW_ContinuousStateType),   intent(in) :: x !< Continuous states
   type(FVW_ConstraintStateType),   intent(in) :: z !< ConstrStates
   type(FVW_InputType),             intent(in) :: u(:) !< Input states
   character(len=*),                intent(in) :: label !< label for print
   integer :: iW
   have_nan=.False.
   do iW = 1,size(p%W)
      if (any(isnan(x%W(iW)%r_NW))) then
         print*,trim(label),'NaN in W(iW)%r_NW'//trim(num2lstr(iW))
         have_nan=.True.
      endif
      if (any(isnan(x%W(iW)%r_FW))) then
         print*,trim(label),'NaN in W(iW)%r_FW'//trim(num2lstr(iW))
         have_nan=.True.
      endif
      if (any(isnan(x%W(iW)%Gamma_NW))) then
         print*,trim(label),'NaN in G_NW'//trim(num2lstr(iW))
         have_nan=.True.
      endif
      if (any(isnan(x%W(iW)%Gamma_FW))) then
         print*,trim(label),'NaN in G_FW'//trim(num2lstr(iW))
         have_nan=.True.
      endif
      if (any(isnan(x%W(iW)%Eps_NW))) then
         print*,trim(label),'NaN in G_FW'//trim(num2lstr(iW))
         have_nan=.True.
      endif
      if (any(isnan(x%W(iW)%Eps_FW))) then
         print*,trim(label),'NaN in G_FW'//trim(num2lstr(iW))
         have_nan=.True.
      endif
      if (any(isnan(z%W(iW)%Gamma_LL))) then
         print*,trim(label),'NaN in G_LL'//trim(num2lstr(iW))
         have_nan=.True.
      endif
   enddo
   do iW=1,size(u)
      if (any(isnan(u(iW)%V_wind))) then
         print*,trim(label),'NaN in Vwind'//trim(num2lstr(iW))
         have_nan=.True.
      endif
   enddo
   if(.false.)print*,m%iStep ! unused var
endfunction
subroutine find_nan_1D(array, varname)
   real(ReKi), dimension(:), intent(in) :: array
   character(len=*), intent(in) :: varname
   logical :: found
   integer :: i, n, tot
   n = size(array)
   found=.false.
   tot=0
   do i =1, n
      if (isnan(array(i)))  then
         if (tot<10) then
            print*,'Position i',i
         endif
         found=.true.
         tot=tot+1
      endif
   enddo
   if (found) then 
      print*,'>>>>>>>>>>>>> NAN ',trim(varname),tot,n
      STOP
   endif
end subroutine
subroutine find_nan_2D(array, varname)
   real(ReKi), dimension(:,:), intent(in) :: array
   character(len=*), intent(in) :: varname
   logical :: found
   integer :: i, n, tot
   n = size(array,2)
   found=.false.
   tot=0
   do i =1, n
      if (any(isnan(array(:,i))))  then
         if (tot<10) then
            print*,'Position i',i
         endif
         found=.true.
         tot=tot+1
      endif
   enddo
   if (found) then 
      print*,'>>>>>>>>>>>>> NAN ',trim(varname),tot,n
      STOP
   endif
end subroutine
subroutine find_nan_3D(array, varname)
   real(ReKi), dimension(:,:,:), intent(in) :: array
   character(len=*), intent(in) :: varname
   logical :: found
   integer :: i, j, n,m,tot
   n = size(array,2)
   m = size(array,3)
   found=.false.
   tot=0
   do i =1, n
      do j =1, m
         if (any(isnan(array(:,i,j))))  then
            if (tot<10) then
               print*,'Position i,j',i,j
            endif
            found=.true.
            tot=tot+1
         endif
      enddo
   enddo
   if (found) then
      print*,'>>>>>>>>>>>>> NAN ',trim(varname), tot,n*m
      STOP
   endif
end subroutine


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
   integer(IntKi) :: iGrid,i,j,k,iW
   real(ReKi) :: xP,yP,zP,dx,dy,dz
   type(GridOutType), pointer :: g

   ! Using array reshaping to ensure a given near or far wake point is always at the same location in the array.
   ! NOTE: Maximum number of points are passed, whether they "exist" or not.
   ! NOTE: InflowWind ignores points at (0,0,0)
   !if (DEV_VERSION) then
   !   ! Removing points that don't exist
   !   !call print_x_NW_FW(p,m,x,'wind befr')
   !   if (m%nNW<=p%nNWMax) then
   !      x%W(iW)%r_NW(1:3, 1:p%W(iW)%nSpan+1, m%nNW+2:p%nNWMax+1, 1:p%nWings) = 0.0_ReKi
   !   endif
   !   if ( ((p%nNWMax<=1) .and. (m%nFW==0)) .or. ((m%nFW>0) .and. (m%nFW<=p%nFWMax))) then
   !      x%W(iW)%r_FW(1:3, 1:FWnSpan+1, m%nFW+2:p%nFWMax+1, 1:p%nWings) = 0.0_ReKi
   !   else
   !      x%W(iW)%r_FW(1:3, 1:FWnSpan+1, m%nFW+1:p%nFWMax+1, 1:p%nWings) = 0.0_ReKi
   !   endif
   !   !call print_x_NW_FW(p,m,x,'wind after')
   !endif

   iP_end=0
   ! --- LL CP
   do iW=1,p%nWings
      iP_start = iP_end+1
      iP_end   = iP_start-1 + p%W(iW)%nSpan
      r_wind(1:3,iP_start:iP_end) = m%W(iW)%CP(1:3,1:p%W(iW)%nSpan)
   enddo
   ! --- NW points
   do iW=1,p%nWings
      iP_start = iP_end+1
      iP_end   = iP_start-1+(p%W(iW)%nSpan+1)*(p%nNWMax+1)
      r_wind(1:3,iP_start:iP_end) = reshape( x%W(iW)%r_NW(1:3,1:p%W(iW)%nSpan+1,1:p%nNWMax+1) , (/ 3, (p%W(iW)%nSpan+1)*(p%nNWMax+1)/))
   enddo
   ! --- FW points
   if (p%nFWMax>0) then
      do iW=1,p%nWings
         iP_start = iP_end+1
         iP_end   = iP_start-1+(FWnSpan+1)*(p%nFWMax+1)
         r_wind(1:3,iP_start:iP_end) = reshape( x%W(iW)%r_FW(1:3,1:FWnSpan+1,1:p%nFWMax+1) , (/ 3, (FWnSpan+1)*(p%nFWMax+1) /))
      enddo
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
   !      x%W(iW)%r_NW(1:3, 1:p%W(iW)%nSpan+1, m%nNW+2:p%nNWMax+1, 1:p%nWings) = -999999.0_ReKi
   !   endif
   !   if ( ((p%nNWMax<=1) .and. (m%nFW==0)) .or. ((m%nFW>0) .and. (m%nFW<=p%nFWMax))) then
   !      x%W(iW)%r_FW(1:3, 1:FWnSpan+1, m%nFW+2:p%nFWMax+1, 1:p%nWings) =-999999.0_ReKi
   !   else
   !      x%W(iW)%r_FW(1:3, 1:FWnSpan+1, m%nFW+1:p%nFWMax+1, 1:p%nWings) =-999999.0_ReKi
   !   endif
   !endif

end subroutine SetRequestedWindPoints


!> Set the requested wind into the correponding misc variables
subroutine DistributeRequestedWind_LL(V_wind, p, m)
   real(ReKi), dimension(:,:),      intent(in   ) :: V_wind  !< Requested wind, packed
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_MiscVarType),           intent(inout) :: m       !< Misc
   !real(ReKi), dimension(:,:,:),    intent(inout) :: Vwnd_LL !< Wind on lifting line
   integer(IntKi) :: iW, iP_start,iP_end   ! Current index of point, start and end of range
   ! Using array reshaping to ensure a given near or far wake point is always at the same location in the array.
   ! NOTE: Maximum number of points are passed, whether they "exist" or not.
   iP_end=0
   ! --- LL CP
   do iW=1,p%nWings
      iP_start = iP_end+1
      iP_end   = iP_start-1 + p%W(iW)%nSpan
      m%W(iW)%Vwnd_CP(1:3,1:p%W(iW)%nSpan) = V_wind(1:3,iP_start:iP_end)
   enddo

   ! TODO TODO LL NODES
   !print*,'TODO transfer of Wind at LL'
   !do iW=1,p%nWings
   !   m%W(iW)%Vwnd_LL(1:3,1:p%W(iW)%nSpan) = m%W(iW)%Vwnd_LL(1:3,1:p%W(iW)%nSpan)
   !   m%W(iW)%Vwnd_LL(1:3,p%W(iW)%nSpan+1) = m%W(iW)%Vwnd_LL(1:3,p%W(iW)%nSpan) ! Last point copy...
   !enddo
end subroutine DistributeRequestedWind_LL

!> Distribute wind onto NW and FW
!! Modifies m%W(:)%Vwind_NW,  m%W(:)%Vwind_FW
subroutine DistributeRequestedWind_NWFW(V_wind, p, m)
   real(ReKi), dimension(:,:),      intent(in   ) :: V_wind  !< Requested wind, packed
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_MiscVarType),           intent(inout) :: m       !< Misc
   !real(ReKi), dimension(:,:,:,:),  intent(inout) :: Vwnd_NW !< Wind on near wake panels
   !real(ReKi), dimension(:,:,:,:),  intent(inout) :: Vwnd_FW !< Wind on near wake panels
   integer(IntKi) :: iW,iP_start,iP_end   ! Current index of point, start and end of range

   iP_end=0
   ! --- LL CP
   do iW=1,p%nWings
      iP_start = iP_end+1
      iP_end   = iP_start-1 + p%W(iW)%nSpan
   enddo
   ! --- NW points
   do iW=1,p%nWings
      iP_start = iP_end+1
      iP_end   = iP_start-1+(p%W(iW)%nSpan+1)*(p%nNWMax+1)
      m%W(iW)%Vwnd_NW(1:3,1:p%W(iW)%nSpan+1,1:p%nNWMax+1) = reshape( V_wind(1:3,iP_start:iP_end),(/ 3, p%W(iW)%nSpan+1, p%nNWMax+1/))
   enddo
   ! --- FW points
   if (p%nFWMax>0) then
      do iW=1,p%nWings
         iP_start = iP_end+1
         iP_end   = iP_start-1+(FWnSpan+1)*(p%nFWMax+1)
         m%W(iW)%Vwnd_FW(1:3,1:FWnSpan+1,1:p%nFWMax+1) = reshape( V_wind(1:3,iP_start:iP_end), (/ 3, FWnSpan+1, p%nFWMax+1 /))
      enddo
   endif
end subroutine DistributeRequestedWind_NWFW

!> Set the requested wind into the correponding misc variables
subroutine DistributeRequestedWind_Grid(V_wind, p, m)
   real(ReKi), dimension(:,:),      intent(in   ) :: V_wind  !< Requested wind, packed
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_MiscVarType), target,   intent(inout) :: m       !< Initial misc/optimization variables
   integer(IntKi)          :: iP_start,iP_end   ! Current index of point, start and end of range
   integer(IntKi) :: iGrid,i,j,k,iW
   type(GridOutType), pointer :: g
   iP_end=0
   ! --- LL CP
   do iW=1,p%nWings
      iP_start = iP_end+1
      iP_end   = iP_start-1 + p%W(iW)%nSpan
   enddo
   ! --- NW points
   do iW=1,p%nWings
      iP_start = iP_end+1
      iP_end   = iP_start-1+(p%W(iW)%nSpan+1)*(p%nNWMax+1)
   enddo
   ! --- FW points
   if (p%nFWMax>0) then
      do iW=1,p%nWings
         iP_start = iP_end+1
         iP_end   = iP_start-1+(FWnSpan+1)*(p%nFWMax+1)
      enddo
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



!> Init States
subroutine FVW_InitStates( x, p, ErrStat, ErrMsg )
   type(FVW_ContinuousStateType),   intent(  out)  :: x              !< States
   type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
   integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   integer(IntKi)          :: ErrStat2       ! temporary error status of the operation
   character(ErrMsgLen)    :: ErrMsg2        ! temporary error message
   character(*), parameter :: RoutineName = 'FVW_InitMiscVars'
   integer :: iW
   ! Initialize ErrStat
   ErrStat = ErrID_None
   ErrMsg  = ""

   allocate(x%W(p%nWings))
   do iW=1,p%nWings
      call AllocAry( x%W(iW)%Gamma_NW,    p%W(iW)%nSpan   , p%nNWMax  , 'NW Panels Circulation', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' ); 
      call AllocAry( x%W(iW)%Gamma_FW,    FWnSpan   , p%nFWMax  , 'FW Panels Circulation', ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' ); 
      call AllocAry( x%W(iW)%Eps_NW  , 3, p%W(iW)%nSpan   , p%nNWMax  , 'NW Panels Reg Param'  , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' );
      call AllocAry( x%W(iW)%Eps_FW  , 3, FWnSpan   , p%nFWMax  , 'FW Panels Reg Param'  , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' );
      ! set x%W(iW)%r_NW and x%W(iW)%r_FW to (0,0,0) so that InflowWind can shortcut the calculations
      call AllocAry( x%W(iW)%r_NW    , 3, p%W(iW)%nSpan+1 , p%nNWMax+1, 'NW Panels Points'     , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' );
      call AllocAry( x%W(iW)%r_FW    , 3, FWnSpan+1 , p%nFWMax+1, 'FW Panels Points'     , ErrStat2, ErrMsg2 );call SetErrStat ( ErrStat2, ErrMsg2, ErrStat,ErrMsg,'FVW_InitStates' );
      if (ErrStat >= AbortErrLev) return
      x%W(iW)%r_NW     = 0.0_ReKi
      x%W(iW)%r_FW     = 0.0_ReKi
      x%W(iW)%Gamma_NW = 0.0_ReKi ! First call of calcoutput, states might not be set 
      x%W(iW)%Gamma_FW = 0.0_ReKi ! NOTE, these values might be mapped from z%W(iW)%Gamma_LL at init
      x%W(iW)%Eps_NW   = 0.001_ReKi 
      x%W(iW)%Eps_FW   = 0.001_ReKi 
   enddo
end subroutine FVW_InitStates

subroutine FVW_InitMiscVarsPostParam( p, m, ErrStat, ErrMsg )
   type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
   type(FVW_MiscVarType),           intent(inout)  :: m              !< Initial misc/optimization variables
   integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   integer(IntKi)          :: ErrStat2       ! temporary error status of the operation
   character(ErrMsgLen)    :: ErrMsg2        ! temporary error message
   character(*), parameter :: RoutineName = 'FVW_InitMiscVarsPostParam'
   integer(IntKi) :: nSeg, nSegP, nSegNW  !< Total number of segments after packing
   integer(IntKi) :: nPart                !< Total number of particles after packing
   integer(IntKi) :: nCPs                 !< Total number of control points
   logical :: bMirror
   logical :: bLLNeedsPart, bWakeNeedsPart
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! --- Counting maximum number of segments and Control Points expected for the whole simulation
   call CountSegments(p, p%nNWMax, p%nFWMax, 1, nSeg, nSegP, nSegNW)
   nCPs = CountCPs(p, p%nNWMax, p%nFWFree)

   bMirror = p%ShearModel==idShearMirror ! Whether or not we mirror the vorticity wrt ground
   if (bMirror) then
      nSeg  = nSeg*2
      nSegP = nSegP*2
   endif
   call AllocAry( m%Sgmt%Connct, 4, nSeg , 'SegConnct' , ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%Sgmt%Connct = -999;
   call AllocAry( m%Sgmt%Points, 3, nSegP, 'SegPoints' , ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%Sgmt%Points = -999999_ReKi;
   call AllocAry( m%Sgmt%Gamma ,    nSeg,  'SegGamma'  , ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%Sgmt%Gamma  = -999999_ReKi;
   call AllocAry( m%Sgmt%Epsilon,   nSeg,  'SegEpsilon', ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%Sgmt%Epsilon= -999999_ReKi;
   m%Sgmt%nAct        = -1  ! Active segments
   m%Sgmt%nActP       = -1
   m%Sgmt%RegFunction = p%RegFunction

   bWakeNeedsPart = p%VelocityMethod(1)==idVelocityPart .or.p%VelocityMethod(1)==idVelocityTreePart
   bLLNeedsPart   = p%VelocityMethod(2)==idVelocityPart .or. p%VelocityMethod(2)==idVelocityTreePart
   if (bLLNeedsPart .or. bWakeNeedsPart) then
      nPart = 0 
      if (bWakeNeedsPart) nPart = max(nPart, nSeg * p%PartPerSegment(1))
      if (bLLNeedsPart)   nPart = max(nPart, nSeg * p%PartPerSegment(2))
      call AllocAry( m%Part%P     , 3, nPart, 'PartP'      , ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%Part%P = -999999_ReKi;
      call AllocAry( m%Part%Alpha , 3, nPart, 'PartAlpha'  , ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%Part%Alpha  = -999999_ReKi;
      call AllocAry( m%Part%RegParam,   nPart, 'PartEpsilon', ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%Part%RegParam= -999999_ReKi;
      m%Part%nAct        = -1  ! Active particles
      m%Part%RegFunction = p%RegFunction
   endif

   ! TODO Figure out Uind, CPs needed for grid
   call AllocAry( m%CPs      , 3,  nCPs, 'CPs'       , ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%CPs= -999999_ReKi;
   call AllocAry( m%Uind     , 3,  nCPs, 'Uind'      , ErrStat2, ErrMsg2 );call SetErrStat(ErrStat2, ErrMsg2, ErrStat,ErrMsg,RoutineName); m%Uind= -999999_ReKi;

end subroutine FVW_InitMiscVarsPostParam

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
   integer :: iW
   ! If the FW contains Shed vorticity, we include the last shed vorticity from the NW, otherwise, we don't!
   ! It's important not to include it, otherwise a strong vortex will be present there with no compensating vorticity from the FW
   LastNWShed = (p%FWShedVorticity ) .or. ((.not.p%FWShedVorticity) .and. (nNW<p%nNWMax))
   ! --- Counting total number of segments
   nSegP=0; nSeg=0; nSegNW=0
   ! NW segments
   if ((nNW-iDepthStart)>=0) then
      do iW=1,p%nWings
         nSegP  = nSegP +         (  (p%W(iW)%nSpan+1)*(nNW-iDepthStart+2)            )
         nSegNW = nSegNW +        (2*(p%W(iW)%nSpan+1)*(nNW-iDepthStart+2)-(p%W(iW)%nSpan+1)-(nNW-iDepthStart+1+1))
         if (.not.LastNWShed) then
            nSegNW =   nSegNW -            (p%W(iW)%nSpan) ! Removing last set of shed segments
         endif
      enddo
   endif
   nSeg=nSegNW
   ! FW segments
   if (nFW>0) then
      do iW=1,p%nWings
         nSegP  = nSegP +            (  (FWnSpan+1)*(nFW+1) )
         if (p%FWShedVorticity) then
            nSeg = nSeg +            (2*(FWnSpan+1)*(nFW+1)-(FWnSpan+1)-(nFW+1))
         else
            nSeg = nSeg +            (  (FWnSpan+1)*(nFW)                    )   ! No Shed vorticity
         endif
      enddo
   endif
end subroutine CountSegments

!> Count how many control points are convecting (needed to compute the wake convection)
pure integer(IntKi) function CountCPs(p, nNW, nFWEff) result(nCPs)
   type(FVW_ParameterType), intent(in   ) :: p       !< Parameters
   integer(IntKi),          intent(in   ) :: nNW     !< Number of NW panels
   integer(IntKi),          intent(in   ) :: nFWEff  !< Number of effective (ie. convecting) FW panels
   integer :: iW
   nCPs=0
   do iW=1,p%nWings
      nCPs = nCPs + (p%W(iW)%nSpan+1)*(nNW+1)
      if (nFWEff>0)  nCPs = nCPs + (FWnSpan+1)*(nFWEff+1)
   enddo
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
            call LatticeToSegments(x%W(iW)%r_NW(1:3,:,1:nNW+1), x%W(iW)%Gamma_NW(:,1:nNW), x%W(iW)%Eps_NW(1:3,:,1:nNW), iDepthStart, SegPoints, SegConnct, SegGamma, SegEpsilon, iHeadP, iHeadC, .True., LastNWShed, .false.)
         enddo
      endif
      if (nFW>0) then
         iHeadC_bkp = iHeadC
         do iW=1,p%nWings
            call LatticeToSegments(x%W(iW)%r_FW(1:3,:,1:nFW+1), x%W(iW)%Gamma_FW(:,1:nFW), x%W(iW)%Eps_FW(1:3,:,1:nFW), 1, SegPoints, SegConnct, SegGamma, SegEpsilon, iHeadP, iHeadC , p%FWShedVorticity, p%FWShedVorticity, .true.)
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
         if (any(SegPoints(3,:)<-999._ReKi)) then
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

      if (DEV_VERSION) then
         call find_nan_2D(SegPoints(:,1:nSegP), 'PackPanelsToSegments SegPoints')
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
   real(ReKi) :: Span !< Wing "span"/length (taken as curvilinear coordinate)
   integer :: iW, iSpan
   ErrStat = ErrID_None
   ErrMsg  = ""
   do iW=1,size(p%W)
      ! --- Compute min max and mean spanwise section lengths
      ds_min  = minval(p%W(iW)%s_LL(2:p%W(iW)%nSpan+1)-p%W(iW)%s_LL(1:p%W(iW)%nSpan))
      ds_max  = maxval(p%W(iW)%s_LL(2:p%W(iW)%nSpan+1)-p%W(iW)%s_LL(1:p%W(iW)%nSpan))
      ds_mean = sum(p%W(iW)%s_LL(2:p%W(iW)%nSpan+1)-p%W(iW)%s_LL(1:p%W(iW)%nSpan))/(p%W(iW)%nSpan+1)
      c_min  = minval(p%W(iW)%chord_LL(:))
      c_max  = maxval(p%W(iW)%chord_LL(:))
      c_mean = sum   (p%W(iW)%chord_LL(:))/(p%W(iW)%nSpan+1)
      d_min  = minval(m%W(iW)%diag_LL(:))
      d_max  = maxval(m%W(iW)%diag_LL(:))
      d_mean = sum   (m%W(iW)%diag_LL(:))/(p%W(iW)%nSpan+1)
      Span    = p%W(iW)%s_LL(p%W(iW)%nSpan+1)-p%W(iW)%s_LL(1)
      RegParam = ds_mean*2

      ! Default init of reg param
      x%W(iW)%Eps_NW(1:3,:,:) = 0.001_ReKi
      x%W(iW)%Eps_FW(1:3,:,:) = 0.001_ReKi
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
            p%WingRegParam=max(0.01_ReKi, p%WingRegParam)
            p%WakeRegParam=max(0.01_ReKi, p%WakeRegParam)
         endif

         ! Set reg param on wing and first NW
         ! NOTE: setting the same in all three directions for now, TODO!
         x%W(iW)%Eps_NW(1:3,:,1) = p%WingRegParam ! First age is always WingRegParam (LL)
         if (p%nNWMax>1) then
            x%W(iW)%Eps_NW(1:3,:,2) = p%WakeRegParam ! Second age is always WakeRegParam
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
         write(*,'(A,1F8.4)')   'WingRegParam      : ', p%WingRegParam
         write(*,'(A,1F9.4)')   'CoreSpreadEddyVisc: ', p%CoreSpreadEddyVisc
      ! Set reg param on wing and first NW
      ! NOTE: setting the same in all three directions for now, TODO!
      x%W(iW)%Eps_NW(1:3,:,1) = p%WingRegParam ! First age is always WingRegParam (LL)
      if (p%nNWMax>1) then
         x%W(iW)%Eps_NW(1:3,:,2) = p%WakeRegParam ! Second age is always WakeRegParam
      endif

      else if (p%RegDeterMethod==idRegDeterChord) then
         ! Using chord to scale the reg param
         do iSpan=1,p%W(iW)%nSpan
            x%W(iW)%Eps_NW(1:3, iSpan, 1) = p%WingRegParam * p%W(iW)%chord_CP(iSpan)
            if (p%nNWMax>1) then
               x%W(iW)%Eps_NW(1:3, iSpan, 2) = p%WakeRegParam * p%W(iW)%chord_CP(iSpan)
            endif
         enddo

      else if (p%RegDeterMethod==idRegDeterSpan) then
         ! Using dr to scale the reg param
         do iSpan=1,p%W(iW)%nSpan
            ds = p%W(iW)%s_LL(iSpan+1)-p%W(iW)%s_LL(iSpan)
            x%W(iW)%Eps_NW(1:3, iSpan, 1) = p%WingRegParam * ds
            if (p%nNWMax>1) then
               x%W(iW)%Eps_NW(1:3, iSpan, 2) = p%WakeRegParam * ds
            endif
         enddo
      else ! Should never happen (caught earlier)
         ErrStat = ErrID_Fatal
         ErrMsg ='Regularization determination method not implemented'
      endif

      if (iW==1) then
      call WrScr(' - OLAF regularization parameters (for wing 1):')
         write(*,'(A,2F8.4)') '    WingReg (min/max) : ', minval(x%W(iW)%Eps_NW(:, :, 1)), maxval(x%W(iW)%Eps_NW(:, :, 1))
         if (p%nNWMax>1) then
            write(*,'(A,2F8.4)')    '    WakeReg (min/max) : ', minval(x%W(iW)%Eps_NW(:,:, 2)), maxval(x%W(iW)%Eps_NW(:,:, 2))
         endif
         write(*,'(A,2F8.4)') '    k = alpha delta nu: ', CoreSpreadAlpha * p%CoreSpreadEddyVisc * p%KinVisc
      endif
   enddo ! Loop on wings

end subroutine FVW_InitRegularization



!> Compute induced velocities from all vortex elements onto nPoints
!! In : x, x%W(iW)%r_NW, x%W(iW)%r_FW, x%W(iW)%Gamma_NW, x%W(iW)%Gamma_FW
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
   real(ReKi), dimension(:,:), allocatable :: CPs  ! TODO get rid of me with dedicated functions
   real(ReKi), dimension(:,:), allocatable :: Uind ! TODO get rid of me with dedicated functions
   ErrStat= ErrID_None
   ErrMsg =''

   ! --- Packing control points
   nCPs = g%nx * g%ny * g%nz
   allocate(CPs(3, nCPs), stat=ErrStat)
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
   allocate(Uind(3, nCPs), stat=ErrStat); Uind=0.0_ReKi
   iHeadP=1
   call FlattenValues(g%uGrid, Uind, iHeadP); ! NOTE: Uind contains uGrid now (Uwnd)

   ! --- Compute induced velocity
   ! Convert Panels to segments, segments to particles, particles to tree
   call InducedVelocitiesAll_Init(p, x, m, m%Sgmt, m%Part, Tree, ErrStat, ErrMsg, allocPart=.false.)
   call InducedVelocitiesAll_Calc(CPs, nCPs, Uind, p, m%Sgmt, m%Part, Tree, ErrStat, ErrMsg)
   call InducedVelocitiesAll_End(p, Tree, m%Part, ErrStat, ErrMsg, deallocPart=.false.)

   ! --- Unpacking induced velocity points
   iHeadP=1
   call DeflateValues(Uind, g%uGrid, iHeadP)

   if(allocated(CPs )) deallocate(CPs , stat=ErrStat)
   if(allocated(Uind)) deallocate(Uind, stat=ErrStat)

end subroutine InducedVelocitiesAll_OnGrid

!> Wrapper to setup part from set of segments
subroutine SegmentsToPartWrap(Sgmt, nSeg, PartPerSegment, RegFunction, Part, allocPart)
   type(T_Sgmt),                    intent(in   ) :: Sgmt  !< Segments
   integer(IntKi),                  intent(in   ) :: nSeg  !< Number of segments to use (might not use all of them)
   integer(IntKi),                  intent(in   ) :: PartPerSegment !< Number of particles per segment
   integer(IntKi),                  intent(in   ) :: RegFunction    !< Regularization function
   type(T_Part),                    intent(inout) :: Part  !< Particles
   logical,                         intent(in   ) :: allocPart !< allocate particles
   integer(IntKi) :: iHeadP
   integer(IntKi) :: nPart
   logical, parameter :: alloc  =.true.  !< Should we allocate the particles?
   iHeadP=1
   nPart = PartPerSegment * nSeg
   ! --- Allocate
   if (allocPart)  then
      if (allocated(Part%P))        deallocate(Part%P)
      if (allocated(Part%Alpha))    deallocate(Part%Alpha)
      if (allocated(Part%RegParam)) deallocate(Part%RegParam)
      allocate(Part%P(3,nPart), Part%Alpha(3,nPart), Part%RegParam(nPart)) ! NOTE: remember to deallocate
   else
      ! check that we have enough space
      if (.not. allocated(Part%P)) then
          print*,'>>> PartP not allocated'; 
          STOP
      endif
      if (size(Part%P,2)<nPart) then
           print*,'>>> PartP storage too small';
           STOP
      endif
   endif
   Part%P(:,:)      = -99999.99_ReKi
   Part%Alpha(:,:)  = -99999.99_ReKi
   Part%RegParam(:) = -99999.99_ReKi
   Part%nAct = nPart ! TODO add iHeadPart  if particles already present

   call SegmentsToPart(Sgmt%Points, Sgmt%Connct, Sgmt%Gamma, Sgmt%Epsilon, 1, nSeg, PartPerSegment, Part%P, Part%Alpha, Part%RegParam, iHeadP)
   if (RegFunction/=idRegNone) then
      Part%RegFunction = idRegExp ! TODO need to find a good equivalence and potentially adapt Epsilon in SegmentsToPart
   endif
   if (DEV_VERSION) then
      call find_nan_2D(Part%P    , 'SegmentsToPartWrap Part%P')
      call find_nan_2D(Part%Alpha, 'SegmentsToPartWrap Part%Alpha')
      if (any(Part%RegParam(:)<-9999.99_ReKi)) then
         print*,'Error in Segment to part conversion'
         STOP
      endif
   endif
end subroutine SegmentsToPartWrap

!> Perform initialization steps before requesting induced velocities from All vortex elements
!! In : x%W(iW)%r_NW, x%W(iW)%r_FW, x%W(iW)%Gamma_NW, x%W(iW)%Gamma_FW
!! Out: Tree, Part, m
subroutine InducedVelocitiesAll_Init(p, x, m, Sgmt, Part, Tree,  ErrStat, ErrMsg, allocPart)
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_ContinuousStateType),   intent(in   ) :: x       !< States
   type(FVW_MiscVarType),           intent(in   ) :: m       !< Misc
   type(T_Sgmt),                    intent(inout) :: Sgmt    !< Segments
   type(T_Part),                    intent(inout) :: Part    !< Particle storage if needed
   type(T_Tree),                    intent(out)   :: Tree    !< Tree of particles if needed
   integer(IntKi),                  intent(  out) :: ErrStat !< Error status of the operation
   character(*),                    intent(  out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
   logical,                         intent(in   ) :: allocPart !< allocate particles
   integer, parameter :: iVel = 1
   ! Local variables
   integer(IntKi) :: nSeg, nSegP
   logical        :: bMirror ! True if we mirror the vorticity wrt ground
   ErrStat= ErrID_None
   ErrMsg =''

   bMirror = p%ShearModel==idShearMirror ! Whether or not we mirror the vorticity wrt ground

   ! --- Packing all vortex elements into a list of segments
   call PackPanelsToSegments(p, x, 1, bMirror, m%nNW, m%nFW, Sgmt%Connct, Sgmt%Points, Sgmt%Gamma, Sgmt%Epsilon, nSeg, nSegP)
   Sgmt%RegFunction=p%RegFunction
   Sgmt%nAct  = nSeg
   Sgmt%nActP = nSegP

   ! --- Convert to particles if needed
   if ((p%VelocityMethod(iVel)==idVelocityTreePart) .or. (p%VelocityMethod(iVel)==idVelocityPart)) then
      call SegmentsToPartWrap(Sgmt, nSeg, p%PartPerSegment(iVel), p%RegFunction, Part, allocPart=allocPart)
   endif

   ! --- Grow tree if needed
   if (p%VelocityMethod(iVel)==idVelocityTreePart) then
      call grow_tree_part(Tree, Part%nAct, Part%P, Part%Alpha, Part%RegFunction, Part%RegParam, 0)

   elseif (p%VelocityMethod(iVel)==idVelocityTreeSeg) then
      call grow_tree_segment(Tree, nSeg, Sgmt%Points, Sgmt%Connct(:,1:nSeg), Sgmt%Gamma(1:nSeg), p%RegFunction, Sgmt%Epsilon(1:nSeg), 0)
   endif

end subroutine InducedVelocitiesAll_Init

!> Compute induced velocity on flat CPs
subroutine InducedVelocitiesAll_Calc(CPs, nCPs, Uind, p, Sgmt, Part, Tree, ErrStat, ErrMsg)
   real(ReKi), dimension(:,:),      intent(in)    :: CPs     !< Control points (3 x nCPs++)
   integer(IntKi)                 , intent(in)    :: nCPs    !< Number of control points on which to compute (nCPs <= size(CPs,2))
   real(ReKi), dimension(:,: )    , intent(inout) :: Uind    !< Induced velocity vector - Side effects!!! (3 x nCPs++)
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(T_Sgmt),                    intent(in   ) :: Sgmt    !< Segments
   type(T_Part),                    intent(in   ) :: Part    !< Particle storage if needed
   type(T_Tree),                    intent(inout) :: Tree    !< Tree of particles if needed
   integer(IntKi),                  intent(  out) :: ErrStat !< Error status of the operation
   character(*),                    intent(  out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
   integer, parameter :: iVel = 1
   ! Local variables
   ErrStat= ErrID_None
   ErrMsg =''

   if (p%VelocityMethod(iVel)==idVelocityBasic) then
      call ui_seg( 1, nCPs, CPs, 1, Sgmt%nAct, Sgmt%Points, Sgmt%Connct, Sgmt%Gamma, Sgmt%RegFunction, Sgmt%Epsilon, Uind)

   elseif (p%VelocityMethod(iVel)==idVelocityTreePart) then
      ! Tree has already been grown with InducedVelocitiesAll_Init
      !call print_tree(Tree)
      call ui_tree_part(Tree, nCPs, CPs, p%TreeBranchFactor(iVel), Tree%DistanceDirect, Uind, ErrStat, ErrMsg)

   elseif (p%VelocityMethod(iVel)==idVelocityPart) then
      call ui_part_nograd(nCPs, CPs, Part%nAct, Part%P, Part%Alpha, Part%RegFunction, Part%RegParam, Uind)

   elseif (p%VelocityMethod(iVel)==idVelocityTreeSeg) then
      call ui_tree_segment(Tree, CPs, nCPs, p%TreeBranchFactor(iVel), Tree%DistanceDirect, Uind, ErrStat, ErrMsg)
   endif
end subroutine InducedVelocitiesAll_Calc


!> Perform termination steps after velocity was requested from all vortex elements
!! InOut: Tree, Part, m
subroutine InducedVelocitiesAll_End(p, Tree, Part, ErrStat, ErrMsg, deallocPart)
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(T_Tree),                    intent(inout) :: Tree    !< Tree of particles if needed
   type(T_Part),                    intent(inout) :: Part    !< Particle storage if needed
   integer(IntKi),                  intent(  out) :: ErrStat !< Error status of the operation
   character(*),                    intent(  out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
   logical,                         intent(in   ) :: deallocPart
   integer, parameter :: iVel = 1
   ! Local variables
   ErrStat= ErrID_None
   ErrMsg =''

   if (p%VelocityMethod(iVel)==idVelocityBasic) then
      ! Nothing

   elseif (p%VelocityMethod(iVel)==idVelocityTreePart) then
      if (deallocPart) deallocate(Part%P, Part%Alpha, Part%RegParam)
      call cut_tree(Tree)

   elseif (p%VelocityMethod(iVel)==idVelocityPart) then
      if (deallocPart) deallocate(Part%P, Part%Alpha, Part%RegParam)

   elseif (p%VelocityMethod(iVel)==idVelocityTreeSeg) then
      call cut_tree(Tree) ! We do not deallocate segment
   endif

end subroutine InducedVelocitiesAll_End




!> Compute induced velocities from all vortex elements onto all the vortex elements
!! In : x%W(iW)%r_NW, x%W(iW)%r_FW, x%W(iW)%Gamma_NW, x%W(iW)%Gamma_FW
!! Out: m%W(iW)%Vind_NW, m%Vind_FW
subroutine WakeInducedVelocities(p, x, m, ErrStat, ErrMsg)
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_ContinuousStateType),   intent(in   ) :: x       !< States
   type(FVW_MiscVarType),           intent(inout) :: m       !< Initial misc/optimization variables
   integer(IntKi),                  intent(  out) :: ErrStat !< Error status of the operation
   character(*),                    intent(  out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi) :: iW, nCPs, iHeadP
   integer(IntKi) :: nFWEff  ! Number of farwake panels that are free at current time step
   integer(IntKi) :: nNWEff  ! Number of nearwake panels that are free at current time step
   type(T_Tree)   :: Tree
   if (OLAF_PROFILING) call tic('WakeInduced Calc')
   ErrStat= ErrID_None
   ErrMsg =''

   nFWEff = min(m%nFW, p%nFWFree)
   nNWEff = min(m%nNW, p%nNWFree)

   ! --- Pack control points
   call PackConvectingPoints() ! m%CPs

   ! --- Compute induced velocity
   ! Convert Panels to segments, segments to particles, particles to tree
   m%Uind=0.0_ReKi ! very important due to side effects of ui_* methods
   m%Uind(:,nCPs+1:)=1000.0_ReKi ! TODO For debugging only
   call InducedVelocitiesAll_Init(p, x, m, m%Sgmt, m%Part, Tree, ErrStat, ErrMsg, allocPart=.false.)
   call InducedVelocitiesAll_Calc(m%CPs, nCPs, m%Uind, p, m%Sgmt, m%Part, Tree, ErrStat, ErrMsg)
   call InducedVelocitiesAll_End(p, Tree, m%Part, ErrStat, ErrMsg, deallocPart=.false.)
   call UnPackInducedVelocity()

   if (DEV_VERSION) then
      print'(A,I0,A,I0,A,I0)','Convection - nSeg:',m%Sgmt%nAct,' - nSegP:',m%Sgmt%nActP, ' - nCPs:',nCPs
   endif
   if (OLAF_PROFILING) call toc()
contains
   !> Pack all the points that convect
   subroutine PackConvectingPoints()
      ! Counting total number of control points that convects
      nCPs = CountCPs(p, nNWEff, nFWEff)
      m%CPs=-999.9_ReKi
      ! Packing
      iHeadP=1
      do iW=1,p%nWings
         CALL LatticeToPoints(x%W(iW)%r_NW(1:3,:,1:nNWEff+1), 1, m%CPs, iHeadP)
      enddo
      if (nFWEff>0) then
         do iW=1,p%nWings
            CALL LatticeToPoints(x%W(iW)%r_FW(1:3,:,1:nFWEff+1), 1, m%CPs, iHeadP)
         enddo
      endif
      if (DEV_VERSION) then
         ! Additional checks
         if (any(m%CPs(1,1:nCPs)<=-99)) then
            call print_x_NW_FW(p,m,x,'pack')
            ErrMsg='PackConvectingPoints: Problem in Control points'; ErrStat=ErrID_Fatal; return
         endif
         if (p%nNWMax==p%nNWFree) then
            ! Number of CP should be number of SegP
            if ((iHeadP-1)/=nCPs) then
               print*,'PackConvectingPoints: Number of points wrongly estimated',nCPs, iHeadP-1
               STOP ! Keep me. The check will be removed once the code is well established
               ErrMsg='PackConvectingPoints: Number of points wrongly estimated '; ErrStat=ErrID_Fatal; return
            endif
         endif
         call find_nan_2D(m%CPs(:,1:nCPs), 'WakeInducedVel CPs')
      endif
   end subroutine
   !> Distribute the induced velocity to the proper location
   subroutine UnPackInducedVelocity()
      do iW=1,p%nWings
         m%W(iW)%Vind_NW = -9999._ReKi !< Safety
         m%W(iW)%Vind_FW = -9999._ReKi !< Safety
         m%W(iW)%Vind_NW(:,:,p%nNWFree+1:) = 2222._ReKi !< Safety
      enddo
      iHeadP=1
      do iW=1,p%nWings
         CALL VecToLattice(m%Uind, 1, m%W(iW)%Vind_NW(:,:,1:nNWEff+1), iHeadP)
      enddo
      if (nFWEff>0) then
         do iW=1,p%nWings
            CALL VecToLattice(m%Uind, 1, m%W(iW)%Vind_FW(1:3,1:FWnSpan+1,1:nFWEff+1), iHeadP)
         enddo
         if (DEV_VERSION) then
            do iW=1,p%nWings
               if (any(m%W(iW)%Vind_FW(1:3,1:FWnSpan+1,1:nFWEff+1)<-99)) then
                  ErrMsg='UnPackInducedVelocity: Problem in FW induced velocity on FW points'; ErrStat=ErrID_Fatal; return
               endif
            enddo
         endif
      endif
      if (DEV_VERSION) then
         if (p%nNWMax==p%nNWFree) then
            ! Number of CP should be number of SegP
            if ((iHeadP-1)/=nCPs) then
               print*,'UnPackInducedVelocity: Number of points wrongly estimated',nCPs, iHeadP-1
               STOP ! Keep me. The check will be removed once the code is well established
               ErrMsg='UnPackInducedVelocity: Number of points wrongly estimated'; ErrStat=ErrID_Fatal; return
            endif
         endif
         call find_nan_2D(m%Uind(:,:), 'WakeInducedVel Uind')
      endif
   end subroutine

end subroutine WakeInducedVelocities

!> Compute induced velocities from all vortex elements onto the lifting line control points
!! In : x%W(iW)%r_NW, x%W(iW)%r_FW, x%W(iW)%Gamma_NW, x%W(iW)%Gamma_FW
!! Out: m%W(iW)%Vind_CP
subroutine LiftingLineInducedVelocities(p, x, InductionAtCP, iDepthStart, m, ErrStat, ErrMsg)
   !real(ReKi), dimension(:,:,:),    intent(in   ) :: CP   !< Control points where velocity is to be evaluated
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_ContinuousStateType),   intent(in   ) :: x       !< States
   logical,                         intent(in   ) :: InductionAtCP !< Compute induction at CP or on LL nodes
   integer(IntKi),                  intent(in   ) :: iDepthStart !< Index where we start packing for NW panels
   type(FVW_MiscVarType),           intent(inout) :: m       !< Initial misc/optimization variables
   integer(IntKi),                  intent(  out) :: ErrStat    !< Error status of the operation
   character(*),                    intent(  out) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi) :: iW, nSeg, nSegP, nCPs, iHeadP
   real(ReKi) :: MaxWingLength, DistanceDirect !< Maximum wing length, used to determined distance for direct evaluation of tree
   real(ReKi),    dimension(:,:), allocatable :: CPs   !< ControlPoints
   real(ReKi),    dimension(:,:), allocatable :: Uind  !< Induced velocity
   type(T_Tree) :: Tree !< Tree of particles/segment if needed
   integer, parameter :: iVel = 2
   logical      :: bMirror
   if (OLAF_PROFILING) call tic('LiftingLine UI Calc')
   ErrStat = ErrID_None
   ErrMsg  = ""

   do iW=1,p%nWings
      m%W(iW)%Vind_CP = -9999._ReKi !< Safety
      m%W(iW)%Vind_LL = -9999._ReKi !< Safety
   enddo
   bMirror = p%ShearModel==idShearMirror ! Whether or not we mirror the vorticity wrt ground

   ! --- Packing all vortex elements into a list of segments
   call PackPanelsToSegments(p, x, iDepthStart, bMirror, m%nNW, m%nFW, m%Sgmt%Connct, m%Sgmt%Points, m%Sgmt%Gamma, m%Sgmt%Epsilon, nSeg, nSegP)
   m%Sgmt%RegFunction=p%RegFunction
   m%Sgmt%nAct  = nSeg
   m%Sgmt%nActP = nSegP

   ! --- Computing induced velocity
   if (nSegP==0) then
      nCPs=0
      do iW=1,p%nWings
         m%W(iW)%Vind_CP = 0.0_ReKi !< Safety
         m%W(iW)%Vind_LL = 0.0_ReKi !< Safety
      enddo
      if (DEV_VERSION) then
         print'(A,I0,A,I0,A,I0,A)','Induction -  nSeg:',nSeg,' - nSegP:',nSegP, ' - nCPs:',nCPs, ' -> No induction'
      endif
   else
      nCPs=0
      if (InductionAtCP) then
         do iW=1,p%nWings
            nCPs = nCPs + p%W(iW)%nSpan
         enddo
      else
         do iW=1,p%nWings
            nCPs = nCPs + p%W(iW)%nSpan+1
         enddo
      endif

      allocate(CPs (1:3,1:nCPs)) ! NOTE: here we do allocate CPs and Uind insteadof using Misc
      allocate(Uind(1:3,1:nCPs)) !       The size is reasonably small, and m%Uind then stay filled with "rollup velocities" (for export)
      Uind=0.0_ReKi !< important due to side effects of ui_seg

      ! --- Pack
      call PackLiftingLinePoints()
      if (DEV_VERSION) then
         print'(A,I0,A,I0,A,I0)','Induction -  nSeg:',nSeg,' - nSegP:',nSegP, ' - nCPs:',nCPs
      endif

      ! --- Compute maximum wing length
      MaxWingLength = 0.0_ReKi
      do iW=1,p%nWings
         MaxWingLength = max(MaxWingLength,  p%W(iW)%s_LL(p%W(iW)%nSpan+1)-p%W(iW)%s_LL(1)) ! Using curvilinear variable for length...
      enddo
      DistanceDirect = MaxWingLength*2.2_ReKi ! Using ~2*R+margin so that an entire rotor will be part of a direct evaluation

      ! --- Compute velocity on LL
      ! TreeSeg is faster but introduce some noise, so we keep this open for the user to choose
      if (p%VelocityMethod(iVel) == idVelocityBasic) then 
         call ui_seg( 1, nCPs, CPs, 1, nSeg, m%Sgmt%Points, m%Sgmt%Connct, m%Sgmt%Gamma, m%Sgmt%RegFunction, m%Sgmt%Epsilon, Uind)

      else if (p%VelocityMethod(iVel) == idVelocityPart) then 
         call SegmentsToPartWrap(m%Sgmt, nSeg, p%PartPerSegment(iVel), p%RegFunction, m%Part, allocPart=.false.)
         call ui_part_nograd(nCPs, CPs, m%Part%nAct, m%Part%P, m%Part%Alpha, m%Part%RegFunction, m%Part%RegParam, Uind)
         !deallocate(Part%P, Part%Alpha, Part%RegParam)

      else if (p%VelocityMethod(iVel) == idVelocityTreeSeg) then 
         call grow_tree_segment(Tree, nSeg, m%Sgmt%Points, m%Sgmt%Connct(:,1:nSeg), m%Sgmt%Gamma(1:nSeg), m%Sgmt%RegFunction, m%Sgmt%Epsilon(1:nSeg), 0)
         call ui_tree_segment(Tree, CPs, nCPs, p%TreeBranchFactor(iVel), DistanceDirect, Uind, ErrStat, ErrMsg)
         call cut_tree(Tree)

      else if (p%VelocityMethod(iVel) == idVelocityTreePart) then 
         call SegmentsToPartWrap(m%Sgmt, nSeg, p%PartPerSegment(iVel), p%RegFunction, m%Part, allocPart=.false.)
         call grow_tree_part(Tree, m%Part%nAct, m%Part%P, m%Part%Alpha, m%Part%RegFunction, m%Part%RegParam, 0)
         call ui_tree_part(Tree, nCPs, CPs, p%TreeBranchFactor(iVel), DistanceDirect, Uind, ErrStat, ErrMsg)
         !deallocate(Part%P, Part%Alpha, Part%RegParam)
         call cut_tree(Tree)
      endif

      ! --- Unpack
      call UnPackLiftingLineVelocities()

      deallocate(Uind)
      deallocate(CPs)
   endif
   if (OLAF_PROFILING) call toc()
contains
   !> Pack all the control points
   subroutine PackLiftingLinePoints()
      iHeadP=1
      if (InductionAtCP) then
         do iW=1,p%nWings
            call LatticeToPoints2D(m%W(iW)%CP(1:3,:), CPs, iHeadP)
         enddo
      else
         do iW=1,p%nWings
            call LatticeToPoints2D(m%W(iW)%r_LL(1:3,:,1), CPs, iHeadP)
         enddo
      endif
      if (DEV_VERSION) then
         if ((iHeadP-1)/=size(CPs,2)) then
            print*,'PackLLPoints: Number of points wrongly estimated',size(CPs,2), iHeadP-1
            STOP ! Keep me. The check will be removed once the code is well established
         endif
         call find_nan_2D(CPs, 'LiftingLineInducedVel CPs')
      endif
      nCPs=iHeadP-1
   end subroutine

   !> Distribute the induced velocity to the proper location
   subroutine UnPackLiftingLineVelocities()
      integer :: iSpan
      iHeadP=1
      if (InductionAtCP) then
         do iW=1,p%nWings
            call VecToLattice2D(Uind, m%W(iW)%Vind_CP(1:3,:), iHeadP)
         enddo
         ! --- Transfer CP to LL (Linear interpolation for interior points and extrapolations at boundaries)
         do iW=1,p%nWings
            call interpextrap_cp2node(p%W(iW)%s_CP(:), m%W(iW)%Vind_CP(1,:), p%W(iW)%s_LL(:), m%W(iW)%Vind_LL(1,:))
            call interpextrap_cp2node(p%W(iW)%s_CP(:), m%W(iW)%Vind_CP(2,:), p%W(iW)%s_LL(:), m%W(iW)%Vind_LL(2,:))
            call interpextrap_cp2node(p%W(iW)%s_CP(:), m%W(iW)%Vind_CP(3,:), p%W(iW)%s_LL(:), m%W(iW)%Vind_LL(3,:))
         enddo
      else
         do iW=1,p%nWings
            call VecToLattice2D(Uind, m%W(iW)%Vind_LL(1:3,:), iHeadP)
         enddo
         ! --- Transfer LL to CP. TODO instead of mean should use weigthed average based on distance to nodes
         do iW=1,p%nWings
            do iSpan=1,p%W(iW)%nSpan
               m%W(iW)%Vind_CP(1:3,iSpan)= (m%W(iW)%Vind_LL(1:3,iSpan)+m%W(iW)%Vind_LL(1:3,iSpan+1))*0.5_ReKi
            enddo
         enddo
      endif
      if (DEV_VERSION) then
         if ((iHeadP-1)/=size(Uind,2)) then
            print*,'UnPackLiftingLineVelocities: Number of points wrongly estimated',size(Uind,2), iHeadP-1
            STOP ! Keep me. The check will be removed once the code is well established
         endif
         call find_nan_2D(Uind, 'LiftingLineInducedVel Uind')
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
   integer(IntKi) :: iAge, iW, iSpan
   integer(IntKi) :: nBelow
   real(ReKi) :: GROUND
   real(ReKi) :: ABOVE_GROUND
   ErrStat = ErrID_None
   ErrMsg  = ""

   if ( p%MHK == 1 .or. p%MHK == 2 ) then
      GROUND         = 1.e-4_ReKi - p%WtrDpth
      ABOVE_GROUND   = 0.1_ReKi - p%WtrDpth
   else
      GROUND         = 1.e-4_ReKi
      ABOVE_GROUND   = 0.1_ReKi
   endif

   nBelow=0
   do iW = 1,p%nWings
      do iAge = 1,m%nNW+1
         do iSpan = 1,p%W(iW)%nSpan+1
            if (x%W(iW)%r_NW(3, iSpan, iAge) < GROUND) then
               x%W(iW)%r_NW(3, iSpan, iAge) = ABOVE_GROUND ! could use m%dxdt
               nBelow=nBelow+1
            endif
         enddo
      enddo
   enddo
   if (m%nFW>0) then
      do iW = 1,p%nWings
         do iAge = 1,m%nFW+1
            do iSpan = 1,FWnSpan
               if (x%W(iW)%r_FW(3, iSpan, iAge) < GROUND) then
                  x%W(iW)%r_FW(3, iSpan, iAge) = ABOVE_GROUND ! could use m%dxdt
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
   real(R8Ki),             intent(in   )  :: M_sg(3,3)               ! m%WithoutSweepPitchTwist                               global  coord to "section" coord
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
   if (EqualRealNos(Urel_s(1),0.0_ReKi)) then
      AxInd = 0.0_ReKi
   else
      AxInd  = -Vind_s(1)/Urel_s(1)
   endif
   if (EqualRealNos(Urel_s(2),0.0_ReKi)) then
      TanInd = 0.0_ReKi
   else
      TanInd =  Vind_s(2)/Urel_s(2)
   end if
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
