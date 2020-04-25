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
   integer(IntKi), parameter, dimension(1) :: idIntMethodVALID      = (/idEuler1 /)
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
   integer(IntKi), parameter :: idRegDeterManual  = 0
   integer(IntKi), parameter :: idRegDeterAuto    = 1
   integer(IntKi), parameter, dimension(2) :: idRegDeterVALID      = (/idRegDeterManual, idRegDeterAuto /)
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
   integer(IntKi)      :: nLines
   integer(IntKi)      :: i
   integer(IntKi)      :: iStat
   integer(IntKi)      :: iUnit
   character(len=1054) :: line
   real(ReKi), dimension(:), allocatable :: sPrescr, GammaPrescr !< Radius
   ErrStat = ErrID_None
   ErrMsg  = ''
   ! --- 
   call GetNewUnit(iUnit)
   open(unit = iUnit, file = CirculationFileName)
   nLines=line_count(iUnit)-1
   ! Read Header
   read(iUnit,*, iostat=istat) line 
   ! Read table:  s/L [-], GammaPresc [m^2/s]
   allocate(sPrescr(1:nLines), GammaPrescr(1:nLines))
   do i=1,nLines
      read(iUnit,*, iostat=istat) sPrescr(i), GammaPrescr(i)
      sPrescr(i)     =   sPrescr(i) * L
      GammaPrescr(i) =   GammaPrescr(i) 
   enddo
   close(iUnit)
   if (istat/=0) then
      ErrStat=ErrID_Fatal
      ErrMsg='Error occured while reading Circulation file: '//trim(CirculationFileName)
      return
   endif
   ! NOTE: TODO TODO TODO THIS ROUTINE PERFORMS NASTY EXTRAPOLATION, SHOULD BE PLATEAUED
   Gamma_CP_LL =  interpolation_array( sPrescr, GammaPrescr, s_CP_LL, size(s_CP_LL), nLines )
contains

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
   type(FVW_MiscVarType),           intent(in   )  :: m              !< Initial misc/optimization variables
   type(FVW_ConstraintStateType),   intent(in   )  :: z              !< Constraints states
   type(FVW_ContinuousStateType),   intent(inout)  :: x              !< Continuous states
   integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
   character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
   integer(IntKi)            :: iW, iRoot
   real(ReKi), dimension(p%nWings) :: FWGamma
   integer(IntKi), parameter :: iAgeFW=1   !< we update the first FW panel
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   ! First Panel of Farwake has coordinates of last panel of near wake always
   if (p%nFWMax>0) then
      FWGamma(:)=0.0_ReKi
      if (m%nNW==p%nNWMax) then
         ! First circulation of Farwake is taken as the max circulation of last NW column
         do iW=1,p%nWings
            !FWGamma = sum(x%Gamma_NW(:,p%nNWMax,iW))/p%nSpan
            FWGamma(iW) = maxval(x%Gamma_NW(:,p%nNWMax,iW))
            x%Gamma_FW(1:FWnSpan,iAgeFW,iW) = FWGamma(iW)
         enddo
      endif

      do iW=1,p%nWings
         ! Find first point (in half span) where circulation is more than 0.1% of MaxGamma, call it the root
         iRoot=1
         ! NOTE: this below won't work for a wing
         ! Need to go from maxgamma location, and integrate spanwise position on both side to find location of tip and root vortex
         !do while ((iRoot<int(p%nSpan/2)) .and. (x%Gamma_NW(iRoot, p%nNWMax,iW)< 0.001*FWGamma(iW) ))
         !   iRoot=iRoot+1
         !enddo

         x%r_FW(1:3,1        ,iAgeFW,iW) =  x%r_NW(1:3,iRoot     ,p%nNWMax+1,iW) ! Point 1 (root)
         x%r_FW(1:3,FWnSpan+1,iAgeFW,iW) =  x%r_NW(1:3,p%nSpan+1 ,p%nNWMax+1,iW) ! Point FWnSpan (tip)
         if ((FWnSpan==2)) then
            ! in between point
            x%r_FW(1:3,2,iAgeFW,iW) =  x%r_NW(1:3,int(p%nSpan+1)/4 ,p%nNWMax+1,iW) ! Point (mid)
         else if ((FWnSpan>2)) then
            ErrMsg='Error: FWnSpan>2 not implemented.'
            ErrStat=ErrID_Fatal
            return
         endif
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
            enddo
         enddo
         x%Gamma_FW(1,1:FWnSpan-1,iW) = -999.9_ReKi ! Nullified
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
            m%dxdt_FW(1:3,iSpan,iAge,iW) = m%dxdt_FW(1:3,iSpan,iAge-1,iW)
         enddo
      enddo
      !m%dxdt_FW(1:3,1:FWnSpan+1,1,iW) = -999999_ReKi ! Important not nullified. The best would be to map the last NW convection velocity for this first row.
   enddo
   do iW=1,p%nWings
      do iAge=p%nNWMax+1,iNWStart+1,-1 
         do iSpan=1,p%nSpan+1
            m%dxdt_NW(1:3,iSpan,iAge,iW) = m%dxdt_NW(1:3,iSpan,iAge-1,iW)
         enddo
      enddo
      m%dxdt_NW(1:3,:,1:iNWStart,iW) = 0.0_ReKi ! Nullified, wing do no convect, handled by LL,NW mapping
   enddo

   if (.false.) print*,m%nNW,z%Gamma_LL(1,1) ! Just to avoid unused var warning
end subroutine PropagateWake


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
   enddo
   print'(A,I0)','FW <<<<<<<<<<<<<<<<<<<< nFW:',m%nFW
   do iAge=1,p%nFWMax+1
      flag='X'
      if ((iAge)<= m%nFW+1) flag='.'
      print'(A,A,I0,A)',flag,'iAge ',iAge,'      Root              Tip'
      print*,trim(label)//'x', x%r_FW(1, 1, iAge,1), x%r_FW(1, FWnSpan+1, iAge,1)
      print*,trim(label)//'y', x%r_FW(2, 1, iAge,1), x%r_FW(2, FWnSpan+1, iAge,1)
      print*,trim(label)//'z', x%r_FW(3, 1, iAge,1), x%r_FW(3, FWnSpan+1, iAge,1)
   enddo
   !print'(A,I0,A,I0)','dxdt NW .....................iNWStart:',iNWStart,' nNW:',m%nNW
   !do iAge=1,p%nNWMax+1
   !   flag='X'
   !   if ((iAge)<= m%nNW+1) flag='.'
   !   print'(A,A,I0,A)',flag,'iAge ',iAge,'      Root              Tip'
   !   print*,trim(label)//'x', m%dxdt_NW(1, 1, iAge,1), m%dxdt_NW(1, p%nSpan+1, iAge,1)
   !   print*,trim(label)//'y', m%dxdt_NW(2, 1, iAge,1), m%dxdt_NW(2, p%nSpan+1, iAge,1)
   !   print*,trim(label)//'z', m%dxdt_NW(3, 1, iAge,1), m%dxdt_NW(3, p%nSpan+1, iAge,1)
   !enddo
   !print'(A,I0)','dxdt FW <<<<<<<<<<<<<<<<<<<< nFW:',m%nFW
   !do iAge=1,p%nFWMax+1
   !   flag='X'
   !   if ((iAge)<= m%nFW+1) flag='.'
   !   print'(A,A,I0,A)',flag,'iAge ',iAge,'      Root              Tip'
   !   print*,trim(label)//'x', m%dxdt_FW(1, 1, iAge,1), m%dxdt_FW(1, FWnSpan+1, iAge,1)
   !   print*,trim(label)//'y', m%dxdt_FW(2, 1, iAge,1), m%dxdt_FW(2, FWnSpan+1, iAge,1)
   !   print*,trim(label)//'z', m%dxdt_FW(3, 1, iAge,1), m%dxdt_FW(3, FWnSpan+1, iAge,1)
   !enddo
endsubroutine


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
   type(FVW_MiscVarType),           intent(in   )              :: m       !< Initial misc/optimization variables
   integer(IntKi)          :: iP_start,iP_end   ! Current index of point, start and end of range

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
subroutine DistributeRequestedWind(V_wind, p, m)
   real(ReKi), dimension(:,:),      intent(in   ) :: V_wind  !< Position where wind is requested
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_MiscVarType),           intent(inout) :: m       !< Initial misc/optimization variables
   integer(IntKi)          :: iP_start,iP_end   ! Current index of point, start and end of range

   ! Using array reshaping to ensure a given near or far wake point is always at the same location in the array.
   ! NOTE: Maximum number of points are passed, whether they "exist" or not. 
   ! --- LL CP
   iP_start=1
   iP_end=p%nWings*p%nSpan
   m%Vwnd_LL(1:3,1:p%nSpan,1:p%nWings) = reshape( V_wind(1:3,iP_start:iP_end), (/ 3, p%nSpan, p%nWings /))
   ! --- NW points
   iP_start=iP_end+1
   iP_end=iP_start-1+(p%nSpan+1)*(p%nNWMax+1)*p%nWings
   m%Vwnd_NW(1:3,1:p%nSpan+1,1:p%nNWMax+1,1:p%nWings) = reshape( V_wind(1:3,iP_start:iP_end), (/ 3, p%nSpan+1, p%nNWMax+1, p%nWings/))
   ! --- FW points
   if (p%nFWMax>0) then
      iP_start=iP_end+1
      iP_end=iP_start-1+(FWnSpan+1)*(p%nFWMax+1)*p%nWings
      m%Vwnd_FW(1:3,1:FWnSpan+1,1:p%nFWMax+1,1:p%nWings) = reshape( V_wind(1:3,iP_start:iP_end), (/ 3, FWnSpan+1, p%nFWMax+1, p%nWings /))
   endif

end subroutine DistributeRequestedWind


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


subroutine PackPanelsToSegments(p, m, x, iDepthStart, bMirror, SegConnct, SegPoints, SegGamma, nSeg, nSegP)
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_MiscVarType),           intent(in   ) :: m       !< Initial misc/optimization variables
   type(FVW_ContinuousStateType),   intent(in   ) :: x       !< States
   integer(IntKi),                  intent(in   ) :: iDepthStart !< Index where we start packing for NW panels
   logical,                         intent(in   ) :: bMirror !< Mirror the vorticity wrt the ground
   integer(IntKi),dimension(:,:), intent(inout) :: SegConnct !< Segment connectivity
   real(ReKi),    dimension(:,:), intent(inout) :: SegPoints !< Segment Points
   real(ReKi),    dimension(:)  , intent(inout) :: SegGamma  !< Segment Circulation
   integer(IntKi), intent(out)                :: nSeg      !< Total number of segments after packing
   integer(IntKi), intent(out)                :: nSegP     !< Total number of segments points after packing
   ! Local
   integer(IntKi) :: iHeadC, iHeadP, nC, nCNW, nP, iW, iHeadC_bkp, i, iMirror
   logical        :: LastNWShed

   ! If the FW contains Shed vorticity, we include the last shed vorticity form the NW, orhtwerise, we don't!
   ! It's important not to include it, otherwise a strong vortex will be present there with no compensating vorticity from the FW
   LastNWShed = (p%FWShedVorticity ) .or. ((.not.p%FWShedVorticity) .and. (m%nNW<p%nNWMax))

   ! Counting total number of segments
   ! Returns nC, nP, nCNW, number of segments (without accounting for mirroring)
   call CountSegments(p, m%nNW, m%nFW, iDepthStart, nC, nP, nCNW)

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
            CALL LatticeToSegments(x%r_NW(1:3,:,1:m%nNW+1,iW), x%Gamma_NW(:,1:m%nNW,iW), iDepthStart, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC, .True., LastNWShed )
         enddo
      endif
      if (m%nFW>0) then
         iHeadC_bkp = iHeadC
         do iW=1,p%nWings
            CALL LatticeToSegments(x%r_FW(1:3,:,1:m%nFW+1,iW), x%Gamma_FW(:,1:m%nFW,iW), 1, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC , p%FWShedVorticity, p%FWShedVorticity)
         enddo
         SegConnct(3,iHeadC_bkp:) = SegConnct(3,iHeadC_bkp:) + m%nNW ! Increasing iDepth (or age) to account for NW
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
            call print_x_NW_FW(p,m,x,'pack')
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
!! NOTE: this should preferably be done at the "panel"/vortex sheet level
subroutine FVW_InitRegularization(p, m, ErrStat, ErrMsg)
   type(FVW_ParameterType),         intent(inout) :: p       !< Parameters
   type(FVW_MiscVarType),           intent(inout) :: m       !< Initial misc/optimization variables
   integer(IntKi),                  intent(  out) :: ErrStat    !< Error status of the operation
   character(*),                    intent(  out) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   ! Local variables
   real(ReKi) :: ds_min, ds_max, ds_mean !< min,max and mean of spanwise sections
   real(ReKi) :: c_min, c_max, c_mean !< min,max and mean of chord
   real(ReKi) :: d_min, d_max, d_mean !< min,max and mean of panel diagonal
   real(ReKi) :: RegParam
   real(ReKi) :: Span !< "Blade span"
   integer :: iW
   ErrStat = ErrID_None
   ErrMsg  = ""
   ! --- Compute min max and mean spanwise section lengths
   iW =1
   ds_min  = minval(m%s_ll(2:p%nSpan+1,iW)-m%s_ll(1:p%nSpan,iW))
   ds_max  = maxval(m%s_ll(2:p%nSpan+1,iW)-m%s_ll(1:p%nSpan,iW))
   ds_mean = sum(m%s_ll(2:p%nSpan+1,iW)-m%s_ll(1:p%nSpan,iW))/(p%nSpan+1)
   c_min  = minval(m%chord_LL(:,iW))
   c_max  = maxval(m%chord_LL(:,iW))
   c_mean = sum   (m%chord_LL(:,iW))/(p%nSpan+1)
   d_min  = minval(m%diag_LL(:,iW))
   d_max  = maxval(m%diag_LL(:,iW))
   d_mean = sum   (m%diag_LL(:,iW))/(p%nSpan+1)
   Span    = m%s_ll(p%nSpan+1,iW)-m%s_ll(1,iW)
   RegParam = ds_mean*2
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
   if (p%RegDeterMethod==idRegDeterAuto) then
      ! TODO this is beta
      print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print*,'!!! NOTE: using optmized wake regularization parameters is still a beta feature!'
      print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      p%WakeRegMethod      = idRegConstant
      p%RegFunction        = idRegVatistas
      p%WakeRegParam       = RegParam
      p%WingRegParam       = RegParam
      p%CoreSpreadEddyVisc = 100
      write(*,'(A)'   )   'The following regularization parameters will be used:'
      write(*,'(A,I0)'   )   'WakeRegMethod     : ', p%WakeRegMethod
      write(*,'(A,I0)'   )   'RegFunction       : ', p%RegFunction
      write(*,'(A,1F8.4)')   'WakeRegParam      : ', p%WakeRegParam
      write(*,'(A,1F8.4)')   'WingRegParam      : ', p%WingRegParam
      write(*,'(A,1F8.4)')   'CoreSpreadEddyVisc: ', p%CoreSpreadEddyVisc
   endif
   ! KEEP ME: potentially perform pre-computation here
   !if (p%WakeRegMethod==idRegConstant) then
   !else if (p%WakeRegMethod==idRegStretching) then
   !else if (p%WakeRegMethod==idRegAge) then
   !else
   !   ErrStat = ErrID_Fatal
   !   ErrMsg ='Regularization method not implemented'
   !endif
end subroutine FVW_InitRegularization


!> Set up regularization parameter based on diffusion method and regularization method
!! NOTE: this should preferably be done at the "panel"/vortex sheet level
subroutine WakeRegularization(p, x, m, SegConnct, SegPoints, SegGamma, SegEpsilon, ErrStat, ErrMsg)
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_ContinuousStateType),   intent(in   ) :: x       !< States
   type(FVW_MiscVarType),           intent(in   ) :: m       !< Initial misc/optimization variables
   integer(IntKi),dimension(:,:)  , intent(in   ) :: SegConnct  !< Segment connectivity
   real(ReKi),    dimension(:,:)  , intent(in   ) :: SegPoints  !< Segment Points
   real(ReKi),    dimension(:)    , intent(in   ) :: SegGamma   !< Segment Circulation
   real(ReKi),    dimension(:)    , intent(  out) :: SegEpsilon !< Segment regularization parameter
   integer(IntKi),                  intent(  out) :: ErrStat    !< Error status of the operation
   character(*),                    intent(  out) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   ! Local variables
   integer(IntKi) :: iSeg
   real(ReKi) :: time
   ErrStat = ErrID_None
   ErrMsg  = ""

   ! 
   if (p%WakeRegMethod==idRegConstant) then
      SegEpsilon=p%WakeRegParam 

   else if (p%WakeRegMethod==idRegStretching) then
      ! TODO
      ErrStat = ErrID_Fatal
      ErrMsg ='Regularization method not implemented'
      if (.false.) print*,m%nNW,x%r_NW(1,1,1,1),SegPoints(1,1),SegGamma(1) ! Needed in the future, Just to avoid unused var warning

   else if (p%WakeRegMethod==idRegAge) then
      do iSeg=1,size(SegEpsilon,1) ! loop on segments
         time = (SegConnct(3, iSeg)-1) * p%DTfvw ! column 3 contains "iDepth", or "iAge", from 1 to nSteps
         SegEpsilon(iSeg) = sqrt( 4._ReKi * CoreSpreadAlpha * p%CoreSpreadEddyVisc * p%KinVisc* time  + p%WakeRegParam**2 )
      enddo

   else
      ErrStat = ErrID_Fatal
      ErrMsg ='Regularization method not implemented'
   endif

end subroutine WakeRegularization


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
   integer(IntKi) :: iW, nSeg, nSegP, nCPs, iHeadP
   integer(IntKi) :: nFWEff  ! Number of farwake panels that are free at current tmie step
   logical        :: bMirror ! True if we mirror the vorticity wrt ground
   ! TODO new options
   integer(IntKi) :: RegFunctionPart
   integer(IntKi) :: nPart
   real(ReKi)     :: DistanceDirect ! Distance under which direct evaluation of the Biot-Savart should be done for tree
   type(T_Tree)   :: Tree
   real(ReKi), dimension(:,:), allocatable :: PartPoints !< Particle points
   real(ReKi), dimension(:,:), allocatable :: PartAlpha  !< Particle circulation
   real(ReKi), dimension(:)  , allocatable :: PartEpsilon !< Regularization parameter
   ErrStat= ErrID_None
   ErrMsg =''

   nFWEff = min(m%nFW, p%nFWFree)
   bMirror = p%ShearModel==idShearMirror ! Whether or not we mirror the vorticity wrt ground

   m%Vind_NW = -9999._ReKi !< Safety
   m%Vind_FW = -9999._ReKi !< Safety

   ! --- Packing all vortex elements into a list of segments
   ! NOTE: modifies m%Seg* 
   call PackPanelsToSegments(p, m, x, 1, bMirror, m%SegConnct, m%SegPoints, m%SegGamma, nSeg, nSegP)

   ! --- Setting up regularization SegEpsilon
   call WakeRegularization(p, x, m, m%SegConnct, m%SegPoints, m%SegGamma, m%SegEpsilon(1:nSeg), ErrStat, ErrMsg)

   ! --- Computing induced velocity
   call PackConvectingPoints()
   if (DEV_VERSION) then
      print'(A,I0,A,I0,A,I0)','Convection - nSeg:',nSeg,' - nSegP:',nSegP, ' - nCPs:',nCPs
   endif

   ! --- Converting to particles
   if ((p%VelocityMethod==idVelocityTree) .or. (p%VelocityMethod==idVelocityPart)) then
      iHeadP=1
      nPart = p%PartPerSegment * nSeg 
      allocate(PartPoints(3,nPart), PartAlpha(3,nPart), PartEpsilon(nPart))
      PartAlpha(:,:)  = -99999.99_ReKi
      PartPoints(:,:) = -99999.99_ReKi
      PartEpsilon(:)  = -99999.99_ReKi
      call SegmentsToPart(m%SegPoints, m%SegConnct, m%SegGamma, m%SegEpsilon, 1, nSeg, p%PartPerSegment, PartPoints, PartAlpha, PartEpsilon, iHeadP)
      if (p%RegFunction/=idRegNone) then
         RegFunctionPart = idRegExp ! TODO need to find a good equivalence and potentially adapt Epsilon in SegmentsToPart
      endif
      if (any(PartEpsilon(:)<-9999.99_ReKi)) then
         print*,'Error in Segment to part conversion'
         STOP
      endif
   endif

   ! --- Getting induced velocity
   m%Uind=0.0_ReKi ! very important due to side effects of ui_* methods
   if (p%VelocityMethod==idVelocityBasic) then
      call ui_seg( 1, nCPs, m%CPs, 1, nSeg, nSeg, nSegP, m%SegPoints, m%SegConnct, m%SegGamma, p%RegFunction, m%SegEpsilon, m%Uind)

   elseif (p%VelocityMethod==idVelocityTree) then

      DistanceDirect = 2*sum(PartEpsilon)/size(PartEpsilon) ! 2*mean(eps), below that distance eps has a strong effect
      call grow_tree(Tree, PartPoints, PartAlpha, RegFunctionPart, PartEpsilon, 0)
      !call print_tree(Tree)
      call ui_tree(Tree, m%CPs, 0, 1, nCPs, p%TreeBranchFactor, DistanceDirect, m%Uind, ErrStat, ErrMsg)
      call cut_tree(Tree)
      deallocate(PartPoints, PartAlpha, PartEpsilon)

   elseif (p%VelocityMethod==idVelocityPart) then
      call ui_part_nograd(m%CPs ,PartPoints, PartAlpha, RegFunctionPart, PartEpsilon, m%Uind, nCPs, nPart)
      deallocate(PartPoints, PartAlpha, PartEpsilon)
   endif
   call UnPackInducedVelocity()

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

end subroutine

!> Compute induced velocities from all vortex elements onto the lifting line control points
!! In : x%r_NW, x%r_FW, x%Gamma_NW, x%Gamma_FW
!! Out: m%Vind_LL
subroutine LiftingLineInducedVelocities(p, x, iDepthStart, m, ErrStat, ErrMsg)
   type(FVW_ParameterType),         intent(in   ) :: p       !< Parameters
   type(FVW_ContinuousStateType),   intent(in   ) :: x       !< States
   integer(IntKi),                  intent(in   ) :: iDepthStart !< Index where we start packing for NW panels
   type(FVW_MiscVarType),           intent(inout) :: m       !< Initial misc/optimization variables
   ! Local variables
   integer(IntKi) :: iW, nSeg, nSegP, nCPs, iHeadP
   real(ReKi),    dimension(:,:), allocatable :: CPs   !< ControlPoints
   real(ReKi),    dimension(:,:), allocatable :: Uind  !< Induced velocity
   integer(IntKi),              intent(  out) :: ErrStat    !< Error status of the operation
   character(*),                intent(  out) :: ErrMsg     !< Error message if ErrStat /= ErrID_None
   logical ::  bMirror 
   ErrStat = ErrID_None
   ErrMsg  = ""
   m%Vind_LL = -9999._ReKi !< Safety
   bMirror = p%ShearModel==idShearMirror ! Whether or not we mirror the vorticity wrt ground

   ! --- Packing all vortex elements into a list of segments
   call PackPanelsToSegments(p, m, x, iDepthStart, bMirror, m%SegConnct, m%SegPoints, m%SegGamma, nSeg, nSegP)

   ! --- Computing induced velocity
   if (nSegP==0) then
      nCPs=0
      m%Vind_LL = 0.0_ReKi
      if (DEV_VERSION) then
         print'(A,I0,A,I0,A,I0,A)','Induction -  nSeg:',nSeg,' - nSegP:',nSegP, ' - nCPs:',nCPs, ' -> No induction'
      endif
   else
      ! --- Setting up regularization
      call WakeRegularization(p, x, m, m%SegConnct(:,1:nSeg), m%SegPoints(:,1:nSegP), m%SegGamma(1:nSeg), m%SegEpsilon(1:nSeg), ErrStat, ErrMsg)

      nCPs=p%nWings * p%nSpan
      allocate(CPs (1:3,1:nCPs)) ! NOTE: here we do allocate CPs and Uind insteadof using Misc 
      allocate(Uind(1:3,1:nCPs)) !       The size is reasonably small, and m%Uind then stay filled with "rollup velocities" (for export)
      Uind=0.0_ReKi !< important due to side effects of ui_seg
      ! ---
      call PackLiftingLinePoints()
      if (DEV_VERSION) then
         print'(A,I0,A,I0,A,I0)','Induction -  nSeg:',nSeg,' - nSegP:',nSegP, ' - nCPs:',nCPs
      endif
      call ui_seg( 1, nCPs, CPs, 1, nSeg, nSeg, nSegP, m%SegPoints, m%SegConnct, m%SegGamma, p%RegFunction, m%SegEpsilon, Uind)
      call UnPackLiftingLineVelocities()

      deallocate(Uind)
      deallocate(CPs)
   endif
contains
   !> Pack all the control points
   subroutine PackLiftingLinePoints()
      iHeadP=1
      do iW=1,p%nWings
         CALL LatticeToPoints(m%CP_LL(1:3,:,iW:iW), 1, CPs, iHeadP)
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
         CALL VecToLattice(Uind, 1, m%Vind_LL(1:3,:,iW:iW), iHeadP)
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
   Re        = Chord * Vrel_norm / KinVisc / 1.0E6

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
   Re        = Chord * Vrel_norm / KinVisc / 1.0E6
end subroutine AlphaVrel_Generic


end module FVW_Subs
