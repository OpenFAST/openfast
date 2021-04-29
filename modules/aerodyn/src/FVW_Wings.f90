module FVW_Wings 

   use NWTC_Library
   use FVW_Types
   use FVW_Subs
   use AirFoilInfo,   only : AFI_ComputeAirfoilCoefs
 
   implicit none

contains
  
   !> Meshing function. Create a 1D mesh of size `n`, based on a `method` and an input vector `x`
   subroutine Meshing(method, x, n, y)
      ! Arguments declarations 
      character(len=*), intent(in)          :: method !< String defining the method used
      integer(IntKi), intent(in)            :: n      !< size of vector y
      real(ReKi), dimension(:),intent(in)   :: x      !< input vector of nodes or (/ min, max/)
      real(ReKi), dimension(:), intent(out) :: y      !< output vector with meshing values
      ! Variable declarations 
      real(ReKi), dimension(:),allocatable :: dx  !<  
      integer::jr
      y = 0.0_ReKi
      select case (method) !
      case ('middle') !
         allocate(dx(1:n))
         dx=diff(x) ! dx is the width of each panel 
         y(1:n)=x(1:n)+dx(1:n)/2._ReKi
         deallocate(dx)

      case ('fullcosineapprox') !
         ! x is assumed to be of size n+1
         if (n==1) then
            y(1)=(x(1)+x(2))/2._ReKi ! middle
            return
         else
            allocate(dx(1:n))
            dx=diff(x) ! dx is the width of each panel 
            y(1) = x(1)+(dx(1)  /(dx(1)  +dx(2)))*dx(1)
            y(n) = x(n)+(dx(n-1)/(dx(n-1)+dx(n)))*dx(n)
            do jr=2,n-1 
               y(jr)=x(jr)+0.25_ReKi*(dx(jr-1)/(dx(jr-1)+dx(jr)) + dx(jr)/(dx(jr)+dx(jr+1))+1 )*dx(jr)
            end do 
            deallocate(dx)
         endif
      end select 

   contains
      !> Compute: x(2:n)-x(1:n-1)
      function diff(d)
         real(ReKi),dimension(:),intent(in) ::d
         real(ReKi),dimension(size(d)-1) ::diff
         integer::i
         do i=1,size(d)-1
            diff(i)=d(i+1)-d(i)
         enddo
      end function
   end subroutine Meshing


   !----------------------------------------------------------------------------------------------------------------------------------
   !> Based on an input mesh, sets the following:
   !!  - s_LL       : Dimensionless spanwise coordinate of LL    
   !!  - s_CP_LL    : Dimensionless spanwise coordinate of LL CP 
   !!  - chord_LL   : chord on LL 
   !!  - chord_LL_CP: chord on LL cp  
   subroutine Wings_Panelling_Init(Meshes, p, m, ErrStat, ErrMsg )
      type(MeshType), dimension(:),    intent(in   )  :: Meshes         !< Wings mesh
      type(FVW_ParameterType),         intent(inout)  :: p              !< Parameters
      type(FVW_MiscVarType),           intent(inout)  :: m              !< Initial misc/optimization variables
      integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
      character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
      ! Local
      !integer(IntKi)          :: ErrStat2       ! temporary error status of the operation
      !character(ErrMsgLen)    :: ErrMsg2        ! temporary error message
      integer(IntKi) :: iW, iSpan
      real(ReKi), dimension(3) :: DP
      real(ReKi), dimension(:),allocatable :: s_in !< Dimensionless spanwise coordinate of input

      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

      ! --- Meshing
      do iW = 1,p%nWings
         if (allocated(s_in)) deallocate(s_in)
         allocate(s_in(1:Meshes(iW)%nNodes))
         ! --- Computing spanwise coordinate of input mesh normalized from 0 to 1
!Note: this info also exists in InitInp%zLocal 
         s_in(:) = -999
         s_in(1) = 0
         do iSpan = 2, Meshes(iW)%nNodes
            DP          = Meshes(iW)%Position(1:3, iSpan) - Meshes(iW)%Position(1:3, iSpan-1)
            s_in(iSpan) = s_in(iSpan-1) + TwoNorm(DP)
         enddo

         ! --- Setting up Lifting line variables based on input  and a "meshing" method (TODO)
         if (Meshes(iW)%nNodes /= p%nSpan+1) then
            ! TODO Possibly interpolate based on FVW meshing
            ! NOTE: p%chord is copied from the InitInput
            ErrMsg ='TODO different discretization InputMesh / vortex code'; ErrStat=ErrID_Fatal; return
         endif
         do iSpan = 1, p%nSpan+1
            p%s_LL    (iSpan, iW) = s_in(iSpan)
            p%chord_LL(iSpan, iW) = p%chord(iSpan,iW)
         enddo
         ! --- Control points spanwise location
         ! NOTE: we use the cos approximation of VanGarrel. For equispacing, it returns mid point
         !       otherwise, points are slightly closer to panels that are shorter
         !call Meshing('middle'           , p%s_LL(:,iW), p%nSpan, p%s_CP_LL(:,iW))
         call Meshing('fullcosineapprox' , p%s_LL(:,iW), p%nSpan, p%s_CP_LL(:,iW))
         call InterpArray(p%s_LL(:,iW), p%chord_LL(:,iW), p%s_CP_LL(:,iW), p%chord_CP_LL(:,iW))

         deallocate(s_in)
      enddo

   end subroutine Wings_Panelling_Init
   !----------------------------------------------------------------------------------------------------------------------------------
   !> Based on an input mesh, sets the following:
   !!  - LE      : Leading edge points                 (3 x nSpan+1 x nWings)
   !!  - TE      : Trailing edge points                (3 x nSpan+1 x nWings)
   !!  - CP_LL   : Coordinates of LL CP"              (3 x nSpan x nWings)
   !!  - Tang    : Unit Tangential vector on LL CP" -
   !!  - Norm    : Unit Normal vector on LL CP    " -
   !!  - Orth    : Unit Orthogonal vector on LL CP" -
   !!  - Vstr_LL : Structural velocity on LL CP" m/s
   subroutine Wings_Panelling(Meshes, p, m, ErrStat, ErrMsg )
      type(MeshType), dimension(:),    intent(in   )  :: Meshes         !< Wings mesh
      type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
      type(FVW_MiscVarType),           intent(inout)  :: m              !< Initial misc/optimization variables
      integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
      character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
      ! Local
      !integer(IntKi)          :: ErrStat2       ! temporary error status of the operation
      !character(ErrMsgLen)    :: ErrMsg2        ! temporary error message
      integer(IntKi) ::iSpan , iW
      real(ReKi), dimension(3) :: P_ref         ! Reference point of Input Mesh (e.g. AeroDynamic Center?)
      real(ReKi), dimension(3) :: DP_LE ! Distance between reference point and Leading edge
      real(ReKi), dimension(3) :: DP_TE ! Distance between reference point and trailing edge
      real(ReKi), dimension(3) :: P1,P2,P3,P4,P5,P7,P8,P6,P9,P10
      real(ReKi), dimension(3) :: DP1, DP2, DP3
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""
      ! --- Position of leading edge (LE) and trailing edge (TE)
      ! NOTE, this assumes one to one between InputMesh and FVW Mesh
      do iW = 1,p%nWings
         do iSpan = 1,p%nSpan+1
            P_ref = Meshes(iW)%Position(1:3, iSpan )+Meshes(iW)%TranslationDisp(1:3, iSpan)
            DP_LE(1:3) =  0.0
            DP_LE(1)   = -p%chord_LL(iSpan,iW)/4.
            DP_TE(1:3) =  0.0
            DP_TE(1)   = +3.*p%chord_LL(iSpan,iW)/4. 
            m%LE(1:3, iSpan, iW) = P_ref + DP_LE(1)*Meshes(iW)%Orientation(2,1:3,iSpan)
            m%TE(1:3, iSpan, iW) = P_ref + DP_TE(1)*Meshes(iW)%Orientation(2,1:3,iSpan)
         enddo         
      enddo
      ! --- Generic code below to compute normal/tangential vectors of a lifting line panel
      ! Notations follow vanGarrel [ECN-C--03-079, Development of a wind turbine aerodynamics simulation module,2003]
      !
      !   P4 -P10---P7------ P3
      !   |
      !   P8                 P6
      !   |        
      !   P1 -P9----P5------ P2
      !
      do iW = 1,p%nWings
         do iSpan = 1,p%nSpan
            P1                    = m%LE(:,iSpan  , iw)
            P4                    = m%LE(:,iSpan+1, iw)
            P3                    = m%TE(:,iSpan+1, iw)
            P2                    = m%TE(:,iSpan  , iw)
            P8                    = (P1+P4)/2
            P6                    = (P2+P3)/2
            P5                    = (P1+P2)/2
            P7                    = (P4+P3)/2
            P9                    = 0.75_ReKi*P1+0.25_ReKi*P2
            P10                   = 0.75_ReKi*P4+0.25_ReKi*P3
            DP1                   = P6-P8
            DP2                   = P10-P9
            DP3                   = P7-P5

            ! NOTE: The denominator below has seg-faulted with Intel Fortran in Release mode, possibly due to nuances in copmiler optimizations.
            !       This code was previously after the m%Norm calculations, but moving it up "fixes" the bug.
            m%Tang(1:3,iSpan,iW)  = (DP1)/TwoNorm(DP1)                       ! tangential unit vector, along chord

            m%Norm(1:3,iSpan,iW)  = cross_product(DP1,DP2)
            m%Norm(1:3,iSpan,iW)  = m%Norm(1:3,iSpan,iW)/TwoNorm(m%Norm(1:3,iSpan,iW))
            ! m%Tscoord(1:3,iSpan) = (DP3)/norm2(DP3)                      ! tangential unit vector, along span, follows ref line
            m%dl  (1:3,iSpan,iW)  = DP2
            m%Orth(1:3,iSpan,iW)  = cross_product(m%Norm(1:3,iSpan,iW),m%Tang(1:3,iSpan,iW)) ! orthogonal vector to N and T
            m%Area(iSpan, iW) = TwoNorm(cross_product(DP1,DP3))
            DP3 = P1-P3
            m%diag_LL(iSpan, iW) = TwoNorm(DP3)
         end do
      enddo
      ! --- Lifting Line/ Bound Circulation panel
      ! For now: goes from 1/4 chord to TE
      ! More panelling options may be considered in the future
      do iW = 1,p%nWings
         do iSpan = 1,p%nSpan+1
            m%r_LL(1:3,iSpan,1,iW)= m%TE(1:3,iSpan,iW)*0.25_ReKi+m%LE(1:3,iSpan,iW)*0.75_ReKi  ! 1/4 chord
            m%r_LL(1:3,iSpan,2,iW)= m%TE(1:3,iSpan,iW)                                         ! TE
         enddo
      enddo

      ! --- Position of control points CP_LL
      ! For now: placed exactly on the LL panel
      ! NOTE: separated from other loops just in case a special discretization is used
      do iW = 1,p%nWings
         call InterpArray(p%s_LL(:,iW), m%r_LL(1,:,1,iW), p%s_CP_LL(:,iW), m%CP_LL(1,:,iW))
         call InterpArray(p%s_LL(:,iW), m%r_LL(2,:,1,iW), p%s_CP_LL(:,iW), m%CP_LL(2,:,iW))
         call InterpArray(p%s_LL(:,iW), m%r_LL(3,:,1,iW), p%s_CP_LL(:,iW), m%CP_LL(3,:,iW))
      enddo

      ! --- Structural velocity on LL
      ! TODO: difference meshes in/LL
      do iW = 1,p%nWings
         call InterpArray(p%s_LL(:,iW), Meshes(iW)%TranslationVel(1,:) ,p%s_CP_LL(:,iW), m%Vstr_LL(1,:,iW))
         call InterpArray(p%s_LL(:,iW), Meshes(iW)%TranslationVel(2,:) ,p%s_CP_LL(:,iW), m%Vstr_LL(2,:,iW))
         call InterpArray(p%s_LL(:,iW), Meshes(iW)%TranslationVel(3,:) ,p%s_CP_LL(:,iW), m%Vstr_LL(3,:,iW))
      enddo
   end subroutine Wings_Panelling

   !----------------------------------------------------------------------------------------------------------------------------------
   !>
   subroutine Wings_ComputeCirculation(t, Gamma_LL, Gamma_LL_prev, u, p, x, m, AFInfo, ErrStat, ErrMsg, iLabel)
      real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
      real(ReKi), dimension(:,:),      intent(inout)  :: Gamma_LL       !< Circulation on all the lifting lines
      real(ReKi), dimension(:,:),      intent(in   )  :: Gamma_LL_prev  !< Previous/Guessed circulation
      type(FVW_InputType),             intent(in   )  :: u              !< Parameters
      type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
      type(FVW_ContinuousStateType),   intent(in   )  :: x              !< Parameters
      type(FVW_MiscVarType),           intent(inout)  :: m              !< Initial misc/optimization variables
      type(AFI_ParameterType),         intent(in   )  :: AFInfo(:)      !< The airfoil parameter data
      integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
      character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
      integer(IntKi), intent(in) :: iLabel
      ! Local
      integer(IntKi) :: iW
      real(DbKi) :: s, RealAxis
      real(ReKi) :: GammaScale
      ErrStat = ErrID_None
      ErrMsg  = ""

      if (t<p%FullCirculationStart) then
         ! The circulation is ramped up progressively, starting from 0 
         if (t<=0.0_DbKi) then
            GammaScale=0.0_ReKi
         else
            s=t/p%FullCirculationStart
            ! If we have at least 10 points we use a smooth Heavyside, otherwise we use a simple linear scaling
            if (p%FullCirculationStart/p%DTfvw >= 9) then
               ! Smooth approximations of the Heavyside function
               ! Example 1: 1/2 (1+2/pi arctan(k x) )  x \in ]-infty,+infty [
               ! Example 2: 1/(1+exp(k x) )            x \in ]-infty,+infty [
               ! Example 3: sin(pi/2*s)**2             s \in [0,1]
               RealAxis = (1-2*s)/(s*(s-1._DbKi)-1.0e-02_DbKi) ! Scaling from 0-1 to real axis
               GammaScale = 1._ReKi- 1._ReKi/(1._ReKi+exp(real(RealAxis,ReKi)))
            else
               GammaScale = s  ! Using a linear scaling
            endif
         endif
      else
         GammaScale=1.0_ReKi
      endif

      if (p%CirculationMethod==idCircPrescribed) then 
         do iW = 1, p%nWings !Loop over lifting lines
            Gamma_LL(1:p%nSpan,iW) = p%PrescribedCirculation(1:p%nSpan)
         enddo
         m%Vind_LL=-9999._ReKi !< Safety 
         m%Vtot_LL=-9999._ReKi !< Safety 

      else if (p%CirculationMethod==idCircPolarData) then 
         ! ---  Solve for circulation using polar data
         CALL Wings_ComputeCirculationPolarData(Gamma_LL, Gamma_LL_prev, p, x, m, AFInfo, GammaScale, ErrStat, ErrMsg, iLabel)

      else if (p%CirculationMethod==idCircNoFlowThrough) then 
         ! ---  Solve for circulation using the no-flow through condition
         ErrMsg='Circulation method nor implemented'; ErrStat=ErrID_Fatal; return ! should never happen
      else
         ErrMsg='Circulation method nor implemented'; ErrStat=ErrID_Fatal; return ! should never happen
      endif

      ! Scale circulation (for initial transient)
      Gamma_LL = Gamma_LL * GammaScale


   endsubroutine Wings_ComputeCirculation

   !----------------------------------------------------------------------------------------------------------------------------------
   !>
   subroutine Wings_ComputeCirculationPolarData(Gamma_LL, Gamma_LL_prev, p, x, m, AFInfo, GammaScale, ErrStat, ErrMsg, iLabel)
      real(ReKi), dimension(:,:),      intent(inout)  :: Gamma_LL       !< Circulation on all the lifting lines
      real(ReKi), dimension(:,:),      intent(in   )  :: Gamma_LL_prev  !< Previous/Guessed circulation
      type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
      type(FVW_ContinuousStateType),   intent(in   )  :: x              !< Parameters
      type(FVW_MiscVarType),           intent(inout)  :: m              !< Initial misc/optimization variables
      real(ReKi),                      intent(in   )  :: GammaScale     !< Scaling factor used at init
      type(AFI_ParameterType),         intent(in   )  :: AFInfo(:)      !< The airfoil parameter data
      integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
      character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
      integer(IntKi), intent(in) :: iLabel
      ! Local
      real(ReKi), dimension(:,:), allocatable :: DGamma        !< 
      real(ReKi), dimension(:,:), allocatable :: GammaLastIter !< 
      logical                                 :: bConverged    !< 
      integer(IntKi)                          :: iIter         !< iteration step number
      real(ReKi)                              :: MeanGamma
      real(ReKi), dimension(:,:,:), allocatable :: Vcst !< Part of the velocity that is constant 
      real(ReKi), dimension(:,:,:), allocatable :: Vvar !< Part of the velocity that is varies due to the solve
      integer(IntKi) :: iW, iSpan, iDepth, iWCP, nCPs
      real(ReKi), dimension(3) :: P1, P2, P3, P4
      real(ReKi) :: Gamm
      ! Error handling
      integer(IntKi)           :: ErrStat2
      character(ErrMsgLen)     :: ErrMsg2

      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

      !print*,'Parameters for circulation solv: ',p%CircSolvConvCrit ,p%CircSolvRelaxation ,p%CircSolvMaxIter   

      allocate(DGamma       (1:p%nSpan,1:p%nWings))
      allocate(GammaLastIter(1:p%nSpan,1:p%nWings))

      ! --- Last iteration circulation
      if (m%FirstCall) then
         ! We find a guess by looking simply at the Wind and Elasticity velocity
         m%Vtot_ll = m%Vwnd_LL - m%Vstr_ll
         call CirculationFromPolarData(GammaLastIter, p, m, AFInfo,ErrStat2,ErrMsg2);  if(Failed()) return;
      else
         ! NOTE: we need to inverse the scaling to speed up the convergence
         if (.not. EqualRealNos(GammaScale, 0.0_ReKi)) then
            GammaLastIter(1:p%nSpan,1:p%nWings) = Gamma_LL_prev(1:p%nSpan,1:p%nWings) / GammaScale 
         else
            GammaLastIter(1:p%nSpan,1:p%nWings) = Gamma_LL_prev(1:p%nSpan,1:p%nWings)
         endif
      endif

      if (any(x%r_NW(1,:,1:m%nNW+1,:)<-999)) then
         ErrMsg='Wings_ComputeCirculationPolarData: Problem in input NW points'; ErrStat=ErrID_Fatal; return
      endif


      ! --- Setting up Vcst: part of the velocity that is constant withing the iteration loop
      !   Vrel_ll_cst = U_u0 - U_body 
      call AllocAry(Vvar,  3, p%nSpan, p%nWings, 'Vvar',  ErrStat2, ErrMsg2);  if(Failed()) return;
      call AllocAry(Vcst,  3, p%nSpan, p%nWings, 'Vcst',  ErrStat2, ErrMsg2);  if(Failed()) return;

      ! Set m%Vind_LL Induced velocity from Known wake only (after iNWStart+1)
      call LiftingLineInducedVelocities(m%CP_LL, p, x, iNWStart+1, m, m%Vind_LL, ErrStat2, ErrMsg2);  if(Failed()) return;

      Vcst = m%Vind_LL + m%Vwnd_LL - m%Vstr_ll

      if (any(m%Vind_LL(1:3,:,:)<-99)) then
         ErrMsg='Wings_ComputeCirculationPolarData: Problem in induced velocity on LL points'; ErrStat=ErrID_Fatal; return
      endif
      if (any(m%Vwnd_LL(1:3,:,:)<-99)) then
         ErrMsg='Wings_ComputeCirculationPolarData: Problem in wind velocity on LL points'; ErrStat=ErrID_Fatal; return
      endif

      ! --- Convergence loop until near wake gives induction coherent with circulation
      bConverged=.false.
      iIter=0
      do while (.not.(bConverged) .and. iIter<p%CircSolvMaxIter) 
          !print*,'------- ITERATION',iIter
          ! --- The induced velocity from the profiles is different at each iteration:
          Vvar=0 
          nCPs=p%nSpan
          do iW=1,p%nWings
             do iSpan=1,p%nSpan
                do iDepth=1,iNWStart ! Two first panels
                   P1=x%r_NW(1:3,iSpan  ,iDepth  ,iW)
                   P2=x%r_NW(1:3,iSpan+1,iDepth  ,iW)
                   P3=x%r_NW(1:3,iSpan+1,iDepth+1,iW)
                   P4=x%r_NW(1:3,iSpan  ,iDepth+1,iW)
                   Gamm=GammaLastIter(iSpan, iW)
                   do iWCP=1,p%nWings
                      call ui_quad_n1(m%CP_LL(1:3,1:p%nSpan,iWCP), nCPs, P1, P2, P3, P4, Gamm, p%RegFunction, x%Eps_NW(1,iSpan,iDepth,iW), Vvar(1:3,1:p%nSpan,iWCP))
                   enddo
                enddo
             enddo
          enddo
          ! Total velocity on the lifting line
          m%Vtot_ll = Vcst + Vvar
          ! --- Computing circulation based on Vtot_LL
          call CirculationFromPolarData(Gamma_LL, p, m, AFInfo,ErrStat2,ErrMsg2);  if(Failed()) return;

          ! --------------------------------------------- 
          ! Differences between iterations and relaxation
          DGamma=Gamma_LL-GammaLastIter 
          GammaLastIter=GammaLastIter+p%CircSolvRelaxation*DGamma

          iIter=iIter+1
          MeanGamma  = sum(abs(GammaLastIter))/(p%nWings*p%nSpan)
          bConverged = maxval(abs(DGamma))/(MeanGamma)<p%CircSolvConvCrit

      end do ! convergence loop
      if (iIter==p%CircSolvMaxIter) then
         if (DEV_VERSION) then
            print'(A,I0,A,I0,A)','Circulation solve, call ',iLabel,', done after ........................ nIter: ', iIter, ' <<< Max reached'
         endif
         Gamma_LL=GammaLastIter ! returning relaxed value if not converged
      else
         if (DEV_VERSION) then
            print'(A,I0,A,I0)','Circulation solve, call ',iLabel,', done after ........................ nIter: ', iIter
         endif
         ! We return Gamma_LL
      endif
      ! KEEP ME: --- ADP: removed m%iStep in favor of m%VTKstep
      !iW=1
      !call Output_Gamma(m%CP_ll(1:3,:,iW), Gamma_LL(:,iW), iW, m%iStep, iLabel, iIter)
      !call print_mean_3d( m%Vwnd_LL(:,:,:), 'Mean wind    vel. LL (cst)')
      !call print_mean_3d( m%Vstr_LL(:,:,:), 'Mean struct  vel. LL (cst)')
      !call print_mean_3d( m%Vind_LL(:,:,:), 'Mean induced vel. LL (cst)')
      !call print_mean_3d( Vvar(:,:,:)     , 'Mean induced vel. LL (var)')
      !call print_mean_3d( Vvar+m%Vind_LL(:,:,:), 'Mean induced vel. LL (tot)')
      !call print_mean_3d( m%Vtot_LL(:,:,:), 'Mean relativevel. LL (tot)')
      !print*,'m%Vind_LL',m%Vind_LL(1,:,:)
      !print*,'m%Vwnd_LL',m%Vwnd_LL(1,:,:)
      !print*,'m%Vcst_LL',Vcst(1,:,:)
      !print*,'Gamm: ',Gamma_LL(1, 1), Gamma_LL(p%nSpan,1)
      m%Vind_LL=-9999._ReKi !< Safety (the induction above was not the true one)
      m%Vtot_LL=-9999._ReKi !< Safety 
      call CleanUp()
   contains

      logical function Failed()
         call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'Wings_ComputeCirculationPolarData')
         Failed =  ErrStat >= AbortErrLev
         if (Failed) call CleanUp()
      end function Failed
      subroutine CleanUp()
         if(allocated(DGamma       ))       deallocate(DGamma       )
         if(allocated(GammaLastIter))       deallocate(GammaLastIter)
         if(allocated(Vcst))                deallocate(Vcst)
         if(allocated(Vvar))                deallocate(Vvar)
      end subroutine
   end subroutine Wings_ComputeCirculationPolarData


   !>  Compute circulation based on polar data
   !! Uses m%Vtot_ll to compute Gamma_ll
   subroutine CirculationFromPolarData(Gamma_LL, p, m, AFInfo, ErrStat, ErrMsg)
      real(ReKi), dimension(:,:),      intent(inout)  :: Gamma_LL       !< Circulation on all the lifting lines
      type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
      type(FVW_MiscVarType),           intent(inout)  :: m              !< Initial misc/optimization variables
      type(AFI_ParameterType),         intent(in   )  :: AFInfo(:)      !< The airfoil parameter data
      integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
      character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None

      ! Local
      integer(IntKi) :: iW, iCP  !< Index on wings and spanwise control points
      real(ReKi), dimension(3) :: N, Tc      !<  Normal and Tangent vector
      real(ReKi), dimension(3) :: Vrel, Vrel_orth, Vjouk, Vjouk_orth
      real(ReKi)               :: Vrel_orth_norm, Vjouk_orth_norm, Vrel_norm
      real(ReKi)               :: alpha, Re, Cl, Cd, Cm
      type(AFI_OutputType)     :: AFI_interp
      integer(IntKi)           :: ErrStat2
      character(ErrMsgLen)     :: ErrMsg2

      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

      do iW=1,p%nWings 
         do icp=1,p%nSpan
            ! Aliases to shorten notations
            N    = m%Norm(1:3, icp, iW) 
            Tc   = m%Tang(1:3, icp, iW)
            Vrel = m%Vtot_LL(1:3,icp,iW)
            ! "Orth": cross sectional plane of the lifting line 
            Vrel_orth(1:3)  = dot_product(Vrel,N)*N + dot_product(Vrel,Tc)*Tc
            Vrel_orth_norm  = TwoNorm(Vrel_orth(1:3))
            Vjouk(1:3)      = cross_product(Vrel,m%dl(1:3,icp,iW))
            Vjouk_orth(1:3) = dot_product(Vjouk,N)*N + dot_product(Vjouk,Tc)*Tc
            Vjouk_orth_norm = TwoNorm(Vjouk_orth)
            Vrel_norm = TwoNorm(Vrel)

            alpha = atan2(dot_product(Vrel,N) , dot_product(Vrel,Tc) ) ! [rad]  
            Re = p%Chord(icp,iW) * Vrel_norm  / p%KinVisc  ! Reynolds number (not in Million)

            !if (p%CircSolvPolar==idPolarAeroDyn) then
               ! compute steady Airfoil Coefs      ! NOTE: UserProp set to 0.0_ReKi (no idea what it does).  Also, note this assumes airfoils at nodes.
!TODO: AFindx is on the nodes, not control points.
            call AFI_ComputeAirfoilCoefs( alpha, Re, 0.0_ReKi,  AFInfo(p%AFindx(icp,iW)), AFI_interp, ErrStat2, ErrMsg2 ); if(Failed()) return;
            Cl = AFI_interp%Cl
            Cd = AFI_interp%Cd
            Cm = AFI_interp%Cm
            ! Simple method:
            !    Gamma_LL=(0.5 * Cl * Vrel_orth_norm*chord)
            ! VanGarrel's method:
            Gamma_LL(icp,iW) =(0.5_ReKi * Cl * Vrel_orth_norm**2*m%Area(icp,iW)/(Vjouk_orth_norm))
            ! Convenient storage
            m%alpha_LL(icp, iW) = alpha ! [rad]
            m%Vreln_LL(icp, iW) = Vrel_norm
         enddo
      enddo
   contains
      logical function Failed()
         character(25)              :: NodeText
         if (ErrStat2 /= ErrID_None) then
            NodeText = '(node '//trim(num2lstr(icp))//', blade '//trim(num2lstr(iW))//')'
            call SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'CirculationFromPolarData'//trim(NodeText))
         end if
         Failed =  ErrStat >= AbortErrLev
      end function Failed
   end subroutine CirculationFromPolarData

end module FVW_Wings
