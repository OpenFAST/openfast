module FVW_Wings 

   use NWTC_Library
   use FVW_Types
   use FVW_Subs

   implicit none

contains

   !----------------------------------------------------------------------------------------------------------------------------------
   !> Based on an input mesh, sets the following:
   !!  - s_LL       : Dimensionless spanwise coordinate of LL    
   !!  - s_CP_LL    : Dimensionless spanwise coordinate of LL CP 
   !!  - chord_LL   : chord on LL 
   !!  - chord_LL_CP: chord on LL cp  
   subroutine Wings_Panelling_Init(Meshes, r, chord, p, m, ErrStat, ErrMsg )
      type(MeshType), dimension(:),    intent(in   )  :: Meshes         !< Wings mesh
      real(ReKi), dimension(:),        intent(in   )  :: r              !< 
      real(ReKi), dimension(:),        intent(in   )  :: chord          !< 
      type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
      type(FVW_MiscVarType),           intent(inout)  :: m              !< Initial misc/optimization variables
      integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
      character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
      ! Local
      integer(IntKi)          :: ErrStat2       ! temporary error status of the operation
      character(ErrMsgLen)    :: ErrMsg2        ! temporary error message
      integer(IntKi) :: iW, iSpan
      real(ReKi), dimension(3) :: First, Last, P1, P2, Pmid, DP
      real(ReKi) :: ds, length
      real(ReKi) :: c1,c2
      real(ReKi), dimension(:),allocatable :: s_in !< Dimensionless spanwise coordinate of input

      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

      ! --- Meshing
      do iW = 1,p%nWings
         if (allocated(s_in)) deallocate(s_in)
         allocate(s_in(1:Meshes(iW)%nNodes))
         ! --- Computing spanwise coordinate of input mesh normalized from 0 to 1
         s_in(:) = -999
         First  = Meshes(iW)%Position(1:3,1        )
         Last   = Meshes(iW)%Position(1:3,p%nSpan+1)
         DP     = Last - First
         length = TwoNorm(DP)
         do iSpan = 1, Meshes(iW)%nNodes
            P1          = Meshes(iW)%Position(1:3, iSpan  )
            DP          = P1-First
            s_in(iSpan) = TwoNorm(DP) / length
         enddo

         ! --- Setting up Lifting line variables based on input  and a "meshing" method (TODO)
         if (Meshes(iW)%nNodes /= p%nSpan+1) then
            ! TODO Possibly interpolate based on FVW meshing
            print*,'TODO different discretization InputMesh / vortex code'
            STOP
         endif
         print*,'Input mesh size',Meshes(iW)%nNodes,' Number of vortex element', p%nSpan
         do iSpan = 1, p%nSpan+1
            m%s_LL    (iSpan, iW) = s_in(iSpan)
            m%chord_LL(iSpan, iW) = chord(iSpan)
         enddo
         ! --- Control points
         ! TODO possibly Control points are not exactly at the middle depending on "meshing" method
         do iSpan = 1, p%nSpan
            m%s_CP_LL (iSpan, iW) = (m%s_LL    (iSpan,iW)+ m%s_LL    (iSpan+1,iW))/2
            m%chord_LL(iSpan, iW) = (m%chord_LL(iSpan,iW)+ m%chord_LL(iSpan+1,iW))/2
         enddo
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
      use Interpolation, only: interp_lin
      type(MeshType), dimension(:),    intent(in   )  :: Meshes         !< Wings mesh
      type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
      type(FVW_MiscVarType),           intent(inout)  :: m              !< Initial misc/optimization variables
      integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
      character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
      ! Local
      integer(IntKi)          :: ErrStat2       ! temporary error status of the operation
      character(ErrMsgLen)    :: ErrMsg2        ! temporary error message
      integer(IntKi) ::iSpan , iW
      real(ReKi), dimension(3) :: P_ref         ! Reference point of Input Mesh (e.g. AeroDynamic Center?)
      real(ReKi), dimension(3) :: DP_LE ! Distance between reference point and Leading edge
      real(ReKi), dimension(3) :: DP_TE ! Distance between reference point and trailing edge
      real(ReKi), dimension(3) :: P1,P2,P3,P4,P5,P7,P8,P6,P9,P10
      real(ReKi), dimension(3) :: DP1, DP2, DP3
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""
      ! --- Position of leading edge and trailing edge
      ! TODO, this assumes one to one between InputMesh and FVW Mesh
      !
      do iW = 1,p%nWings
         do iSpan = 1,p%nSpan+1
            P_ref = Meshes(iW)%Position(1:3, iSpan )  
            DP_LE(1:3) =  0.0
            DP_LE(1)   = -m%chord_LL(iSpan,iW)/2  ! TODO TODO TODO Use orientation and might not be c/2
            DP_TE(1:3) =  0.0
            DP_TE(1)   = +m%chord_LL(iSpan,iW)/2  ! TODO TODO TODO Use orientation and might not be c/2
            m%LE(1:3, iSpan, iW) = P_ref + DP_LE
            m%TE(1:3, iSpan, iW) = P_ref + DP_TE
         enddo         
      enddo
      ! --- Generic code below to compute normal/tangential vectors of a lifting line panel
      ! Notations follow vanGarrel [TODO REF]
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
            m%Norm(1:3,iSpan,iW)  = cross_product(DP1,DP2)
            m%Norm(1:3,iSpan,iW)  = m%Norm(1:3,iSpan,iW)/norm2(m%Norm(1:3,iSpan,iW))
            m%Tang(1:3,iSpan,iW)  = (DP1)/norm2(DP1)                       ! tangential unit vector, along chord
            ! m%Tscoord(1:3,iSpan) = (DP3)/norm2(DP3)                      ! tangential unit vector, along span, follows ref line
            m%dl  (1:3,iSpan,iW)  = DP2
            m%Orth(1:3,iSpan,iW)  = cross_product(m%Norm(1:3,iSpan,iW),m%Tang(1:3,iSpan,iW)) ! orthogonal vector to N and T
            m%Area(iSpan, iW) = norm2(cross_product(DP1,DP3));
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
      ! For now: placed at the "chordwise" middle of the LL panel
      ! NOTE: separated from other loops just in case a special discretization is used
      do iW = 1,p%nWings
         call interp_lin(m%s_LL(:,iW), m%r_LL(1,:,1,iW)*0.5_ReKi+m%r_LL(1,:,2,iW)*0.5_ReKi ,m%s_CP_LL(:,iW), m%CP_LL(1,:,iW))
         call interp_lin(m%s_LL(:,iW), m%r_LL(2,:,1,iW)*0.5_ReKi+m%r_LL(2,:,2,iW)*0.5_ReKi ,m%s_CP_LL(:,iW), m%CP_LL(2,:,iW))
         call interp_lin(m%s_LL(:,iW), m%r_LL(3,:,1,iW)*0.5_ReKi+m%r_LL(3,:,2,iW)*0.5_ReKi ,m%s_CP_LL(:,iW), m%CP_LL(3,:,iW))
      enddo

   end subroutine Wings_Panelling




   !----------------------------------------------------------------------------------------------------------------------------------
   !>
   subroutine Wings_ComputeCirculation(t, Gamma_LL, Gamma_LL_prev, u, p, x, m, ErrStat, ErrMsg)
      real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
      real(ReKi), dimension(:,:),      intent(inout)  :: Gamma_LL       !< Circulation on all the lifting lines
      real(ReKi), dimension(:,:),      intent(in   )  :: Gamma_LL_prev  !< Previous/Guessed circulation
      type(FVW_InputType),             intent(in   )  :: u              !< Parameters
      type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
      type(FVW_ContinuousStateType),   intent(in   )  :: x              !< Parameters
      type(FVW_MiscVarType),           intent(in   )  :: m              !< Initial misc/optimization variables
      integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
      character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
      ! Local
      integer(IntKi) :: iW
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

      if (p%CirculationMethod==idCircPrescribed) then 
         !print*,'>>>Prescribing circulation'
         do iW = 1, p%nWings !Loop over lifting lines
            Gamma_LL(1:p%nSpan,iW) = p%PrescribedCirculation(1:p%nSpan)
         enddo

      else if (p%CirculationMethod==idCircPolarData) then 
         ! ---  Solve for circulation using polar data
         ! TODO
         print*,'>>>Circulation solving with polar data'
         CALL Wings_ComputeCirculationPolarData(t, Gamma_LL, Gamma_LL_prev, u, p, x, m, ErrStat, ErrMsg)

      else if (p%CirculationMethod==idCircNoFlowThrough) then 
         ! ---  Solve for circulation using the no-flow through condition
         ! TODO
         print*,'Circulation method nor implemented', p%CirculationMethod
         STOP
      else
         print*,'Circulation method nor implemented', p%CirculationMethod ! Will never happen
         STOP
      endif

      if (t<p%FullCirculationStart) then
         ! The circulation is ramped up progressively, starting from 0 
         ! TODO use a smooth approximation of HeavySide function instead of linear
         print*,'Slow start'
         Gamma_LL = (t/p%FullCirculationStart)*Gamma_LL
      endif

   endsubroutine Wings_ComputeCirculation

   !----------------------------------------------------------------------------------------------------------------------------------
   !>
   subroutine Wings_ComputeCirculationPolarData(t, Gamma_LL, Gamma_LL_prev, u, p, x, m, ErrStat, ErrMsg)
      real(DbKi),                      intent(in   )  :: t           !< Current simulation time in seconds
      real(ReKi), dimension(:,:),      intent(inout)  :: Gamma_LL       !< Circulation on all the lifting lines
      real(ReKi), dimension(:,:),      intent(in   )  :: Gamma_LL_prev  !< Previous/Guessed circulation
      type(FVW_InputType),             intent(in   )  :: u              !< Parameters
      type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
      type(FVW_ContinuousStateType),   intent(in   )  :: x              !< Parameters
      type(FVW_MiscVarType),           intent(in   )  :: m              !< Initial misc/optimization variables
      integer(IntKi),                  intent(  out)  :: ErrStat        !< Error status of the operation
      character(*),                    intent(  out)  :: ErrMsg         !< Error message if ErrStat /= ErrID_None
      ! Local
      real(ReKi), dimension(:,:), allocatable :: DGamma        !< 
      real(ReKi), dimension(:,:), allocatable :: GammaLastIter !< 
      logical                                 :: bConverged    !< 
      integer(IntKi)                          :: iIter         !< iteration step number
      real(ReKi)                              :: MeanGamma
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""

      print*,'Parameters for circulation solv: ',p%CircSolvConvCrit ,p%CircSolvRelaxation ,p%CircSolvMaxIter   

      allocate(DGamma       (1:p%nSPan,1:p%nWings))
      allocate(GammaLastIter(1:p%nSPan,1:p%nWings))
      !
      GammaLastIter = Gamma_LL_prev

      ! Building Vrel_cst This part do not change wihtin the iteration loop
      ! Remember: uiu0 contains U0 and vorticity (free and prescribed)
      !do icp=1,SW%ncp_ll_tot 
      ! TODO
      !    SW%Vrel_ll_cst(1:3,icp)  = SW%U_uiu0(1:3,icp) - SW%U_body(1:3,icp) + SW%U_solv(1:3,icp)
      !end do 

      ! --- Convergence loop until near wake gives induction coherent with circulation
      bConverged=.false.
      iIter=0
      do while (.not.(bConverged) .and. iIter<p%CircSolvMaxIter) 
          !! Setting up wings and nw panels intensities
          !if (pAlgo%EmissionMethod==idEMethCreateTEBEforeSolving) then 
          !    ! if CreateTEBeforeSolving, then need to get influece of the TE panel 
          !    SW%GammaPanels(1:SW%ncp_ll_tot)                 = SW%GammaBoundLast
          !    SW%GammaPanels(1+SW%ncp_ll_tot:2*SW%ncp_ll_tot) = SW%GammaBoundLast ! The NW panels
          !else
          !    SW%GammaPanels=SW%GammaBoundLast
          !end if 
          !! The induced velocity from the profiles is different at each iteration:
          !call fUi_VoRing_IR(Panl,SW%CP_ll,SW%ncp_ll_tot,Panl%IUseProfilesWithNW,SW%GammaPanels,SW%U_prof)  ! no side effects
          !! Building Vrel ( taking advantage of the computation of the constant part)
          !do icp=1,SW%ncp_ll_tot 
          !    m%Vtot_ll(1:3,icp) = SW%Vrel_ll_cst(1:3,icp) + SW%U_prof(1:3,icp)
          !end do 
          ! --- Computing circulation based on Vtot_LL
          call CirculationFromPolarData(Gamma_LL, p, m)

          ! --------------------------------------------- 
          ! Differences between iterations and relaxation
          DGamma=Gamma_LL-GammaLastIter 
          GammaLastIter=GammaLastIter+p%CircSolvRelaxation*DGamma

          iIter=iIter+1
          MeanGamma  = sum(abs(GammaLastIter))/(p%nWings*p%nSpan)
          bConverged = maxval(abs(DGamma))/(MeanGamma)<p%CircSolvConvCrit
      end do ! convergence loop
      if (iIter==p%CircSolvMaxIter) then
          print*,'Maximum number of iterations reached: ',iIter
          Gamma_LL=GammaLastIter ! returning relaxed value if not converged
       else
          print*,'Circulation solve done after:', iIter,' iterations'
          ! We return Gamma_LL
       endif

      deallocate(DGamma       )
      deallocate(GammaLastIter)
      STOP

   end subroutine


   !>  Compute circulation based on polar data
   !! Uses m%Vtot_ll to compute Gamma_ll
   subroutine CirculationFromPolarData(Gamma_LL, p, m)
      real(ReKi), dimension(:,:),      intent(inout)  :: Gamma_LL       !< Circulation on all the lifting lines
      type(FVW_ParameterType),         intent(in   )  :: p              !< Parameters
      type(FVW_MiscVarType),           intent(in   )  :: m              !< Initial misc/optimization variables
      ! Local
      integer(IntKi) :: iW, iCP  !< Index on wings and spanwise control points
      real(ReKi), dimension(3) :: N, Tc      !<  Normal and Tangent vector
      real(ReKi), dimension(3) :: Vrel, Vrel_orth, Vjouk, Vjouk_orth
      real(ReKi)               :: Vrel_orth_norm, Vjouk_orth_norm
      real(ReKi)               :: alpha, Re, Cl

      do iW=1,p%nWings 
         do icp=1,p%nSpan
            ! Aliases to shorten notations
            N    = m%Norm(1:3, icp, iW) 
            Tc   = m%Tang(1:3, icp, iW)
            Vrel = m%Vtot_LL(1:3,icp,iW)
            ! "Orth": cross sectional plane of the lifting line 
            Vrel_orth(1:3)  = dot_product(Vrel,N)*N + dot_product(Vrel,Tc)*Tc
            Vrel_orth_norm  = norm2(Vrel_orth(1:3))
            Vjouk(1:3)      = cross_product(Vrel,m%dl(1:3,icp,iW))
            Vjouk_orth(1:3) = dot_product(Vjouk,N)*N + dot_product(Vjouk,Tc)*Tc
            Vjouk_orth_norm = norm2(Vjouk_orth)

            alpha = atan2(dot_product(Vrel,N) , dot_product(Vrel,Tc) ) ! [rad]  
            !Re    = LL%Vrel_orth_norm(icp)*LL%chord(icp)/KinVisc/(1.E6_MK) ! TODO TODO TODO KinVisc

            if (p%PrescribedPolar==idPolarAeroDyn) then
               print*,'TODO TODO TODO Get Cl, Cd, Cm from alpha, Re and AirfoilInfo'
               STOP
            else if (p%PrescribedPolar==idPolar2PiAlpha) then
               Cl=TwoPi*alpha
            else if (p%PrescribedPolar==idPolar2PiSinAlpha) then
               Cl=TwoPi*sin(alpha)
            else
               print*,'Unknown PrescribedPolar value'
               STOP
            endif
            ! Simple method:
            !    Gamma_LL=(0.5 * Cl * Vrel_orth_norm*chord)
            ! VanGarrel's method:
            Gamma_LL(icp,iW) =(0.5_ReKi * Cl * Vrel_orth_norm**2*m%Area(icp,iW)/(Vjouk_orth_norm))
         enddo
      enddo
   end subroutine CirculationFromPolarData



end module FVW_Wings
