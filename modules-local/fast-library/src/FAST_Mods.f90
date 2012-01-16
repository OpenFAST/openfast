!=======================================================================
MODULE ADAMSInput


   ! This MODULE stores FAST-to-ADAMS, ADAMS-specifc input parameters.


USE                             Precision


REAL(ReKi)                   :: BoomRad                                         ! Radius of the tail boom used for tail boom GRAPHICS.
REAL(ReKi)                   :: BPActrDmp                                       ! Blade pitch actuator damping          constant, (N-m/rad/s).
REAL(ReKi)                   :: BPActrSpr                                       ! Blade pitch actuator spring stiffness constant, (N-m/rad).
REAL(ReKi)                   :: CRatioBEA                                       ! The ratio of CMatrix to KMatrix for the blade extensional deflection.
REAL(ReKi)                   :: CRatioBGJ                                       ! The ratio of CMatrix to KMatrix for the blade torsion     deflection.
REAL(ReKi)                   :: CRatioTEA                                       ! The ratio of CMatrix to KMatrix for the tower extensional deflection.
REAL(ReKi)                   :: CRatioTGJ                                       ! The ratio of CMatrix to KMatrix for the tower torsion     deflection.
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for modeling the hydrodynamic loading and
!jmj   mooring system dynamics for floating wind turbines.  Do this by allowing
!jmj   a keyword in place of the integers 0 or 1 in input PtfmLdMod when
!jmj   PtfmModel = 3:
!jmj Also, add an undocumented feature for modeling the hydrodynamic loading on
!jmj   a monopile.  Do this by reading in addition inputs from the platform
!jmj   file if they exist:
REAL(ReKi), PARAMETER        :: FrSrfcSpc =  5.0                                ! Distance between points on the still water level plane along the incident wave propogation heading direction for depicting the free surface where the elevation of the incident waves will be computed used for free surface GRAPHICS. (meters)  !JASON: MAKE THIS AN ACTUAL INPUT TO THE PROGRAM IN ADAMSFile WHEN YOU DOCUMENT THESE ROUTINES!!!!!
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
REAL(ReKi)                   :: GBoxLength                                      ! Length, width, height of the gearbox for gearbox GRAPHICS.
REAL(ReKi)                   :: GenLength                                       ! Length of the generator used for gen. GRAPHICS.
REAL(ReKi)                   :: GenRad                                          ! Radius of the generator used for gen. GRAPHICS.
REAL(ReKi)                   :: HubCylRad                                       ! Radius of hub cylincder used for hub GRAPHICS.
REAL(ReKi)                   :: HSSLength                                       ! Length of high-speed shaft for HSS GRAPHICS.
REAL(ReKi)                   :: HSSRad                                          ! Radius of the high-speed shaft used for HSS GRAPHICS.
REAL(ReKi)                   :: LSSLength                                       ! Length of low-speed shaft for LSS GRAPHICS.
REAL(ReKi)                   :: LSSRad                                          ! Radius of the low-speed shaft used for LSS GRAPHICS.
REAL(ReKi)                   :: NacLength                                       ! Length of nacelle used for the nacelle GRAPHICS.
REAL(ReKi)                   :: NacRadBot                                       ! Bottom radius of nacelle FRUSTUM used for the nacelle GRAPHICS.
REAL(ReKi)                   :: NacRadTop                                       ! Top    radius of nacelle FRUSTUM used for the nacelle GRAPHICS.
REAL(ReKi)                   :: ThkOvrChrd                                      ! Ratio of blade thickness to blade chord used for blade element GRAPHICS.
REAL(ReKi)                   :: TwrBaseRad                                      ! Tower base radius used for linearly tapered tower GRAPHICS.
REAL(ReKi)                   :: TwrTopRad                                       ! Tower top  radius used for linearly tapered tower GRAPHICS.

!jmj Start of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Simplify the SFORCE used to generate free surface GRAPHICS in the
!jmj   FAST-to-ADAMS preprocessor.  Also, eliminate the free surface DOFs
!jmj   during a linearization analysis:
!jmj Also, replace the hard-coded mooring line restoring calculation with a
!jmj   general purpose, quasi-static solution based on the analytical catenary
!jmj   cable equations with seabed interaction:
INTEGER(4)                   :: NFreeSrfc = -1                                  ! Number of points on free surface (not including the zero'th point) where the elevation of the incident waves will be computed (computed every FrSrfcSpc meters along the incident wave propogation heading direction for a length of the rotor diameter).
INTEGER(4)                   :: NLnNodes  = 10                                  ! Number of nodes per line for mooring line GRAPHICS.  !JASON: MAKE THIS AN ACTUAL INPUT TO THE PROGRAM IN ADAMSFile WHEN YOU DOCUMENT THESE ROUTINES!!!!!
!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.
INTEGER(4)                   :: NSides                                          ! The number of sides used in GRAPHICS CYLINDER and FRUSTUM statements.

!bjj chg: LOGICAL(1)                   :: MakeLINacf                                      ! Switch for making an ADAMS/LINEAR control command file.  To prevent an ADAMS/LINEAR control command file to be made, and to not include the RESULTS statement in the ADAMS dataset, set to .FALSE.
LOGICAL                      :: MakeLINacf                                      ! Switch for making an ADAMS/LINEAR control command file.  To prevent an ADAMS/LINEAR control command file to be made, and to not include the RESULTS statement in the ADAMS dataset, set to .FALSE.
!bjj chg: LOGICAL(1)                   :: SaveGrphcs                                      ! Switch to determine whether or note GRAPHICS output is saved in an ADAMS analysis.
LOGICAL                      :: SaveGrphcs                                      ! Switch to determine whether or note GRAPHICS output is saved in an ADAMS analysis.



END MODULE ADAMSInput
!=======================================================================
MODULE AeroElem


   ! This MODULE stores FAST/AeroDyn interface variables.


USE                             Precision
!bjj start of proposed change
USE                             AeroDyn  ! for type;  Precision is also included so the previous line could be removed, too.
!bjj end of proposed change


!bjj rm:REAL(ReKi)                   :: ElAeroLoc (3)                                   ! Vector location of the current aerodynamic element relative to the undeflected tower top.
!bjj start of proposed change
!rmREAL(ReKi)                   :: VES       (3)                                   ! Absolute velocity of blade point S in frame E.
!bjj end of proposed change
!BJJ RM:REAL(ReKi)                   :: ElRad                                           ! Radius of the curent element measured radially from the center of rotor rotation.

!bjj start of proposed change
TYPE(AllAeroMarkers)          :: ADAeroMarkers
!TYPE(Marker)                   :: ADCurrentMarker
TYPE(AeroLoadsOptions)        :: ADIntrfaceOptions
!TYPE(CalcOutput),ALLOCATABLE  :: ADCurrentOutputs(:,:)
TYPE(AllAeroLoads)            :: ADAeroLoads
TYPE(AeroConfig)              :: ADInterfaceComponents                        ! The configuration markers that make up the bodies where aerodynamic calculations will be needed

!REAL(ReKi), ALLOCATABLE, SAVE  :: ElAeroFrc(:,:,:) ! Array of DFN, DFT, PMA aero forces and moment on elements


INTEGER                       :: NumADBldNodes = 0                               ! Number of blade nodes in AeroDyn
!bjj rm:LOGICAL                       :: ADFirstLoop
!bjj rm:REAL(ReKi)                    :: Prev_Aero_t = 0.0   ! Last time aero forces were updated - set to -1 so aero calcs are done at t = 0

!bjj end of proposed change


END MODULE AeroElem
!=======================================================================
MODULE Blades


   ! This MODULE stores input variables for the blades.


USE                             Precision


REAL(ReKi)                   :: AdjBlMs                                         ! Factor to adjust blade mass density.
REAL(ReKi)                   :: AdjEdSt                                         ! Factor to adjust edge stiffness.
REAL(ReKi)                   :: AdjFlSt                                         ! Factor to adjust flap stiffness.
REAL(ReKi), ALLOCATABLE      :: AerCen    (:)                                   ! Aerodynamic center for distributed input data.
REAL(ReKi), ALLOCATABLE      :: AeroCent  (:,:)                                 ! Aerodynamic center for analysis nodes.
REAL(ReKi), ALLOCATABLE      :: AeroTwst  (:)                                   ! Aerodynamic twist of the blade at the analysis nodes.
REAL(ReKi), ALLOCATABLE      :: Alpha     (:)                                   ! Blade coupling coefficient between flap and twist for a given input station.
REAL(ReKi), ALLOCATABLE      :: AxRedBld  (:,:,:,:)                             ! The axial-reduction terms of the blade shape function.
REAL(ReKi), ALLOCATABLE      :: BAlpha    (:,:)                                 ! Interpolated blade coupling coefficient between flap and twist.
REAL(ReKi), ALLOCATABLE      :: BldEDamp  (:,:)                                 ! Blade edgewise damping coefficients.
REAL(ReKi)                   :: BldEdDmp  (1)                                   ! Blade structural damping ratios in edgewise direction.
REAL(ReKi), ALLOCATABLE      :: BldFDamp  (:,:)                                 ! Blade flapwise damping coefficients.
REAL(ReKi)                   :: BldFlDmp  (2)                                   ! Blade structural damping ratios in flapwise direction.
REAL(ReKi)                   :: BldFlexL                                        ! Flexible blade length.
REAL(ReKi), ALLOCATABLE      :: BlFract   (:)                                   ! Blade fractional radius for distributed input data.
REAL(ReKi), ALLOCATABLE      :: BMassDen  (:)                                   ! Blade mass density for distributed input data.
REAL(ReKi), ALLOCATABLE      :: CAeroTwst (:)                                   ! Cosine of the aerodynamic twist of the blade at the analysis nodes.
REAL(ReKi), ALLOCATABLE      :: CBE       (:,:,:)                               ! Generalized edgewise damping of the blades.
REAL(ReKi), ALLOCATABLE      :: CBF       (:,:,:)                               ! Generalized flapwise damping of the blades.
REAL(ReKi), ALLOCATABLE      :: cgOffBEdg (:,:)                                 ! Interpolated blade edge (along local aerodynamic yb-axis) mass cg offset.
REAL(ReKi), ALLOCATABLE      :: cgOffBFlp (:,:)                                 ! Interpolated blade flap (along local aerodynamic xb-axis) mass cg offset.
REAL(ReKi), ALLOCATABLE      :: Chord     (:)                                   ! Chord of the blade at the analysis nodes.
REAL(ReKi), ALLOCATABLE      :: CThetaS   (:,:)                                 ! COS( ThetaS )
REAL(ReKi), ALLOCATABLE      :: DRNodes   (:)                                   ! Length of variable-spaced blade elements.
REAL(ReKi), ALLOCATABLE      :: EAOffBEdg (:,:)                                 ! Interpolated blade edge (along local aerodynamic yb-axis) elastic axis offset.
REAL(ReKi), ALLOCATABLE      :: EAOffBFlp (:,:)                                 ! Interpolated blade flap (along local aerodynamic xb-axis) elastic axis offset.
REAL(ReKi), ALLOCATABLE      :: EAStff    (:)                                   ! Blade extensional stiffness for a given input station.
REAL(ReKi), ALLOCATABLE      :: EdgcgOf   (:)                                   ! Blade edge (along local aerodynamic yb-axis) mass cg offset for a given input station.
REAL(ReKi), ALLOCATABLE      :: EdgEAOf   (:)                                   ! Blade edge (along local aerodynamic yb-axis) elastic axis offset for a given input station.
REAL(ReKi), ALLOCATABLE      :: EdgIner   (:)                                   ! Blade edge (about local structural xb-axis) mass inertia per unit length for a given input station.
REAL(ReKi), ALLOCATABLE      :: EdgStff   (:)                                   ! Blade edge stiffness for distributed input data.
REAL(ReKi), ALLOCATABLE      :: FlpcgOf   (:)                                   ! Blade flap (along local aerodynamic xb-axis) mass cg offset for a given input station.
REAL(ReKi), ALLOCATABLE      :: FlpEAOf   (:)                                   ! Blade flap (along local aerodynamic xb-axis) elastic axis offset for a given input station.
REAL(ReKi), ALLOCATABLE      :: FlpIner   (:)                                   ! Blade flap (about local structural yb-axis) mass inertia per unit length for a given input station.
REAL(ReKi), ALLOCATABLE      :: FlpStff   (:)                                   ! Blade flap stiffness for distributed input data.
REAL(ReKi), ALLOCATABLE      :: FStTunr   (:,:)                                 ! Blade flapwise modal stiffness tuners (stored for all blades).
REAL(ReKi)                   :: FlStTunr  (2)                                   ! Blade flapwise modal stiffness tuners (input).
REAL(ReKi), ALLOCATABLE      :: GJStff    (:)                                   ! Blade torsional stiffness for a given input station.
REAL(ReKi), ALLOCATABLE      :: InerBEdg  (:,:)                                 ! Interpolated blade edge (about local structural xb-axis) mass inertia per unit length.
REAL(ReKi), ALLOCATABLE      :: InerBFlp  (:,:)                                 ! Interpolated blade flap (about local structural yb-axis) mass inertia per unit length.
REAL(ReKi), ALLOCATABLE      :: KBE       (:,:,:)                               ! Generalized edgewise stiffness of the blades.
REAL(ReKi), ALLOCATABLE      :: KBF       (:,:,:)                               ! Generalized flapwise stiffness of the blades.
REAL(ReKi), ALLOCATABLE      :: MassB     (:,:)                                 ! Interpolated lineal blade mass density.
REAL(ReKi), ALLOCATABLE      :: PrecrvRef (:)                                   ! Offset for defining the reference axis from the pitch axis for precurved blades at a given input station.
REAL(ReKi), ALLOCATABLE      :: PreswpRef (:)                                   ! Offset for defining the reference axis from the pitch axis for preswept  blades at a given input station.
REAL(ReKi), ALLOCATABLE      :: RefAxisxb (:,:)                                 ! Interpolated Offset for defining the reference axis from the pitch axis for precurved blades at a given input station (along xb-axis).
REAL(ReKi), ALLOCATABLE      :: RefAxisyb (:,:)                                 ! Interpolated Offset for defining the reference axis from the pitch axis for preswept  blades at a given input station (along yb-axis).
REAL(ReKi), ALLOCATABLE      :: RNodes    (:)                                   ! Radius to analysis nodes relative to hub ( 0 < RNodes(:) < BldFlexL )
REAL(ReKi), ALLOCATABLE      :: RNodesNorm(:)                                   ! Normalized radius to analysis nodes relative to hub ( 0 < RNodesNorm(:) < 1 )
REAL(ReKi), ALLOCATABLE      :: rSAerCenn1(:,:)                                 ! Distance from point S on a blade to the aerodynamic center in the n1 direction (m).
REAL(ReKi), ALLOCATABLE      :: rSAerCenn2(:,:)                                 ! Distance from point S on a blade to the aerodynamic center in the n2 direction (m).
REAL(ReKi), ALLOCATABLE      :: SAeroTwst (:)                                   ! Sine of the aerodynamic twist of the blade at the analysis nodes.
REAL(ReKi), ALLOCATABLE      :: StiffBE   (:,:)                                 ! Interpolated edgewise blade stiffness.
REAL(ReKi), ALLOCATABLE      :: StiffBEA  (:,:)                                 ! Interpolated blade extensional stiffness.
REAL(ReKi), ALLOCATABLE      :: StiffBF   (:,:)                                 ! Interpolated flapwise blade stiffness.
REAL(ReKi), ALLOCATABLE      :: StiffBGJ  (:,:)                                 ! Interpolated blade torsional stiffness.
REAL(ReKi), ALLOCATABLE      :: SThetaS   (:,:)                                 ! SIN( ThetaS )
REAL(ReKi), ALLOCATABLE      :: StrcTwst  (:)                                   ! Structural twist for distributed input data.
REAL(ReKi), ALLOCATABLE      :: ThetaS    (:,:)                                 ! Structural twist for analysis nodes.
REAL(ReKi), ALLOCATABLE      :: TwistedSF (:,:,:,:,:)                           ! Interpolated lineal blade mass density.

INTEGER(4)                   :: BldNodes                                        ! Number of blade nodes used in the analysis.
INTEGER(4)                   :: NBlInpSt                                        ! Number of blade input stations.
INTEGER(4)                   :: TipNode                                         ! Index of the additional node located at the blade tip = BldNodes + 1


END MODULE Blades
!=======================================================================
MODULE CoordSys


   ! This MODULE stores coordinate sytems used internally by FAST.  The 3
   !   components of each vector correspond to the z1, z2, and z3 components
   !   of the individual vectors.
   ! NOTE: the orientations of most of these coordinate systems will change
   !   every time step.


USE                             Precision


REAL(ReKi)                   :: a1       (3)                                    ! Vector / direction a1 (=  xt from the IEC coord. system).
REAL(ReKi)                   :: a2       (3)                                    ! Vector / direction a2 (=  zt from the IEC coord. system).
REAL(ReKi)                   :: a3       (3)                                    ! Vector / direction a3 (= -yt from the IEC coord. system).
REAL(ReKi)                   :: b1       (3)                                    ! Vector / direction b1 (=  xp from the IEC coord. system).
REAL(ReKi)                   :: b2       (3)                                    ! Vector / direction b2 (=  zp from the IEC coord. system).
REAL(ReKi)                   :: b3       (3)                                    ! Vector / direction b3 (= -yp from the IEC coord. system).
REAL(ReKi)                   :: c1       (3)                                    ! Vector / direction c1 (=  xs from the IEC coord. system).
REAL(ReKi)                   :: c2       (3)                                    ! Vector / direction c2 (=  zs from the IEC coord. system).
REAL(ReKi)                   :: c3       (3)                                    ! Vector / direction c3 (= -ys from the IEC coord. system).
REAL(ReKi)                   :: d1       (3)                                    ! Vector / direction d1 (=  xn from the IEC coord. system).
REAL(ReKi)                   :: d2       (3)                                    ! Vector / direction d2 (=  zn from the IEC coord. system).
REAL(ReKi)                   :: d3       (3)                                    ! Vector / direction d3 (= -yn from the IEC coord. system).
REAL(ReKi)                   :: e1       (3)                                    ! Vector / direction e1 (=  xa from the IEC coord. system).
REAL(ReKi)                   :: e2       (3)                                    ! Vector / direction e2 (=  ya from the IEC coord. system).
REAL(ReKi)                   :: e3       (3)                                    ! Vector / direction e3 (=  za from the IEC coord. system).
REAL(ReKi)                   :: f1       (3)                                    ! Vector / direction f1.
REAL(ReKi)                   :: f2       (3)                                    ! Vector / direction f2.
REAL(ReKi)                   :: f3       (3)                                    ! Vector / direction f3.
REAL(ReKi)                   :: g1       (3)                                    ! Vector / direction g1 (=  xh from the IEC coord. system).
REAL(ReKi)                   :: g2       (3)                                    ! Vector / direction g2 (=  yh from the IEC coord. system).
REAL(ReKi)                   :: g3       (3)                                    ! Vector / direction g3 (=  zh from the IEC coord. system).
REAL(ReKi), ALLOCATABLE      :: i1       (:,:)                                  ! i1(K,:) = vector / direction i1 for blade K (=  xcK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE      :: i2       (:,:)                                  ! i2(K,:) = vector / direction i2 for blade K (=  ycK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE      :: i3       (:,:)                                  ! i3(K,:) = vector / direction i3 for blade K (=  zcK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE      :: j1       (:,:)                                  ! j1(K,:) = vector / direction j1 for blade K (=  xbK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE      :: j2       (:,:)                                  ! j2(K,:) = vector / direction j2 for blade K (=  ybK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE      :: j3       (:,:)                                  ! j3(K,:) = vector / direction j3 for blade K (=  zbK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE      :: m1       (:,:,:)                                ! m1(K,J,:) = vector / direction m1 for node J of blade K (used to calc. and return aerodynamic loads from AeroDyn).
REAL(ReKi), ALLOCATABLE      :: m2       (:,:,:)                                ! m2(K,J,:) = vector / direction m2 for node J of blade K (used to calc. and return aerodynamic loads from AeroDyn).
REAL(ReKi), ALLOCATABLE      :: m3       (:,:,:)                                ! m3(K,J,:) = vector / direction m3 for node J of blade K (used to calc. and return aerodynamic loads from AeroDyn).
REAL(ReKi), ALLOCATABLE      :: n1       (:,:,:)                                ! n1(K,J,:) = vector / direction n1 for node J of blade K (= LxbK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE      :: n2       (:,:,:)                                ! n2(K,J,:) = vector / direction n2 for node J of blade K (= LybK from the IEC coord. system).
REAL(ReKi), ALLOCATABLE      :: n3       (:,:,:)                                ! n3(K,J,:) = vector / direction n3 for node J of blade K (= LzbK from the IEC coord. system).
REAL(ReKi)                   :: p1       (3)                                    ! Vector / direction p1 (used to calc. and return tail aerodynamic loads from AeroDyn).
REAL(ReKi)                   :: p2       (3)                                    ! Vector / direction p2 (used to calc. and return tail aerodynamic loads from AeroDyn).
REAL(ReKi)                   :: p3       (3)                                    ! Vector / direction p3 (used to calc. and return tail aerodynamic loads from AeroDyn).
REAL(ReKi)                   :: rf1      (3)                                    ! Vector / direction rf1 (rotor-furl coordinate system = d1 when rotor-furl angle = 0).
REAL(ReKi)                   :: rf2      (3)                                    ! Vector / direction rf2 (rotor-furl coordinate system = d2 when rotor-furl angle = 0).
REAL(ReKi)                   :: rf3      (3)                                    ! Vector / direction rf3 (rotor-furl coordinate system = d3 when rotor-furl angle = 0).
REAL(ReKi)                   :: rfa      (3)                                    ! Vector / direction of the rotor-furl axis.
REAL(ReKi), ALLOCATABLE      :: t1       (:,:)                                  ! Vector / direction t1 for tower node J (=  Lxt from the IEC coord. system).
REAL(ReKi), ALLOCATABLE      :: t2       (:,:)                                  ! Vector / direction t2 for tower node J (=  Lzt from the IEC coord. system).
REAL(ReKi), ALLOCATABLE      :: t3       (:,:)                                  ! Vector / direction t3 for tower node J (= -Lyt from the IEC coord. system).
REAL(ReKi), ALLOCATABLE      :: te1      (:,:,:)                                ! te1(K,J,:) = vector / direction te1 for node J of blade K (used to calc. noise).
REAL(ReKi), ALLOCATABLE      :: te2      (:,:,:)                                ! te2(K,J,:) = vector / direction te2 for node J of blade K (used to calc. noise).
REAL(ReKi), ALLOCATABLE      :: te3      (:,:,:)                                ! te3(K,J,:) = vector / direction te3 for node J of blade K (used to calc. noise).
REAL(ReKi)                   :: tf1      (3)                                    ! Vector / direction tf1 (tail-furl coordinate system = d1 when rotor-furl angle = 0).
REAL(ReKi)                   :: tf2      (3)                                    ! Vector / direction tf2 (tail-furl coordinate system = d2 when rotor-furl angle = 0).
REAL(ReKi)                   :: tf3      (3)                                    ! Vector / direction tf3 (tail-furl coordinate system = d3 when rotor-furl angle = 0).
REAL(ReKi)                   :: tfa      (3)                                    ! Vector / direction of the tail-furl axis.
REAL(ReKi)                   :: z1       (3)                                    ! Vector / direction z1 (=  xi from the IEC coord. system).
REAL(ReKi)                   :: z2       (3)                                    ! Vector / direction z2 (=  zi from the IEC coord. system).
REAL(ReKi)                   :: z3       (3)                                    ! Vector / direction z3 (= -yi from the IEC coord. system).


END MODULE CoordSys
!=======================================================================
MODULE Constants


  ! This MODULE stores various constants.


USE                             Precision


!rmREAL(ReKi), PARAMETER        :: D2R      =  0.017453293                         ! Factor to convert degrees to radians.
REAL(ReKi), PARAMETER        :: Inv2Pi   =  0.15915494                          ! 0.5/Pi.
!rmREAL(ReKi), PARAMETER        :: Pi       =  3.1415927                           ! Ratio of a circle's circumference to its diameter.
!rmREAL(ReKi), PARAMETER        :: PiOvr2   =  1.5707963                           ! Pi/2.
!rmREAL(ReKi), PARAMETER        :: R2D      = 57.295780                            ! Factor to convert radians to degrees.
!rmREAL(ReKi), PARAMETER        :: RPM2RPS  =  0.10471976                          ! Factor to convert revolutions per minute to radians per second.
!rmREAL(ReKi), PARAMETER        :: RPS2RPM  =  9.5492966                           ! Factor to convert radians per second to revolutions per minute.
!rmREAL(ReKi), PARAMETER        :: TwoPi    =  6.2831853                           ! 2*Pi.
REAL(ReKi)                   :: TwoPiNB                                         ! 2*Pi/NumBl.  This constant is calculated in fast_io.f90/Inputs()
!rm
!rm
END MODULE Constants
!=======================================================================
MODULE DOFs


   ! This MODULE stores variables related to degrees of freedom.


USE                             Precision


REAL(ReKi), ALLOCATABLE      :: Q        (:,:)                                  ! Displacement matrix.
REAL(ReKi), ALLOCATABLE      :: QD       (:,:)                                  ! Velocity matrix.
REAL(ReKi), ALLOCATABLE      :: QD2      (:,:)                                  ! Acceleration matrix.

INTEGER(4), ALLOCATABLE      :: Diag     (:)                                    ! Array containing the indices of SrtPS() associated with each enabled DOF; that is, SrtPS(Diag(I)) = I.
INTEGER(4), ALLOCATABLE      :: DOF_BE   (:,:)                                  ! DOF indices for blade edge.
INTEGER(4), ALLOCATABLE      :: DOF_BF   (:,:)                                  ! DOF indices for blade flap.
INTEGER(4), PARAMETER        :: DOF_DrTr = 14                                   ! DOF index for drivetrain rotational-flexibility.
INTEGER(4), PARAMETER        :: DOF_GeAz = 13                                   ! DOF index for the generator azimuth.
INTEGER(4), PARAMETER        :: DOF_Hv   =  3                                   ! DOF index for platform heave.
INTEGER(4), PARAMETER        :: DOF_P    =  5                                   ! DOF index for platform pitch.
INTEGER(4), PARAMETER        :: DOF_R    =  4                                   ! DOF index for platform roll.
INTEGER(4), PARAMETER        :: DOF_RFrl = 12                                   ! DOF index for rotor-furl.
INTEGER(4), PARAMETER        :: DOF_Sg   =  1                                   ! DOF index for platform surge.
INTEGER(4), PARAMETER        :: DOF_Sw   =  2                                   ! DOF index for platform sway.
INTEGER(4), PARAMETER        :: DOF_Teet = 22                                   ! DOF index for rotor-teeter.
INTEGER(4), PARAMETER        :: DOF_TFA1 =  7                                   ! DOF index for 1st tower fore-aft mode.
INTEGER(4), PARAMETER        :: DOF_TFA2 =  9                                   ! DOF index for 2nd tower fore-aft mode.
INTEGER(4), PARAMETER        :: DOF_TFrl = 15                                   ! DOF index for tail-furl.
INTEGER(4), PARAMETER        :: DOF_TSS1 =  8                                   ! DOF index for 1st tower side-to-side mode.
INTEGER(4), PARAMETER        :: DOF_TSS2 = 10                                   ! DOF index for 2nd tower side-to-side mode.
INTEGER(4), PARAMETER        :: DOF_Y    =  6                                   ! DOF index for platform yaw.
INTEGER(4), PARAMETER        :: DOF_Yaw  = 11                                   ! DOF index for nacelle-yaw.
INTEGER(4), ALLOCATABLE      :: IC       (:)                                    ! Array which stores pointers to predictor-corrector results.
INTEGER(4)                   :: NActvDOF                                        ! The number of active (enabled) DOFs in the model.
INTEGER(4)                   :: NAug                                            ! Dimension of augmented solution matrix.
INTEGER(4)                   :: NDOF                                            ! Number of total DOFs.
INTEGER(4), PARAMETER        :: NMX      =  9                                   ! Used in updating predictor-corrector values.
INTEGER(4)                   :: NPA                                             ! Number of DOFs                  that contribute to the angular velocity of the tail                                                      (body A) in the inertia frame.
INTEGER(4)                   :: NPB                                             ! Number of DOFs                  that contribute to the angular velocity of the tower top / baseplate                                     (body B) in the inertia frame.
INTEGER(4)                   :: NPCE                                            ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the hub center of mass                                                              (point C) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4)                   :: NPDE                                            ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the center of mass of the structure that furls with the rotor (not including rotor) (point D) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4)                   :: NPF                                             ! Number of DOFs                  that contribute to the angular velocity of the tower elements                                            (body F) in the inertia frame.
INTEGER(4)                   :: NPG                                             ! Number of DOFs                  that contribute to the angular velocity of the generator                                                 (body G) in the inertia frame.
INTEGER(4)                   :: NPH                                             ! Number of DOFs                  that contribute to the angular velocity of the hub                                                       (body H) in the inertia frame.
INTEGER(4)                   :: NPIE                                            ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the tail boom center of mass                                                        (point I) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4)                   :: NPL                                             ! Number of DOFs                  that contribute to the angular velocity of the low-speed shaft                                           (body L) in the inertia frame.
INTEGER(4)                   :: NPM                                             ! Number of DOFs                  that contribute to the angular velocity of the blade elements                                            (body M) in the inertia frame.
INTEGER(4)                   :: NPN                                             ! Number of DOFs                  that contribute to the angular velocity of the nacelle                                                   (body N) in the inertia frame.
INTEGER(4)                   :: NPTE                                            ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the tower nodes                                                                     (point T) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4)                   :: NPTTE                                           ! Number of tower DOFs            that contribute to the QD2T-related linear accelerations of the tower nodes                                                                     (point T) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4)                   :: NPR                                             ! Number of DOFs                  that contribute to the angular velocity of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame.
INTEGER(4), ALLOCATABLE      :: NPSBE    (:)                                    ! Number of blade DOFs            that contribute to the QD2T-related linear accelerations of the blade nodes                                                                     (point S) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE      :: NPSE     (:)                                    ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the blade nodes                                                                     (point S) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4)                   :: NPUE                                            ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the nacelle center of mass                                                          (point U) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4)                   :: NPX                                             ! Number of DOFs                  that contribute to the angular velocity of the platform                                                  (body X) in the inertia frame.
INTEGER(4)                   :: NPYE                                            ! Number of DOFs                  that contribute to the QD2T-related linear accelerations of the platform center of mass                                                         (point Y) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE      :: PA       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the tail                                                      (body A) in the inertia frame.
INTEGER(4), ALLOCATABLE      :: PB       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the tower top / baseplate                                     (body B) in the inertia frame.
INTEGER(4), ALLOCATABLE      :: PCE      (:)                                    ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the hub center of mass                                                              (point C) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE      :: PDE      (:)                                    ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the center of mass of the structure that furls with the rotor (not including rotor) (point D) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE      :: PF       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the tower elements                                            (body F) in the inertia frame.
INTEGER(4), ALLOCATABLE      :: PG       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the generator                                                 (body G) in the inertia frame.
INTEGER(4), ALLOCATABLE      :: PH       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the hub                                                       (body H) in the inertia frame.
INTEGER(4), ALLOCATABLE      :: PIE      (:)                                    ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the tail boom center of mass                                                        (point I) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE      :: PL       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the low-speed shaft                                           (body L) in the inertia frame.
INTEGER(4), ALLOCATABLE      :: PM       (:,:)                                  ! Array of DOF indices (pointers) that contribute to the angular velocity of the blade elements                                            (body M) in the inertia frame.
INTEGER(4), ALLOCATABLE      :: PN       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the nacelle                                                   (body N) in the inertia frame.
INTEGER(4), ALLOCATABLE      :: PTE      (:)                                    ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the tower nodes                                                                     (point T) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE      :: PTTE     (:)                                    ! Array of tower DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the tower nodes                                                               (point T) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE      :: PR       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame.
INTEGER(4), ALLOCATABLE      :: PS       (:)                                    ! Array of DOF indices (pointers) to the active (enabled) DOFs/states.
INTEGER(4), ALLOCATABLE      :: PSBE     (:,:)                                  ! Array of blade DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the blade nodes                                                               (point S) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE      :: PSE      (:,:)                                  ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the blade nodes                                                                     (point S) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE      :: PUE      (:)                                    ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the nacelle center of mass                                                          (point U) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE      :: PX       (:)                                    ! Array of DOF indices (pointers) that contribute to the angular velocity of the platform                                                  (body X) in the inertia frame.
INTEGER(4), ALLOCATABLE      :: PYE      (:)                                    ! Array of DOF indices (pointers) that contribute to the QD2T-related linear accelerations of the platform center of mass                                                         (point Y) in the inertia frame, based on which DOFs are presently enabled.
INTEGER(4), ALLOCATABLE      :: SrtPS    (:)                                    ! Sorted (from smallest to largest DOF index) version of PS().
INTEGER(4), ALLOCATABLE      :: SrtPSNAUG(:)                                    ! SrtPS() with the additional value of NAUG.

!bjj chg: LOGICAL(1), ALLOCATABLE      :: DOF_Flag (:)                                    ! Array which stores values of the feature flags for each DOF.
LOGICAL,    ALLOCATABLE      :: DOF_Flag (:)                                    ! Array which stores values of the feature flags for each DOF.

CHARACTER(99), ALLOCATABLE   :: DOF_Desc (:)                                    ! Array which stores descriptions of each DOF.


END MODULE DOFs
!=======================================================================
MODULE DriveTrain


   ! This MODULE stores variables for the drivetrain.


USE                             Precision


REAL(ReKi)                   :: DTTorDmp                                        ! Drivetrain torsional damper
REAL(ReKi)                   :: DTTorSpr                                        ! Drivetrain torsional spring
REAL(ReKi)                   :: ElecPwr                                         ! Electrical power, W.
REAL(ReKi)                   :: GBRatio                                         ! Gearbox ratio
REAL(ReKi)                   :: GBoxEff                                         ! Gearbox efficiency.
REAL(ReKi)                   :: GenCTrq                                         ! Constant generator torque.
REAL(ReKi)                   :: GenEff                                          ! Generator efficiency
REAL(ReKi)                   :: GenSpRZT                                        ! Difference between rated and zero-torque generator speeds for SIG.
REAL(ReKi)                   :: GenSpRat                                        ! Rated generator speed.
REAL(ReKi)                   :: GenSpZT                                         ! Zero-torque generator speed.
REAL(ReKi)                   :: GenTrq                                          ! Electrical generator torque.
REAL(ReKi)                   :: HSSBrDT                                         ! Time it takes for HSS brake to reach full deployment once deployed.
REAL(ReKi)                   :: HSSBrTqF                                        ! Fully deployed HSS brake torque
REAL(ReKi)                   :: HSSBrTrq                                        ! Instantaneous HSS brake torque
REAL(ReKi)                   :: HSSBrTrqC                                       ! A copy of the value of HSSBrTrq calculated in SUBROUTINE DrvTrTrq().
REAL(ReKi)                   :: SIG_PORt                                        ! Pull-out ratio (Tpullout/Trated).
REAL(ReKi)                   :: SIG_POSl                                        ! Pullout slip.
REAL(ReKi)                   :: SIG_POTq                                        ! Pullout torque.
REAL(ReKi)                   :: SIG_RtSp                                        ! Rated speed.
REAL(ReKi)                   :: SIG_RtTq                                        ! Rated torque.
REAL(ReKi)                   :: SIG_SlPc                                        ! Rated generator slip percentage.
REAL(ReKi)                   :: SIG_Slop                                        ! Torque/Speed slope for simple induction generator.
REAL(ReKi)                   :: SIG_SySp                                        ! Synchronous (zero-torque) generator speed.
REAL(ReKi)                   :: TEC_A0                                          ! A0 term for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_C0                                          ! C0 term for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_C1                                          ! C1 term for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_C2                                          ! C2 term for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_Freq                                        ! Line frequency for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_K1                                          ! K1 term for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_K2                                          ! K2 term for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_MR                                          ! Magnetizing reactance for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_Re1                                         ! Thevenin's equivalent stator resistance (ohms)
REAL(ReKi)                   :: TEC_RLR                                         ! Rotor leakage reactance for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_RRes                                        ! Rotor resistance for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_SLR                                         ! Stator leakage reactance for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_SRes                                        ! Stator resistance for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_SySp                                        ! Synchronous speed for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_V1a                                         ! Source voltage for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_VLL                                         ! Line-to-line RMS voltage for Thevenin-equivalent circuit.
REAL(ReKi)                   :: TEC_Xe1                                         ! Thevenin's equivalent stator leakage reactance (ohms)

INTEGER(4)                   :: GenDir                                          ! Direction of the generator = +/- 1 (+ 1 = same direction as LSS; -1 = opposite direction of LSS).
INTEGER(4)                   :: TEC_NPol                                        ! Number of poles for Thevenin-equivalent circuit.

!bjj chg: LOGICAL(1)                   :: GBRevers                                        ! Gearbox reversal flag.
LOGICAL                      :: GBRevers                                        ! Gearbox reversal flag.


END MODULE DriveTrain
!=======================================================================
MODULE EnvCond


   ! This MODULE stores input variables for environmental conditions.


USE                             Precision


REAL(ReKi)                   :: AirDens                                         ! Air density = RHO.
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for modeling the hydrodynamic loading and
!jmj   mooring system dynamics for floating wind turbines.  Do this by allowing
!jmj   a keyword in place of the integers 0 or 1 in input PtfmLdMod when
!jmj   PtfmModel = 3:
!jmj Also, add an undocumented feature for modeling the hydrodynamic loading on
!jmj   a monopile.  Do this by reading in addition inputs from the platform
!jmj   file if they exist:
REAL(ReKi)                   :: CurrDIDir = 0.0                                 ! Depth-independent current heading direction.
REAL(ReKi)                   :: CurrDIV   = 0.0                                 ! Depth-independent current velocity.
REAL(ReKi)                   :: CurrNSDir = 0.0                                 ! Near-surface current heading direction.
REAL(ReKi)                   :: CurrNSRef = 0.0                                 ! Near-surface current reference depth.
REAL(ReKi)                   :: CurrNSV0  = 0.0                                 ! Near-surface current velocity at still water level.
REAL(ReKi)                   :: CurrSSDir = 0.0                                 ! Sub-surface current heading direction.
REAL(ReKi)                   :: CurrSSV0  = 0.0                                 ! Sub-surface current velocity at still water level.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
REAL(ReKi)                   :: Gravity                                         ! Gravitational acceleration.

!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for modeling the hydrodynamic loading and
!jmj   mooring system dynamics for floating wind turbines.  Do this by allowing
!jmj   a keyword in place of the integers 0 or 1 in input PtfmLdMod when
!jmj   PtfmModel = 3:
!jmj Also, add an undocumented feature for modeling the hydrodynamic loading on
!jmj   a monopile.  Do this by reading in addition inputs from the platform
!jmj   file if they exist:
REAL(ReKi)                   :: WaveDir   = 0.0                                 ! Wave heading direction.
REAL(ReKi)                   :: WaveDT    = 0.0                                 ! Time step for incident wave calculations.
REAL(ReKi)                   :: WaveHs    = 0.0                                 ! Significant wave height.
REAL(ReKi)                   :: WavePkShp = 1.0                                 ! Peak shape parameter of incident wave spectrum.
REAL(ReKi)                   :: WaveTMax  = 0.0                                 ! Analysis time for incident wave calculations.
REAL(ReKi)                   :: WaveTp    = 0.0                                 ! Peak spectral period.
REAL(ReKi)                   :: WtrDens                                         ! Water density.
REAL(ReKi)                   :: WtrDpth                                         ! Water depth.

INTEGER(4)                   :: CurrMod                                         ! Current profile model switch.
INTEGER(4)                   :: WaveStMod = 0                                   ! Model switch for stretching incident wave kinematics to instantaneous free surface.
INTEGER(4)                   :: WaveMod   = 0                                   ! Incident wave kinematics model switch.
INTEGER(4)                   :: WaveSeed (2) = 0                                ! Random seeds of incident waves.

!bjj start of proposed change v6.02d-bjj
!rmCHARACTER(99)                :: GHWvFile  = ''                                  ! The root name of GH Bladed files containing wave data.
CHARACTER(1024)              :: GHWvFile  = ''                                  ! The root name of GH Bladed files containing wave data.
!bjj end of proposed change v6.02d-bjj
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.


END MODULE EnvCond
!=======================================================================
MODULE Features


   ! This MODULE stores input variables for feature switches.

!bjj chg: changed all of these from LOGICAL(1) to LOGICAL
LOGICAL                   :: CompAero                                        ! Compute aerodynamic forces switch.
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for modeling the hydrodynamic loading and
!jmj   mooring system dynamics for floating wind turbines.  Do this by allowing
!jmj   a keyword in place of the integers 0 or 1 in input PtfmLdMod when
!jmj   PtfmModel = 3:
!jmj Also, add an undocumented feature for modeling the hydrodynamic loading on
!jmj   a monopile.  Do this by reading in addition inputs from the platform
!jmj   file if they exist:
LOGICAL                   :: CompHydro = .FALSE.                             ! Compute hydrodynamic forces switch.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.

LOGICAL                   :: CompNoise                                       ! Compute aerodynamic noise  switch.
LOGICAL                   :: DrTrDOF                                         ! Drivetrain rotational-flexibility DOF.
LOGICAL                   :: EdgeDOF                                         ! Edgewise blade mode DOF.
LOGICAL                   :: FlapDOF1                                        ! First flapwise blade mode DOF.
LOGICAL                   :: FlapDOF2                                        ! Second flapwise blade mode DOF.
LOGICAL                   :: GenDOF                                          ! Generator DOF.
LOGICAL                   :: PtfmHvDOF = .FALSE.                             ! Platform vertical heave translation DOF. (Initialized to .FALSE. b/c not all models will read in PtfmFile)
LOGICAL                   :: PtfmPDOF  = .FALSE.                             ! Platform pitch tilt rotation DOF. (Initialized to .FALSE. b/c not all models will read in PtfmFile)
LOGICAL                   :: PtfmRDOF  = .FALSE.                             ! Platform roll tilt rotation DOF. (Initialized to .FALSE. b/c not all models will read in PtfmFile)
LOGICAL                   :: PtfmSgDOF = .FALSE.                             ! Platform horizontal surge translation DOF. (Initialized to .FALSE. b/c not all models will read in PtfmFile)
LOGICAL                   :: PtfmSwDOF = .FALSE.                             ! Platform horizontal sway translation DOF. (Initialized to .FALSE. b/c not all models will read in PtfmFile)
LOGICAL                   :: PtfmYDOF  = .FALSE.                             ! Platform yaw rotation DOF. (Initialized to .FALSE. b/c not all models will read in PtfmFile)
LOGICAL                   :: RFrlDOF   = .FALSE.                             ! Rotor-furl DOF. (Initialized to .FALSE. b/c not all models read in FurlFile)
LOGICAL                   :: TeetDOF   = .FALSE.                             ! Rotor-teeter DOF. (Initialized to .FALSE. b/c the 3-blader requires it to be .FALSE.)
LOGICAL                   :: TFrlDOF   = .FALSE.                             ! Tail-furl DOF. (Initialized to .FALSE. b/c not all models read in FurlFile)
LOGICAL                   :: TwFADOF1                                        ! First tower fore-aft bending-mode DOF.
LOGICAL                   :: TwFADOF2                                        ! Second tower fore-aft bending-mode DOF.
LOGICAL                   :: TwSSDOF1                                        ! First tower side-to-side bending-mode DOF.
LOGICAL                   :: TwSSDOF2                                        ! Second tower side-to-side bending-mode DOF.
LOGICAL                   :: YawDOF                                          ! Nacelle-yaw DOF.


END MODULE Features
!=======================================================================
MODULE General


   ! This MODULE stores input variables for general program control.


INTEGER(4)                   :: ADAMSPrep                                       ! ADAMS preprocessor mode {1: Run FAST, 2: use FAST as a preprocessor to create equivalent ADAMS model, 3: do both} (switch).
INTEGER(4)                   :: AnalMode                                        ! FAST analysis mode {1: Run a time-marching simulation, 2: create a periodic linearized model} (switch).
!bjj rmINTEGER(4), PARAMETER        :: FlagType  = 1                                   ! Switch for telling if a variable is a flag.
!bjj rmINTEGER(4), PARAMETER        :: Numeric   = 2                                   ! Switch for telling if a variable is a number.
INTEGER(4)                   :: PtfmModel                                       ! Platform model {0: none, 1: onshore, 2: fixed bottom offshore, 3: floating offshore} (switch).
!bjj rmINTEGER(4), PARAMETER        :: String    = 3                                   ! Switch for telling if a variable is a string.
INTEGER(4)                   :: StrtTime (8)                                    ! Start time of simulation.
INTEGER(4)                   :: UnAC      = 24                                  ! I/O unit number for the ADAMS control output file (.acf) useful for an ADAMS SIMULATE analysis.
INTEGER(4)                   :: UnAD      = 23                                  ! I/O unit number for the ADAMS dataset output file (.adm).
INTEGER(4)                   :: UnAL      = 25                                  ! I/O unit number for the ADAMS control output file (.acf) useful for an ADAMS LINEAR analysis.
!bjj rm: INTEGER(4)                   :: UnEc      = 19                                  ! I/O unit number for the echo file. -- bjj: in NWTC_Library now
INTEGER(4)                   :: UnIn      = 20                                  ! I/O unit number for the input files.
INTEGER(4)                   :: UnLn      = 26                                  ! I/O unit number for the FAST linear output file (.lin).
INTEGER(4)                   :: UnNoSpec  = 27                                  ! I/O unit number for the noise spectr output file.
INTEGER(4)                   :: UnNoSPL   = 28                                  ! I/O unit number for the noise SPL output file.
INTEGER(4)                   :: UnOu      = 21                                  ! I/O unit number for the tabular output file.
INTEGER(4)                   :: UnSu      = 22                                  ! I/O unit number for the summary output file.

!bjj chg: changed all 4 from LOGICAL(1) to LOGICAL
LOGICAL                      :: Cmpl4SFun = .FALSE.                             ! Is FAST being compiled as an S-Function for Simulink?
LOGICAL                      :: Furling                                         ! Read in additional model properties for furling turbine?
LOGICAL                      :: SumDisp                                         ! Display summary data on screen?
LOGICAL                      :: SumPrint                                        ! Print summary data to "*.fsm"?

!bjj start of proposed change v6.02d-bjj
!rmCHARACTER(99)                :: ADAMSFile                                       ! The name of the file containing ADAMS-specific data inputs.
!rmCHARACTER(99)                :: ADFile                                          ! The name of the AeroDyn input file.
!rmCHARACTER(99), ALLOCATABLE   :: BldFile  (:)                                    ! The names of the blade-data input files.
!rmCHARACTER(1024)              :: DirRoot                                         ! The name of the root file including the full path to the current working directory.
!rmCHARACTER(99)                :: DynBrkFi                                        ! The name of the dynamic generator brake input file.
!rmCHARACTER(99)                :: FTitle                                          ! The title line from the primary input file.
!rmCHARACTER(99)                :: FurlFile                                        ! The name of the furling-data input file.
!rmCHARACTER(99)                :: LinFile                                         ! The name of the file containing FAST linearization control input parameters.
!rmCHARACTER(99)                :: NoiseFile                                       ! The name of the file containing aerodynamic noise input parameters.
!rm!bjj rmCHARACTER(99)                :: OutFile                                         ! The name of the output file.
!rmCHARACTER(99)                :: PriFile   = 'primary.fst'                       ! The name of the primary input file.  Can be overwritten on command line.
!rm!bjj Start of proposed change vXX NWTC_Lib
!rm!rmCHARACTER( 4), PARAMETER     :: ProgName  = 'FAST'                              ! The name of this program.
!rm!rmCHARACTER(62)                :: ProgVer                                         ! The version of this program.
!rm!bjj End of proposed change vXX NWTC_Lib
!rmCHARACTER(99)                :: PtfmFile                                        ! The name of the platform-data input file.
!rmCHARACTER(99)                :: RootName                                        ! The root name of the input and output files.
!rmCHARACTER(99)                :: TwrFile                                         ! The name of the tower-data input file.
CHARACTER(1024)              :: ADAMSFile                                       ! The name of the file containing ADAMS-specific data inputs.
CHARACTER(1024)              :: ADFile                                          ! The name of the AeroDyn input file.
CHARACTER(1024), ALLOCATABLE :: BldFile  (:)                                    ! The names of the blade-data input files.
CHARACTER(1024)              :: DirRoot                                         ! The name of the root file including the full path to the current working directory.
CHARACTER(1024)              :: DynBrkFi                                        ! The name of the dynamic generator brake input file.
CHARACTER(1024)              :: FTitle                                          ! The title line from the primary input file.
CHARACTER(1024)              :: FurlFile                                        ! The name of the furling-data input file.
CHARACTER(1024)              :: LinFile                                         ! The name of the file containing FAST linearization control input parameters.
CHARACTER(1024)              :: NoiseFile                                       ! The name of the file containing aerodynamic noise input parameters.
CHARACTER(1024)              :: PriFile   = 'primary.fst'                       ! The name of the primary input file.  Can be overwritten on command line.
CHARACTER(1024)              :: PtfmFile                                        ! The name of the platform-data input file.
CHARACTER(1024)              :: RootName                                        ! The root name of the input and output files.
CHARACTER(1024)              :: TwrFile                                         ! The name of the tower-data input file.
!bjj end of proposed change v6.02d-bjj


END MODULE General
!=======================================================================
MODULE InitCond


   ! This MODULE stores input variables for initial conditions.


USE                             Precision


REAL(ReKi)                   :: Azimuth                                         ! Initial azimuth angle for blade 1.
REAL(ReKi), ALLOCATABLE      :: BlPitchInit(:)                                  ! Initial blade pitch angles at the start of the simulation.
REAL(ReKi)                   :: IPDefl                                          ! Initial in-plane blade-tip deflection.
REAL(ReKi)                   :: NacYaw                                          ! Initial or fixed nacelle-yaw angle.
REAL(ReKi)                   :: OoPDefl                                         ! Initial out-of-plane blade-tip displacement.
REAL(ReKi)                   :: PtfmHeave = 0.0                                 ! Initial or fixed vertical heave translational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: PtfmPitch = 0.0                                 ! Initial or fixed pitch tilt rotational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: PtfmRoll  = 0.0                                 ! Initial or fixed roll tilt rotational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: PtfmSurge = 0.0                                 ! Initial or fixed horizontal surge translational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: PtfmSway  = 0.0                                 ! Initial or fixed horizontal sway translational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: PtfmYaw   = 0.0                                 ! Initial or fixed yaw rotational displacement of platform. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: QAzimInit                                       ! Initial value of the internal generator azimuth DOF (Q(DOF_GeAz)).
REAL(ReKi)                   :: RotFurl   = 0.0                                 ! Initial or fixed rotor-furl angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RotSpeed                                        ! Initial or fixed rotor speed.
REAL(ReKi)                   :: TailFurl  = 0.0                                 ! Initial or fixed tail-furl angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TTDspFA                                         ! Initial fore-aft tower-top displacement.
REAL(ReKi)                   :: TTDspSS                                         ! Initial side-to-side tower-top displacement.
REAL(ReKi)                   :: TeetDefl  = 0.0                                 ! Initial or fixed teeter angle. (Initialized to zero b/c the 3-blader requires it to be zero)

!bjj chg: LOGICAL(1), ALLOCATABLE      :: DOF_FlagInit(:)                                 ! Array which stores initial values of the feature flags for each DOF (at the start of the simulation).
LOGICAL,    ALLOCATABLE      :: DOF_FlagInit(:)                                 ! Array which stores initial values of the feature flags for each DOF (at the start of the simulation).


END MODULE InitCond
!=======================================================================
MODULE Linear


   ! This MODULE stores variables for a FAST linearization analysis.


USE                             Precision


REAL(ReKi)                   :: AbsQDNorm = 0.0                                 ! 2-norm of the absolute difference between the velocites     of two consecutive periods.
REAL(ReKi)                   :: AbsQNorm  = 0.0                                 ! 2-norm of the absolute difference between the displacements of two consecutive periods.
REAL(ReKi)                   :: DelGenTrq = 0.0                                 ! Pertubation in generator torque using during FAST linearization (zero otherwise).
REAL(ReKi)                   :: DispTol                                         ! Convergence tolerance for the 2-norm of the absolute difference between the displacements of two consecutive periods (rad).
REAL(ReKi)                   :: Period                                          ! Steady state period of solution.
REAL(ReKi), ALLOCATABLE      :: QD2op    (:,:)                                  ! Periodic steady state operating accelerations.
REAL(ReKi), ALLOCATABLE      :: QDop     (:,:)                                  ! Periodic steady state operating velocities.
REAL(ReKi), ALLOCATABLE      :: Qop      (:,:)                                  ! Periodic steady state operating displacements.
REAL(ReKi)                   :: VelTol                                          ! Convergence tolerance for the 2-norm of the absolute difference between the velocities    of two consecutive periods (rad/s).

INTEGER(4)                   :: CntrlInpt(7)                                    ! List   of control inputs [1 to NInputs] {1: nacelle yaw angle, 2: nacelle yaw rate, 3: generator torque, 4: collective blade pitch, 5: individual pitch of blade 1, 6: individual pitch of blade 2, 7: individual pitch of blade 3 [unavailable for 2-bladed turbines]} (-) [unused if NInputs=0]
INTEGER(4)                   :: Disturbnc(7)                                    ! List   of input wind disturbances [1 to NDisturbs] {1: horizontal hub-height wind speed, 2: horizontal wind direction, 3: vertical wind speed, 4: horizontal wind shear, 5: vertical power law wind shear, 6: linear vertical wind shear, 7: horizontal hub-height wind gust} (-) [unused if NDisturbs=0]
INTEGER(4)                   :: Iteration = 0                                   ! Current iteration (number of periods to convergence)
INTEGER(4)                   :: MdlOrder                                        ! Order of output linearized model (1: 1st order A, B, Bd; 2: 2nd order M, C, K, F, Fd) (switch)
INTEGER(4)                   :: NAzimStep                                       ! Number of azimuth steps in periodic linearized model (-).
INTEGER(4)                   :: NDisturbs                                       ! Number of wind disturbances [0 to 7] (-)
INTEGER(4)                   :: NInputs                                         ! Number of control inputs [0 (none) or 1 to 4+NumBl] (-)
INTEGER(4)                   :: NStep                                           ! Number of time steps in one Period.
INTEGER(4)                   :: TrimCase                                        ! Trim case {1: find nacelle yaw, 2: find generator torque, 3: find collective blade pitch} (switch) [used only when CalcStdy=True and GenDOF=True]

!bjj chg: LOGICAL(1)                   :: CalcStdy                                        ! Calculate periodic steady state condition (False: linearize about zero) (switch).
LOGICAL                      :: CalcStdy                                        ! Calculate periodic steady state condition (False: linearize about zero) (switch).
!bjj chg: LOGICAL(1)                   :: IgnoreMOD = .FALSE.                             ! Ignore the use of function MOD in SUBROUTINE CalcOuts()?
LOGICAL                      :: IgnoreMOD = .FALSE.                             ! Ignore the use of function MOD in SUBROUTINE CalcOuts()?


END MODULE Linear
!=======================================================================
MODULE MassInert


   ! This MODULE stores input variables for turbine mass and inertias.


USE                             Precision


REAL(ReKi)                   :: AtfaIner                                        ! Inertia of tail boom about the tail-furl axis whose origin is the tail boom center of mass.
REAL(ReKi), ALLOCATABLE      :: BldCG     (:)                                   ! Blade center of mass wrt the blade root.
REAL(ReKi), ALLOCATABLE      :: BldMass   (:)                                   ! Blade masses
REAL(ReKi)                   :: BoomMass  = 0.0                                 ! Tail boom mass. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi), ALLOCATABLE      :: FirstMom  (:)                                   ! First mass moment of inertia of blades wrt the root.
REAL(ReKi)                   :: GenIner                                         ! Generator inertia about HSS.
REAL(ReKi)                   :: Hubg1Iner                                       ! Inertia of hub about g1-axis (rotor centerline).
REAL(ReKi)                   :: Hubg2Iner                                       ! Inertia of hub about g2-axis (transverse to the cyclinder and passing through its c.g.).
REAL(ReKi)                   :: HubIner                                         ! Hub inertia about teeter axis (2-blader) or rotor axis (3-blader).
REAL(ReKi)                   :: HubMass                                         ! Hub mass.
REAL(ReKi)                   :: Nacd2Iner                                       ! Inertia of nacelle about the d2-axis whose origin is the nacelle center of mass.
REAL(ReKi)                   :: NacMass                                         ! Nacelle mass.
REAL(ReKi)                   :: NacYIner                                        ! Nacelle yaw inertia.
REAL(ReKi)                   :: PtfmMass  = 0.0                                 ! Platform mass. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: PtfmPIner = 0.0                                 ! Platform inertia for pitch tilt rotation about the platform CM. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: PtfmRIner = 0.0                                 ! Platform inertia for roll tilt rotation about the platform CM. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: PtfmYIner = 0.0                                 ! Platform inertia for yaw rotation about the platform CM. (Initialized to zero b/c not all models will read in PtfmFile)
REAL(ReKi)                   :: RFrlIner  = 0.0                                 ! Rotor-furl inertia about rotor-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlMass  = 0.0                                 ! Rotor-furl mass. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RotIner                                         ! Inertia of rotor about its centerline.
REAL(ReKi)                   :: RotMass                                         ! Rotor mass (blades, tips, and hub)
REAL(ReKi)                   :: RrfaIner                                        ! Inertia of structure that furls with the rotor (not including rotor) about the rotor-furl axis whose origin is the center of mass of the structure that furls with the rotor (not including rotor).
REAL(ReKi), ALLOCATABLE      :: SecondMom (:)                                   ! Second mass moment of inertia of blades wrt the root.
REAL(ReKi)                   :: TFinMass  = 0.0                                 ! Tail fin mass. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlIner  = 0.0                                 ! Tail boom inertia about tail-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi), ALLOCATABLE      :: TipMass   (:)                                   ! Tip-brake masses.
REAL(ReKi)                   :: TotalMass                                       ! Mass of turbine + platform.
REAL(ReKi)                   :: TurbMass                                        ! Mass of turbine (tower + rotor + nacelle).
REAL(ReKi)                   :: TwrMass                                         ! Mass of tower.
REAL(ReKi)                   :: TwrTpMass                                       ! Tower-top mass (rotor + nacelle).
REAL(ReKi)                   :: YawBrMass                                       ! Yaw bearing mass.


END MODULE MassInert
!=======================================================================
MODULE Modes


   ! This MODULE stores variables for mode shapes and CONTAINs FUNCTION SHP.


USE                             Precision


!eab Start of proposed change.  v6.10b-eab  24-Jul-2008.
!eab Allow the mode shape to be described by a polynomial greater than sixth 
!eab   order. Do this by setting the order of the polynomial to a parameter.
INTEGER(4), PARAMETER        :: PolyOrd = 6                                     ! Order of the polynomial describing the mode shape.

!eab End of proposed change.  v6.10b-eab  24-Jul-2008.
REAL(ReKi), ALLOCATABLE      :: BldEdgSh (:,:)                                  ! Blade-edge-mode shape coefficients.
REAL(ReKi), ALLOCATABLE      :: BldFl1Sh (:,:)                                  ! Blade-flap-mode-1 shape coefficients.
REAL(ReKi), ALLOCATABLE      :: BldFl2Sh (:,:)                                  ! Blade-flap-mode-2 shape coefficients.
REAL(ReKi), ALLOCATABLE      :: FreqBE   (:,:,:)                                ! Blade edgewise natural frequencies (both w/ and w/o centrifugal stiffening)
REAL(ReKi), ALLOCATABLE      :: FreqBF   (:,:,:)                                ! Blade flapwise natural frequencies (both w/ and w/o centrifugal stiffening)
REAL(ReKi)                   :: FreqTFA  (2,2)                                  ! Computed fore-aft tower natural frequencies.
REAL(ReKi)                   :: FreqTSS  (2,2)                                  ! Computed side-to-side tower natural frequencies.
!eab Start of proposed change.  v6.10b-eab  24-Jul-2008.
!eab Allow shape coefficient arrays to exceed the sixth position so as to
!eab   accommodate all coefficients of a polynomial of any order.
!remove6.10b   REAL(ReKi)                   :: TwFAM1Sh (2:6)                                  ! Tower fore-aft mode-1 shape coefficients.
!remove6.10b   REAL(ReKi)                   :: TwFAM2Sh (2:6)                                  ! Tower fore-aft mode-2 shape coefficients.
!remove6.10b   REAL(ReKi)                   :: TwSSM1Sh (2:6)                                  ! Tower side-to-side mode-1 shape coefficients.
!remove6.10b   REAL(ReKi)                   :: TwSSM2Sh (2:6)                                  ! Tower side-to-side mode-2 shape coefficients.
REAL(ReKi)                   :: TwFAM1Sh (2:PolyOrd)                            ! Tower fore-aft mode-1 shape coefficients.
REAL(ReKi)                   :: TwFAM2Sh (2:PolyOrd)                            ! Tower fore-aft mode-2 shape coefficients.    
REAL(ReKi)                   :: TwSSM1Sh (2:PolyOrd)                            ! Tower side-to-side mode-1 shape coefficients.
REAL(ReKi)                   :: TwSSM2Sh (2:PolyOrd)                            ! Tower side-to-side mode-2 shape coefficients.
!eab End of proposed change.  v6.10b-eab  24-Jul-2008.

!bjj chg: LOGICAL(1)                   :: CalcBMode                                       ! T: calculate blade mode shapes internally, F: use blade mode shapes from the blade file.
LOGICAL                      :: CalcBMode                                       ! T: calculate blade mode shapes internally, F: use blade mode shapes from the blade file.
!bjj chg: LOGICAL(1), ALLOCATABLE      :: CalcBModes(:)                                   ! Holds CalcBMode for all of the blades.
LOGICAL,    ALLOCATABLE      :: CalcBModes(:)                                   ! Holds CalcBMode for all of the blades.
!bjj chg: LOGICAL(1)                   :: CalcTMode                                       ! T: calculate tower mode shapes internally, F: use tower mode shapes from the tower file.
LOGICAL                      :: CalcTMode                                       ! T: calculate tower mode shapes internally, F: use tower mode shapes from the tower file.


CONTAINS
!=======================================================================

   FUNCTION SHP(Fract, FlexL, ModShpAry, Deriv)

   ! SHP calculates the Derive-derivative of the shape function ModShpAry at Fract.
   ! NOTE: This function only works for Deriv = 0, 1, or 2.

      USE                           Precision


      IMPLICIT                      NONE


   ! Passed variables:

      REAL(ReKi), INTENT(IN )    :: FlexL          ! Length of flexible beam, (m)
      REAL(ReKi), INTENT(IN )    :: Fract          ! Fractional distance along flexible beam, 0<=Frac<=1
!eab Start of proposed change.  v6.10b-eab  24-Jul-2008.
!eab Allow the shape coefficient array to exceed the sixth position so as to
!eab   accommodate all coefficients of a polynomial of any order.
!remove6.10b      REAL(ReKi), INTENT(IN )    :: ModShpAry(2:6)                              ! Array holding mode shape coefficients
      REAL(ReKi), INTENT(IN )    :: ModShpAry(2:PolyOrd)                        ! Array holding mode shape coefficients
!eab End of proposed change.  v6.10b-eab  24-Jul-2008.
      REAL(ReKi)                 :: SHP            ! The shape function returned by this function.

      INTEGER(4), INTENT(IN )    :: Deriv          ! Which derivative to compute Deriv = 0 (regular function SHP), 1 (D(SHP)/DZ), 2 (D2(SHP)/DZ2)


   ! Lccal variables:

      INTEGER(4)                 :: CoefTmp        !Temporary coefficient
      INTEGER(4)                 :: I              !Counts through polynomial array.
      INTEGER(4)                 :: Swtch(0:2)     !Corresponds to which derivative to compute.  Sets all portions of the coefficient = 0 except those that are relevant.


      Swtch        = 0 !Initialize Swtch(:) to 0
      Swtch(Deriv) = 1
      SHP          = 0.0

!eab Start of proposed change.  v6.10b-eab  24-Jul-2008.
!eab Loop through all mode shape terms.
!remove6.10b      DO I = 2,6
      DO I = 2,PolyOrd
!eab End of proposed change.  v6.10b-eab  24-Jul-2008.
         CoefTmp = Swtch(0) + ( Swtch(1)*I ) + ( Swtch(2)*I*( I - 1 ) )

         IF ( (I == 2) .AND. (Deriv == 2) ) THEN
            SHP = ModShpAry(I)*CoefTmp/( FlexL**Deriv )
         ELSE
            SHP = SHP + ModShpAry(I)*CoefTmp*( Fract**( I - Deriv ) )/( FlexL**Deriv )
         ENDIF
      ENDDO !I

      RETURN
   END FUNCTION SHP
!=======================================================================
END MODULE Modes
!=======================================================================
MODULE NacelleYaw


   ! This MODULE stores variables for nacelle yaw.


USE                             Precision


REAL(ReKi)                   :: YawSpr                                          ! Nacelle-yaw spring constant.
REAL(ReKi)                   :: YawDamp                                         ! Nacelle-yaw constant.
REAL(ReKi)                   :: YawNeut                                         ! Neutral yaw position.
REAL(ReKi)                   :: YawRateNeut = 0.0                               ! Neutral yaw rate.


END MODULE NacelleYaw
!=======================================================================
MODULE Output


   ! This MODULE stores variables used for output.


USE                             Precision


   ! Defined TYPEs:

TYPE OutPar                                                                     ! User-defined type for output parameters
   CHARACTER(10)             :: Name                                            ! Name of the output parameter.
   CHARACTER(10)             :: Units                                           ! Units corresponding to the output parameter.
END TYPE OutPar


   ! Parameters:

   ! Indices for computing output channels:
   ! NOTE: These variables are declared in numerical order, instead of
   !       alphabetical order.
   ! NOTE: If an index for computing a NEW output channel is ever designated to
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Increase the upper limit for the number of blade and tower strain gage
!jmj   locations from 5 to 9 and add new output parameters for the local loads
!jmj   and motions at the additional strain gage locations:
!jmj Also, add an undocumented feature for outputting the incident wave
!jmj   elevation at the platform reference point and the incident wave
!jmj   kinematics at up to 9 nodes along the undeflected tower [not floating]
!jmj   or undisplaced platform [floating]:
!remove6.02a   !       be greater than 286, reDIMENSION array AllOuts() accordingly.
!jmj Start of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Add blade strain gage output parameters for the local loads and motions of
!jmj   blades 2 and 3:
!jmj Also, replace the hard-coded mooring line restoring calculation with a
!jmj   general purpose, quasi-static solution based on the analytical catenary
!jmj   cable equations with seabed interaction:
!remove6.02b   !       be greater than 389, reDIMENSION array AllOuts() accordingly.
   !       be greater than 533, reDIMENSION array AllOuts() accordingly.
!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Start of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Add blade strain gage output parameters for the local loads and motions of
!jmj   blades 2 and 3:
!jmj Also, replace the hard-coded mooring line restoring calculation with a
!jmj   general purpose, quasi-static solution based on the analytical catenary
!jmj   cable equations with seabed interaction:
!remove6.02b   ! NOTE: If an index ever becomes greater or equal to 500, the logic to
   ! NOTE: If an index ever becomes greater or equal to 1000, the logic to
!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.
   !       create ARRAY/1 in the FAST-to-ADAMS preprocessor will have to be
   !       changed.

   ! Time:

INTEGER(4), PARAMETER        :: Time                 =   0


   ! Blade 1 Tip Motions:

INTEGER(4), PARAMETER        :: TipDxc1              =   1
INTEGER(4), PARAMETER        :: TipDyc1              =   2
INTEGER(4), PARAMETER        :: TipDzc1              =   3
INTEGER(4), PARAMETER        :: TipDxb1              =   4
INTEGER(4), PARAMETER        :: TipDyb1              =   5
INTEGER(4), PARAMETER        :: TipALxb1             =   6
INTEGER(4), PARAMETER        :: TipALyb1             =   7
INTEGER(4), PARAMETER        :: TipALzb1             =   8
INTEGER(4), PARAMETER        :: TipRDxb1             =   9
INTEGER(4), PARAMETER        :: TipRDyb1             =  10
INTEGER(4), PARAMETER        :: TipRDzc1             =  11
INTEGER(4), PARAMETER        :: TipClrnc1            =  12


   ! Blade 1 Local Span Motions:

INTEGER(4), PARAMETER        :: Spn1ALxb1            =  13
INTEGER(4), PARAMETER        :: Spn1ALyb1            =  14
INTEGER(4), PARAMETER        :: Spn1ALzb1            =  15
INTEGER(4), PARAMETER        :: Spn2ALxb1            =  16
INTEGER(4), PARAMETER        :: Spn2ALyb1            =  17
INTEGER(4), PARAMETER        :: Spn2ALzb1            =  18
INTEGER(4), PARAMETER        :: Spn3ALxb1            =  19
INTEGER(4), PARAMETER        :: Spn3ALyb1            =  20
INTEGER(4), PARAMETER        :: Spn3ALzb1            =  21
INTEGER(4), PARAMETER        :: Spn4ALxb1            =  22
INTEGER(4), PARAMETER        :: Spn4ALyb1            =  23
INTEGER(4), PARAMETER        :: Spn4ALzb1            =  24
INTEGER(4), PARAMETER        :: Spn5ALxb1            =  25
INTEGER(4), PARAMETER        :: Spn5ALyb1            =  26
INTEGER(4), PARAMETER        :: Spn5ALzb1            =  27
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Increase the upper limit for the number of blade and tower strain gage
!jmj   locations from 5 to 9 and add new output parameters for the local loads
!jmj   and motions at the additional strain gage locations.  While doing this,
!jmj   renumber the remaining output PARAMETERs accordingly, without further
!jmj   documentation:
INTEGER(4), PARAMETER        :: Spn6ALxb1            =  28
INTEGER(4), PARAMETER        :: Spn6ALyb1            =  29
INTEGER(4), PARAMETER        :: Spn6ALzb1            =  30
INTEGER(4), PARAMETER        :: Spn7ALxb1            =  31
INTEGER(4), PARAMETER        :: Spn7ALyb1            =  32
INTEGER(4), PARAMETER        :: Spn7ALzb1            =  33
INTEGER(4), PARAMETER        :: Spn8ALxb1            =  34
INTEGER(4), PARAMETER        :: Spn8ALyb1            =  35
INTEGER(4), PARAMETER        :: Spn8ALzb1            =  36
INTEGER(4), PARAMETER        :: Spn9ALxb1            =  37
INTEGER(4), PARAMETER        :: Spn9ALyb1            =  38
INTEGER(4), PARAMETER        :: Spn9ALzb1            =  39
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.

   ! Blade 2 Tip Motions:

INTEGER(4), PARAMETER        :: TipDxc2              =  40
INTEGER(4), PARAMETER        :: TipDyc2              =  41
INTEGER(4), PARAMETER        :: TipDzc2              =  42
INTEGER(4), PARAMETER        :: TipDxb2              =  43
INTEGER(4), PARAMETER        :: TipDyb2              =  44
INTEGER(4), PARAMETER        :: TipALxb2             =  45
INTEGER(4), PARAMETER        :: TipALyb2             =  46
INTEGER(4), PARAMETER        :: TipALzb2             =  47
INTEGER(4), PARAMETER        :: TipRDxb2             =  48
INTEGER(4), PARAMETER        :: TipRDyb2             =  49
INTEGER(4), PARAMETER        :: TipRDzc2             =  50
INTEGER(4), PARAMETER        :: TipClrnc2            =  51


!jmj Start of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Add blade strain gage output parameters for the local loads and motions of
!jmj   blades 2 and 3:
!jmj While doing this, renumber the remaining output PARAMETERs accordingly,
!jmj   without further documentation:
   ! Blade 2 Local Span Motions:

INTEGER(4), PARAMETER        :: Spn1ALxb2            =  52
INTEGER(4), PARAMETER        :: Spn1ALyb2            =  53
INTEGER(4), PARAMETER        :: Spn1ALzb2            =  54
INTEGER(4), PARAMETER        :: Spn2ALxb2            =  55
INTEGER(4), PARAMETER        :: Spn2ALyb2            =  56
INTEGER(4), PARAMETER        :: Spn2ALzb2            =  57
INTEGER(4), PARAMETER        :: Spn3ALxb2            =  58
INTEGER(4), PARAMETER        :: Spn3ALyb2            =  59
INTEGER(4), PARAMETER        :: Spn3ALzb2            =  60
INTEGER(4), PARAMETER        :: Spn4ALxb2            =  61
INTEGER(4), PARAMETER        :: Spn4ALyb2            =  62
INTEGER(4), PARAMETER        :: Spn4ALzb2            =  63
INTEGER(4), PARAMETER        :: Spn5ALxb2            =  64
INTEGER(4), PARAMETER        :: Spn5ALyb2            =  65
INTEGER(4), PARAMETER        :: Spn5ALzb2            =  66
INTEGER(4), PARAMETER        :: Spn6ALxb2            =  67
INTEGER(4), PARAMETER        :: Spn6ALyb2            =  68
INTEGER(4), PARAMETER        :: Spn6ALzb2            =  69
INTEGER(4), PARAMETER        :: Spn7ALxb2            =  70
INTEGER(4), PARAMETER        :: Spn7ALyb2            =  71
INTEGER(4), PARAMETER        :: Spn7ALzb2            =  72
INTEGER(4), PARAMETER        :: Spn8ALxb2            =  73
INTEGER(4), PARAMETER        :: Spn8ALyb2            =  74
INTEGER(4), PARAMETER        :: Spn8ALzb2            =  75
INTEGER(4), PARAMETER        :: Spn9ALxb2            =  76
INTEGER(4), PARAMETER        :: Spn9ALyb2            =  77
INTEGER(4), PARAMETER        :: Spn9ALzb2            =  78


!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.
   ! Blade 3 Tip Motions:

INTEGER(4), PARAMETER        :: TipDxc3              =  79
INTEGER(4), PARAMETER        :: TipDyc3              =  80
INTEGER(4), PARAMETER        :: TipDzc3              =  81
INTEGER(4), PARAMETER        :: TipDxb3              =  82
INTEGER(4), PARAMETER        :: TipDyb3              =  83
INTEGER(4), PARAMETER        :: TipALxb3             =  84
INTEGER(4), PARAMETER        :: TipALyb3             =  85
INTEGER(4), PARAMETER        :: TipALzb3             =  86
INTEGER(4), PARAMETER        :: TipRDxb3             =  87
INTEGER(4), PARAMETER        :: TipRDyb3             =  88
INTEGER(4), PARAMETER        :: TipRDzc3             =  89
INTEGER(4), PARAMETER        :: TipClrnc3            =  90


!jmj Start of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Add blade strain gage output parameters for the local loads and motions of
!jmj   blades 2 and 3:
   ! Blade 3 Local Span Motions:

INTEGER(4), PARAMETER        :: Spn1ALxb3            =  91
INTEGER(4), PARAMETER        :: Spn1ALyb3            =  92
INTEGER(4), PARAMETER        :: Spn1ALzb3            =  93
INTEGER(4), PARAMETER        :: Spn2ALxb3            =  94
INTEGER(4), PARAMETER        :: Spn2ALyb3            =  95
INTEGER(4), PARAMETER        :: Spn2ALzb3            =  96
INTEGER(4), PARAMETER        :: Spn3ALxb3            =  97
INTEGER(4), PARAMETER        :: Spn3ALyb3            =  98
INTEGER(4), PARAMETER        :: Spn3ALzb3            =  99
INTEGER(4), PARAMETER        :: Spn4ALxb3            = 100
INTEGER(4), PARAMETER        :: Spn4ALyb3            = 101
INTEGER(4), PARAMETER        :: Spn4ALzb3            = 102
INTEGER(4), PARAMETER        :: Spn5ALxb3            = 103
INTEGER(4), PARAMETER        :: Spn5ALyb3            = 104
INTEGER(4), PARAMETER        :: Spn5ALzb3            = 105
INTEGER(4), PARAMETER        :: Spn6ALxb3            = 106
INTEGER(4), PARAMETER        :: Spn6ALyb3            = 107
INTEGER(4), PARAMETER        :: Spn6ALzb3            = 108
INTEGER(4), PARAMETER        :: Spn7ALxb3            = 109
INTEGER(4), PARAMETER        :: Spn7ALyb3            = 110
INTEGER(4), PARAMETER        :: Spn7ALzb3            = 111
INTEGER(4), PARAMETER        :: Spn8ALxb3            = 112
INTEGER(4), PARAMETER        :: Spn8ALyb3            = 113
INTEGER(4), PARAMETER        :: Spn8ALzb3            = 114
INTEGER(4), PARAMETER        :: Spn9ALxb3            = 115
INTEGER(4), PARAMETER        :: Spn9ALyb3            = 116
INTEGER(4), PARAMETER        :: Spn9ALzb3            = 117


!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.
   ! Blade Pitch Motions:

INTEGER(4), PARAMETER        :: PtchPMzc1            = 118
INTEGER(4), PARAMETER        :: PtchPMzc2            = 119
INTEGER(4), PARAMETER        :: PtchPMzc3            = 120


   ! Teeter Motions:

INTEGER(4), PARAMETER        :: TeetPya              = 121
INTEGER(4), PARAMETER        :: TeetVya              = 122
INTEGER(4), PARAMETER        :: TeetAya              = 123


   ! Shaft Motions:

INTEGER(4), PARAMETER        :: LSSTipPxa            = 124
INTEGER(4), PARAMETER        :: LSSTipVxa            = 125
INTEGER(4), PARAMETER        :: LSSTipAxa            = 126
INTEGER(4), PARAMETER        :: LSSGagPxa            = 127
INTEGER(4), PARAMETER        :: LSSGagVxa            = 128
INTEGER(4), PARAMETER        :: LSSGagAxa            = 129
INTEGER(4), PARAMETER        :: HSShftV              = 130
INTEGER(4), PARAMETER        :: HSShftA              = 131
INTEGER(4), PARAMETER        :: TipSpdRat            = 132


   ! Nacelle IMU Motions:

INTEGER(4), PARAMETER        :: NcIMUTVxs            = 133
INTEGER(4), PARAMETER        :: NcIMUTVys            = 134
INTEGER(4), PARAMETER        :: NcIMUTVzs            = 135
INTEGER(4), PARAMETER        :: NcIMUTAxs            = 136
INTEGER(4), PARAMETER        :: NcIMUTAys            = 137
INTEGER(4), PARAMETER        :: NcIMUTAzs            = 138
INTEGER(4), PARAMETER        :: NcIMURVxs            = 139
INTEGER(4), PARAMETER        :: NcIMURVys            = 140
INTEGER(4), PARAMETER        :: NcIMURVzs            = 141
INTEGER(4), PARAMETER        :: NcIMURAxs            = 142
INTEGER(4), PARAMETER        :: NcIMURAys            = 143
INTEGER(4), PARAMETER        :: NcIMURAzs            = 144


   ! Rotor-Furl Motions:

INTEGER(4), PARAMETER        :: RotFurlP             = 145
INTEGER(4), PARAMETER        :: RotFurlV             = 146
INTEGER(4), PARAMETER        :: RotFurlA             = 147


   ! Yaw Motions:

INTEGER(4), PARAMETER        :: YawPzn               = 148
INTEGER(4), PARAMETER        :: YawVzn               = 149
INTEGER(4), PARAMETER        :: YawAzn               = 150
INTEGER(4), PARAMETER        :: NacYawErr            = 151


   ! Tower-Top / Yaw Bearing Motions:

INTEGER(4), PARAMETER        :: YawBrTDxp            = 152
INTEGER(4), PARAMETER        :: YawBrTDyp            = 153
INTEGER(4), PARAMETER        :: YawBrTDzp            = 154
INTEGER(4), PARAMETER        :: YawBrTDxt            = 155
INTEGER(4), PARAMETER        :: YawBrTDyt            = 156
INTEGER(4), PARAMETER        :: YawBrTDzt            = 157
INTEGER(4), PARAMETER        :: YawBrTAxp            = 158
INTEGER(4), PARAMETER        :: YawBrTAyp            = 159
INTEGER(4), PARAMETER        :: YawBrTAzp            = 160
INTEGER(4), PARAMETER        :: YawBrRDxt            = 161
INTEGER(4), PARAMETER        :: YawBrRDyt            = 162
INTEGER(4), PARAMETER        :: YawBrRDzt            = 163
INTEGER(4), PARAMETER        :: YawBrRVxp            = 164
INTEGER(4), PARAMETER        :: YawBrRVyp            = 165
INTEGER(4), PARAMETER        :: YawBrRVzp            = 166
INTEGER(4), PARAMETER        :: YawBrRAxp            = 167
INTEGER(4), PARAMETER        :: YawBrRAyp            = 168
INTEGER(4), PARAMETER        :: YawBrRAzp            = 169


   ! Tail-Furl Motions:

INTEGER(4), PARAMETER        :: TailFurlP            = 170
INTEGER(4), PARAMETER        :: TailFurlV            = 171
INTEGER(4), PARAMETER        :: TailFurlA            = 172


   ! Local Tower Motions:

INTEGER(4), PARAMETER        :: TwHt1ALxt            = 173
INTEGER(4), PARAMETER        :: TwHt1ALyt            = 174
INTEGER(4), PARAMETER        :: TwHt1ALzt            = 175
INTEGER(4), PARAMETER        :: TwHt2ALxt            = 176
INTEGER(4), PARAMETER        :: TwHt2ALyt            = 177
INTEGER(4), PARAMETER        :: TwHt2ALzt            = 178
INTEGER(4), PARAMETER        :: TwHt3ALxt            = 179
INTEGER(4), PARAMETER        :: TwHt3ALyt            = 180
INTEGER(4), PARAMETER        :: TwHt3ALzt            = 181
INTEGER(4), PARAMETER        :: TwHt4ALxt            = 182
INTEGER(4), PARAMETER        :: TwHt4ALyt            = 183
INTEGER(4), PARAMETER        :: TwHt4ALzt            = 184
INTEGER(4), PARAMETER        :: TwHt5ALxt            = 185
INTEGER(4), PARAMETER        :: TwHt5ALyt            = 186
INTEGER(4), PARAMETER        :: TwHt5ALzt            = 187
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Increase the upper limit for the number of blade and tower strain gage
!jmj   locations from 5 to 9 and add new output parameters for the local loads
!jmj   and motions at the additional strain gage locations:
INTEGER(4), PARAMETER        :: TwHt6ALxt            = 188
INTEGER(4), PARAMETER        :: TwHt6ALyt            = 189
INTEGER(4), PARAMETER        :: TwHt6ALzt            = 190
INTEGER(4), PARAMETER        :: TwHt7ALxt            = 191
INTEGER(4), PARAMETER        :: TwHt7ALyt            = 192
INTEGER(4), PARAMETER        :: TwHt7ALzt            = 193
INTEGER(4), PARAMETER        :: TwHt8ALxt            = 194
INTEGER(4), PARAMETER        :: TwHt8ALyt            = 195
INTEGER(4), PARAMETER        :: TwHt8ALzt            = 196
INTEGER(4), PARAMETER        :: TwHt9ALxt            = 197
INTEGER(4), PARAMETER        :: TwHt9ALyt            = 198
INTEGER(4), PARAMETER        :: TwHt9ALzt            = 199
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.


   ! Platform Motions:

INTEGER(4), PARAMETER        :: PtfmTDxt             = 200
INTEGER(4), PARAMETER        :: PtfmTDyt             = 201
INTEGER(4), PARAMETER        :: PtfmTDzt             = 202
INTEGER(4), PARAMETER        :: PtfmTDxi             = 203
INTEGER(4), PARAMETER        :: PtfmTDyi             = 204
INTEGER(4), PARAMETER        :: PtfmTDzi             = 205
INTEGER(4), PARAMETER        :: PtfmTVxt             = 206
INTEGER(4), PARAMETER        :: PtfmTVyt             = 207
INTEGER(4), PARAMETER        :: PtfmTVzt             = 208
INTEGER(4), PARAMETER        :: PtfmTVxi             = 209
INTEGER(4), PARAMETER        :: PtfmTVyi             = 210
INTEGER(4), PARAMETER        :: PtfmTVzi             = 211
INTEGER(4), PARAMETER        :: PtfmTAxt             = 212
INTEGER(4), PARAMETER        :: PtfmTAyt             = 213
INTEGER(4), PARAMETER        :: PtfmTAzt             = 214
INTEGER(4), PARAMETER        :: PtfmTAxi             = 215
INTEGER(4), PARAMETER        :: PtfmTAyi             = 216
INTEGER(4), PARAMETER        :: PtfmTAzi             = 217
INTEGER(4), PARAMETER        :: PtfmRDxi             = 218
INTEGER(4), PARAMETER        :: PtfmRDyi             = 219
INTEGER(4), PARAMETER        :: PtfmRDzi             = 220
INTEGER(4), PARAMETER        :: PtfmRVxt             = 221
INTEGER(4), PARAMETER        :: PtfmRVyt             = 222
INTEGER(4), PARAMETER        :: PtfmRVzt             = 223
INTEGER(4), PARAMETER        :: PtfmRVxi             = 224
INTEGER(4), PARAMETER        :: PtfmRVyi             = 225
INTEGER(4), PARAMETER        :: PtfmRVzi             = 226
INTEGER(4), PARAMETER        :: PtfmRAxt             = 227
INTEGER(4), PARAMETER        :: PtfmRAyt             = 228
INTEGER(4), PARAMETER        :: PtfmRAzt             = 229
INTEGER(4), PARAMETER        :: PtfmRAxi             = 230
INTEGER(4), PARAMETER        :: PtfmRAyi             = 231
INTEGER(4), PARAMETER        :: PtfmRAzi             = 232



   ! Blade 1 Root Loads:

INTEGER(4), PARAMETER        :: RootFxc1             = 233
INTEGER(4), PARAMETER        :: RootFyc1             = 234
INTEGER(4), PARAMETER        :: RootFzc1             = 235
INTEGER(4), PARAMETER        :: RootFxb1             = 236
INTEGER(4), PARAMETER        :: RootFyb1             = 237
INTEGER(4), PARAMETER        :: RootMxc1             = 238
INTEGER(4), PARAMETER        :: RootMyc1             = 239
INTEGER(4), PARAMETER        :: RootMzc1             = 240
INTEGER(4), PARAMETER        :: RootMxb1             = 241
INTEGER(4), PARAMETER        :: RootMyb1             = 242


   ! Blade 1 Local Span Loads:

INTEGER(4), PARAMETER        :: Spn1MLxb1            = 243
INTEGER(4), PARAMETER        :: Spn1MLyb1            = 244
INTEGER(4), PARAMETER        :: Spn1MLzb1            = 245
INTEGER(4), PARAMETER        :: Spn2MLxb1            = 246
INTEGER(4), PARAMETER        :: Spn2MLyb1            = 247
INTEGER(4), PARAMETER        :: Spn2MLzb1            = 248
INTEGER(4), PARAMETER        :: Spn3MLxb1            = 249
INTEGER(4), PARAMETER        :: Spn3MLyb1            = 250
INTEGER(4), PARAMETER        :: Spn3MLzb1            = 251
INTEGER(4), PARAMETER        :: Spn4MLxb1            = 252
INTEGER(4), PARAMETER        :: Spn4MLyb1            = 253
INTEGER(4), PARAMETER        :: Spn4MLzb1            = 254
INTEGER(4), PARAMETER        :: Spn5MLxb1            = 255
INTEGER(4), PARAMETER        :: Spn5MLyb1            = 256
INTEGER(4), PARAMETER        :: Spn5MLzb1            = 257
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Increase the upper limit for the number of blade and tower strain gage
!jmj   locations from 5 to 9 and add new output parameters for the local loads
!jmj   and motions at the additional strain gage locations:
INTEGER(4), PARAMETER        :: Spn6MLxb1            = 258
INTEGER(4), PARAMETER        :: Spn6MLyb1            = 259
INTEGER(4), PARAMETER        :: Spn6MLzb1            = 260
INTEGER(4), PARAMETER        :: Spn7MLxb1            = 261
INTEGER(4), PARAMETER        :: Spn7MLyb1            = 262
INTEGER(4), PARAMETER        :: Spn7MLzb1            = 263
INTEGER(4), PARAMETER        :: Spn8MLxb1            = 264
INTEGER(4), PARAMETER        :: Spn8MLyb1            = 265
INTEGER(4), PARAMETER        :: Spn8MLzb1            = 266
INTEGER(4), PARAMETER        :: Spn9MLxb1            = 267
INTEGER(4), PARAMETER        :: Spn9MLyb1            = 268
INTEGER(4), PARAMETER        :: Spn9MLzb1            = 269
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.


   ! Blade 2 Root Loads:

INTEGER(4), PARAMETER        :: RootFxc2             = 270
INTEGER(4), PARAMETER        :: RootFyc2             = 271
INTEGER(4), PARAMETER        :: RootFzc2             = 272
INTEGER(4), PARAMETER        :: RootFxb2             = 273
INTEGER(4), PARAMETER        :: RootFyb2             = 274
INTEGER(4), PARAMETER        :: RootMxc2             = 275
INTEGER(4), PARAMETER        :: RootMyc2             = 276
INTEGER(4), PARAMETER        :: RootMzc2             = 277
INTEGER(4), PARAMETER        :: RootMxb2             = 278
INTEGER(4), PARAMETER        :: RootMyb2             = 279


!jmj Start of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Add blade strain gage output parameters for the local loads and motions of
!jmj   blades 2 and 3:
   ! Blade 2 Local Span Loads:

INTEGER(4), PARAMETER        :: Spn1MLxb2            = 280
INTEGER(4), PARAMETER        :: Spn1MLyb2            = 281
INTEGER(4), PARAMETER        :: Spn1MLzb2            = 282
INTEGER(4), PARAMETER        :: Spn2MLxb2            = 283
INTEGER(4), PARAMETER        :: Spn2MLyb2            = 284
INTEGER(4), PARAMETER        :: Spn2MLzb2            = 285
INTEGER(4), PARAMETER        :: Spn3MLxb2            = 286
INTEGER(4), PARAMETER        :: Spn3MLyb2            = 287
INTEGER(4), PARAMETER        :: Spn3MLzb2            = 288
INTEGER(4), PARAMETER        :: Spn4MLxb2            = 289
INTEGER(4), PARAMETER        :: Spn4MLyb2            = 290
INTEGER(4), PARAMETER        :: Spn4MLzb2            = 291
INTEGER(4), PARAMETER        :: Spn5MLxb2            = 292
INTEGER(4), PARAMETER        :: Spn5MLyb2            = 293
INTEGER(4), PARAMETER        :: Spn5MLzb2            = 294
INTEGER(4), PARAMETER        :: Spn6MLxb2            = 295
INTEGER(4), PARAMETER        :: Spn6MLyb2            = 296
INTEGER(4), PARAMETER        :: Spn6MLzb2            = 297
INTEGER(4), PARAMETER        :: Spn7MLxb2            = 298
INTEGER(4), PARAMETER        :: Spn7MLyb2            = 299
INTEGER(4), PARAMETER        :: Spn7MLzb2            = 300
INTEGER(4), PARAMETER        :: Spn8MLxb2            = 301
INTEGER(4), PARAMETER        :: Spn8MLyb2            = 302
INTEGER(4), PARAMETER        :: Spn8MLzb2            = 303
INTEGER(4), PARAMETER        :: Spn9MLxb2            = 304
INTEGER(4), PARAMETER        :: Spn9MLyb2            = 305
INTEGER(4), PARAMETER        :: Spn9MLzb2            = 306


!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.
   ! Blade 3 Root Loads:

INTEGER(4), PARAMETER        :: RootFxc3             = 307
INTEGER(4), PARAMETER        :: RootFyc3             = 308
INTEGER(4), PARAMETER        :: RootFzc3             = 309
INTEGER(4), PARAMETER        :: RootFxb3             = 310
INTEGER(4), PARAMETER        :: RootFyb3             = 311
INTEGER(4), PARAMETER        :: RootMxc3             = 312
INTEGER(4), PARAMETER        :: RootMyc3             = 313
INTEGER(4), PARAMETER        :: RootMzc3             = 314
INTEGER(4), PARAMETER        :: RootMxb3             = 315
INTEGER(4), PARAMETER        :: RootMyb3             = 316


!jmj Start of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Add blade strain gage output parameters for the local loads and motions of
!jmj   blades 2 and 3:
   ! Blade 3 Local Span Loads:

INTEGER(4), PARAMETER        :: Spn1MLxb3            = 317
INTEGER(4), PARAMETER        :: Spn1MLyb3            = 318
INTEGER(4), PARAMETER        :: Spn1MLzb3            = 319
INTEGER(4), PARAMETER        :: Spn2MLxb3            = 320
INTEGER(4), PARAMETER        :: Spn2MLyb3            = 321
INTEGER(4), PARAMETER        :: Spn2MLzb3            = 322
INTEGER(4), PARAMETER        :: Spn3MLxb3            = 323
INTEGER(4), PARAMETER        :: Spn3MLyb3            = 324
INTEGER(4), PARAMETER        :: Spn3MLzb3            = 325
INTEGER(4), PARAMETER        :: Spn4MLxb3            = 326
INTEGER(4), PARAMETER        :: Spn4MLyb3            = 327
INTEGER(4), PARAMETER        :: Spn4MLzb3            = 328
INTEGER(4), PARAMETER        :: Spn5MLxb3            = 329
INTEGER(4), PARAMETER        :: Spn5MLyb3            = 330
INTEGER(4), PARAMETER        :: Spn5MLzb3            = 331
INTEGER(4), PARAMETER        :: Spn6MLxb3            = 332
INTEGER(4), PARAMETER        :: Spn6MLyb3            = 333
INTEGER(4), PARAMETER        :: Spn6MLzb3            = 334
INTEGER(4), PARAMETER        :: Spn7MLxb3            = 335
INTEGER(4), PARAMETER        :: Spn7MLyb3            = 336
INTEGER(4), PARAMETER        :: Spn7MLzb3            = 337
INTEGER(4), PARAMETER        :: Spn8MLxb3            = 338
INTEGER(4), PARAMETER        :: Spn8MLyb3            = 339
INTEGER(4), PARAMETER        :: Spn8MLzb3            = 340
INTEGER(4), PARAMETER        :: Spn9MLxb3            = 341
INTEGER(4), PARAMETER        :: Spn9MLyb3            = 342
INTEGER(4), PARAMETER        :: Spn9MLzb3            = 343


!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.
   ! Hub and Rotor Loads:

INTEGER(4), PARAMETER        :: LSShftFxa            = 344
INTEGER(4), PARAMETER        :: LSShftFya            = 345
INTEGER(4), PARAMETER        :: LSShftFza            = 346
INTEGER(4), PARAMETER        :: LSShftFys            = 347
INTEGER(4), PARAMETER        :: LSShftFzs            = 348
INTEGER(4), PARAMETER        :: LSShftMxa            = 349
INTEGER(4), PARAMETER        :: LSSTipMya            = 350
INTEGER(4), PARAMETER        :: LSSTipMza            = 351
INTEGER(4), PARAMETER        :: LSSTipMys            = 352
INTEGER(4), PARAMETER        :: LSSTipMzs            = 353
INTEGER(4), PARAMETER        :: CThrstAzm            = 354
INTEGER(4), PARAMETER        :: CThrstRad            = 355
INTEGER(4), PARAMETER        :: RotPwr               = 356
INTEGER(4), PARAMETER        :: RotCq                = 357
INTEGER(4), PARAMETER        :: RotCp                = 358
INTEGER(4), PARAMETER        :: RotCt                = 359


   ! Shaft Strain Gage Loads:

INTEGER(4), PARAMETER        :: LSSGagMya            = 360
INTEGER(4), PARAMETER        :: LSSGagMza            = 361
INTEGER(4), PARAMETER        :: LSSGagMys            = 362
INTEGER(4), PARAMETER        :: LSSGagMzs            = 363


   ! Generator and High-Speed Shaft Loads:

INTEGER(4), PARAMETER        :: HSShftTq             = 364
INTEGER(4), PARAMETER        :: HSShftPwr            = 365
INTEGER(4), PARAMETER        :: HSShftCq             = 366
INTEGER(4), PARAMETER        :: HSShftCp             = 367
INTEGER(4), PARAMETER        :: HSSBrTq              = 368
INTEGER(4), PARAMETER        :: GenTq                = 369
INTEGER(4), PARAMETER        :: GenPwr               = 370
INTEGER(4), PARAMETER        :: GenCq                = 371
INTEGER(4), PARAMETER        :: GenCp                = 372


   ! Rotor-Furl Bearing Loads:

INTEGER(4), PARAMETER        :: RFrlBrM              = 373


   ! Tower-Top / Yaw Bearing Loads:

INTEGER(4), PARAMETER        :: YawBrFxn             = 374
INTEGER(4), PARAMETER        :: YawBrFyn             = 375
INTEGER(4), PARAMETER        :: YawBrFzn             = 376
INTEGER(4), PARAMETER        :: YawBrFxp             = 377
INTEGER(4), PARAMETER        :: YawBrFyp             = 378
INTEGER(4), PARAMETER        :: YawBrMxn             = 379
INTEGER(4), PARAMETER        :: YawBrMyn             = 380
INTEGER(4), PARAMETER        :: YawBrMzn             = 381
INTEGER(4), PARAMETER        :: YawBrMxp             = 382
INTEGER(4), PARAMETER        :: YawBrMyp             = 383


   ! Tower Base Loads:

INTEGER(4), PARAMETER        :: TwrBsFxt             = 384
INTEGER(4), PARAMETER        :: TwrBsFyt             = 385
INTEGER(4), PARAMETER        :: TwrBsFzt             = 386
INTEGER(4), PARAMETER        :: TwrBsMxt             = 387
INTEGER(4), PARAMETER        :: TwrBsMyt             = 388
INTEGER(4), PARAMETER        :: TwrBsMzt             = 389


   ! Tail-Furl Bearing Loads:

INTEGER(4), PARAMETER        :: TFrlBrM              = 390


   ! Local Tower Loads:

INTEGER(4), PARAMETER        :: TwHt1MLxt            = 391
INTEGER(4), PARAMETER        :: TwHt1MLyt            = 392
INTEGER(4), PARAMETER        :: TwHt1MLzt            = 393
INTEGER(4), PARAMETER        :: TwHt2MLxt            = 394
INTEGER(4), PARAMETER        :: TwHt2MLyt            = 395
INTEGER(4), PARAMETER        :: TwHt2MLzt            = 396
INTEGER(4), PARAMETER        :: TwHt3MLxt            = 397
INTEGER(4), PARAMETER        :: TwHt3MLyt            = 398
INTEGER(4), PARAMETER        :: TwHt3MLzt            = 399
INTEGER(4), PARAMETER        :: TwHt4MLxt            = 400
INTEGER(4), PARAMETER        :: TwHt4MLyt            = 401
INTEGER(4), PARAMETER        :: TwHt4MLzt            = 402
INTEGER(4), PARAMETER        :: TwHt5MLxt            = 403
INTEGER(4), PARAMETER        :: TwHt5MLyt            = 404
INTEGER(4), PARAMETER        :: TwHt5MLzt            = 405
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Increase the upper limit for the number of blade and tower strain gage
!jmj   locations from 5 to 9 and add new output parameters for the local loads
!jmj   and motions at the additional strain gage locations:
INTEGER(4), PARAMETER        :: TwHt6MLxt            = 406
INTEGER(4), PARAMETER        :: TwHt6MLyt            = 407
INTEGER(4), PARAMETER        :: TwHt6MLzt            = 408
INTEGER(4), PARAMETER        :: TwHt7MLxt            = 409
INTEGER(4), PARAMETER        :: TwHt7MLyt            = 410
INTEGER(4), PARAMETER        :: TwHt7MLzt            = 411
INTEGER(4), PARAMETER        :: TwHt8MLxt            = 412
INTEGER(4), PARAMETER        :: TwHt8MLyt            = 413
INTEGER(4), PARAMETER        :: TwHt8MLzt            = 414
INTEGER(4), PARAMETER        :: TwHt9MLxt            = 415
INTEGER(4), PARAMETER        :: TwHt9MLyt            = 416
INTEGER(4), PARAMETER        :: TwHt9MLzt            = 417
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.


   ! Platform Loads:

INTEGER(4), PARAMETER        :: PtfmFxt              = 418
INTEGER(4), PARAMETER        :: PtfmFyt              = 419
INTEGER(4), PARAMETER        :: PtfmFzt              = 420
INTEGER(4), PARAMETER        :: PtfmFxi              = 421
INTEGER(4), PARAMETER        :: PtfmFyi              = 422
INTEGER(4), PARAMETER        :: PtfmFzi              = 423
INTEGER(4), PARAMETER        :: PtfmMxt              = 424
INTEGER(4), PARAMETER        :: PtfmMyt              = 425
INTEGER(4), PARAMETER        :: PtfmMzt              = 426
INTEGER(4), PARAMETER        :: PtfmMxi              = 427
INTEGER(4), PARAMETER        :: PtfmMyi              = 428
INTEGER(4), PARAMETER        :: PtfmMzi              = 429


!jmj Start of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Replace the hard-coded mooring line restoring calculation with a general
!jmj   purpose, quasi-static solution based on the analytical catenary cable
!jmj   equations with seabed interaction:
   ! Mooring Line Loads:

INTEGER(4), PARAMETER        :: Fair1Ten             = 430
INTEGER(4), PARAMETER        :: Fair1Ang             = 431
INTEGER(4), PARAMETER        :: Anch1Ten             = 432
INTEGER(4), PARAMETER        :: Anch1Ang             = 433
INTEGER(4), PARAMETER        :: Fair2Ten             = 434
INTEGER(4), PARAMETER        :: Fair2Ang             = 435
INTEGER(4), PARAMETER        :: Anch2Ten             = 436
INTEGER(4), PARAMETER        :: Anch2Ang             = 437
INTEGER(4), PARAMETER        :: Fair3Ten             = 438
INTEGER(4), PARAMETER        :: Fair3Ang             = 439
INTEGER(4), PARAMETER        :: Anch3Ten             = 440
INTEGER(4), PARAMETER        :: Anch3Ang             = 441
INTEGER(4), PARAMETER        :: Fair4Ten             = 442
INTEGER(4), PARAMETER        :: Fair4Ang             = 443
INTEGER(4), PARAMETER        :: Anch4Ten             = 444
INTEGER(4), PARAMETER        :: Anch4Ang             = 445
INTEGER(4), PARAMETER        :: Fair5Ten             = 446
INTEGER(4), PARAMETER        :: Fair5Ang             = 447
INTEGER(4), PARAMETER        :: Anch5Ten             = 448
INTEGER(4), PARAMETER        :: Anch5Ang             = 449
INTEGER(4), PARAMETER        :: Fair6Ten             = 450
INTEGER(4), PARAMETER        :: Fair6Ang             = 451
INTEGER(4), PARAMETER        :: Anch6Ten             = 452
INTEGER(4), PARAMETER        :: Anch6Ang             = 453
INTEGER(4), PARAMETER        :: Fair7Ten             = 454
INTEGER(4), PARAMETER        :: Fair7Ang             = 455
INTEGER(4), PARAMETER        :: Anch7Ten             = 456
INTEGER(4), PARAMETER        :: Anch7Ang             = 457
INTEGER(4), PARAMETER        :: Fair8Ten             = 458
INTEGER(4), PARAMETER        :: Fair8Ang             = 459
INTEGER(4), PARAMETER        :: Anch8Ten             = 460
INTEGER(4), PARAMETER        :: Anch8Ang             = 461
INTEGER(4), PARAMETER        :: Fair9Ten             = 462
INTEGER(4), PARAMETER        :: Fair9Ang             = 463
INTEGER(4), PARAMETER        :: Anch9Ten             = 464
INTEGER(4), PARAMETER        :: Anch9Ang             = 465


!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.

   ! Wind Motions:

INTEGER(4), PARAMETER        :: WindVxi              = 466
INTEGER(4), PARAMETER        :: WindVyi              = 467
INTEGER(4), PARAMETER        :: WindVzi              = 468
INTEGER(4), PARAMETER        :: TotWindV             = 469
INTEGER(4), PARAMETER        :: HorWindV             = 470
INTEGER(4), PARAMETER        :: HorWndDir            = 471
INTEGER(4), PARAMETER        :: VerWndDir            = 472



   ! Tail Fin Element Aerodynamics:

INTEGER(4), PARAMETER        :: TFinAlpha            = 473
INTEGER(4), PARAMETER        :: TFinCLift            = 474
INTEGER(4), PARAMETER        :: TFinCDrag            = 475
INTEGER(4), PARAMETER        :: TFinDnPrs            = 476
INTEGER(4), PARAMETER        :: TFinCPFx             = 477
INTEGER(4), PARAMETER        :: TFinCPFy             = 478



!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for outputting the incident wave elevation at
!jmj   the platform reference point and the incident wave kinematics at up to 9
!jmj   nodes along the undeflected tower [not floating] or undisplaced platform
!jmj   [floating]:
   ! Wave Motions:

INTEGER(4), PARAMETER        :: WaveElev             = 479

INTEGER(4), PARAMETER        :: Wave1Vxi             = 480
INTEGER(4), PARAMETER        :: Wave1Vyi             = 481
INTEGER(4), PARAMETER        :: Wave1Vzi             = 482
INTEGER(4), PARAMETER        :: Wave1Axi             = 483
INTEGER(4), PARAMETER        :: Wave1Ayi             = 484
INTEGER(4), PARAMETER        :: Wave1Azi             = 485
INTEGER(4), PARAMETER        :: Wave2Vxi             = 486
INTEGER(4), PARAMETER        :: Wave2Vyi             = 487
INTEGER(4), PARAMETER        :: Wave2Vzi             = 488
INTEGER(4), PARAMETER        :: Wave2Axi             = 489
INTEGER(4), PARAMETER        :: Wave2Ayi             = 490
INTEGER(4), PARAMETER        :: Wave2Azi             = 491
INTEGER(4), PARAMETER        :: Wave3Vxi             = 492
INTEGER(4), PARAMETER        :: Wave3Vyi             = 493
INTEGER(4), PARAMETER        :: Wave3Vzi             = 494
INTEGER(4), PARAMETER        :: Wave3Axi             = 495
INTEGER(4), PARAMETER        :: Wave3Ayi             = 496
INTEGER(4), PARAMETER        :: Wave3Azi             = 497
INTEGER(4), PARAMETER        :: Wave4Vxi             = 498
INTEGER(4), PARAMETER        :: Wave4Vyi             = 499
INTEGER(4), PARAMETER        :: Wave4Vzi             = 500
INTEGER(4), PARAMETER        :: Wave4Axi             = 501
INTEGER(4), PARAMETER        :: Wave4Ayi             = 502
INTEGER(4), PARAMETER        :: Wave4Azi             = 503
INTEGER(4), PARAMETER        :: Wave5Vxi             = 504
INTEGER(4), PARAMETER        :: Wave5Vyi             = 505
INTEGER(4), PARAMETER        :: Wave5Vzi             = 506
INTEGER(4), PARAMETER        :: Wave5Axi             = 507
INTEGER(4), PARAMETER        :: Wave5Ayi             = 508
INTEGER(4), PARAMETER        :: Wave5Azi             = 509
INTEGER(4), PARAMETER        :: Wave6Vxi             = 510
INTEGER(4), PARAMETER        :: Wave6Vyi             = 511
INTEGER(4), PARAMETER        :: Wave6Vzi             = 512
INTEGER(4), PARAMETER        :: Wave6Axi             = 513
INTEGER(4), PARAMETER        :: Wave6Ayi             = 514
INTEGER(4), PARAMETER        :: Wave6Azi             = 515
INTEGER(4), PARAMETER        :: Wave7Vxi             = 516
INTEGER(4), PARAMETER        :: Wave7Vyi             = 517
INTEGER(4), PARAMETER        :: Wave7Vzi             = 518
INTEGER(4), PARAMETER        :: Wave7Axi             = 519
INTEGER(4), PARAMETER        :: Wave7Ayi             = 520
INTEGER(4), PARAMETER        :: Wave7Azi             = 521
INTEGER(4), PARAMETER        :: Wave8Vxi             = 522
INTEGER(4), PARAMETER        :: Wave8Vyi             = 523
INTEGER(4), PARAMETER        :: Wave8Vzi             = 524
INTEGER(4), PARAMETER        :: Wave8Axi             = 525
INTEGER(4), PARAMETER        :: Wave8Ayi             = 526
INTEGER(4), PARAMETER        :: Wave8Azi             = 527
INTEGER(4), PARAMETER        :: Wave9Vxi             = 528
INTEGER(4), PARAMETER        :: Wave9Vyi             = 529
INTEGER(4), PARAMETER        :: Wave9Vzi             = 530
INTEGER(4), PARAMETER        :: Wave9Axi             = 531
INTEGER(4), PARAMETER        :: Wave9Ayi             = 532
INTEGER(4), PARAMETER        :: Wave9Azi             = 533



!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.


INTEGER(4), PARAMETER        :: MaxOutPts            = 200                      ! The maximum number of output channels which can be outputted by the code.


!bjj rm NWTC_Lib: CHARACTER(1), PARAMETER      :: Tab                  = CHAR( 9 )                ! The tab character.


   ! Regular variables:

! SEE NOTE ABOVE FOR SIZE (DIMENSION) OF THE VARIABLE BELOW:
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Increase the upper limit for the number of blade and tower strain gage
!jmj   locations from 5 to 9 and add new output parameters for the local loads
!jmj   and motions at the additional strain gage locations:
!jmj Also, add an undocumented feature for outputting the incident wave
!jmj   elevation at the platform reference point and the incident wave
!jmj   kinematics at up to 9 nodes along the undeflected tower [not floating]
!jmj   or undisplaced platform [floating]:
!remove6.02aREAL(ReKi)                   :: AllOuts  (0:286)                                ! An array holding the value of all of the calculated (not selected) output channels.
!jmj Start of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Add blade strain gage output parameters for the local loads and motions of
!jmj   blades 2 and 3:
!jmj Also, replace the hard-coded mooring line restoring calculation with a
!jmj   general purpose, quasi-static solution based on the analytical catenary
!jmj   cable equations with seabed interaction:
!remove6.02bREAL(ReKi)                   :: AllOuts  (0:389)                                ! An array holding the value of all of the calculated (not selected) output channels.
REAL(ReKi)                   :: AllOuts  (0:533)                                ! An array holding the value of all of the calculated (not selected) output channels.
!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
! SEE NOTE ABOVE FOR SIZE (DIMENSION) OF THE PREVIOUS VARIABLE:
REAL(ReKi), ALLOCATABLE      :: LinAccES (:,:,:)                                ! Total linear acceleration of a point on a   blade (point S) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: LinAccET (:,:)                                  ! Total linear acceleration of a point on the tower (point T) in the inertia frame (body E for earth).
!jmj Start of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Replace the hard-coded mooring line restoring calculation with a general
!jmj   purpose, quasi-static solution based on the analytical catenary cable
!jmj   equations with seabed interaction:
REAL(ReKi), ALLOCATABLE      :: LSNodes  (:,:)                                  ! Unstretched arc distance along mooring line from anchor to each node where the line position and tension can be output (meters).
!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.
REAL(ReKi), ALLOCATABLE      :: FrcS0B   (:,:)                                  ! Total force at the blade root (point S(0)) due to the blade.
REAL(ReKi), ALLOCATABLE      :: FTHydro  (:,:)                                  ! Total hydrodynamic force per unit length acting on the tower at point T.
REAL(ReKi), ALLOCATABLE      :: MFHydro  (:,:)                                  ! Total hydrodynamic moment per unit length acting on a tower element (body F) at point T.
REAL(ReKi), ALLOCATABLE      :: MomH0B   (:,:)                                  ! Total moment at the hub (body H) / blade root (point S(0)) due to the blade.
REAL(ReKi)                   :: NcIMUxn                                         ! Downwind distance from the tower-top to the nacelle IMU.
REAL(ReKi)                   :: NcIMUyn                                         ! Lateral  distance from the tower-top to the nacelle IMU.
REAL(ReKi)                   :: NcIMUzn                                         ! Vertical distance from the tower-top to the nacelle IMU.
REAL(ReKi), ALLOCATABLE      :: OutData  (:)                                    ! Array to contain the output data.
REAL(ReKi)                   :: ShftGagL                                        ! Distance from hub or teeter pin to shaft strain gages.
REAL(ReKi)                   :: SttsTime                                        ! Amount of time between screen status messages (sec).
REAL(ReKi)                   :: TStart                                          ! Time to begin tabular output.
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for outputting the incident wave elevation at
!jmj   the platform reference point and the incident wave kinematics at up to 9
!jmj   nodes along the undeflected tower [not floating] or undisplaced platform
!jmj   [floating]:
REAL(ReKi)                   :: WaveElevxi(1) = (/ 0.0 /)                       ! xi-coordinates for points where the incident wave elevation can be output (meters).
REAL(ReKi)                   :: WaveElevyi(1) = (/ 0.0 /)                       ! yi-coordinates for points where the incident wave elevation can be output (meters).
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.

!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Increase the upper limit for the number of blade and tower strain gage
!jmj   locations from 5 to 9 and add new output parameters for the local loads
!jmj   and motions at the additional strain gage locations:
!remove6.02aINTEGER(4)                   :: BldGagNd (5)                                    ! Nodes closest to the blade strain gages.
INTEGER(4)                   :: BldGagNd (9)                                    ! Nodes closest to the blade strain gages.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
INTEGER(4)                   :: DecFact                                         ! Decimation factor for tabular output.
!jmj Start of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Replace the hard-coded mooring line restoring calculation with a general
!jmj   purpose, quasi-static solution based on the analytical catenary cable
!jmj   equations with seabed interaction:
!JASON: ADD OUTPUTS FOR THE MOORING LINE POSITION AND EFFECTIVE TENSION AT EACH NODE.  USE NAMES SUCH AS: Ln#Nd#Pxi, Ln#Nd#Pyi, Ln#Nd#Pzi, Ln#Nd#Ten WHERE # REPRESENTS THE LINE NUMBER OR NODE NUMBER!!!
INTEGER(4)                   :: LineNodes    = 0                                ! Number of nodes per line where the mooring line position and tension can be output (-).
!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.
INTEGER(4)                   :: NBlGages                                        ! Number of blade strain gages.
INTEGER(4)                   :: NTwGages                                        ! Number of tower strain gages.
INTEGER(4)                   :: NumOuts                                         ! Number of parameters in the output list.
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for outputting the incident wave elevation at
!jmj   the platform reference point and the incident wave kinematics at up to 9
!jmj   nodes along the undeflected tower [not floating] or undisplaced platform
!jmj   [floating]:
INTEGER(4)                   :: NWaveElev    = 1                                ! Number of points where the incident wave elevation  can be output (-).
INTEGER(4)                   :: NWaveKin     = 0                                ! Number of points where the incident wave kinematics can be output (-).
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
INTEGER(4), ALLOCATABLE      :: OutInd   (:)                                    ! Array of indices for OutData(:)
INTEGER(4), ALLOCATABLE      :: OutSign  (:)                                    ! Sign of output channel (+1 = normal, -1 = reversed).
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Increase the upper limit for the number of blade and tower strain gage
!jmj   locations from 5 to 9 and add new output parameters for the local loads
!jmj   and motions at the additional strain gage locations:
!jmj Also, add an undocumented feature for outputting the incident wave
!jmj   elevation at the platform reference point and the incident wave
!jmj   kinematics at up to 9 nodes along the undeflected tower [not floating]
!jmj   or undisplaced platform [floating]:
!remove6.02aINTEGER(4)                   :: TwrGagNd (5)                                    ! Nodes closest to the tower strain gages.
INTEGER(4)                   :: TwrGagNd (9)                                    ! Nodes closest to the tower strain gages.
INTEGER(4)                   :: WaveKinNd(9)                                    ! List of tower [not floating] or platform [floating] nodes that have wave kinematics sensors.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.

!bjj chg: LOGICAL(1)                   :: Echo                                            ! Flag for echoing input to the echo file.
!bjj rm NWTC_Library: LOGICAL                      :: Echo                                            ! Flag for echoing input to the echo file.
LOGICAL                      :: WrEcho
!bjj chg: LOGICAL(1)                   :: TabDelim                                        ! Flag to cause tab-delimited output.
LOGICAL                      :: TabDelim                                        ! Flag to cause tab-delimited output.

CHARACTER(20)                :: OutFmt                                          ! Output format for tabular data.
CHARACTER(10)                :: OutList  (MaxOutPts)                            ! List of output parameters.

TYPE (OutPar), ALLOCATABLE   :: OutParam (:)                                    ! Names and units of all output parameters.


END MODULE Output
!=======================================================================
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Rename MODULE PlatformLd() to Platform():
!remove6.02aMODULE PlatformLd
MODULE Platform
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.


   ! This MODULE stores input variables for platform loading.


USE                             Precision

!jmj Start of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Replace the hard-coded mooring line restoring calculation with a general
!jmj   purpose, quasi-static solution based on the analytical catenary cable
!jmj   equations with seabed interaction:
REAL(ReKi), ALLOCATABLE      :: LAnchxi  (:)                                    ! xi-coordinate of each anchor   in the inertial frame        coordinate system.
REAL(ReKi), ALLOCATABLE      :: LAnchyi  (:)                                    ! yi-coordinate of each anchor   in the inertial frame        coordinate system.
REAL(ReKi), ALLOCATABLE      :: LAnchzi  (:)                                    ! zi-coordinate of each anchor   in the inertial frame        coordinate system.
REAL(ReKi), ALLOCATABLE      :: LDiam    (:)                                    ! Effective diameter of each mooring line for calculation of the line buoyancy.
REAL(ReKi), ALLOCATABLE      :: LEAStff  (:)                                    ! Extensional stiffness of each mooring line.
REAL(ReKi), ALLOCATABLE      :: LFairxt  (:)                                    ! xt-coordinate of each fairlead in the tower base / platform coordinate system.
REAL(ReKi), ALLOCATABLE      :: LFairyt  (:)                                    ! yt-coordinate of each fairlead in the tower base / platform coordinate system.
REAL(ReKi), ALLOCATABLE      :: LFairzt  (:)                                    ! zt-coordinate of each fairlead in the tower base / platform coordinate system.
REAL(ReKi), ALLOCATABLE      :: LMassDen (:)                                    ! Mass density of each mooring line.
REAL(ReKi), ALLOCATABLE      :: LSeabedCD(:)                                    ! Coefficient of seabed static friction drag of each mooring line (a negative value indicates no seabed).
REAL(ReKi), ALLOCATABLE      :: LTenTol  (:)                                    ! Convergence tolerance within Newton-Raphson iteration of each mooring line specified as a fraction of tension.
REAL(ReKi), ALLOCATABLE      :: LUnstrLen(:)                                    ! Unstretched length of each mooring line.
!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Start of proposed change.  v6.10a-jmj  21-Feb-2007.
!jmj Make sure the seabed GRAPHICS in FAST-to-ADAMS has a radius at least as
!jmj   large as the maximum mooring line anchor radius:
REAL(ReKi)                   :: MaxLRadAnch  = 0.0                              ! Maximum value of input array LRadAnch.
!jmj End of proposed change.  v6.10a-jmj  21-Feb-2007.

REAL(ReKi)                   :: PtfmAM (6,6) = 0.0                              ! Platform added mass matrix.
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for modeling the hydrodynamic loading and
!jmj   mooring system dynamics for floating wind turbines.  Do this by allowing
!jmj   a keyword in place of the integers 0 or 1 in input PtfmLdMod when
!jmj   PtfmModel = 3:
REAL(ReKi)                   :: PtfmCD                                          ! Effective platform normalized hydrodynamic viscous drag coefficient in calculation of viscous drag term from Morison's equation.
REAL(ReKi)                   :: PtfmDiam                                        ! Effective platform diameter in calculation of viscous drag term from Morison's equation.
REAL(ReKi)                   :: PtfmDraft                                       ! Effective platform draft    in calculation of viscous drag term from Morison's equation.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.

REAL(ReKi)                   :: PtfmFt   (6) = 0.0                              ! The surge/xi (1), sway/yi (2), and heave/zi (3)-components of the portion of the platform force at the platform reference (point Z) and the roll/xi (4), pitch/yi (5), and yaw/zi (6)-components of the portion of the platform moment acting at the platform (body X) / platform reference (point Z) associated with everything but the QD2T()'s.

!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for modeling the hydrodynamic loading and
!jmj   mooring system dynamics for floating wind turbines.  Do this by allowing
!jmj   a keyword in place of the integers 0 or 1 in input PtfmLdMod when
!jmj   PtfmModel = 3:
REAL(ReKi)                   :: PtfmVol0                                        ! Displaced volume of water when the platform is in its undisplaced position.
REAL(ReKi)                   :: RdtnDT                                          ! Time step for wave radiation kernel calculations.
!jmj Start of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Put in some logic to ensure that the hydrodynamic loads are time invariant
!jmj   when linearizing a model:
!remove6.02bREAL(ReKi)                   :: RdtnTMax                                        ! Analysis time for wave radiation kernel calculations.
REAL(ReKi)                   :: RdtnTMax     = 0.0                              ! Analysis time for wave radiation kernel calculations.
!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.

!jmj Start of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Replace the hard-coded mooring line restoring calculation with a general
!jmj   purpose, quasi-static solution based on the analytical catenary cable
!jmj   equations with seabed interaction:
INTEGER(4)                   :: LineMod                                         ! Mooring line model switch.
INTEGER(4)                   :: NumLines     = 0                                ! Number of mooring lines.
!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.

INTEGER(4)                   :: PtfmLdMod    = 0                                ! Platform loading model switch. (Initialized to zero b/c not all models read in PtfmFile)
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for modeling the hydrodynamic loading and
!jmj   mooring system dynamics for floating wind turbines.  Do this by allowing
!jmj   a keyword in place of the integers 0 or 1 in input PtfmLdMod when
!jmj   PtfmModel = 3:
INTEGER(4)                   :: PtfmNodes                                       ! Number of platform nodes used in calculation of viscous drag term from Morison's equation.

!bjj start of proposed change v6.02d-bjj
!rmCHARACTER(99)                :: WAMITFile                                       ! Root name of WAMIT output files containing the linear, nondimensionalized, hydrostatic restoring matrix (.hst extension), frequency-dependent hydrodynamic added mass matrix and damping matrix (.1 extension), and frequency- and direction-dependent wave excitation force vector per unit wave amplitude (.3 extension).
CHARACTER(1024)              :: WAMITFile                                       ! Root name of WAMIT output files containing the linear, nondimensionalized, hydrostatic restoring matrix (.hst extension), frequency-dependent hydrodynamic added mass matrix and damping matrix (.1 extension), and frequency- and direction-dependent wave excitation force vector per unit wave amplitude (.3 extension).
!bjj end of proposed change v6.02d-bjj
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.


!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Rename MODULE PlatformLd() to Platform():
!remove6.02aEND MODULE PlatformLd
END MODULE Platform
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
!=======================================================================
MODULE RotorFurling


   ! This MODULE stores input variables for rotor-furling.


USE                             Precision


REAL(ReKi)                   :: RFrlCDmp  = 0.0                                 ! Rotor-furl rate-independent Coulomb-damping moment. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlDmp   = 0.0                                 ! Rotor-furl damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlDSDmp = 0.0                                 ! Rotor-furl down-stop damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlDSDP  = 0.0                                 ! Rotor-furl down-stop damper position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlDSSP  = 0.0                                 ! Rotor-furl down-stop spring position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlDSSpr = 0.0                                 ! Rotor-furl down-stop spring constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlSpr   = 0.0                                 ! Rotor-furl spring constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlUSDmp = 0.0                                 ! Rotor-furl up-stop damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlUSDP  = 0.0                                 ! Rotor-furl up-stop damper position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlUSSP  = 0.0                                 ! Rotor-furl up-stop spring position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlUSSpr = 0.0                                 ! Rotor-furl up-stop spring constant. (Initialized to zero b/c not all models read in FurlFile)

INTEGER(4)                   :: RFrlMod   = 0                                   ! Rotor-furl spring/damper model switch. (Initialized to zero b/c not all models read in FurlFile)


END MODULE RotorFurling
!=======================================================================
MODULE RtHndSid


   ! This MODULE stores variables used in RtHS.


USE                             Precision


REAL(ReKi)                   :: AngAccEBt(3)                                    ! Portion of the angular acceleration of the base plate                                                (body B) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi)                   :: AngAccERt(3)                                    ! Portion of the angular acceleration of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi)                   :: AngAccEXt(3)                                    ! Portion of the angular acceleration of the platform                                                  (body X) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for modeling the hydrodynamic loading on a
!jmj   monopile.  Do this by reading in addition inputs from the platform file
!jmj   if they exist:
REAL(ReKi), ALLOCATABLE      :: AngAccEFt(:,:)                                  ! Portion of the angular acceleration of tower element J                                               (body F) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.

REAL(ReKi), ALLOCATABLE      :: AngPosHM (:,:,:)                                ! Angular position of eleMent J of blade K                                          (body M) in the hub           (body H          ).
REAL(ReKi)                   :: AngPosXB (3)                                    ! Angular position of the base plate                                                (body B) in the platform      (body X          ).
REAL(ReKi)                   :: AngVelEB (3)                                    ! Angular velocity of the base plate                                                (body B) in the inertia frame (body E for earth).
REAL(ReKi)                   :: AngVelER (3)                                    ! Angular velocity of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame (body E for earth).
REAL(ReKi)                   :: AngVelEX (3)                                    ! Angular velocity of the platform                                                  (body X) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: AugMat   (:,:)                                  ! The augmented matrix used for the solution of the QD2T()s.
REAL(ReKi)                   :: FKAero   (3)                                    ! The tail fin aerodynamic force acting at point K, the center-of-pressure of the tail fin.
REAL(ReKi)                   :: FrcONcRtt(3)                                    ! Portion of the force at yaw bearing         (point O   ) due to the nacelle, generator, and rotor associated with everything but the QD2T()'s.
REAL(ReKi)                   :: FrcPRott (3)                                    ! Portion of the force at the teeter pin      (point P   ) due to the rotor associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE      :: FrcS0Bt  (:,:)                                  ! Portion of the force at the blade root      (point S(0)) due to the blade associated with everything but the QD2T()'s.
REAL(ReKi)                   :: FrcT0Trbt(3)                                    ! Portion of the force at tower base          (point T(0)) due to the turbine associated with everything but the QD2T()'s.
REAL(ReKi)                   :: FrcVGnRtt(3)                                    ! Portion of the force at the rotor-furl axis (point V   ) due to the structure that furls with the rotor, generator, and rotor associated with everything but the QD2T()'s.
REAL(ReKi)                   :: FrcWTailt(3)                                    ! Portion of the force at the  tail-furl axis (point W   ) due to the tail associated with everything but the QD2T()'s.
REAL(ReKi)                   :: FrcZAllt (3)                                    ! Portion of the force at platform reference  (point Z   ) due to everything associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE      :: FSAero   (:,:,:)                                ! The aerodynamic force per unit span acting on a blade at point S.
REAL(ReKi), ALLOCATABLE      :: FSTipDrag(:,:)                                  ! The aerodynamic force at a blade tip resulting from tip drag.
REAL(ReKi), ALLOCATABLE      :: FTAero   (:,:)                                  ! The aerodynamic force per unit length acting on the tower at point T.
REAL(ReKi), ALLOCATABLE      :: FTHydrot (:,:)                                  ! Portion of the hydrodynamic force per unit length acting on the tower at point T associated with everything but the QD2T()'s.
REAL(ReKi)                   :: FZHydrot (3)                                    ! Portion of the platform hydrodynamic force at the platform reference (point Z) associated with everything but the QD2T()'s.
REAL(ReKi)                   :: GBoxEffFac                                      ! The factor used to apply the gearbox efficiency effects to the equation associated with the generator DOF.
REAL(ReKi)                   :: LinAccEIMUt(3)                                  ! Portion of the linear acceleration of the nacelle IMU      (point IMU) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi)                   :: LinAccEOt(3)                                    ! Portion of the linear acceleration of the base plate         (point O) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE      :: LinAccESt(:,:,:)                                ! Portion of the linear acceleration of a point on a blade     (point S) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE      :: LinAccETt(:,:)                                  ! Portion of the linear acceleration of a point on the tower   (point T) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi)                   :: LinAccEZt (3)                                   ! Portion of the linear acceleration of the platform reference (point Z) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.
REAL(ReKi)                   :: LinVelEIMU(3)                                   ! Linear velocity of the nacelle IMU  (point IMU) in the inertia frame.
REAL(ReKi)                   :: LinVelEZ  (3)                                   ! Linear velocity of platform reference (point Z) in the inertia frame.
!bjj start of propsoed change
REAL(ReKi), ALLOCATABLE      :: LinVelESm2 (:)                                  ! The m2-component (closest to tip) of LinVelES.
!bjj end of proposed change
REAL(ReKi)                   :: MAAero   (3)                                    ! The tail fin aerodynamic moment acting at point K, the center-of-pressure of the tail fin.
REAL(ReKi), ALLOCATABLE      :: MFAero   (:,:)                                  ! The aerodynamic moment per unit length acting on the tower at point T.
REAL(ReKi), ALLOCATABLE      :: MFHydrot (:,:)                                  ! Portion of the hydrodynamic moment per unit length acting on a tower element (body F) at point T associated with everything but the QD2T()'s.
REAL(ReKi)                   :: MomBNcRtt(3)                                    ! Portion of the moment at the base plate (body B) / yaw bearing                       (point O   ) due to the nacelle, generator, and rotor associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE      :: MomH0Bt  (:,:)                                  ! Portion of the moment at the hub        (body H) / blade root                        (point S(0)) due to the blade associated with everything but the QD2T()'s.
REAL(ReKi)                   :: MomLPRott(3)                                    ! Portion of the moment at the teeter pin (point P) on the low-speed shaft             (body L    ) due to the rotor associated with everything but the QD2T()'s.
REAL(ReKi)                   :: MomNGnRtt(3)                                    ! Portion of the moment at the nacelle    (body N) / selected point on rotor-furl axis (point V   ) due the structure that furls with the rotor, generator, and rotor associated with everything but the QD2T()'s.
REAL(ReKi)                   :: MomNTailt(3)                                    ! Portion of the moment at the nacelle    (body N) / selected point on  tail-furl axis (point W   ) due the tail associated with everything but the QD2T()'s.
REAL(ReKi)                   :: MomX0Trbt(3)                                    ! Portion of the moment at the platform   (body X) / tower base                        (point T(0)) due to the turbine associated with everything but the QD2T()'s.
REAL(ReKi)                   :: MomXAllt (3)                                    ! Portion of the moment at the platform   (body X) / platform reference                (point Z   ) due to everything associated with everything but the QD2T()'s.
REAL(ReKi), ALLOCATABLE      :: MMAero   (:,:,:)                                ! The aerodynamic moment per unit span acting on a blade at point S.
REAL(ReKi)                   :: MXHydrot (3)                                    ! Portion of the platform hydrodynamic moment acting at the platform (body X) / platform reference (point Z) associated with everything but the  QD2T()'s.
REAL(ReKi), ALLOCATABLE      :: PAngVelEA(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the tail                                                      (body A) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PAngVelEB(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the base plate                                                (body B) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PAngVelEF(:,:,:,:)                              ! Partial angular velocity (and its 1st time derivative) of tower element J                                               (body F) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PAngVelEG(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the generator                                                 (body G) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PAngVelEH(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the hub                                                       (body H) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PAngVelEL(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the low-speed shaft                                           (body L) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PAngVelEM(:,:,:,:,:)                            ! Partial angular velocity (and its 1st time derivative) of eleMent J of blade K                                          (body M) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PAngVelEN(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the nacelle                                                   (body N) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PAngVelER(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the structure that furls with the rotor (not including rotor) (body R) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PAngVelEX(:,:,:)                                ! Partial angular velocity (and its 1st time derivative) of the platform                                                (body B) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PFrcONcRt(:,:)                                  ! Partial force at the yaw bearing     (point O   ) due to the nacelle, generator, and rotor.
REAL(ReKi), ALLOCATABLE      :: PFrcPRot (:,:)                                  ! Partial force at the teeter pin      (point P   ) due to the rotor.
REAL(ReKi), ALLOCATABLE      :: PFrcS0B  (:,:,:)                                ! Partial force at the blade root      (point S(0)) due to the blade.
REAL(ReKi), ALLOCATABLE      :: PFrcT0Trb(:,:)                                  ! Partial force at the tower base      (point T(0)) due to the turbine.
REAL(ReKi), ALLOCATABLE      :: PFrcVGnRt(:,:)                                  ! Partial force at the rotor-furl axis (point V   ) due to the structure that furls with the rotor, generator, and rotor.
REAL(ReKi), ALLOCATABLE      :: PFrcWTail(:,:)                                  ! Partial force at the  tail-furl axis (point W   ) due to the tail.
REAL(ReKi), ALLOCATABLE      :: PFrcZAll (:,:)                                  ! Partial force at the platform reference (point Z) due to everything.
REAL(ReKi), ALLOCATABLE      :: PFTHydro (:,:,:)                                ! Partial hydrodynamic force per unit length acting on the tower at point T.
REAL(ReKi)                   :: PFZHydro (6,3)                                  ! Partial platform hydrodynamic force at the platform reference (point Z).
REAL(ReKi), ALLOCATABLE      :: PLinVelEC(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the hub center of mass            (point C) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PLinVelED(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the center of mass of the structure that furls with the rotor (not including rotor) (point D) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PLinVelEI(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the tail boom center of mass                                                        (point I) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PLinVelEIMU(:,:,:)                              ! Partial linear velocity (and its 1st time derivative) of the nacelle IMU                                                                   (point IMU) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PLinVelEJ(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the tail fin  center of mass                                                        (point J) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PLinVelEK(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the tail fin  center of pressure                                                    (point K) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PLinVelEO(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the base plate                                                                      (point O) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PLinVelEP(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the teeter pin                                                                      (point P) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PLinVelEQ(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the apex of rotation                                                                (point Q) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PLinVelES(:,:,:,:,:)                            ! Partial linear velocity (and its 1st time derivative) of a point on a blade                                                                  (point S) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PLinVelET(:,:,:,:)                              ! Partial linear velocity (and its 1st time derivative) of a point on the tower                                                                (point T) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PLinVelEU(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the nacelle center of mass                                                          (point U) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PLinVelEV(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the selected point on the rotor-furl axis                                           (point V) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PLinVelEW(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the selected point on the  tail-furl axis                                           (point W) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PLinVelEY(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the platform mass center                                                            (point Y) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PLinVelEZ(:,:,:)                                ! Partial linear velocity (and its 1st time derivative) of the platform reference point                                                        (point Z) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: PMFHydro (:,:,:)                                ! Partial hydrodynamic moment per unit length acting on a tower element (body F) at point T.
REAL(ReKi), ALLOCATABLE      :: PMomBNcRt(:,:)                                  ! Partial moment at the base plate (body B) / yaw bearing                       (point O   ) due the nacelle, generator, and rotor.
REAL(ReKi), ALLOCATABLE      :: PMomH0B  (:,:,:)                                ! Partial moment at the hub        (body H) / blade root                        (point S(0)) due to the blade.
REAL(ReKi), ALLOCATABLE      :: PMomLPRot(:,:)                                  ! Partial moment at the teeter pin (point P) on the low-speed shaft             (body L    ) due to the rotor.
REAL(ReKi), ALLOCATABLE      :: PMomNGnRt(:,:)                                  ! Partial moment at the nacelle    (body N) / selected point on rotor-furl axis (point V   ) due the structure that furls with the rotor, generator, and rotor.
REAL(ReKi), ALLOCATABLE      :: PMomNTail(:,:)                                  ! Partial moment at the nacelle    (body N) / selected point on  tail-furl axis (point W   ) due the tail.
REAL(ReKi), ALLOCATABLE      :: PMomX0Trb(:,:)                                  ! Partial moment at the platform   (body X) / tower base                        (point T(0)) due to the turbine.
REAL(ReKi), ALLOCATABLE      :: PMomXAll (:,:)                                  ! Partial moment at the platform   (body X) / platform reference                (point Z   ) due to the everything.
REAL(ReKi)                   :: PMXHydro (6,3)                                  ! Partial platform hydrodynamic moment at the platform (body X) / platform reference (point Z).
REAL(ReKi), ALLOCATABLE      :: QDT      (:)                                    ! Current estimate of QD.
REAL(ReKi), ALLOCATABLE      :: QD2T     (:)                                    ! Solution (acceleration) vector.
REAL(ReKi), ALLOCATABLE      :: QD2TC    (:)                                    ! A copy of the value of QD2T used in SUBROUTINE FixHSSBrTq().
REAL(ReKi), ALLOCATABLE      :: OgnlGeAzRo(:)                                   ! The original elements of AugMat that formed the DOF_GeAz equation before application of known initial conditions.
REAL(ReKi), ALLOCATABLE      :: QT       (:)                                    ! Current estimate of Q for each degree of freedom.
REAL(ReKi)                   :: rO       (3)                                    ! Position vector from inertial frame origin             to tower-top / base plate (point O).
REAL(ReKi), ALLOCATABLE      :: rQS      (:,:,:)                                ! Position vector from the apex of rotation (point Q   ) to a point on a blade (point S).
REAL(ReKi), ALLOCATABLE      :: rS       (:,:,:)                                ! Position vector from inertial frame origin             to a point on a blade (point S).
REAL(ReKi), ALLOCATABLE      :: rS0S     (:,:,:)                                ! Position vector from the blade root       (point S(0)) to a point on a blade (point S).
REAL(ReKi)                   :: rT0O     (3)                                    ! Position vector from the tower base       (point T(0)) to tower-top / base plate (point O).
REAL(ReKi), ALLOCATABLE      :: rT0T     (:,:)                                  ! Position vector from a height of TwrRBHt (base of flexible portion of tower) (point T(0)) to a point on the tower (point T).
REAL(ReKi)                   :: rZ       (3)                                    ! Position vector from inertia frame origin to platform reference (point Z).
REAL(ReKi)                   :: rZO      (3)                                    ! Position vector from platform reference   (point Z   ) to tower-top / base plate (point O).
REAL(ReKi), ALLOCATABLE      :: rZT      (:,:)                                  ! Position vector from platform reference   (point Z   ) to a point on a tower     (point T).
REAL(ReKi), ALLOCATABLE      :: SolnVec  (:)                                    ! Solution vector found by solving the equations of motion.
REAL(ReKi)                   :: TeetAng                                         ! Current teeter angle = QT(DOF_Teet) for 2-blader or 0 for 3-blader (this is used in place of QT(DOF_Teet) throughout RtHS().
REAL(ReKi)                   :: TeetAngVel                                      ! Angular velocity of the teeter motion.


END MODULE RtHndSid
!=======================================================================
MODULE SimCont


   ! This MODULE stores variables for simulation control.


USE                             Precision


REAL(ReKi)                   :: DT                                              ! Integration time step.
REAL(ReKi)                   :: DT24                                            ! DT/24.
REAL(ReKi)                   :: TMax                                            ! Total run time.
REAL(ReKi)                   :: ZTime    = 0.0                                  ! Current simulation time.

REAL(4)                      :: UsrTime1                                        ! User CPU time for simulation initialization.
!!bjj start of proposed change
!REAL(ReKi)                   :: UsrTime0                                        ! Initial time
!!bjj end of proposed change

INTEGER(4)                   :: Step     = 0                                    ! Current simulation time step.


END MODULE SimCont
!=======================================================================
MODULE TailAero


   ! This MODULE stores input variables for tail fin aerodynamics.


USE                             Precision


REAL(ReKi)                   :: SQRTTFinA = 0.0                                 ! = SQRT( TFinArea )
REAL(ReKi)                   :: TFinAOA   = 0.0                                 ! Angle-of-attack between the relative wind velocity and tail fin chordline
REAL(ReKi)                   :: TFinArea  = 0.0                                 ! Tail fin planform area.
REAL(ReKi)                   :: TFinCD    = 0.0                                 ! Tail fin drag            coefficient resulting from current TFinAOA
REAL(ReKi)                   :: TFinCL    = 0.0                                 ! Tail fin lift            coefficient resulting from current TFinAOA
REAL(ReKi)                   :: TFinCM    = 0.0                                 ! Tail fin pitching moment coefficient resulting from current TFinAOA
REAL(ReKi)                   :: TFinKFx   = 0.0                                 ! Aerodynamic force  at the tail fin center-of-pressure (point K) along tail fin chordline pointing toward tail fin trailing edge (N)
REAL(ReKi)                   :: TFinKFy   = 0.0                                 ! Aerodynamic force  at the tail fin center-of-pressure (point K) normal to plane of tail fin pointing towards suction surface    (N)
REAL(ReKi)                   :: TFinKMz   = 0.0                                 ! Aerodynamic moment at the tail fin center-of-pressure (point K) in plane of tail fin normal to chordline and nominally upward   (N-m)
REAL(ReKi)                   :: TFinQ     = 0.0                                 ! Dynamic pressure of the relative wind velocity

INTEGER(4)                   :: TFinMod   = 0                                   ! Tail fin aerodynamics model switch. (Initialized to zero b/c not all models read in FurlFile)
INTEGER(4)                   :: TFinNFoil = 1                                   ! Tail fin airfoil number. (iniated to first airfoil number)

!bjj chg: LOGICAL(1)                   :: SubAxInd  = .FALSE.                             ! Subtract average rotor axial induction when computing relative wind-inflow at tail fin?
LOGICAL                      :: SubAxInd  = .FALSE.                             ! Subtract average rotor axial induction when computing relative wind-inflow at tail fin?


END MODULE TailAero
!=======================================================================
MODULE TailFurling


   ! This MODULE stores input variables for tail-furling.


USE                             Precision


REAL(ReKi)                   :: TFrlCDmp  = 0.0                                 ! Tail-furl rate-independent Coulomb-damping moment. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlDmp   = 0.0                                 ! Tail-furl damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlDSDmp = 0.0                                 ! Tail-furl down-stop damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlDSDP  = 0.0                                 ! Tail-furl down-stop damper position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlDSSP  = 0.0                                 ! Tail-furl down-stop spring position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlDSSpr = 0.0                                 ! Tail-furl down-stop spring constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlSpr   = 0.0                                 ! Tail-furl spring constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlUSDmp = 0.0                                 ! Tail-furl up-stop damping constant. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlUSDP  = 0.0                                 ! Tail-furl up-stop damper position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlUSSP  = 0.0                                 ! Tail-furl up-stop spring position. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlUSSpr = 0.0                                 ! Tail-furl up-stop spring constant. (Initialized to zero b/c not all models read in FurlFile)

INTEGER(4)                   :: TFrlMod   = 0                                   ! Tail-furl spring/damper model switch. (Initialized to zero b/c not all models read in FurlFile)


END MODULE TailFurling
!=======================================================================
MODULE TeeterVars


   ! This MODULE stores input variables for rotor teeter.


USE                             Precision


REAL(ReKi)                   :: TeetCDmp = 0.0                                  ! Rotor-teeter rate-independent Coulomb-damping. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                   :: TeetDmp  = 0.0                                  ! Rotor-teeter damping constant. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                   :: TeetDmpP = 0.0                                  ! Rotor-teeter damper position. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                   :: TeetHSSp = 0.0                                  ! Rotor-teeter hard-stop linear-spring constant. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                   :: TeetHStP = 0.0                                  ! Rotor-teeter hard-stop position. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                   :: TeetSSSp = 0.0                                  ! Rotor-teeter soft-stop linear-spring constant. (initiated to zero b/c the 3-blader requires it to be zero)
REAL(ReKi)                   :: TeetSStP = 0.0                                  ! Rotor-teeter soft-stop position. (initiated to zero b/c the 3-blader requires it to be zero)

INTEGER(4)                   :: TeetMod  = 0                                    ! Rotor-teeter spring/damper model switch. (initiated to zero b/c the 3-blader requires it to be zero)


END MODULE TeeterVars
!=======================================================================
MODULE TipBrakes


   ! This MODULE stores input variables for tip brakes.


USE                             Precision


REAL(ReKi)                   :: TBDrCon                                         ! Instantaneous tip-brake drag constant, Cd*Area.
REAL(ReKi)                   :: TBDrConD                                        ! Tip-brake drag constant during fully-deployed operation, Cd*Area.
REAL(ReKi)                   :: TBDrConN                                        ! Tip-brake drag constant during normal operation, Cd*Area.
REAL(ReKi)                   :: TpBrDT                                          ! Time for tip-brake to reach full deployment once released (sec).


END MODULE TipBrakes
!=======================================================================
MODULE Tower


   ! This MODULE stores variables for the tower.


USE                             Precision


REAL(ReKi)                   :: AdjFASt                                          ! Factor to adjust tower fore-aft stiffness.
REAL(ReKi)                   :: AdjSSSt                                          ! Factor to adjust tower side-to-side stiffness.
REAL(ReKi)                   :: AdjTwMa                                          ! Factor to adjust tower mass density.
REAL(ReKi), ALLOCATABLE      :: AxRedTFA  (:,:,:)                                ! The axial-reduction terms for the fore-aft tower mode shapes.
REAL(ReKi), ALLOCATABLE      :: AxRedTSS  (:,:,:)                                ! The axial-reduction terms for the side-to-side tower mode shapes.
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for modeling the hydrodynamic loading on a
!jmj   monopile.  Do this by reading in addition inputs from the platform file
!jmj   if they exist:
REAL(ReKi), ALLOCATABLE      :: CAT       (:)                                   ! Interpolated, normalized hydrodynamic added mass   coefficient in Morison's equation.
REAL(ReKi), ALLOCATABLE      :: CDT       (:)                                   ! Interpolated, normalized hydrodynamic viscous drag coefficient in Morison's equation.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
REAL(ReKi), ALLOCATABLE      :: cgOffTFA  (:)                                    ! Interpolated tower fore-aft mass cg offset.
REAL(ReKi), ALLOCATABLE      :: cgOffTSS  (:)                                    ! Interpolated tower side-to-side mass cg offset.
REAL(ReKi)                   :: CTFA      (2,2)                                  ! Generalized damping of tower in fore-aft direction.
REAL(ReKi)                   :: CTSS      (2,2)                                  ! Generalized damping of tower in side-to-side direction.
REAL(ReKi), ALLOCATABLE      :: DHNodes   (:)                                    ! Length of variable-length tower elements
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for modeling the hydrodynamic loading on a
!jmj   monopile.  Do this by reading in addition inputs from the platform file
!jmj   if they exist:
REAL(ReKi), ALLOCATABLE      :: DiamT     (:)                                   ! Interpolated tower diameter in Morison's equation.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
REAL(ReKi)                   :: FAStTunr  (2)                                    ! Tower fore-aft modal stiffness tuners.
REAL(ReKi), ALLOCATABLE      :: HNodes    (:)                                    ! Location of variable-spaced tower nodes (relative to the tower rigid base height)
REAL(ReKi), ALLOCATABLE      :: HNodesNorm(:)                                    ! Normalized location of variable-spaced tower nodes (relative to the tower rigid base height) (0 < HNodesNorm(:) < 1)
REAL(ReKi), ALLOCATABLE      :: HtFract   (:)                                    ! Fractional height of the flexible portion of tower for a given input station.
REAL(ReKi), ALLOCATABLE      :: InerTFA   (:)                                    ! Interpolated tower fore-aft (about yt-axis) mass inertia per unit length.
REAL(ReKi), ALLOCATABLE      :: InerTSS   (:)                                    ! Interpolated tower side-to-side (about xt-axis) mass inertia per unit length.
REAL(ReKi)                   :: KTFA      (2,2)                                  ! Generalized stiffness of tower in fore-aft direction.
REAL(ReKi)                   :: KTSS      (2,2)                                  ! Generalized stiffness of tower in side-to-side direction.
REAL(ReKi), ALLOCATABLE      :: MassT     (:)                                    ! Interpolated lineal mass density of tower.
REAL(ReKi)                   :: SSStTunr  (2)                                    ! Tower side-to-side modal stiffness tuners.
REAL(ReKi), ALLOCATABLE      :: StiffTEA  (:)                                    ! Interpolated tower extensional stiffness.
REAL(ReKi), ALLOCATABLE      :: StiffTFA  (:)                                    ! Interpolated fore-aft tower stiffness.
REAL(ReKi), ALLOCATABLE      :: StiffTGJ  (:)                                    ! Interpolated tower torsional stiffness.
REAL(ReKi), ALLOCATABLE      :: StiffTSS  (:)                                    ! Interpolated side-side tower stiffness.
REAL(ReKi), ALLOCATABLE      :: TMassDen  (:)                                    ! Tower mass density for a given input station.
REAL(ReKi), ALLOCATABLE      :: TwEAStif  (:)                                    ! Tower extensional stiffness for a given input station.
REAL(ReKi), ALLOCATABLE      :: TwFAcgOf  (:)                                    ! Tower fore-aft (along the xt-axis) mass cg offset for a given input station.
REAL(ReKi), ALLOCATABLE      :: TwFAIner  (:)                                    ! Tower fore-aft (about yt-axis) mass inertia per unit length for a given input station.
REAL(ReKi), ALLOCATABLE      :: TwFAStif  (:)                                    ! Tower fore-aft stiffness for a given input station.
REAL(ReKi), ALLOCATABLE      :: TwGJStif  (:)                                    ! Tower torsional stiffness for a given input station.
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for modeling the hydrodynamic loading on a
!jmj   monopile.  Do this by reading in addition inputs from the platform file
!jmj   if they exist:
REAL(ReKi)                   :: TwrAM     (6,6) = 0.0                           ! Added mass matrix of the current tower element per unit length.
REAL(ReKi)                   :: TwrCA    = 0.0                                  ! Normalized hydrodynamic added mass   coefficient in Morison's equation.
REAL(ReKi)                   :: TwrCD    = 0.0                                  ! Normalized hydrodynamic viscous drag coefficient in Morison's equation.
REAL(ReKi)                   :: TwrDiam  = 0.0                                  ! Tower diameter in Morison's equation.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
REAL(ReKi)                   :: TwrFADmp  (2)                                    ! Tower fore-aft structural damping ratios.
REAL(ReKi), ALLOCATABLE      :: TwrFASF   (:,:,:)                                ! Tower fore-aft shape functions.
REAL(ReKi)                   :: TwrFlexL                                         ! Height / length of the flexible portion of the tower.
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for modeling the hydrodynamic loading on a
!jmj   monopile.  Do this by reading in addition inputs from the platform file
!jmj   if they exist:
REAL(ReKi)                   :: TwrFt     (6)   = 0.0                            ! The surge/xi (1), sway/yi (2), and heave/zi (3)-components of the portion of the tower force at the current tower element (point T) and the roll/xi (4), pitch/yi (5), and yaw/zi (6)-components of the portion of the tower moment acting at the current tower element (body F) / (point T) per unit length associated with everything but the QD2T()'s.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.

REAL(ReKi)                   :: TwrSSDmp  (2)                                    ! Tower side-to-side structural damping ratios.
REAL(ReKi), ALLOCATABLE      :: TwrSSSF   (:,:,:)                                ! Tower side-to-side shape functions.
REAL(ReKi), ALLOCATABLE      :: TwSScgOf  (:)                                    ! Tower fore-aft (along the yt-axis) mass cg offset for a given input station.
REAL(ReKi), ALLOCATABLE      :: TwSSIner  (:)                                    ! Tower side-to-side (about xt-axis) mass inertia per unit length for a given input station.
REAL(ReKi), ALLOCATABLE      :: TwSSStif  (:)                                    ! Tower side-to-side stiffness for a given input station.

INTEGER(4)                   :: NTwInpSt                                        ! Number of tower input stations.
INTEGER(4)                   :: TTopNode                                        ! Index of the additional node located at the tower-top = TwrNodes + 1
!jmj Start of proposed change.  v6.02a-jmj  25-Aug-2006.
!jmj Add an undocumented feature for modeling the hydrodynamic loading on a
!jmj   monopile.  Do this by reading in addition inputs from the platform file
!jmj   if they exist:
INTEGER(4)                   :: TwrLdMod = 0                                    ! Tower loading model switch.
!jmj End of proposed change.  v6.02a-jmj  25-Aug-2006.
INTEGER(4)                   :: TwrNodes                                        ! Number of tower nodes used in the analysis.


END MODULE Tower
!=======================================================================
MODULE TurbConf


   ! This MODULE stores variables for turbine configuration.


USE                             Precision


REAL(ReKi)                   :: AvgNrmTpRd                                      ! Average tip radius normal to the saft.
REAL(ReKi)                   :: AzimB1Up                                        ! Azimuth value to use for I/O when blade 1 points up.
REAL(ReKi)                   :: BoomCMxn  = 0.0                                 ! Downwind distance from tower-top to tail boom CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: BoomCMyn  = 0.0                                 ! Lateral  distance from tower-top to tail boom CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: BoomCMzn  = 0.0                                 ! Vertical distance from tower-top to tail boom CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: CosDel3   = 1.0                                 ! Cosine of the Delta-3 angle for teetering rotors.
REAL(ReKi), ALLOCATABLE      :: CosPreC   (:)                                   ! Cosines of the precone angles.
REAL(ReKi)                   :: CRFrlSkew = 1.0                                 ! Cosine of the rotor-furl axis skew angle.
REAL(ReKi)                   :: CRFrlSkw2 = 1.0                                 ! Cosine-squared of the rotor-furl axis skew angle.
REAL(ReKi)                   :: CRFrlTilt = 1.0                                 ! Cosine of the rotor-furl axis tilt angle.
REAL(ReKi)                   :: CRFrlTlt2 = 1.0                                 ! Cosine-squared of the rotor-furl axis tilt angle.
REAL(ReKi)                   :: CShftSkew = 1.0                                 ! Cosine of the shaft skew angle.
REAL(ReKi)                   :: CShftTilt                                       ! Cosine of the shaft tilt angle.
REAL(ReKi)                   :: CSRFrlSkw = 0.0                                 ! Cosine*Sine of the rotor-furl axis skew angle.
REAL(ReKi)                   :: CSRFrlTlt = 0.0                                 ! Cosine*Sine of the rotor-furl axis tilt angle.
REAL(ReKi)                   :: CSTFrlSkw = 0.0                                 ! Cosine*Sine of the tail-furl axis skew angle.
REAL(ReKi)                   :: CSTFrlTlt = 0.0                                 ! Cosine*Sine of the tail-furl axis tilt angle.
REAL(ReKi)                   :: CTFinBank = 1.0                                 ! Cosine of the tail fin planform  bank angle.
REAL(ReKi)                   :: CTFinSkew = 1.0                                 ! Cosine of the tail fin chordline skew angle.
REAL(ReKi)                   :: CTFinTilt = 1.0                                 ! Cosine of the tail fin chordline tilt angle.
REAL(ReKi)                   :: CTFrlSkew = 1.0                                 ! Cosine of the tail-furl axis skew angle.
REAL(ReKi)                   :: CTFrlSkw2 = 1.0                                 ! Cosine-squared of the tail-furl axis skew angle.
REAL(ReKi)                   :: CTFrlTilt = 1.0                                 ! Cosine of the tail-furl axis tilt angle.
REAL(ReKi)                   :: CTFrlTlt2 = 1.0                                 ! Cosine-squared of the tail-furl axis tilt angle.
REAL(ReKi)                   :: Delta3    = 0.0                                 ! Delta-3 angle for teetering rotors.
REAL(ReKi)                   :: FASTHH                                          ! Hub-height as computed using FAST inputs [= TowerHt + Twr2Shft + OverHang*SIN( ShftTilt ) ].
REAL(ReKi)                   :: HubCM                                           ! Distance from rotor apex to hub mass.
REAL(ReKi)                   :: HubRad                                          ! Preconed hub radius.
REAL(ReKi)                   :: NacCMxn                                         ! Downwind distance from tower-top to nacelle CM.
REAL(ReKi)                   :: NacCMyn                                         ! Lateral  distance from tower-top to nacelle CM.
REAL(ReKi)                   :: NacCMzn                                         ! Vertical distance from tower-top to nacelle CM.
REAL(ReKi)                   :: OverHang                                        ! Distance from yaw axis to rotor apex or teeter pin.
REAL(ReKi), ALLOCATABLE      :: PreCone   (:)                                   ! Rotor precone angle.
REAL(ReKi)                   :: ProjArea                                        ! Swept area of the rotor projected onto the rotor plane (the plane normal to the low-speed shaft).
REAL(ReKi)                   :: PtfmCM    = 0.0                                 ! Downward distance from the ground [onshore] or MSL [offshore] to the platform CM. (Initialized to zero b/c not all models read in PtfmFile)
REAL(ReKi)                   :: PtfmRef   = 0.0                                 ! Downward distance from the ground [onshore] or MSL [offshore] to the platform reference point. (Initialized to zero b/c not all models read in PtfmFile)
!bjj start of proposed change
!rmREAL(ReKi)                   :: RefHH                                           ! Vertical distance between AeroDyn's wind reference (hub) height (variable HH     ) and FAST's inertia frame reference point (variable PtfmRef); that is, RefHH    = HH      + PtfmRef.
!bjj end of proposed change
REAL(ReKi)                   :: RefTwrHt                                        ! Vertical distance between FAST's undisplaced tower       height (variable TowerHt) and FAST's inertia frame reference point (variable PtfmRef); that is, RefTwrHt = TowerHt + PtfmRef.
REAL(ReKi)                   :: RFrlCMxn  = 0.0                                 ! Downwind distance from tower-top to rotor-furl CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlCMyn  = 0.0                                 ! Lateral  distance from tower-top to rotor-furl CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlCMzn  = 0.0                                 ! Vertical distance from tower-top to rotor-furl CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlPntxn = 0.0                                 ! Downwind distance from tower-top to arbitrary point on rotor-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlPntyn = 0.0                                 ! Lateral  distance from tower-top to arbitrary point on rotor-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlPntzn = 0.0                                 ! Vertical distance from tower-top to arbitrary point on rotor-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlSkew  = 0.0                                 ! Rotor-furl axis skew angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: RFrlTilt  = 0.0                                 ! Rotor-furl axis tilt angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: rVDxn     = 0.0                                 ! xn-component of position vector rVD.
REAL(ReKi)                   :: rVDyn     = 0.0                                 ! yn-component of position vector rVD.
REAL(ReKi)                   :: rVDzn     = 0.0                                 ! zn-component of position vector rVD.
REAL(ReKi)                   :: rVIMUxn                                         ! xn-component of position vector rVIMU.
REAL(ReKi)                   :: rVIMUyn                                         ! yn-component of position vector rVIMU.
REAL(ReKi)                   :: rVIMUzn                                         ! zn-component of position vector rVIMU.
REAL(ReKi)                   :: rVPxn     = 0.0                                 ! xn-component of position vector rVP.
REAL(ReKi)                   :: rVPyn     = 0.0                                 ! yn-component of position vector rVP.
REAL(ReKi)                   :: rVPzn                                           ! zn-component of position vector rVP (need not be initialized to zero).
REAL(ReKi)                   :: rWIxn     = 0.0                                 ! xn-component of position vector rWI.
REAL(ReKi)                   :: rWIyn     = 0.0                                 ! yn-component of position vector rWI.
REAL(ReKi)                   :: rWIzn     = 0.0                                 ! zn-component of position vector rWI.
REAL(ReKi)                   :: rWJxn     = 0.0                                 ! xn-component of position vector rWJ.
REAL(ReKi)                   :: rWJyn     = 0.0                                 ! yn-component of position vector rWJ.
REAL(ReKi)                   :: rWJzn     = 0.0                                 ! zn-component of position vector rWJ.
REAL(ReKi)                   :: rWKxn     = 0.0                                 ! xn-component of position vector rWK.
REAL(ReKi)                   :: rWKyn     = 0.0                                 ! yn-component of position vector rWK.
REAL(ReKi)                   :: rWKzn     = 0.0                                 ! zn-component of position vector rWK.
REAL(ReKi)                   :: rZT0zt                                          ! zt-component of position vector rZT0.
REAL(ReKi)                   :: rZYzt     = 0.0                                 ! zt-component of position vector rZY.
REAL(ReKi)                   :: ShftSkew  = 0.0                                 ! Rotor shaft skew angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: ShftTilt                                        ! Rotor shaft tilt angle.
REAL(ReKi)                   :: SinDel3   = 0.0                                 ! Sine of the Delta-3 angle for teetering rotors.
REAL(ReKi), ALLOCATABLE      :: SinPreC   (:)                                   ! Sines of the precone angles.
REAL(ReKi)                   :: SRFrlSkew = 0.0                                 ! Sine of the rotor-furl axis skew angle.
REAL(ReKi)                   :: SRFrlSkw2 = 0.0                                 ! Sine-squared of the rotor-furl axis skew angle.
REAL(ReKi)                   :: SRFrlTilt = 0.0                                 ! Sine of the rotor-furl axis tilt angle.
REAL(ReKi)                   :: SRFrlTlt2 = 0.0                                 ! Sine-squared of the rotor-furl axis tilt angle.
REAL(ReKi)                   :: SShftSkew = 0.0                                 ! Sine of the shaft skew angle.
REAL(ReKi)                   :: SShftTilt                                       ! Sine of the shaft tilt angle.
REAL(ReKi)                   :: STFinBank = 0.0                                 ! Sine of the tail fin planform  bank angle.
REAL(ReKi)                   :: STFinSkew = 0.0                                 ! Sine of the tail fin chordline skew angle.
REAL(ReKi)                   :: STFinTilt = 0.0                                 ! Sine of the tail fin chordline tilt angle.
REAL(ReKi)                   :: STFrlSkew = 0.0                                 ! Sine of the tail-furl axis skew angle.
REAL(ReKi)                   :: STFrlSkw2 = 0.0                                 ! Sine-squared of the tail-furl axis skew angle.
REAL(ReKi)                   :: STFrlTilt = 0.0                                 ! Sine of the tail-furl axis tilt angle.
REAL(ReKi)                   :: STFrlTlt2 = 0.0                                 ! Sine-squared of the tail-furl axis tilt angle.
REAL(ReKi)                   :: TFrlPntxn = 0.0                                 ! Downwind distance from tower-top to arbitrary point on tail-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlPntyn = 0.0                                 ! Lateral  distance from tower-top to arbitrary point on tail-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlPntzn = 0.0                                 ! Vertical distance from tower-top to arbitrary point on tail-furl axis. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlSkew  = 0.0                                 ! Rotor-furl axis skew angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFrlTilt  = 0.0                                 ! Rotor-furl axis tilt angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFinBank  = 0.0                                 ! Tail fin planform  bank angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFinCMxn  = 0.0                                 ! Downwind distance from tower-top to tail fin CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFinCMyn  = 0.0                                 ! Lateral  distance from tower-top to tail fin CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFinCMzn  = 0.0                                 ! Vertical distance from tower-top to tail fin CM. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFinCPxn  = 0.0                                 ! Downwind distance from tower-top to tail fin CP. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFinCPyn  = 0.0                                 ! Lateral  distance from tower-top to tail fin CP. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFinCPzn  = 0.0                                 ! Vertical distance from tower-top to tail fin CP. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFinSkew  = 0.0                                 ! Tail fin chordline skew angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TFinTilt  = 0.0                                 ! Tail fin chordline tilt angle. (Initialized to zero b/c not all models read in FurlFile)
REAL(ReKi)                   :: TipRad                                          ! Preconed blade-tip radius.
REAL(ReKi)                   :: TowerHt                                         ! Height of tower above ground level.
REAL(ReKi)                   :: Twr2Shft                                        ! Vertical distance from the tower-top to the rotor shaft.
REAL(ReKi)                   :: TwrDraft  = 0.0                                 ! Downward distance from the ground [onshore] or MSL [offshore] to the tower base platform connection. (Initialized to zero b/c not all models read in PtfmFile)
REAL(ReKi)                   :: TwrRBHt                                         ! Tower rigid base height.
REAL(ReKi)                   :: UndSling                                        ! Undersling length.
REAL(ReKi)                   :: Yaw2Shft  = 0.0                                 ! Lateral distance from the yaw axis to the rotor shaft. (Initialized to zero b/c not all models read in FurlFile)

INTEGER(4)                   :: NumBl                                           ! Number of blades.
INTEGER(4)                   :: PSpnElN   = 1                                   ! Number of the innermost blade element which is still part of the pitchable portion of the blade for partial-span pitch control.


END MODULE TurbConf
!=======================================================================
MODULE TurbCont


   ! This MODULE stores input variables for turbine control.


USE                             Precision


REAL(ReKi), ALLOCATABLE      :: BlPitch  (:)                                    ! Initial and current blade pitch angles.
REAL(ReKi), ALLOCATABLE      :: BlPitchCom(:)                                   ! Commanded blade pitch angles.
REAL(ReKi), ALLOCATABLE      :: BlPitchF (:)                                    ! Final blade pitch.
REAL(ReKi), ALLOCATABLE      :: BlPitchFrct(:)                                  ! Blade pitch angle fractions used for the override pitch maneuver calculation.
REAL(ReKi), ALLOCATABLE      :: BlPitchI (:)                                    ! Initial blade pitch angles at the start of the override pitch maneuver.
REAL(ReKi)                   :: NacYawF                                         ! Final yaw angle after override yaw maneuver.
REAL(ReKi)                   :: SpdGenOn                                        ! Generator speed to turn on the generator for a startup.
REAL(ReKi), ALLOCATABLE      :: TBDepISp (:)                                    ! Deployment-initiation speed for the tip brakes.
REAL(ReKi)                   :: THSSBrDp                                        ! Time to initiate deployment of the shaft brake.
REAL(ReKi)                   :: THSSBrFl                                        ! Time at which shaft brake is fully deployed.
REAL(ReKi)                   :: TiDynBrk                                        ! Time to initiate deployment of the dynamic generator brake.
REAL(ReKi)                   :: TimGenOf                                        ! Time to turn off generator for braking or modeling a run-away.
REAL(ReKi)                   :: TimGenOn                                        ! Time to turn on generator for startup.
REAL(ReKi)                   :: TPCOn                                           ! Time to enable active pitch control.
REAL(ReKi), ALLOCATABLE      :: TPitManE (:)                                    ! Time to end pitch maneuvers for each blade.
REAL(ReKi), ALLOCATABLE      :: TPitManS (:)                                    ! Time to start pitch maneuvers for each blade.
REAL(ReKi), ALLOCATABLE      :: TTpBrDp  (:)                                    ! Times to initiate deployment of tip brakes.
REAL(ReKi), ALLOCATABLE      :: TTpBrFl  (:)                                    ! Times at which tip brakes are fully deployed.
REAL(ReKi)                   :: TYawManE                                        ! Time to end override yaw maneuver.
REAL(ReKi)                   :: TYawManS                                        ! Time to start override yaw maneuver.
REAL(ReKi)                   :: TYCOn                                           ! Time to enable active yaw control.
REAL(ReKi)                   :: VS_Rgn2K                                        ! Generator torque constant in Region 2 (HSS side), N-m/rpm^2.
REAL(ReKi)                   :: VS_RtGnSp                                       ! Rated generator speed (HSS side), rpm.
REAL(ReKi)                   :: VS_RtTq                                         ! Rated generator torque/constant generator torque in Region 3 (HSS side), N-m.
REAL(ReKi)                   :: VS_Slope                                        ! Torque/speed slope of region 2 1/2 induction generator.
REAL(ReKi)                   :: VS_SlPc                                         ! Rated generator slip percentage in Region 2 1/2, %.
REAL(ReKi)                   :: VS_SySp                                         ! Synchronous speed of region 2 1/2 induction generator.
REAL(ReKi)                   :: VS_TrGnSp                                       ! Transitional generator speed between regions 2 and 2 1/2.
REAL(ReKi)                   :: YawPosCom                                       ! Commanded yaw angle from user-defined routines, rad.
REAL(ReKi)                   :: YawRateCom                                      ! Commanded yaw rate  from user-defined routines, rad/s.

INTEGER(4)                   :: GenModel                                        ! Generator model
INTEGER(4)                   :: HSSBrMode                                       ! HSS brake model.
INTEGER(4)                   :: PCMode                                          ! Pitch control mode
INTEGER(4)                   :: VSContrl                                        ! Variable-speed-generator control switch.
INTEGER(4)                   :: YCMode                                          ! Yaw control mode

!bjj chg: LOGICAL(1), ALLOCATABLE      :: BegPitMan(:)                                    ! .TRUE. before the override pitch manuever has begun (begin pitch manuever).
LOGICAL,    ALLOCATABLE      :: BegPitMan(:)                                    ! .TRUE. before the override pitch manuever has begun (begin pitch manuever).
!bjj chg: LOGICAL(1)                   :: GenTiStp                                        ! Stop generator based upon T: time or F: generator power = 0.
LOGICAL                      :: GenTiStp                                        ! Stop generator based upon T: time or F: generator power = 0.
!bjj chg: LOGICAL(1)                   :: GenTiStr                                        ! Start generator based upon T: time or F: generator speed.
LOGICAL                      :: GenTiStr                                        ! Start generator based upon T: time or F: generator speed.


END MODULE TurbCont
!=======================================================================
