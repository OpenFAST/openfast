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
REAL(ReKi), PARAMETER        :: FrSrfcSpc =  5.0                                ! Distance between points on the still water level plane along the incident wave propogation heading direction for depicting the free surface where the elevation of the incident waves will be computed used for free surface GRAPHICS. (meters)  !JASON: MAKE THIS AN ACTUAL INPUT TO THE PROGRAM IN ADAMSFile WHEN YOU DOCUMENT THESE ROUTINES!!!!!
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

INTEGER(4)                   :: NFreeSrfc = -1                                  ! Number of points on free surface (not including the zero'th point) where the elevation of the incident waves will be computed (computed every FrSrfcSpc meters along the incident wave propogation heading direction for a length of the rotor diameter).
INTEGER(4)                   :: NLnNodes  = 10                                  ! Number of nodes per line for mooring line GRAPHICS.  !JASON: MAKE THIS AN ACTUAL INPUT TO THE PROGRAM IN ADAMSFile WHEN YOU DOCUMENT THESE ROUTINES!!!!!
INTEGER(4)                   :: NSides                                          ! The number of sides used in GRAPHICS CYLINDER and FRUSTUM statements.

LOGICAL                      :: MakeLINacf                                      ! Switch for making an ADAMS/LINEAR control command file.  To prevent an ADAMS/LINEAR control command file to be made, and to not include the RESULTS statement in the ADAMS dataset, set to .FALSE.
LOGICAL                      :: SaveGrphcs                                      ! Switch to determine whether or note GRAPHICS output is saved in an ADAMS analysis.



END MODULE ADAMSInput
!=======================================================================
MODULE AeroElem


   ! This MODULE stores FAST/AeroDyn interface variables.


USE                             Precision
USE                             AeroDyn  ! for type;  Precision is also included so the previous line could be removed, too.


TYPE(AllAeroMarkers)          :: ADAeroMarkers
TYPE(AeroLoadsOptions)        :: ADIntrfaceOptions
TYPE(AllAeroLoads)            :: ADAeroLoads
TYPE(AeroConfig)              :: ADInterfaceComponents                        ! The configuration markers that make up the bodies where aerodynamic calculations will be needed

INTEGER                       :: NumADBldNodes = 0                               ! Number of blade nodes in AeroDyn


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


REAL(ReKi), PARAMETER        :: Inv2Pi   =  0.15915494                          ! 0.5/Pi.
REAL(ReKi)                   :: TwoPiNB                                         ! 2*Pi/NumBl.  This constant is calculated in fast_io.f90/Inputs()

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
REAL(ReKi)                   :: HSSBrFrac                                       ! Fraction of full braking torque: 0 (off) <= HSSBrFrac <= 1 (full), (-). !bjj: used to be local variable in FAST.f90/Subroutine DrvTrTrq()
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

LOGICAL                      :: GBRevers                                        ! Gearbox reversal flag.


END MODULE DriveTrain
!=======================================================================
MODULE EnvCond


   ! This MODULE stores input variables for environmental conditions.


USE                             Precision


REAL(ReKi)                   :: AirDens                                         ! Air density = RHO.
REAL(ReKi)                   :: CurrDIDir = 0.0                                 ! Depth-independent current heading direction.
REAL(ReKi)                   :: CurrDIV   = 0.0                                 ! Depth-independent current velocity.
REAL(ReKi)                   :: CurrNSDir = 0.0                                 ! Near-surface current heading direction.
REAL(ReKi)                   :: CurrNSRef = 0.0                                 ! Near-surface current reference depth.
REAL(ReKi)                   :: CurrNSV0  = 0.0                                 ! Near-surface current velocity at still water level.
REAL(ReKi)                   :: CurrSSDir = 0.0                                 ! Sub-surface current heading direction.
REAL(ReKi)                   :: CurrSSV0  = 0.0                                 ! Sub-surface current velocity at still water level.
REAL(ReKi)                   :: Gravity                                         ! Gravitational acceleration.

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

CHARACTER(1024)              :: GHWvFile  = ''                                  ! The root name of GH Bladed files containing wave data.


END MODULE EnvCond
!=======================================================================
MODULE Features


   ! This MODULE stores input variables for feature switches.

LOGICAL                   :: CompAero                                        ! Compute aerodynamic forces switch.
LOGICAL                   :: CompHydro = .FALSE.                             ! Compute hydrodynamic forces switch.

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
INTEGER(4)                   :: PtfmModel                                       ! Platform model {0: none, 1: onshore, 2: fixed bottom offshore, 3: floating offshore} (switch).
INTEGER(4)                   :: StrtTime (8)                                    ! Start time of simulation.
INTEGER(4)                   :: UnAC      = 24                                  ! I/O unit number for the ADAMS control output file (.acf) useful for an ADAMS SIMULATE analysis.
INTEGER(4)                   :: UnAD      = 23                                  ! I/O unit number for the ADAMS dataset output file (.adm).
INTEGER(4)                   :: UnAL      = 25                                  ! I/O unit number for the ADAMS control output file (.acf) useful for an ADAMS LINEAR analysis.
INTEGER(4)                   :: UnIn      = 20                                  ! I/O unit number for the input files.
INTEGER(4)                   :: UnLn      = 26                                  ! I/O unit number for the FAST linear output file (.lin).
INTEGER(4)                   :: UnNoSpec  = 27                                  ! I/O unit number for the noise spectr output file.
INTEGER(4)                   :: UnNoSPL   = 28                                  ! I/O unit number for the noise SPL output file.
INTEGER(4)                   :: UnOu      = 21                                  ! I/O unit number for the tabular output file.
INTEGER(4)                   :: UnOuBin   = 29                                  ! I/O unit number for the binary output file.
INTEGER(4)                   :: UnSu      = 22                                  ! I/O unit number for the summary output file.

LOGICAL                      :: Cmpl4SFun  = .FALSE.                            ! Is FAST being compiled as an S-Function for Simulink?
LOGICAL                      :: Cmpl4LV    = .FALSE.                            ! Is FAST being compiled for Labview?
LOGICAL                      :: Furling                                         ! Read in additional model properties for furling turbine?
LOGICAL                      :: SumDisp                                         ! Display summary data on screen?
LOGICAL                      :: SumPrint                                        ! Print summary data to "*.fsm"?

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

LOGICAL                      :: CalcStdy                                        ! Calculate periodic steady state condition (False: linearize about zero) (switch).
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


INTEGER(4), PARAMETER        :: PolyOrd = 6                                     ! Order of the polynomial describing the mode shape.

REAL(ReKi), ALLOCATABLE      :: BldEdgSh (:,:)                                  ! Blade-edge-mode shape coefficients.
REAL(ReKi), ALLOCATABLE      :: BldFl1Sh (:,:)                                  ! Blade-flap-mode-1 shape coefficients.
REAL(ReKi), ALLOCATABLE      :: BldFl2Sh (:,:)                                  ! Blade-flap-mode-2 shape coefficients.
REAL(ReKi), ALLOCATABLE      :: FreqBE   (:,:,:)                                ! Blade edgewise natural frequencies (both w/ and w/o centrifugal stiffening)
REAL(ReKi), ALLOCATABLE      :: FreqBF   (:,:,:)                                ! Blade flapwise natural frequencies (both w/ and w/o centrifugal stiffening)
REAL(ReKi)                   :: FreqTFA  (2,2)                                  ! Computed fore-aft tower natural frequencies.
REAL(ReKi)                   :: FreqTSS  (2,2)                                  ! Computed side-to-side tower natural frequencies.
REAL(ReKi)                   :: TwFAM1Sh (2:PolyOrd)                            ! Tower fore-aft mode-1 shape coefficients.
REAL(ReKi)                   :: TwFAM2Sh (2:PolyOrd)                            ! Tower fore-aft mode-2 shape coefficients.    
REAL(ReKi)                   :: TwSSM1Sh (2:PolyOrd)                            ! Tower side-to-side mode-1 shape coefficients.
REAL(ReKi)                   :: TwSSM2Sh (2:PolyOrd)                            ! Tower side-to-side mode-2 shape coefficients.

LOGICAL                      :: CalcBMode                                       ! T: calculate blade mode shapes internally, F: use blade mode shapes from the blade file.
LOGICAL,    ALLOCATABLE      :: CalcBModes(:)                                   ! Holds CalcBMode for all of the blades.
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
      REAL(ReKi), INTENT(IN )    :: ModShpAry(2:PolyOrd)                        ! Array holding mode shape coefficients
      REAL(ReKi)                 :: SHP            ! The shape function returned by this function.

      INTEGER(4), INTENT(IN )    :: Deriv          ! Which derivative to compute Deriv = 0 (regular function SHP), 1 (D(SHP)/DZ), 2 (D2(SHP)/DZ2)


   ! Lccal variables:

      INTEGER(4)                 :: CoefTmp        !Temporary coefficient
      INTEGER(4)                 :: I              !Counts through polynomial array.
      INTEGER(4)                 :: Swtch(0:2)     !Corresponds to which derivative to compute.  Sets all portions of the coefficient = 0 except those that are relevant.


      Swtch        = 0 !Initialize Swtch(:) to 0
      Swtch(Deriv) = 1
      SHP          = 0.0

      DO I = 2,PolyOrd
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


USE                             NWTC_Library

! ==================================================================================================="
! NOTE: The following lines of code were generated by a Matlab script called "Write_ChckOutLst.m"
!      using the parameters listed in the "OutListParameters.xlsx" Excel file. Any changes to these 
!      lines should be modified in the Matlab script and/or Excel worksheet as necessary. 
! ==================================================================================================="
! This code was generated by Write_ChckOutLst.m at 30-Jan-2012 12:30:08.
   TYPE OutParmType                                                               ! User-defined type for output parameters
      INTEGER                    :: Indx                                          ! Index into AllOuts output array
      CHARACTER(10)              :: Name                                          ! Name of the output parameter.
      CHARACTER(10)              :: Units                                         ! Units corresponding to the output parameter.
      INTEGER                    :: SignM                                         ! sign (output multiplier).
   END TYPE OutParmType


     ! Parameters:


     ! Indices for computing output channels:
     ! NOTES: 
     !    (1) These parameters are in the order stored in "OutListParameters.xlsx"
     !    (2) Array AllOuts() must be dimensioned to the value of the largest output parameter
     !    (3) If an index (MaxOutPts) ever becomes greater or equal to 1000, the logic to create ARRAY/1 in the FAST-to-ADAMS preprocessor will have to be changed.

     !  Time: 

   INTEGER, PARAMETER             :: Time      =    0


  ! Wind Motions:

   INTEGER, PARAMETER             :: WindVxi   =    1
   INTEGER, PARAMETER             :: WindVyi   =    2
   INTEGER, PARAMETER             :: WindVzi   =    3
   INTEGER, PARAMETER             :: TotWindV  =    4
   INTEGER, PARAMETER             :: HorWindV  =    5
   INTEGER, PARAMETER             :: HorWndDir =    6
   INTEGER, PARAMETER             :: VerWndDir =    7


  ! Blade 1 Tip Motions:

   INTEGER, PARAMETER             :: TipDxc1   =    8
   INTEGER, PARAMETER             :: TipDyc1   =    9
   INTEGER, PARAMETER             :: TipDzc1   =   10
   INTEGER, PARAMETER             :: TipDxb1   =   11
   INTEGER, PARAMETER             :: TipDyb1   =   12
   INTEGER, PARAMETER             :: TipALxb1  =   13
   INTEGER, PARAMETER             :: TipALyb1  =   14
   INTEGER, PARAMETER             :: TipALzb1  =   15
   INTEGER, PARAMETER             :: TipRDxb1  =   16
   INTEGER, PARAMETER             :: TipRDyb1  =   17
   INTEGER, PARAMETER             :: TipRDzc1  =   18
   INTEGER, PARAMETER             :: TipClrnc1 =   19


  ! Blade 2 Tip Motions:

   INTEGER, PARAMETER             :: TipDxc2   =   20
   INTEGER, PARAMETER             :: TipDyc2   =   21
   INTEGER, PARAMETER             :: TipDzc2   =   22
   INTEGER, PARAMETER             :: TipDxb2   =   23
   INTEGER, PARAMETER             :: TipDyb2   =   24
   INTEGER, PARAMETER             :: TipALxb2  =   25
   INTEGER, PARAMETER             :: TipALyb2  =   26
   INTEGER, PARAMETER             :: TipALzb2  =   27
   INTEGER, PARAMETER             :: TipRDxb2  =   28
   INTEGER, PARAMETER             :: TipRDyb2  =   29
   INTEGER, PARAMETER             :: TipRDzc2  =   30
   INTEGER, PARAMETER             :: TipClrnc2 =   31


  ! Blade 3 Tip Motions:

   INTEGER, PARAMETER             :: TipDxc3   =   32
   INTEGER, PARAMETER             :: TipDyc3   =   33
   INTEGER, PARAMETER             :: TipDzc3   =   34
   INTEGER, PARAMETER             :: TipDxb3   =   35
   INTEGER, PARAMETER             :: TipDyb3   =   36
   INTEGER, PARAMETER             :: TipALxb3  =   37
   INTEGER, PARAMETER             :: TipALyb3  =   38
   INTEGER, PARAMETER             :: TipALzb3  =   39
   INTEGER, PARAMETER             :: TipRDxb3  =   40
   INTEGER, PARAMETER             :: TipRDyb3  =   41
   INTEGER, PARAMETER             :: TipRDzc3  =   42
   INTEGER, PARAMETER             :: TipClrnc3 =   43


  ! Blade 1 Local Span Motions:

   INTEGER, PARAMETER             :: Spn1ALxb1 =   44
   INTEGER, PARAMETER             :: Spn1ALyb1 =   45
   INTEGER, PARAMETER             :: Spn1ALzb1 =   46
   INTEGER, PARAMETER             :: Spn2ALxb1 =   47
   INTEGER, PARAMETER             :: Spn2ALyb1 =   48
   INTEGER, PARAMETER             :: Spn2ALzb1 =   49
   INTEGER, PARAMETER             :: Spn3ALxb1 =   50
   INTEGER, PARAMETER             :: Spn3ALyb1 =   51
   INTEGER, PARAMETER             :: Spn3ALzb1 =   52
   INTEGER, PARAMETER             :: Spn4ALxb1 =   53
   INTEGER, PARAMETER             :: Spn4ALyb1 =   54
   INTEGER, PARAMETER             :: Spn4ALzb1 =   55
   INTEGER, PARAMETER             :: Spn5ALxb1 =   56
   INTEGER, PARAMETER             :: Spn5ALyb1 =   57
   INTEGER, PARAMETER             :: Spn5ALzb1 =   58
   INTEGER, PARAMETER             :: Spn6ALxb1 =   59
   INTEGER, PARAMETER             :: Spn6ALyb1 =   60
   INTEGER, PARAMETER             :: Spn6ALzb1 =   61
   INTEGER, PARAMETER             :: Spn7ALxb1 =   62
   INTEGER, PARAMETER             :: Spn7ALyb1 =   63
   INTEGER, PARAMETER             :: Spn7ALzb1 =   64
   INTEGER, PARAMETER             :: Spn8ALxb1 =   65
   INTEGER, PARAMETER             :: Spn8ALyb1 =   66
   INTEGER, PARAMETER             :: Spn8ALzb1 =   67
   INTEGER, PARAMETER             :: Spn9ALxb1 =   68
   INTEGER, PARAMETER             :: Spn9ALyb1 =   69
   INTEGER, PARAMETER             :: Spn9ALzb1 =   70
   INTEGER, PARAMETER             :: Spn1TDxb1 =   71
   INTEGER, PARAMETER             :: Spn1TDyb1 =   72
   INTEGER, PARAMETER             :: Spn1TDzb1 =   73
   INTEGER, PARAMETER             :: Spn2TDxb1 =   74
   INTEGER, PARAMETER             :: Spn2TDyb1 =   75
   INTEGER, PARAMETER             :: Spn2TDzb1 =   76
   INTEGER, PARAMETER             :: Spn3TDxb1 =   77
   INTEGER, PARAMETER             :: Spn3TDyb1 =   78
   INTEGER, PARAMETER             :: Spn3TDzb1 =   79
   INTEGER, PARAMETER             :: Spn4TDxb1 =   80
   INTEGER, PARAMETER             :: Spn4TDyb1 =   81
   INTEGER, PARAMETER             :: Spn4TDzb1 =   82
   INTEGER, PARAMETER             :: Spn5TDxb1 =   83
   INTEGER, PARAMETER             :: Spn5TDyb1 =   84
   INTEGER, PARAMETER             :: Spn5TDzb1 =   85
   INTEGER, PARAMETER             :: Spn6TDxb1 =   86
   INTEGER, PARAMETER             :: Spn6TDyb1 =   87
   INTEGER, PARAMETER             :: Spn6TDzb1 =   88
   INTEGER, PARAMETER             :: Spn7TDxb1 =   89
   INTEGER, PARAMETER             :: Spn7TDyb1 =   90
   INTEGER, PARAMETER             :: Spn7TDzb1 =   91
   INTEGER, PARAMETER             :: Spn8TDxb1 =   92
   INTEGER, PARAMETER             :: Spn8TDyb1 =   93
   INTEGER, PARAMETER             :: Spn8TDzb1 =   94
   INTEGER, PARAMETER             :: Spn9TDxb1 =   95
   INTEGER, PARAMETER             :: Spn9TDyb1 =   96
   INTEGER, PARAMETER             :: Spn9TDzb1 =   97
   INTEGER, PARAMETER             :: Spn1RDxb1 =   98
   INTEGER, PARAMETER             :: Spn1RDyb1 =   99
   INTEGER, PARAMETER             :: Spn1RDzb1 =  100
   INTEGER, PARAMETER             :: Spn2RDxb1 =  101
   INTEGER, PARAMETER             :: Spn2RDyb1 =  102
   INTEGER, PARAMETER             :: Spn2RDzb1 =  103
   INTEGER, PARAMETER             :: Spn3RDxb1 =  104
   INTEGER, PARAMETER             :: Spn3RDyb1 =  105
   INTEGER, PARAMETER             :: Spn3RDzb1 =  106
   INTEGER, PARAMETER             :: Spn4RDxb1 =  107
   INTEGER, PARAMETER             :: Spn4RDyb1 =  108
   INTEGER, PARAMETER             :: Spn4RDzb1 =  109
   INTEGER, PARAMETER             :: Spn5RDxb1 =  110
   INTEGER, PARAMETER             :: Spn5RDyb1 =  111
   INTEGER, PARAMETER             :: Spn5RDzb1 =  112
   INTEGER, PARAMETER             :: Spn6RDxb1 =  113
   INTEGER, PARAMETER             :: Spn6RDyb1 =  114
   INTEGER, PARAMETER             :: Spn6RDzb1 =  115
   INTEGER, PARAMETER             :: Spn7RDxb1 =  116
   INTEGER, PARAMETER             :: Spn7RDyb1 =  117
   INTEGER, PARAMETER             :: Spn7RDzb1 =  118
   INTEGER, PARAMETER             :: Spn8RDxb1 =  119
   INTEGER, PARAMETER             :: Spn8RDyb1 =  120
   INTEGER, PARAMETER             :: Spn8RDzb1 =  121
   INTEGER, PARAMETER             :: Spn9RDxb1 =  122
   INTEGER, PARAMETER             :: Spn9RDyb1 =  123
   INTEGER, PARAMETER             :: Spn9RDzb1 =  124


  ! Blade 2 Local Span Motions:

   INTEGER, PARAMETER             :: Spn1ALxb2 =  125
   INTEGER, PARAMETER             :: Spn1ALyb2 =  126
   INTEGER, PARAMETER             :: Spn1ALzb2 =  127
   INTEGER, PARAMETER             :: Spn2ALxb2 =  128
   INTEGER, PARAMETER             :: Spn2ALyb2 =  129
   INTEGER, PARAMETER             :: Spn2ALzb2 =  130
   INTEGER, PARAMETER             :: Spn3ALxb2 =  131
   INTEGER, PARAMETER             :: Spn3ALyb2 =  132
   INTEGER, PARAMETER             :: Spn3ALzb2 =  133
   INTEGER, PARAMETER             :: Spn4ALxb2 =  134
   INTEGER, PARAMETER             :: Spn4ALyb2 =  135
   INTEGER, PARAMETER             :: Spn4ALzb2 =  136
   INTEGER, PARAMETER             :: Spn5ALxb2 =  137
   INTEGER, PARAMETER             :: Spn5ALyb2 =  138
   INTEGER, PARAMETER             :: Spn5ALzb2 =  139
   INTEGER, PARAMETER             :: Spn6ALxb2 =  140
   INTEGER, PARAMETER             :: Spn6ALyb2 =  141
   INTEGER, PARAMETER             :: Spn6ALzb2 =  142
   INTEGER, PARAMETER             :: Spn7ALxb2 =  143
   INTEGER, PARAMETER             :: Spn7ALyb2 =  144
   INTEGER, PARAMETER             :: Spn7ALzb2 =  145
   INTEGER, PARAMETER             :: Spn8ALxb2 =  146
   INTEGER, PARAMETER             :: Spn8ALyb2 =  147
   INTEGER, PARAMETER             :: Spn8ALzb2 =  148
   INTEGER, PARAMETER             :: Spn9ALxb2 =  149
   INTEGER, PARAMETER             :: Spn9ALyb2 =  150
   INTEGER, PARAMETER             :: Spn9ALzb2 =  151
   INTEGER, PARAMETER             :: Spn1TDxb2 =  152
   INTEGER, PARAMETER             :: Spn1TDyb2 =  153
   INTEGER, PARAMETER             :: Spn1TDzb2 =  154
   INTEGER, PARAMETER             :: Spn2TDxb2 =  155
   INTEGER, PARAMETER             :: Spn2TDyb2 =  156
   INTEGER, PARAMETER             :: Spn2TDzb2 =  157
   INTEGER, PARAMETER             :: Spn3TDxb2 =  158
   INTEGER, PARAMETER             :: Spn3TDyb2 =  159
   INTEGER, PARAMETER             :: Spn3TDzb2 =  160
   INTEGER, PARAMETER             :: Spn4TDxb2 =  161
   INTEGER, PARAMETER             :: Spn4TDyb2 =  162
   INTEGER, PARAMETER             :: Spn4TDzb2 =  163
   INTEGER, PARAMETER             :: Spn5TDxb2 =  164
   INTEGER, PARAMETER             :: Spn5TDyb2 =  165
   INTEGER, PARAMETER             :: Spn5TDzb2 =  166
   INTEGER, PARAMETER             :: Spn6TDxb2 =  167
   INTEGER, PARAMETER             :: Spn6TDyb2 =  168
   INTEGER, PARAMETER             :: Spn6TDzb2 =  169
   INTEGER, PARAMETER             :: Spn7TDxb2 =  170
   INTEGER, PARAMETER             :: Spn7TDyb2 =  171
   INTEGER, PARAMETER             :: Spn7TDzb2 =  172
   INTEGER, PARAMETER             :: Spn8TDxb2 =  173
   INTEGER, PARAMETER             :: Spn8TDyb2 =  174
   INTEGER, PARAMETER             :: Spn8TDzb2 =  175
   INTEGER, PARAMETER             :: Spn9TDxb2 =  176
   INTEGER, PARAMETER             :: Spn9TDyb2 =  177
   INTEGER, PARAMETER             :: Spn9TDzb2 =  178
   INTEGER, PARAMETER             :: Spn1RDxb2 =  179
   INTEGER, PARAMETER             :: Spn1RDyb2 =  180
   INTEGER, PARAMETER             :: Spn1RDzb2 =  181
   INTEGER, PARAMETER             :: Spn2RDxb2 =  182
   INTEGER, PARAMETER             :: Spn2RDyb2 =  183
   INTEGER, PARAMETER             :: Spn2RDzb2 =  184
   INTEGER, PARAMETER             :: Spn3RDxb2 =  185
   INTEGER, PARAMETER             :: Spn3RDyb2 =  186
   INTEGER, PARAMETER             :: Spn3RDzb2 =  187
   INTEGER, PARAMETER             :: Spn4RDxb2 =  188
   INTEGER, PARAMETER             :: Spn4RDyb2 =  189
   INTEGER, PARAMETER             :: Spn4RDzb2 =  190
   INTEGER, PARAMETER             :: Spn5RDxb2 =  191
   INTEGER, PARAMETER             :: Spn5RDyb2 =  192
   INTEGER, PARAMETER             :: Spn5RDzb2 =  193
   INTEGER, PARAMETER             :: Spn6RDxb2 =  194
   INTEGER, PARAMETER             :: Spn6RDyb2 =  195
   INTEGER, PARAMETER             :: Spn6RDzb2 =  196
   INTEGER, PARAMETER             :: Spn7RDxb2 =  197
   INTEGER, PARAMETER             :: Spn7RDyb2 =  198
   INTEGER, PARAMETER             :: Spn7RDzb2 =  199
   INTEGER, PARAMETER             :: Spn8RDxb2 =  200
   INTEGER, PARAMETER             :: Spn8RDyb2 =  201
   INTEGER, PARAMETER             :: Spn8RDzb2 =  202
   INTEGER, PARAMETER             :: Spn9RDxb2 =  203
   INTEGER, PARAMETER             :: Spn9RDyb2 =  204
   INTEGER, PARAMETER             :: Spn9RDzb2 =  205


  ! Blade 3 Local Span Motions:

   INTEGER, PARAMETER             :: Spn1ALxb3 =  206
   INTEGER, PARAMETER             :: Spn1ALyb3 =  207
   INTEGER, PARAMETER             :: Spn1ALzb3 =  208
   INTEGER, PARAMETER             :: Spn2ALxb3 =  209
   INTEGER, PARAMETER             :: Spn2ALyb3 =  210
   INTEGER, PARAMETER             :: Spn2ALzb3 =  211
   INTEGER, PARAMETER             :: Spn3ALxb3 =  212
   INTEGER, PARAMETER             :: Spn3ALyb3 =  213
   INTEGER, PARAMETER             :: Spn3ALzb3 =  214
   INTEGER, PARAMETER             :: Spn4ALxb3 =  215
   INTEGER, PARAMETER             :: Spn4ALyb3 =  216
   INTEGER, PARAMETER             :: Spn4ALzb3 =  217
   INTEGER, PARAMETER             :: Spn5ALxb3 =  218
   INTEGER, PARAMETER             :: Spn5ALyb3 =  219
   INTEGER, PARAMETER             :: Spn5ALzb3 =  220
   INTEGER, PARAMETER             :: Spn6ALxb3 =  221
   INTEGER, PARAMETER             :: Spn6ALyb3 =  222
   INTEGER, PARAMETER             :: Spn6ALzb3 =  223
   INTEGER, PARAMETER             :: Spn7ALxb3 =  224
   INTEGER, PARAMETER             :: Spn7ALyb3 =  225
   INTEGER, PARAMETER             :: Spn7ALzb3 =  226
   INTEGER, PARAMETER             :: Spn8ALxb3 =  227
   INTEGER, PARAMETER             :: Spn8ALyb3 =  228
   INTEGER, PARAMETER             :: Spn8ALzb3 =  229
   INTEGER, PARAMETER             :: Spn9ALxb3 =  230
   INTEGER, PARAMETER             :: Spn9ALyb3 =  231
   INTEGER, PARAMETER             :: Spn9ALzb3 =  232
   INTEGER, PARAMETER             :: Spn1TDxb3 =  233
   INTEGER, PARAMETER             :: Spn1TDyb3 =  234
   INTEGER, PARAMETER             :: Spn1TDzb3 =  235
   INTEGER, PARAMETER             :: Spn2TDxb3 =  236
   INTEGER, PARAMETER             :: Spn2TDyb3 =  237
   INTEGER, PARAMETER             :: Spn2TDzb3 =  238
   INTEGER, PARAMETER             :: Spn3TDxb3 =  239
   INTEGER, PARAMETER             :: Spn3TDyb3 =  240
   INTEGER, PARAMETER             :: Spn3TDzb3 =  241
   INTEGER, PARAMETER             :: Spn4TDxb3 =  242
   INTEGER, PARAMETER             :: Spn4TDyb3 =  243
   INTEGER, PARAMETER             :: Spn4TDzb3 =  244
   INTEGER, PARAMETER             :: Spn5TDxb3 =  245
   INTEGER, PARAMETER             :: Spn5TDyb3 =  246
   INTEGER, PARAMETER             :: Spn5TDzb3 =  247
   INTEGER, PARAMETER             :: Spn6TDxb3 =  248
   INTEGER, PARAMETER             :: Spn6TDyb3 =  249
   INTEGER, PARAMETER             :: Spn6TDzb3 =  250
   INTEGER, PARAMETER             :: Spn7TDxb3 =  251
   INTEGER, PARAMETER             :: Spn7TDyb3 =  252
   INTEGER, PARAMETER             :: Spn7TDzb3 =  253
   INTEGER, PARAMETER             :: Spn8TDxb3 =  254
   INTEGER, PARAMETER             :: Spn8TDyb3 =  255
   INTEGER, PARAMETER             :: Spn8TDzb3 =  256
   INTEGER, PARAMETER             :: Spn9TDxb3 =  257
   INTEGER, PARAMETER             :: Spn9TDyb3 =  258
   INTEGER, PARAMETER             :: Spn9TDzb3 =  259
   INTEGER, PARAMETER             :: Spn1RDxb3 =  260
   INTEGER, PARAMETER             :: Spn1RDyb3 =  261
   INTEGER, PARAMETER             :: Spn1RDzb3 =  262
   INTEGER, PARAMETER             :: Spn2RDxb3 =  263
   INTEGER, PARAMETER             :: Spn2RDyb3 =  264
   INTEGER, PARAMETER             :: Spn2RDzb3 =  265
   INTEGER, PARAMETER             :: Spn3RDxb3 =  266
   INTEGER, PARAMETER             :: Spn3RDyb3 =  267
   INTEGER, PARAMETER             :: Spn3RDzb3 =  268
   INTEGER, PARAMETER             :: Spn4RDxb3 =  269
   INTEGER, PARAMETER             :: Spn4RDyb3 =  270
   INTEGER, PARAMETER             :: Spn4RDzb3 =  271
   INTEGER, PARAMETER             :: Spn5RDxb3 =  272
   INTEGER, PARAMETER             :: Spn5RDyb3 =  273
   INTEGER, PARAMETER             :: Spn5RDzb3 =  274
   INTEGER, PARAMETER             :: Spn6RDxb3 =  275
   INTEGER, PARAMETER             :: Spn6RDyb3 =  276
   INTEGER, PARAMETER             :: Spn6RDzb3 =  277
   INTEGER, PARAMETER             :: Spn7RDxb3 =  278
   INTEGER, PARAMETER             :: Spn7RDyb3 =  279
   INTEGER, PARAMETER             :: Spn7RDzb3 =  280
   INTEGER, PARAMETER             :: Spn8RDxb3 =  281
   INTEGER, PARAMETER             :: Spn8RDyb3 =  282
   INTEGER, PARAMETER             :: Spn8RDzb3 =  283
   INTEGER, PARAMETER             :: Spn9RDxb3 =  284
   INTEGER, PARAMETER             :: Spn9RDyb3 =  285
   INTEGER, PARAMETER             :: Spn9RDzb3 =  286


  ! Blade Pitch Motions:

   INTEGER, PARAMETER             :: PtchPMzc1 =  287
   INTEGER, PARAMETER             :: PtchPMzc2 =  288
   INTEGER, PARAMETER             :: PtchPMzc3 =  289


  ! Teeter Motions:

   INTEGER, PARAMETER             :: TeetPya   =  290
   INTEGER, PARAMETER             :: TeetVya   =  291
   INTEGER, PARAMETER             :: TeetAya   =  292


  ! Shaft Motions:

   INTEGER, PARAMETER             :: LSSTipPxa =  293
   INTEGER, PARAMETER             :: LSSTipVxa =  294
   INTEGER, PARAMETER             :: LSSTipAxa =  295
   INTEGER, PARAMETER             :: LSSGagPxa =  296
   INTEGER, PARAMETER             :: LSSGagVxa =  297
   INTEGER, PARAMETER             :: LSSGagAxa =  298
   INTEGER, PARAMETER             :: HSShftV   =  299
   INTEGER, PARAMETER             :: HSShftA   =  300
   INTEGER, PARAMETER             :: TipSpdRat =  301


  ! Nacelle IMU Motions:

   INTEGER, PARAMETER             :: NcIMUTVxs =  302
   INTEGER, PARAMETER             :: NcIMUTVys =  303
   INTEGER, PARAMETER             :: NcIMUTVzs =  304
   INTEGER, PARAMETER             :: NcIMUTAxs =  305
   INTEGER, PARAMETER             :: NcIMUTAys =  306
   INTEGER, PARAMETER             :: NcIMUTAzs =  307
   INTEGER, PARAMETER             :: NcIMURVxs =  308
   INTEGER, PARAMETER             :: NcIMURVys =  309
   INTEGER, PARAMETER             :: NcIMURVzs =  310
   INTEGER, PARAMETER             :: NcIMURAxs =  311
   INTEGER, PARAMETER             :: NcIMURAys =  312
   INTEGER, PARAMETER             :: NcIMURAzs =  313


  ! Rotor-Furl Motions:

   INTEGER, PARAMETER             :: RotFurlP  =  314
   INTEGER, PARAMETER             :: RotFurlV  =  315
   INTEGER, PARAMETER             :: RotFurlA  =  316


  ! Tail-Furl Motions:

   INTEGER, PARAMETER             :: TailFurlP =  317
   INTEGER, PARAMETER             :: TailFurlV =  318
   INTEGER, PARAMETER             :: TailFurlA =  319


  ! Nacelle Yaw Motions:

   INTEGER, PARAMETER             :: YawPzn    =  320
   INTEGER, PARAMETER             :: YawVzn    =  321
   INTEGER, PARAMETER             :: YawAzn    =  322
   INTEGER, PARAMETER             :: NacYawErr =  323


  ! Tower-Top / Yaw Bearing Motions:

   INTEGER, PARAMETER             :: YawBrTDxp =  324
   INTEGER, PARAMETER             :: YawBrTDyp =  325
   INTEGER, PARAMETER             :: YawBrTDzp =  326
   INTEGER, PARAMETER             :: YawBrTDxt =  327
   INTEGER, PARAMETER             :: YawBrTDyt =  328
   INTEGER, PARAMETER             :: YawBrTDzt =  329
   INTEGER, PARAMETER             :: YawBrTAxp =  330
   INTEGER, PARAMETER             :: YawBrTAyp =  331
   INTEGER, PARAMETER             :: YawBrTAzp =  332
   INTEGER, PARAMETER             :: YawBrRDxt =  333
   INTEGER, PARAMETER             :: YawBrRDyt =  334
   INTEGER, PARAMETER             :: YawBrRDzt =  335
   INTEGER, PARAMETER             :: YawBrRVxp =  336
   INTEGER, PARAMETER             :: YawBrRVyp =  337
   INTEGER, PARAMETER             :: YawBrRVzp =  338
   INTEGER, PARAMETER             :: YawBrRAxp =  339
   INTEGER, PARAMETER             :: YawBrRAyp =  340
   INTEGER, PARAMETER             :: YawBrRAzp =  341


  ! Local Tower Motions:

   INTEGER, PARAMETER             :: TwHt1ALxt =  342
   INTEGER, PARAMETER             :: TwHt1ALyt =  343
   INTEGER, PARAMETER             :: TwHt1ALzt =  344
   INTEGER, PARAMETER             :: TwHt2ALxt =  345
   INTEGER, PARAMETER             :: TwHt2ALyt =  346
   INTEGER, PARAMETER             :: TwHt2ALzt =  347
   INTEGER, PARAMETER             :: TwHt3ALxt =  348
   INTEGER, PARAMETER             :: TwHt3ALyt =  349
   INTEGER, PARAMETER             :: TwHt3ALzt =  350
   INTEGER, PARAMETER             :: TwHt4ALxt =  351
   INTEGER, PARAMETER             :: TwHt4ALyt =  352
   INTEGER, PARAMETER             :: TwHt4ALzt =  353
   INTEGER, PARAMETER             :: TwHt5ALxt =  354
   INTEGER, PARAMETER             :: TwHt5ALyt =  355
   INTEGER, PARAMETER             :: TwHt5ALzt =  356
   INTEGER, PARAMETER             :: TwHt6ALxt =  357
   INTEGER, PARAMETER             :: TwHt6ALyt =  358
   INTEGER, PARAMETER             :: TwHt6ALzt =  359
   INTEGER, PARAMETER             :: TwHt7ALxt =  360
   INTEGER, PARAMETER             :: TwHt7ALyt =  361
   INTEGER, PARAMETER             :: TwHt7ALzt =  362
   INTEGER, PARAMETER             :: TwHt8ALxt =  363
   INTEGER, PARAMETER             :: TwHt8ALyt =  364
   INTEGER, PARAMETER             :: TwHt8ALzt =  365
   INTEGER, PARAMETER             :: TwHt9ALxt =  366
   INTEGER, PARAMETER             :: TwHt9ALyt =  367
   INTEGER, PARAMETER             :: TwHt9ALzt =  368
   INTEGER, PARAMETER             :: TwHt1TDxt =  369
   INTEGER, PARAMETER             :: TwHt1TDyt =  370
   INTEGER, PARAMETER             :: TwHt1TDzt =  371
   INTEGER, PARAMETER             :: TwHt2TDxt =  372
   INTEGER, PARAMETER             :: TwHt2TDyt =  373
   INTEGER, PARAMETER             :: TwHt2TDzt =  374
   INTEGER, PARAMETER             :: TwHt3TDxt =  375
   INTEGER, PARAMETER             :: TwHt3TDyt =  376
   INTEGER, PARAMETER             :: TwHt3TDzt =  377
   INTEGER, PARAMETER             :: TwHt4TDxt =  378
   INTEGER, PARAMETER             :: TwHt4TDyt =  379
   INTEGER, PARAMETER             :: TwHt4TDzt =  380
   INTEGER, PARAMETER             :: TwHt5TDxt =  381
   INTEGER, PARAMETER             :: TwHt5TDyt =  382
   INTEGER, PARAMETER             :: TwHt5TDzt =  383
   INTEGER, PARAMETER             :: TwHt6TDxt =  384
   INTEGER, PARAMETER             :: TwHt6TDyt =  385
   INTEGER, PARAMETER             :: TwHt6TDzt =  386
   INTEGER, PARAMETER             :: TwHt7TDxt =  387
   INTEGER, PARAMETER             :: TwHt7TDyt =  388
   INTEGER, PARAMETER             :: TwHt7TDzt =  389
   INTEGER, PARAMETER             :: TwHt8TDxt =  390
   INTEGER, PARAMETER             :: TwHt8TDyt =  391
   INTEGER, PARAMETER             :: TwHt8TDzt =  392
   INTEGER, PARAMETER             :: TwHt9TDxt =  393
   INTEGER, PARAMETER             :: TwHt9TDyt =  394
   INTEGER, PARAMETER             :: TwHt9TDzt =  395
   INTEGER, PARAMETER             :: TwHt1RDxt =  396
   INTEGER, PARAMETER             :: TwHt1RDyt =  397
   INTEGER, PARAMETER             :: TwHt1RDzt =  398
   INTEGER, PARAMETER             :: TwHt2RDxt =  399
   INTEGER, PARAMETER             :: TwHt2RDyt =  400
   INTEGER, PARAMETER             :: TwHt2RDzt =  401
   INTEGER, PARAMETER             :: TwHt3RDxt =  402
   INTEGER, PARAMETER             :: TwHt3RDyt =  403
   INTEGER, PARAMETER             :: TwHt3RDzt =  404
   INTEGER, PARAMETER             :: TwHt4RDxt =  405
   INTEGER, PARAMETER             :: TwHt4RDyt =  406
   INTEGER, PARAMETER             :: TwHt4RDzt =  407
   INTEGER, PARAMETER             :: TwHt5RDxt =  408
   INTEGER, PARAMETER             :: TwHt5RDyt =  409
   INTEGER, PARAMETER             :: TwHt5RDzt =  410
   INTEGER, PARAMETER             :: TwHt6RDxt =  411
   INTEGER, PARAMETER             :: TwHt6RDyt =  412
   INTEGER, PARAMETER             :: TwHt6RDzt =  413
   INTEGER, PARAMETER             :: TwHt7RDxt =  414
   INTEGER, PARAMETER             :: TwHt7RDyt =  415
   INTEGER, PARAMETER             :: TwHt7RDzt =  416
   INTEGER, PARAMETER             :: TwHt8RDxt =  417
   INTEGER, PARAMETER             :: TwHt8RDyt =  418
   INTEGER, PARAMETER             :: TwHt8RDzt =  419
   INTEGER, PARAMETER             :: TwHt9RDxt =  420
   INTEGER, PARAMETER             :: TwHt9RDyt =  421
   INTEGER, PARAMETER             :: TwHt9RDzt =  422
   INTEGER, PARAMETER             :: TwHt1TPxi =  423
   INTEGER, PARAMETER             :: TwHt1TPyi =  424
   INTEGER, PARAMETER             :: TwHt1TPzi =  425
   INTEGER, PARAMETER             :: TwHt2TPxi =  426
   INTEGER, PARAMETER             :: TwHt2TPyi =  427
   INTEGER, PARAMETER             :: TwHt2TPzi =  428
   INTEGER, PARAMETER             :: TwHt3TPxi =  429
   INTEGER, PARAMETER             :: TwHt3TPyi =  430
   INTEGER, PARAMETER             :: TwHt3TPzi =  431
   INTEGER, PARAMETER             :: TwHt4TPxi =  432
   INTEGER, PARAMETER             :: TwHt4TPyi =  433
   INTEGER, PARAMETER             :: TwHt4TPzi =  434
   INTEGER, PARAMETER             :: TwHt5TPxi =  435
   INTEGER, PARAMETER             :: TwHt5TPyi =  436
   INTEGER, PARAMETER             :: TwHt5TPzi =  437
   INTEGER, PARAMETER             :: TwHt6TPxi =  438
   INTEGER, PARAMETER             :: TwHt6TPyi =  439
   INTEGER, PARAMETER             :: TwHt6TPzi =  440
   INTEGER, PARAMETER             :: TwHt7TPxi =  441
   INTEGER, PARAMETER             :: TwHt7TPyi =  442
   INTEGER, PARAMETER             :: TwHt7TPzi =  443
   INTEGER, PARAMETER             :: TwHt8TPxi =  444
   INTEGER, PARAMETER             :: TwHt8TPyi =  445
   INTEGER, PARAMETER             :: TwHt8TPzi =  446
   INTEGER, PARAMETER             :: TwHt9TPxi =  447
   INTEGER, PARAMETER             :: TwHt9TPyi =  448
   INTEGER, PARAMETER             :: TwHt9TPzi =  449
   INTEGER, PARAMETER             :: TwHt1RPxi =  450
   INTEGER, PARAMETER             :: TwHt1RPyi =  451
   INTEGER, PARAMETER             :: TwHt1RPzi =  452
   INTEGER, PARAMETER             :: TwHt2RPxi =  453
   INTEGER, PARAMETER             :: TwHt2RPyi =  454
   INTEGER, PARAMETER             :: TwHt2RPzi =  455
   INTEGER, PARAMETER             :: TwHt3RPxi =  456
   INTEGER, PARAMETER             :: TwHt3RPyi =  457
   INTEGER, PARAMETER             :: TwHt3RPzi =  458
   INTEGER, PARAMETER             :: TwHt4RPxi =  459
   INTEGER, PARAMETER             :: TwHt4RPyi =  460
   INTEGER, PARAMETER             :: TwHt4RPzi =  461
   INTEGER, PARAMETER             :: TwHt5RPxi =  462
   INTEGER, PARAMETER             :: TwHt5RPyi =  463
   INTEGER, PARAMETER             :: TwHt5RPzi =  464
   INTEGER, PARAMETER             :: TwHt6RPxi =  465
   INTEGER, PARAMETER             :: TwHt6RPyi =  466
   INTEGER, PARAMETER             :: TwHt6RPzi =  467
   INTEGER, PARAMETER             :: TwHt7RPxi =  468
   INTEGER, PARAMETER             :: TwHt7RPyi =  469
   INTEGER, PARAMETER             :: TwHt7RPzi =  470
   INTEGER, PARAMETER             :: TwHt8RPxi =  471
   INTEGER, PARAMETER             :: TwHt8RPyi =  472
   INTEGER, PARAMETER             :: TwHt8RPzi =  473
   INTEGER, PARAMETER             :: TwHt9RPxi =  474
   INTEGER, PARAMETER             :: TwHt9RPyi =  475
   INTEGER, PARAMETER             :: TwHt9RPzi =  476


  ! Platform Motions:

   INTEGER, PARAMETER             :: PtfmTDxt  =  477
   INTEGER, PARAMETER             :: PtfmTDyt  =  478
   INTEGER, PARAMETER             :: PtfmTDzt  =  479
   INTEGER, PARAMETER             :: PtfmTDxi  =  480
   INTEGER, PARAMETER             :: PtfmTDyi  =  481
   INTEGER, PARAMETER             :: PtfmTDzi  =  482
   INTEGER, PARAMETER             :: PtfmTVxt  =  483
   INTEGER, PARAMETER             :: PtfmTVyt  =  484
   INTEGER, PARAMETER             :: PtfmTVzt  =  485
   INTEGER, PARAMETER             :: PtfmTVxi  =  486
   INTEGER, PARAMETER             :: PtfmTVyi  =  487
   INTEGER, PARAMETER             :: PtfmTVzi  =  488
   INTEGER, PARAMETER             :: PtfmTAxt  =  489
   INTEGER, PARAMETER             :: PtfmTAyt  =  490
   INTEGER, PARAMETER             :: PtfmTAzt  =  491
   INTEGER, PARAMETER             :: PtfmTAxi  =  492
   INTEGER, PARAMETER             :: PtfmTAyi  =  493
   INTEGER, PARAMETER             :: PtfmTAzi  =  494
   INTEGER, PARAMETER             :: PtfmRDxi  =  495
   INTEGER, PARAMETER             :: PtfmRDyi  =  496
   INTEGER, PARAMETER             :: PtfmRDzi  =  497
   INTEGER, PARAMETER             :: PtfmRVxt  =  498
   INTEGER, PARAMETER             :: PtfmRVyt  =  499
   INTEGER, PARAMETER             :: PtfmRVzt  =  500
   INTEGER, PARAMETER             :: PtfmRVxi  =  501
   INTEGER, PARAMETER             :: PtfmRVyi  =  502
   INTEGER, PARAMETER             :: PtfmRVzi  =  503
   INTEGER, PARAMETER             :: PtfmRAxt  =  504
   INTEGER, PARAMETER             :: PtfmRAyt  =  505
   INTEGER, PARAMETER             :: PtfmRAzt  =  506
   INTEGER, PARAMETER             :: PtfmRAxi  =  507
   INTEGER, PARAMETER             :: PtfmRAyi  =  508
   INTEGER, PARAMETER             :: PtfmRAzi  =  509


  ! Blade 1 Root Loads:

   INTEGER, PARAMETER             :: RootFxc1  =  510
   INTEGER, PARAMETER             :: RootFyc1  =  511
   INTEGER, PARAMETER             :: RootFzc1  =  512
   INTEGER, PARAMETER             :: RootFxb1  =  513
   INTEGER, PARAMETER             :: RootFyb1  =  514
   INTEGER, PARAMETER             :: RootMxc1  =  515
   INTEGER, PARAMETER             :: RootMyc1  =  516
   INTEGER, PARAMETER             :: RootMzc1  =  517
   INTEGER, PARAMETER             :: RootMxb1  =  518
   INTEGER, PARAMETER             :: RootMyb1  =  519


  ! Blade 2 Root Loads:

   INTEGER, PARAMETER             :: RootFxc2  =  520
   INTEGER, PARAMETER             :: RootFyc2  =  521
   INTEGER, PARAMETER             :: RootFzc2  =  522
   INTEGER, PARAMETER             :: RootFxb2  =  523
   INTEGER, PARAMETER             :: RootFyb2  =  524
   INTEGER, PARAMETER             :: RootMxc2  =  525
   INTEGER, PARAMETER             :: RootMyc2  =  526
   INTEGER, PARAMETER             :: RootMzc2  =  527
   INTEGER, PARAMETER             :: RootMxb2  =  528
   INTEGER, PARAMETER             :: RootMyb2  =  529


  ! Blade 3 Root Loads:

   INTEGER, PARAMETER             :: RootFxc3  =  530
   INTEGER, PARAMETER             :: RootFyc3  =  531
   INTEGER, PARAMETER             :: RootFzc3  =  532
   INTEGER, PARAMETER             :: RootFxb3  =  533
   INTEGER, PARAMETER             :: RootFyb3  =  534
   INTEGER, PARAMETER             :: RootMxc3  =  535
   INTEGER, PARAMETER             :: RootMyc3  =  536
   INTEGER, PARAMETER             :: RootMzc3  =  537
   INTEGER, PARAMETER             :: RootMxb3  =  538
   INTEGER, PARAMETER             :: RootMyb3  =  539


  ! Blade 1 Local Span Loads:

   INTEGER, PARAMETER             :: Spn1MLxb1 =  540
   INTEGER, PARAMETER             :: Spn1MLyb1 =  541
   INTEGER, PARAMETER             :: Spn1MLzb1 =  542
   INTEGER, PARAMETER             :: Spn2MLxb1 =  543
   INTEGER, PARAMETER             :: Spn2MLyb1 =  544
   INTEGER, PARAMETER             :: Spn2MLzb1 =  545
   INTEGER, PARAMETER             :: Spn3MLxb1 =  546
   INTEGER, PARAMETER             :: Spn3MLyb1 =  547
   INTEGER, PARAMETER             :: Spn3MLzb1 =  548
   INTEGER, PARAMETER             :: Spn4MLxb1 =  549
   INTEGER, PARAMETER             :: Spn4MLyb1 =  550
   INTEGER, PARAMETER             :: Spn4MLzb1 =  551
   INTEGER, PARAMETER             :: Spn5MLxb1 =  552
   INTEGER, PARAMETER             :: Spn5MLyb1 =  553
   INTEGER, PARAMETER             :: Spn5MLzb1 =  554
   INTEGER, PARAMETER             :: Spn6MLxb1 =  555
   INTEGER, PARAMETER             :: Spn6MLyb1 =  556
   INTEGER, PARAMETER             :: Spn6MLzb1 =  557
   INTEGER, PARAMETER             :: Spn7MLxb1 =  558
   INTEGER, PARAMETER             :: Spn7MLyb1 =  559
   INTEGER, PARAMETER             :: Spn7MLzb1 =  560
   INTEGER, PARAMETER             :: Spn8MLxb1 =  561
   INTEGER, PARAMETER             :: Spn8MLyb1 =  562
   INTEGER, PARAMETER             :: Spn8MLzb1 =  563
   INTEGER, PARAMETER             :: Spn9MLxb1 =  564
   INTEGER, PARAMETER             :: Spn9MLyb1 =  565
   INTEGER, PARAMETER             :: Spn9MLzb1 =  566
   INTEGER, PARAMETER             :: Spn1FLxb1 =  567
   INTEGER, PARAMETER             :: Spn1FLyb1 =  568
   INTEGER, PARAMETER             :: Spn1FLzb1 =  569
   INTEGER, PARAMETER             :: Spn2FLxb1 =  570
   INTEGER, PARAMETER             :: Spn2FLyb1 =  571
   INTEGER, PARAMETER             :: Spn2FLzb1 =  572
   INTEGER, PARAMETER             :: Spn3FLxb1 =  573
   INTEGER, PARAMETER             :: Spn3FLyb1 =  574
   INTEGER, PARAMETER             :: Spn3FLzb1 =  575
   INTEGER, PARAMETER             :: Spn4FLxb1 =  576
   INTEGER, PARAMETER             :: Spn4FLyb1 =  577
   INTEGER, PARAMETER             :: Spn4FLzb1 =  578
   INTEGER, PARAMETER             :: Spn5FLxb1 =  579
   INTEGER, PARAMETER             :: Spn5FLyb1 =  580
   INTEGER, PARAMETER             :: Spn5FLzb1 =  581
   INTEGER, PARAMETER             :: Spn6FLxb1 =  582
   INTEGER, PARAMETER             :: Spn6FLyb1 =  583
   INTEGER, PARAMETER             :: Spn6FLzb1 =  584
   INTEGER, PARAMETER             :: Spn7FLxb1 =  585
   INTEGER, PARAMETER             :: Spn7FLyb1 =  586
   INTEGER, PARAMETER             :: Spn7FLzb1 =  587
   INTEGER, PARAMETER             :: Spn8FLxb1 =  588
   INTEGER, PARAMETER             :: Spn8FLyb1 =  589
   INTEGER, PARAMETER             :: Spn8FLzb1 =  590
   INTEGER, PARAMETER             :: Spn9FLxb1 =  591
   INTEGER, PARAMETER             :: Spn9FLyb1 =  592
   INTEGER, PARAMETER             :: Spn9FLzb1 =  593


  ! Blade 2 Local Span Loads:

   INTEGER, PARAMETER             :: Spn1MLxb2 =  594
   INTEGER, PARAMETER             :: Spn1MLyb2 =  595
   INTEGER, PARAMETER             :: Spn1MLzb2 =  596
   INTEGER, PARAMETER             :: Spn2MLxb2 =  597
   INTEGER, PARAMETER             :: Spn2MLyb2 =  598
   INTEGER, PARAMETER             :: Spn2MLzb2 =  599
   INTEGER, PARAMETER             :: Spn3MLxb2 =  600
   INTEGER, PARAMETER             :: Spn3MLyb2 =  601
   INTEGER, PARAMETER             :: Spn3MLzb2 =  602
   INTEGER, PARAMETER             :: Spn4MLxb2 =  603
   INTEGER, PARAMETER             :: Spn4MLyb2 =  604
   INTEGER, PARAMETER             :: Spn4MLzb2 =  605
   INTEGER, PARAMETER             :: Spn5MLxb2 =  606
   INTEGER, PARAMETER             :: Spn5MLyb2 =  607
   INTEGER, PARAMETER             :: Spn5MLzb2 =  608
   INTEGER, PARAMETER             :: Spn6MLxb2 =  609
   INTEGER, PARAMETER             :: Spn6MLyb2 =  610
   INTEGER, PARAMETER             :: Spn6MLzb2 =  611
   INTEGER, PARAMETER             :: Spn7MLxb2 =  612
   INTEGER, PARAMETER             :: Spn7MLyb2 =  613
   INTEGER, PARAMETER             :: Spn7MLzb2 =  614
   INTEGER, PARAMETER             :: Spn8MLxb2 =  615
   INTEGER, PARAMETER             :: Spn8MLyb2 =  616
   INTEGER, PARAMETER             :: Spn8MLzb2 =  617
   INTEGER, PARAMETER             :: Spn9MLxb2 =  618
   INTEGER, PARAMETER             :: Spn9MLyb2 =  619
   INTEGER, PARAMETER             :: Spn9MLzb2 =  620
   INTEGER, PARAMETER             :: Spn1FLxb2 =  621
   INTEGER, PARAMETER             :: Spn1FLyb2 =  622
   INTEGER, PARAMETER             :: Spn1FLzb2 =  623
   INTEGER, PARAMETER             :: Spn2FLxb2 =  624
   INTEGER, PARAMETER             :: Spn2FLyb2 =  625
   INTEGER, PARAMETER             :: Spn2FLzb2 =  626
   INTEGER, PARAMETER             :: Spn3FLxb2 =  627
   INTEGER, PARAMETER             :: Spn3FLyb2 =  628
   INTEGER, PARAMETER             :: Spn3FLzb2 =  629
   INTEGER, PARAMETER             :: Spn4FLxb2 =  630
   INTEGER, PARAMETER             :: Spn4FLyb2 =  631
   INTEGER, PARAMETER             :: Spn4FLzb2 =  632
   INTEGER, PARAMETER             :: Spn5FLxb2 =  633
   INTEGER, PARAMETER             :: Spn5FLyb2 =  634
   INTEGER, PARAMETER             :: Spn5FLzb2 =  635
   INTEGER, PARAMETER             :: Spn6FLxb2 =  636
   INTEGER, PARAMETER             :: Spn6FLyb2 =  637
   INTEGER, PARAMETER             :: Spn6FLzb2 =  638
   INTEGER, PARAMETER             :: Spn7FLxb2 =  639
   INTEGER, PARAMETER             :: Spn7FLyb2 =  640
   INTEGER, PARAMETER             :: Spn7FLzb2 =  641
   INTEGER, PARAMETER             :: Spn8FLxb2 =  642
   INTEGER, PARAMETER             :: Spn8FLyb2 =  643
   INTEGER, PARAMETER             :: Spn8FLzb2 =  644
   INTEGER, PARAMETER             :: Spn9FLxb2 =  645
   INTEGER, PARAMETER             :: Spn9FLyb2 =  646
   INTEGER, PARAMETER             :: Spn9FLzb2 =  647


  ! Blade 3 Local Span Loads:

   INTEGER, PARAMETER             :: Spn1MLxb3 =  648
   INTEGER, PARAMETER             :: Spn1MLyb3 =  649
   INTEGER, PARAMETER             :: Spn1MLzb3 =  650
   INTEGER, PARAMETER             :: Spn2MLxb3 =  651
   INTEGER, PARAMETER             :: Spn2MLyb3 =  652
   INTEGER, PARAMETER             :: Spn2MLzb3 =  653
   INTEGER, PARAMETER             :: Spn3MLxb3 =  654
   INTEGER, PARAMETER             :: Spn3MLyb3 =  655
   INTEGER, PARAMETER             :: Spn3MLzb3 =  656
   INTEGER, PARAMETER             :: Spn4MLxb3 =  657
   INTEGER, PARAMETER             :: Spn4MLyb3 =  658
   INTEGER, PARAMETER             :: Spn4MLzb3 =  659
   INTEGER, PARAMETER             :: Spn5MLxb3 =  660
   INTEGER, PARAMETER             :: Spn5MLyb3 =  661
   INTEGER, PARAMETER             :: Spn5MLzb3 =  662
   INTEGER, PARAMETER             :: Spn6MLxb3 =  663
   INTEGER, PARAMETER             :: Spn6MLyb3 =  664
   INTEGER, PARAMETER             :: Spn6MLzb3 =  665
   INTEGER, PARAMETER             :: Spn7MLxb3 =  666
   INTEGER, PARAMETER             :: Spn7MLyb3 =  667
   INTEGER, PARAMETER             :: Spn7MLzb3 =  668
   INTEGER, PARAMETER             :: Spn8MLxb3 =  669
   INTEGER, PARAMETER             :: Spn8MLyb3 =  670
   INTEGER, PARAMETER             :: Spn8MLzb3 =  671
   INTEGER, PARAMETER             :: Spn9MLxb3 =  672
   INTEGER, PARAMETER             :: Spn9MLyb3 =  673
   INTEGER, PARAMETER             :: Spn9MLzb3 =  674
   INTEGER, PARAMETER             :: Spn1FLxb3 =  675
   INTEGER, PARAMETER             :: Spn1FLyb3 =  676
   INTEGER, PARAMETER             :: Spn1FLzb3 =  677
   INTEGER, PARAMETER             :: Spn2FLxb3 =  678
   INTEGER, PARAMETER             :: Spn2FLyb3 =  679
   INTEGER, PARAMETER             :: Spn2FLzb3 =  680
   INTEGER, PARAMETER             :: Spn3FLxb3 =  681
   INTEGER, PARAMETER             :: Spn3FLyb3 =  682
   INTEGER, PARAMETER             :: Spn3FLzb3 =  683
   INTEGER, PARAMETER             :: Spn4FLxb3 =  684
   INTEGER, PARAMETER             :: Spn4FLyb3 =  685
   INTEGER, PARAMETER             :: Spn4FLzb3 =  686
   INTEGER, PARAMETER             :: Spn5FLxb3 =  687
   INTEGER, PARAMETER             :: Spn5FLyb3 =  688
   INTEGER, PARAMETER             :: Spn5FLzb3 =  689
   INTEGER, PARAMETER             :: Spn6FLxb3 =  690
   INTEGER, PARAMETER             :: Spn6FLyb3 =  691
   INTEGER, PARAMETER             :: Spn6FLzb3 =  692
   INTEGER, PARAMETER             :: Spn7FLxb3 =  693
   INTEGER, PARAMETER             :: Spn7FLyb3 =  694
   INTEGER, PARAMETER             :: Spn7FLzb3 =  695
   INTEGER, PARAMETER             :: Spn8FLxb3 =  696
   INTEGER, PARAMETER             :: Spn8FLyb3 =  697
   INTEGER, PARAMETER             :: Spn8FLzb3 =  698
   INTEGER, PARAMETER             :: Spn9FLxb3 =  699
   INTEGER, PARAMETER             :: Spn9FLyb3 =  700
   INTEGER, PARAMETER             :: Spn9FLzb3 =  701


  ! Hub and Rotor Loads:

   INTEGER, PARAMETER             :: LSShftFxa =  702
   INTEGER, PARAMETER             :: LSShftFya =  703
   INTEGER, PARAMETER             :: LSShftFza =  704
   INTEGER, PARAMETER             :: LSShftFys =  705
   INTEGER, PARAMETER             :: LSShftFzs =  706
   INTEGER, PARAMETER             :: LSShftMxa =  707
   INTEGER, PARAMETER             :: LSSTipMya =  708
   INTEGER, PARAMETER             :: LSSTipMza =  709
   INTEGER, PARAMETER             :: LSSTipMys =  710
   INTEGER, PARAMETER             :: LSSTipMzs =  711
   INTEGER, PARAMETER             :: CThrstAzm =  712
   INTEGER, PARAMETER             :: CThrstRad =  713
   INTEGER, PARAMETER             :: RotPwr    =  714
   INTEGER, PARAMETER             :: RotCq     =  715
   INTEGER, PARAMETER             :: RotCp     =  716
   INTEGER, PARAMETER             :: RotCt     =  717


  ! Shaft Strain Gage Loads:

   INTEGER, PARAMETER             :: LSSGagMya =  718
   INTEGER, PARAMETER             :: LSSGagMza =  719
   INTEGER, PARAMETER             :: LSSGagMys =  720
   INTEGER, PARAMETER             :: LSSGagMzs =  721


  ! Generator and High-Speed Shaft Loads:

   INTEGER, PARAMETER             :: HSShftTq  =  722
   INTEGER, PARAMETER             :: HSShftPwr =  723
   INTEGER, PARAMETER             :: HSShftCq  =  724
   INTEGER, PARAMETER             :: HSShftCp  =  725
   INTEGER, PARAMETER             :: GenTq     =  726
   INTEGER, PARAMETER             :: GenPwr    =  727
   INTEGER, PARAMETER             :: GenCq     =  728
   INTEGER, PARAMETER             :: GenCp     =  729
   INTEGER, PARAMETER             :: HSSBrTq   =  730


  ! Rotor-Furl Bearing Loads:

   INTEGER, PARAMETER             :: RFrlBrM   =  731


  ! Tail-Furl Bearing Loads:

   INTEGER, PARAMETER             :: TFrlBrM   =  732


  ! Tail Fin Aerodynamic Loads:

   INTEGER, PARAMETER             :: TFinAlpha =  733
   INTEGER, PARAMETER             :: TFinCLift =  734
   INTEGER, PARAMETER             :: TFinCDrag =  735
   INTEGER, PARAMETER             :: TFinDnPrs =  736
   INTEGER, PARAMETER             :: TFinCPFx  =  737
   INTEGER, PARAMETER             :: TFinCPFy  =  738


  ! Tower-Top / Yaw Bearing Loads:

   INTEGER, PARAMETER             :: YawBrFxn  =  739
   INTEGER, PARAMETER             :: YawBrFyn  =  740
   INTEGER, PARAMETER             :: YawBrFzn  =  741
   INTEGER, PARAMETER             :: YawBrFxp  =  742
   INTEGER, PARAMETER             :: YawBrFyp  =  743
   INTEGER, PARAMETER             :: YawBrMxn  =  744
   INTEGER, PARAMETER             :: YawBrMyn  =  745
   INTEGER, PARAMETER             :: YawBrMzn  =  746
   INTEGER, PARAMETER             :: YawBrMxp  =  747
   INTEGER, PARAMETER             :: YawBrMyp  =  748


  ! Tower Base Loads:

   INTEGER, PARAMETER             :: TwrBsFxt  =  749
   INTEGER, PARAMETER             :: TwrBsFyt  =  750
   INTEGER, PARAMETER             :: TwrBsFzt  =  751
   INTEGER, PARAMETER             :: TwrBsMxt  =  752
   INTEGER, PARAMETER             :: TwrBsMyt  =  753
   INTEGER, PARAMETER             :: TwrBsMzt  =  754


  ! Local Tower Loads:

   INTEGER, PARAMETER             :: TwHt1MLxt =  755
   INTEGER, PARAMETER             :: TwHt1MLyt =  756
   INTEGER, PARAMETER             :: TwHt1MLzt =  757
   INTEGER, PARAMETER             :: TwHt2MLxt =  758
   INTEGER, PARAMETER             :: TwHt2MLyt =  759
   INTEGER, PARAMETER             :: TwHt2MLzt =  760
   INTEGER, PARAMETER             :: TwHt3MLxt =  761
   INTEGER, PARAMETER             :: TwHt3MLyt =  762
   INTEGER, PARAMETER             :: TwHt3MLzt =  763
   INTEGER, PARAMETER             :: TwHt4MLxt =  764
   INTEGER, PARAMETER             :: TwHt4MLyt =  765
   INTEGER, PARAMETER             :: TwHt4MLzt =  766
   INTEGER, PARAMETER             :: TwHt5MLxt =  767
   INTEGER, PARAMETER             :: TwHt5MLyt =  768
   INTEGER, PARAMETER             :: TwHt5MLzt =  769
   INTEGER, PARAMETER             :: TwHt6MLxt =  770
   INTEGER, PARAMETER             :: TwHt6MLyt =  771
   INTEGER, PARAMETER             :: TwHt6MLzt =  772
   INTEGER, PARAMETER             :: TwHt7MLxt =  773
   INTEGER, PARAMETER             :: TwHt7MLyt =  774
   INTEGER, PARAMETER             :: TwHt7MLzt =  775
   INTEGER, PARAMETER             :: TwHt8MLxt =  776
   INTEGER, PARAMETER             :: TwHt8MLyt =  777
   INTEGER, PARAMETER             :: TwHt8MLzt =  778
   INTEGER, PARAMETER             :: TwHt9MLxt =  779
   INTEGER, PARAMETER             :: TwHt9MLyt =  780
   INTEGER, PARAMETER             :: TwHt9MLzt =  781
   INTEGER, PARAMETER             :: TwHt1FLxt =  782
   INTEGER, PARAMETER             :: TwHt1FLyt =  783
   INTEGER, PARAMETER             :: TwHt1FLzt =  784
   INTEGER, PARAMETER             :: TwHt2FLxt =  785
   INTEGER, PARAMETER             :: TwHt2FLyt =  786
   INTEGER, PARAMETER             :: TwHt2FLzt =  787
   INTEGER, PARAMETER             :: TwHt3FLxt =  788
   INTEGER, PARAMETER             :: TwHt3FLyt =  789
   INTEGER, PARAMETER             :: TwHt3FLzt =  790
   INTEGER, PARAMETER             :: TwHt4FLxt =  791
   INTEGER, PARAMETER             :: TwHt4FLyt =  792
   INTEGER, PARAMETER             :: TwHt4FLzt =  793
   INTEGER, PARAMETER             :: TwHt5FLxt =  794
   INTEGER, PARAMETER             :: TwHt5FLyt =  795
   INTEGER, PARAMETER             :: TwHt5FLzt =  796
   INTEGER, PARAMETER             :: TwHt6FLxt =  797
   INTEGER, PARAMETER             :: TwHt6FLyt =  798
   INTEGER, PARAMETER             :: TwHt6FLzt =  799
   INTEGER, PARAMETER             :: TwHt7FLxt =  800
   INTEGER, PARAMETER             :: TwHt7FLyt =  801
   INTEGER, PARAMETER             :: TwHt7FLzt =  802
   INTEGER, PARAMETER             :: TwHt8FLxt =  803
   INTEGER, PARAMETER             :: TwHt8FLyt =  804
   INTEGER, PARAMETER             :: TwHt8FLzt =  805
   INTEGER, PARAMETER             :: TwHt9FLxt =  806
   INTEGER, PARAMETER             :: TwHt9FLyt =  807
   INTEGER, PARAMETER             :: TwHt9FLzt =  808


  ! Platform Loads:

   INTEGER, PARAMETER             :: PtfmFxt   =  809
   INTEGER, PARAMETER             :: PtfmFyt   =  810
   INTEGER, PARAMETER             :: PtfmFzt   =  811
   INTEGER, PARAMETER             :: PtfmFxi   =  812
   INTEGER, PARAMETER             :: PtfmFyi   =  813
   INTEGER, PARAMETER             :: PtfmFzi   =  814
   INTEGER, PARAMETER             :: PtfmMxt   =  815
   INTEGER, PARAMETER             :: PtfmMyt   =  816
   INTEGER, PARAMETER             :: PtfmMzt   =  817
   INTEGER, PARAMETER             :: PtfmMxi   =  818
   INTEGER, PARAMETER             :: PtfmMyi   =  819
   INTEGER, PARAMETER             :: PtfmMzi   =  820


  ! Mooring Line Loads:

   INTEGER, PARAMETER             :: Fair1Ten  =  821
   INTEGER, PARAMETER             :: Fair1Ang  =  822
   INTEGER, PARAMETER             :: Anch1Ten  =  823
   INTEGER, PARAMETER             :: Anch1Ang  =  824
   INTEGER, PARAMETER             :: Fair2Ten  =  825
   INTEGER, PARAMETER             :: Fair2Ang  =  826
   INTEGER, PARAMETER             :: Anch2Ten  =  827
   INTEGER, PARAMETER             :: Anch2Ang  =  828
   INTEGER, PARAMETER             :: Fair3Ten  =  829
   INTEGER, PARAMETER             :: Fair3Ang  =  830
   INTEGER, PARAMETER             :: Anch3Ten  =  831
   INTEGER, PARAMETER             :: Anch3Ang  =  832
   INTEGER, PARAMETER             :: Fair4Ten  =  833
   INTEGER, PARAMETER             :: Fair4Ang  =  834
   INTEGER, PARAMETER             :: Anch4Ten  =  835
   INTEGER, PARAMETER             :: Anch4Ang  =  836
   INTEGER, PARAMETER             :: Fair5Ten  =  837
   INTEGER, PARAMETER             :: Fair5Ang  =  838
   INTEGER, PARAMETER             :: Anch5Ten  =  839
   INTEGER, PARAMETER             :: Anch5Ang  =  840
   INTEGER, PARAMETER             :: Fair6Ten  =  841
   INTEGER, PARAMETER             :: Fair6Ang  =  842
   INTEGER, PARAMETER             :: Anch6Ten  =  843
   INTEGER, PARAMETER             :: Anch6Ang  =  844
   INTEGER, PARAMETER             :: Fair7Ten  =  845
   INTEGER, PARAMETER             :: Fair7Ang  =  846
   INTEGER, PARAMETER             :: Anch7Ten  =  847
   INTEGER, PARAMETER             :: Anch7Ang  =  848
   INTEGER, PARAMETER             :: Fair8Ten  =  849
   INTEGER, PARAMETER             :: Fair8Ang  =  850
   INTEGER, PARAMETER             :: Anch8Ten  =  851
   INTEGER, PARAMETER             :: Anch8Ang  =  852
   INTEGER, PARAMETER             :: Fair9Ten  =  853
   INTEGER, PARAMETER             :: Fair9Ang  =  854
   INTEGER, PARAMETER             :: Anch9Ten  =  855
   INTEGER, PARAMETER             :: Anch9Ang  =  856


  ! Wave Motions:

   INTEGER, PARAMETER             :: WaveElev  =  857
   INTEGER, PARAMETER             :: Wave1Vxi  =  858
   INTEGER, PARAMETER             :: Wave1Vyi  =  859
   INTEGER, PARAMETER             :: Wave1Vzi  =  860
   INTEGER, PARAMETER             :: Wave1Axi  =  861
   INTEGER, PARAMETER             :: Wave1Ayi  =  862
   INTEGER, PARAMETER             :: Wave1Azi  =  863
   INTEGER, PARAMETER             :: Wave2Vxi  =  864
   INTEGER, PARAMETER             :: Wave2Vyi  =  865
   INTEGER, PARAMETER             :: Wave2Vzi  =  866
   INTEGER, PARAMETER             :: Wave2Axi  =  867
   INTEGER, PARAMETER             :: Wave2Ayi  =  868
   INTEGER, PARAMETER             :: Wave2Azi  =  869
   INTEGER, PARAMETER             :: Wave3Vxi  =  870
   INTEGER, PARAMETER             :: Wave3Vyi  =  871
   INTEGER, PARAMETER             :: Wave3Vzi  =  872
   INTEGER, PARAMETER             :: Wave3Axi  =  873
   INTEGER, PARAMETER             :: Wave3Ayi  =  874
   INTEGER, PARAMETER             :: Wave3Azi  =  875
   INTEGER, PARAMETER             :: Wave4Vxi  =  876
   INTEGER, PARAMETER             :: Wave4Vyi  =  877
   INTEGER, PARAMETER             :: Wave4Vzi  =  878
   INTEGER, PARAMETER             :: Wave4Axi  =  879
   INTEGER, PARAMETER             :: Wave4Ayi  =  880
   INTEGER, PARAMETER             :: Wave4Azi  =  881
   INTEGER, PARAMETER             :: Wave5Vxi  =  882
   INTEGER, PARAMETER             :: Wave5Vyi  =  883
   INTEGER, PARAMETER             :: Wave5Vzi  =  884
   INTEGER, PARAMETER             :: Wave5Axi  =  885
   INTEGER, PARAMETER             :: Wave5Ayi  =  886
   INTEGER, PARAMETER             :: Wave5Azi  =  887
   INTEGER, PARAMETER             :: Wave6Vxi  =  888
   INTEGER, PARAMETER             :: Wave6Vyi  =  889
   INTEGER, PARAMETER             :: Wave6Vzi  =  890
   INTEGER, PARAMETER             :: Wave6Axi  =  891
   INTEGER, PARAMETER             :: Wave6Ayi  =  892
   INTEGER, PARAMETER             :: Wave6Azi  =  893
   INTEGER, PARAMETER             :: Wave7Vxi  =  894
   INTEGER, PARAMETER             :: Wave7Vyi  =  895
   INTEGER, PARAMETER             :: Wave7Vzi  =  896
   INTEGER, PARAMETER             :: Wave7Axi  =  897
   INTEGER, PARAMETER             :: Wave7Ayi  =  898
   INTEGER, PARAMETER             :: Wave7Azi  =  899
   INTEGER, PARAMETER             :: Wave8Vxi  =  900
   INTEGER, PARAMETER             :: Wave8Vyi  =  901
   INTEGER, PARAMETER             :: Wave8Vzi  =  902
   INTEGER, PARAMETER             :: Wave8Axi  =  903
   INTEGER, PARAMETER             :: Wave8Ayi  =  904
   INTEGER, PARAMETER             :: Wave8Azi  =  905
   INTEGER, PARAMETER             :: Wave9Vxi  =  906
   INTEGER, PARAMETER             :: Wave9Vyi  =  907
   INTEGER, PARAMETER             :: Wave9Vzi  =  908
   INTEGER, PARAMETER             :: Wave9Axi  =  909
   INTEGER, PARAMETER             :: Wave9Ayi  =  910
   INTEGER, PARAMETER             :: Wave9Azi  =  911


  ! Internal Degrees of Freedom:

   INTEGER, PARAMETER             :: Q_B1E1    =  912
   INTEGER, PARAMETER             :: Q_B2E1    =  913
   INTEGER, PARAMETER             :: Q_B3E1    =  914
   INTEGER, PARAMETER             :: Q_B1F1    =  915
   INTEGER, PARAMETER             :: Q_B2F1    =  916
   INTEGER, PARAMETER             :: Q_B3F1    =  917
   INTEGER, PARAMETER             :: Q_B1F2    =  918
   INTEGER, PARAMETER             :: Q_B2F2    =  919
   INTEGER, PARAMETER             :: Q_B3F2    =  920
   INTEGER, PARAMETER             :: Q_Teet    =  921
   INTEGER, PARAMETER             :: Q_DrTr    =  922
   INTEGER, PARAMETER             :: Q_GeAz    =  923
   INTEGER, PARAMETER             :: Q_RFrl    =  924
   INTEGER, PARAMETER             :: Q_TFrl    =  925
   INTEGER, PARAMETER             :: Q_Yaw     =  926
   INTEGER, PARAMETER             :: Q_TFA1    =  927
   INTEGER, PARAMETER             :: Q_TSS1    =  928
   INTEGER, PARAMETER             :: Q_TFA2    =  929
   INTEGER, PARAMETER             :: Q_TSS2    =  930
   INTEGER, PARAMETER             :: Q_Sg      =  931
   INTEGER, PARAMETER             :: Q_Sw      =  932
   INTEGER, PARAMETER             :: Q_Hv      =  933
   INTEGER, PARAMETER             :: Q_R       =  934
   INTEGER, PARAMETER             :: Q_P       =  935
   INTEGER, PARAMETER             :: Q_Y       =  936
   INTEGER, PARAMETER             :: QD_B1E1   =  937
   INTEGER, PARAMETER             :: QD_B2E1   =  938
   INTEGER, PARAMETER             :: QD_B3E1   =  939
   INTEGER, PARAMETER             :: QD_B1F1   =  940
   INTEGER, PARAMETER             :: QD_B2F1   =  941
   INTEGER, PARAMETER             :: QD_B3F1   =  942
   INTEGER, PARAMETER             :: QD_B1F2   =  943
   INTEGER, PARAMETER             :: QD_B2F2   =  944
   INTEGER, PARAMETER             :: QD_B3F2   =  945
   INTEGER, PARAMETER             :: QD_Teet   =  946
   INTEGER, PARAMETER             :: QD_DrTr   =  947
   INTEGER, PARAMETER             :: QD_GeAz   =  948
   INTEGER, PARAMETER             :: QD_RFrl   =  949
   INTEGER, PARAMETER             :: QD_TFrl   =  950
   INTEGER, PARAMETER             :: QD_Yaw    =  951
   INTEGER, PARAMETER             :: QD_TFA1   =  952
   INTEGER, PARAMETER             :: QD_TSS1   =  953
   INTEGER, PARAMETER             :: QD_TFA2   =  954
   INTEGER, PARAMETER             :: QD_TSS2   =  955
   INTEGER, PARAMETER             :: QD_Sg     =  956
   INTEGER, PARAMETER             :: QD_Sw     =  957
   INTEGER, PARAMETER             :: QD_Hv     =  958
   INTEGER, PARAMETER             :: QD_R      =  959
   INTEGER, PARAMETER             :: QD_P      =  960
   INTEGER, PARAMETER             :: QD_Y      =  961
   INTEGER, PARAMETER             :: QD2_B1E1  =  962
   INTEGER, PARAMETER             :: QD2_B2E1  =  963
   INTEGER, PARAMETER             :: QD2_B3E1  =  964
   INTEGER, PARAMETER             :: QD2_B1F1  =  965
   INTEGER, PARAMETER             :: QD2_B2F1  =  966
   INTEGER, PARAMETER             :: QD2_B3F1  =  967
   INTEGER, PARAMETER             :: QD2_B1F2  =  968
   INTEGER, PARAMETER             :: QD2_B2F2  =  969
   INTEGER, PARAMETER             :: QD2_B3F2  =  970
   INTEGER, PARAMETER             :: QD2_Teet  =  971
   INTEGER, PARAMETER             :: QD2_DrTr  =  972
   INTEGER, PARAMETER             :: QD2_GeAz  =  973
   INTEGER, PARAMETER             :: QD2_RFrl  =  974
   INTEGER, PARAMETER             :: QD2_TFrl  =  975
   INTEGER, PARAMETER             :: QD2_Yaw   =  976
   INTEGER, PARAMETER             :: QD2_TFA1  =  977
   INTEGER, PARAMETER             :: QD2_TSS1  =  978
   INTEGER, PARAMETER             :: QD2_TFA2  =  979
   INTEGER, PARAMETER             :: QD2_TSS2  =  980
   INTEGER, PARAMETER             :: QD2_Sg    =  981
   INTEGER, PARAMETER             :: QD2_Sw    =  982
   INTEGER, PARAMETER             :: QD2_Hv    =  983
   INTEGER, PARAMETER             :: QD2_R     =  984
   INTEGER, PARAMETER             :: QD2_P     =  985
   INTEGER, PARAMETER             :: QD2_Y     =  986


     ! The maximum number of output channels which can be output by the code.
   INTEGER, PARAMETER             :: MaxOutPts =  986

!End of code generated by Matlab script



   ! Regular variables:

! SEE NOTE ABOVE FOR SIZE (DIMENSION) OF THE VARIABLE BELOW:
REAL(ReKi)                   :: AllOuts  (0:MaxOutPts)                          ! An array holding the value of all of the calculated (not only selected) output channels.
! SEE NOTE ABOVE FOR SIZE (DIMENSION) OF THE PREVIOUS VARIABLE:
REAL(DbKi), ALLOCATABLE      :: TimeData (:)                                    ! Array to contain the time output data for the binary file (first output time and a time [fixed] increment)
REAL(ReKi), ALLOCATABLE      :: AllOutData (:,:)                                ! Array to contain all the output data (time history of all outputs); Index 1 is NumOuts, Index 2 is Time step.
REAL(ReKi), ALLOCATABLE      :: LinAccES (:,:,:)                                ! Total linear acceleration of a point on a   blade (point S) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: LinAccET (:,:)                                  ! Total linear acceleration of a point on the tower (point T) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: LSNodes  (:,:)                                  ! Unstretched arc distance along mooring line from anchor to each node where the line position and tension can be output (meters).
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
REAL(ReKi)                   :: WaveElevxi(1) = (/ 0.0 /)                       ! xi-coordinates for points where the incident wave elevation can be output (meters).
REAL(ReKi)                   :: WaveElevyi(1) = (/ 0.0 /)                       ! yi-coordinates for points where the incident wave elevation can be output (meters).

INTEGER(4)                   :: BldGagNd (9)                                    ! Nodes closest to the blade strain gages.
INTEGER                      :: CurrOutStep                                     ! Time index into the AllOutData arrat
INTEGER(4)                   :: DecFact                                         ! Decimation factor for tabular output.
!JASON: ADD OUTPUTS FOR THE MOORING LINE POSITION AND EFFECTIVE TENSION AT EACH NODE.  USE NAMES SUCH AS: Ln#Nd#Pxi, Ln#Nd#Pyi, Ln#Nd#Pzi, Ln#Nd#Ten WHERE # REPRESENTS THE LINE NUMBER OR NODE NUMBER!!!
INTEGER(4)                   :: LineNodes    = 0                                ! Number of nodes per line where the mooring line position and tension can be output (-).
INTEGER(4)                   :: NBlGages                                        ! Number of blade strain gages.
INTEGER(4)                   :: NTwGages                                        ! Number of tower strain gages.
INTEGER(4)                   :: NumOuts                                         ! Number of parameters in the output list.
INTEGER(IntKi)               :: NOutSteps                                       ! Maximum number of output steps
INTEGER(4)                   :: NWaveElev    = 1                                ! Number of points where the incident wave elevation  can be output (-).
INTEGER(4)                   :: NWaveKin     = 0                                ! Number of points where the incident wave kinematics can be output (-).
INTEGER(4)                   :: TwrGagNd (9)                                    ! Nodes closest to the tower strain gages.
INTEGER(4)                   :: WaveKinNd(9)                                    ! List of tower [not floating] or platform [floating] nodes that have wave kinematics sensors.

INTEGER(B2Ki), PARAMETER     :: FileFmtID_WithTime    = 1                       ! ID for OutputFileFmtID to specify that the time channel should be included in the output file (use if the output can occur at variable times)
INTEGER(B2Ki), PARAMETER     :: FileFmtID_WithoutTime = 2                       ! ID for OutputFileFmtID to specify that the time channel does not need to be included in the output file (used only with constant time-step output)
INTEGER(B2Ki)                :: OutputFileFmtID = FileFmtID_WithoutTime         ! A format specifier for the binary output file format (1=include time channel as packed 32-bit binary; 2=don't include time channel)

LOGICAL                      :: WrEcho
LOGICAL                      :: WrBinOutFile  = .true.                          ! Write a binary output file? (.outb)
LOGICAL                      :: WrTxtOutFile  = .true.                          ! Write a text (formatted) output file? (.out)
LOGICAL                      :: TabDelim                                        ! Flag to cause tab-delimited output.

CHARACTER(20)                :: OutFmt                                          ! Output format for tabular data.
CHARACTER(10)                :: OutList  (MaxOutPts)                            ! List of output parameters.
CHARACTER(1024)              :: FileDesc                                        ! Description of run to include in binary output file

TYPE(OutParmType),ALLOCATABLE:: OutParam (:)                                    ! Names and units of all output parameters.

   ! Let's make these parameters into arrays so we can loop through them
   
INTEGER,  PARAMETER          :: WaveVxi(9) = (/Wave1Vxi,Wave2Vxi,Wave3Vxi,Wave4Vxi,Wave5Vxi,Wave6Vxi,Wave7Vxi,Wave8Vxi,Wave9Vxi/)
INTEGER,  PARAMETER          :: WaveVyi(9) = (/Wave1Vyi,Wave2Vyi,Wave3Vyi,Wave4Vyi,Wave5Vyi,Wave6Vyi,Wave7Vyi,Wave8Vyi,Wave9Vyi/)
INTEGER,  PARAMETER          :: WaveVzi(9) = (/Wave1Vzi,Wave2Vzi,Wave3Vzi,Wave4Vzi,Wave5Vzi,Wave6Vzi,Wave7Vzi,Wave8Vzi,Wave9Vzi/)
INTEGER,  PARAMETER          :: WaveAxi(9) = (/Wave1Axi,Wave2Axi,Wave3Axi,Wave4Axi,Wave5Axi,Wave6Axi,Wave7Axi,Wave8Axi,Wave9Axi/)
INTEGER,  PARAMETER          :: WaveAyi(9) = (/Wave1Ayi,Wave2Ayi,Wave3Ayi,Wave4Ayi,Wave5Ayi,Wave6Ayi,Wave7Ayi,Wave8Ayi,Wave9Ayi/)
INTEGER,  PARAMETER          :: WaveAzi(9) = (/Wave1Azi,Wave2Azi,Wave3Azi,Wave4Azi,Wave5Azi,Wave6Azi,Wave7Azi,Wave8Azi,Wave9Azi/)
INTEGER,  PARAMETER          :: FairTen(9) = (/Fair1Ten,Fair2Ten,Fair3Ten,Fair4Ten,Fair5Ten,Fair6Ten,Fair7Ten,Fair8Ten,Fair9Ten/)
INTEGER,  PARAMETER          :: FairAng(9) = (/Fair1Ang,Fair2Ang,Fair3Ang,Fair4Ang,Fair5Ang,Fair6Ang,Fair7Ang,Fair8Ang,Fair9Ang/)
INTEGER,  PARAMETER          :: AnchTen(9) = (/Anch1Ten,Anch2Ten,Anch3Ten,Anch4Ten,Anch5Ten,Anch6Ten,Anch7Ten,Anch8Ten,Anch9Ten/)
INTEGER,  PARAMETER          :: AnchAng(9) = (/Anch1Ang,Anch2Ang,Anch3Ang,Anch4Ang,Anch5Ang,Anch6Ang,Anch7Ang,Anch8Ang,Anch9Ang/)

INTEGER,  PARAMETER          :: TipDxc( 3)  = (/TipDxc1,  TipDxc2,  TipDxc3/)
INTEGER,  PARAMETER          :: TipDyc( 3)  = (/TipDyc1,  TipDyc2,  TipDyc3/)
INTEGER,  PARAMETER          :: TipDzc( 3)  = (/TipDzc1,  TipDzc2,  TipDzc3/)
INTEGER,  PARAMETER          :: TipDxb( 3)  = (/TipDxb1,  TipDxb2,  TipDxb3/)
INTEGER,  PARAMETER          :: TipDyb( 3)  = (/TipDyb1,  TipDyb2,  TipDyb3/)
INTEGER,  PARAMETER          :: TipALxb(3)  = (/TipALxb1, TipALxb2, TipALxb3/)
INTEGER,  PARAMETER          :: TipALyb(3)  = (/TipALyb1, TipALyb2, TipALyb3/)
INTEGER,  PARAMETER          :: TipALzb(3)  = (/TipALzb1, TipALzb2, TipALzb3/)
INTEGER,  PARAMETER          :: TipRDxb(3)  = (/TipRDxb1, TipRDxb2, TipRDxb3/)
INTEGER,  PARAMETER          :: TipRDyb(3)  = (/TipRDyb1, TipRDyb2, TipRDyb3/)
INTEGER,  PARAMETER          :: TipRDzc(3)  = (/TipRDzc1, TipRDzc2, TipRDzc3/)
INTEGER,  PARAMETER          :: TipClrnc(3) = (/TipClrnc1,TipClrnc2,TipClrnc3/)
INTEGER,  PARAMETER          :: PtchPMzc(3) = (/PtchPMzc1,PtchPMzc2,PtchPMzc3/)

INTEGER,  PARAMETER          :: RootFxc(3) = (/ RootFxc1,RootFxc2,RootFxc3 /)
INTEGER,  PARAMETER          :: RootFyc(3) = (/ RootFyc1,RootFyc2,RootFyc3 /)
INTEGER,  PARAMETER          :: RootFzc(3) = (/ RootFzc1,RootFzc2,RootFzc3 /)
INTEGER,  PARAMETER          :: RootFxb(3) = (/ RootFxb1,RootFxb2,RootFxb3 /)
INTEGER,  PARAMETER          :: RootFyb(3) = (/ RootFyb1,RootFyb2,RootFyb3 /)
INTEGER,  PARAMETER          :: RootMxc(3) = (/ RootMxc1,RootMxc2,RootMxc3 /)
INTEGER,  PARAMETER          :: RootMyc(3) = (/ RootMyc1,RootMyc2,RootMyc3 /)
INTEGER,  PARAMETER          :: RootMzc(3) = (/ RootMzc1,RootMzc2,RootMzc3 /)
INTEGER,  PARAMETER          :: RootMxb(3) = (/ RootMxb1,RootMxb2,RootMxb3 /)
INTEGER,  PARAMETER          :: RootMyb(3) = (/ RootMyb1,RootMyb2,RootMyb3 /)

INTEGER,  PARAMETER          :: SpnALxb(9, 3) = RESHAPE( (/ &
                                    Spn1ALxb1,Spn2ALxb1,Spn3ALxb1,Spn4ALxb1,Spn5ALxb1,Spn6ALxb1,Spn7ALxb1,Spn8ALxb1,Spn9ALxb1, &
                                    Spn1ALxb2,Spn2ALxb2,Spn3ALxb2,Spn4ALxb2,Spn5ALxb2,Spn6ALxb2,Spn7ALxb2,Spn8ALxb2,Spn9ALxb2, &
                                    Spn1ALxb3,Spn2ALxb3,Spn3ALxb3,Spn4ALxb3,Spn5ALxb3,Spn6ALxb3,Spn7ALxb3,Spn8ALxb3,Spn9ALxb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnALyb(9, 3) = RESHAPE( (/ &
                                    Spn1ALyb1,Spn2ALyb1,Spn3ALyb1,Spn4ALyb1,Spn5ALyb1,Spn6ALyb1,Spn7ALyb1,Spn8ALyb1,Spn9ALyb1, &
                                    Spn1ALyb2,Spn2ALyb2,Spn3ALyb2,Spn4ALyb2,Spn5ALyb2,Spn6ALyb2,Spn7ALyb2,Spn8ALyb2,Spn9ALyb2, &
                                    Spn1ALyb3,Spn2ALyb3,Spn3ALyb3,Spn4ALyb3,Spn5ALyb3,Spn6ALyb3,Spn7ALyb3,Spn8ALyb3,Spn9ALyb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnALzb(9, 3) = RESHAPE( (/ &
                                    Spn1ALzb1,Spn2ALzb1,Spn3ALzb1,Spn4ALzb1,Spn5ALzb1,Spn6ALzb1,Spn7ALzb1,Spn8ALzb1,Spn9ALzb1, &
                                    Spn1ALzb2,Spn2ALzb2,Spn3ALzb2,Spn4ALzb2,Spn5ALzb2,Spn6ALzb2,Spn7ALzb2,Spn8ALzb2,Spn9ALzb2, &
                                    Spn1ALzb3,Spn2ALzb3,Spn3ALzb3,Spn4ALzb3,Spn5ALzb3,Spn6ALzb3,Spn7ALzb3,Spn8ALzb3,Spn9ALzb3  &
                                /), (/9, 3/) )

INTEGER,  PARAMETER          :: SpnFLxb(9,3) = RESHAPE( (/ &
                                    Spn1FLxb1,Spn2FLxb1,Spn3FLxb1,Spn4FLxb1,Spn5FLxb1,Spn6FLxb1,Spn7FLxb1,Spn8FLxb1,Spn9FLxb1, &
                                    Spn1FLxb2,Spn2FLxb2,Spn3FLxb2,Spn4FLxb2,Spn5FLxb2,Spn6FLxb2,Spn7FLxb2,Spn8FLxb2,Spn9FLxb2, &
                                    Spn1FLxb3,Spn2FLxb3,Spn3FLxb3,Spn4FLxb3,Spn5FLxb3,Spn6FLxb3,Spn7FLxb3,Spn8FLxb3,Spn9FLxb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnFLyb(9,3) = RESHAPE( (/ &
                                    Spn1FLyb1,Spn2FLyb1,Spn3FLyb1,Spn4FLyb1,Spn5FLyb1,Spn6FLyb1,Spn7FLyb1,Spn8FLyb1,Spn9FLyb1, &
                                    Spn1FLyb2,Spn2FLyb2,Spn3FLyb2,Spn4FLyb2,Spn5FLyb2,Spn6FLyb2,Spn7FLyb2,Spn8FLyb2,Spn9FLyb2, &
                                    Spn1FLyb3,Spn2FLyb3,Spn3FLyb3,Spn4FLyb3,Spn5FLyb3,Spn6FLyb3,Spn7FLyb3,Spn8FLyb3,Spn9FLyb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnFLzb(9,3) = RESHAPE( (/ &
                                    Spn1FLzb1,Spn2FLzb1,Spn3FLzb1,Spn4FLzb1,Spn5FLzb1,Spn6FLzb1,Spn7FLzb1,Spn8FLzb1,Spn9FLzb1, &
                                    Spn1FLzb2,Spn2FLzb2,Spn3FLzb2,Spn4FLzb2,Spn5FLzb2,Spn6FLzb2,Spn7FLzb2,Spn8FLzb2,Spn9FLzb2, &
                                    Spn1FLzb3,Spn2FLzb3,Spn3FLzb3,Spn4FLzb3,Spn5FLzb3,Spn6FLzb3,Spn7FLzb3,Spn8FLzb3,Spn9FLzb3  &
                                /), (/9, 3/) )
                                
INTEGER,  PARAMETER          :: SpnMLxb(9,3) = RESHAPE( (/ &
                                    Spn1MLxb1,Spn2MLxb1,Spn3MLxb1,Spn4MLxb1,Spn5MLxb1,Spn6MLxb1,Spn7MLxb1,Spn8MLxb1,Spn9MLxb1, &
                                    Spn1MLxb2,Spn2MLxb2,Spn3MLxb2,Spn4MLxb2,Spn5MLxb2,Spn6MLxb2,Spn7MLxb2,Spn8MLxb2,Spn9MLxb2, &
                                    Spn1MLxb3,Spn2MLxb3,Spn3MLxb3,Spn4MLxb3,Spn5MLxb3,Spn6MLxb3,Spn7MLxb3,Spn8MLxb3,Spn9MLxb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnMLyb(9,3) = RESHAPE( (/ &
                                    Spn1MLyb1,Spn2MLyb1,Spn3MLyb1,Spn4MLyb1,Spn5MLyb1,Spn6MLyb1,Spn7MLyb1,Spn8MLyb1,Spn9MLyb1, &
                                    Spn1MLyb2,Spn2MLyb2,Spn3MLyb2,Spn4MLyb2,Spn5MLyb2,Spn6MLyb2,Spn7MLyb2,Spn8MLyb2,Spn9MLyb2, &
                                    Spn1MLyb3,Spn2MLyb3,Spn3MLyb3,Spn4MLyb3,Spn5MLyb3,Spn6MLyb3,Spn7MLyb3,Spn8MLyb3,Spn9MLyb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnMLzb(9,3) = RESHAPE( (/ &
                                    Spn1MLzb1,Spn2MLzb1,Spn3MLzb1,Spn4MLzb1,Spn5MLzb1,Spn6MLzb1,Spn7MLzb1,Spn8MLzb1,Spn9MLzb1, &
                                    Spn1MLzb2,Spn2MLzb2,Spn3MLzb2,Spn4MLzb2,Spn5MLzb2,Spn6MLzb2,Spn7MLzb2,Spn8MLzb2,Spn9MLzb2, &
                                    Spn1MLzb3,Spn2MLzb3,Spn3MLzb3,Spn4MLzb3,Spn5MLzb3,Spn6MLzb3,Spn7MLzb3,Spn8MLzb3,Spn9MLzb3  &
                                /), (/9, 3/) )
                                
INTEGER,  PARAMETER          :: SpnTDxb(9,3) = RESHAPE( (/ &
                                    Spn1TDxb1,Spn2TDxb1,Spn3TDxb1,Spn4TDxb1,Spn5TDxb1,Spn6TDxb1,Spn7TDxb1,Spn8TDxb1,Spn9TDxb1, &
                                    Spn1TDxb2,Spn2TDxb2,Spn3TDxb2,Spn4TDxb2,Spn5TDxb2,Spn6TDxb2,Spn7TDxb2,Spn8TDxb2,Spn9TDxb2, &
                                    Spn1TDxb3,Spn2TDxb3,Spn3TDxb3,Spn4TDxb3,Spn5TDxb3,Spn6TDxb3,Spn7TDxb3,Spn8TDxb3,Spn9TDxb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnTDyb(9,3) = RESHAPE( (/ &
                                    Spn1TDyb1,Spn2TDyb1,Spn3TDyb1,Spn4TDyb1,Spn5TDyb1,Spn6TDyb1,Spn7TDyb1,Spn8TDyb1,Spn9TDyb1, &
                                    Spn1TDyb2,Spn2TDyb2,Spn3TDyb2,Spn4TDyb2,Spn5TDyb2,Spn6TDyb2,Spn7TDyb2,Spn8TDyb2,Spn9TDyb2, &
                                    Spn1TDyb3,Spn2TDyb3,Spn3TDyb3,Spn4TDyb3,Spn5TDyb3,Spn6TDyb3,Spn7TDyb3,Spn8TDyb3,Spn9TDyb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnTDzb(9,3) = RESHAPE( (/ &
                                    Spn1TDzb1,Spn2TDzb1,Spn3TDzb1,Spn4TDzb1,Spn5TDzb1,Spn6TDzb1,Spn7TDzb1,Spn8TDzb1,Spn9TDzb1, &
                                    Spn1TDzb2,Spn2TDzb2,Spn3TDzb2,Spn4TDzb2,Spn5TDzb2,Spn6TDzb2,Spn7TDzb2,Spn8TDzb2,Spn9TDzb2, &
                                    Spn1TDzb3,Spn2TDzb3,Spn3TDzb3,Spn4TDzb3,Spn5TDzb3,Spn6TDzb3,Spn7TDzb3,Spn8TDzb3,Spn9TDzb3  &
                                /), (/9, 3/) )

INTEGER,  PARAMETER          :: SpnRDxb(9,3) = RESHAPE( (/ &
                                    Spn1RDxb1,Spn2RDxb1,Spn3RDxb1,Spn4RDxb1,Spn5RDxb1,Spn6RDxb1,Spn7RDxb1,Spn8RDxb1,Spn9RDxb1, &
                                    Spn1RDxb2,Spn2RDxb2,Spn3RDxb2,Spn4RDxb2,Spn5RDxb2,Spn6RDxb2,Spn7RDxb2,Spn8RDxb2,Spn9RDxb2, &
                                    Spn1RDxb3,Spn2RDxb3,Spn3RDxb3,Spn4RDxb3,Spn5RDxb3,Spn6RDxb3,Spn7RDxb3,Spn8RDxb3,Spn9RDxb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnRDyb(9,3) = RESHAPE( (/ &
                                    Spn1RDyb1,Spn2RDyb1,Spn3RDyb1,Spn4RDyb1,Spn5RDyb1,Spn6RDyb1,Spn7RDyb1,Spn8RDyb1,Spn9RDyb1, &
                                    Spn1RDyb2,Spn2RDyb2,Spn3RDyb2,Spn4RDyb2,Spn5RDyb2,Spn6RDyb2,Spn7RDyb2,Spn8RDyb2,Spn9RDyb2, &
                                    Spn1RDyb3,Spn2RDyb3,Spn3RDyb3,Spn4RDyb3,Spn5RDyb3,Spn6RDyb3,Spn7RDyb3,Spn8RDyb3,Spn9RDyb3  &
                                /), (/9, 3/) )
INTEGER,  PARAMETER          :: SpnRDzb(9,3) = RESHAPE( (/ &
                                    Spn1RDzb1,Spn2RDzb1,Spn3RDzb1,Spn4RDzb1,Spn5RDzb1,Spn6RDzb1,Spn7RDzb1,Spn8RDzb1,Spn9RDzb1, &
                                    Spn1RDzb2,Spn2RDzb2,Spn3RDzb2,Spn4RDzb2,Spn5RDzb2,Spn6RDzb2,Spn7RDzb2,Spn8RDzb2,Spn9RDzb2, &
                                    Spn1RDzb3,Spn2RDzb3,Spn3RDzb3,Spn4RDzb3,Spn5RDzb3,Spn6RDzb3,Spn7RDzb3,Spn8RDzb3,Spn9RDzb3  &
                                /), (/9, 3/) )

                                
INTEGER,  PARAMETER          :: TwHtALxt(9) = (/ &
                                    TwHt1ALxt,TwHt2ALxt,TwHt3ALxt,TwHt4ALxt,TwHt5ALxt,TwHt6ALxt,TwHt7ALxt,TwHt8ALxt,TwHt9ALxt /)
INTEGER,  PARAMETER          :: TwHtALyt(9) = (/ &
                                    TwHt1ALyt,TwHt2ALyt,TwHt3ALyt,TwHt4ALyt,TwHt5ALyt,TwHt6ALyt,TwHt7ALyt,TwHt8ALyt,TwHt9ALyt /)
INTEGER,  PARAMETER          :: TwHtALzt(9) = (/ &
                                    TwHt1ALzt,TwHt2ALzt,TwHt3ALzt,TwHt4ALzt,TwHt5ALzt,TwHt6ALzt,TwHt7ALzt,TwHt8ALzt,TwHt9ALzt /)

INTEGER,  PARAMETER          :: TwHtMLxt(9) = (/ &
                                    TwHt1MLxt,TwHt2MLxt,TwHt3MLxt,TwHt4MLxt,TwHt5MLxt,TwHt6MLxt,TwHt7MLxt,TwHt8MLxt,TwHt9MLxt /)
INTEGER,  PARAMETER          :: TwHtMLyt(9) = (/ &
                                    TwHt1MLyt,TwHt2MLyt,TwHt3MLyt,TwHt4MLyt,TwHt5MLyt,TwHt6MLyt,TwHt7MLyt,TwHt8MLyt,TwHt9MLyt /)
INTEGER,  PARAMETER          :: TwHtMLzt(9) = (/ &
                                    TwHt1MLzt,TwHt2MLzt,TwHt3MLzt,TwHt4MLzt,TwHt5MLzt,TwHt6MLzt,TwHt7MLzt,TwHt8MLzt,TwHt9MLzt /)

INTEGER,  PARAMETER          :: TwHtFLxt(9) = (/ &
                                    TwHt1FLxt,TwHt2FLxt,TwHt3FLxt,TwHt4FLxt,TwHt5FLxt,TwHt6FLxt,TwHt7FLxt,TwHt8FLxt,TwHt9FLxt /)
INTEGER,  PARAMETER          :: TwHtFLyt(9) = (/ &
                                    TwHt1FLyt,TwHt2FLyt,TwHt3FLyt,TwHt4FLyt,TwHt5FLyt,TwHt6FLyt,TwHt7FLyt,TwHt8FLyt,TwHt9FLyt /)
INTEGER,  PARAMETER          :: TwHtFLzt(9) = (/ &
                                    TwHt1FLzt,TwHt2FLzt,TwHt3FLzt,TwHt4FLzt,TwHt5FLzt,TwHt6FLzt,TwHt7FLzt,TwHt8FLzt,TwHt9FLzt /)

INTEGER,  PARAMETER          :: TwHtTDxt(9) = (/ &
                                    TwHt1TDxt,TwHt2TDxt,TwHt3TDxt,TwHt4TDxt,TwHt5TDxt,TwHt6TDxt,TwHt7TDxt,TwHt8TDxt,TwHt9TDxt /)
INTEGER,  PARAMETER          :: TwHtTDyt(9) = (/ &
                                    TwHt1TDyt,TwHt2TDyt,TwHt3TDyt,TwHt4TDyt,TwHt5TDyt,TwHt6TDyt,TwHt7TDyt,TwHt8TDyt,TwHt9TDyt /)
INTEGER,  PARAMETER          :: TwHtTDzt(9) = (/ &
                                    TwHt1TDzt,TwHt2TDzt,TwHt3TDzt,TwHt4TDzt,TwHt5TDzt,TwHt6TDzt,TwHt7TDzt,TwHt8TDzt,TwHt9TDzt /)

INTEGER,  PARAMETER          :: TwHtRDxt(9) = (/ &
                                    TwHt1RDxt,TwHt2RDxt,TwHt3RDxt,TwHt4RDxt,TwHt5RDxt,TwHt6RDxt,TwHt7RDxt,TwHt8RDxt,TwHt9RDxt /)
INTEGER,  PARAMETER          :: TwHtRDyt(9) = (/ &
                                    TwHt1RDyt,TwHt2RDyt,TwHt3RDyt,TwHt4RDyt,TwHt5RDyt,TwHt6RDyt,TwHt7RDyt,TwHt8RDyt,TwHt9RDyt /)
INTEGER,  PARAMETER          :: TwHtRDzt(9) = (/ &
                                    TwHt1RDzt,TwHt2RDzt,TwHt3RDzt,TwHt4RDzt,TwHt5RDzt,TwHt6RDzt,TwHt7RDzt,TwHt8RDzt,TwHt9RDzt /)

INTEGER,  PARAMETER          :: TwHtTPxi(9) = (/ &
                                    TwHt1TPxi,TwHt2TPxi,TwHt3TPxi,TwHt4TPxi,TwHt5TPxi,TwHt6TPxi,TwHt7TPxi,TwHt8TPxi,TwHt9TPxi /)
INTEGER,  PARAMETER          :: TwHtTPyi(9) = (/ &
                                    TwHt1TPyi,TwHt2TPyi,TwHt3TPyi,TwHt4TPyi,TwHt5TPyi,TwHt6TPyi,TwHt7TPyi,TwHt8TPyi,TwHt9TPyi /)
INTEGER,  PARAMETER          :: TwHtTPzi(9) = (/ &
                                    TwHt1TPzi,TwHt2TPzi,TwHt3TPzi,TwHt4TPzi,TwHt5TPzi,TwHt6TPzi,TwHt7TPzi,TwHt8TPzi,TwHt9TPzi /)

INTEGER,  PARAMETER          :: TwHtRPxi(9) = (/ &
                                    TwHt1RPxi,TwHt2RPxi,TwHt3RPxi,TwHt4RPxi,TwHt5RPxi,TwHt6RPxi,TwHt7RPxi,TwHt8RPxi,TwHt9RPxi /)
INTEGER,  PARAMETER          :: TwHtRPyi(9) = (/ &
                                    TwHt1RPyi,TwHt2RPyi,TwHt3RPyi,TwHt4RPyi,TwHt5RPyi,TwHt6RPyi,TwHt7RPyi,TwHt8RPyi,TwHt9RPyi /)
INTEGER,  PARAMETER          :: TwHtRPzi(9) = (/ &
                                    TwHt1RPzi,TwHt2RPzi,TwHt3RPzi,TwHt4RPzi,TwHt5RPzi,TwHt6RPzi,TwHt7RPzi,TwHt8RPzi,TwHt9RPzi /)
   
END MODULE Output
!=======================================================================
MODULE Platform


   ! This MODULE stores input variables for platform loading.


USE                             Precision

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
REAL(ReKi)                   :: MaxLRadAnch  = 0.0                              ! Maximum value of input array LRadAnch.

REAL(ReKi)                   :: PtfmAM (6,6) = 0.0                              ! Platform added mass matrix.
REAL(ReKi)                   :: PtfmCD                                          ! Effective platform normalized hydrodynamic viscous drag coefficient in calculation of viscous drag term from Morison's equation.
REAL(ReKi)                   :: PtfmDiam                                        ! Effective platform diameter in calculation of viscous drag term from Morison's equation.
REAL(ReKi)                   :: PtfmDraft                                       ! Effective platform draft    in calculation of viscous drag term from Morison's equation.

REAL(ReKi)                   :: PtfmFt   (6) = 0.0                              ! The surge/xi (1), sway/yi (2), and heave/zi (3)-components of the portion of the platform force at the platform reference (point Z) and the roll/xi (4), pitch/yi (5), and yaw/zi (6)-components of the portion of the platform moment acting at the platform (body X) / platform reference (point Z) associated with everything but the QD2T()'s.

REAL(ReKi)                   :: PtfmVol0                                        ! Displaced volume of water when the platform is in its undisplaced position.
REAL(ReKi)                   :: RdtnDT                                          ! Time step for wave radiation kernel calculations.
REAL(ReKi)                   :: RdtnTMax     = 0.0                              ! Analysis time for wave radiation kernel calculations.

INTEGER(4)                   :: LineMod                                         ! Mooring line model switch.
INTEGER(4)                   :: NumLines     = 0                                ! Number of mooring lines.

INTEGER(4)                   :: PtfmLdMod    = 0                                ! Platform loading model switch. (Initialized to zero b/c not all models read in PtfmFile)
INTEGER(4)                   :: PtfmNodes                                       ! Number of platform nodes used in calculation of viscous drag term from Morison's equation.

CHARACTER(1024)              :: WAMITFile                                       ! Root name of WAMIT output files containing the linear, nondimensionalized, hydrostatic restoring matrix (.hst extension), frequency-dependent hydrodynamic added mass matrix and damping matrix (.1 extension), and frequency- and direction-dependent wave excitation force vector per unit wave amplitude (.3 extension).


END MODULE Platform
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
REAL(ReKi), ALLOCATABLE      :: AngAccEFt(:,:)                                  ! Portion of the angular acceleration of tower element J                                               (body F) in the inertia frame (body E for earth) associated with everything but the QD2T()'s.

REAL(ReKi), ALLOCATABLE      :: AngPosEF (:,:)                                  ! Angular position of the current point on the tower                                (body F) in the inertial frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: AngPosXF (:,:)                                  ! Angular position of the current point on the tower                                (body F) in the platform      (body X          ).
REAL(ReKi), ALLOCATABLE      :: AngPosHM (:,:,:)                                ! Angular position of eleMent J of blade K                                          (body M) in the hub           (body H          ).
REAL(ReKi)                   :: AngPosXB (3)                                    ! Angular position of the base plate                                                (body B) in the platform      (body X          ).
REAL(ReKi)                   :: AngVelEB (3)                                    ! Angular velocity of the base plate                                                (body B) in the inertia frame (body E for earth).
REAL(ReKi), ALLOCATABLE      :: AngVelEF  (:,:)                                 ! Angular velocity of the current point on the tower                                (body F) in the inertia frame (body E for earth).
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
REAL(ReKi), ALLOCATABLE      :: LinVelET  (:,:)                                 ! Linear velocity of current point on the tower         (point T) in the inertia frame.
REAL(ReKi), ALLOCATABLE      :: LinVelESm2 (:)                                  ! The m2-component (closest to tip) of LinVelES.
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
REAL(ReKi), ALLOCATABLE      :: rT       (:,:)                                  ! Position vector from inertial frame origin to the current node (point T(HNodes(J)).
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
!bjj: these variables should be initialized in an intialization subroutine (for Simulink)


REAL(ReKi)                   :: DT                                              ! Integration time step.
REAL(ReKi)                   :: DT24                                            ! DT/24.
REAL(ReKi)                   :: TMax                                            ! Total run time.
REAL(ReKi)                   :: ZTime    = 0.0                                  ! Current simulation time.

REAL(4)                      :: UsrTime1                                        ! User CPU time for simulation initialization.

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
REAL(ReKi), ALLOCATABLE      :: CAT       (:)                                   ! Interpolated, normalized hydrodynamic added mass   coefficient in Morison's equation.
REAL(ReKi), ALLOCATABLE      :: CDT       (:)                                   ! Interpolated, normalized hydrodynamic viscous drag coefficient in Morison's equation.
REAL(ReKi), ALLOCATABLE      :: cgOffTFA  (:)                                    ! Interpolated tower fore-aft mass cg offset.
REAL(ReKi), ALLOCATABLE      :: cgOffTSS  (:)                                    ! Interpolated tower side-to-side mass cg offset.
REAL(ReKi)                   :: CTFA      (2,2)                                  ! Generalized damping of tower in fore-aft direction.
REAL(ReKi)                   :: CTSS      (2,2)                                  ! Generalized damping of tower in side-to-side direction.
REAL(ReKi), ALLOCATABLE      :: DHNodes   (:)                                    ! Length of variable-length tower elements
REAL(ReKi), ALLOCATABLE      :: DiamT     (:)                                   ! Interpolated tower diameter in Morison's equation.
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
REAL(ReKi)                   :: TwrAM     (6,6) = 0.0                           ! Added mass matrix of the current tower element per unit length.
REAL(ReKi)                   :: TwrCA    = 0.0                                  ! Normalized hydrodynamic added mass   coefficient in Morison's equation.
REAL(ReKi)                   :: TwrCD    = 0.0                                  ! Normalized hydrodynamic viscous drag coefficient in Morison's equation.
REAL(ReKi)                   :: TwrDiam  = 0.0                                  ! Tower diameter in Morison's equation.
REAL(ReKi)                   :: TwrFADmp  (2)                                    ! Tower fore-aft structural damping ratios.
REAL(ReKi), ALLOCATABLE      :: TwrFASF   (:,:,:)                                ! Tower fore-aft shape functions.
REAL(ReKi)                   :: TwrFlexL                                         ! Height / length of the flexible portion of the tower.
REAL(ReKi)                   :: TwrFt     (6)   = 0.0                            ! The surge/xi (1), sway/yi (2), and heave/zi (3)-components of the portion of the tower force at the current tower element (point T) and the roll/xi (4), pitch/yi (5), and yaw/zi (6)-components of the portion of the tower moment acting at the current tower element (body F) / (point T) per unit length associated with everything but the QD2T()'s.

REAL(ReKi)                   :: TwrSSDmp  (2)                                    ! Tower side-to-side structural damping ratios.
REAL(ReKi), ALLOCATABLE      :: TwrSSSF   (:,:,:)                                ! Tower side-to-side shape functions.
REAL(ReKi), ALLOCATABLE      :: TwSScgOf  (:)                                    ! Tower fore-aft (along the yt-axis) mass cg offset for a given input station.
REAL(ReKi), ALLOCATABLE      :: TwSSIner  (:)                                    ! Tower side-to-side (about xt-axis) mass inertia per unit length for a given input station.
REAL(ReKi), ALLOCATABLE      :: TwSSStif  (:)                                    ! Tower side-to-side stiffness for a given input station.

INTEGER(4)                   :: NTwInpSt                                        ! Number of tower input stations.
INTEGER(4)                   :: TTopNode                                        ! Index of the additional node located at the tower-top = TwrNodes + 1
INTEGER(4)                   :: TwrLdMod = 0                                    ! Tower loading model switch.
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

LOGICAL,    ALLOCATABLE      :: BegPitMan(:)                                    ! .TRUE. before the override pitch manuever has begun (begin pitch manuever).
LOGICAL                      :: GenTiStp                                        ! Stop generator based upon T: time or F: generator power = 0.
LOGICAL                      :: GenTiStr                                        ! Start generator based upon T: time or F: generator speed.


END MODULE TurbCont
!=======================================================================
