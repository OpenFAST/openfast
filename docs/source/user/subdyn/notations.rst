
.. contents::
   :depth: 3
..

**NOTICE**

This report was prepared as an account of work sponsored by an agency of
the United States government. Neither the United States government nor
any agency thereof, nor any of their employees, makes any warranty,
express or implied, or assumes any legal liability or responsibility for
the accuracy, completeness, or usefulness of any information, apparatus,
product, or process disclosed, or represents that its use would not
infringe privately owned rights. Reference herein to any specific
commercial product, process, or service by trade name, trademark,
manufacturer, or otherwise does not necessarily constitute or imply its
endorsement, recommendation, or favoring by the United States government
or any agency thereof. The views and opinions of authors expressed
herein do not necessarily state or reflect those of the United States
government or any agency thereof.

   | This report is available at no cost from the National Renewable
     Energy
   | Laboratory (NREL) at www.nrel.gov/publications.

   Available electronically at http://www.osti.gov/scitech

   | Available for a processing fee to U.S. Department of Energy
   | and its contractors, in paper, from:

   | U.S. Department of Energy
   | Office of Scientific and Technical Information

   | P.O. Box 62
   | Oak Ridge, TN 37831-0062
   | phone: 865.576.8401
   | fax: 865.576.5728
   | email:
     `mailto:reports@adonis.osti.gov <mailto:reports@adonis.osti.gov>`__

   Available for sale to the public, in paper, from:

   | U.S. Department of Commerce
   | National Technical Information Service
   | 5285 Port Royal Road
   | Springfield, VA 22161
   | phone: 800.553.6847
   | fax: 703.605.6900
   | email: orders@ntis.fedworld.gov
   | online ordering: http://www.ntis.gov/help/ordermethods.aspx

   *Cover Photos: (left to right) photo by Pat Corkery, NREL 16416,
   photo from SunEdison, NREL 17423, photo by Pat Corkery, NREL 16560,
   photo by Dennis Schroeder, NREL 17613, photo by Dean Armstrong, NREL
   17436, photo by Pat Corkery, NREL 17721.*

|Recycle Logo| Printed on paper containing at least 50% wastepaper,
including 10% post consumer waste.

Preface

This manual offers both a quick reference guide and a more in-depth
theory guide for the SubDyn software program. It is intended to be used
by the general user in combination with the FAST and HydroDyn manuals.
The manual will be updated as new releases are issued and as needed to
provide further information on advancements or modifications to the
software.

**
**

Acknowledgments

The authors would like to acknowledge Huimin Song, a post-doctoral
researcher at the National Renewable Energy Laboratory, for her
contribution to the development of SubDyn. Amy Robertson and Fabian
Wendt provided helpful feedback on the document, and Braulio Barahona
produced some of the graphics. We also are grateful to the U.S.
Department of Energy Wind and Water Power Program for supporting the
development of this software.

Nomenclature

:math:`\zeta` Damping ratio (percent of critical) or diagonal damping
matrix

:math:`\omega` Eigenfrequency

:math:`\rho` Material density

:math:`\Phi_{L}` *L*\ ×\ *L* matrix representing the C-B eigenmodes

:math:`\Phi_{m}` *L*\ ×\ *m* matrix representing the retained C-B
eigenmodes

:math:`\Phi_{R}` *L*\ ×\ *R* matrix representing the interior nodal
displacements for a static, rigid-body motion of the boundary nodes

:math:`{\bar{\Phi}}_{R}` Matrix representing the interior nodal
displacements for a static, rigid-body motion of the boundary nodes,
after the restrained degrees of freedoms are removed

:math:`\Omega_{L}` *L*\ ×\ *L* diagonal matrix representing the
eigenfrequencies of all C-B eigenmodes

:math:`\Omega_{m}` *m*\ ×\ *m* diagonal matrix representing the
eigenfrequencies of the retained C-B eigenmodes

ν Poisson’s ratio

*A* State matrix in state-space state equation

AF Apparent fixity

:math:`A_{\text{sx}}` Member shear area along local *x*-axis

:math:`A_{\text{sy}}` Member shear area along local *y*-axis

:math:`A_{z}` Member cross-sectional area, normal to local *z*-axis

*B* Input matrix in state-space state equation

C-B Craig Bampton

:math:`C` System damping matrix

:math:`C_{1}` State matrix in state-space output equation for ElastoDyn

:math:`C_{2}` State matrix in state-space output equation for HydroDyn

:math:`C_{\text{LL}}` Damping matrix partition referred to the *L*
interior degrees of freedom (*L*\ ×\ *L*)

:math:`C_{\text{LR}}` Damping matrix partition referred to the *L* and
*R* degrees of freedom (*L*\ ×\ *R*)

:math:`C_{\text{RL}}` Damping matrix partition referred to the *R* and
*L* boundary degrees of freedom (*R*\ ×\ *L*)

:math:`C_{\text{RR}}` Damping matrix partition referred to the *R*
boundary degrees of freedom (*R*\ ×\ *R*)

:math:`C_{33}^{\left( \text{HD} \right)}` HydroDyn’s linear damping
matrix heave element

:math:`\text{dt}` Time step in the integration algorithm

:math:`dt_{\max}` Maximum time-step recommended

*D\ i* Member inner diameter

*D\ o* Member outer diameter

:math:`D_{1}` Input transmission matrix in state-space output equation
for ElastoDyn

:math:`D_{2}` Input transmission matrix in state-space output equation
for HydroDyn

DOF Degree of freedom

:math:`\left\lbrack D_{C} \right\rbrack` Member direction cosine 3x3
matrix, transforming local coordinates to global coordinates

*E* Young’s modulus

*f* Frequency

:math:`f_{\max}` Maximum frequency

*F* Vector of external forces and moments

:math:`F_{G}` Array of self-weight forces and moments

:math:`F_{L}` Array of external (hydrodynamic) forces and moments at the
interior nodes

:math:`{\widetilde{F}}_{L}` Interior forces and moments (*L*\ ×1) array
accounting for all interior modes
=\ :math:`\Phi_{L}^{T}\left( F_{L} + F_{\text{Lg}} \right)`

:math:`F_{\text{Lg}}` Array of self-weight gravity forces and moments at
the interior nodes

:math:`{\widetilde{F}}_{m}` Interior forces and moments (*L*\ ×1) array
accounting for only the retained modes
=\ :math:`\Phi_{m}^{T}\left( F_{L} + F_{\text{Lg}} \right)`

:math:`F_{R}` Array of external (hydrodynamic) forces and moments at the
boundary nodes

:math:`{\bar{F}}_{R}` Array of external (hydrodynamic) forces and
moments at the interface nodes

:math:`F_{\text{Re}a\text{ct}}` Substructure base reaction forces and
moments, applied to the substructure by the restraint

:math:`F_{\text{Rg}}` Array of self-weight forces and moments at the
boundary nodes

:math:`{\bar{F}}_{\text{Rg}}` Array of self-weight forces and moments at
the interface nodes

:math:`F_{\text{TP}}` TP reaction forces and moments, applied to the
substructure by the remainder of the turbine system

:math:`{\widetilde{F}}_{\text{TP}}` TP forces and moments after C-B
reduction, =
:math:`F_{\text{TP}} + T_{I}^{T}{\bar{F}}_{\text{Rg}} + T_{I}^{T}{\bar{\Phi}}_{R}^{T}\left( F_{L} + F_{\text{Lg}} \right)`

*F\ X* Substructure base reaction force along global *X*-axis; or

forcing vector in state-space state equation

*F\ Y* Substructure base reaction force along global *Y*-axis

*F\ Y1* Forcing vector in state-space output equation for ElastoDyn

*F\ Y2* Forcing vector in state-space output equation for HydroDyn

*F\ Z* Substructure base reaction force along global *Z*-axis

:math:`F_{I}^{e}` Element inertia forces

:math:`F_{S}^{e}` Element static forces

FEA Finite-element analysis

FEM Finite-element method

:math:`g` Gravity, unsigned magnitude

*G* Shear modulus

:math:`i` Member or element index

:math:`{\widehat{i}}_{e}` Unit vector along the element local *x*-axis

:math:`I` Identity matrix

:math:`\widehat{I}` Unit vector along the global *X*-axis

:math:`{\widehat{j}}_{e}` Unit vector along the element local *y*-axis

:math:`J` Generic second area moment of inertia

:math:`\widehat{J}` Unit vector along the global *Y*-axis

:math:`J_{x}` Second area moment of inertia about the local, principal
*x*-axis

:math:`J_{\text{xx}}` Second area moment of inertia about the local
*x*-axis

:math:`k` Element stiffness matrix (12x12) in global coordinate system

:math:`k_{e}` Element stiffness matrix (12x12)

:math:`{\widehat{k}}_{e}` Unit vector along the element local *z*-axis

:math:`k_{\text{ax}}` Shear area factor along local *x*-axis

:math:`k_{\text{ay}}` Shear area factor along local *y*-axis

:math:`K` System stiffness matrix

:math:`\widehat{K}` Unit vector along the global *Z*-axis

:math:`K_{\text{BB}}` Matrix partition after C-B system reduction =
:math:`K_{\text{RR}} + K_{\text{RL}}\Phi_{R}`

:math:`{\bar{K}}_{\text{BB}}` :math:`K_{\text{BB}}` after removal of
constrained DOF rows and columns

:math:`{\widetilde{K}}_{\text{BB}}` Substructure equivalent stiffness
matrix referred to the TP reference point, =
:math:`T_{I}^{T}{\bar{K}}_{\text{BB}}T_{I}`

:math:`K_{\text{LL}}` Stiffness matrix partition referred to the *L*
interior DOFs (*L*\ ×\ *L*)

:math:`K_{\text{LR}}` Stiffness matrix partition referred to the *L* and
*R* DOFs (*L*\ ×\ *R*)

:math:`K_{\text{RL}}` Stiffness matrix partition referred to the *R* and
*L* boundary DOFs (*R*\ ×\ *L*)

:math:`K_{\text{RR}}` Stiffness matrix partition referred to the R
boundary DOFs (*R*\ ×\ *R*)

:math:`K_{\text{sx}}` Shear correction factor along local *x*-axis

:math:`K_{\text{sy}}` Shear correction factor along local *y*-axis

:math:`K_{33}^{\left( \text{SD} \right)}` Substructure equivalent heave
stiffness

*L* Interior nodes’ DOFs

LFEB Linear frame finite-element beam model

:math:`L_{e}` Length of element

:math:`L_{\text{exy}}` Length of element projection in the global *XY*
plane

:math:`m_{e}` Element stiffness matrix (12x12)

:math:`m` Element stiffness matrix (12x12) in global coordinate system;
or

number of retained C-B modes

:math:`M` System mass matrix

MSL Mean sea level

:math:`M_{\text{BB}}` Matrix partition after C-B system reduction =
:math:`M_{\text{RR}} + M_{\text{RL}}\Phi_{R} + \Phi_{R}^{T}M_{\text{LR}} + \Phi_{R}^{T}M_{\text{LL}}\Phi_{R}`

:math:`{\bar{M}}_{\text{BB}}` :math:`M_{\text{BB}}` after removal of
constrained DOF rows and columns

:math:`{\widetilde{M}}_{\text{BB}}` Substructure equivalent mass matrix
referred to the TP reference point =
:math:`T_{I}^{T}{\bar{M}}_{\text{BB}}T_{I}`

:math:`M_{\text{Bm}}` Matrix partition after C-B system reduction =
:math:`M_{\text{mB}}^{T}`

:math:`{\bar{M}}_{\text{Bm}}` :math:`M_{\text{Bm}}` after removal of
constrained DOF rows and columns

:math:`{\widetilde{M}}_{\text{Bm}}` Matrix partition =
:math:`T_{I}^{T}{\bar{M}}_{\text{Bm}}`

:math:`M_{\text{LL}}` Mass matrix partition referred to the *L* interior
DOFs (*L*\ ×\ *L*)

:math:`M_{\text{LR}}` Mass matrix partition referred to the *L* and *R*
DOFs (*L*\ ×\ *R*)

:math:`M_{\text{mB}}` Matrix partition after C-B system reduction =
:math:`\Phi_{m}^{T}M_{\text{LR}} + \Phi_{m}^{T}M_{\text{LL}}\Phi_{R}`

:math:`{\widetilde{M}}_{\text{mB}}` Matrix partition =
:math:`{\widetilde{M}}_{\text{Bm}}^{T}`

:math:`M_{\text{RL}}` Mass matrix partition referred to the *R* and *L*
boundary DOFs (*R*\ ×\ *L*)

:math:`M_{\text{RR}}` Mass matrix partition referred to the *R* boundary
DOFs (*R*\ ×\ *R*)

*M\ X* Substructure base reaction moment along global *X*-axis

*M\ Y* Substructure base reaction moment along global *Y*-axis

*M\ Z* Substructure base reaction moment along global *Z*-axis

:math:`M^{\left( \text{ED} \right)}` 6x6 mass matrix from ElastoDyn

:math:`M_{33}^{(SD)}` Substructure equivalent heave mass

:math:`n` The n\ :sup:`th` time step

:math:`\text{NIN}` Number of interface nodes

:math:`N_{\text{react}}` Number of restrained nodes

*p* State-space parameters

:math:`q_{L}` Modal coefficients for all interior nodes’ DOF modes

:math:`q_{L0}` Modal coefficients for all interior nodes’ DOF modes
assumed operating in static fashion

:math:`q_{m}` Modal coefficients for the retained modes

:math:`q_{m0}` Modal coefficients for the retained modes assumed
operating in static fashion

*R* Boundary nodes’ DOFs

:math:`\overrightarrow{R}` Reaction forces at the base of the
substructure

*SSI* Soil-structure interaction

:math:`t` time; or

thickness

TP Transition piece

:math:`T_{I}` Matrix to transform interface nodes’ DOFs to TP DOFs;
:math:`\left( 6 \cdot \text{NIN} \right) \times 6` matrix

:math:`T_{\text{react}}` Auxiliary matrix
(:math:`6x(6 \cdot N_{\text{react}})`) to link restrained nodes’ forces
and moments to the substructure base reaction vector

*u, u\ i* State-space formulation inputs, generic i-th input

:math:`U` Vector of nodal displacements

:math:`U_{e}` Vector of element nodes’ displacements (DOFs)

:math:`U_{L}` Vector of interior nodes’ displacements (DOFs)

:math:`{\widehat{U}}_{L}` Time-varying components of the interior nodes’
displacements (DOFs)

:math:`U_{L0}` Static components of the interior nodes’ displacements
(DOFs) (*L*\ ×1)

:math:`U_{L0m}` Static components of the interior nodes’ displacements
(DOFs) (*L*\ ×1), but obtained considering the first *m* C-B eigenmodes
only

:math:`U_{R}` Vector of boundary nodes’ displacements (DOFs)

:math:`{\bar{U}}_{R}` Vector of interface nodes’ displacements (DOFs)

:math:`U^{e}` Element nodes’ displacements (DOFs)

:math:`{\widehat{U}}^{e}` Time-varying components of the element nodes’
displacements (DOFs)

:math:`U_{L0}^{e}` Static components of the element nodes’ displacements
(DOFs)

:math:`U_{L0m}^{e}` Static components of the element nodes’
displacements (DOFs), but obtained considering the first *m* C-B
eigenmodes only

:math:`U_{\text{TP}}` TP reference point displacements (DOFs)

:math:`{\widehat{U}}_{\text{TP}}` Time-varying components of the TP
reference point displacements (DOFs)

:math:`U_{TP0}` Static components of the TP reference point
displacements (DOFs)

*x, x*\ :sub:`i` State-space formulation states, generic i-th state; or

generic local *x*-coordinate

:math:`x_{e}` Element local *x*-coordinates

*X* Global or substructure coordinate; or

state-space state equation(s)

:math:`X_{E}` Member end node *X*-coordinate in global coordinate system

:math:`X_{S}` Member start node *X*-coordinate in global coordinate
system

:math:`X_{\text{INi}}` *X*-coordinate in global coordinate system of the
generic interface node

*X\ SS* Global or substructure coordinate

:math:`X_{\text{TP}}` TP reference point *X*-coordinate in global
coordinate system

*y, y\ i* State-space formulation outputs, generic i-th output; or

generic *y*-coordinate

:math:`y_{e}` Element local *y*-coordinates

*Y* Global or substructure coordinates; or

state-space output equation(s)

:math:`Y_{E}` Member end node *Y*-coordinate in global coordinate system

:math:`Y_{\text{INi}}` *Y*-coordinate in global coordinate system of the
generic interface node

:math:`Y_{S}` Member start node *Y*-coordinate in global coordinate
system

*Y\ SS* Global or substructure coordinates

:math:`Y_{\text{TP}}` TP reference point *Y*-coordinate in global
coordinate system

*Y\ 1* State-space output equation for ElastoDyn

*Y\ 2* State-space output equation for HydroDyn

:math:`z_{e}` Element local *z*-coordinate

*Z* Global or substructure coordinate

:math:`Z_{E}` Member end node *Z*-coordinate in global coordinate system

:math:`Z_{\text{INi}}` *Z*-coordinate in global coordinate system of the
generic interface node

:math:`Z_{S}` Member start node *Z*-coordinate in global coordinate
system

*Z\ SS* Global or substructure coordinate

:math:`Z_{\text{TP}}` TP reference point *Z*-coordinate in global
coordinate system

Table of Contents

`Preface 4 <#_Toc413741142>`__

`Acknowledgments 5 <#_Toc413741143>`__

`Nomenclature 6 <#_Toc413741144>`__

`1 Introduction 14 <#_Toc401565516>`__

`2 Running SubDyn 17 <#_Ref391919518>`__

`2.1 Downloading the SubDyn Software 17 <#_Toc381965029>`__

`2.1.1 Stand-Alone SubDyn Archive 17 <#_Toc401565519>`__

`2.1.2 FAST Archive 17 <#_Toc401565520>`__

`2.2 Running SubDyn 18 <#_Toc413741150>`__

`2.2.1 Running the Stand-Alone SubDyn Program 18 <#_Toc413741151>`__

`2.2.2 Running SubDyn Coupled to FAST 18 <#_Toc382336008>`__

`3 Input Files 19 <#_Ref391919668>`__

`3.1 Units 19 <#_Toc401565524>`__

`3.2 SubDyn Driver Input File 19 <#_Toc381965070>`__

`3.3 SubDyn Primary Input File 20 <#_Ref382059081>`__

`3.3.1 Simulation Controls 21 <#_Toc401565527>`__

`3.3.2 FEA and Craig-Bampton Parameters 21 <#_Toc401565528>`__

`3.3.3 Structure Joints 22 <#_Toc401565529>`__

`3.3.4 Base Reaction Joints 22 <#_Toc401565530>`__

`3.3.5 Interface Joints 23 <#_Toc401565531>`__

`3.3.6 Members 23 <#_Ref391919995>`__

`3.3.7 Member Cross-Section Properties 23 <#_Toc519676721>`__

`3.3.8 Member Cosine Matrices COSM (i,j) 24 <#_Ref393270930>`__

`3.3.9 Joint Additional Concentrated Masses 24 <#_Ref393289423>`__

`3.3.10 Output: Summary and Outfile 24 <#_Toc401565536>`__

`3.3.11 Member Output List 25 <#_Ref393272400>`__

`3.3.12 Output Channels- SDOutList Section 25 <#_Toc401565538>`__

`4 Output Files 26 <#_Ref391919685>`__

`4.1 Echo Files 26 <#_Toc401565540>`__

`4.2 Summary File 26 <#_Ref391920238>`__

`4.3 Results File 27 <#_Toc401565542>`__

`5 Modeling Considerations 28 <#_Toc381965090>`__

`5.1 Model Discretization 28 <#_Toc401565544>`__

`5.2 Foundations 29 <#_Toc413741175>`__

`5.3 Member Overlap 29 <#_Toc413741176>`__

`5.4 Substructure Tower/Turbine Coupling 29 <#_Ref413052060>`__

`5.5 Self-Weight Calculations 31 <#_Toc401565546>`__

`5.6 Note On Other Load Calculations 31 <#_Toc413741179>`__

`5.7 Craig-Bampton Guidelines 31 <#_Ref391971454>`__

`5.8 Integration Time Step Guidelines 32 <#_Ref399231319>`__

`6 SubDyn Theory 33 <#_Ref394401650>`__

`6.1 Overview 33 <#_Toc401565550>`__

`6.2 Integration with the FAST Modularization Framework
34 <#_Toc401565551>`__

`6.3 Coordinate Systems 35 <#_Toc401565552>`__

`6.3.1 Global and Substructure Coordinate System: or ( Figure 4)
36 <#_Toc401565553>`__

`6.3.2 Member or Element Local Coordinate System (Figure 5)
36 <#_Toc401565554>`__

`6.3.3 Local to Global Transformation 37 <#_Toc401565555>`__

`6.4 Linear Finite-Element Beam Model 38 <#_Toc401565556>`__

`6.4.1 Element Formulation 39 <#_Toc401565557>`__

`6.4.2 Self-Weight Loads 40 <#_Toc401565558>`__

`6.5 Dynamic System of Equations and C-B Reduction
42 <#_Toc401565559>`__

`6.5.1 State-Space Formulation 46 <#_Ref394398934>`__

`6.5.2 Member Force Calculation 48 <#_Toc401565561>`__

`6.5.3 Reaction Calculation 48 <#_Toc401565562>`__

`6.5.4 Time Integration 49 <#_Ref394398970>`__

`6.5.5 Static-Improvement Method 49 <#_Toc401565564>`__

`7 Known Limitations and Future Work 53 <#_Toc401565565>`__

`8 References 54 <#_Ref391919755>`__

`Appendix A. OC4 Jacket Input File (CertTest Test04)
55 <#_Toc401565567>`__

`Appendix B. OC4 Jacket Driver File 60 <#_Toc401565568>`__

`Appendix C. List of Output Channels 61 <#_Toc401565569>`__

`Appendix D. Compiling Stand-Alone SubDyn 63 <#_Toc401565570>`__

`Appendix E. Major Changes in SubDyn 64 <#_Toc401565571>`__

List of Figures

`Figure 1. SubDyn, HydroDyn, and FAST 8 coupled interaction
15 <#_Ref408208410>`__

`Figure 2. WinZip Self-Extractor main window 17 <#_Toc413741206>`__

`Figure 3. SubDyn layout within the modularization framework
34 <#_Ref393719536>`__

`Figure 4. Global (coincident with the substructure) coordinate system.
Also shown are the DOFs associated with the TP reference point.
35 <#_Ref407623061>`__

`Figure 5. The element coordinate system. The sketched member contains
four elements, and the second element is called out with nodes S and E.
37 <#_Ref408218733>`__

List of Tables

`Table 1. TP Reference Point Inputs Time-Series Data File Contents
20 <#_Toc401063869>`__

`Table C-1. List of Output Channels. 61 <#_Toc417640980>`__


