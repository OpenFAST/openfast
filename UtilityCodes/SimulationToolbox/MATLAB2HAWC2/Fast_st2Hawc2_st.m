% Function for converting structural data from FAST to Hawc2 format (blade
% and tower) and for defining dummy data for shaft and hub.
%
% In:   FastPar     -   Fast input parameters in MATLAB structure extracted
%                       using Fast2Matlab function
%       TowerData   -   FAST Tower Data in matlab structure
%       BladeData   -   FAST Blade Data in matlab structure
%       Hawc2Set    -   Hawc2 additional settings
%
% Out:  OutData     -   Structure with data for all main bodies in Hawc2
%                       format
%
% Knud A. Kragh

function [OutData]=Fast_st2Hawc2_st(FastPar,TowerData,BladeData,Hawc2Set)

%% Tower

Ix=TowerData.TowProp(:,3)*Hawc2Set.Tower.FlapStiffTune./Hawc2Set.Tower.E;
Iy=TowerData.TowProp(:,4)*Hawc2Set.Tower.EdgeStiffTune./Hawc2Set.Tower.E;

if sum(TowerData.TowProp(:,5))<1e-5
    J=Hawc2Set.Tower.J;
else
    J=TowerData.TowProp(:,5)/Hawc2Set.TowProp.G;
end

A=TowerData.TowProp(:,2)./Hawc2Set.Tower.rho;

OutData.Tower(:,1)=(GetFastPar(FastPar,'TowerHt')).*TowerData.TowProp(:,1);              % Length
OutData.Tower(:,2)=TowerData.TowProp(:,2);              % Mass per unit length
OutData.Tower(:,3)=TowerData.TowProp(:,10);              % xm, xc2-coordinate from C1/2 to mass center [m]
OutData.Tower(:,4)=TowerData.TowProp(:,9);              % ym, xc2-coordinate from C1/2 to mass center [m]
OutData.Tower(:,5)=sqrt(Ix./A);              % rix, radius of inertia related to elastic center.
OutData.Tower(:,6)=sqrt(Iy./A);              % riy, radius of inertia related to elastic center.
OutData.Tower(:,7)=TowerData.TowProp(:,10);              % xs, xc2-coordinate from C1/2 to shear center [m]
OutData.Tower(:,8)=TowerData.TowProp(:,9);              % ys, xc2-coordinate from C1/2 to shear center [m]
OutData.Tower(:,9)=Hawc2Set.Tower.E;              % E
OutData.Tower(:,10)=Hawc2Set.Tower.G;              % G
OutData.Tower(:,11)=Ix;              % Ix, area moment of inertia with respect to principal bending xe axis [m4]
OutData.Tower(:,12)=Iy;              % Iy, area moment of inertia with respect to principal bending ye axis [m4]
OutData.Tower(:,13)=J;              % K, torsional stiffness CORRECT TO GET FAST INP.
OutData.Tower(:,14)=Hawc2Set.Tower.kx;              % kx shear factor for force in principal bending xe direction [-]
OutData.Tower(:,15)=Hawc2Set.Tower.ky;              % ky shear factor for force in principal bending ye direction [-]
OutData.Tower(:,16)=A;              % A
OutData.Tower(:,17)=0.0;              % Struct pitch
OutData.Tower(:,18)=TowerData.TowProp(:,10);              % xe, xc2-coordinate from C1/2 to center of elasticity [m]
OutData.Tower(:,19)=TowerData.TowProp(:,9);              % xe, yc2-coordinate from C1/2 to center of elasticity [m]

%% Adding top section to tower

OutData.Tower(end+1:end+2,1)=[GetFastPar(FastPar,'TowerHt')+0.1 GetFastPar(FastPar,'TowerHt')+GetFastPar(FastPar,'Twr2Shft')];              % Length
OutData.Tower(end-1:end,2)=1;              % Mass per unit length
OutData.Tower(end-1:end,3)=0;              % xm, xc2-coordinate from C1/2 to mass center [m]
OutData.Tower(end-1:end,4)=0;              % ym, xc2-coordinate from C1/2 to mass center [m]
OutData.Tower(end-1:end,5)=1;              % rix, radius of inertia related to elastic center.
OutData.Tower(end-1:end,6)=1;              % riy, radius of inertia related to elastic center.
OutData.Tower(end-1:end,7)=0;              % xs, xc2-coordinate from C1/2 to shear center [m]
OutData.Tower(end-1:end,8)=0;              % ys, xc2-coordinate from C1/2 to shear center [m]
OutData.Tower(end-1:end,9)=2.11e11;              % E
OutData.Tower(end-1:end,10)=8e10;              % G
OutData.Tower(end-1:end,11)=1e4;              % Ix, area moment of inertia with respect to principal bending xe axis [m4]
OutData.Tower(end-1:end,12)=1e4;              % Iy, area moment of inertia with respect to principal bending ye axis [m4]
OutData.Tower(end-1:end,13)=5;              % K, torsional stiffness CORRECT TO GET FAST INP.
OutData.Tower(end-1:end,14)=0.5;              % kx shear factor for force in principal bending xe direction [-]
OutData.Tower(end-1:end,15)=0.5;              % ky shear factor for force in principal bending ye direction [-]
OutData.Tower(end-1:end,16)=10;              % A
OutData.Tower(end-1:end,17)=0.0;              % Struct pitch
OutData.Tower(end-1:end,18)=0;              % xe, xc2-coordinate from C1/2 to center of elasticity [m]
OutData.Tower(end-1:end,19)=0;              % xe, yc2-coordinate from C1/2 to center of elasticity [m]


%% Shaft
Ix=Hawc2Set.Shaft.Ix;
Iy=Hawc2Set.Shaft.Iy;
A=Hawc2Set.Shaft.A;

OutData.Shaft=zeros(3,19);
OutData.Shaft(:,1)=[0 abs(GetFastPar(FastPar,'OverHang'))/2 abs(GetFastPar(FastPar,'OverHang'))];              % Length
OutData.Shaft(:,2)=1;              % Mass per unit length (dummy)
OutData.Shaft(:,3)=0.0;              % xm, xc2-coordinate from C1/2 to mass center [m]
OutData.Shaft(:,4)=0.0;              % ym, xc2-coordinate from C1/2 to mass center [m]
OutData.Shaft(:,5)=sqrt(Ix./A);              % rix, radius of inertia related to elastic center.
OutData.Shaft(:,6)=sqrt(Iy./A);               % riy, radius of inertia related to elastic center.
OutData.Shaft(:,7)=0.0;              % xs, xc2-coordinate from C1/2 to shear center [m]
OutData.Shaft(:,8)=0.0;              % ys, xc2-coordinate from C1/2 to shear center [m]
OutData.Shaft(:,9)=Hawc2Set.Shaft.E;              % E
OutData.Shaft(:,10)=Hawc2Set.Shaft.G;              % G
OutData.Shaft(:,11)=Ix;              % Ix, area moment of inertia with respect to principal bending xe axis [m4]
OutData.Shaft(:,12)=Iy;              % Iy, area moment of inertia with respect to principal bending ye axis [m4]
OutData.Shaft(:,13)=GetFastPar(FastPar,'DTTorSpr')./OutData.Shaft(:,10)*-GetFastPar(FastPar,'OverHang');              % K, torsional stiffness
OutData.Shaft(:,14)=Hawc2Set.Shaft.kx;              % kx shear factor for force in principal bending xe direction [-]
OutData.Shaft(:,15)=Hawc2Set.Shaft.ky;              % ky shear factor for force in principal bending ye direction [-]
OutData.Shaft(:,16)=A;              % A
OutData.Shaft(:,17)=0.0;              % Struct pitch
OutData.Shaft(:,18)=0.0;              % xe, xc2-coordinate from C1/2 to center of elasticity [m]
OutData.Shaft(:,19)=0.0;              % xe, yc2-coordinate from C1/2 to center of elasticity [m]

%% Hub
OutData.Hub=zeros(2,19);
OutData.Hub(:,1)=[0 GetFastPar(FastPar,'HubRad')];              % Length
OutData.Hub(:,2)=1;              % Mass per unit length (dummy)
OutData.Hub(:,3)=0.0;              % xm, xc2-coordinate from C1/2 to mass center [m]
OutData.Hub(:,4)=0.0;              % ym, xc2-coordinate from C1/2 to mass center [m]
OutData.Hub(:,5)=1;              % rix, radius of inertia related to elastic center.
OutData.Hub(:,6)=1;              % riy, radius of inertia related to elastic center.
OutData.Hub(:,7)=0.0;              % xs, xc2-coordinate from C1/2 to shear center [m]
OutData.Hub(:,8)=0.0;              % ys, xc2-coordinate from C1/2 to shear center [m]
OutData.Hub(:,9)=2.11e11;              % E
OutData.Hub(:,10)=8e10;              % G
OutData.Hub(:,11)=129;              % Ix, area moment of inertia with respect to principal bending xe axis [m4]
OutData.Hub(:,12)=129;              % Iy, area moment of inertia with respect to principal bending ye axis [m4]
OutData.Hub(:,13)=50;              % K, torsional stiffness
OutData.Hub(:,14)=5;              % kx shear factor for force in principal bending xe direction [-]
OutData.Hub(:,15)=5;              % ky shear factor for force in principal bending ye direction [-]
OutData.Hub(:,16)=1;              % A
OutData.Hub(:,17)=0.0;              % Struct pitch
OutData.Hub(:,18)=0.0;              % xe, xc2-coordinate from C1/2 to center of elasticity [m]
OutData.Hub(:,19)=0.0;              % xe, yc2-coordinate from C1/2 to center of elasticity [m]

%% Blade

Ix=BladeData.BldProp(:,5)*Hawc2Set.Blade.FlapStiffTune./Hawc2Set.Blade.E;
Iy=BladeData.BldProp(:,6)*Hawc2Set.Blade.EdgeStiffTune./Hawc2Set.Blade.E;

if sum(BladeData.BldProp(:,7))<1e-5
    J=Hawc2Set.Blade.J;
else
    J=BladeData.BldProp(:,7)/Hawc2Set.Blade.G;
end
A=BladeData.BldProp(:,4)./Hawc2Set.Blade.rho;

OutData.Blade(:,1)=round(100*((GetFastPar(FastPar,'TipRad')-GetFastPar(FastPar,'HubRad'))*BladeData.BldProp(:,1)))/100;              % Length
OutData.Blade(:,2)=BladeData.BldProp(:,4);              % Mass per unit length
OutData.Blade(:,3)=BladeData.BldProp(:,15);              % xm, xc2-coordinate from C1/2 to mass center [m]
OutData.Blade(:,4)=BladeData.BldProp(:,14);              % ym, xc2-coordinate from C1/2 to mass center [m]
OutData.Blade(:,5)=sqrt(Ix./A);             % rix, radius of inertia related to elastic center.
OutData.Blade(:,6)=sqrt(Iy./A);             % riy, radius of inertia related to elastic center.
OutData.Blade(:,7)=BladeData.BldProp(:,17);              % xs, xc2-coordinate from C1/2 to shear center [m]
OutData.Blade(:,8)=BladeData.BldProp(:,16);              % ys, xc2-coordinate from C1/2 to shear center [m]
OutData.Blade(:,9)=Hawc2Set.Blade.E;              % E
OutData.Blade(:,10)=Hawc2Set.Blade.G;              % G
OutData.Blade(:,11)=Ix;              % Ix, area moment of inertia with respect to principal bending xe axis [m4]
OutData.Blade(:,12)=Iy;              % Iy, area moment of inertia with respect to principal bending ye axis [m4]
OutData.Blade(:,13)=J;              % K, torsional stiffness
OutData.Blade(:,14)=Hawc2Set.Blade.kx;              % kx shear factor for force in principal bending xe direction [-]
OutData.Blade(:,15)=Hawc2Set.Blade.ky;              % ky shear factor for force in principal bending ye direction [-]
OutData.Blade(:,16)=A;              % A
OutData.Blade(:,17)=0.0;              % Struct pitch
OutData.Blade(:,18)=BladeData.BldProp(:,17);              % xe, xc2-coordinate from C1/2 to center of elasticity [m]
OutData.Blade(:,19)=BladeData.BldProp(:,16);              % ye, yc2-coordinate from C1/2 to center of elasticity [m]