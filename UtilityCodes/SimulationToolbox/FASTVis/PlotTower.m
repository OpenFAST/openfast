% Function for plotting tower
%
%
%
% Knud A. Kragh

function PlotTower(TowerData,FastPar,R)

height=GetFastPar(FastPar,'TowerHt')+GetFastPar(FastPar,'Twr2Shft');

[X,Y,Z]=cylinder(R);
Z=Z*height;

surf(X,Y,Z)

