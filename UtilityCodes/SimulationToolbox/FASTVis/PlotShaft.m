% Function for plotting shaft
%
%
% Knud A. Kragh

function PlotShaft(FastPar,R)

ShaftHeight=GetFastPar(FastPar,'TowerHt')+GetFastPar(FastPar,'Twr2Shft');
ShaftLength=GetFastPar(FastPar,'OverHang');

[X,Y,Z]=cylinder(R);
Z=Z*ShaftLength;

surf(X,Z,Y+ShaftHeight)
