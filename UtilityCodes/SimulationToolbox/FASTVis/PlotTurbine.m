% Function for plotting the entire turbine
%
%
%
% Knud A. Kragh

function PlotTurbine(FastPar,TowerData,Chord_front1,C2_def1,Chord_back1,Chord_front2,C2_def2,Chord_back2,Chord_front3,C2_def3,Chord_back3,R,R2)

PlotTower(TowerData,FastPar,R)
hold on

PlotShaft(FastPar,R2)
PlotRotor(Chord_front1,C2_def1,Chord_back1,Chord_front2,C2_def2,Chord_back2,Chord_front3,C2_def3,Chord_back3,GetFastPar(FastPar,'OverHang'),GetFastPar(FastPar,'TowerHt')+GetFastPar(FastPar,'Twr2Shft'))

