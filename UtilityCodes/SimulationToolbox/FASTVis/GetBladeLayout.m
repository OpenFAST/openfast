% Function for Calculating blade layout of blade modeled in FAST
%
%
%
% Knud A. Kragh

function [Chord_front C2_def Chord_back]=GetBladeLayout(BladeData,Hawc2_Ae,FastPar)

%C2_def=[-BladeData.Data(:,13) BladeData.Data(:,12) round(100*((GetFastPar(FastPar,'TipRad')-GetFastPar(FastPar,'HubRad'))*BladeData.Data(:,1)))/100];
R=round(100*((GetFastPar(FastPar,'TipRad')-GetFastPar(FastPar,'HubRad'))*BladeData.BldProp(:,1)))/100;
C2_def=[zeros(size(R)) zeros(size(R)) R];
twist=-BladeData.BldProp(:,3);
L_chord=interp1(Hawc2_Ae(:,1),Hawc2_Ae(:,2),C2_def(:,3));

% Modify C2_def to have straight leading edge
InitRad=L_chord(1)/2;
C2_def(2:end,1)=C2_def(2:end,1)-(L_chord(2:end)/2-InitRad);

% Apply pre bent and pre sweep to c2_def
C2_def(:,1)=C2_def(:,1)-BladeData.BldProp(:,13);
C2_def(:,2)=C2_def(:,2)+BladeData.BldProp(:,12);

% Chord coord. UN TWISTED
Chord_front(:,1)=C2_def(:,1)+L_chord/2;
Chord_front(:,2)=0;
Chord_front(:,3)=C2_def(:,3);

Chord_back(:,1)=C2_def(:,1)-L_chord/2;
Chord_back(:,2)=0;
Chord_back(:,3)=C2_def(:,3);

% Applying Twist
Chord_front(:,1)=Chord_front(:,1).*cos(twist*pi/180);
Chord_front(:,2)=L_chord/2.*sin(twist*pi/180);

Chord_back(:,1)=Chord_back(:,1).*cos(twist*pi/180);
Chord_back(:,2)=-L_chord/2.*sin(twist*pi/180);




