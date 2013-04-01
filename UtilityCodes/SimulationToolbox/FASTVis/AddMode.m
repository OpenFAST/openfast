% Adding modeshape to blade or tower
%
%
%
% Knud A. Kragh

function [Chord_front C2_def Chord_back]=AddMode(Chord_frontOri,C2_defOri,Chord_backOri,dirMode,dirLengthAxis,coefs,Amp)

Chord_front=Chord_frontOri;
C2_def=C2_defOri;
Chord_back=Chord_backOri;

NormLength=abs(C2_defOri(:,dirLengthAxis))/max(abs(C2_defOri(:,dirLengthAxis)));

for i=1:length(NormLength)
    Mode(i)=coefs(1)*NormLength(i).^2+coefs(2)*NormLength(i).^3+coefs(3)*NormLength(i).^4+coefs(4)*NormLength(i).^5+coefs(5)*NormLength(i).^6;
end
Mode=Mode*Amp;

Chord_front(:,dirMode)=Chord_front(:,dirMode)+Mode';
C2_def(:,dirMode)=C2_def(:,dirMode)+Mode';
Chord_back(:,dirMode)=Chord_back(:,dirMode)+Mode';

% [NormLength Mode']
% 
% figure(10)
% plot(NormLength,Mode)