% Function for plotting rotor
%
%
%
% Knud A. Kragh


function PlotRotor(Chord_front1,C2_def1,Chord_back1,Chord_front2,C2_def2,Chord_back2,Chord_front3,C2_def3,Chord_back3,Y_offset,Z_offset)

RotAng2=120;
RotAng2=RotAng2*pi/180;

RotAng3=240;
RotAng3=RotAng3*pi/180;

TransMat2=[cos(RotAng2)   0   -sin(RotAng2)
    0             1   0
    sin(RotAng2)   0   cos(RotAng2)];

TransMat3=[cos(RotAng3)   0   -sin(RotAng3)
    0             1   0
    sin(RotAng3)   0   cos(RotAng3)];

for i=1:size(Chord_front1,1)
    temp=TransMat2*Chord_front2(i,:)';
    Chord_front2(i,:)=temp';
    temp=TransMat3*Chord_front3(i,:)';
    Chord_front3(i,:)=temp';
    
    
    temp=TransMat2*C2_def2(i,:)';
    C2_def2(i,:)=temp';
    temp=TransMat3*C2_def3(i,:)';
    C2_def3(i,:)=temp';
    
    temp=TransMat2*Chord_back2(i,:)';
    Chord_back2(i,:)=temp';
    temp=TransMat3*Chord_back3(i,:)';
    Chord_back3(i,:)=temp';
end

Chord_front1(:,2)=Chord_front1(:,2)+Y_offset;
Chord_front2(:,2)=Chord_front2(:,2)+Y_offset;
Chord_front3(:,2)=Chord_front3(:,2)+Y_offset;

C2_def1(:,2)=C2_def1(:,2)+Y_offset;
C2_def2(:,2)=C2_def2(:,2)+Y_offset;
C2_def3(:,2)=C2_def3(:,2)+Y_offset;

Chord_back1(:,2)=Chord_back1(:,2)+Y_offset;
Chord_back2(:,2)=Chord_back2(:,2)+Y_offset;
Chord_back3(:,2)=Chord_back3(:,2)+Y_offset;

Chord_front1(:,3)=Chord_front1(:,3)+Z_offset;
Chord_front2(:,3)=Chord_front2(:,3)+Z_offset;
Chord_front3(:,3)=Chord_front3(:,3)+Z_offset;

C2_def1(:,3)=C2_def1(:,3)+Z_offset;
C2_def2(:,3)=C2_def2(:,3)+Z_offset;
C2_def3(:,3)=C2_def3(:,3)+Z_offset;

Chord_back1(:,3)=Chord_back1(:,3)+Z_offset;
Chord_back2(:,3)=Chord_back2(:,3)+Z_offset;
Chord_back3(:,3)=Chord_back3(:,3)+Z_offset;

PlotBlade(Chord_front1,C2_def1,Chord_back1)
hold on
PlotBlade(Chord_front2,C2_def2,Chord_back2)
hold on
PlotBlade(Chord_front3,C2_def3,Chord_back3)

% str1(1) = {'Fast Viz. by K. A. Kragh'};
% text(10,10,str1)
title('Fast Vis. by K. A. Kragh','fontsize',14)

