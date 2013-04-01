% Function for Plotting Blade layout of blade modeled in FAST
%
%
%
% Knud A. Kragh

function PlotBlade(Chord_front,C2_def,Chord_back)

plot3(C2_def(:,1),C2_def(:,2),C2_def(:,3),'r--')
hold on
plot3(Chord_front(:,1),Chord_front(:,2),Chord_front(:,3))
plot3(Chord_back(:,1),Chord_back(:,2),Chord_back(:,3))
axis equal


% Connecting Leading and trailing edge
for i=1:size(Chord_front,1)
Connect(1,:)=[Chord_back(i,1) Chord_back(i,2) Chord_back(i,3)];
Connect(2,:)=[Chord_front(i,1) Chord_front(i,2) Chord_front(i,3)];
plot3(Connect(:,1),Connect(:,2),Connect(:,3))
end
hold off
title('Fast Vis. by K. A. Kragh','fontsize',14)

