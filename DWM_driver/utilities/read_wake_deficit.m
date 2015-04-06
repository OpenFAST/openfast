clear all
%% read the wake deficit file
% Read Wind File, wind file is in .txt format, this is a example
speed_profile=textread('D:\V80_FAST8\DWM-results\WakeU_Turbine_3.txt'); %

ppr      = 50;
domain_r = 5;
domain_x = 36;

for i=1:ppr*domain_r;
    WakeU_profile(i,:)=speed_profile((i-1)*domain_x/2*ppr+1:i*domain_x/2*ppr);
end

%% plot the wake deficit profile at certain downstream location

% U_profile(:,3*ppr+1) means the axisymmetric wake profile is
% at the 3D downstream position

figure (1)   
hold all
plot(0:0.02:size(WakeU_profile(:,3*ppr+1))/50-0.02,WakeU_profile(:,3*ppr+1),'b--','linewidth',3);  % @3D
plot(0:0.02:size(WakeU_profile(:,4*ppr+1))/50-0.02,WakeU_profile(:,4*ppr+1),'g','linewidth',3);    % @4D
plot(0:0.02:size(WakeU_profile(:,5*ppr+1))/50-0.02,WakeU_profile(:,5*ppr+1),'r','linewidth',3);    % @5D
plot(0:0.02:size(WakeU_profile(:,6*ppr+1))/50-0.02,WakeU_profile(:,6*ppr+1),'c','linewidth',3);    % @6D
plot(0:0.02:size(WakeU_profile(:,7*ppr+1))/50-0.02,WakeU_profile(:,7*ppr+1),'m','linewidth',3);    % @7D
plot(0:0.02:size(WakeU_profile(:,8*ppr+1))/50-0.02,WakeU_profile(:,8*ppr+1),'k','linewidth',3);    % @8D
axis([0 2.0 0.4 1.1])
grid on
title('Wake Deficit Profile')
ylabel('Normalized Axial velocity [/]')
xlabel('Radial position [R]')
legend('@ downstream 3D','@ downstream 4D','@ downstream 5D','@ downstream 6D','@ downstream 7D','@ downstream 8D')







