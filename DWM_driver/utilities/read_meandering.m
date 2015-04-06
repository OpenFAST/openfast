clear all

Meandering_simulation_time_length = 300;  % Meandering_simulation_time_length
Meandering_Moving_time            = 60;   % Meandering_Moving_time
TurbSimHubHt                      = 205;  % TurbsimRefHt
HubHt                             = 70;   % HubHt


%% Read Menadering Center Position at each situation
center_position=textread('D:\V80_FAST8\DWM-results\WC_Turbine_1.txt'); % output .txt files from the DWM model

for i=1:3;
    wake_center_1(:,i) = center_position((i-1)*(Meandering_Moving_time+1)*Meandering_simulation_time_length+1 : i*(Meandering_Moving_time+1)*Meandering_simulation_time_length);
    for j=1:51;
        wake_center_2(:,j,i) = wake_center_1((j-1)*Meandering_simulation_time_length+1:j*Meandering_simulation_time_length,i);
        for k=1:80;
            wake_center(k,j,i) = wake_center_2(k,j,i);
        end
    end
end

wake_center(:,:,3) = wake_center(:,:,3) - TurbSimHubHt + HubHt;

%% plot the meandered wake center profile at each time step

% The wake meandering center is plotted with respect to "time steps". 

% ex. For flying time = 49(time steps)  ( => release time = 1, flying time = 49 ,
% => the wake center profile consists of the the cross-planes from wake_center(1,50,:)... to ....(50,1,:)


for n=1:41  % for flying time 40
    x_axis_40(n)=wake_center(42-n,n,1);    % the first cross scetion is @ turbine plane when n=1
    y_axis_40(n)=wake_center(42-n,n,2);
    z_axis_40(n)=wake_center(42-n,n,3);  % shift the height of the hub height to 90m
end

for n=1:31  % for flying time 30
    x_axis_30(n)=wake_center(32-n,n,1);    % the first cross scetion is @ turbine plane when n=1
    y_axis_30(n)=wake_center(32-n,n,2);
    z_axis_30(n)=wake_center(32-n,n,3); 
end

for n=1:11  % for flying time 10
    x_axis_10(n)=wake_center(12-n,n,1);    % the first cross scetion is @ turbine plane when n=1
    y_axis_10(n)=wake_center(12-n,n,2);
    z_axis_10(n)=wake_center(12-n,n,3); 
end

figure (1)
hold all
plot (x_axis_40(:),y_axis_40(:),'r','linewidth',2);  % wake center position at y coordinate
plot (x_axis_30(:),y_axis_30(:),'k','linewidth',2);  
plot (x_axis_10(:),y_axis_10(:),'b','linewidth',2);  
grid on
legend('Lateral Position @ flying time = 40','Lateral Position @ flying time = 30','Lateral Position @ flying time = 10');
xlabel('Downstream Distance [m]');
title('Wake Center Lateral Position @ Flying Time 40,30,10')

figure (2)
hold all
plot (x_axis_40(:),z_axis_40(:),'r','linewidth',2);  % wake center position at z coordinate
plot (x_axis_30(:),z_axis_30(:),'k','linewidth',2);
plot (x_axis_10(:),z_axis_10(:),'b','linewidth',2);
grid on
legend('Vertical Position @ flying time = 40','Vertical Position @ flying time = 30','Vertical Position @ flying time = 10');
xlabel('Downstream Distance [m]');
title('Wake Center Vertical Position @ Flying Time 40,30,10')