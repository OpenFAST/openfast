%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2016 AeroDyn standalone performance code wrapper
%  Updated by Robynne Murray, NREL, February 2017
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

NumCases=3;

TimeSteps=80;

%% Run the AeroDyn model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Running AeroDyn model.....');
display('  ');
command='AeroDyn_Driver_Win32 NRELOffshrBsline5MW_Onshore_Test01.dvr' ;   % Runs the driver in a command prompt window
status = dos(command);  %Zero means it worked!


%% Obtain performance results from output files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust these depending on what you want to output from the code
for i=1:NumCases
    fid = fopen( ['AOC_Test01.' num2str(i) '.out'] );
    Read_OutData = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines',9);
    OutData=Read_OutData(1,:);
    OutData= cell2mat(OutData);
    
    Time_curve=(OutData(:,1));
    Rt_Speed_curve=(OutData(:,2));
    B1N1Phi_curve=(OutData(:,5));
    B1N2Phi_curve=(OutData(:,6));
    
% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust these based on what outputs you want to look at 
figure(1)
plot(Time_curve,B1N1Phi_curve,'-k')
grid on
hold on
plot(Time_curve,B1N2Phi_curve,'-r')
title('time-phi')
xlabel('Time')
ylabel('Phi')
legend('Blade 1 Node 1 Phi','Blade 1 Node 2 Phi')



end


