%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2016 AeroDyn standalone performance code wrapper
%  Robynne Murray, NREL
%
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

NumCases=16;

%% Run the AeroDyn model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('Running AeroDyn model.....');
display('  ');
command='AeroDyn_Driver_Win32 Bahaj_MHK_driver.dvr' ;   % Runs the driver in a command prompt window
status = dos(command);  %Zero means it worked!


%% Obtain performance results from output files %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjust these depending on what you want to output from the code
for i=1:NumCases
    fid = fopen( ['Bahaj_MHK_driver.' num2str(i) '.out'] );      %Open the output files in a loop depending on number of cases
    Read_OutData = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines',9);
    OutData=Read_OutData(1,:);
    OutData= cell2mat(OutData);
    
    TSR=mean(OutData(:,3));   %Take the mean at each case number (TSR) to compare rotor power and thrust for each TSR
    Cp=mean(OutData(:,6));
    Ct=mean(OutData(:,7));
    
    TSR_curve(i)=TSR;
    Cp_curve(i)=Cp;
    Ct_curve(i)=Ct;   %Save the mean values at each TSR in an array 
    
end

%% Experimental data to used as a comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bahaj experimental data U = 1.73m/s, pitch = 0 deg (Bahaj et al.,2007)
Cp_dataTSR173 = [4.17061600000000 4.42338100000000 4.66034800000000 4.89731400000000 5.13428100000000 5.37124800000000 5.37124800000000 5.59241700000000 5.84518200000000 6.08214800000000 6.30331800000000 6.54028400000000 6.77725100000000 7.01421800000000 7.21958900000000 7.44075800000000 7.69352300000000 8.15165900000000 8.60979500000000 9.06793000000000 9.54186400000000 9.87361800000000;];
Cp_data173 = [0.413793000000000 0.430885000000000 0.437181000000000 0.446177000000000 0.445277000000000 0.454273000000000 0.457871000000000 0.452474000000000 0.454273000000000 0.452474000000000 0.449775000000000 0.447976000000000 0.441679000000000 0.435382000000000 0.429085000000000 0.412894000000000 0.409295000000000 0.389505000000000 0.361619000000000 0.331034000000000 0.287856000000000 0.258171000000000;];
Ct_dataTSR173 = [4.18495300000000 4.43573700000000 4.67084600000000 4.90595600000000 5.14106600000000 5.39185000000000 5.37617600000000 5.61128500000000 5.86206900000000 5.87774300000000 6.09717900000000 6.09717900000000 6.31661400000000 6.55172400000000 6.78683400000000 7.02194400000000 7.25705300000000 7.46081500000000 7.71159900000000 8.16614400000000 8.63636400000000 9.09090900000000 9.54545500000000 9.89028200000000;];
Ct_data173 = [0.644510000000000 0.672997000000000 0.701484000000000 0.726409000000000 0.740653000000000 0.753116000000000 0.772700000000000 0.778042000000000 0.795846000000000 0.802967000000000 0.811869000000000 0.817211000000000 0.820772000000000 0.836795000000000 0.851039000000000 0.859941000000000 0.874184000000000 0.884866000000000 0.890208000000000 0.902671000000000 0.915134000000000 0.925816000000000 0.925816000000000 0.938279000000000;];

%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(TSR_curve,Cp_curve,'-k')
grid on
hold on
scatter(Cp_dataTSR173,Cp_data173,'k','d', 'LineWidth',1.)
% axis([ 3 10 0 1])
plot(TSR_curve,Ct_curve,'--k')
scatter(Ct_dataTSR173,Ct_data173,'k', 'LineWidth',1.)
title('Cp and Ct-TSR for Bahaj Turbine ','FontSize',12,'FontName','Times')
xlabel('TSR','FontSize',14,'FontName','Times')
ylabel('C_P  and  C_T','FontSize',14,'FontName','Times')
legend('AeroDyn C_P','Experiment C_P','AeroDyn C_T','Experiment C_T')
