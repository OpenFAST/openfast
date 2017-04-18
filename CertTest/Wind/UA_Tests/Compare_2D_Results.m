function Compare_2D_Results(localFolder, NRELFolder, ExpFolder, legendStr)
% Compare NREL's UnsteadyAero Models with experimental data for the:
%   NACA, LS417, and S809 airfoils.  See the following reference for
%   details:
%
% Rick Damiani, Greg Hayman, and Jason M. Jonkman. 
% "Development and Validation of a New Unsteady Airfoil Aerodynamics Model Within AeroDyn", 
% 34th Wind Energy Symposium, AIAA SciTech, (AIAA 2016-1007). 
%
if (nargin < 4)
   legendStr = {'Experimental','Local-UAMod2','Local-UAMod3','NREL-UAMod2','NREL-UAMod3'};
end
if (nargin < 1)
    localFolder = '.\Results';
end 
if (nargin < 2)
    NRELFolder = '.\NREL_Results';
end
if (nargin < 3)
    ExpFolder = '.\Exp_Results';
end

fileExp = [ExpFolder '\05014051_Validacion.txt'];
fileSim2 = [localFolder '\05014051_NACA_UAMod2.out'];
fileSimNREL2 = [NRELFolder '\05014051_NACA_UAMod2.out'];
fileSim3 = [localFolder '\05014051_NACA_UAMod3.out'];
fileSimNREL3 = [NRELFolder '\05014051_NACA_UAMod3.out'];
plot2D_results(fileExp, fileSim2, fileSimNREL2, fileSim3, fileSimNREL3,[.2,2.2],[-.1,1.2],[-.5,.2],legendStr, 'NACA 0015, k=0.077, Re=1.48e6, Mean AoA = 15 deg');

fileExp = [ExpFolder '\05014061_Validacion.txt'];
fileSim2 = [localFolder '\05014061_NACA_UAMod2.out'];
fileSimNREL2 = [NRELFolder '\05014061_NACA_UAMod2.out'];
fileSim3 = [localFolder '\05014061_NACA_UAMod3.out'];
fileSimNREL3 = [NRELFolder '\05014061_NACA_UAMod3.out'];
plot2D_results(fileExp, fileSim2, fileSimNREL2, fileSim3, fileSimNREL3,[.2,2.2],[-.1,1.2],[-.5,.2],legendStr, 'NACA 0015, k=0.077, Re=1.48e6, Mean AoA = 20 deg');

%LS
fileExp = [ExpFolder '\C10m100AM14_l417_Validacion.txt'];
fileSim2 = [localFolder '\C10m100AM14_LS417_UAMod2.out'];
fileSimNREL2 = [NRELFolder '\C10m100AM14_LS417_UAMod2.out'];
fileSim3 = [localFolder '\C10m100AM14_LS417_UAMod3.out'];
fileSimNREL3 = [NRELFolder '\C10m100AM14_LS417_UAMod3.out'];
plot2D_results(fileExp, fileSim2, fileSimNREL2, fileSim3, fileSimNREL3,[0.5,2.2],[-.05,0.45],[-.35,.05],legendStr,'LS(1)-0417MOD, k=0.077, Re=1.0e6');

fileExp = [ExpFolder '\C10m150AM14_l417_Validacion.txt'];
fileSim2 = [localFolder '\C10m150AM14_LS417_UAMod2.out'];
fileSimNREL2 = [NRELFolder '\C10m150AM14_LS417_UAMod2.out'];
fileSim3 = [localFolder '\C10m150AM14_LS417_UAMod3.out'];
fileSimNREL3 = [NRELFolder '\C10m150AM14_LS417_UAMod3.out'];
plot2D_results(fileExp, fileSim2, fileSimNREL2, fileSim3, fileSimNREL3,[0.5,2.2],[-.05,0.45],[-.35,.05],legendStr,'LS(1)-0417MOD, k=0.077, Re=1.5e6');

%S809
fileExp = [ExpFolder '\C10l100AM20_s809_Validacion.txt'];
fileSim2 = [localFolder '\C10l100AM20_s809_UAMod2.out'];
fileSimNREL2 = [NRELFolder '\C10l100AM20_s809_UAMod2.out'];
fileSim3 = [localFolder '\C10l100AM20_s809_UAMod3.out'];
fileSimNREL3 = [NRELFolder '\C10l100AM20_s809_UAMod3.out'];
plot2D_results(fileExp, fileSim2, fileSimNREL2, fileSim3, fileSimNREL3,[0.4,2.0],[-.05,0.9],[-.4,.15],legendStr, 'S809, k=0.025, Re=1.0e6');

fileExp = [ExpFolder '\C10m100AM20_s809_Validacion.txt'];
fileSim2 = [localFolder '\C10m100AM20_s809_UAMod2.out'];
fileSimNREL2 = [NRELFolder '\C10m100AM20_s809_UAMod2.out'];
fileSim3 = [localFolder '\C10m100AM20_s809_UAMod3.out'];
fileSimNREL3 = [NRELFolder '\C10m100AM20_s809_UAMod3.out'];
plot2D_results(fileExp, fileSim2, fileSimNREL2, fileSim3, fileSimNREL3,[0.4,2.0],[-.05,0.9],[-.4,.15],legendStr, 'S809, k=0.05, Re=1.0e6');

return;

end

function plot2D_results(fileExp, fileSim2, fileSimNREL2, fileSim3, fileSimNREL3, yLimitsCl, yLimitsCd, yLimitsCm, legendStr, titleStr) %, fileSimADv15)

   % Experimental Data
   [dataFile] = importdata(fileExp);
   data = dataFile.data;

   aoa_exp  = data(:,1);
   cl_exp   = data(:,2);
   cd_exp   = data(:,3);
   cm_exp   = data(:,4);

   % Local UA Stand-alone
   [dataFile] = importdata(fileSim2);
   data = dataFile.data;
   time = data(:,1);
   aoa_local2  = data(:,2);
   cl_local2   = data(:,6);
   cd_local2   = data(:,7);
   cm_local2   = data(:,8);
   
   [dataFile] = importdata(fileSimNREL2);
   data = dataFile.data;
   aoa_NREL2  = data(:,2);
   cl_NREL2   = data(:,6);
   cd_NREL2   = data(:,7);
   cm_NREL2   = data(:,8);
   
   [dataFile] = importdata(fileSim3);
   data = dataFile.data;
   time = data(:,1);
   aoa_local3  = data(:,2);
   cl_local3   = data(:,6);
   cd_local3   = data(:,7);
   cm_local3   = data(:,8);
   
   [dataFile] = importdata(fileSimNREL3);
   data = dataFile.data;
   aoa_NREL3  = data(:,2);
   cl_NREL3   = data(:,6);
   cd_NREL3   = data(:,7);
   cm_NREL3   = data(:,8);
    
   figure
   plot(aoa_exp,cl_exp,'k-','Marker','.');
   hold on;
   plot(aoa_local2(500:end),cl_local2(500:end),'m','LineWidth',2)
   plot(aoa_local3(500:end),cl_local3(500:end),'b','LineWidth',2)
   grid;
   ylabel('Cl');    
   xlabel('AOA');
   title(['Cl vs AOA : ' titleStr]);
   plot(aoa_NREL2(500:end),cl_NREL2(500:end),'r','LineWidth',2) 
   plot(aoa_NREL3(500:end),cl_NREL3(500:end),'c','LineWidth',2)
   legend(legendStr,'Location','northwest');   
   ylim(yLimitsCl);
   
   figure
   plot(aoa_exp,cd_exp,'k-','Marker','.');
   hold on;
   plot(aoa_local2(500:end),cd_local2(500:end),'m','LineWidth',2)
   plot(aoa_local3(500:end),cd_local3(500:end),'b','LineWidth',2)
   grid;
   ylabel('Cd');    
   xlabel('AOA');
   title(['Cd vs AOA : ' titleStr]);
   plot(aoa_NREL2(500:end),cd_NREL2(500:end),'r','LineWidth',2) 
   plot(aoa_NREL3(500:end),cd_NREL3(500:end),'c','LineWidth',2)
   legend(legendStr,'Location','northwest');   
   ylim(yLimitsCd);
   
   figure
   plot(aoa_exp,cm_exp,'k-','Marker','.');
   hold on;
   plot(aoa_local2(500:end),cm_local2(500:end),'m','LineWidth',2)
   plot(aoa_local3(500:end),cm_local3(500:end),'b','LineWidth',2)
   grid;
   ylabel('Cm');    
   xlabel('AOA');
   title(['Cm vs AOA : ' titleStr]);
   plot(aoa_NREL2(500:end),cm_NREL2(500:end),'r','LineWidth',2) 
   plot(aoa_NREL3(500:end),cm_NREL3(500:end),'c','LineWidth',2)
   legend(legendStr,'Location','northwest');   
   ylim(yLimitsCm);
   
   
   return;
end
