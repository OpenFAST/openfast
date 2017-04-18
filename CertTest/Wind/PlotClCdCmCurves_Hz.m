
function PlotClCdCmCurves_Hz(Index, baseName, baseFolder)
if (nargin < 2)
  Index = 2;  % Valid values are 1, 2, 3
end 

if (nargin < 2)
   baseName = 'AOC';
   %baseName = 'NRELOffshrBsline5MW_Onshore';
end
if (nargin < 3)
   baseFolder = '.\NREL_Results';
   %baseFolder = '.\Results';
end
LineColors     = {[0 0 0],[0 1 1],[1 0 1],[0 1 0],[0 0 1],[1 0 0]};
Markers        = {'o'    ,'s'    ,'d'    ,'v'    ,'^'    ,'.'    };
delim='';
HeaderRows= 8;
NameLine= 7;
UnitsLine= 8;
legendStr = {'1000 Hz', '20 Hz', '10 Hz'};
   % Loop over different timestep sizes
for iCase=3:-1:1
   FileName = [baseFolder '\' baseName '_Test0' num2str(iCase) '.' num2str(Index) '.out'];
   [Channels, ChanName, ChanUnit,DescStr] = ReadFASTtext(FileName, delim, HeaderRows, NameLine, UnitsLine );
   
   if ( iCase == 3 ) 
      AoACols = find(cellfun('length',regexp(ChanName,'Alpha')) == 1);
      VIndxCols = find(cellfun('length',regexp(ChanName,'VIndx')) == 1);
      VIndyCols = find(cellfun('length',regexp(ChanName,'VIndy')) == 1);
      PhiCols = find(cellfun('length',regexp(ChanName,'Phi')) == 1);
      AxIndCols = find(cellfun('length',regexp(ChanName,'AxInd')) == 1);
      TnIndCols = find(cellfun('length',regexp(ChanName,'TnInd')) == 1);
      ClCols  = find(cellfun('length',regexp(ChanName,'Cl')) == 1);
      CdCols  = find(cellfun('length',regexp(ChanName,'Cd')) == 1);
      CmCols  = find(cellfun('length',regexp(ChanName,'Cm')) == 1);
      numNodes = length(AoACols); 
   end
         %numSteps = size(Channels,1);
         %chanPerNode = size(Channels,2)/numNodes;
         %AoAData = zeros(numSteps, numNodes, numCases);
         %titleText = DescStr;
   for iNode = 1:numNodes
      if ( iCase == 3 )
         
         ClFigs(iNode) = figure();
         grid;
         ylabel('Cl');    
         xlabel('AOA');
         title(['Node Output #' num2str(iNode)]);
         
         CdFigs(iNode) = figure();
         grid;
         ylabel('Cd');    
         xlabel('AOA');
         title(['Node Output #' num2str(iNode)]);
         
         CmFigs(iNode) = figure();
         grid;
         ylabel('Cm');    
         xlabel('AOA');
         title(['Node Output #' num2str(iNode)]);
         
      end
      figure(ClFigs(iNode));
      hold on;
      plot(Channels(2:end,AoACols(iNode)),Channels(2:end,ClCols(iNode)),'Color',LineColors{iCase},'LineWidth',2);
     % plot(Channels(500:end,ClCols(iNode)),'Color',LineColors{iCase},'LineWidth',2);
      figure(CdFigs(iNode));
      hold on;
      plot(Channels(2:end,AoACols(iNode)),Channels(2:end,CdCols(iNode)),'Color',LineColors{iCase},'LineWidth',2);
      figure(CmFigs(iNode));
      hold on;
      plot(Channels(2:end,AoACols(iNode)),Channels(2:end,CmCols(iNode)),'Color',LineColors{iCase},'LineWidth',2);
  
   end
end
 % Legends
 for iNode = 1:numNodes  
    figure(ClFigs(iNode));
    legend(legendStr);
    figure(CdFigs(iNode));
    legend(legendStr);
    figure(CmFigs(iNode));
    legend(legendStr);
 end
end

% for i=2:numNodes
%    figure;
%    plot(Channels(5:end,AoACols(i)),Channels(5:end,ClCols(i)),'r','LineWidth',2) 
%    % hold on;
%    % plot(aoa_local2(500:end),cl_local2(500:end),'m','LineWidth',2)
%    % plot(aoa_local3(500:end),cl_local3(500:end),'b','LineWidth',2)
%    grid;
%    ylabel('Cl');    
%    xlabel('AOA');
%    % plot(AoACols(500:end),ClCols(500:end),'r','LineWidth',2) 
%    % 
%    % legend(legendStr,'Location','northwest');   
%    title( ['Node' num2str(i) ': ' titleText])    
% end
% end