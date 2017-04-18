
function PlotClCdCmCurves_DT(Index, baseName, baseFolder )

if (nargin < 2)
Index = 2;  % Valid values are 1, 2, 3
end 

if (nargin < 2)
   %baseName = 'AOC';
   baseName = 'NRELOffshrBsline5MW_Onshore';
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
%legendStr = {'0.000138 sec', '0.00138 sec', '0.0138 sec','0.1.38 sec'};
legendStr = {'0.00138 sec', '0.0138 sec','0.1.38 sec'};
   % Loop over different timestep sizes
for iCase=3:-1:1
   FileName = [baseFolder '\' baseName '_Test0' num2str(Index) '.' num2str(iCase) '.out'];
   [Channels, ChanName, ChanUnit,DescStr] = ReadFASTtext(FileName, delim, HeaderRows, NameLine, UnitsLine );
   
   if ( iCase == 3 ) 
      AoACols = find(cellfun('length',regexp(ChanName,'Alpha')) == 1);
      ClCols  = find(cellfun('length',regexp(ChanName,'Cl')) == 1);
      CdCols  = find(cellfun('length',regexp(ChanName,'Cd')) == 1);
      CmCols  = find(cellfun('length',regexp(ChanName,'Cm')) == 1);
      numNodes = length(AoACols); 
   end
         %numSteps = size(Channels,1);
         %chanPerNode = size(Channels,2)/numNodes;
         %AoAData = zeros(numSteps, numNodes, numCases);
         %titleText = DescStr;
   for iNode = 2:numNodes   
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
      plot(Channels(5:end,AoACols(iNode)),Channels(5:end,ClCols(iNode)),'Color',LineColors{iCase},'LineWidth',2);
      figure(CdFigs(iNode));
      hold on;
      plot(Channels(5:end,AoACols(iNode)),Channels(5:end,CdCols(iNode)),'Color',LineColors{iCase},'LineWidth',2);
      figure(CmFigs(iNode));
      hold on;
      plot(Channels(5:end,AoACols(iNode)),Channels(5:end,CmCols(iNode)),'Color',LineColors{iCase},'LineWidth',2);
  
   end
end
 % Legends
 for iNode = 2:numNodes  
    figure(ClFigs(iNode));
    legend(legendStr);
    figure(CdFigs(iNode));
    legend(legendStr);
    figure(CmFigs(iNode));
    legend(legendStr);
 end
end

