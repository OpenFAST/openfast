function [ADPar, ADBladeRootname] = newInputs_AD_v15(ADPar, AOC_AeroDyn_Test01)
%[ADPar] = newInputs_AD_v15(ADPar, ADrootname)
% ADPar is the data structure containing already-filled parameters for
%       AeroDyn, which will be modified for AeroDyn v15.
% ADrootname is the base file name (without path or extension) of the
%       AeroDyn input file [to set name of blade files]

%% ........................................................................
%% Environmental conditions:
Patm		  = 101325     ; %value from AD 14
Pvap	          = 2000        ; %value from AD 14
FluidDepth          = .6; %value from AD 14

%% ........................................................................
%% BEMT options:
[IndModel] = GetFASTPar(ADPar,'IndModel');        %Induction-factor model [NONE or WAKE or SWIRL] (unquoted string)
[InfModel] = GetFASTPar(ADPar,'InfModel');        %Inflow model [DYNIN or EQUIL] (unquoted string)
[TLModel ] = GetFASTPar(ADPar,'TLModel');         %Tip-loss model (EQUIL only) [PRANDtl, GTECH, or NONE] (unquoted string)
[HLModel ] = GetFASTPar(ADPar,'HLModel');         %Hub-loss model (EQUIL only) [PRANdtl or NONE] (unquoted string)
[AToler  ] = GetFASTPar(ADPar,'AToler');          %Induction-factor tolerance (convergence criteria) (-)

if (strcmpi(IndModel,'"none"'))
    WakeMod = 0;
    TanInd  = 'False'; % not used
else
    WakeMod = 1;    
    if strcmpi(InfModel,'"dynin"')
        disp('Warning: GDW not supported; will use BEMT instead.')        
    end

    if strcmpi( IndModel, '"SWIRL"' )
        TanInd  = 'true';
    else %strcmpi( IndModel, '"WAKE"' )
        TanInd  = 'False';
    end

end

SkewMod = 2;
if strcmpi(InfModel,'"equildb"')
    AIDrag = 'True';
    TIDrag = 'True';
elseif strcmpi(InfModel,'"equilda"')
    AIDrag = 'True';
    TIDrag = 'False';
elseif strcmpi(InfModel,'"equildt"')
    AIDrag = 'False';
    TIDrag = 'True';
else
    AIDrag = 'False';
    TIDrag = 'False';    
end
%IndToler = AToler/50;
IndToler = '"Default"';
MaxIter  = 100;
    
if strcmpi(TLModel,'"none"') 
    TipLoss = 'False';
else
    TipLoss = 'True';
    if strcmpi(TLModel,'"GTECH"') 
        disp('Warning: GTECH model not supported in AeroDyn 15. Using PRANDtl instead.')
    end
end

if strcmpi(HLModel,'"none"') 
    HubLoss = 'False';
else
    HubLoss = 'True';
end

%% ........................................................................
%% Unsteady Aero options:
[StallMod] = GetFASTPar(ADPar,'StallMod');        %Dynamic stall included [BEDDOES or STEADY] (unquoted string)
if strcmpi(StallMod,'"beddoes"')
    AFAeroMod = 2;
else %steady
    AFAeroMod = 1;
end

UAMod = 3; %2 is more stable at this point
FLookup = 'True'; 


%% ........................................................................
%% Airfoil Info:
[NumFoil ] = GetFASTPar(ADPar,'NumFoil');         %Number of airfoil files (-)

InCol_Alfa  = 1;
InCol_Cl    = 2;
InCol_Cd    = 3;
InCol_Cm    = 4;
InCol_Cpmin = 0;

NumAFfiles = NumFoil;
% AFNames = FoilNm ; I don't think this is necessary because I already
% modified the Matlab2FAST and FAST2Matlab code to do this.
% However, someone will have to change the airfoil files to the new format.
%--------------

%% ........................................................................
%% rotor/blade properties:
[UseCm   ] = GetFASTPar(ADPar,'UseCm');           %Use aerodynamic pitching moment model? [USE_CM or NO_CM] (unquoted string)
[BldNodes] = GetFASTPar(ADPar,'BldNodes');        %Number of blade nodes used for analysis (-)

if strcmpi(UseCm,'"USE_CM"')
    UseBlCm = 'True';
else
    UseBlCm = 'False';
end

% ADBlFile = '"AeroDyn_Blade.dat"';

% BldNodeHeaders = ADPar.BldNodesHdr;

NumBlNds = BldNodes + 2;
ADPar.BldNodesHdr = {'BlSpn','BlCrvAC','BlSwpAC','BlCrvAng','BlTwist','BlChord','BlAFID'};
BldNodeTable = ADPar.BldNodes; %old values; 
%we'll overwrite BldNodes table here, with new columns and 2 new nodes:

rHub = BldNodeTable{1,1}-0.5*BldNodeTable{1,3}; %RNodes(1) - 0.5*DRNodes(1)

% node at the blade root:
ADPar.BldNodes{1,1} = 0.0;                              %BlSpn
ADPar.BldNodes{1,2} = 0.0;                              %BlCrvAC
ADPar.BldNodes{1,3} = 0.0;                              %BlSwpAC
ADPar.BldNodes{1,4} = 0.0;                              %BlCrvAng
ADPar.BldNodes{1,5} = BldNodeTable{1,2};                %BlTwist
ADPar.BldNodes{1,6} = BldNodeTable{1,4};                %BlChord
ADPar.BldNodes{1,7} = BldNodeTable{1,5};                %BlAFID

for i=1:BldNodes
    ADPar.BldNodes{i+1,1} = BldNodeTable{i,1} - rHub;   %BlSpn
    ADPar.BldNodes{i+1,2} = 0.0;                        %BlCrvAC
    ADPar.BldNodes{i+1,3} = 0.0;                        %BlSwpAC
    ADPar.BldNodes{i+1,4} = 0.0;                        %BlCrvAng
    ADPar.BldNodes{i+1,5} = BldNodeTable{i,2};          %BlTwist
    ADPar.BldNodes{i+1,6} = BldNodeTable{i,4};          %BlChord
    ADPar.BldNodes{i+1,7} = BldNodeTable{i,5};          %BlAFID
end

% node at the blade tip:
ADPar.BldNodes{NumBlNds,1} = BldNodeTable{BldNodes,1} - rHub + 0.5*BldNodeTable{BldNodes,3}; %BlSpn = RNodes(end) + 0.5*DRNodes(end)
ADPar.BldNodes{NumBlNds,2} = 0.0;                       %BlCrvAC
ADPar.BldNodes{NumBlNds,3} = 0.0;                       %BlSwpAC
ADPar.BldNodes{NumBlNds,4} = 0.0;                       %BlCrvAng
ADPar.BldNodes{NumBlNds,5} = BldNodeTable{BldNodes,2};  %BlTwist
ADPar.BldNodes{NumBlNds,6} = BldNodeTable{BldNodes,4};  %BlChord
ADPar.BldNodes{NumBlNds,7} = BldNodeTable{BldNodes,5};  %BlAFID


%% ........................................................................
%% tower shadow/aero:
[TwrShad,  err1] = GetFASTPar(ADPar,'TwrShad');         %Tower-shadow velocity deficit (-)
if ~err1 && strcmpi(TwrShad, '"newtower"' )
% if TwrShad is "NEWTOWER" we can convert tower influence

    % TwrShadow exists already
    [TwrPotent] = GetFASTPar(ADPar,'TwrPotent');     %Calculate tower potential flow (flag)
    if strcmpi(TwrPotent,'false')
        TwrPotent = 0;
    else
        TwrPotent = 1; % ignore the tower file, so that section will need to be filled out        
    end
    [ADPar] = SetFASTPar(ADPar,'TwrPotent',TwrPotent);     %Calculate tower potential flow (flag)
    
    
    [CalcTwrAero] = GetFASTPar(ADPar,'CalcTwrAero'); %Calculate aerodynamic drag of the tower at the ElastoDyn nodes.
    TwrAero = CalcTwrAero;
    
elseif ~err1 && TwrShad > 0
    
    % create TwrDiam and TwrCd values to get tower shadow results
    % comparable to AD14.
   [ShadHWid] = GetFASTPar(ADPar,'ShadHWid');           
   [T_Shad_Refpt] = GetFASTPar(ADPar,'T_Shad_Refpt'); 

   uref = TwrShad;
   bref = ShadHWid;
   lref = T_Shad_Refpt;

   TwrDiam = 2 * bref.^2 / lref;  % note this isn't necessarially the tower diameter, but will give similar tower shadow to AD14
   TwrCd   = 2 * uref * bref ./ TwrDiam;

   % need TwrElev from ED
   %add this to the tower table!!! FIX ME                    
   fprintf( 'TwrDiam = %f\n', TwrDiam);
   fprintf( 'TwrCd   = %f\n', TwrCd);
   
    TwrAero = 'False';        

else    
    
% if it's not "NEWTOWER", we're going to ignore TwrShad, ShadHWid, and T_Shad_Refpt
    disp('Warning: AeroDyn v15 tower shadow model is not compatible with AeroDyn v14. Ignoring tower shadow.');
    TwrAero   = 'False';        
end
    
%% ........................................................................
%% add new variables to the ADPar structure:    


% table: RNodes, AeroTwst, DRNodes, Chord, NFoil,     
    %----------------------------------------------------------------------
    % Create new fields for AeroDyn v15.00.x:
    %----------------------------------------------------------------------      
    
    %ADBlFile(1), ADBlFile(2), ADBlFile(3)
        
    newVars = { 'Echo', 'WakeMod', 'AFAeroMod', 'TwrAero', ...
                'SpdSound', ...
                'SkewMod', 'TipLoss', 'HubLoss', 'TanInd', 'AIDrag', 'TIDrag', 'IndToler', 'MaxIter', ...
                'UAMod', 'FLookup', ...
                'InCol_Alfa', 'InCol_Cl', 'InCol_Cd', 'InCol_Cm', 'InCol_Cpmin', 'NumAFfiles', ...   
                'UseBlCm', 'NumBlNds', ...
                'SumPrint', 'NBlOuts', 'NTwOuts'};
                
    newVals = { 'False', WakeMod, AFAeroMod, TwrAero, ...
                SpdSound, ...
                SkewMod', TipLoss, HubLoss, TanInd, AIDrag, TIDrag, IndToler, MaxIter, ...
                UAMod, FLookup, ...
                InCol_Alfa, InCol_Cl, InCol_Cd, InCol_Cm, InCol_Cpmin, NumAFfiles, ...   
                UseBlCm, NumBlNds, ...
                'False',      0,       0       };
           
            
    n = length(ADPar.Label);
    for i = 1:length(newVars);      
        n = n + 1;
        ADPar.Label{n} = newVars{i};
        ADPar.Val{n}   = newVals{i};
    end 

    
    [~, err1] = GetFASTPar(ADPar,'TwrShadow');     %Calculate tower shadow (flag)
    if err1
        n = n + 1;
        ADPar.Label{n} = 'TwrShadow';
        ADPar.Val{n}   = 'False';
    end 
    
    [~, err1] = GetFASTPar(ADPar,'TwrPotent');     %Calculate tower shadow (flag)
    if err1
        n = n + 1;
        ADPar.Label{n} = 'TwrPotent';
        ADPar.Val{n}   = 0;
    end 
    
    ADBladeRootname = [ AOC_AeroDyn_Test01 '_blade.dat' ];
    for i=1:3
        n = n + 1;
        ADPar.Label{n} = ['ADBlFile(' num2str(i) ')'];
        ADPar.Val{n}   = [ '"' ADBladeRootname '"' ];
    end

return;
