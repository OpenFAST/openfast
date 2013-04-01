function [LinResults]=get_Lin4GUI(FASTLinName)
% Function formed from Eigenanalysis to parse linearization results.
%
% Input: FASTLinName = fast lin output filename including its extension (e.g. '.lin')
% Output:
%     LinResults.DescStates=DescStates;
%     LinResults.DescCntrlInpt=DescCntrlInpt;
%     LinResults.DescDisturbnc=DescDisturbnc;
%     LinResults.DescOutput=DescOutput;
%     LinResults.OutName=OutName;
%     LinResults.Azimuth=Azimuth;
%     LinResults.Amat=AMat;
%     LinResults.Bmat=BMat;
%     LinResults.Bdmat=Bdmat;
%     LinResults.Cmat=CMat;
%     LinResults.Dmat=DMat;
%     LinResults.Ddmat=DdMat;
%     LinResults.NRotTripletStates=NRotTripletStates;
%     LinResults.RotTripletIndicesCntrlInpt=RotTripletIndicesCntrlInpt;
%     LinResults.NRotTripletCntrlInpt=NRotTripletCntrlInpt;
%     LinResults.RotTripletIndicesCntrlInpt=RotTripletIndicesCntrlInpt;
%     LinResults.NRotTripletOutput=NRotTripletOutput;
%     LinResults.RotTripletIndicesOutput=RotTripletIndicesOutput;
%
% If the linearization did not use ctrl inputs (dist, outputs, etc.)
% the corresponding entry in LinResults will still exist, but contain an
% empty member.

% Open the FAST linear file:
FileID = fopen( FASTLinName, 'r' );

% Abort if file is not found:
if ( FileID == -1 )
    error( ['FAST linearization output file "', FASTLinName, '" not found.  Aborting.'] );
    clear ClearRootName RootName FASTLinName FileID;
    return;
end

% Read in (and ignore) the first portion of the file:
for Row = 1:10
    line   = fgetl( FileID );
end

% Read in the azimuth-average rotor speed, the steady state period of solution,
%   then ignore the next three lines:
line = fgetl( FileID );
RotSpeed   = str2num( line(55:68) );    % in (rad/s)
line = fgetl( FileID );
if ( RotSpeed > 0 )
    Period = str2num( line(55:68) );    % in (sec)
else
    clear Period;
end
line = fgetl( FileID );
line = fgetl( FileID );
line = fgetl( FileID );

% Read in the number of equally spaced azimuth steps, the model order, the
%   number of active (enabled) DOFs, the number of control inputs, the
%   number of wind input disturbances, and the number of output measurements,
%   then ignore the next four lines:
line = fgetl( FileID );
NAzimStep  = str2num( line(55:58) );
line = fgetl( FileID );
MdlOrder   = str2num( line(55:58) );
line = fgetl( FileID );
NActvDOF   = str2num( line(55:58) );
line = fgetl( FileID );
NInputs    = str2num( line(55:58) );
line = fgetl( FileID );
NDisturbs  = str2num( line(55:58) );
line = fgetl( FileID );
NumOuts    = str2num( line(55:58) );
line = fgetl( FileID );
line = fgetl( FileID );
line = fgetl( FileID );
line = fgetl( FileID );

% Define the number of states:
N = 2*NActvDOF;                     % the number of states in the linearized model {x}

% Initialize the matrices before reading them in:
DescStates        = cell (         NActvDOF ,1        );
Azimuth           = zeros(         NAzimStep,1        );
xdop              = zeros(N                 ,NAzimStep);
xop               = zeros(N                 ,NAzimStep);
AMat              = zeros(N       ,N        ,NAzimStep);
    BMat          = zeros(N       ,NInputs  ,NAzimStep);
    BdMat         = zeros(N       ,NDisturbs,NAzimStep);
    CMat          = zeros(NumOuts ,N        ,NAzimStep);
    DMat      = zeros(NumOuts ,NInputs  ,NAzimStep);
    DdMat     = zeros(NumOuts ,NDisturbs,NAzimStep);
if ( NInputs > 0 )
    DescCntrlInpt = cell (         NInputs  ,1        );
else
    DescCntrlInpt = {};
end
if ( NDisturbs > 0 )
    DescDisturbnc = cell (         NDisturbs,1        );
else
    DescDisturbnc = {};
end
if ( NumOuts > 0 )
    DescOutput    = cell (NumOuts           ,1        );
    OutName       = cell (NumOuts           ,1        );
    yop           = zeros(NumOuts           ,NAzimStep);
else
    DescOutput    = {};
    OutName       = {};
    yop           = zeros(NumOuts           ,NAzimStep);
end
if ( MdlOrder == 2 )
    qd2op         = zeros(NActvDOF          ,NAzimStep);
    qdop          = zeros(NActvDOF          ,NAzimStep);
    qop           = zeros(NActvDOF          ,NAzimStep);
    MassMat       = zeros(NActvDOF,NActvDOF ,NAzimStep);
    DampMat       = zeros(NActvDOF,NActvDOF ,NAzimStep);
    StffMat       = zeros(NActvDOF,NActvDOF ,NAzimStep);
    if ( NInputs > 0 )
        FMat      = zeros(NActvDOF,NInputs  ,NAzimStep);
    end
    if ( NDisturbs > 0 )
        FdMat     = zeros(NActvDOF,NDisturbs,NAzimStep);
    end
    if ( NumOuts > 0 )
        VelCMat   = zeros(NumOuts ,NActvDOF ,NAzimStep);
        DspCMat   = zeros(NumOuts ,NActvDOF ,NAzimStep);
    end
end

% Read in the description of the model states, control inputs, input wind
%   disturbances, and output measurements:
if ( MdlOrder == 1 )
    for Row = 1:NActvDOF
        DescStates   {Row} = strtrim( fgetl( FileID ) );
    end
    line = fgetl( FileID );
else    % ( MdlOrder == 2 )
    for Row = 1:NActvDOF
        DescStates   {Row} = strtrim( fgetl( FileID ) );
    end
end
line = fgetl( FileID );
line = fgetl( FileID );
line = fgetl( FileID );
line = fgetl( FileID );
if ( NInputs > 0 )
    for Row = 1:NInputs
        DescCntrlInpt{Row} = strtrim( fgetl( FileID ) );
    end
else
    line = fgetl( FileID );
end
line = fgetl( FileID );
line = fgetl( FileID );
line = fgetl( FileID );
line = fgetl( FileID );
if ( NDisturbs > 0 )
    for Row = 1:NDisturbs
        DescDisturbnc{Row} = strtrim( fgetl( FileID ) );
    end
else
    line = fgetl( FileID );
end
line = fgetl( FileID );
line = fgetl( FileID );
line = fgetl( FileID );
line = fgetl( FileID );
if ( NumOuts > 0 )
    for Row = 1:NumOuts
        DescOutput   {Row} = strtrim( fgetl( FileID ) );
        EqlCol             = min(strfind(DescOutput{Row},'=')) + 1;               % find the index following the index of the first instance of string '='
        OutName      {Row} = strtok(DescOutput{Row}(EqlCol:length(DescOutput{Row})));
    end
else
    line = fgetl( FileID );
end
line = fgetl( FileID );
line = fgetl( FileID );
line = fgetl( FileID );
line = fgetl( FileID );

% Loop through all azimuth steps: 
for L = 1:NAzimStep
    
% Read in the azimuth step and reference azimuth angle from the azimuth title,
%   then ignore the two column headers:
    line = fgetl( FileID );
    Azimuth(L) = str2num( line(21:26) );   % in (deg)
    AzimB1Up   = str2num( line(60:65) );   % in (deg)
    line = fgetl( FileID );
    line = fgetl( FileID );
    
% Read in the vectors and matrices:
    if ( MdlOrder == 1 )

        for Row = 1:N
            line = fgetl( FileID );
        
            xdop(Row,L) = str2num( line( 1:10) );
            xop (Row,L) = str2num( line(14:23) );
        
            for Col = 1:N
                AMat         (Row,Col,L) = str2num( line(( 27     + 11*(               Col - 1 ) ):( 36     + 11*(               Col - 1 ) )) );
            end    
            for Col = 1:NInputs
                BMat         (Row,Col,L) = str2num( line(( 27 + 2 + 11*(           N + Col - 1 ) ):( 36 + 2 + 11*(           N + Col - 1 ) )) ); 
            end
            if ( NInputs > 0 )      % control inputs exist
                for Col = 1:NDisturbs
                    BdMat    (Row,Col,L) = str2num( line(( 27 + 4 + 11*( NInputs + N + Col - 1 ) ):( 36 + 4 + 11*( NInputs + N + Col - 1 ) )) ); 
                end
            else                    % no control inputs selected
                for Col = 1:NDisturbs
                    BdMat    (Row,Col,L) = str2num( line(( 27 + 2 + 11*(           N + Col - 1 ) ):( 36 + 2 + 11*(           N + Col - 1 ) )) ); 
                end   
            end
        end

        if ( NumOuts > 0 )
            line = fgetl( FileID ); % ignore the header
            line = fgetl( FileID ); % for the output
            line = fgetl( FileID ); % matrices
  
            for Row = 1:NumOuts
                line = fgetl( FileID );
        
                yop (Row,L) = str2num( line( 1:10) );
        
                for Col = 1:N
                    CMat     (Row,Col,L) = str2num( line(( 27     + 11*(               Col - 1 ) ):( 36     + 11*(               Col - 1 ) )) );
                end    
                for Col = 1:NInputs
                    DMat     (Row,Col,L) = str2num( line(( 27 + 2 + 11*(           N + Col - 1 ) ):( 36 + 2 + 11*(           N + Col - 1 ) )) ); 
                end   
                if ( NInputs > 0 )      % control inputs exist
                    for Col = 1:NDisturbs
                        DdMat(Row,Col,L) = str2num( line(( 27 + 4 + 11*( NInputs + N + Col - 1 ) ):( 36 + 4 + 11*( NInputs + N + Col - 1 ) )) ); 
                    end
                else                    % no control inputs selected
                    for Col = 1:NDisturbs
                        DdMat(Row,Col,L) = str2num( line(( 27 + 2 + 11*(           N + Col - 1 ) ):( 36 + 2 + 11*(           N + Col - 1 ) )) ); 
                    end   
                end
            end
        end
                
    else    % ( MdlOrder == 2 )
        
        for Row = 1:NActvDOF
            line = fgetl( FileID );
        
            qd2op(Row,L) = str2num( line( 1:10) );
            qdop (Row,L) = str2num( line(14:23) );
            qop  (Row,L) = str2num( line(27:36) );
        
            for Col = 1:NActvDOF
                MassMat      (Row,Col,L) = str2num( line(( 40     + 11*(                        Col - 1 ) ):( 49     + 11*(                        Col - 1 ) )) );
                DampMat      (Row,Col,L) = str2num( line(( 40 + 2 + 11*(             NActvDOF + Col - 1 ) ):( 49 + 2 + 11*(             NActvDOF + Col - 1 ) )) );
                StffMat      (Row,Col,L) = str2num( line(( 40 + 4 + 11*(           2*NActvDOF + Col - 1 ) ):( 49 + 4 + 11*(           2*NActvDOF + Col - 1 ) )) );   
            end    
            for Col = 1:NInputs
                FMat         (Row,Col,L) = str2num( line(( 40 + 6 + 11*(           3*NActvDOF + Col - 1 ) ):( 49 + 6 + 11*(           3*NActvDOF + Col - 1 ) )) );
            end    
            if ( NInputs > 0 )      % control inputs exist
                for Col = 1:NDisturbs
                    FdMat    (Row,Col,L) = str2num( line(( 40 + 8 + 11*( NInputs + 3*NActvDOF + Col - 1 ) ):( 49 + 8 + 11*( NInputs + 3*NActvDOF + Col - 1 ) )) ); 
                end
            else                    % no control inputs selected
                for Col = 1:NDisturbs
                    FdMat    (Row,Col,L) = str2num( line(( 40 + 6 + 11*(           3*NActvDOF + Col - 1 ) ):( 49 + 6 + 11*(           3*NActvDOF + Col - 1 ) )) ); 
                end   
            end
        end
        
        if ( NumOuts > 0 )
            line = fgetl( FileID ); % ignore the header
            line = fgetl( FileID ); % for the output
            line = fgetl( FileID ); % matrices
  
            for Row = 1:NumOuts
                line = fgetl( FileID );
        
                yop (Row,L) = str2num( line( 1:10) );
        
                for Col = 1:NActvDOF
                    VelCMat  (Row,Col,L) = str2num( line(( 40 + 2 + 11*(             NActvDOF + Col - 1 ) ):( 49 + 2 + 11*(             NActvDOF + Col - 1 ) )) );
                    DspCMat  (Row,Col,L) = str2num( line(( 40 + 4 + 11*(           2*NActvDOF + Col - 1 ) ):( 49 + 4 + 11*(           2*NActvDOF + Col - 1 ) )) );   
                end    
                for Col = 1:NInputs
                    DMat     (Row,Col,L) = str2num( line(( 40 + 6 + 11*(           3*NActvDOF + Col - 1 ) ):( 49 + 6 + 11*(           3*NActvDOF + Col - 1 ) )) );
                end    
                if ( NInputs > 0 )      % control inputs exist
                    for Col = 1:NDisturbs
                        DdMat(Row,Col,L) = str2num( line(( 40 + 8 + 11*( NInputs + 3*NActvDOF + Col - 1 ) ):( 49 + 8 + 11*( NInputs + 3*NActvDOF + Col - 1 ) )) ); 
                    end
                else                    % no control inputs selected
                    for Col = 1:NDisturbs
                        DdMat(Row,Col,L) = str2num( line(( 40 + 6 + 11*(           3*NActvDOF + Col - 1 ) ):( 49 + 6 + 11*(           3*NActvDOF + Col - 1 ) )) ); 
                    end   
                end  
            end
        end

    end
    
% Ignore the blank lines:
    line = fgetl( FileID );
    line = fgetl( FileID );
  
end

% Close the FAST linear file:
fclose( FileID );


% ----------- Assemble a 1st order model from the 2nd order model ---------
% ----------- Azimuth average the matrices --------------------------------

if ( MdlOrder == 2 )

% Form the [A], [B], and [Bd] matrices (so that {xdot} = [A]{x} + [B]{u} +
%   [Bd]{ud}) at each azimuth step using the 2nd order state matrices:
% Also, form the 1st order [C] matrix using [VelC] and [DspC]:
    for L = 1:NAzimStep
        xdop     ( 1:N                                   , L ) = [ qdop( 1:NActvDOF, L ); qd2op( 1:NActvDOF, L ) ];
        xop      ( 1:N                                   , L ) = [ qop(  1:NActvDOF, L ); qdop(  1:NActvDOF, L ) ];
        AMat     ( 1:NActvDOF        , ( NActvDOF + 1 ):N, L ) =  eye(NActvDOF);
        AMat     ( ( NActvDOF + 1 ):N, 1:NActvDOF        , L ) = -inv(MassMat(:,:,L))*StffMat(:,:,L);
        AMat     ( ( NActvDOF + 1 ):N, ( NActvDOF + 1 ):N, L ) = -inv(MassMat(:,:,L))*DampMat(:,:,L);
        if ( NInputs > 0 )
            BMat ( ( NActvDOF + 1 ):N, 1:NInputs         , L ) =  inv(MassMat(:,:,L))*FMat   (:,:,L);
        end
        if ( NDisturbs > 0 )
            BdMat( ( NActvDOF + 1 ):N, 1:NDisturbs       , L ) =  inv(MassMat(:,:,L))*FdMat  (:,:,L);
        end
        if ( NumOuts > 0 )
            CMat ( 1:NumOuts         , 1:N               , L ) = [ DspCMat(1:NumOuts,1:NActvDOF,L),  VelCMat(1:NumOuts,1:NActvDOF,L) ];   
        end 
    end
end;

% ----------- Find multi-blade coordinate (MBC) transformation indices ----
% Find the number of, and indices for, state triplets in the rotating
%   frame:
NRotTripletStates = 0;                  % first initialize to zero
if ( exist('RotTripletIndicesStates') )
    clear RotTripletIndicesStates;      % start from scratch if necessary
end
for I = 1:NActvDOF	% loop through all active (enabled) DOFs
    BldCol = strfind(DescStates{I},'blade');                    % find the starting index of the string 'blade'
    if ( ~isempty(BldCol) )     % true if the DescStates{I} contains the string 'blade'
        EqlCol = min(strfind(DescStates{I},'=')) + 1;           % find the index following the index of the first instance of string '='
        [ K, OK ]     = str2num(DescStates{I}(BldCol+6));       % save the blade number  for the initial   blade
        if ( OK && isreal(K) )  % true if the 2nd character after the string 'blade' is a real number (i.e., a blade number)
            Tmp              = zeros(1,3);                      % first initialize to zero
            Tmp(        1,K) = I;                               % save the index         for the initial   blade
            for J = (I+1):NActvDOF  % loop through all remaining active (enabled) DOFs
                if ( ( length(DescStates{J}) >= BldCol ) && ( strcmpi(DescStates{J}(EqlCol:BldCol),DescStates{I}(EqlCol:BldCol)) ) )    % true if we have the same state from a different blade
                    K = str2num(DescStates{J}(BldCol+6));       % save the blade numbers for the remaining blades
                    Tmp(1,K) = J;                               % save the indices       for the remaining blades
                end
            end	                    % J - all remaining active (enabled) DOFs
            if ( all(Tmp) )     % true if all the elements of Tmp are nonzero; thus, we found a triplet of rotating indices
                NRotTripletStates                            = NRotTripletStates + 1;   % this  is  the number  of  state triplets in the rotating frame
                RotTripletIndicesStates(NRotTripletStates,:) = Tmp;                     % these are the indices for state triplets in the rotating frame
            end
        end
    end
end                 % I - all active (enabled) DOFs

% Find the number of, and indices for, control input triplets in the
%   rotating frame:
NRotTripletCntrlInpt = 0;               % first initialize to zero
for I = 1:NInputs	% loop through all control inputs
    BldCol = strfind(DescCntrlInpt{I},'blade');                 % find the starting index of the string 'blade'
    if ( ~isempty(BldCol) )     % true if the DescCntrlInpt{I} contains the string 'blade'
        EqlCol = min(strfind(DescCntrlInpt{I},'=')) + 1;        % find the index following the index of the first instance of string '='
        [ K, OK ]     = str2num(DescCntrlInpt{I}(BldCol+6));    % save the blade number  for the initial   blade
        if ( OK && isreal(K) )  % true if the 2nd character after the string 'blade' is a real number (i.e., a blade number)
            Tmp              = zeros(1,3);                      % first initialize to zero
            Tmp(        1,K) = I;                               % save the index         for the initial   blade
            for J = (I+1):NInputs	% loop through all remaining control inputs
                if ( ( length(DescCntrlInpt{J}) >= BldCol ) && ( strcmpi(DescCntrlInpt{J}(EqlCol:BldCol),DescCntrlInpt{I}(EqlCol:BldCol)) ) )	% true if we have the same control input from a different blade
                    K = str2num(DescCntrlInpt{J}(BldCol+6));    % save the blade numbers for the remaining blades
                    Tmp(1,K) = J;                               % save the indices       for the remaining blades
                end
            end                     % J - all remaining active control inputs
            if ( all(Tmp) )     % true if all the elements of Tmp are nonzero; thus, we found a triplet of rotating indices
                NRotTripletCntrlInpt                               = NRotTripletCntrlInpt + 1;  % this  is  the number  of  control input triplets in the rotating frame
                RotTripletIndicesCntrlInpt(NRotTripletCntrlInpt,:) = Tmp;                       % these are the indices for control input triplets in the rotating frame
            end
        end
    end
end                 % I - all control inputs

% Find the number of, and indices for, output measurement triplets in the
%   rotating frame:
NRotTripletOutput = 0;                  % first initialize to zero
if ( exist('RotTripletIndicesOutput') )
    clear RotTripletIndicesOutput;      % start from scratch if necessary
end
for I = 1:NumOuts	% loop through all output measurements
    [ K, OK ]     = str2num(OutName{I}(length(OutName{I})));	% save the blade number  for the initial   blade
    if ( OK && isreal(K) )  % true if the last character of string OutName{I}' is a real number (i.e., a blade number)
        Tmp              = zeros(1,3);                          % first initialize to zero
        Tmp(        1,K) = I;                                   % save the index         for the initial   blade
        for J = (I+1):NumOuts	% loop through all remaining output measurements
            if ( ( length(OutName{J}) == length(OutName{I}) ) && ( strcmpi(OutName{J}(1:length(OutName{J})-1),OutName{I}(1:length(OutName{I})-1)) ) )	% true if we have the same output measurement from a different blade
                K = str2num(OutName{J}(length(OutName{J})));    % save the blade numbers for the remaining blades
                Tmp(1,K) = J;                                   % save the indices       for the remaining blades
            end
        end                     % J - all remaining active output measurements
        if ( all(Tmp) )     % true if all the elements of Tmp are nonzero; thus, we found a triplet of rotating indices
            NRotTripletOutput                            = NRotTripletOutput + 1;	% this  is  the number  of  output measurement triplets in the rotating frame
            RotTripletIndicesOutput(NRotTripletOutput,:) = Tmp;                     % these are the indices for output measurement triplets in the rotating frame
        end
    end
end                 % I - all output measurements

%Pack results into a data structure
LinResults.DescStates=DescStates;
LinResults.DescCntrlInpt=DescCntrlInpt;
LinResults.DescDisturbnc=DescDisturbnc;
LinResults.DescOutput=DescOutput;
LinResults.OutName=OutName;
LinResults.Azimuth=Azimuth;
LinResults.Amat=AMat;
LinResults.Bmat=BMat;
LinResults.Bdmat=BdMat;
LinResults.Cmat=CMat;
LinResults.Dmat=DMat;
LinResults.Ddmat=DdMat;
LinResults.NRotTripletStates=NRotTripletStates;
LinResults.RotTripletIndicesCntrlInpt=RotTripletIndicesCntrlInpt;
LinResults.NRotTripletCntrlInpt=NRotTripletCntrlInpt;
LinResults.RotTripletIndicesCntrlInpt=RotTripletIndicesCntrlInpt;
LinResults.NRotTripletOutput=NRotTripletOutput;
LinResults.RotTripletIndicesOutput=RotTripletIndicesOutput;



