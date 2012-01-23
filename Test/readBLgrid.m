function [velocity, y, z, nz, ny, dz, dy, dt, zHub, z1, SummVars] = readBLgrid(FileName)
%[velocity, y, z, nz, ny, dz, dy, dt, zHub, z1, SummVars] = readBLgrid(FileName)
%
% Input:
% FileName       - string, containing file name to open (.wnd extension is optional)
%
% Output:
%  velocity      - 4-D vector: time, velocity component, iy, iz 
%  y             - 1-D vector: horizontal locations y(iy)
%  z             - 1-D vector: vertical locations z(iz)
%  nz, ny        - scalars: number of points in the vertical/horizontal
%                  direction of the grid
%  dz, dy, dt    - scalars: distance between two points in the vertical [m]/
%                  horizontal [m]/time [s] dimension
% zHub           - hub height [m]
% z1             - vertical location of bottom of grid [m above ground level]
% SumVars        - variables from the summary file (zHub, Clockwise, UBAR, TI_u, TI_v, TI_w)


len    = length(FileName);
ending = FileName(len-3:len);

if strcmpi( ending, '.wnd' )
    FileName = FileName(1:len-4);
end

%-------------------------------------------------------------

    % initialize variables
fileFmt  = 'int16';
ConvFact = 1.0; %results in meters and seconds

str      = {'HUB HEIGHT','CLOCKWISE','UBAR','TI(U','TI(V','TI(W'};  %MUST be in UPPER case
numVars  = length(str);   
SummVars = zeros(numVars, 1);

%% -----------------------------------------
%  READ THE HEADER OF THE BINARY FILE 
%  ----------------------------------------- 
fid_wnd   = fopen( [ FileName '.wnd' ] );
if ( fid_wnd <= 0 )
   error( 'Wind file could not be opened.' );
end
 
nffc  = fread( fid_wnd, 1, 'int16' );                     % number of components
      
if nffc ~= -99  % AN OLD-STYLE AERODYN WIND FILE
    dz      = fread( fid_wnd, 1, 'int16' );               % delta z in mm
    dy      = fread( fid_wnd, 1, 'int16' );               % delta y in mm
    dx      = fread( fid_wnd, 1, 'int16' );               % delta x (actually t in this case) in mm
    nt      = fread( fid_wnd, 1, 'int16' );               % half number of time steps
    MFFWS   = fread( fid_wnd, 1, 'int16' );               % 10 times mean FF wind speed, should be equal to MWS
              fread( fid_wnd, 5, 'int16' );               % unnecessary lines
    nz      = fread( fid_wnd, 1, 'int16' );               % 1000 times number of points in vertical direction, max 32
    ny      = fread( fid_wnd, 1, 'int16' );               % 1000 times the number of points in horizontal direction, max 32
              fread( fid_wnd, 3*(-nffc-1), 'int16' );

        % convert the integers to real numbers 
    nffc     = -nffc;
    dz       = 0.001*ConvFact*dz;
    dy       = 0.001*ConvFact*dy;
    dx       = 0.001*ConvFact*dx;
    MFFWS    = 0.1*ConvFact*MFFWS;
    nz       = fix( mod(nz,2^16) / 1000 );                % the mod 2^16 is a work around for somewhat larger grids
    ny       = fix( mod(ny,2^16) / 1000 );                % the mod 2^16 is a work around for somewhat larger grids
        
else %== -99, THE NEWER-STYLE AERODYN WIND FILE
    fc       = fread( fid_wnd, 1, 'int16' );              % should be 4 to allow turbulence intensity to be stored in the header

    if fc == 4
        nffc     = fread( fid_wnd, 1, 'int32' );              % number of components (should be 3)
        lat      = fread( fid_wnd, 1, 'float32' );            % latitude (deg)
        z0       = fread( fid_wnd, 1, 'float32' );            % Roughness length (m)
        zOffset  = fread( fid_wnd, 1, 'float32' );            % Reference height (m) = Z(1) + GridHeight / 2.0
        TI_U     = fread( fid_wnd, 1, 'float32' );            % Turbulence Intensity of u component (%)
        TI_V     = fread( fid_wnd, 1, 'float32' );            % Turbulence Intensity of v component (%)
        TI_W     = fread( fid_wnd, 1, 'float32' );            % Turbulence Intensity of w component (%)
    else
        if fc > 2
            nffc = 3;
        else
            nffc = 1;
        end
        TI_U  = 1;
        TI_V  = 1;
        TI_W  = 1;
        
        if fc == 8 %MANN model
            HeadRec = fread(fid_wnd,1,'int32');
            tmp     = fread(fid_wnd,1,'int32');  %nffc?
        end
        
    end %fc == 4
        
    dz       = fread( fid_wnd, 1, 'float32' );            % delta z in m 
    dy       = fread( fid_wnd, 1, 'float32' );            % delta y in m
    dx       = fread( fid_wnd, 1, 'float32' );            % delta x in m           
    nt       = fread( fid_wnd, 1, 'int32' );              % half the number of time steps
    MFFWS    = fread( fid_wnd, 1, 'float32');             % mean full-field wind speed

               fread( fid_wnd, 3, 'float32' );            % zLu, yLu, xLu: unused variables (for BLADED)
               fread( fid_wnd, 2, 'int32' );              % unused variables (for BLADED)
    nz       = fread( fid_wnd, 1, 'int32' );              % number of points in vertical direction
    ny       = fread( fid_wnd, 1, 'int32' );              % number of points in horizontal direction
               fread( fid_wnd, 3*(nffc-1), 'int32' );     % other length scales: unused variables (for BLADED)                

%     SummVars{numVars-3} = MFFWS;
%     SummVars{numVars-2} = TI_U;
%     SummVars{numVars-1} = TI_V;
%     SummVars{numVars}   = TI_W;
    SummVars(3:6) = [MFFWS, TI_U, TI_V, TI_W];
    if fc == 8        
        gamma  = fread(fid_wnd,1,'float32');               % MANN model shear parameter
        Scale  = fread(fid_wnd,1,'float32');               % MANN model scale length
        
    end
    
end % old or new bladed styles

nt     = max([nt*2,1]);
dt     = dx/MFFWS;
                
%% -----------------------------------------
%  READ THE SUMMARY FILE FOR SCALING FACTORS
%  -----------------------------------------                   
disp('Reading the summary file....');

indx     = SummVars;
fid_sum  = fopen( [ FileName '.sum' ] );

if ( fid_sum <= 0 )
    fclose(fid_sum);
    fclose(fid_wnd);
    
    error(['Could not open the summary file: ' FileName '.sum']);
    return;
end

while ( any( indx == 0 ) )  %MFFWS and the TIs should not be zero
    line  = fgetl(fid_sum);
    if ~ischar(line) 
        % We reached the end of the file

        fclose(fid_sum);
        fclose(fid_wnd);
        
        error('Reached the end of summary file without all necessary data.');            
    end

    line  = upper(line);
    findx = strfind(line,'=')+1;        %first index
    if isempty(findx)
        findx = 1;
    end    
    lindx = length(line);               %last  index

    i = 1;
    while i <= numVars
        if indx(i)==0

            k = findstr(line, str{i} );                
            if ~isempty(k)              % we found a string we're looking for                
                indx(i) = k;
                k=strfind(line,'%');
                if ~isempty(k)
                    lindx = max(findx,k-1);
                end
                
                tmp = strtok(line(findx:lindx));
                if isempty( str2num(tmp) ) %bjj: was str2double, but that returns NaN for 'T'/'F'
                    if strcmpi( tmp(1), 'T' )
                        SummVars(i) = 1;
                    else
                        SummVars(i) = -1;  %use this for false instead of zero.
                    end
                else
                    SummVars(i) = str2double( tmp );
                    break;
                end %if isempty( str2num(tmp) )
            end % ~isempty(k)
        end % indx(i)==0
        i = i + 1;
    end %while i
            
end %while any(indx==0)
   
% read the rest of the file to get the grid height offset, if it's there
ZGoffset = 0.0;

while ( true )
    line  = fgetl(fid_sum);

    if ~ischar(line)
        break;
    end

    line  = upper(line);
    findx = strfind(line,'HEIGHT OFFSET');

    if ~isempty(findx)
        lindx = length(line);
        findx = strfind(line,'=')+1;
        ZGoffset = str2double( strtok(line(findx:lindx))); %z grid offset
        break;
    end            
end  
fclose(fid_sum);      

%-----------------------------------------
%READ THE GRID DATA FROM THE BINARY FILE
%-----------------------------------------                   
disp('Reading and scaling the grid data...');

% nffc     = 3;
nv       = nffc*ny*nz;               % the size of one time step
Scale    = 0.00001*SummVars(3)*SummVars(4:6);
Offset   = [SummVars(3) 0 0];

velocity = zeros(nt,nffc,ny,nz);

%bjj: disp( [nv nffc ny nz nt] );

if SummVars(2) > 0 %clockwise rotation
    %flip the y direction....
    %let's change the dimension of velocity so that it's 4-d instead of 3-d   
    y_ix = ny:-1:1;
else
    y_ix = 1:ny;    
end

% [v cnt] = fread( fid_wnd, nv, fileFmt );
% if cnt < nv
%     error(['Could not read entire file: at grid record ' num2str( (it-1)*nv+cnt2 ) ' of ' num2str(nrecs)]);
% end
% disp('Scaling the grid data...');

for it = 1:nt
    
    [v cnt] = fread( fid_wnd, nv, fileFmt );
    if cnt < nv
        error(['Could not read entire file: at grid record ' num2str( (it-1)*nv+cnt ) ' of ' num2str(nv*nt)]);
    end
    
    cnt2 = 1;
    for iz = 1:nz
        for iy = y_ix
            for k=1:nffc
                velocity(it,k,iy,iz) = v(cnt2)*Scale(k) + Offset(k);
                cnt2 = cnt2 + 1;
            end %for k
        end %iy
    end % iz   
    
end %it

fclose(fid_wnd);


y    = (0:ny-1)*dy - dy*(ny-1)/2;

zHub = SummVars(1);
z1   = zHub - ZGoffset - dz*(nz-1)/2;  %this is the bottom of the grid
z    = (0:nz-1)*dz + z1;

disp('Finished.');
disp('');

return;