function [v, y, z, dt, ny, nz] = loadFFtxt(fname)
% load formatted full-field turbulence files   
% v(time,iy,iz)

%line  2   description
%line  4   header
%line  5   grid size/resolution data
%line  7   z coordinate description
%line  8   z coordinates
%line 10   y coordinate description
%line 11   y coordinates
%line 13   time, hub wind
%line 14ff grid wind speeds

fid = fopen(fname);
if fid > 0
   
    for i=1:5
        line = fgetl(fid);
    end
    tmp   = sscanf(line,'%f', 7);
    ny    = tmp(1);
    nz    = tmp(2);
    dy    = tmp(3);
    dz    = tmp(4);
    dt    = tmp(5);
    zhh   = tmp(6);
    mffws = tmp(7);
    
    for i=6:8
        line = fgetl(fid);
    end
    z = sscanf(line,'%f',nz);
    z = z + zhh;  %z coordinates are now height above ground, instead of relative to hub.

    for i=9:11
        line = fgetl(fid);
    end
    y = sscanf(line,'%f',ny);      
    
    nt        = 0;
    gridSize  = ny*nz;
    
    [tmp1, cnt] = fscanf(fid,'%f',2);
    while cnt == 2
        
        nt     = nt + 1;
        t(nt)  = tmp1(1);
        hh(nt) = tmp1(2);
        
        [tmp2, cnt] = fscanf(fid,'%f',[ny,nz]);               
        v(nt,:,:)   = fliplr(tmp2);
        
            % read the next grouping, if it exists
        [tmp1, cnt] = fscanf(fid,'%f',2);        
    end    %while
    
    fclose(fid);
    
else
    error([ 'Error opening ' fname '.']);
end    

return;