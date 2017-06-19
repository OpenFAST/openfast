function writeFARMOutputs()

fid=fopen('OutListParameters.csv','wt');
x = {'Category','Name', 'Other Name(s)', 'Description', 'Convention', 'Units','Invalid Channel Criteria' };
csvFun = @(str)sprintf('%s,',str);
xchar = cellfun(csvFun, x, 'UniformOutput', false);
xchar = strcat(xchar{:});
xchar = strcat(xchar(1:end-1),'\n');
fprintf(fid,xchar);


% Super Controller

   % Global (turbine-independent) super controller input
fprintf(fid,'Global Super Controller Input\n');
for beta = 1:9
   x = {'',['SCGblIn' num2str(beta,'%1d')], ' ', ['Global (turbine-independent) super controller input ' num2str(beta,'%1d')], ' ', '(user)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);
end
  
   % Turbine-dependent super controller input ? for turbine ?
fprintf(fid,'Turbine-dependent Super Controller Input\n');
for alpha = 1:9
   for beta = 1:9
      x = {'',['SCT' num2str(alpha,'%1d') 'In' num2str(beta,'%1d')], ' ', ['Turbine-dependent super controller input ' num2str(beta,'%1d') ' for turbine ' num2str(alpha,'%1d')], ' ', '(user)' };
      csvFun = @(str)sprintf('%s,',str);
      xchar = cellfun(csvFun, x, 'UniformOutput', false);
      xchar = strcat(xchar{:});
      xchar = strcat(xchar(1:end-1),'\n');
      fprintf(fid,xchar);
   end
end

   % Global (turbine-independent) super controller output
fprintf(fid,'Global Super Controller Output\n');
for beta = 1:9
   x = {'',['SCGblOt' num2str(beta,'%1d')], ' ', ['Global (turbine-independent) super controller output ' num2str(beta,'%1d')], ' ', '(user)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);
end

   % Turbine-dependent super controller output ? for turbine ?
fprintf(fid,'Turbine-dependent Super Controller Output\n');
for alpha = 1:9
   for beta = 1:9
      x = {'',['SCT' num2str(alpha,'%1d') 'Ot' num2str(beta,'%1d')], ' ', ['Turbine-dependent super controller output ' num2str(beta,'%1d') ' for turbine ' num2str(alpha,'%1d')], ' ', '(user)' };
      csvFun = @(str)sprintf('%s,',str);
      xchar = cellfun(csvFun, x, 'UniformOutput', false);
      xchar = strcat(xchar{:});
      xchar = strcat(xchar(1:end-1),'\n');
      fprintf(fid,xchar);
   end
end

% Wind Turbine and its Inflow
fprintf(fid,'Rotor Centerline Orientation\n');
   % Orientation of the rotor centerline for turbine ? in the global coordinate system
for alpha = 1:9
   x = {'',['RtAxsXT' num2str(alpha,'%1d')], ' ', ['X-component of the rotor centerline orientation for turbine ' num2str(alpha,'%1d') ' in the global coordinate system'], ' ', '(-)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);
end

for alpha = 1:9
   x = {'',['RtAxsYT' num2str(alpha,'%1d')], ' ', ['Y-component of the rotor centerline orientation for turbine ' num2str(alpha,'%1d') ' in the global coordinate system'], ' ', '(-)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);
end

for alpha = 1:9
   x = {'',['RtAxsZT' num2str(alpha,'%1d')], ' ', ['Z-component of the rotor centerline orientation for turbine ' num2str(alpha,'%1d') ' in the global coordinate system'], ' ', '(-)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);
end

   % Position of the rotor (hub) center for turbine ? in the global coordinate system
fprintf(fid,'Position of the Rotor (Hub) Center\n');
for alpha = 1:9
   x = {'',['RtPosXT' num2str(alpha,'%1d')], ' ', ['X-component of the position of the rotor (hub) center for turbine ' num2str(alpha,'%1d') ' in the global coordinate system'], ' ', '(m)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);
end

for alpha = 1:9
   x = {'',['RtPosYT' num2str(alpha,'%1d')], ' ', ['Y-component of the position of the rotor (hub) center for turbine ' num2str(alpha,'%1d') ' in the global coordinate system'], ' ', '(m)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);
end

for alpha = 1:9
   x = {'',['RtPosZT' num2str(alpha,'%1d')], ' ', ['Z-component of the position of the rotor (hub) center for turbine ' num2str(alpha,'%1d') ' in the global coordinate system'], ' ', '(m)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);
end

   % Rotor diameter for turbine ?
fprintf(fid,'Rotor Diamete\n');
for alpha = 1:9
   x = {'',['RtDiamT' num2str(alpha,'%1d')], ' ', ['Rotor diameter for turbine ' num2str(alpha,'%1d')], ' ', '(m)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);
end


   % Nacelle-yaw error for turbine ?
fprintf(fid,'Nacelle-Yaw Error\n');
for alpha = 1:9
   x = {'',['YawErrT' num2str(alpha,'%1d')], ' ', ['Nacelle-yaw error for turbine ' num2str(alpha,'%1d')], ' ', '(deg)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);
end

   % Ambient turbulence intensity  of the wind at the rotor disk for turbine ?
fprintf(fid,'Ambient Turbulence Intensity of the Wind at the Rotor Disk\n');
for alpha = 1:9
   x = {'',['TIAmbT' num2str(alpha,'%1d')], ' ', ['Ambient turbulence intensity of the wind at the rotor disk for turbine ' num2str(alpha,'%1d')], ' ', '(percent)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);
end

   % Rotor-disk-averaged ambient wind speed (normal to disk, not including structural motion, local induction or wakes from upstream turbines) for turbine ?
fprintf(fid,'Rotor-Disk-Averaged Ambient Wind Speed\n');
for alpha = 1:9
   x = {'',['RtVAmbT' num2str(alpha,'%1d')], ' ', ['Rotor-disk-averaged ambient wind speed (normal to disk: not including structural motion: local induction or wakes from upstream turbines) for turbine ' num2str(alpha,'%1d')], ' ', '(m/s)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);
end

   % Rotor-disk-averaged relative wind speed (normal to disk, including structural motion and wakes from upstream turbines, but not including local induction) for turbine ?
fprintf(fid,'Rotor-Disk-Averaged Relative Wind Speed\n');
for alpha = 1:9
   x = {'',['RtVRelT' num2str(alpha,'%1d')], ' ', ['Rotor-disk-averaged relative wind speed (normal to disk: including structural motion and wakes from upstream turbines: but not including local induction) for turbine ' num2str(alpha,'%1d')], ' ', '(m/s)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);
end

   % Azimuthally averaged thrust force coefficient (normal to disk) for radial node ? of turbine ?
fprintf(fid,'Azimuthally Averaged Thrust Force Coefficient\n');
for alpha = 1:9
   for beta = 1:20
      x = {'',['CtT' num2str(alpha,'%1d') 'N' num2str(beta,'%02d')], ' ', ['Azimuthally averaged thrust force coefficient (normal to disk) for radial node ' num2str(beta,'%02d') ' of turbine ' num2str(alpha,'%1d')], ' ', '(-)' };
      csvFun = @(str)sprintf('%s,',str);
      xchar = cellfun(csvFun, x, 'UniformOutput', false);
      xchar = strcat(xchar{:});
      xchar = strcat(xchar(1:end-1),'\n');
      fprintf(fid,xchar);
   end
end

% Wake (for an Individual Rotor)

   % Orientation of the wake centerline for downstream distance ? of turbine ? in the global coordinate system
fprintf(fid,'Wake Centerline Orientation\n');
for alpha = 1:9
   for gamma = 1:9
      x = {'',['WkAxsXT' num2str(alpha,'%1d') 'D' num2str(gamma,'%1d')], ' ', ['X-component of the wake centerline orientation for downstream distance ' num2str(gamma,'%1d') ' of turbine ' num2str(alpha,'%1d') ' in the global coordinate system'], ' ', '(-)' };
      csvFun = @(str)sprintf('%s,',str);
      xchar = cellfun(csvFun, x, 'UniformOutput', false);
      xchar = strcat(xchar{:});
      xchar = strcat(xchar(1:end-1),'\n');
      fprintf(fid,xchar);
   end
end  
for alpha = 1:9
   for gamma = 1:9
      x = {'',['WkAxsYT' num2str(alpha,'%1d') 'D' num2str(gamma,'%1d')], ' ', ['Y-component of the wake centerline orientation for downstream distance ' num2str(gamma,'%1d') ' of turbine ' num2str(alpha,'%1d') ' in the global coordinate system'], ' ', '(-)' };
      csvFun = @(str)sprintf('%s,',str);
      xchar = cellfun(csvFun, x, 'UniformOutput', false);
      xchar = strcat(xchar{:});
      xchar = strcat(xchar(1:end-1),'\n');
      fprintf(fid,xchar);
   end
 end  
for alpha = 1:9
   for gamma = 1:9
      x = {'',['WkAxsZT' num2str(alpha,'%1d') 'D' num2str(gamma,'%1d')], ' ', ['Z-component of the wake centerline orientation for downstream distance ' num2str(gamma,'%1d') ' of turbine ' num2str(alpha,'%1d') ' in the global coordinate system'], ' ', '(-)' };
      csvFun = @(str)sprintf('%s,',str);
      xchar = cellfun(csvFun, x, 'UniformOutput', false);
      xchar = strcat(xchar{:});
      xchar = strcat(xchar(1:end-1),'\n');
      fprintf(fid,xchar);
   end
end  
   
   % Center position of the wake centerline for downstream distance ? of turbine ? in the global coordinate system
fprintf(fid,'Center Position of Wake Centerline\n');
for alpha = 1:9
   for gamma = 1:9
      x = {'',['WkPosXT' num2str(alpha,'%1d') 'D' num2str(gamma,'%1d')], ' ', ['X-component of the center position of the wake centerline for downstream distance ' num2str(gamma,'%1d') ' of turbine ' num2str(alpha,'%1d') ' in the global coordinate system'], ' ', '(m)' };
      csvFun = @(str)sprintf('%s,',str);
      xchar = cellfun(csvFun, x, 'UniformOutput', false);
      xchar = strcat(xchar{:});
      xchar = strcat(xchar(1:end-1),'\n');
      fprintf(fid,xchar);
   end
end  
for alpha = 1:9
   for gamma = 1:9
      x = {'',['WkPosYT' num2str(alpha,'%1d') 'D' num2str(gamma,'%1d')], ' ', ['Y-component of the center position of the wake centerline for downstream distance ' num2str(gamma,'%1d') ' of turbine ' num2str(alpha,'%1d') ' in the global coordinate system'], ' ', '(m)' };
      csvFun = @(str)sprintf('%s,',str);
      xchar = cellfun(csvFun, x, 'UniformOutput', false);
      xchar = strcat(xchar{:});
      xchar = strcat(xchar(1:end-1),'\n');
      fprintf(fid,xchar);
   end
 end  
for alpha = 1:9
   for gamma = 1:9
      x = {'',['WkPosZT' num2str(alpha,'%1d') 'D' num2str(gamma,'%1d')], ' ', ['Z-component of the center position of the wake centerline for downstream distance ' num2str(gamma,'%1d') ' of turbine ' num2str(alpha,'%1d') ' in the global coordinate system'], ' ', '(m)' };
      csvFun = @(str)sprintf('%s,',str);
      xchar = cellfun(csvFun, x, 'UniformOutput', false);
      xchar = strcat(xchar{:});
      xchar = strcat(xchar(1:end-1),'\n');
      fprintf(fid,xchar);
   end
end  
   
   % Advection, deflection, and meandering velocity (not including the horizontal wake-deflection correction) of the wake for downstream distance ? of turbine ? in the global coordinate system
fprintf(fid,'Advection: Deflection: and Meandering Velocity\n');
for alpha = 1:9
   for gamma = 1:9
      x = {'',['WkVelXT' num2str(alpha,'%1d') 'D' num2str(gamma,'%1d')], ' ', ['X-component of the Advection: deflection: and meandering velocity (not including the horizontal wake-deflection correction) of the wake for downstream distance ' num2str(gamma,'%1d') ' of turbine ' num2str(alpha,'%1d') ' in the global coordinate system'], ' ', '(m/s)' };
      csvFun = @(str)sprintf('%s,',str);
      xchar = cellfun(csvFun, x, 'UniformOutput', false);
      xchar = strcat(xchar{:});
      xchar = strcat(xchar(1:end-1),'\n');
      fprintf(fid,xchar);
   end
end  
for alpha = 1:9
   for gamma = 1:9
      x = {'',['WkVelYT' num2str(alpha,'%1d') 'D' num2str(gamma,'%1d')], ' ', ['Y-component of the Advection: deflection: and meandering velocity (not including the horizontal wake-deflection correction) of the wake for downstream distance ' num2str(gamma,'%1d') ' of turbine ' num2str(alpha,'%1d') ' in the global coordinate system'], ' ', '(m/s)' };
      csvFun = @(str)sprintf('%s,',str);
      xchar = cellfun(csvFun, x, 'UniformOutput', false);
      xchar = strcat(xchar{:});
      xchar = strcat(xchar(1:end-1),'\n');
      fprintf(fid,xchar);
   end
 end  
for alpha = 1:9
   for gamma = 1:9
      x = {'',['WkVelZT' num2str(alpha,'%1d') 'D' num2str(gamma,'%1d')], ' ', ['Z-component of the Advection: deflection: and meandering velocity (not including the horizontal wake-deflection correction) of the wake for downstream distance ' num2str(gamma,'%1d') ' of turbine ' num2str(alpha,'%1d') ' in the global coordinate system'], ' ', '(m/s)' };
      csvFun = @(str)sprintf('%s,',str);
      xchar = cellfun(csvFun, x, 'UniformOutput', false);
      xchar = strcat(xchar{:});
      xchar = strcat(xchar(1:end-1),'\n');
      fprintf(fid,xchar);
   end
end    

   % Wake diameter for downstream distance ? of turbine ?
fprintf(fid,'Wake Diameter\n');
for alpha = 1:9
   for gamma = 1:9
      x = {'',['WkDiamT' num2str(alpha,'%1d') 'D' num2str(gamma,'%1d')], ' ', ['Wake diameter for downstream distance ' num2str(gamma,'%1d') ' of turbine ' num2str(alpha,'%1d')], ' ', '(m)' };
      csvFun = @(str)sprintf('%s,',str);
      xchar = cellfun(csvFun, x, 'UniformOutput', false);
      xchar = strcat(xchar{:});
      xchar = strcat(xchar(1:end-1),'\n');
      fprintf(fid,xchar);
   end
end    

   % Axial and radial wake velocity deficits for radial node ? and downstream distance ? of turbine ?
fprintf(fid,'Axial and Radial Wake Velocity Deficits\n');
for alpha = 1:9
   for beta = 1:20
      for gamma = 1:9
         x = {'',['WkDfVxT' num2str(alpha,'%1d') 'N' num2str(beta,'%02d') 'D' num2str(gamma,'%1d')], ' ', ['Axial wake velocity deficits for radial node ' num2str(beta,'%02d') ' and downstream distance ' num2str(gamma,'%1d') ' of turbine ' num2str(alpha,'%1d')], ' ', '(m/s)' };
         csvFun = @(str)sprintf('%s,',str);
         xchar = cellfun(csvFun, x, 'UniformOutput', false);
         xchar = strcat(xchar{:});
         xchar = strcat(xchar(1:end-1),'\n');
         fprintf(fid,xchar);
      end
   end
end  
for alpha = 1:9
   for beta = 1:20
      for gamma = 1:9
         x = {'',['WkDfVrT' num2str(alpha,'%1d') 'N' num2str(beta,'%02d') 'D' num2str(gamma,'%1d')], ' ', ['Radial wake velocity deficits for radial node ' num2str(beta,'%02d') ' and downstream distance ' num2str(gamma,'%1d') ' of turbine ' num2str(alpha,'%1d')], ' ', '(m/s)' };
         csvFun = @(str)sprintf('%s,',str);
         xchar = cellfun(csvFun, x, 'UniformOutput', false);
         xchar = strcat(xchar{:});
         xchar = strcat(xchar(1:end-1),'\n');
         fprintf(fid,xchar);
      end
   end
end  

   % Total eddy viscosity, and individual contributions to the eddy viscosity from ambient turbulence and the shear layer, for radial node ? and downstream distance ? of turbine ?
fprintf(fid,'Total Eddy Viscosity and Individual Contributions\n');
for alpha = 1:9
   for beta = 1:20
      for gamma = 1:9
         x = {'',['EddVisT' num2str(alpha,'%1d') 'N' num2str(beta,'%02d') 'D' num2str(gamma,'%1d')], ' ', ['Total eddy viscosity for radial node ' num2str(beta,'%02d') ' and downstream distance ' num2str(gamma,'%1d') ' of turbine ' num2str(alpha,'%1d')], ' ', '(m^2/s)' };
         csvFun = @(str)sprintf('%s,',str);
         xchar = cellfun(csvFun, x, 'UniformOutput', false);
         xchar = strcat(xchar{:});
         xchar = strcat(xchar(1:end-1),'\n');
         fprintf(fid,xchar);
      end
   end
end  
for alpha = 1:9
   for beta = 1:20
      for gamma = 1:9
         x = {'',['EddAmbT' num2str(alpha,'%1d') 'N' num2str(beta,'%02d') 'D' num2str(gamma,'%1d')], ' ', ['Contribution to the eddy viscosity from ambient turbulence for radial node ' num2str(beta,'%02d') ' and downstream distance ' num2str(gamma,'%1d') ' of turbine ' num2str(alpha,'%1d')], ' ', '(m^2/s)' };
         csvFun = @(str)sprintf('%s,',str);
         xchar = cellfun(csvFun, x, 'UniformOutput', false);
         xchar = strcat(xchar{:});
         xchar = strcat(xchar(1:end-1),'\n');
         fprintf(fid,xchar);
      end
   end
end  
for alpha = 1:9
   for beta = 1:20
      for gamma = 1:9
         x = {'',['EddShrT' num2str(alpha,'%1d') 'N' num2str(beta,'%02d') 'D' num2str(gamma,'%1d')], ' ', ['Contribution to the eddy viscosity from the shear layer for radial node ' num2str(beta,'%02d') ' and downstream distance ' num2str(gamma,'%1d') ' of turbine ' num2str(alpha,'%1d')], ' ', '(m^2/s)' };
         csvFun = @(str)sprintf('%s,',str);
         xchar = cellfun(csvFun, x, 'UniformOutput', false);
         xchar = strcat(xchar{:});
         xchar = strcat(xchar(1:end-1),'\n');
         fprintf(fid,xchar);
      end
   end
end

% Ambient Wind and Array effects

   % Ambient wind velocity (not including wakes) for point ? in global coordinates (from the low-resolution domain)
fprintf(fid,'Ambient Wind Velocity from Low-resolution Domain\n');
for beta = 1:9
   x = {'',['W' num2str(beta,'%02d') 'VAmbX'], ' ', ['X-component of the ambient wind velocity (not including wakes) for point ' num2str(beta,'%02d') ' in global coordinates (from the low-resolution domain)'], ' ', '(m/s)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);   
end

for beta = 1:9
   x = {'',['W' num2str(beta,'%02d') 'VAmbY'], ' ', ['Y-component of the ambient wind velocity (not including wakes) for point ' num2str(beta,'%02d') ' in global coordinates (from the low-resolution domain)'], ' ', '(m/s)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);   
end   

for beta = 1:9
   x = {'',['W' num2str(beta,'%02d') 'VAmbZ'], ' ', ['Z-component of the ambient wind velocity (not including wakes) for point ' num2str(beta,'%02d') ' in global coordinates (from the low-resolution domain)'], ' ', '(m/s)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);   
end   
   
   
   % Disturbed wind velocity (including wakes) for point ? in the global coordinate system (from the low-resolution domain)
fprintf(fid,'Disturbed Wind Velocity from Low-resolution Domain\n');
for beta = 1:9
   x = {'',['W' num2str(beta,'%02d') 'VDisX'], ' ', ['X-component of the disturbed wind velocity (including wakes) for point ' num2str(beta,'%02d') ' in global coordinates (from the low-resolution domain)'], ' ', '(m/s)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);   
end

for beta = 1:9
   x = {'',['W' num2str(beta,'%02d') 'VDisY'], ' ', ['Y-component of the disturbed wind velocity (including wakes) for point ' num2str(beta,'%02d') ' in global coordinates (from the low-resolution domain)'], ' ', '(m/s)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);   
end   

for beta = 1:9
   x = {'',['W' num2str(beta,'%02d') 'VDisZ'], ' ', ['Z-component of the disturbed wind velocity (including wakes) for point ' num2str(beta,'%02d') ' in global coordinates (from the low-resolution domain)'], ' ', '(m/s)' };
   csvFun = @(str)sprintf('%s,',str);
   xchar = cellfun(csvFun, x, 'UniformOutput', false);
   xchar = strcat(xchar{:});
   xchar = strcat(xchar(1:end-1),'\n');
   fprintf(fid,xchar);   
end     
   
   
   
fclose(fid);
end
