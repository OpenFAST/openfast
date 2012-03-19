function CompareWNDfiles(FileName1,FileName2)

%-----------------------------------
% TurbSim/BLADED binary form (*.wnd)
%-----------------------------------
    comp = 'UVW';   

    fprintf( '\n\n\n%s\n\n', 'Comparing TurbSim binary *.wnd files...');
    fprintf( '%2s%10s = %s\n', '','File 1 ', FileName1);
    fprintf( '%2s%10s = %s\n', '','File 2 ', FileName2);
    fprintf( '\n' );
    
    [velocity_1, y_1, z_1, nz_1, ny_1, dz_1, dy_1, dt_1] = readBLgrid(FileName1);
    [velocity_2, y_2, z_2, nz_2, ny_2, dz_2, dy_2, dt_2] = readBLgrid(FileName2);
        
    PlotColors = jet(ny_1*nz_1 + ny_2*nz_2);
    
        % Calculate the differences....
   
    fprintf( '\n' );        
    if ny_1 ~= ny_2
        fprintf( '%2s%s\n', '','Horizontal grid dimensions not the same: ');
        fprintf( '%5s%10s = %3.0f, %10s = %3.0f\n', '','File 1', ny_1, 'File 2', ny_2);
    else
        fprintf( '%2s%50s %10.5g m\n', '','Maximum difference in horizontal grid locations: ', norm(y_1-y_2,inf) );
    end
    
    if nz_1 ~= nz_2
        fprintf( '%2s%s\n', '','Vertical grid dimensions not the same: ');
        fprintf( '%5s%10s = %3.0f, %10s = %3.0f\n', '','File 1', nz_1, 'File 2', nz_2);
    else
        fprintf( '%2s%50s %10.5g m\n', '','Maximum difference in vertical grid locations: ', norm(z_1-z_2,inf) );
    end
    
    if dt_1 ~= dt_2
        fprintf( '%2s%s\n', '','Time steps not the same: ');
        fprintf( '%5s%10s = %6.4f, %10s = %6.4f\n', '','File 1', dt_1, 'File 2', dt_2);
    end
      
        % Velocity(timestep,component,iy,iz)    
    if any( size(velocity_1) ~= size(velocity_2) )
        fprintf( '%2s%s\n', '','Dimensions of velocity matrix not the same: ');
        fprintf( '%5s%10s = (%8.0f,%2.0f,%3.0f,%3.0f)', '','File 1', size(velocity_1) );
        fprintf(   ' %10s = (%8.0f,%2.0f,%3.0f,%3.0f)\n',  'File 2', size(velocity_2) );
        
        if ny_1 == ny_2 && nz_1 == nz_2
            t1 = min(size(velocity_1,1),size(velocity_2,1));
            t2 = min(size(velocity_1,2),size(velocity_2,2));
            dv = velocity_1(1:t1,1:t2,:,:) - velocity_2(1:t1,1:t2,:,:);
            myTxt = sprintf('Maximum difference in velocity in first %s time steps of %s components:', num2str(t1), num2str(t2) );
            fprintf( '%5s%50s %10.5g m/s\n', '', myTxt, norm(dv(:),inf) );
        end
    else
        dv = velocity_1-velocity_2;
        fprintf( '%2s%50s %10.5g m/s\n', '','Maximum difference in velocity: ', norm(dv(:),inf) );
    end
    
%%
% getMeanPSD( {velocity_1,velocity_2}, {dt_1, dt_2} );
    
%%    
%--------------------------------------------------        
% Plot time series
%--------------------------------------------------        

    figure;
    x_1 = [0:size(velocity_1,1)-1]*dt_1;
    x_2 = [0:size(velocity_2,1)-1]*dt_2;
    
    for ic=1:3
            % Plot time series
        subplot(3,2,2*ic-1)
        hold on;
        ylabel([ comp(ic) '-component velocities (m/s)'])
        xlabel(  'time (s)' )
        if ic==1
            title( ['TurbSim Binary WND Full-Field File Comparisons: ' FileName1 ' and ' FileName2], 'Interpreter','none' );
        end
        
        i_color = 1;
                                
        for iy = 1:ny_1
            for iz = 1:nz_1
                plot(x_1, squeeze( velocity_1(:,ic,iy,iz) ), 'Color', PlotColors(i_color,:) );
                i_color = i_color + 1;
            end
        end
                
        for iy = 1:ny_2
            for iz = 1:nz_2
                plot(x_2, squeeze( velocity_2(:,ic,iy,iz) ), 'Color', PlotColors(i_color,:) );
                i_color = i_color + 1;
            end
        end                   
        
            % Plot differences
        subplot(3,2,2*ic)
        hold on;
        ylabel([ 'Differences in ' comp(ic) '-component velocities (m/s)'])
        xlabel(  'time (s)' )
        if ic==1
            title( ['TurbSim Binary Full-Field File Comparisons: ' FileName1 ' and ' FileName2], 'Interpreter','none' );
        end

        i_color = 1;
    
        if exist('dv','var')
            for iy = 1:size(dv,3)
                for iz = 1:size(dv,4)
                    plot(x_1, squeeze( dv(:,ic,iy,iz) ), 'Color', PlotColors(i_color,:) );
                    i_color = i_color + 2;
                end
            end
        end
  
%         if exist('dz','var')
%             for iz = 1:size(dz,3)
%                 plot(x_1, squeeze( dz(:,ic,iz) ), 'Color', PlotColors(i_color,:) );
%                 i_color = i_color + 2;
%             end
%         end
        
    end %for ic
    
%--------------------------------------------------    
% plot the mean wind velocity & direction profiles
%--------------------------------------------------    
    figure;
     
    p1 = zeros(nz_1,2);
    z1 = z_1;
    
    p2 = zeros(nz_2,2);
    z2 = z_2;
    
%     hubY_1 = [floor(ny_1/2) ceil(ny_1/2)];
%     hubY_2 = [floor(ny_2/2) ceil(ny_2/2)];
    hubY_1 = 1:ny_1;
    hubY_2 = 1:ny_2;

    for i = 1:nz_1
        tmp     = sqrt(velocity_1(:,1,hubY_1,i).^2 + velocity_1(:,2,hubY_1,i).^2 + velocity_1(:,3,hubY_1,i).^2);
        p1(i,1) = mean( tmp(:) );
        tmp     = atan2(velocity_1(:,2,hubY_1,i),velocity_1(:,1,hubY_1,i));
        p1(i,2) = mean( tmp(:) )* 180./pi;        
    end        
   
    for i = 1:nz_2
        tmp     = sqrt(velocity_2(:,1,hubY_2,i).^2 + velocity_2(:,2,hubY_2,i).^2 + velocity_2(:,3,hubY_2,i).^2);
        p2(i,1) = mean( tmp(:) );
        tmp     = atan2(velocity_2(:,2,hubY_2,i),velocity_2(:,1,hubY_2,i));
        p2(i,2) = mean( tmp(:) )* 180./pi;        
    end

    subplot(1,2,1)
    plot(p1(:,1),z1, '-b+', p2(:,1),z2, ':ro');

    title( ['TurbSim Binary Full-Field File Comparisons: ' FileName1 ' and ' FileName2], 'Interpreter','none' );

    xlabel('Mean Total Wind Speed (m/s)');
    ylabel('Height above ground (m)');
    h=legend(FileName1, FileName2);  
    set(h,'Interpreter','none','FontSize',7);   %note: title() cannot follow this line, or the interpreter has trouble!

   
    subplot(1,2,2)
    plot(p1(:,2),z1,'-b+',p2(:,2),z2,':ro');
    title( ['TurbSim Binary Full-Field File Comparisons: ' FileName1 ' and ' FileName2], 'Interpreter','none' );

    xlabel('Mean Horizontal Wind Angle (deg, ccw from +x)');
    ylabel('Height above ground (m)');
    h=legend(FileName1, FileName2);  
    set(h,'Interpreter','none','FontSize',7);  %note: title() cannot follow this line, or the interpreter has trouble!
    
    
return;    