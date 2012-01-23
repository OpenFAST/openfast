function CompareFFtxtFiles(FileName1,FileName2)

%-----------------------------------
% Full-field time series (*.u .v .w)
%-----------------------------------
    comp = FileName1(length(FileName1):length(FileName1));

    fprintf( '\n\n\n%s\n\n', ['Comparing TurbSim formatted full-field *.' comp ' files...']);
    fprintf( '%2s%10s = %s\n', '','File 1 ', FileName1);
    fprintf( '%2s%10s = %s\n', '','File 2 ', FileName2);
    fprintf( '\n' );    
    
    [velocity_1, y_1, z_1, dt_1, ny_1, nz_1] = loadFFtxt(FileName1);
    [velocity_2, y_2, z_2, dt_2, ny_2, nz_2] = loadFFtxt(FileName2);
    
    PlotColors = jet(ny_1*nz_1 + ny_2*nz_2);
    nt = size(velocity_1,1);
    
        % Calculate the differences....
   
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
        fprintf( '%5s%10s = %3.0f, %10s = %3.0f\n', '','File 1', dt_1, 'File 2', dt_2);
    end
      
        % Velocity(timestep,iy,iz)    
    if any( size(velocity_1) ~= size(velocity_2) )
        fprintf( '%2s%s\n', '','Dimensions of velocity matrix not the same: ');
        fprintf( '%5s%10s = (%8.0f,%3.0f,%3.0f)', '','File 1', size(velocity_1) );
        fprintf(   ' %10s = (%8.0f,%3.0f,%3.0f)\n',  'File 2', size(velocity_2) );
        
        if ny_1 == ny_2 && nz_1 == nz_2
            nt = min(nt,size(velocity_2,1));
            dv = velocity_1(1:nt,:,:) - velocity_2(1:nt,:,:);
            myTxt = sprintf('Maximum difference in velocity in first %s time steps of component:', num2str(nt) );
            fprintf( '%5s%50s %10.5g m/s\n', '', myTxt, norm(dv(:),inf) );
        end
    else
        dv = velocity_1-velocity_2;        
        fprintf( '%2s%50s %10.5g m/s\n', '','Maximum difference in velocity: ', norm(dv(:),inf) );
    end
    
%%
% getMeanPSD( {velocity_1,velocity_2}, {dt_1, dt_2} );

%%
    figure;
    x_1 = [0:size(velocity_1,1)-1]*dt_1;
    x_2 = [0:size(velocity_2,1)-1]*dt_2;
    
        % Plot time series
    subplot(2,1,1)
    hold on;
    ylabel([ comp '-component velocities (m/s)'])
    xlabel([ 'Time (s)'])
    title( ['TurbSim Formatted Full-Field File Comparison: ' FileName1 ' and ' FileName2], 'Interpreter','none' );


    i_color = 1;
    for iy = 1:ny_1
        for iz = 1:nz_1
            plot(x_1, squeeze( velocity_1(:,iy,iz) ), 'Color', PlotColors(i_color,:) );
            i_color = i_color + 1;
        end
    end

    for iy = 1:ny_2
        for iz = 1:nz_2
            plot(x_2, squeeze( velocity_2(:,iy,iz) ), 'Color', PlotColors(i_color,:) );
            i_color = i_color + 1;
        end
    end                   

        % Plot differences
    subplot(2,1,2)
    hold on;
    ylabel({'Differences in ';[ comp '-component velocities (m/s)']})
    xlabel([ 'Time (s)'])


    i_color = 1;
    if exist('dv', 'var')
        x_1 = [0:size(dv,1)-1]*dt_1;
        for iy = 1:size(dv,2)
            for iz = 1:size(dv,3)
                plot(x_1, squeeze( dv(:,iy,iz) ), 'Color', PlotColors(i_color,:) );
                i_color = i_color + 2;
            end
        end
    end
          
    
return;    