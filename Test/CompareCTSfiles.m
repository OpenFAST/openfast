function CompareCTSfiles(FileName1,FileName2)


fntsz = 8;
    
    
    fprintf( '\n\n\n%45s\n\n', 'Comparing TurbSim Coherent Turbulence *.cts files...');                   
    
    fprintf( '%2s%10s = %s\n', '','File 1 ', FileName1);
    fprintf( '%2s%10s = %s\n', '','File 2 ', FileName2);
    fprintf( '\n' );

    
    [time_1,cts_1,data_1,lab_1] = readCTSfile(FileName1);
    [time_2,cts_2,data_2      ] = readCTSfile(FileName2);

    diff2 = data_1 - data_2;  %difference in header scaling data
    
    nl = length(lab_1);
    fprintf( '%2s%s\n', '','Differences between files: ');
    fprintf( '%4s ', ' ');
    for i = 1:nl
        fprintf( '%12s ', lab_1{i});
    end
    fprintf( '%12s ', 'Time', 'CTS File');
    fprintf( '\n%4s', ' ');
    for i = 1:nl
        fprintf( '%12s ', '(?)');
    end
    fprintf( '%12s ', '(s)', '(#)');
    
    fprintf( '\n%4s ', 'Max ' );
    for i=1:nl
        fprintf( '%12.4f ', diff2(i) );           
    end
    

    if any(size(cts_1) ~= size(cts_2))
        fprintf('\n');
        fprintf( '%2s%s\n', '','Dimensions of data not the same: ');
        fprintf( '%5s%10s = (%8.0f,%3.0f)', '','File 1', size(cts_1) );
        fprintf(   ' %10s = (%8.0f,%3.0f)\n',  'File 2', size(cts_2) );                  
    else
        nr = length(time_1);        
        diff = [time_1 cts_1] - [time_2 cts_2];                

        for i=1:2
            fprintf( '%12.4f ', norm(diff(:,i),inf) );           
        end
        fprintf( '\n%4s ', 'RMS ' );
        for i=1:nl
            fprintf( '%12.4f ', diff2(i) );           
        end
        for i=1:2
            fprintf( '%12.4f ', norm(diff(:,i),2)/sqrt(nr) );           
        end
        fprintf( '\n%4s ', 'Min ');
        for i=1:nl
            fprintf( '%12.4f ', diff2(i) );           
        end
        for i=1:2
            fprintf( '%12.4f ', norm(diff(:,i),-inf) );           
        end
    end
    fprintf( '\n' );        
       
    figure;
    
    % plot the file # vs. time 
    subplot(5,1,1:2);
    
    plot(time_1,cts_1,'bo', time_2,cts_2,'r+');
    set(gca,'FontSize',fntsz);
    xlabel('Time (s)');
    ylabel('CTS file (#)');
    h=legend(FileName1,FileName2,'location','northwest');
    set(h,'Interpreter','none','FontSize',fntsz);
    
    if length(time_1) == length(time_2)
        subplot(5,1,3);
        plot(time_1, diff);        
        legend({'Time (s)','File (#)'},'fontsize',fntsz);
        xlabel('Time (s) [from File 1]');
        ylabel('Differences');
    end

    subplot(5,1,4:5)
    bar([data_1 data_2 diff2]);
    set(gca,'xticklabel',lab_1);
    h=legend(FileName1,FileName2,'Difference','location','northwest');
    set(h,'Interpreter','none','FontSize',fntsz);
    
    
return;