function CompareHHfiles(FileName1,FileName2)

    hhFile  = 1;
    datFile = 2;
    binFile = 3;

    len1 = length(deblank(FileName1));
    len2 = length(deblank(FileName2));

    if strcmpi(FileName1(len1-3:len1),'.dat') && strcmpi(FileName2(len2-3:len2),'.dat') 
        fileType = datFile;
    elseif strcmpi(FileName1(len1-3:len1),'.bin') && strcmpi(FileName2(len2-3:len2),'.bin') 
        fileType = binFile;
    else
        fileType = hhFile;            
    end

    npr = 5;
    npc = 1;
    np  = npr*npc;
    
    
    switch fileType
        case hhFile
            nHeader = 8;
            nUnits  = 6;
            
            fprintf( '\n\n\n%45s\n\n', 'Comparing TurbSim hub-height wind speed *.hh files...');
        case datFile
            nHeader = 4;
            nUnits  = 3;  
            
            fprintf( '\n\n\n%45s\n\n', 'Comparing TurbSim hub-height formatted *.dat files...');
            
        case binFile            
            fprintf( '\n\n\n%45s\n\n', 'Comparing TurbSim hub-height binary *.bin files...');
            
    end
        
    
    fprintf( '%2s%10s = %s\n', '','File 1 ', FileName1);
    fprintf( '%2s%10s = %s\n', '','File 2 ', FileName2);
    fprintf( '\n' );

    if fileType == binFile
        [data_1, time_1, txt_1, dat_units] = readHHbin(FileName1);
        [data_2, time_2, txt_2, dat_units] = readHHbin(FileName2);
    else
        [data_1, time_1, txt_1, dat_units] = loadColumnData(FileName1, 'all', '',[nHeader,0],nUnits);
        [data_2, time_2, txt_2, dat_units] = loadColumnData(FileName2, 'all', '',[nHeader,0],nUnits);
   
        if fileType == datFile
            dat_units{14} = '(s)'; %the first column is time & stored at the end
            for i=2:9
                dat_units{i-1} = '(m/s)';
            end
            for i=10:14
                dat_units{i-1} = '(m/s)^2';
            end
        end
    end
    
    
    if any(size(data_1) ~= size(data_2))
        fprintf( '%2s%43s\n', '','Dimensions of data not the same: ');
        fprintf( '%5s%10s = (%8.0f,%3.0f)', '','File 1', size(data_1) );
        fprintf(   ' %10s = (%8.0f,%3.0f)\n',  'File 2', size(data_2) );  
        
        nc = min( size(data_1,2), size(data_2,2) );
        
    else
        diff = data_1 - data_2;
        [nr,nc]   = size(diff);
        fprintf( '%2s%s\n', '','Differences between files: ');
        fprintf( '%4s', ' ');
        for i=[nc+1, 1:nc]
            fprintf( '%7s ', txt_1{i});
        end
        fprintf( '\n%4s', ' ');
        for i=[nc+1, 1:nc]
            fprintf( '%7s ', dat_units{i});
        end

        fprintf( '\n%4s%7.4f ', 'Max ', norm(time_1-time_2,inf) );
        for i=1:nc
            fprintf( '%7.4f ', norm(diff(:,i),inf) );           
        end
        fprintf( '\n%4s%7.4f ', 'RMS ', norm(time_1-time_2,2)/sqrt(nr) );
        for i=1:nc
            fprintf( '%7.4f ', norm(diff(:,i),2)/sqrt(nr) );           
        end
        fprintf( '\n%4s%7.4f ', 'Min ', norm(time_1-time_2,-inf) );
        for i=1:nc
            fprintf( '%7.4f ', norm(diff(:,i),-inf) );           
        end
        
        fprintf( '\n' );
    end
       
    i_plt = 1;
    for i=1:nc
        if i_plt == 1
            figure;
            lgnd = true;
        else
            lgnd = false;
        end
        subplot(npr,npc,i_plt)
%bjj         subplot(npr,npc,i_plt,'v6'); % using v6 is necessary to align the two axes in plotyy.  Otherwise they do not overlap in R2006b.
        i_plt = mod(i_plt,np)+1; 
        
        if length(time_1) == length(time_2)
            [ax, h1 ] = plotyy( [time_1 time_2], [data_1(:,i) data_2(:,i)], time_1, diff(:,i) );
            set(ax(2),'FontSize',8);
            set(get(ax(2),'Ylabel'),'String',['Difference',         ' ', dat_units{i} ]);
        else
            h1 = plot( time_1, data_1(:,i), time_2, data_2(:,i) );
            ax = gca;
        end
        
        set(ax(1),'FontSize',8);
        set(get(ax(1),'Ylabel'),'String',[ txt_1{i}, ' ', dat_units{i} ]);
        
%         plot(time_1, diff(:,i));
        xlabel([ txt_1{nc+1} ' ' dat_units{nc+1}], 'FontSize',7);
        if lgnd
            title( ['TurbSim Hub-Height File Comparison: ' FileName1 ' and ' FileName2], 'Interpreter','none' );
            
%             h=legend(FileName1, FileName2);
%             set(h,'Interpreter','none','FontSize',7)
            legend(h1,{FileName1, FileName2}, 'Interpreter','none','FontSize',7); %must specifiy h1 for R2008a
         
% If the title is called after legend(), the "interpreter" property doesn't
% get set correctly.  So, I moved it earlier.
%             title( ['TurbSim Hub-Height File Comparison: ' FileName1 ' and ' FileName2], 'Interpreter','none' );
        end

    end
    
return;