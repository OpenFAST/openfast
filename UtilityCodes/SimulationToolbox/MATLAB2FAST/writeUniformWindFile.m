function writeUniformWindFile(wind,filename)
Windheader={'!Wind file for sheared 18 m/s wind with 30 degree direction.'
            '!       Wind    Wind    Vert.   Horiz.  Vert.   LinV    Gust'
            '!Time   Speed   Dir     Speed   Shear   Shear   Shear   Speed'};
%convert cell array to char array
Windheader=strvcat(Windheader);

%set up numeric format string
numformat=repmat('%5.2f\t',1,8);
numformat=[numformat,'\n'];

%add any missing columns from the wind input:
[n,m]=size(wind);
if m==1 %assume wind is sampled at 1 Hz
    wind=[(0:(length(wind)-1))',wind,zeros(length(wind(:,1)),6)];
elseif m==2 %assume time and uniform wind
    wind=[wind,zeros(length(wind(:,1)),6)];    
elseif m==4 %assume time, wind speed, h-shear, and l-v-shear
    wind=[wind(:,1:2),zeros(n,2),wind(:,3),zeros(n,1),wind(:,4),zeros(n,1)];
% % elseif m>4 %assume all components are specified

end;

windf=fopen(filename,'w');
if windf>=0
    %print wind header lines 1-by-1
    for index=1:length(Windheader(:,1))
        fprintf(windf,'%s',Windheader(index,:));
        fprintf(windf,'\n');
    end;
    
    %use format string to print wind data in one shot
    fprintf(windf,numformat,wind');

    fclose(windf);
else
    error('Could not open file for writing');
end
