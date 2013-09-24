clear all;
close all;
clc;

[Channels, ChanName, ChanUnit, FileID] = ReadFASTbinary('Test20.outb');

for i = 1:length(ChanName)
    if strcmp(ChanName(i),'Fair1Ten')
        index_fair1ten = i;
    end
    if strcmp(ChanName(i),'Fair2Ten')
        index_fair2ten = i;
    end
    if strcmp(ChanName(i),'PtfmSurge')
        index_surge = i;
    end
    if strcmp(ChanName(i),'PtfmSway')
        index_sway = i;
    end
    if strcmp(ChanName(i),'PtfmHeave')
        index_heave = i;
    end
    if strcmp(ChanName(i),'PtfmRoll')
        index_roll = i;
    end
    if strcmp(ChanName(i),'PtfmPitch')
        index_pitch = i;
    end
    if strcmp(ChanName(i),'PtfmYaw')
        index_yaw = i;
    end
end

array = [Channels(:,1),Channels(:,index_fair1ten),Channels(:,index_fair2ten),Channels(:,index_surge),Channels(:,index_sway),Channels(:,index_heave),Channels(:,index_roll),Channels(:,index_pitch),Channels(:,index_yaw)];

fid = fopen('test20.out','w');
fprintf(fid,'%s\r\n','');
dlmwrite('test20.out',array,...
            '-append',...  %# Print the matrix
            'delimiter','\t',...
            'newline','pc');
     
fclose(fid);