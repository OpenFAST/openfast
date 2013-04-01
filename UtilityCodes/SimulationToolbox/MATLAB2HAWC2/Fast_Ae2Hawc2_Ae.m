% Fucntion for converting Fast AE and to Hawc2 ae
%
% In:   OutPath         -       Path to folder holding hawc2 files
%       OutFileNames    -       Names of output files .ae and .pc
%       AeroDyn_AE      -       Aerodyne data extracted using
%       AeroDyn2Matlab.m
%       PC_Data         -       Profile coefficient extracted using
%       Fast_PC2Matlab.m
%       FastPar         -       Extracted using Fast2Matlab.m
%
% Knud A. Kragh

function Hawc2_Ae=Fast_Ae2Hawc2_Ae(AeroDyn_AE)

Hawc2_Ae=[];
PC_active=unique(AeroDyn_AE.BldNodes(:,end));

Hawc2_Ae(1,1:4)=[0 AeroDyn_AE.BldNodes(1,4) round(AeroDyn_AE.BldNodes(1,5)) 1];

for i=1:size(AeroDyn_AE.BldNodes,1)
    Hawc2_Ae(i+1,1:4)=[Hawc2_Ae(i,1)+AeroDyn_AE.BldNodes(i,3) AeroDyn_AE.BldNodes(i,4) round(AeroDyn_AE.BldNodes(i,5)) 1];
end

i=1;
while i<size(Hawc2_Ae,1)
    if Hawc2_Ae(i,3)~=Hawc2_Ae(i+1,3)
        Hawc2_Ae(end+1,:)=0;
        Hawc2_Ae(i+1:end,:)=Hawc2_Ae(i:end-1,:);
        Hawc2_Ae(i+1,1)=Hawc2_Ae(i,1)+0.01;
        Hawc2_Ae(i+1,3)=Hawc2_Ae(i+2,3);
        i=i+1;
    end
    i=i+1;
end