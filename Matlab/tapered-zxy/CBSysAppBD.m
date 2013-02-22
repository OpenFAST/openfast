function [MMCB, KKCB, CBrDOF, IdRCB, FCB] = CBSysAppBD(IdR, CstrNodes, MMCB, KKCB, nCstrNodes, FCB)
%CBSYSAPPBD Summary of this function goes here
%   Detailed explanation goes here


nCdof = nCstrNodes*6;
bdID = zeros(1,nCdof);
nIdR = max(size(IdR));
CBrDOF = max(size(KKCB)) - nCdof;
IdRCB = IdR;

trgtdof = zeros(1, nCdof);
for i = 1:nCstrNodes
    NodeNum = CstrNodes(i);
    trgtdof(i*6-5 : i*6) = NodeNum*6-5 : NodeNum*6;
end

for i = 1:nCdof
    tdof = trgtdof(i);

    for j = 1:nIdR
        if( IdR(j) == tdof )
            bdID(i) = j;
        end
    end
end

IdRCB(bdID) = [];
MMCB(bdID, :) = [];
MMCB(:, bdID) = [];

KKCB(bdID, :) = [];
KKCB(:, bdID) = [];

FCB(bdID) = [];


end

