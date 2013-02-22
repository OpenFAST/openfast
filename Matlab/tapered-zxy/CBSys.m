 function [MM, KK] = CBSys(MBB, MmB, KBB, OmegaL, IdR, NumRdModes)

%SYSCB Summary of this function goes here
%   Detailed explanation goes here

nIdR = max(size(IdR));
CBDOF = nIdR + NumRdModes;

MM = zeros(CBDOF, CBDOF);
KK = zeros(CBDOF, CBDOF);

nBB = 1:nIdR;
nmm = (nIdR + 1) : CBDOF;

MM(nBB, nBB) = MBB;
MM(nBB, nmm) = MmB';
MM(nmm, nBB) = MmB;
MM(nmm, nmm) = eye(NumRdModes);

KK(nBB, nBB) = KBB;

for i = 1:NumRdModes
    KK(nIdR+i, nIdR+i) = OmegaL(i)^2;
end

end % function

