function [MBB, MmB, KBB, PhiL, OmegaL, PhiR] = CBMatrix(KRR, KRL, KLL, MLL, MRR, MRL, NumRdModes)
%CBMATRIX Summary of this function goes here
%   Detailed explanation goes here

MLR = MRL';
KLR = KRL';

PhiR = -inv(KLL)*KLR;

[Phi, Omega] = eig(KLL, MLL);

LDOF = max(size(KLL));

mu = Phi'*MLL*Phi;
mu2 = eye(LDOF);

for i = 1:LDOF
    mu2(i, i) = mu(i, i)^-0.5;
end

Phi = Phi*mu2;


for i = 1: LDOF
    a(i) = Omega(i,i);
end

b = Phi;
c = [a;b];
d = sortrows(c', 1);
d = d';

OmegaL = sqrt(d(1,1:NumRdModes));
PhiL = d((2:LDOF+1), 1:NumRdModes);


MRLPhiR = MRL*PhiR;
MBB = MRR + MRLPhiR + MRLPhiR' + PhiR'*MLL*PhiR;
MmB = PhiL'*MLR + PhiL'*MLL*PhiR;
KBB = KRR + KRL*PhiR;

end

