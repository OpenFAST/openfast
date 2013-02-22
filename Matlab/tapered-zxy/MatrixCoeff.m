function [AA, BB, CC, DD, EE, FF, GG, HH, PP ] ...
         = MatrixCoeff(a0, a1, a2, b0, b1, b2, b3, b4, E, G, kappa, rho)
%MATRIXCOEFF Summary of this function goes here
%   Detailed explanation goes here
aa = [a0, a1, a2];
bb = [b0, b1, b2, b3, b4];

AA = E*aa;
BB = G*kappa*aa;
CC = G*kappa*aa;
DD = G*2*bb;
EE = E*bb;
FF = E*bb;

GG = rho*aa;
HH = 2*rho*bb;
PP = rho*bb;

end

