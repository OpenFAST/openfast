function [ F ] = ForceVector( NFc, Fcs, NNode )
%FORCEVECTOR Summary of this function goes here
%   Detailed explanation goes here

F = zeros(NNode*6, 1);

for i = 1: NFc
    n = Fcs(i, 1);
    F(n*6-5) = Fcs(i, 2);
    F(n*6-4) = Fcs(i, 3);
    F(n*6-3) = Fcs(i, 4);
    F(n*6-2) = Fcs(i, 5);
    F(n*6-1) = Fcs(i, 6);
    F(n*6-0) = Fcs(i, 7);
end



end

