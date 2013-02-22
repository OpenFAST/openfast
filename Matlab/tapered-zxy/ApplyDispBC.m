function [ K, M, F] = ApplyDispBC( K,M, BCNodes, F )
%APPLYDISPBC Summary of this function goes here
%   Detailed explanation goes here

[m, n] = size(BCNodes);
n = max(m, n);

BCs = ones(n, 7);
BCs(:, 1) = BCNodes;

for i = 1:n
    nodenumber = BCs(i, 1);
    for j = 1:6
        tii = (nodenumber-1)*6 + j; % target row/column
        
        K(:, tii) = 0;
        K(tii, :) = 0;
        K(tii, tii) = 1;

        M(:, tii) = 0;
        M(tii, :) = 0;
        M(tii, tii) = 0;

        F(tii) = 0;
    end


end







end % function

