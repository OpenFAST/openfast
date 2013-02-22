function [KT, MT] = ReducedKM(KT, MT, CstrNodes)
%REDUCEDKM Summary of this function goes here
%   Detailed explanation goes here

[m, n]=size(CstrNodes);
n = max(m,n);
indx = zeros(1, 6*n);

for i = 1:n
    bn = CstrNodes(i);
    a = (bn*6-5):(bn*6);
    
    indx(i*6-5 : i*6) = a;
end

KT(indx, :) = [];
KT(:, indx) = [];

MT(indx, :) = [];
MT(:, indx) = [];


end

