function [ IdRjBar, PhiRjBar, IdBC ] = GetPhiRbar( IdR, PhiR, CstrNodes )
%GETPHIRBAR Summary of this function goes here
%   Detailed explanation goes here

n = max(size(IdR));

nc = max(size(CstrNodes));

innc = zeros(1, nc*6);

for i = 1: nc
    for j = 1:6
        nj = (CstrNodes(i)-1)*6 + j;
        
        for k = 1:n
            if(IdR(k) == nj)
                innc( (i-1)*6 + j ) = k;
            end
        end
        
    end
end

IdBC = IdR(innc);

IdR(innc) = [];
PhiR(:, innc) = [];

IdRjBar = IdR;
PhiRjBar = PhiR;

end

