function [KTCR, MTCR] = RcvKMfCR(iK, jK, KCR, MCR, TDOF)

%RCVKMFCR Summary of this function goes here
%   Detailed explanation goes here

KTCR = zeros(TDOF, TDOF);
MTCR = zeros(TDOF, TDOF);

for i = 1:TDOF
    j1 = iK(i);
    j2 = iK(i+1) - 1;
    
    for j = j1:j2
        k = jK(j);
        KTCR(i, k) = KCR(j);
        KTCR(k, i) = KTCR(i, k);
        MTCR(i, k) = MCR(j);
        MTCR(k, i) = MTCR(i, k);
    end
    
end


end

