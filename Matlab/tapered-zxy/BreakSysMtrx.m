function [KRR, KLL, KRL, MRR, MLL, MRL, IdR, IdL, FR, FL] = BreakSysMtrx(KT, MT, rhs, RNodes, TDOF)

%CB Summary of this function goes here
%   Detailed explanation goes here

nr = max(size(RNodes));
IdTR = zeros(1, TDOF);

for i = 1:nr
    inr = RNodes(i);
    IdTR(inr*6-5 : inr*6) = inr*6-5 : inr*6;
end

IdT = 1:TDOF;

IdTT = [IdTR; IdT];
IdTT = sortrows(IdTT', -1);
IdTT = IdTT';
IdT = IdTT(2,:);

IdR = IdT(1:nr*6);
IdR = sort(IdR);
IdL = IdT(nr*6+1 : TDOF);
IdL = sort(IdL);

KRR = KT(IdR, IdR);
KRL = KT(IdR, IdL);
KLL = KT(IdL, IdL);

MRR = MT(IdR, IdR);
MRL = MT(IdR, IdL);
MLL = MT(IdL, IdL);

FR = rhs(IdR);
FL = rhs(IdL);


end

