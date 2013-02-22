function [MMCBp, KKCBp, CBpDOF, FCBp, T] ...
        = ReduIntFaceToPoint(MMCBr, KKCBr, TPnode, IdRCB, CBrDOF, Nodes, FCBr)

%REDUINTFACETOPOINT Summary of this function goes here
%   Detailed explanation goes here

InterfaceDOF = max(size(IdRCB));

T = zeros(InterfaceDOF, 6);

for i = 1: InterfaceDOF
    a = IdRCB(i);
    rem = mod(a, 6);
    n = ceil(a/6);
    
    x = Nodes(n, 2);
    y = Nodes(n, 3);
    z = Nodes(n, 4);
    
    dx = x - TPnode(1);
    dy = y - TPnode(2);
    dz = z - TPnode(3);
    
    switch rem
        case(1)
            T(i, :) = [1 0 0 0 dz -dy];
        case(2)
            T(i, :) = [0 1 0 -dz 0 dx];
        case(3)
            T(i, :) = [0 0 1 dy -dx 0];
        case(4)
            T(i, :) = [0 0 0 1 0 0];
        case(5)
            T(i, :) = [0 0 0 0 1 0];
        case(0)
            T(i, :) = [0 0 0 0 0 1];
        otherwise
            disp('unknown case number');
    end
end

CBpDOF = CBrDOF-InterfaceDOF + 6;
MMCBp = zeros(CBpDOF, CBpDOF);
MMCBp(1:6, 1:6) = T'*MMCBr(1:InterfaceDOF, 1:InterfaceDOF)*T;
MMCBp(1:6, 7:CBpDOF) = T'*MMCBr(1:InterfaceDOF, InterfaceDOF +1: CBrDOF);
MMCBp(7:CBpDOF, 1:6 ) = MMCBr(InterfaceDOF +1: CBrDOF, 1:InterfaceDOF )*T;
MMCBp(7:CBpDOF, 7:CBpDOF ) = MMCBr(InterfaceDOF +1: CBrDOF,InterfaceDOF +1: CBrDOF);

KKCBp = zeros(CBpDOF, CBpDOF);
KKCBp(1:6, 1:6) = T'*KKCBr(1:InterfaceDOF, 1:InterfaceDOF)*T;
KKCBp(1:6, 7:CBpDOF) = T'*KKCBr(1:InterfaceDOF, InterfaceDOF +1: CBrDOF);
KKCBp(7:CBpDOF, 1:6 ) = KKCBr(InterfaceDOF +1: CBrDOF, 1:InterfaceDOF )*T;
KKCBp(7:CBpDOF, 7:CBpDOF ) = KKCBr(InterfaceDOF +1: CBrDOF,InterfaceDOF +1: CBrDOF);

FCBp = zeros(CBpDOF, 1);
FCBp(1:6) = T'*FCBr(1:InterfaceDOF); 
FCBp(7:CBpDOF) = FCBr(InterfaceDOF +1: CBrDOF);


end

