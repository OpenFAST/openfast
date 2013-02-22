function [ K, M, F, rTDOF, TDOF, IdR, IdL, PhiR, PhiL, Ti ] = BFEM( filename, Div, NumRdModes, CstrNodes,  IntFcNodes, TypeFlag , TPnode)
%BEAMFEM Summary of this function goes here
%   Detailed explanation goes here
% dynamic analysis
% calculate eigen values and eigen vectors

% filename = 'TaperedBeam2.txt';
% filename = 'JacketData.txt';
% filename = 'test1.txt';
% filename = 'jacketonly.txt';
% filename = 'TripodOnly.txt';

% filename = 'OC3Monopile.txt';
% filename = 'OC3Tripod0.txt';
% filename = 'ThreeLegFullLength.txt';
% filename = 'MMI7MWd.txt';

global g;
 g = -9.80665;
% g = 0;

global kappa;
kappa = 0.53;

global NFigure;
NFigure = 0;

t = cputime;
disp('Program start at: ');
disp(t);

pi = acos(-1);
ElemType = 1; % 1: 2-node uniform; 2: 2-node taper; 3: 3-node taper
shear = 1;
%Div = 2;
%NumRdModes = 8;

% CstrNodes = [61, 62, 63, 64];
% CstrNodes = [57, 58, 59, 60];
% CstrNodes = [65,66,67];
% IntFcNodes = [53, 54, 55, 56];
% IntFcNodes = [1];

[NNode, NElem, NProp, Nodes, Elems, Props, PMs, NPM, NFc, Fcs] = ReadData(filename, Div, ElemType);


[JointsConnE, JointsConnN] = JointConnect( NNode, NElem, Nodes, Elems );

t = cputime;
[ iK, jK, NNZ ] = InitIkJk( NNode, JointsConnN ); % for two-node element
t = cputime - t;
display(['InitIkJk ',num2str(t)]);


TDOF = NNode*6;

t = cputime;
[KT, MT, KCR, MCR, Gv] = AssembleKM( NElem, Nodes, Elems, Props, TDOF, ...
    pi , shear, ElemType, PMs, NPM, jK, iK, NNZ);
t = cputime - t;
display(['AssembleKM ',num2str(t)]);

t = cputime;
% from compressed row to full matrix
[KTCR, MTCR] = RcvKMfCR(iK, jK, KCR, MCR, TDOF);
t = cputime - t;
display(['Recover K and M ',num2str(t)]);

errK = 0;
errM = 0;
for i = 1:TDOF
    for j = 1:TDOF
        errK = errK + ( KT(i, j)-KTCR(i,j) )^2;
        errM = errM + ( MT(i, j)-MTCR(i,j) )^2;
    end
end
errK = sqrt(errK);
errM = sqrt(errM);

 nCstrNodes = max(size(CstrNodes));
 
  rTDOF = TDOF - 6*nCstrNodes;

%    KT = KTCR;
%    MT = MTCR;

% [Kr, Mr] = ReducedKM(KT, MT, CstrNodes);

% applied force 
Fv = ForceVector(NFc, Fcs, NNode);
rhs = Fv + Gv;


if (nCstrNodes ~= 0)

    [Kr, Mr, rhsr] = ApplyDispBC(KT, MT, CstrNodes, rhs);
    rTDOF = TDOF;

    if (TypeFlag == 0)
        K = Kr;
        M = Mr;
        F = rhsr;
        IdR = 0;
        IdL = 0;
        PhiR = 0;
        PhiL = 0;
        Ti = 0;
        
        return;
    end
else
    if (TypeFlag == 0)
        K = KT;
        M = MT;
        F = rhs;
        IdR = 0;
        IdL = 0;
        PhiR = 0;
        PhiL = 0;
        Ti = 0;
        
        return;
    end
  
    
end
%--------------------------------------------------------------------------
% C-B reduction

%STOP here;

t = cputime;
disp('Begin C-B: break sys matrices');
disp(t);

RNodes = [CstrNodes, IntFcNodes];

%[KRR, KLL, KRL, MRR, MLL, MRL, IdR, IdL] = BreakSysMtrx(KT, MT, RNodes, TDOF);

[KRR, KLL, KRL, MRR, MLL, MRL, IdR, IdL, FR, FL] = BreakSysMtrx(KT, MT, rhs, RNodes, TDOF);


t = cputime;

disp('Begin C-B: CB matrices');

[MBB, MmB, KBB, PhiL, OmegaL, PhiR] = CBMatrix(KRR, KRL, KLL, MLL, MRR, MRL, NumRdModes);

t = cputime - t;
disp(t);

t = cputime;
disp('Begin C-B: CB sys matrix');
disp(t);

[MMCB, KKCB] = CBSys(MBB, MmB, KBB, OmegaL, IdR, NumRdModes);

FCB = [FR + PhiR'*FL; PhiL'*FL]; 

[MMCBr, KKCBr, CBrDOF, IdRCB, FCBr] = CBSysAppBD(IdR, CstrNodes, MMCB, KKCB, nCstrNodes, FCB);

% t = cputime;
% disp('C-B eigen solve: ')
% [VCB, DCB]= eig(KKCBr, MMCBr);
% 
% dt = cputime -t;
% disp(dt);
% 
% 
% for i = 1:CBrDOF
%     OmegaCB(i) = sqrt(DCB(i, i))/2/pi;
% end
% 
% OmegaCB = sort(OmegaCB);
% omega2 = OmegaCB(1:CBrDOF)';
% % omega1 = omega(1:min(rTDOF,100))';
% 
% t = cputime;
% disp('Program ends at: ')
% disp(t);
% 

% reduction of a redundant interface to a point

[MMCBp, KKCBp, CBpDOF, FCBp, Ti] ...
        = ReduIntFaceToPoint(MMCBr, KKCBr, TPnode, IdRCB, CBrDOF, Nodes, FCBr);

 
if (TypeFlag == 1 )
    K = KKCBp;
    M = MMCBp;
    F = FCBp;
    rTDOF = CBpDOF;
end
    
%     
%     
% t = cputime;
% disp('C-B eigen solve: ')
% [VCB, DCB]= eig(KKCBp, MMCBp);
% 
% dt = cputime -t;
% disp(dt);
% 
% 
% for i = 1:CBpDOF
%     OmegaCB(i) = sqrt(DCB(i, i))/2/pi;
% end
% 
% OmegaCB = sort(OmegaCB);
% omega2 = OmegaCB(1:CBpDOF)';
% 
% t = cputime;
% disp('Program ends at: ')
% disp(t);
% 


%--------------------------------------------------------------------------
















end

