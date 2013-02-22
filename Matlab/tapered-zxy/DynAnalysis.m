% dynamic analysis
% calculate eigen values and eigen vectors

clc;
clear all;

% filename = 'TaperedBeam2.txt';
% filename = 'JacketData.txt';
% filename = 'test1.txt';
 filename = 'jacketonly.txt';
% filename = 'TripodOnly.txt';
% filename = 'OC4JacketTowerRNA.txt';

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
Div = 2;
NumRdModes = 1008;

 CstrNodes = [61, 62, 63, 64];
% CstrNodes = [81, 82];
% CstrNodes = [65,66,67];
 IntFcNodes = [53, 54, 55, 56];
 nCstrNodes = 4;
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

% t = cputime;
% % from compressed row to full matrix
% [KTCR, MTCR] = RcvKMfCR(iK, jK, KCR, MCR, TDOF);
% t = cputime - t;
% display(['Recover K and M ',num2str(t)]);
% 
% errK = 0;
% errM = 0;
% for i = 1:TDOF
%     for j = 1:TDOF
%         errK = errK + ( KT(i, j)-KTCR(i,j) )^2;
%         errM = errM + ( MT(i, j)-MTCR(i,j) )^2;
%     end
% end
% errK = sqrt(errK);
% errM = sqrt(errM);
% 
%  nCstrNodes = max(size(CstrNodes));
%  rTDOF = TDOF - 6*nCstrNodes;
% 
% %    KT = KTCR;
% %    MT = MTCR;
% 
% % [Kr, Mr] = ReducedKM(KT, MT, CstrNodes);
% 
% applied force 
Fv = ForceVector(NFc, Fcs, NNode);
rhs = Fv + Gv;
% 
% 
% [Kr, Mr, rhsr] = ApplyDispBC(KT, MT, CstrNodes, rhs);
% rTDOF = TDOF;
% 
% 
% % rhs = Fv;
% %-----------------------------------------------------------------------
% %   static solve for full fem 
% %----------------------------------------------------------------------- 
% t = cputime;
% disp('Full fem static solve: ')
% 
% disp_u = Kr\rhsr;
% 
% dt = cputime -t;
% disp(dt);
% 
% 
%  
% %-----------------------------------------------------------------------
% %   eigen solve for full fem 
% %----------------------------------------------------------------------- 
% % 
% t = cputime;
% disp('Full fem eigen solve: ')
% [V, D]= eig(Kr, Mr);
% dt = cputime -t;
% disp(dt);
% 
% 
% for i = 1: rTDOF
%     omega(i) = sqrt(D(i, i))/2/pi;
% end
% 
% 
% OO=sort(omega)';
% 
% 
% 
% % sort frequencies and mode shapes
% a = omega;
% b = V;
% c = [a;b];
% d = sortrows(c', 1);
% d = d';
% 
% omega = d(1,:);
% 
% V1 = d((2:rTDOF+1),:);
% 
% 
% % add boundary nodal displacements to mode shapes
% bdnodes = sort(CstrNodes);
% Vx = zeros(TDOF, rTDOF);
% Vx(1:rTDOF, 1:rTDOF) = V1;
% bdblock = zeros(6, rTDOF);
% 
% for i = 1:nCstrNodes
%     bdn = bdnodes(i);
%     start_insert = bdn*6-5;
%     end_insert = bdn*6;
%     if(start_insert < rTDOF)
%         block1 = Vx(1:start_insert-1, :);
%         block2 = Vx(start_insert:TDOF, :);
%         Vnew = [block1; bdblock; block2];
%         Vx = Vnew(1:TDOF, 1:rTDOF);
%     end
% 
% end
% 
% % plot the mode shapes
%  maxn = min(20,rTDOF);
%  omega1 = omega(1:maxn)';
%  maxnn = 5;
% %V2 = Vx(:, 1:maxnn);
% V2 = V1(:, 1:maxnn);
% 
% 
% fac = 1;
% %PlotModalResultsJackets(maxnn, NElem, Nodes, Elems, V2, omega1, fac, ElemType);
% 
% 



% 
% 
% filename = 'C:\Documents and Settings\hsong\My Documents\Visual Studio 2008\Projects\BeamFEM\BeamFEM\testframe_eigen_results.txt';
% fpr = fopen(filename, 'r');
% 
% %temp = fscanf(fpr, '%g',[1,1]  );
% Ka = fscanf(fpr, '%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g',[24, 24]);
% %temp = fscanf(fpr, '%s' );
% Ma = fscanf(fpr, '%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g',[24, 24]);
% 
% 
% [Va, Da]= eig(Ka, Ma);
% 
% for i = 1:24
%     wa(i) = sqrt(Da(i, i))/2/pi;
% end
% 
% wa = sort(wa)';
% 
% fclose(fpr);
% 
% 
% errKa = Ka - Kr;
% errMa = Ma - Mr;
% 
% 
% 










% TowerTopNode = 132;
% TowerBtmNode = 41;
% 
% nn1 = TowerTopNode*6-5;
% nn2 = TowerBtmNode*6-5;
% 
% u_tower_top = sqrt( disp_u(nn1)^2+disp_u(nn1+1)^2); 
% u_tower_bottom = sqrt( disp_u(nn2)^2+disp_u(nn2+1)^2); 
% 


%--------------------------------------------------------------------------
% C-B reduction

%STOP here;

t = cputime;
disp('Begin C-B: break sys matrices');
disp(t);

RNodes = [CstrNodes, IntFcNodes];

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
TPnode = [0, 0, 20.15];

[MMCBp, KKCBp, CBpDOF, FCBp] ...
        = ReduIntFaceToPoint(MMCBr, KKCBr, TPnode, IdRCB, CBrDOF, Nodes, FCBr);

 
    
KKCB = [];
MMCB = [];
MMCBr=[];
KKCBr=[];
KLL=[];
MLL=[];
KT=[];
MT=[];
PhiL = [];    
    
t = cputime;
disp('C-B eigen solve: ')
[VCB, DCB]= eig(KKCBp, MMCBp);

dt = cputime -t;
disp(dt);


for i = 1:CBpDOF
    OmegaCB(i) = sqrt(DCB(i, i))/2/pi;
end

OmegaCB = sort(OmegaCB);
omega2 = OmegaCB(1:CBpDOF)';

t = cputime;
disp('Program ends at: ')
disp(t);



%--------------------------------------------------------------------------














