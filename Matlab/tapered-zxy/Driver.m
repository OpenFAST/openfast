
clc;
clear all;



% % tower and RNA
% filename = 'OC4JacketTowerRNA.txt';
filename = 'TestAssembleTop.txt';
Div = 1;
NumRdModes = 1;
nCstrNodes = 0;
CstrNodes = [];
IntFcNodes = [];
TypeFlag = 0; % Full FEM
TPnode = [0,0,18.15];

[ Krna, Mrna, Frna, rTDOFrna, TDOFrna, IdRrna, IdLrna, PhiRrna, PhiLrna, Tirna ]...
    = BFEM( filename, Div, NumRdModes, CstrNodes,  IntFcNodes, TypeFlag , TPnode);
% 

% % jacket only Full FEM
% filename = 'jacketonly.txt';
% Div = 2;
% NumRdModes = 1;
% nCstrNodes = 4;
% CstrNodes = [61, 62, 63, 64];
% IntFcNodes = [1];
% TypeFlag = 0; % Full FEM
% TPnode = [0,0,18.15];
% 
% [ K, M, F, rTDOF, TDOF ] = BFEM( filename, Div, NumRdModes, CstrNodes,  IntFcNodes, TypeFlag , TPnode);


filename = 'TestAssembleBase.txt';
Div = 2;
NumRdModes = 1;
nCstrNodes = 1;
CstrNodesj = [1];
IntFcNodesj = [2];
TypeFlag = 1; % C-B reduction
TPnode = [0,0,45];

% % jacket only C-B
% filename = 'OC4JacketSub.txt';
% Div = 1;
% NumRdModes = 10;
% nCstrNodes = 4;
% CstrNodesj = [61, 62, 63, 64];
% IntFcNodesj = [53, 54, 55, 56];
% TypeFlag = 1; % C-B reduction
% TPnode = [0,0,18.15];
% 
[ Kj, Mj, Fj, rTDOFj, TDOFj, IdRj, IdLj, PhiRj, PhiLj, Tij ]...
    = BFEM( filename, Div, NumRdModes, CstrNodesj,  IntFcNodesj, TypeFlag , TPnode);



% assemble rna and jacket
TDOF = rTDOFrna + rTDOFj - 6;

KT1 = zeros(TDOF, TDOF);
MT1 = zeros(TDOF, TDOF);

KT2 = zeros(TDOF, TDOF);
MT2 = zeros(TDOF, TDOF);

FT1 = zeros(TDOF, 1);
FT2 = zeros(TDOF, 1);

KT1(1:rTDOFrna, 1:rTDOFrna) = Krna;
KT2(rTDOFrna-5 : TDOF,  rTDOFrna-5 : TDOF) = Kj;
KT = KT1+KT2;

MT1(1:rTDOFrna, 1:rTDOFrna) = Mrna;
MT2(rTDOFrna-5 : TDOF,  rTDOFrna-5 : TDOF) = Mj;
MT = MT1+MT2;

FT1(1:rTDOFrna) = Frna;
FT2(rTDOFrna-5 : TDOF) = Fj;
FT = FT1 + FT2; 

rTDOF = TDOF;


%-----------------------------------------------------------------------
%   static solve for full fem 
%----------------------------------------------------------------------- 
disp('Full fem static solve: ')

Frhs = zeros(TDOF, 1);
Frhs(1) = 4000e3;


disp_u = KT\Frhs;

disp(disp_u(1));

disp(disp_u(TDOF-11-NumRdModes));


%-----------------------------------------------------------------------
%   Recover displacement of the internal DOFs of Jacket
%----------------------------------------------------------------------- 


disp_j = disp_u(TDOF-5-NumRdModes : TDOF);
disp_TP = disp_j(1:6);

[IdRjBar, PhiRjBar, IdBC] = GetPhiRbar(IdRj, PhiRj, CstrNodesj );

UL = [PhiRjBar*Tij, PhiLj]*disp_j;
URbar = Tij*disp_TP;

UBc = zeros(size(IdBC));

IdJ = [IdBC, IdRjBar, IdLj]';
UJ = [UBc'; URbar; UL];

U_J = [IdJ, UJ];

U_J = sort(U_J, 1);

% idn = 1:64;
% U_X = U_J(idn*6-5, :);

%-----------------------------------------------------------------------
%   Caculate member forces
%----------------------------------------------------------------------- 



%%%%%%%%%%%
   Stop
%%%%%%%%%%%

%-----------------------------------------------------------------------
%   eigen solve for full fem 
%----------------------------------------------------------------------- 
% 
% rTDOF = rTDOFj;
% TDOF = TDOFj;
% KT = Kj;
% MT = Mj;

t = cputime;
disp('Full fem eigen solve: ')
[V, D]= eig(KT, MT);
dt = cputime -t;
disp(dt);


for i = 1: rTDOF
    omega(i) = sqrt(D(i, i))/2/pi;
end


OO=sort(omega)';



% sort frequencies and mode shapes
a = omega;
b = V;
c = [a;b];
d = sortrows(c', 1);
d = d';

omega = d(1,:);

V1 = d((2:rTDOF+1),:);


% add boundary nodal displacements to mode shapes
bdnodes = sort(CstrNodes);
Vx = zeros(TDOF, rTDOF);
Vx(1:rTDOF, 1:rTDOF) = V1;
bdblock = zeros(6, rTDOF);

for i = 1:nCstrNodes
    bdn = bdnodes(i);
    start_insert = bdn*6-5;
    end_insert = bdn*6;
    if(start_insert < rTDOF)
        block1 = Vx(1:start_insert-1, :);
        block2 = Vx(start_insert:TDOF, :);
        Vnew = [block1; bdblock; block2];
        Vx = Vnew(1:TDOF, 1:rTDOF);
    end

end

% plot the mode shapes
 maxn = min(20,rTDOF);
 omega1 = omega(1:maxn)';
 maxnn = 5;
%V2 = Vx(:, 1:maxnn);
V2 = V1(:, 1:maxnn);


fac = 1;
%PlotModalResultsJackets(maxnn, NElem, Nodes, Elems, V2, omega1, fac, ElemType);












%-----------------------------------------------------
%  check the frequency of a single beam structure
%-----------------------------------------------------

I = 0.25*pi*0.5^4;
A = pi*0.5^2;
m = 7850*A;
L = 90;

o = 22.0345*sqrt(2.1e11*I/m/L^4);
o/pi*0.5

