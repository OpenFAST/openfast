function [KT, MT, KCR, MCR, Gv] = AssembleKM( NElem, Nodes, Elems, Props, TDOF,...
    pi, shear, ElemType, PMs, NPM, jK, iK, NNZ)
%ASSEMBLEK Summary of this function goes here
%   Detailed explanation goes here

global kappa;

KT = zeros(TDOF, TDOF);
MT = zeros(TDOF, TDOF);
Gv = zeros(TDOF, 1);

KCR = zeros(NNZ, 1);
MCR = zeros(NNZ, 1);

TotalMass = 0;
TowerMass = 0;
JacketMass = 0;
TPMass = 0;
PileMass = 0;

for i = 1:NElem
    if ElemType ==1
        nnie = 2;
    else
        nnie = ElemType ; % total number of nodes in one element
    end
    
    for j = 1:nnie
        nn(j) = Elems(i, j+1);
    end
    
    p1 = Elems(i, nnie+2);
    p2 = Elems(i, nnie+3);
    
    d1 = Props(p1, 5);
    d2 = Props(p2, 5);
    t1 = Props(p1, 6);
    t2 = Props(p2, 6);
    E = Props(p1, 2);
    G = Props(p1, 3);
    rho = Props(p1, 4);
    
    x1 = Nodes(nn(1), 2);
    y1 = Nodes(nn(1), 3);
    z1 = Nodes(nn(1), 4);

    x2 = Nodes(nn(nnie), 2);
    y2 = Nodes(nn(nnie), 3);
    z2 = Nodes(nn(nnie), 4);

    
    DirCos = DirCosMtrx( x1, y1, z1, x2, y2, z2 );
    L = norm([x1-x2, y1-y2, z1-z2]);

    if (ElemType == 1) % two node classical element
        r1 = 0.25*(d1+d2);
        t = 0.5*(t1+t2);

        if( t == 0 )
            r2 = 0;
        else
            r2 = r1 - t;
        end

        A = pi*( r1^2 - r2^2 );
        I = 0.25*pi*( r1^4 - r2^4 );
        ka = kappa;
        J = 0.5*pi*( r1^4 - r2^4 );

        Ki = ElemK(A, L, I, J, ka, E, G, DirCos, shear );    
        Mi = ElemM(A, L, I, J, rho, DirCos);    
        Gi = ElemG(A, L, rho, DirCos);

        
        ElemMass = L*A*rho;
%         ElemMass = pi*L/3*rho*(d2^2+d1^2+d2*d1)/4 - ...
%             pi*L/3*rho*((d2*0.5-t2)^2+(d1*0.5-t1)^2+(d2*0.5-t2)*(d1*0.5-t1));

        
%         a = (n1*6-5): (n1*6);
%         b = (n2*6-5): (n2*6);
% 
%         KT(a, a) = KT(a, a) + Ki(1:6, 1:6);
%         KT(a, b) = KT(a, b) + Ki(1:6, 7:12);
%         KT(b, a) = KT(b, a) + Ki(7:12, 1:6);
%         KT(b, b) = KT(b, b) + Ki(7:12, 7:12);
% 
% 
%         MT(a, a) = MT(a, a) + Mi(1:6, 1:6);
%         MT(a, b) = MT(a, b) + Mi(1:6, 7:12);
%         MT(b, a) = MT(b, a) + Mi(7:12, 1:6);
%         MT(b, b) = MT(b, b) + Mi(7:12, 7:12);
        
    elseif(ElemType == 2)    % two node tapered element
    
        R1 = 0.5*d1;
        r1 = R1-t1;
        R2 = 0.5*d2;
        r2 = R2-t2;

        [a0, a1, a2, b0, b1, b2, b3, b4] = CrsSecCoeff(R1, r1, R2, r2, L, pi);

        [AA, BB, CC, DD, EE, FF, GG, HH, PP ] ...
         = MatrixCoeff(a0, a1, a2, b0, b1, b2, b3, b4, E, G, kappa, rho);

     
         Ki = ElemK1(AA, BB, CC, DD, EE, FF, L, DirCos);
         Mi = ElemM1(GG, HH, PP, L, DirCos);       
 
         A = pi*( r1^2 - r2^2 );

         Gi = ElemG(A, L, rho, DirCos);
    
    elseif(ElemType == 3)   % three node tapered element 
    
        R1 = 0.5*d1;
        r1 = R1-t1;
        R2 = 0.5*d2;
        r2 = R2-t2;

        [a0, a1, a2, b0, b1, b2, b3, b4] = CrsSecCoeff(R1, r1, R2, r2, L, pi);

        [AA, BB, CC, DD, EE, FF, GG, HH, PP ] ...
         = MatrixCoeff(a0, a1, a2, b0, b1, b2, b3, b4, E, G, kappa, rho);
 
         Ki = ElemK2(AA, BB, CC, DD, EE, FF, L, DirCos);
         Mi = ElemM2(GG, HH, PP, L, DirCos);       

    end % end if element type

    for j = 1:nnie % assemble to the system K and M
        
        jn = nn(j);
        ind_j = jn*6-5: jn*6;
        je = j*6-5 : j*6;
        
        Gv(ind_j) = Gv(ind_j) + Gi(je);
        
        for k = 1:nnie
            kn = nn(k);
            ind_k = kn*6-5:kn*6;
            ke = k*6-5:k*6;
            
            KT(ind_j, ind_k) = KT(ind_j, ind_k) + Ki(je, ke); 
            MT(ind_j, ind_k) = MT(ind_j, ind_k) + Mi(je, ke); 

        end

    end
    
    % total mass
    
%    TotalMass = TotalMass + ElemMass;
%     
%    %if (i>=2 && i<= 3 )
%    % if ( i>=54 && i<= 83 )     % OC3 Tripod
%    % if( ( i>=130 && i<=209 ) ) % OC4 jacket
%     if( ( i>=109 && i<=198 ) ) % MMI
%         TowerMass = TowerMass + ElemMass;
%     else
%         JacketMass = JacketMass + ElemMass;
%     end
%     
%     if ( p1 == 95 )
%         TPMass = TPMass + ElemMass;
%     end
    
%      if ( p1 == 5 || p1 == 6 ) % OC4 jacket
%          PileMass = PileMass + ElemMass;
%      end
    
    %-----------------------------------------------------------------
    for j = 1:nnie % assemble to system K and M in compressed row format
                   % for two-node element
        
        jn = nn(j);
        ind_j = jn*6-5: jn*6; % target rows
        je = j*6-5 : j*6;     % rows in element K
        for k = 1:nnie
            kn = nn(k);
            ind_k = kn*6-5:kn*6; % target columns
            ke = k*6-5:k*6;      % columns in element K
            for iin = 1:6 % 
                iim = je(iin);   % index i in Ke
                ti = ind_j(iin); % target i in KT
                for jjn = 1:6
                    jjm = ke(jjn);   % index j in Ke
                    tj = ind_k(jjn); % target j in KT

                    if (ti<=tj)
                        beg_jA = iK(ti); % find KT(i, j) index in KCR
                        end_jA = iK(ti+1) - 1;

                        tjCR = 0;
                        for ri = beg_jA:end_jA
                            if(jK(ri) == tj)
                                tjCR = ri;
                            end
                        end

                        if (tjCR == 0)
                           display(['K(',num2str(ti),',',num2str(tj), ') not found!']);
                        end
                        
                        KCR(tjCR) = KCR(tjCR) + 0.5*(Ki(iim, jjm) + Ki(jjm, iim));
                        MCR(tjCR) = MCR(tjCR) + Mi(iim, jjm);

                    end
                end

            end
                
          

        end

    end % end assemble KCR, MCR
    %-----------------------------------------------------------------

    
  
end



disp(' ')
disp('Total mass: ');
disp(TotalMass);

disp(' ')
disp('Tower: ');
disp(TowerMass);

disp(' ')
disp('TP: ');
disp(TPMass - 8786.035-4571.653);
%disp(TPMass - 17425.48705);
%disp(TPMass);

disp(' ')
disp('Pile: ');
disp(PileMass);


disp(' ')
disp('Jacket: ');
%disp(JacketMass - 666e3 - PileMass);
disp(TotalMass - TowerMass);

% add point mass 

for i = 1: NPM
    n = PMs(i, 1);
    nn = 6*n-5;
    
    MT(nn  , nn  ) = MT(nn  , nn  ) + PMs(i, 2);
    MT(nn+1, nn+1) = MT(nn+1, nn+1) + PMs(i, 3);
    MT(nn+2, nn+2) = MT(nn+2, nn+2) + PMs(i, 4);
    MT(nn+3, nn+3) = MT(nn+3, nn+3) + PMs(i, 5);
    MT(nn+4, nn+4) = MT(nn+4, nn+4) + PMs(i, 6);
    MT(nn+5, nn+5) = MT(nn+5, nn+5) + PMs(i, 7);
    
    
end % end point mass

% add point mass to MCR
for i = 1: NPM
    n = PMs(i, 1);
    ind_i = n*6-5:n*6;
    
    for j = 1:6
        tj = ind_i(j); % target row/colum in MT
        
        tjCR = iK(tj);
        MCR(tjCR) = MCR(tjCR) + PMs(i, j+1);
    end
    
end % end point mass


end

