function PlotModalResultsJackets(maxn, NElem, Nodes, Elems, V, omega, p, ElemType)


factors(1:maxn) = 10^p;
%factors(1:maxn) = 0;
for nn = 1: maxn
    
    

    MN = nn;   % mode number
    
    fct = factors(nn);
    
    
    figure(MN)

    for i = 1: NElem

          for j = 1:ElemType+1
                mm = Elems(i, 1+j);
                X(j) = Nodes(mm, 2) + fct*V(mm*6-5, nn);
                Y(j) = Nodes(mm, 3) + fct*V(mm*6-4, nn);
                Z(j) = Nodes(mm, 4) + fct*V(mm*6-3, nn);

          end
    


%         N0 = Elems(i, 2);
%         N1 = Elems(i, 3);
% 
% 
%         x0 = Nodes(N0, 2);
%         y0 = Nodes(N0, 3);
%         z0 = Nodes(N0, 4);
% 
%         dx0 = V(N0*6-5, nn);
%         dy0 = V(N0*6-4, nn);
%         dz0 = V(N0*6-3, nn);
% 
%         x0 = x0 + dx0*fct;
%         y0 = y0 + dy0*fct;
%         z0 = z0 + dz0*fct;
% 
%         x1 = Nodes(N1, 2);
%         y1 = Nodes(N1, 3);
%         z1 = Nodes(N1, 4);
% 
%         dx1 = V(N1*6-5, nn);
%         dy1 = V(N1*6-4, nn);
%         dz1 = V(N1*6-3, nn);
% 
%         x1 = x1 + dx1*fct;
%         y1 = y1 + dy1*fct;
%         z1 = z1 + dz1*fct;
% 
%         X = [x0 x1]; Y = [y0 y1]; Z = [z0, z1];
        plot3(X, Y, Z, '.-k', 'linewidth',1.5); 
        hold on;
    end
    axis off;
    title(['\fontsize{16}Mode # ',num2str(MN),': \omega = ',num2str(omega(MN))]);

    axis equal;
    hold off;
    
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
































