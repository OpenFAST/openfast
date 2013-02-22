% plot the jacket structure
%function PlotStructure

function PlotStructure(file)

global NFigure;

fpr = fopen(file, 'r');


% read joint and member data
NumData = fscanf(fpr, '%g %g %g %g',[1, 4]);
NNode = NumData(1);
NElem = NumData(2);
NProp = NumData(3);
ElemType = NumData(4);

Nodes = fscanf(fpr, '%g %g %g %g',[4, NNode]);
Nodes = Nodes';

if(ElemType == 1)
    nne = ElemType;
    Elems = fscanf(fpr, '%g %g %g %g %g',[5, NElem]);
end
if(ElemType == 2)
    nne = ElemType -1;
    Elems = fscanf(fpr, '%g %g %g %g %g ',[5, NElem]);
end
if(ElemType == 3)
    nne = ElemType -1;
    Elems = fscanf(fpr, '%g %g %g %g %g %g',[6, NElem]);
end

Elems = Elems';

Props = fscanf(fpr, '%g %g %g ',[3, NProp]);
Props = Props';

OldNodeNo = Nodes(:,1);
NewNodeNo = 1:1:NNode;
Nodes(:, 1) = NewNodeNo';   % renumber node number
OldElems = Elems;


% change the node number in elements with new numbering
for i = 1: NElem
    node1 = Elems(i, 2);
    node2 = Elems(i, 3);
    
    for j = 1: NNode
        if (node1 == OldNodeNo(j) )
            Elems(i, 2) = NewNodeNo(j);
        end
        
        if (node2 == OldNodeNo(j) )
            Elems(i, 3) = NewNodeNo(j);
        end
    end
end
Elems(:, 1) = [1:1:NElem]';


NFigure = NFigure + 1;
figure (NFigure);

%plot3(Nodes(SampleJoint,2),Nodes(SampleJoint,3),Nodes(SampleJoint,4),'*g', 'linewidth',7);
%hold on;
colors = zeros(NProp, 3);
dc = 1/(NProp-1);
colors(:, 1) = 0.5*1:-0.5*dc:0;
colors(:, 2) = 0:0.99*dc:0.99*1;
colors(:, 3) = 0.5*1:-0.5*dc:0;


for i = 1: NElem
    
    for j = 1:nne+1
        X(j) = Nodes(Elems(i, 1+j), 2);
        Y(j) = Nodes(Elems(i, 1+j), 3);
        Z(j) = Nodes(Elems(i, 1+j), 4);
        
    end
    
    propn = Elems(i, 1+nne+1+1);
%     node1 = Elems(i, 2);
%     node2 = Elems(i, 3);
%     
%     propn = Elems(i, 4);
    
%     for j = 1: NNode
%         nodenumber = Nodes(j, 1);
%         if(nodenumber == node1)
%             x_node1 = Nodes(j, 2);
%             y_node1 = Nodes(j, 3);
%             z_node1 = Nodes(j, 4);
%         end 
%         if(nodenumber == node2)
%             x_node2 = Nodes(j, 2);
%             y_node2 = Nodes(j, 3);
%             z_node2 = Nodes(j, 4);
%         end 
%         
%     end


%     
%         X = [x_node1, x_node2];
%         Y = [y_node1, y_node2];
%         Z = [z_node1, z_node2];
%  




%        if(propn == 1)
%          plot3(X, Y, Z, '.-y', 'linewidth',1.5);
%        
%        elseif(propn == 2)
%          plot3(X, Y, Z, '.-k', 'linewidth',1.5);
%        
%        elseif(propn == 3)
%          plot3(X, Y, Z, '.-b', 'linewidth',1.5);
%          
%        elseif(propn == 4)
%          plot3(X, Y, Z, '.-m', 'linewidth',1.5);
%        
% %        elseif(propn == 5)
% %          plot3(X, Y, Z, '.-r', 'linewidth',1.5);
% %          
% %        elseif(propn == 6)
% %          plot3(X, Y, Z, '.-g', 'linewidth',1.5);
% %          
% %        elseif(propn == 7)
% %          plot3(X, Y, Z, '.-c', 'linewidth',1.5);
%          
% %        elseif(propn == 8)
% %          plot3(X, Y, Z, '--k', 'linewidth',1.5);
% %          
% %        elseif(propn == 9)
% %          plot3(X, Y, Z, '--m', 'linewidth',1.5);
%          
%        else
%          plot3(X, Y, Z, '-.', 'linewidth',1.5);
%          %plot3(X, Y, Z, '.-', 'color', colors(propn, 1:3), 'linewidth',1.5);
%        end

plot3(X, Y, Z, '.-k', 'linewidth',1.5);

        hold on;

end
axis off;

axis equal;
hold off;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% JointsConn = zeros(NNode, MaxMemberAtOneJoint);
% 
% for i = 1: NNode
%     
%     JointsConn(i, 1) = Nodes(i, 1); 
%     
%     k = 0;
%     for j = 1: NElem
%         if ( JointsConn(i, 1)== Elems(j, 2)||JointsConn(i, 1)== Elems(j, 3) )
%             k = k + 1;
%             JointsConn(i, k+2) = Elems(j, 1);
%         end
%     end
%     JointsConn(i, 2) = k;
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

% 
% figure (2)
% 
% 
% %plot3(Nodes(SampleJoint,2),Nodes(SampleJoint,3),Nodes(SampleJoint,4),'og', 'linewidth',5);
% %hold on;
% 
% for i = 1: JointsConn(SampleJoint, 2)
%     node1 = Elems(JointsConn(SampleJoint, i+2), 2);
%     node2 = Elems(JointsConn(SampleJoint, i+2), 3);
%     
%     if(node2 == SampleJoint)
%         node2 = node1;
%         node1 = SampleJoint;
%     end
% 
%     
%     x_node1 = Nodes(node1, 2);
%     y_node1 = Nodes(node1, 3);
%     z_node1 = Nodes(node1, 4);
%     
%     x_node2 = Nodes(node2, 2);
%     y_node2 = Nodes(node2, 3);
%     z_node2 = Nodes(node2, 4);
% 
%     X = [x_node1, x_node2];
%     Y = [y_node1, y_node2];
%     Z = [z_node1, z_node2];
%     
%     plot3(x_node1,y_node1,z_node1,'.g', 'linewidth',1.5);
%     hold on;
%     plot3(x_node2,y_node2,z_node2,'.m', 'linewidth',1.5);
%     hold on;
%     plot3(X, Y, Z, '-', 'linewidth',1);
%     axis equal;
%     hold on;
% end
% 
% axis off;
% hold off;
% 
% % 
% 




end % end function













