% plot the jacket structure
%function PlotStructure

function Plot3DStructure(file)

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




end % end function













