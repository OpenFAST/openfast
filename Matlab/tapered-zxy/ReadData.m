% read data file

function [nnNNode, newNElem, kprp, nnNodes, nnElems, NewProps, ...
    PMs, NPM, NFc, Fcs] = ReadData(filename, Div, ElemType)

% nnNNode, newNElem, pprp, ElemType



% read joints and members data
fpr = fopen(filename, 'r');

NumData = fscanf(fpr, '%g %g %g',[1, 5]);
NNode = NumData(1);
NProp = NumData(3);
NElem = NumData(2);
NPM = NumData(4);
NFc = NumData(5);

Nodes = fscanf(fpr, '%g %g %g %g',[4, NNode]);
Nodes = Nodes';

Elems = fscanf(fpr, '%g %g %g %g %g',[5, NElem]);
Elems = Elems';

Props = fscanf(fpr, '%g %g %g %g %g %g',[6, NProp]);
Props = Props';

PMs = fscanf(fpr, '%g %g %g %g %g %g %g ',[7, NPM]);
PMs = PMs';

Fcs = fscanf(fpr, '%g %g %g %g %g %g %g ',[7, NFc]);
Fcs = Fcs';

% MassProp = fscanf(fpr, '%g %g %g',[3, NMass]);
% MassProp = MassProp';


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


% change the node number in point mass with new numbering
for i = 1: NPM
    node = PMs(i, 1);
    
    for j = 1: NNode
        if (node == OldNodeNo(j) )
            PMs(i, 1) = NewNodeNo(j);
        end
    end
end


% change the node number in force with new numbering
for i = 1: NFc
    node = Fcs(i, 1);
    
    for j = 1: NNode
        if (node == OldNodeNo(j) )
            Fcs(i, 1) = NewNodeNo(j);
        end
    end
end




fclose(fpr);

% discretize the element using same divisions for every member
if(Div > 1)
    newNNode = NNode + (Div - 1)*NElem;
    newNElem = NElem*Div; 
    newNProp = newNElem*2;

    NewNodes = zeros(newNNode, 4);
    NewElems = zeros(newNElem, 5);
    NewProps = zeros(newNProp, 6);

    NewNodes(1:NNode, 1:4) = Nodes;
    NewElems(1:NElem, 1:5) = Elems;
    NewProps(1:NProp, 1:6) = Props;

    kp = NNode;
    ke = 0;
    kprp = NProp;

    for i = 1:NElem
        % creating new nodes
        node1 = Elems(i, 2);
        node2 = Elems(i, 3);
        prop1 = Elems(i, 4);
        prop2 = Elems(i, 5);

        x1 = Nodes(node1, 2);
        y1 = Nodes(node1, 3);
        z1 = Nodes(node1, 4);

        x2 = Nodes(node2, 2);
        y2 = Nodes(node2, 3);
        z2 = Nodes(node2, 4);

        dx = (x2 - x1)/Div;
        dy = (y2 - y1)/Div;
        dz = (z2 - z1)/Div;

        d1 = Props(Elems(i, 4), 5);
        t1 = Props(Elems(i, 4), 6);

        d2 = Props(Elems(i, 5), 5);
        t2 = Props(Elems(i, 5), 6);

        dd = (d2 - d1)/Div;
        dt = (t2 - t1)/Div;

        % first node/element connect to node1
        kp = kp+1;
        NewNodes(kp, 1) = kp;
        NewNodes(kp, 2) = x1 + dx;
        NewNodes(kp, 3) = y1 + dy;
        NewNodes(kp, 4) = z1 + dz;

        if(dd ~=0 || dt ~=0)
            kprp = kprp + 1;
            NewProps(kprp, 1) = kprp;
            NewProps(kprp, 2) = Props(Elems(i, 4), 2);
            NewProps(kprp, 3) = Props(Elems(i, 4), 3);
            NewProps(kprp, 4) = Props(Elems(i, 4), 4);
            NewProps(kprp, 5) = d1 + dd;
            NewProps(kprp, 6) = t1 + dt;

            nprp = kprp;
        else
            nprp = prop1;
        end

        ke = ke+1;
        NewElems(ke, 1) = ke;
        NewElems(ke, 2) = node1;
        NewElems(ke, 3) = kp;
        NewElems(ke, 4) = prop1;
        NewElems(ke, 5) = nprp;

        pprp = nprp;


        % interior nodes and elements

        for j = 2: Div-1
            kp = kp+1;
            NewNodes(kp, 1) = kp;
            NewNodes(kp, 2) = x1 + j*dx;
            NewNodes(kp, 3) = y1 + j*dy;
            NewNodes(kp, 4) = z1 + j*dz;

            if(dd ~=0 || dt ~=0)
                kprp = kprp + 1;
                NewProps(kprp, 1) = kprp;
                NewProps(kprp, 2) = Props(Elems(i, 4), 2);
                NewProps(kprp, 3) = Props(Elems(i, 4), 3);
                NewProps(kprp, 4) = Props(Elems(i, 4), 4);
                NewProps(kprp, 5) = d1 + j*dd;
                NewProps(kprp, 6) = t1 + j*dt;

                nprp = kprp;
            else
                nprp = Elems(i, 4);
            end

            ke = ke+1;
            NewElems(ke, 1) = ke;
            NewElems(ke, 2) = kp - 1;
            NewElems(ke, 3) = kp;
            NewElems(ke, 4) = pprp;
            NewElems(ke, 5) = nprp;

            pprp = nprp;
        end

        % last element: one connect to node2
        ke = ke+1;
        NewElems(ke, 1) = ke;
        NewElems(ke, 2) = kp;
        NewElems(ke, 3) = node2;
        NewElems(ke, 4) = nprp;
        NewElems(ke, 5) = prop2;

        pprp = nprp;

    end
else
    newNNode = NNode;
    newNElem = NElem; 
    kprp = NProp;

    NewNodes(1:NNode, 1:4) = Nodes;
    NewElems(1:NElem, 1:5) = Elems;
    NewProps(1:NProp, 1:6) = Props;


end % if


% for higher order elements, add nodes to elements
if(ElemType > 2)
    nne = ElemType -1;
    nnNNode = newNNode + (nne - 1)*newNElem;

    nnNodes = zeros(nnNNode, 4);
    nnElems = zeros(newNElem, 5+nne-1);

    nnNodes(1:newNNode, 1:4) = NewNodes;
    nnElems(1:newNElem, 1:2) = NewElems(1:newNElem, 1:2);
    nnElems(1:newNElem, 5+nne-3:5+nne-1) = NewElems(1:newNElem, 3:5);

    kp = newNNode;
    
    for i = 1:newNElem
        node1 = NewElems(i, 2);
        node2 = NewElems(i, 3);

        x1 = NewNodes(node1, 2);
        y1 = NewNodes(node1, 3);
        z1 = NewNodes(node1, 4);

        x2 = NewNodes(node2, 2);
        y2 = NewNodes(node2, 3);
        z2 = NewNodes(node2, 4);

        dx = (x2 - x1)/nne;
        dy = (y2 - y1)/nne;
        dz = (z2 - z1)/nne;
        
         for j = 1: nne-1
            kp = kp+1;
            nnNodes(kp, 1) = kp;
            nnNodes(kp, 2) = x1 + j*dx;
            nnNodes(kp, 3) = y1 + j*dy;
            nnNodes(kp, 4) = z1 + j*dz;

            nnElems(i, 2+j) = kp;
         end       

    end
    
else
    nne = ElemType -1;
    nnNNode = newNNode;
    nnNodes = NewNodes;
    nnElems = NewElems;
end



filename = 'newdata.txt';
fpr = fopen(filename, 'w');
fprintf(fpr, '%d %d %d %d\n', nnNNode, newNElem, kprp, ElemType);
for i = 1: nnNNode
    fprintf(fpr, '%d %15.7f %15.7f %15.7f \n', nnNodes(i, 1), nnNodes(i, 2:4));
end
for i = 1: newNElem
    if(ElemType == 1)
        fprintf(fpr, '%d %d %d %d %d   ', nnElems(i, 1:5+ElemType-1));
        fprintf(fpr, '\n');
    end
    
    if(ElemType == 2)
        fprintf(fpr, '%d %d %d %d %d   ', nnElems(i, 1:5+nne-1));
        fprintf(fpr, '\n');
    end
    
    if(ElemType == 3)
        fprintf(fpr, '%d %d %d %d %d %d ', nnElems(i, 1:5+nne-1));
        fprintf(fpr, '\n');
    end
    
end
for i = 1: kprp
    fprintf(fpr, '%d %15.7e %15.7e %15.7e %15.7f %15.7f \n', NewProps(i, 1:6));
end

for i = 1: NPM 
    fprintf(fpr, '%d %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e \n', PMs(i, 1:7));
end


fclose(fpr);


 PlotStructure(filename)


end % end function
