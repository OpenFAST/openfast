function [JointsConnE, JointsConnN] = JointConnect( NNode, NElem, Nodes, Elems )
%JOINTCONNECT Summary of this function goes here
%   Detailed explanation goes here

%----------------  Joint connectivity ------------------------------------
MaxMemberAtOneJoint = 10;
JointsConnE = zeros(NNode, MaxMemberAtOneJoint);
JointsConnN = zeros(NNode, MaxMemberAtOneJoint);

for i = 1: NNode
    
    JointsConnE(i, 1) = Nodes(i, 1); 
    JointsConnN(i, 1) = Nodes(i, 1); 
    
    k = 0;
    for j = 1: NElem
        if ( JointsConnE(i, 1)== Elems(j, 2)||JointsConnE(i, 1)== Elems(j, 3) )
            k = k + 1;
            JointsConnE(i, k+2) = Elems(j, 1);
            if (JointsConnE(i, 1)== Elems(j, 2))
                JointsConnN(i, k+2) = Elems(j, 3);
            else
                JointsConnN(i, k+2) = Elems(j, 2);
            end
        end
    end
    JointsConnE(i, 2) = k;
    JointsConnN(i, 2) = k;
end



end

