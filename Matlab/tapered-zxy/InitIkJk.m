function [ iK, jK, NNZ ] = InitIkJk( NNode, JointsConnN )
%INITIKJK Summary of this function goes here
%   Detailed explanation goes here
NNZ = 0;
iK = zeros(NNode*6+1,1);
for i=1: NNode
    bgn_col = i*6-5:i*6;
    NN = JointsConnN(i, 2);
    ad_col = zeros(1, NN*6);
    ad_node = 0;
    
    
    for j = 1: NN % column numbers ecxept the one same as the row number
        tn = JointsConnN(i, j+2);
        if(tn > i) % only count the node number which is larger than the joint node
            ad_node = ad_node + 1;
            jj = ad_node;
            ad_col(jj*6-5:jj*6) = tn*6-5:tn*6;
        end
    end
    
    total_num_col = 6 + ad_node*6;
    col_array = [bgn_col, ad_col];
    
    for r = 1:6
        iK((i-1)*6+r) = NNZ + 1;

        for k = r: total_num_col
            NNZ = NNZ+1;
            jK(NNZ) = col_array(k);
        end
    end
    
    
end

iK(NNode*6 +1) = iK(1)+ NNZ;
jK = jK';

end
