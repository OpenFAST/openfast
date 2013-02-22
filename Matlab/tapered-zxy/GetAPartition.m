function [ Ap, IAp, JAp, NZAp ] = GetAPartition( A, IA, JA, NA, NZA, Iarr, NIa, Jarr, NJa )
%GETAPARTITION Summary of this function goes here
%   Detailed explanation goes here

NZAp = 0;
Ap = zeros(NIa*NJa,1);
IAp = zeros(Nia+1,1) ;
JAp = zeros(NIa*NJa,1);

for i = 1: NIa
    ia = Iarr(i);
    
    j1 = JA( IA(ia) );
    j2 = JA( IA(ia+1) - 1 );
    
    for j = 1: NJa
       ja = Jarr(j);
       
        for k = j1: j2
            
            
        end
        
    end
    
end



end

