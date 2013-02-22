function Gi = ElemG(A, L, rho, DirCos)
%GRAVITYFORCE Summary of this function goes here
%   Detailed explanation goes here

global g;

Gi = zeros(12, 1);

Gi(3) = 0.5*L*rho*A*g;
Gi(9) = Gi(3);

a = 1.0/12.0*g*L^2*rho*A;

% t3*t4-t1*t6
Gi(4) = a*( DirCos(1, 3)*DirCos(2, 1) - DirCos(1, 1)*DirCos(2, 3) ) ;
Gi(5) = a*( DirCos(1, 3)*DirCos(2, 2) - DirCos(1, 2)*DirCos(2, 3) ) ;

Gi(10) = -Gi(4);
Gi(11) = -Gi(5);


end

