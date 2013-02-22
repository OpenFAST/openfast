function [a0, a1, a2, b0, b1, b2, b3, b4] = CrsSecCoeff(R1, r1, R2, r2, L, pi)
%CRSSECCOEFF Summary of this function goes here
%   Detailed explanation goes here

a0 = R1^2-r1^2;
a1 = 2/L*(r1^2 -r1*r2+R1*(R2-R1));
a2 = 1/L^2*( (R1-R2)^2 - (r1-r2)^2 );

b0 = -r1^4+R1^4;
b1 = ( 4*( r1^4-r1^3*r2+R1^3*(-R1+R2) ) )/L;
b2 = -((6*(r1^4-2*r1^3*r2+r1^2*r2^2-R1^2*(R1-R2)^2))/L^2);
b3 = (4*(r1^4-3*r1^3*r2+3*r1^2*r2^2-r1*r2^3-R1*(R1-R2)^3))/L^3;
b4 = (-r1^4+R1^4+4*r1^3*r2-6*r1^2*r2^2+4*r1*r2^3-r2^4-4*R1^3*R2+6*R1^2*R2^2-4*R1*R2^3+R2^4)/L^4;

a0 = a0*pi;
a1 = a1*pi;
a2 = a2*pi;

b0 = b0*pi*0.25;
b1 = b1*pi*0.25;
b2 = b2*pi*0.25;
b3 = b3*pi*0.25;
b4 = b4*pi*0.25;


end

