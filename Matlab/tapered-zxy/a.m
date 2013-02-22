clc;
clear all;

pi = acos(-1);

filename = 'C:\Documents and Settings\hsong\My Documents\Visual Studio 2008\Projects\BeamFEM\BeamFEM\testframe_eigen_results.txt';
fpr = fopen(filename, 'r');

%temp = fscanf(fpr, '%g',[1,1]  );
K = fscanf(fpr, '%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g',[24, 24]);
%temp = fscanf(fpr, '%s' );
M = fscanf(fpr, '%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g',[24, 24]);


[V, D]= eig(K, M);

for i = 1:24
    w(i) = sqrt(D(i, i))/2/pi;
end

w = sort(w);

EI = 210e9*pi*0.25*(1.5^4-1.4^4);

m = 7850*pi*(1.5^2-1.4^2);

L = 45;

omega= 3.516*sqrt(EI/m/L^4);

fclose(fpr);