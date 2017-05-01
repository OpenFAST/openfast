% Test14.m
% Written by J. Jonkman, NREL
% Last update: 07/27/2016

% This m-file is used to call perform eigenanalysis using file 'Test14.lin' 
%   and echo out the average sorted natural frequencies.

% ----------- Call GetMats_f8.m and mbc3 using Test14.1.lin ---------------
FileNames = {'Test14.1.lin'};
GetMats_f8

mbc3

% ----------- Echo out the average sorted natural frequencies -------------
N
sort( MBC_NaturalFrequencyHz )


% ----------- Exit MATLAB -------------------------------------------------
exit;


