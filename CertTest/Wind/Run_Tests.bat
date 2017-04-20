set DriverProg=..\bin\AeroDyn_Driver_Win32.exe
set Compare=FC /T

%DriverProg% AOC_Test01.dvr
%DriverProg% AOC_Test02.dvr
%DriverProg% AOC_Test03.dvr

%Compare% .\AOC_Test01.1.out .\NREL_Results\AOC_Test01.1.out
%Compare% .\AOC_Test01.2.out .\NREL_Results\AOC_Test01.2.out
%Compare% .\AOC_Test01.3.out .\NREL_Results\AOC_Test01.3.out
%Compare% .\AOC_Test02.1.out .\NREL_Results\AOC_Test02.1.out
%Compare% .\AOC_Test02.2.out .\NREL_Results\AOC_Test02.2.out
%Compare% .\AOC_Test02.3.out .\NREL_Results\AOC_Test02.3.out
%Compare% .\AOC_Test03.1.out .\NREL_Results\AOC_Test03.1.out
%Compare% .\AOC_Test03.2.out .\NREL_Results\AOC_Test03.2.out
%Compare% .\AOC_Test03.3.out .\NREL_Results\AOC_Test03.3.out

:: move AOC_Test*.out .\NREL_Results



%DriverProg% NRELOffshrBsline5MW_Onshore_Test01.dvr
%DriverProg% NRELOffshrBsline5MW_Onshore_Test02.dvr
%DriverProg% NRELOffshrBsline5MW_Onshore_Test03.dvr

%Compare% .\NRELOffshrBsline5MW_Onshore_Test01.1.out .\NREL_Results\NRELOffshrBsline5MW_Onshore_Test01.1.out
%Compare% .\NRELOffshrBsline5MW_Onshore_Test01.2.out .\NREL_Results\NRELOffshrBsline5MW_Onshore_Test01.2.out
%Compare% .\NRELOffshrBsline5MW_Onshore_Test01.3.out .\NREL_Results\NRELOffshrBsline5MW_Onshore_Test01.3.out
%Compare% .\NRELOffshrBsline5MW_Onshore_Test02.1.out .\NREL_Results\NRELOffshrBsline5MW_Onshore_Test02.1.out
%Compare% .\NRELOffshrBsline5MW_Onshore_Test02.2.out .\NREL_Results\NRELOffshrBsline5MW_Onshore_Test02.2.out
%Compare% .\NRELOffshrBsline5MW_Onshore_Test02.3.out .\NREL_Results\NRELOffshrBsline5MW_Onshore_Test02.3.out
%Compare% .\NRELOffshrBsline5MW_Onshore_Test03.1.out .\NREL_Results\NRELOffshrBsline5MW_Onshore_Test03.1.out
%Compare% .\NRELOffshrBsline5MW_Onshore_Test03.2.out .\NREL_Results\NRELOffshrBsline5MW_Onshore_Test03.2.out
%Compare% .\NRELOffshrBsline5MW_Onshore_Test03.3.out .\NREL_Results\NRELOffshrBsline5MW_Onshore_Test03.3.out

:: move NRELOffshrBsline5MW_Onshore_Test*.out .\NREL_Results\

set DriverProg=