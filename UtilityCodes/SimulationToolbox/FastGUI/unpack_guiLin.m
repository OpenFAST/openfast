function [A,B,Bd,C,D,Dd,RotorSpeed,HHWind,Pitch,GenTq...
	,LinearInputDesc,LinearOutputDesc,LinResults_text]=unpack_guiLin(guiLinFile)

try 
	load(guiLinFile); % gui_LinResults_W14_Rspd37p1_P12_Tq3519
catch
	disp('cannot find gui lin file:');
	disp(guiLinFile);
	error('cannot find gui lin file');
end;

RotorSpeed=LinResults.RotorSpeed;
HHWind=LinResults.HHWind;
Pitch=LinResults.trimPitch;
GenTq=LinResults.trimGenTq;

%% make sure i/o descriptions are available
if isfield(FastData,'LinResults_text')
	LinResults_text=FastData.LinResults_text;
else
	error('"LinResults_text" is not part of FastData')
end;
if isfield(FastData,'LinearInputDesc')
	LinearInputDesc=FastData.LinearInputDesc;
else
	error('"LinearInputDesc" is not part of FastData');
end;
if isfield(FastData,'LinearInputDesc')
	LinearInputDesc=FastData.LinearInputDesc;
else
	error('"LinearInputDesc" is not part of FastData');
end;
if isfield(FastData,'FastParams')
	LinearOutputDesc=FastData.FastParams.OutList;
	LinearOutputDesc=strrep(LinearOutputDesc,'"','');
	[n,m]=size(LinearOutputDesc);
	if m>n
		LinearOutputDesc=LinearOutputDesc';
	end;
else
	error('Could not find output desc in FastParams.');
end;
%% Unpack linearization matrices
A=LinResults.Amat;
B=LinResults.Bmat;
Bd=LinResults.Bdmat;
C=LinResults.Cmat;
D=LinResults.Dmat;
Dd=LinResults.Ddmat;
