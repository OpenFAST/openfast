function [Pitch,GenTq]=get_FastLinPitchTorque(fast_text);

%get index to lines starting with data
index=regexp(fast_text,'\n\s*\d+\s+\d+');
index=index(end);   %last line of data
fast_text=fast_text(index:end);

GenTq=sscanf(fast_text,'%f',5); %torqe is 5th number
GenTq=GenTq(end)*1000;
Pitch=sscanf(fast_text,'%f',6); %pitch is 6th number
Pitch=Pitch(end);
