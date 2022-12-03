%% getSeqFrames.m
% Wynter Alford
% November 2022
%
% For a chosen pulse sequence, returns the toggled frames of the Iz operator
% in terms of 1-spin Pauli matrices.

function framesList = getSeqFrames(sequenceName,x,y,z)
    fseq = getSequence(sequenceName,x,y);

    UDxx=expm(-1i*x*pi/2);
    UDyy=expm(-1i*y*pi/2);
    UDxxbar=expm(1i*x*pi/2);
    UDyybar=expm(1i*y*pi/2);

    framesList = cell(1,length(fseq.Pulses)+1);
    for p = 0:length(fseq.Pulses)
        toggle = getURF(p,fseq.Pulses,2,x,y,UDxx,UDyy,UDxxbar,UDyybar);
        framesList{p+1} = toggle'*z*toggle; 
    end
end
   