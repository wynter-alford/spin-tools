% Ben Alford, November 2020
%
% Calculates the toggling unitary for a pulse sequence after the first
% 'frame' pulses.  Assumes pulses have zero length. See the script magnus.m
% for example usage.

function URF = getURF(frame)
    global Pulses dim
    
    if frame < 1
        URF = speye(dim,dim); %returns the identity if frame == 0

    else
        URF = expm(-1i*Pulses{1}*pi/2);
    end
    
    if frame > 1
        for j=2:frame
            URF = expm(-1i*Pulses{j}*pi/2) * URF;
        end
    end
end