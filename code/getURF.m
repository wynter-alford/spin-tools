%% getURF.m
% Wynter Alford
% November 2020
%
% Calculates the toggling unitary for a pulse sequence after the first
% 'frame' pulses.  Assumes pulses are instantaneous, as is standard in AHT
% computations for pulse sequences.
%
% Make sure that the Ux, etc passed are for instantaneous pulses


function URF = getURF(frame,Pulses,X,Y,Ux,Uy,Uxbar,Uybar) 
    global dim
    
    if frame < 1
        URF = speye(dim,dim); %returns the identity if frame == 0
    elseif Pulses{1}==X
            URF = Ux;
    elseif Pulses{1}==Y
        URF = Uy;
    elseif Pulses{1}==-X
        URF = Uxbar;
    elseif Pulses{1}==-Y
        URF = Uybar;
    end
    
    if frame > 1
        for j=2:frame            
            if Pulses{j}==X
                    URF = Ux * URF;
            elseif Pulses{j}==Y
                URF = Uy * URF;
            elseif Pulses{j}==-X
                URF = Uxbar * URF;
            elseif Pulses{j}==-Y
                URF = Uybar * URF;
            end
        end
    end
end
