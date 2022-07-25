%% getSequence.m
% Wynter Alford
% Created February 2022; updated frequently
%
% Add sequences to this file following the format already present to be
% able to test new sequences. The pulse array just needs to list the pulses
% in order.  The Taus array lists the number of tau delays BEFORE each
% pulse.
%
% For example, the sequence WHH is:
% tau, X, tau, -Y, 2*tau, Y, tau, -X, tau
% This gives a Taus array of [1 1 2 1 1] and Pulses array {X, -Y, Y, -X}
%
% For sequences like the YXX sequences that start with a pulse, Taus(1)=0

function sequence = getSequence(sequenceName, X, Y)

    %WAHUHA
    if strcmp(sequenceName, 'WHH')
        sequence.Pulses = {X, -Y, Y, -X};
        sequence.Taus = [1 1 2 1 1];

    %MREV-8
    elseif strcmp(sequenceName, 'MREV8')
        sequence.Pulses = {-X, -Y, Y, X, X, -Y, Y, -X};
        sequence.Taus = [1 1 2 1 2 1 2 1 1];

    %BR-24    
    elseif strcmp(sequenceName, 'BR24')
        sequence.Pulses = {X, -Y, Y, -X, -X, -Y, Y, X, -Y, X, -X, Y, Y, X, -Y, X,    -X, Y, Y, X,   -X, -Y, -X, -Y };
        sequence.Taus = [1 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 1];
         
    %CORY 48
    elseif strcmp(sequenceName, 'CORY48')
        sequence.Pulses = {X, Y, -X, Y, X, Y, X, Y, X, -Y, X, Y, -Y, -X, Y, -X, -Y, -X, -Y, -X, -Y, X, -Y, -X, -X, Y, -X, -Y, -X, Y, X, -Y, -X, -Y, X, -Y, Y, -X, Y, X, Y, -X, -Y, X, Y, X, -Y, X};
        sequence.Taus = [1 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 1];
    
    %YXX-48
    elseif strcmp(sequenceName,'YXX48')
        sequence.Pulses = {Y, -X,-X,Y,-X,-X,-Y,X,X,Y,-X,-X,-Y,X,X,-Y,X,X,Y,-X,-X,Y,-X,-X,-Y,X,X,Y,-X,-X,-Y,X,X,-Y,X,X,Y,-X,-X,-Y,X,X,Y,-X,-X,-Y,X,X};
        sequence.Taus = ones(49,1);
        sequence.Taus(1)=0;
    
    %YXX-24
    elseif strcmp(sequenceName,'YXX24')
        sequence.Pulses = {-Y,X,-X,Y,-X,-X,Y,-X,X,-Y,X,X,Y,-X,X,-Y,X,X,-Y,X,-X,Y,-X,-X};
        sequence.Taus = ones(25,1);
        sequence.Taus(1)=0;
    
    %AZ-48 (OLD VERSION)
    elseif strcmp(sequenceName,'AZ48_OLD')
        sequence.Pulses = {-X,Y,Y,X,Y,Y,-Y,X,X,-Y,X,X,Y,X,X,-Y,X,X,-Y,X,-Y,X,X,-Y,-X,-X,Y,Y,-X,Y,Y,-Y,X,-Y,-Y,X,-Y,X,X,-Y,X,X,-Y,-X,-X,-Y,-X,-X};
        sequence.Taus = ones(49,1);
        sequence.Taus(1)=0;
    
    %Symmetrized 48 from YXX-24
    elseif strcmp(sequenceName,'YXX24S')
        sequence.Pulses = {-Y,X,-X,Y,-X,-X,Y,-X,X,-Y,X,X,Y,-X,X,-Y,X,X,-Y,X,-X,Y,-X,-X,X,X,-Y,X,-X,Y,-X,-X,Y,-X,X,-Y,-X,-X,Y,-X,X,-Y,X,X,-Y,X,-X,Y};
        sequence.Taus = [0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
    
    elseif strcmp(sequenceName,'AS10M')
        sequence.Pulses = {-X,Y,2*Y,-Y,-X,X,Y,2*Y,-Y,X};
        sequence.Taus = ones(11,1);
        sequence.Taus(6) = 2;
        
    elseif strcmp(sequenceName,'AS10E')
        sequence.Pulses = {X, -Y, Y, -X, 2*Y, X, Y, -Y, -X, 2*Y}; %(added pi-pulse with 0tau after to fix basis)
        sequence.Taus = ones(11,1);
        sequence.Taus(11) = 0;
        
    elseif strcmp(sequenceName,'AS10E-2cycle')
        sequence.Pulses = {X, -Y, Y, -X, 2*Y, X, Y, -Y, -X, X, -Y, Y, -X, 2*Y, X, Y, -Y, -X};
        sequence.Taus = ones(19,1);
        
    elseif strcmp(sequenceName,'WHHWHH')
        sequence.Pulses = {X, -Y, Y, -X, X, -Y, Y, -X};
        sequence.Taus = [1 1 2 1 2 1 2 1 1];
            
    elseif strcmp(sequenceName,'ML10')
        sequence.Pulses = {X,X,Y,X,X,-X,-X,-Y,-X,-X};
        sequence.Taus = [1 1 1 1 1 2 1 1 1 1 1];
    
    elseif strcmp(sequenceName,'MG8')
        sequence.Pulses = {X,Y,Y,X,X,Y,Y,X};
        sequence.Taus = [1 1 2 1 2 1 2 1 1];
        
    elseif strcmp(sequenceName,'AZ48') % Designed by Will Kaufman; his best 48-pulse sequence
        sequence.Pulses = {X, X, -Y, X, X, -Y, -X, -Y, -Y, -X, -Y, -Y, -X, -Y, -Y, -X, -Y, -Y, Y, Y, Y, -X, Y, Y, -X, -X, Y, Y, -X, Y, Y, Y, -X, Y, Y, -X, X, -Y, -Y, X, -Y, -Y, -Y, -X, -X, -Y, -X, -X};
        sequence.Taus = ones(49,1);
        sequence.Taus(1)=0;
        
    elseif strcmp(sequenceName,'I24') % one of Owen's test sequences
        sequence.Pulses = {X, Y, Y, X, -Y, -X, -Y, -X, Y, X, Y, X, -Y, -X, -X, -Y, -X, Y, -X, -Y, X, Y, -X, Y};
        sequence.Taus = [1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1];
    
    elseif strcmp(sequenceName,'SED24') % Designed by Owen Eskandari; particularly robust against phase transients
        sequence.Pulses = {X, Y, -X, -Y, Y, -X, Y, -X, Y, X, -X, Y, X, Y, -X, -Y, Y, X, -Y, -X, -Y, -X, X, -Y};
        sequence.Taus = [1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1];
        
    elseif strcmp(sequenceName,'SED48') % Designed by Owen Eskandari; best overall 48-pulse sequence
        sequence.Pulses = {-X, -Y, -X, Y, -X, Y, -Y, -X, Y, -X, -Y, -X, X, Y, X, -Y, X, -Y, Y, X, -Y, X, Y, X, X, Y, Y, X, Y, X, -Y, X, -X, Y, -X, -Y, -X, -Y, -X, Y, -Y, X, -Y, -X, -Y, -X, X, Y};
        sequence.Taus = [1 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 1];

    elseif strcmp(sequenceName,'SED96') % Designed by Owen Eskandari; best overall 96-pulse sequence
        sequence.Pulses = {-X, Y, -Y, -X, Y, -X, X, Y, -Y, X, -X, -Y, X, -Y, X, -Y, Y, -X, -X, -Y, Y, X, Y, X, Y, X, Y, -X, -X, Y, -Y, -X, -Y, X, X, -Y, X, Y, -X, -Y, X, -Y, Y, -X, -Y, -X, X, Y, Y, X, X, Y, -X, Y, X, Y, Y, X, Y, -X, -Y, -X, X, -Y, -X, -Y, X, -Y, -Y, -X, -Y, X, Y, -X, X, Y, X, Y, -X, Y, -Y, -X, -Y, X, -Y, -X, -Y, X, -X, Y, -X, -Y, -X, -Y, Y, X};
        sequence.Taus =  [1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1];
        error("SED96 is currently not functional")
    
    elseif strcmp(sequenceName,'TestXY')
        sequence.Pulses = {X, Y};
        sequence.Taus = [1 1 1];
    else
        error("Invalid Sequence Name. See getSequence.m for all valid sequence names.")    
    
    end
end
