function sequence = getSequence(sequenceName, X, Y)
    %global X Y %#ok<GVMIS> 
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
    
    %AZ-48
    elseif strcmp(sequenceName,'AZ48')
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
    end
end
