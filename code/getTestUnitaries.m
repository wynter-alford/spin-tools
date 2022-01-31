function testUnitaries = getTestUnitaries(sequence,Hsys,pulse,tau,f1)
    global dim 
    
    UtauCell = {speye(dim,dim),expm(-1i*Hsys*2*pi*(tau-pulse)), expm(-1i*Hsys*2*pi*(2*tau-pulse))};
        
    testUnitary = speye(dim,dim);
    deltaUnitary = speye(dim,dim);
    
    for p=1:length(sequence.Pulses)
        Utau = UtauCell{sequence.Taus(p)+1};
        nextUs = getNextUs(sequence,Hsys,pulse,f1,p);
        nextU = nextUs{1};
        nextUD = nextUs{2};
        testUnitary = nextU * Utau * testUnitary;
        deltaUnitary = nextUD *Utau * deltaUnitary;
    end
    testUnitary = UtauCell{sequence.Taus(length(sequence.Taus))+1} * testUnitary;
    deltaUnitary = UtauCell{sequence.Taus(length(sequence.Taus))+1} * deltaUnitary;
    
    testUnitaries = {testUnitary,deltaUnitary};
end

function nextUs = getNextUs(sequence,Hsys,pulse,f1,p)
    nextU = expm(-1i*2*pi*(Hsys+f1*sequence.Pulses{p})*pulse);
    nextUD = expm(-1i*Hsys*2*pi*pulse/2)*expm(-1i*pi*sequence.Pulses{p}/2)*expm(-1i*Hsys*2*pi*pulse/2);  
    nextUs = {nextU,nextUD};
end