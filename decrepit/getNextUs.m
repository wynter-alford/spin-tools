function nextUs = getNextUs(sequence,Hsys,pulse,f1,p)
    nextU = expm(-1i*2*pi*(Hsys+f1*sequence.Pulses{p})*pulse);
    nextUD = expm(-1i*Hsys*2*pi*pulse/2)*expm(-1i*pi*sequence.Pulses{p}/2)*expm(-1i*Hsys*2*pi*pulse/2);  
    nextUs = {nextU,nextUD};
end
