function maxTerm = getMaxTerm(sequenceName, plus)
    seq = getSequence(sequenceName);
    ll = length(seq.Pulses);
    if ll<5
        maxTerm = 14;
    elseif ll<11
        maxTerm = 10;
    elseif ll<25
        maxTerm = 6;
    else
        maxTerm = 2;
    end

    if exist('plus','var')
        maxTerm = maxTerm + plus;
    end
end