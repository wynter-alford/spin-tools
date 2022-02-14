global knownOmegas knownPs knownQs %#ok<NUSED>

if strcmp(mode,'time')
    done = false;
    tic
    mt = 1;
    while ~done
        MagnusTerms{mt} = (1i/tCyc)*Omega(mt,toggledHsys,Taus);
        raw_hsizes(d,c,mt) = matOrder(MagnusTerms{mt})/hsys;
        elapsed = toc;
        
        if elapsed > computationTime
            done = true;
            maxTerm = mt - 1;
            mode = 'max'; % ensures the same number of terms are computed for each Hdip, test param
        else
            mt = mt + 1;
        end
    end
    
elseif strcmp(mode,'max')
    for mt=1:maxTerm+1
        MagnusTerms{mt} = (1i/tCyc)*Omega(mt,toggledHsys,Taus);
        raw_hsizes(d,c,mt) = matOrder(MagnusTerms{mt})/hsys;
    end   
end
