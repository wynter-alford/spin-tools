global knownOmegas knownPs knownQs %#ok<NUSED>

if strcmp(mode,'time')
    done = false;
    tic
    mt = 1;
    while ~done
        MagnusTerms{mt} = (1i/tCyc)*Omega(mt,toggledHsys,Taus);
        trunc = trunc + MagnusTerms{mt};  
        raw_hsizes(d,c,mt) = matOrder(MagnusTerms{mt})/hsys;
        raw_C0(d,c,mt) = matOrder(comm(MagnusTerms{mt},MagnusTerms{1}))/hsys;
        raw_CS(d,c,mt) = matOrder(comm(MagnusTerms{mt},trunc))/hsys;

        mt
        elapsed = toc
        if elapsed > computationTime
            done = true;
        else
            mt = mt + 1;
        end
    end    
elseif strcmp(mode,'max')    
    for mt=1:maxTerm+1
        MagnusTerms{mt} = (1i/tCyc)*Omega(mt,toggledHsys,Taus);
        trunc = trunc + MagnusTerms{mt};
        raw_hsizes(d,c,mt) = matOrder(MagnusTerms{mt})/hsys;
    end   
end
