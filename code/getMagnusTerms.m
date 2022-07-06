%% getMagnusTerms.m
% Wynter Alford
% February 2022
%
% Computes, using the recursive method, the terms of the Magnus expansion
% for a given pulse sequence.  If the 'time' mode is chosen, keeps
% computing additional terms until some time limit is exceeded.  If the
% 'max' mode is chosen, computes up to some specified maxTerm.

global knownOmegas knownPs knownQs %#ok<NUSED>

if strcmp(mode,'time')
    done = false;
    tic
    mt = 1;
    while ~done
        MagnusTerms{mt} = (1i/tCyc)*Omega(mt,toggledHsys,Taus); %#ok<*SAGROW> 
        raw_hsizes(d,c,mt) = log10(specnorm(MagnusTerms{mt}));
        elapsed = toc;
        mt
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
        raw_hsizes(d,c,mt) = log10(specnorm(MagnusTerms{mt}));
    end   
end
