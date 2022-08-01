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
    termCalcInd = 1;
    while ~done
        MagnusTerms{termCalcInd} = (1i/tCyc)*Omega(termCalcInd,toggledHsys,Taus); %#ok<*SAGROW> 
        raw_hsizes(paramValInd,coupMatInd,termCalcInd) = log10(specnorm(MagnusTerms{termCalcInd}));
        elapsed = toc;
        termCalcInd %#ok<NOPTS>
        if elapsed > computationTime
            done = true;
            maxTerm = termCalcInd - 1;
            mode = 'max'; % ensures the same number of terms are computed for each Hdip, test param
        else
            termCalcInd = termCalcInd + 1;
        end
    end
    
elseif strcmp(mode,'max')
    for termCalcInd=1:maxTerm+1
        MagnusTerms{termCalcInd} = (1i/tCyc)*Omega(termCalcInd,toggledHsys,Taus);
        raw_hsizes(paramValInd,coupMatInd,termCalcInd) = log10(specnorm(MagnusTerms{termCalcInd}));
    end   
end
