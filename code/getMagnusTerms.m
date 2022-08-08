%% getMagnusTerms.m
% Wynter Alford
% February 2022
%
% Computes, using the recursive method, the terms of the Magnus expansion
% for a given pulse sequence.  If the 'time' mode is chosen, keeps
% computing additional terms until some time limit is exceeded.  If the
% 'max' mode is chosen, computes up to some specified maxTerm.

global knownOmegas knownPs knownQs %#ok<GVMIS,NUSED>

knownOmegas = {};
knownPs = {};
knownQs = {};        
MagnusTerms = {};

if strcmp(mode,'time')
    done = false;
    tic
    termCalcInd = 1;
    while ~done
        MagnusTerms{termCalcInd} = (1i/tCyc)*Omega(termCalcInd,toggledHsys,Taus); %#ok<*SAGROW> 
        raw_HHT(paramValInd,coupMatInd,termCalcInd) = log10(specnorm(MagnusTerms{termCalcInd}));
        elapsed = toc;
        termCalcInd %#ok<NOPTS>
        if elapsed > computationTime/min(pulseDivs,2)
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
        raw_HHT(paramValInd,coupMatInd,termCalcInd) = log10(specnorm(MagnusTerms{termCalcInd}));
    end   
end

% Finite Pulse AHT
knownOmegas = {};
knownPs = {};
knownQs = {};
for termCalcInd=1:maxTerm+1
    if pulseDivs > 1
        MagnusTermsP{termCalcInd} = (1i/tCyc)*Omega(termCalcInd,toggledHsysP,TausP);
        raw_HHP(paramValInd,coupMatInd,termCalcInd) = log10(specnorm(MagnusTermsP{termCalcInd}));
    else
        MagnusTermsP{termCalcInd} = nan;
        raw_HHP(paramValInd,coupMatInd,termCalcInd) = nan;
    end
end  