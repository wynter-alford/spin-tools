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
    magTermInd = 1;
    while ~done
        MagnusTerms{magTermInd} = (1i/tCyc)*Omega(magTermInd,toggledHsys,Taus); %#ok<*SAGROW> 
        raw_MHHT(paramValInd,coupMatInd,magTermInd) = log10(specnorm(MagnusTerms{magTermInd}));
        elapsed = toc;
        magTermInd %#ok<NOPTS>
        if elapsed > computationTime/min(pulseDivs,2)
            done = true;
            maxTerm = magTermInd - 1;
            mode = 'max'; % ensures the same number of terms are computed for each Hdip, test param
        else
            magTermInd = magTermInd + 1;
        end
    end
    
elseif strcmp(mode,'max')
    for magTermInd=1:maxTerm+1
        MagnusTerms{magTermInd} = (1i/tCyc)*Omega(magTermInd,toggledHsys,Taus);
        raw_MHHT(paramValInd,coupMatInd,magTermInd) = log10(specnorm(MagnusTerms{magTermInd}));
    end   
end

% Finite Pulse AHT
knownOmegas = {};
knownPs = {};
knownQs = {};
for magTermInd=1:maxTerm+1
    if pulseDivs > 1
        MagnusTermsP{magTermInd} = (1i/tCyc)*Omega(magTermInd,toggledHsysP,TausP);
        raw_MHHP(paramValInd,coupMatInd,magTermInd) = log10(specnorm(MagnusTermsP{magTermInd}));
    else
        MagnusTermsP{magTermInd} = nan;
        raw_MHHP(paramValInd,coupMatInd,magTermInd) = nan;
    end
end  