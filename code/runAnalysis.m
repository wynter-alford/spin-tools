%% runAnalysis.m
% Wynter Alford
% Original version from April 2020. Last updated July 2021. 
%
% Analyze pulse sequences in terms of the Magnus Expansion.  Compute terms
% in the Magnus Series.  Compute fidelity metrics describing how well the
% propagator from the computed Magnus terms resembles the actual
% experimental propagator, with either finite length or instantaneous
% pulses. New features added frequently as I continue to research this
% problem.

%% Configuration

% Derived parameters
dim = 2^N;
f1 = 1/(4*pulse); 
w1 = 2*pi*f1;
rotationError = 1+overRotation; %1 is a perfect pulse

% Decrepit
mode = 'max';
computationTime = 1;

%Configure Test Variable
testValueMax = getVarMax(testVarName,varMaxMod);

%% Initialize Quantum Setting

initVars            % pauli matrices
initCollectiveObs   % collective spin operators
initHamiltonians    % Hamiltonians
initResultStorage   % arrays for storing results

testVars = zeros(testValueCount,1);

%% Iterate over different parameter values, compute AHT terms and fidelities

for paramValInd=1:length(testVars)
    
    % Configure Test Variable and Test Sequence
    if strcmp(testVarName,'tau')||strcmp(testVarName,'Tau')
        tau = pulse + paramValInd*(testValueMax/testValueCount);
        testVars(paramValInd)=tau;
    elseif strcmp(testVarName,'tau_lo')
        tau = paramValInd*(testValueMax/testValueCount); % WARNING: May interact poorly with finite-pulse results
        testVars(paramValInd)=tau;
    elseif strcmp(testVarName,'Delta')||strcmp(testVarName,'delta')||strcmp(testVarName,'delta_lo')
        Delta = 2*(paramValInd-(testValueCount/2))*(testValueMax/testValueCount);
        testVars(paramValInd)=Delta;
    elseif strcmp(testVarName,'Coupling')||strcmp(testVarName,'coupling')||strcmp(testVarName,'coupling_lo')
        coupling = paramValInd*(testValueMax/testValueCount);
        testVars(paramValInd)=coupling;
    elseif strcmp(testVarName,'Overrotations')||strcmp(testVarName,'overrotations')||strcmp(testVarName,'overrot_hi')
        overRotation = 2*(paramValInd-(testValueCount/2))*(testValueMax/testValueCount);
        rotationError = 1+overRotation; % rotationError = 1 is a pulse with no rotation error. <1 or >1 is an under/over rotation.
        testVars(paramValInd) = overRotation;
    elseif strcmp(testVarName,'phaseTrans')
        phaseTrans = 2*(paramValInd-(testValueCount/2))*(testValueMax/testValueCount);
        testVars(paramValInd) = phaseTrans;
    elseif strcmp(testVarName,'WT1') % for omega_D * tau troubleshooting
        testVars(paramValInd) = coupling*tau;
    end
    
    sequence = getSequence(sequenceName,X,Y);   
    Pulses = sequence.Pulses;
    Taus = tau*sequence.Taus;        
    tCyc = sum(Taus);
    
    % For each dipolar coupling matrix
    for coupMatInd=1:couplingsCount  
        
        Hdip = Hdips{coupMatInd};
        Hsys = Hdip*coupling + Z*Delta;        
        hsys = specnorm(Hsys);
        
        % obtain experimental unitary
        initPulses
        getTestUnitaries
        
        % Set up "traditional" (HHT) toggled Hsys
        toggledHsys = cell(1,length(Pulses)+1);
        toggledHsysP = cell(1,length(Pulses)*pulseDivs+1);
        for p = 0:length(Pulses)
            toggle = getURF(p,Pulses,dim,X,Y,UDx,UDy,UDxbar,UDybar);
            toggledHsys{p+1} = toggle'*Hsys*toggle; 
        end

        % Set up finite-pulse (HHP) toggled Hsys
        if pulseDivs > 1
            PulsesP = repmat({ones(1)},1,length(Pulses)*pulseDivs);
            for semiPulseInd=1:length(PulsesP)
                PulsesP{semiPulseInd}=PulsesP{semiPulseInd}*Pulses{ceil(semiPulseInd/pulseDivs)};
            end
            for p = 0:length(Pulses)*pulseDivs
                toggle = getURF(p,PulsesP,dim,X,Y,UdivX,UdivY,UdivXbar,UdivYbar);
                toggledHsysP{p+1} = toggle'*Hsys*toggle;
            end

            TausP = Taus - (pulseDivs - 1) * pulse / pulseDivs;
            TausP(1) = tau - pulse/2;
            TausP(end) = tau - (pulse/2 - pulse/pulseDivs);
            tauPInd = 1;
            while tauPInd < length(TausP)  && pulseDivs > 1
                TausP = [TausP(1:tauPInd) (pulse/pulseDivs)*ones(1,pulseDivs-1) TausP(tauPInd+1:end)];
                tauPInd = tauPInd + pulseDivs; 
            end
        end
         
        if do_Magnus
            
            % Compute Magnus Series terms
            getMagnusTerms

            % obtain theoretical unitaries from AHT and compute Fidelity
            AHTUnitaries = {};
            AHTUnitariesP = {};
            sumAH = sparse(dim,dim);
            sumAHP = sparse(dim,dim);

            for termInd=1:maxTerm+1
                sumAH = sumAH + MagnusTerms{termInd};
                sumAHP = sumAHP + MagnusTermsP{termInd};
                AHTUnitaries{termInd} = expm(-1i*sumAH*2*pi*tCyc); %#ok<*SAGROW> 
                AHTUnitariesP{termInd} = expm(-1i*sumAHP*2*pi*tCyc);   


                raw_MOFT(paramValInd,coupMatInd,termInd)=overlap(expUnitary, AHTUnitaries{termInd}, N);  
                raw_MOIT(paramValInd,coupMatInd,termInd) = overlap(deltaUnitary, AHTUnitaries{termInd}, N);
                raw_MOFP(paramValInd,coupMatInd,termInd) = overlap(expUnitary, AHTUnitariesP{termInd}, N);

                raw_MDFT(paramValInd,coupMatInd,termInd)=uDist(expUnitary,AHTUnitaries{termInd},N);
                raw_MDIT(paramValInd,coupMatInd,termInd)=uDist(deltaUnitary,AHTUnitaries{termInd},N);
                raw_MDFP(paramValInd,coupMatInd,termInd) = overlap(expUnitary, AHTUnitariesP{termInd}, N);
            end

            % Time-suspension fidelities
            raw_SOFS(paramValInd,coupMatInd) = overlap(expUnitary,speye(dim,dim),N);
            raw_SOIS(paramValInd,coupMatInd) = overlap(deltaUnitary,speye(dim,dim),N); 

            raw_SDFS(paramValInd,coupMatInd) = uDist(expUnitary,speye(dim,dim),N);
            raw_SDIS(paramValInd,coupMatInd) = uDist(deltaUnitary,speye(dim,dim),N); 
        end
        
        if do_Floquet            
            getFloquetTerms % computes Floquet propagators and computes fidelities.           
        end
        
    end
    
    % Average Raw Results   
    for termAvgInd = 1:maxTerm+1
        results_MHHT(paramValInd,termAvgInd)=mean(raw_MHHT(paramValInd,:,termAvgInd));
        results_MHHP(paramValInd,termAvgInd)=mean(raw_MHHP(paramValInd,:,termAvgInd));
        results_MOFT(paramValInd,termAvgInd)=mean(raw_MOFT(paramValInd,:,termAvgInd));
        results_MOIT(paramValInd,termAvgInd)=mean(raw_MOIT(paramValInd,:,termAvgInd)); 
        results_MOFP(paramValInd,termAvgInd)=mean(raw_MOFP(paramValInd,:,termAvgInd)); 
        results_MDFT(paramValInd,termAvgInd)=mean(raw_MDFT(paramValInd,:,termAvgInd));
        results_MDIT(paramValInd,termAvgInd)=mean(raw_MDIT(paramValInd,:,termAvgInd));
        results_MDFP(paramValInd,termAvgInd)=mean(raw_MDFP(paramValInd,:,termAvgInd));
    end
    for termAvgInd = 1:maxFloqN+1
        results_QOCT(paramValInd,termAvgInd)=mean(raw_QOCT(paramValInd,:,termAvgInd));
        results_QDCT(paramValInd,termAvgInd)=mean(raw_QDCT(paramValInd,:,termAvgInd));
        results_QOCP(paramValInd,termAvgInd)=mean(raw_QOCP(paramValInd,:,termAvgInd));
        results_QDCP(paramValInd,termAvgInd)=mean(raw_QDCP(paramValInd,:,termAvgInd));
        results_QOFT(paramValInd,termAvgInd)=mean(raw_QOFT(paramValInd,:,termAvgInd));
        results_QOIT(paramValInd,termAvgInd)=mean(raw_QOIT(paramValInd,:,termAvgInd));
        results_QOFP(paramValInd,termAvgInd)=mean(raw_QOFP(paramValInd,:,termAvgInd));
        results_QOIP(paramValInd,termAvgInd)=mean(raw_QOIP(paramValInd,:,termAvgInd));
        results_QDFT(paramValInd,termAvgInd)=mean(raw_QDFT(paramValInd,:,termAvgInd));
        results_QDFP(paramValInd,termAvgInd)=mean(raw_QDFP(paramValInd,:,termAvgInd));
        results_QDIT(paramValInd,termAvgInd)=mean(raw_QDIT(paramValInd,:,termAvgInd));
        results_QDIP(paramValInd,termAvgInd)=mean(raw_QDIP(paramValInd,:,termAvgInd));
    end

    results_SOFS(paramValInd)=mean(raw_SOFS(paramValInd,:));
    results_SOIS(paramValInd)=mean(raw_SOIS(paramValInd,:));
    results_SDFS(paramValInd)=mean(raw_SDFS(paramValInd,:));
    results_SDIS(paramValInd)=mean(raw_SDIS(paramValInd,:));
    
    if exist('WTs','var') % for sweeping over the product coupling*tau; off by default
        WTs(WTind) = coupling * tau;       
        WTF0s(WTind) = results_MOFT(paramValInd,termAvgInd);
        WTHsizes(WTind,:)=results_MHHT(paramValInd,:);
        WTtaus(WTind)=tau;
        WTcoups(WTind)=coupling;
        
        WTind = WTind + 1;
    end
    
    % progress tracker for my impatient self
    strcat(testVarName,'_',sequenceName,'_',string(paramValInd))
end

%% Save Result Output
if saveRuns
    its = num2str(tau*10^7);
    tauString = strcat(its(1),',',its(2));
    defColors
    if strcmp(testVarName,'coupling')||strcmp(testVarName,'Coupling')||strcmp(testVarName,'coupling_lo')
        fileDescriptor = strcat(date,'_',sequenceName,'_',testVarName,tauString,']',string(Delta),'_REC_',string(maxTerm));
    elseif strcmp(testVarName,'tau')||strcmp(testVarName,'Tau')||strcmp(testVarName,'tau_lo')
        fileDescriptor = strcat(date,'_',sequenceName,'_',testVarName,string(coupling),']',string(Delta),'_REC_',string(maxTerm));
    end  

    if ~exist('fileNameAdd','var')
        fileNameAdd='';
    end
    fileDescriptor = strcat(fileDescriptor,fileNameAdd);
    
        
    if do_Magnus
        fileDescriptor = strcat(fileDescriptor,"[MAG]");
    end
    if do_Floquet
        fileDescriptor = strcat(fileDescriptor,"[FLQ]");
    end
    
    save(strcat(fileDescriptor,".mat"))
end