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
global dim Pulses knownOmegas knownQs knownPs 

% Derived parameters
dim = 2^N;
f1 = 1/(4*pulse); 
w1 = 2*pi*f1;
rotationError = 1+overRotation; %1 is a perfect pulse

%Configure Test Variable
testValueMax = getVarMax(testVarName,varMaxMod);

%% Initialize Quantum Setting

initVars
initCollectiveObs
initHamiltonians

%% Construct Result Storage Arrays
% Array naming convention is results_ABC with a 3-letter code.
% 
% The first letter is the METRIC in question.  The metric options are:
% O: overlap(Uexp, AHTUnitaries{termInd}) is the propagator overlap of the
% experimental and AHT propagators.
% D: uDist(Uexp, AHTUnitaries{termInd}) is the unitary distance (currently
% using the spectral norm) between the experimental and AHT propagators.
% H: specnorm(MagnusTerms{termInd}) is the norm of the average Hamiltonian
% term.
%
% The second letter is the PULSE TREATMENT for constructing Uexp, the
% EXPERIMENTAL propagator. The options are:
% I: Instantaneous pulses (deltaUnitary)
% F: Finite pulses (expUnitary)
% H: used if the Metric is H, since in that case treatment of Uexp is
% irrelevant.
%
% The third letter is the PULSE METHOD FOR AHT, indicating whether the
% Average Hamiltonian computations used pulse divisions. The options are:
% T: Traditional analysis assuming instantaneous pulses
% P: Pulse length incorporated into the AHT computations.
% S: Uexp compared only to the identity matrix; AHT terms not computed
% (the true Time-Suspension fidelity) (Never occurs with HH).
%
% For example, results_OFT is the overlap between testUnitary (Uexp for
% finite pulses) and AHTUnitaries{termInd} (U_AHT,n not accounting for
% pulse widths in the AHT computations).

testVars = zeros(testValueCount,1);

% Size of Magnus Terms
results_HHT = zeros(length(testVars),maxTerm+1);
results_HHP = zeros(length(testVars),maxTerm+1);

raw_HHT = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_HHP = zeros(length(testVars),couplingsCount,maxTerm+1);

% Hn Overlaps (fn)
results_OIT = zeros(length(testVars),maxTerm+1);
results_OFT = zeros(length(testVars),maxTerm+1);
results_OFP = zeros(length(testVars),maxTerm+1);
results_OFS = zeros(length(testVars),1);
results_OIS = zeros(length(testVars),1);

raw_OIT = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_OFT = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_OFP = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_OFS = zeros(length(testVars),couplingsCount);
raw_OIS = zeros(length(testVars),couplingsCount);

% Hn Propagator Distances (dn)
results_DIT = zeros(length(testVars),maxTerm+1);
results_DFT = zeros(length(testVars),maxTerm+1);
results_DFP = zeros(length(testVars),maxTerm+1);
results_DFS = zeros(length(testVars),1);
results_DIS = zeros(length(testVars),1);

raw_DFS = zeros(length(testVars),couplingsCount);
raw_DIS = zeros(length(testVars),couplingsCount);
raw_DFT = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_DIT = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_DFP = zeros(length(testVars),couplingsCount,maxTerm+1);

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
            toggle = getURF(p,Pulses,X,Y,UDx,UDy,UDxbar,UDybar);
            toggledHsys{p+1} = toggle'*Hsys*toggle; 
        end

        % Set up finite-pulse (HHP) toggled Hsys
        if pulseDivs > 1
            PulsesP = repmat({ones(1)},1,length(Pulses)*pulseDivs);
            for semiPulseInd=1:length(PulsesP)
                PulsesP{semiPulseInd}=PulsesP{semiPulseInd}*Pulses{ceil(semiPulseInd/pulseDivs)};
            end
            for p = 0:length(Pulses)*pulseDivs
                toggle = getURF(p,PulsesP,X,Y,UdivX,UdivY,UdivXbar,UdivYbar);
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
                
        % Compute Magnus Series terms
        getMagnusTerms;

        % obtain theoretical unitaries from AHT and compute Fidelity
        AHTUnitaries = {};
        AHTUnitariesP = {};
        sumAH = sparse(dim,dim);
        sumAHP = sparse(dim,dim);
        
        for termInd=1:maxTerm+1
            sumAH = sumAH + MagnusTerms{termInd};
            sumAHP = sumAHP + MagnusTermsP{termInd};
            AHTUnitaries{termInd} = expm(-1i*sumAH*2*pi*tCyc); %#ok<SAGROW>
            AHTUnitariesP{termInd} = expm(-1i*sumAHP*2*pi*tCyc); %#ok<SAGROW>  

            
            raw_OFT(paramValInd,coupMatInd,termInd)=overlap(expUnitary, AHTUnitaries{termInd}, N);
            raw_OIT(paramValInd,coupMatInd,termInd) = overlap(deltaUnitary, AHTUnitaries{termInd}, N);
            raw_OFP(paramValInd,coupMatInd,termInd) = overlap(expUnitary, AHTUnitariesP{termInd}, N);

            raw_DFT(paramValInd,coupMatInd,termInd)=uDist(expUnitary,AHTUnitaries{termInd},N);
            raw_DIT(paramValInd,coupMatInd,termInd)=uDist(deltaUnitary,AHTUnitaries{termInd},N);
            raw_DFP(paramValInd,coupMatInd,termInd) = overlap(expUnitary, AHTUnitariesP{termInd}, N);
        end
        
        % Time-suspension fidelities
        raw_OFS(paramValInd,coupMatInd) = overlap(expUnitary,speye(dim,dim),N);
        raw_OIS(paramValInd,coupMatInd) = overlap(deltaUnitary,speye(dim,dim),N); 
    
        raw_DFS(paramValInd,coupMatInd) = uDist(expUnitary,speye(dim,dim),N);
        raw_DIS(paramValInd,coupMatInd) = uDist(deltaUnitary,speye(dim,dim),N); 
        
    end
    
    % Average Raw Results   
    for termAvgInd = 1:maxTerm+1
        results_HHT(paramValInd,termAvgInd)=mean(raw_HHT(paramValInd,:,termAvgInd));
        results_HHP(paramValInd,termAvgInd)=mean(raw_HHP(paramValInd,:,termAvgInd));
        results_OFT(paramValInd,termAvgInd)=mean(raw_OFT(paramValInd,:,termAvgInd));
        results_OIT(paramValInd,termAvgInd)=mean(raw_OIT(paramValInd,:,termAvgInd)); 
        results_OFP(paramValInd,termAvgInd)=mean(raw_OFP(paramValInd,:,termAvgInd)); 
        results_DFT(paramValInd,termAvgInd)=mean(raw_DFT(paramValInd,:,termAvgInd));
        results_DIT(paramValInd,termAvgInd)=mean(raw_DIT(paramValInd,:,termAvgInd));
        results_DFP(paramValInd,termAvgInd)=mean(raw_DFP(paramValInd,:,termAvgInd));
    end

    results_OFS(paramValInd)=mean(raw_OFS(paramValInd,:));
    results_OIS(paramValInd)=mean(raw_OIS(paramValInd,:));
    results_DFS(paramValInd)=mean(raw_DFS(paramValInd,:));
    results_DIS(paramValInd)=mean(raw_DIS(paramValInd,:));
    
    if exist('WTs','var') % for sweeping over the product coupling*tau; off by default
        WTs(WTind) = coupling * tau;       
        WTF0s(WTind) = results_OFT(paramValInd,termAvgInd);
        WTHsizes(WTind,:)=results_HHT(paramValInd,:);
        WTtaus(WTind)=tau;
        WTcoups(WTind)=coupling;
        
        WTind = WTind + 1;
    end
    
    % progress tracker for my impatient self
    %strcat(testVarName,'_',sequenceName,'_',string(paramValInd))
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
    fileDescriptor = strcat(fileDescriptor,fileNameAdd,'.mat');
    save(fileDescriptor)
end