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

testVars = zeros(testValueCount,1);

% Size of Magnus Terms
results_hsizes = zeros(length(testVars),maxTerm+1);
raw_hsizes = zeros(length(testVars),couplingsCount,maxTerm+1);

% Hn Overlaps (fn)
results_f = zeros(length(testVars),maxTerm+1);
results_Df = zeros(length(testVars),maxTerm+1);
raw_f = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_Df = zeros(length(testVars),couplingsCount,maxTerm+1);

% Hn Propagator Distances (dn)
results_d = zeros(length(testVars),maxTerm+1);
results_Dd = zeros(length(testVars),maxTerm+1);
raw_d = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_Dd = zeros(length(testVars),couplingsCount,maxTerm+1);

% Time-Suspension Fidelities
results_fTS = zeros(length(testVars),1);
results_DfTS = zeros(length(testVars),1);
raw_fTS = zeros(length(testVars),couplingsCount);
raw_DfTS = zeros(length(testVars),couplingsCount);

results_dTS = zeros(length(testVars),1);
results_DdTS = zeros(length(testVars),1);
raw_dTS = zeros(length(testVars),couplingsCount);
raw_DdTS = zeros(length(testVars),couplingsCount);

% Commutation with Other Terms
% results_C0 = zeros(length(testVars),maxTerm+1); % [Hn,H0]
% results_CS = zeros(length(testVars),maxTerm+1); % [Hn,H0+H1+...+H(n-1)]
% raw_C = zeros(length(testVars),couplingsCount,maxTerm+1);
% raw_CS = zeros(length(testVars),couplingsCount,maxTerm+1);

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
        
        toggledHsys = {};
        for p = 0:length(Pulses)
            toggledHsys{p+1} = getURF(p,X,Y,UDx,UDy,UDxbar,UDybar)'*Hsys*getURF(p,X,Y,UDx,UDy,UDxbar,UDybar); %#ok<*SAGROW> 
        end
                
        %Calculate Magnus terms
        knownOmegas = {};
        knownPs = {};
        knownQs = {};        
        MagnusTerms = {};
        getMagnusTerms;

        % obtain theoretical unitaries from AHT and compute Fidelity
        AHTUnitaries = {};
        sumAH = sparse(dim,dim);
        
        for termInd=1:maxTerm+1
            sumAH = sumAH + MagnusTerms{termInd};
            AHTUnitaries{termInd} = expm(-1i*sumAH*2*pi*tCyc);
            
            raw_f(paramValInd,coupMatInd,termInd)=overlap(expUnitary, AHTUnitaries{termInd}, N);
            raw_Df(paramValInd,coupMatInd,termInd) = overlap(deltaUnitary, AHTUnitaries{termInd}, N);
            
            raw_d(paramValInd,coupMatInd,termInd)=uDist(expUnitary,AHTUnitaries{termInd},N);
            raw_Dd(paramValInd,coupMatInd,termInd)=uDist(expUnitary,AHTUnitaries{termInd},N);

        end
        
        % Time-suspension fidelities
        raw_fTS(paramValInd,coupMatInd) = overlap(expUnitary,speye(dim,dim),N);
        raw_DfTS(paramValInd,coupMatInd) = overlap(deltaUnitary,speye(dim,dim),N); 
    
        raw_dTS(paramValInd,coupMatInd) = uDist(expUnitary,speye(dim,dim),N);
        raw_DdTS(paramValInd,coupMatInd) = uDist(deltaUnitary,speye(dim,dim),N); 
        
    end
    
    % Average Raw Results   
    for termAvgInd = 1:maxTerm+1
        results_hsizes(paramValInd,termAvgInd)=mean(raw_hsizes(paramValInd,:,termAvgInd));
        results_f(paramValInd,termAvgInd)=mean(raw_f(paramValInd,:,termAvgInd));
        results_Df(paramValInd,termAvgInd)=mean(raw_Df(paramValInd,:,termAvgInd)); 
        results_d(paramValInd,termAvgInd)=mean(raw_d(paramValInd,:,termAvgInd));
        results_Dd(paramValInd,termAvgInd)=mean(raw_Dd(paramValInd,:,termAvgInd));
    end

    results_fTS(paramValInd)=mean(raw_fTS(paramValInd,:));
    results_DfTS(paramValInd)=mean(raw_DfTS(paramValInd,:));
    results_dTS(paramValInd)=mean(raw_dTS(paramValInd,:));
    results_DdTS(paramValInd)=mean(raw_DdTS(paramValInd,:));
    
    if exist('WTs','var') % for sweeping over the product coupling*tau; off by default
        WTs(WTind) = coupling * tau;       
        WTF0s(WTind) = results_f(paramValInd,termAvgInd);
        WTHsizes(WTind,:)=results_hsizes(paramValInd,:);
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