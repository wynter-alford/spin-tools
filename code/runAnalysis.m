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
f1 = 1/4/pulse; 
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
results_C0 = zeros(length(testVars),maxTerm+1);
results_CS = zeros(length(testVars),maxTerm+1);
raw_C = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_CS = zeros(length(testVars),couplingsCount,maxTerm+1);

%% Iterate over different parameter values, compute AHT terms and fidelities

for d=1:length(testVars)
    
    % Configure Test Variable and Test Sequence
    if strcmp(testVarName,'tau')||strcmp(testVarName,'Tau')
        tau = pulse + d*(testValueMax/testValueCount);
        %tau = 10^(-d);
        testVars(d)=tau;
    elseif strcmp(testVarName,'tau_lo')
        tau = d*(testValueMax/testValueCount);
        testVars(d)=tau;
    elseif strcmp(testVarName,'Delta')||strcmp(testVarName,'delta')||strcmp(testVarName,'delta_lo')
        Delta = 2*(d-(testValueCount/2))*(testValueMax/testValueCount);
        testVars(d)=Delta;
    elseif strcmp(testVarName,'Coupling')||strcmp(testVarName,'coupling')||strcmp(testVarName,'coupling_lo')
        coupling = d*(testValueMax/testValueCount);
        testVars(d)=coupling;
    elseif strcmp(testVarName,'Overrotations')||strcmp(testVarName,'overrotations')||strcmp(testVarName,'overrot_hi')
        overRotation = 2*(d-(testValueCount/2))*(testValueMax/testValueCount);
        rotationError = 1+overRotation;
        testVars(d) = overRotation;
    elseif strcmp(testVarName,'phaseTrans')
        phaseTrans = 2*(d-(testValueCount/2))*(testValueMax/testValueCount);
        testVars(d) = phaseTrans;
    end
  
    sequence = getSequence(sequenceName,X,Y);
    

    % Add overrotations and phase transients
    for editPulse = 1:length(sequence.Pulses)
        sequence.Pulses{editPulse}=pulseError(sequence.Pulses{editPulse},rotationError,phaseTrans,X,Y);
    end
    
    Pulses = sequence.Pulses;
    Taus = tau * sequence.Taus;        
    tCyc = sum(Taus);
    
    % For each dipolar coupling matrix
    for c=1:couplingsCount  
        
        Hdip = Hdips{c};
        Hsys = Hdip*coupling + Z*Delta;        
        hsys = specnorm(Hsys);
        
        toggledHsys = {};
        for p = 0:length(Pulses)
            toggledHsys{p+1} = getURF(p)'*Hsys*getURF(p); %#ok<*SAGROW> 
        end
        
        %Calculate Magnus terms
        knownOmegas = {};
        knownPs = {};
        knownQs = {};        
        MagnusTerms = {};
        getMagnusTerms;
        
        % obtain experimental unitary
        testUnitaries = getTestUnitaries(sequence,Hsys,pulse,tau,f1);
        expUnitary = testUnitaries{1};
        deltaUnitary = testUnitaries{2};

        % obtain theoretical unitaries from AHT and compute Fidelity
        AHTUnitaries = {};
        sumAH = sparse(dim,dim);
        
        for au=1:maxTerm+1
            sumAH = sumAH + MagnusTerms{au};
            AHTUnitaries{au} = expm(-1i*sumAH*2*pi*tCyc);
            
            raw_f(d,c,au)=overlap(expUnitary, AHTUnitaries{au}, N);
            raw_Df(d,c,au) = overlap(deltaUnitary, AHTUnitaries{au}, N);
            
            raw_d(d,c,au)=uDist(expUnitary,AHTUnitaries{au},N);
            raw_Dd(d,c,au)=uDist(expUnitary,AHTUnitaries{au},N);

        end
        
        % Time-suspension fidelities
        raw_fTS(d,c) = overlap(expUnitary,speye(dim,dim),N);
        raw_DfTS(d,c) = overlap(deltaUnitary,speye(dim,dim),N); 
    
        raw_dTS(d,c) = uDist(expUnitary,speye(dim,dim),N);
        raw_DdTS(d,c) = uDist(deltaUnitary,speye(dim,dim),N); 
        
    end
    
    % Average Raw Results   
    for mt = 1:maxTerm+1
        results_hsizes(d,mt)=mean(raw_hsizes(d,:,mt));
        results_f(d,mt)=mean(raw_f(d,:,mt));
        results_Df(d,mt)=mean(raw_Df(d,:,mt)); 
        results_d(d,mt)=mean(raw_d(d,:,mt));
        results_Dd(d,mt)=mean(raw_Dd(d,:,mt));
    end

    results_fTS(d)=mean(raw_fTS(d,:));
    results_DfTS(d)=mean(raw_DfTS(d,:));
    results_dTS(d)=mean(raw_dTS(d,:));
    results_DdTS(d)=mean(raw_DdTS(d,:));
    
    % progress tracker for my impatient self
    % strcat(testVarName,'_',sequenceName,'_',string(d))
end

%% Make Plotting More Convenient by saving colors now

myColors = {[250 190 212],[230 25 75],[245 130 48],  [210 245 60], [60 180 75], [70 240 240], [0 130 200], [145 30 180], [0 0 0]};
% colors are:        pink          red      orange           lime          green        cyan         blue         violet       black  

for i=1:length(myColors)
    myColors{i} = myColors{i}/255;    
end


%% Save Result Output
if saveRuns
    its = num2str(tau*10^7);
    tauString = strcat(its(1),',',its(2));
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