%% sequenceAnalysis.m
% Wynter Alford
% January 2022
%
% For a given pulse sequence, does the following:
%    - Calculates the first [maxTerm] terms of the magnus
%    expansion for that sequence using the recursive method
%
%    - Calculates the fidelity of that sequence against H=0, H=H0,
%    and H=H0+H1+... up to the highest order term chosen
%
%    - Calculates fidelity when using experimental or instantaneous pulses
%
%    - Can analyze sequence behaviour and robustness with respect to 
%    offsets, tau spacing, pulse length, dipolar coupling strength,
%    overrotations, and phase transients

%% Control Parameters (change these)

global dim pulse f1 Pulses knownOmegas knownQs knownPs

% Valid sequences are: WHH, MREV8, CORY48, YXX48, YXX24, YXX24S, AZ48,
%    AS10M, AS10E, WHHWHH, ML10, WPW9-2Cycle, MG8

sequenceName = 'MREV8' ;  % select sequence to test over.
testVarName = 'coupling_lo'; % select parameter to test over (Delta, Tau, Coupling, coupling_lo, phaseTrans, delta_lo, 'overrot_hi')
varMaxMod = 1; % factor by which to multiply the default testVarMax

mode = 'max'; %can be 'max' to get up to maxTerm terms, or 'time' to compute for a certain amount of time
maxTerm = 10; % highest Magnus series term to compute 
computationTime = 

% Control Parameters
testValueCount = 1;
couplingsCount = 4; % how many different coupling matrices to average over

N = 4;
pulse = 1.4e-6;
tau = 3.4e-6;  % delay spacing
coupling = 5000;
Delta = 00;
overRotation = 0; %0 is a perfect pulse, +/- 0.01 is a 1% over/under rotation, etc.
phaseTrans = 0; %0 is no phase transient, 1 is a pi/2 phase transient

%% Configuration

% Derived parameters
dim = 2^N;
f1 = 1/4/pulse; %Can adjust f1 and w1 by changing 'pulse' variable
w1 = 2*pi*f1;
rotationError = 1+overRotation; %1 is a perfect pulse

%Configure Test Variable
testValueMax = getVarMax(testVarName,varMaxMod);

initVars
initCollectiveObs
% % initvars, initCollectiveObs;
% z=0.5*sparse([1 0; 0 -1]);
% x=0.5*sparse( [ 0 1;1 0]); 
% y=1i*0.5*sparse([0 -1;1 0]);
% 
% Z=sparse(dim,dim);
% X=sparse(dim,dim);
% Y=sparse(dim,dim);
% for k=1:N
%       Z = Z + mykron(speye(2^(k-1)),z,speye(2^(N-k)));
%       X = X + mykron(speye(2^(k-1)),x,speye(2^(N-k)));
%       Y = Y + mykron(speye(2^(k-1)),y,speye(2^(N-k)));
% end

% %% Initialize Hamiltonians (modified code)
% 
%  Hdips = cell(couplingsCount,1);
%  % Generate [couplingsCount] dipole Hamiltonians, each with a different coupling matrix
%  for j=1:couplingsCount
%      dip = abs(randn(N));
%      dip = triu(dip,1) + triu(dip,1)';
%      Hdips{j} = getHdip(N, dim, x, y, z, dip);
%  end
load('Standard_Dipole_Hamiltonians(4).mat','Hdips');

%% Iterate over different parameter values to see how term magnitude changes

testVars = zeros(testValueCount,1);

% Size of Magnus Terms
results_hsizes = zeros(length(testVars),maxTerm+1);
raw_hsizes = zeros(length(testVars),couplingsCount,maxTerm+1);

% Hn Fidelities
results_f = zeros(length(testVars),maxTerm+1);
results_Df = zeros(length(testVars),maxTerm+1);
raw_f = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_Df = zeros(length(testVars),couplingsCount,maxTerm+1);

% Time-Suspension Fidelities
results_fTS = zeros(length(testVars),1);
results_DfTS = zeros(length(testVars),1);
raw_fTS = zeros(length(testVars),couplingsCount);
raw_DfTS = zeros(length(testVars),couplingsCount);

% Commutation with Other Terms
results_C0 = zeros(length(testVars),maxTerm+1);
results_CS = zeros(length(testVars),maxTerm+1);
raw_C = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_CS = zeros(length(testVars),couplingsCount,maxTerm+1);

for d=1:length(testVars)
    
    % Configure Test Variable and Test Sequence
    if strcmp(testVarName,'tau')||strcmp(testVarName,'Tau')
        tau = pulse + d*(testValueMax/testValueCount);
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
        hsys = matOrder(Hsys);
        
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
            raw_f(d,c,au)=metric(expUnitary, AHTUnitaries{au}, N);
            raw_Df(d,c,au) = metric(deltaUnitary, AHTUnitaries{au}, N);
        end
        
        % Time-suspension fidelities
        raw_fTS(d,c) = metric(expUnitary,speye(dim,dim),N);
        raw_DfTS(d,c) = metric(deltaUnitary,speye(dim,dim),N); 

    end
    
    % Average Raw Results   
    for mt = 1:maxTerm+1
        results_hsizes(d,mt)=mean(raw_hsizes(d,:,mt));
        results_f(d,mt)=mean(raw_f(d,:,mt));
        results_Df(d,mt)=mean(raw_Df(d,:,mt)); 
        results_C0(d,mt)=mean(raw_C0(d,:,mt));
        results_CS(d,mt)=mean(raw_CS(d,:,mt));
    end

    results_fTS(d)=mean(raw_fTS(d,:));
    results_DfTS(d)=mean(raw_DfTS(d,:));

    % progress tracker for my impatient self
    strcat(testVarName,'_',sequenceName,'_',string(d))
end


%% Save Result Output
its = num2str(tau*10^7);
tauString = strcat(its(1),',',its(2));
if strcmp(testVarName,'coupling')||strcmp(testVarName,'Coupling')||strcmp(testVarName,'coupling_lo')
    fileDescriptor = strcat(date,'_',sequenceName,'_',testVarName,tauString,']',string(Delta),'_REC_',string(maxTerm),'.mat');
elseif strcmp(testVarName,'tau')||strcmp(testVarName,'Tau')
    fileDescriptor = strcat(date,'_',sequenceName,'_',testVarName,string(coupling),']',string(Delta),'_REC_',string(maxTerm),'.mat');
end
save(fileDescriptor, 'sequenceName','results_hsizes', 'testVars','tau','coupling','Delta','N','couplingsCount','results_f','results_Df','pulse','results_fTS','results_DfTS','Hdips','maxTerm','tCyc','results_C0','results_CS','coupling','Delta')
