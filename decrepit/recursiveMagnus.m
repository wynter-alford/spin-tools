%% recursiveMagnus.m
% Wynter Alford
% January 2021
%
% Calculates terms in the Magnus expansion using the recursive form from
% Klarsfeld & Oteo, 1989

%% Setup
global dim pulse f1 Pulses X Y Z getNextUsM knownOmegas knownPs knownQs %#ok<*GVMIS>


% Valid sequences are: WHH, MREV8, CORY48, YXX48, YXX24, YXX24S, AZ48
sequenceName = 'WHH' ;  % select sequence to test over.
testVarName = 'tau'; % select parameter to test over (Delta, Tau, Coupling, coupling_lo)


% Control Parameters
testValueCount = 20;
couplingsCount = 4; % how many different coupling matrices to average over
maxTerm = 12; % how many orders of AHT to compute

N = 4;
pulse = 1.4e-6;
tau = 4.4e-6;  % delay spacing
coupling = 1500;
Delta = 0;
overRotation = 0; %0 is a perfect pulse, +/- 0.01 is a 1% over/under rotation, etc.


%Derived parameters
dim = 2^N;
f1 = 1/4/pulse; %Can adjust f1 and w1 by changing 'pulse' variable
w1 = 2*pi*f1;
rotationError = 1+overRotation; %1 is a perfect pulse

%Configure Test Variable
if strcmp(testVarName,'tau')||strcmp(testVarName,'Tau')
    %testValueMax = 10e-6;
    testValueMax = 3e-6;
elseif strcmp(testVarName,'Delta')||strcmp(testVarName,'delta')
    testValueMax = 1000;
elseif strcmp(testVarName,'Coupling')||strcmp(testVarName,'coupling')
    testValueMax = 5000;
elseif strcmp(testVarName,'coupling_lo')
    testValueMax = 5000;
elseif strcmp(testVarName,'Overrotations')||strcmp(testVarName,'overrotations')
    testValueMax = 0.04;
end

% initvars, initCollectiveObs;
z=0.5*sparse([1 0; 0 -1]);
x=0.5*sparse( [ 0 1;1 0]); 
y=1i*0.5*sparse([0 -1;1 0]);

Z=sparse(dim,dim);
X=sparse(dim,dim);
Y=sparse(dim,dim);
for k=1:N
      Z = Z + mykron(speye(2^(k-1)),z,speye(2^(N-k)));
      X = X + mykron(speye(2^(k-1)),x,speye(2^(N-k)));
      Y = Y + mykron(speye(2^(k-1)),y,speye(2^(N-k)));
end

getNextUsM = memoize(@getNextUs); % used later but needs to be defined here

%% Initialize Hamiltonians (modified code)

Hdips = cell(couplingsCount,1);
% Generate [couplingsCount] dipole Hamiltonians, each with a different coupling matrix
for j=1:couplingsCount
    dip = abs(randn(N));
    dip = triu(dip,1) + triu(dip,1)';
    Hdips{j} = getHdip(N, dim, x, y, z, dip);
end

%% Commutator Method

testVars = zeros(testValueCount,1);
results = zeros(testValueCount,maxTerm);

for d=1:length(testVars)
    
    % select variable to test over based on input
    
    if strcmp(testVarName,'tau')||strcmp(testVarName,'Tau')
        tau = d*(testValueMax/testValueCount); % ADJUST FOR TEST VAR 
        testVars(d)=tau; % ADJUST FOR TEST VAR
    elseif strcmp(testVarName,'Delta')||strcmp(testVarName,'delta')
        Delta = 2*(d-(testValueCount/2))*(testValueMax/testValueCount);
        testVars(d)=Delta;
    elseif strcmp(testVarName,'Coupling')||strcmp(testVarName,'coupling')||strcmp(testVarName,'coupling_lo')
        coupling = d*(testValueMax/testValueCount);
        testVars(d)=coupling;
    elseif strcmp(testVarName,'Overrotations')||strcmp(testVarName,'overrotations')
        overRotation = 2*(d-(testValueCount/2))*(testValueMax/testValueCount);
        rotationError = 1+overRotation;
        testVars(d) = overRotation;
    end
    
        
    for c=1:couplingsCount
        
        % Get sequence to test
        
        sequence = getSequence(sequenceName);        
        Pulses = sequence.Pulses;
        for editPulse = 1:length(Pulses)
            Pulses{editPulse}=Pulses{editPulse}*rotationError;
        end
        Taus = tau * sequence.Taus;
        tCyc = sum(Taus);

        Hdip = Hdips{c};
        Hsys = Hdip*coupling + Z*Delta;
        
        hsys = matOrder(Hsys);
        
        toggledHsys = {};
        for p = 0:length(Pulses)
            toggledHsys{p+1} = getURF(p)'*Hsys*getURF(p); %#ok<SAGROW> 
        end
        knownOmegas = {};
        knownPs = {};
        knownQs = {};
        for hterm=1:maxTerm
            results(d,hterm)=results(d,hterm)+matOrder((1i/tCyc)*Omega(hterm,toggledHsys,Taus));
        end
        
    end
end

results = results/4;


%% Recursive Method
global knownOmegas knownPs knownQs %#ok<*GVMIS> 
knownOmegas = {};
knownPs = {};
knownQs = {};

rH0 = (1i/tCyc)*Omega(1,toggledHsys,Taus); %should correspond to H0
rH1 = (1i/tCyc)*Omega(2,toggledHsys,Taus); %should correspond to H1
rH2 = (1i/tCyc)*Omega(3,toggledHsys,Taus);
rH3 = (1i/tCyc)*Omega(4,toggledHsys,Taus);

Omega1=Omega(1,toggledHsys,Taus);
Omega2=Omega(2,toggledHsys,Taus);
Omega3=Omega(3,toggledHsys,Taus);
P1=P(1,toggledHsys,Taus);
P2=P(2,toggledHsys,Taus);
P3=P(3,toggledHsys,Taus);

tH1=toggledHsys{1};
tH2=toggledHsys{2};

H0 %#ok<NOPTS>
rH0 %#ok<NOPTS>
H1 %#ok<NOPTS>
rH1 %#ok<NOPTS>
H2 %#ok<NOPTS>
rH2 %#ok<NOPTS>
