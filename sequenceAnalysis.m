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

%% Setup
global dim pulse f1 Pulses X Y Z getNextUsM

% Valid sequences are: WHH, MREV8, CORY48, YXX48, YXX24, YXX24S, AZ48,
%    AS10M, AS10E, WHHWHH, ML10, WPW9-2Cycle

sequenceName = 'MREV8' ;  % select sequence to test over.
testVarName = 'coupling_lo'; % select parameter to test over (Delta, Tau, Coupling, coupling_lo, phaseTrans, delta_lo, 'overrot_hi')
maxTerm = 16; % highest Magnus series term to compute (can do 8+ for WHH; no higher than 4 for 48 pulses)

% Control Parameters
testValueCount = 1;
couplingsCount = 1; % how many different coupling matrices to average over

N = 4;
pulse = 1.4e-6;
tau = 5.4e-6;  % delay spacing
coupling = 5000;
Delta = 500;
overRotation = 0; %0 is a perfect pulse, +/- 0.01 is a 1% over/under rotation, etc.
phaseTrans = 0; %0 is no phase transient, 1 is a pi/2 phase transient

% Derived parameters
dim = 2^N;
f1 = 1/4/pulse; %Can adjust f1 and w1 by changing 'pulse' variable
w1 = 2*pi*f1;
rotationError = 1+overRotation; %1 is a perfect pulse

%Configure Test Variable
if strcmp(testVarName,'tau')||strcmp(testVarName,'Tau')
    testValueMax = 10e-6;
elseif strcmp(testVarName,'Delta')||strcmp(testVarName,'delta')
    testValueMax = 1000;
elseif strcmp(testVarName,'Coupling')||strcmp(testVarName,'coupling')
    testValueMax = 50000;
elseif strcmp(testVarName,'coupling_lo')
    testValueMax = 5000;
elseif strcmp(testVarName,'Overrotations')||strcmp(testVarName,'overrotations')
    testValueMax = 0.03;
elseif strcmp(testVarName,'phaseTrans')
    testValueMax = 0.15;
elseif strcmp(testVarName,'delta_lo')
    testValueMax = 400;
elseif strcmp(testVarName,'overrot_hi')
    testValueMax = 0.2;
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


% %% Initialize Hamiltonians (modified code)
% 
%  Hdips = cell(couplingsCount,1);
%  % Generate [couplingsCount] dipole Hamiltonians, each with a different coupling matrix
%  for j=1:couplingsCount
%      dip = abs(randn(N));
%      dip = triu(dip,1) + triu(dip,1)';
%      Hdips{j} = getHdip(N, dim, x, y, z, dip);
%  end

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
  
    sequence = getSequence(sequenceName);
    
    % Add overrotations and phase transients
    for editPulse = 1:length(sequence.Pulses)
        sequence.Pulses{editPulse}=pulseError(sequence.Pulses{editPulse},rotationError,phaseTrans);
    end
    
    Pulses = sequence.Pulses;
    Taus = tau * sequence.Taus;        
    tCyc = sum(Taus);
    

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
        length(knownOmegas)
        knownPs = {};
        knownQs = {};
        
        MagnusTerms = {};
        trunc = zeros(dim,dim);
        for mt=1:maxTerm+1
            length(knownOmegas)
            MagnusTerms{mt} = (1i/tCyc)*Omega(mt,toggledHsys,Taus);
            trunc = trunc + MagnusTerms{mt};
            
            raw_hsizes(d,c,mt) = matOrder(MagnusTerms{mt});
            raw_C0(d,c,mt) = matOrder(comm(MagnusTerms{mt},MagnusTerms{1}));
            raw_CS(d,c,mt) = matOrder(comm(MagnusTerms{mt},trunc));
        end       
                
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

fileDescriptor = strcat(date,'_',sequenceName,'_',testVarName,string(coupling),'|',string(Delta),'_REC_',string(maxTerm),'.mat');
save(fileDescriptor, 'sequenceName','results_hsizes', 'testVars','tau','coupling','Delta','N','couplingsCount','results_f','results_Df','pulse','results_fTS','results_DfTS','Hdips','maxTerm','tCyc','results_C0','results_CS')















%% FUNCTION DEFINITIONS
function ct = comm(A,B) % calculates the commutator of a pair of matrices
    ct = A*B-B*A;
end

function order = matOrder(A) % calculates the magnitude of a matrix
    order = real(sqrt(trace(A*A)));
end

% Find the frame transformation U_rf to eliminate the RF Hamiltonian
function URF = getURF(frame)
    global dim Pulses
    
    if frame < 1
        URF = speye(dim,dim); %returns the identity if frame == 0

    else
        URF = expm(-1i*Pulses{1}*pi/2);
    end
    
    if frame > 1
        for j=2:frame
            URF = expm(-1i*Pulses{j}*pi/2) * URF;
        end
    end
end

% SEQUENCES
function sequence = getSequence(sequenceName)
    global X Y
    %WAHUHA
    if strcmp(sequenceName, 'WHH')
        sequence.Pulses = {X, -Y, Y, -X};
        sequence.Taus = [1 1 2 1 1];

    %MREV-8
    elseif strcmp(sequenceName, 'MREV8')
        sequence.Pulses = {-X, -Y, Y, X, X, -Y, Y, -X};
        sequence.Taus = [1 1 2 1 2 1 2 1 1];

    %BR-24    
    elseif strcmp(sequenceName, 'BR24')
        sequence.Pulses = {X, -Y, Y, -X, -X, -Y, Y, X, -Y, X, -X, Y, Y, X, -Y, X,    -X, Y, Y, X,   -X, -Y, -X, -Y };
        sequence.Taus = [1 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 1];
         
    %CORY 48
    elseif strcmp(sequenceName, 'CORY48')
        sequence.Pulses = {X, Y, -X, Y, X, Y, X, Y, X, -Y, X, Y, -Y, -X, Y, -X, -Y, -X, -Y, -X, -Y, X, -Y, -X, -X, Y, -X, -Y, -X, Y, X, -Y, -X, -Y, X, -Y, Y, -X, Y, X, Y, -X, -Y, X, Y, X, -Y, X};
        sequence.Taus = [1 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 1];
    
    %YXX-48
    elseif strcmp(sequenceName,'YXX48')
        sequence.Pulses = {Y, -X,-X,Y,-X,-X,-Y,X,X,Y,-X,-X,-Y,X,X,-Y,X,X,Y,-X,-X,Y,-X,-X,-Y,X,X,Y,-X,-X,-Y,X,X,-Y,X,X,Y,-X,-X,-Y,X,X,Y,-X,-X,-Y,X,X};
        sequence.Taus = ones(49,1);
        sequence.Taus(1)=0;
    
    %YXX-24
    elseif strcmp(sequenceName,'YXX24')
        sequence.Pulses = {-Y,X,-X,Y,-X,-X,Y,-X,X,-Y,X,X,Y,-X,X,-Y,X,X,-Y,X,-X,Y,-X,-X};
        sequence.Taus = ones(25,1);
        sequence.Taus(1)=0;
    
    %AZ-48
    elseif strcmp(sequenceName,'AZ48')
        sequence.Pulses = {-X,Y,Y,X,Y,Y,-Y,X,X,-Y,X,X,Y,X,X,-Y,X,X,-Y,X,-Y,X,X,-Y,-X,-X,Y,Y,-X,Y,Y,-Y,X,-Y,-Y,X,-Y,X,X,-Y,X,X,-Y,-X,-X,-Y,-X,-X};
        sequence.Taus = ones(49,1);
        sequence.Taus(1)=0;
    
    %Symmetrized 48 from YXX-24
    elseif strcmp(sequenceName,'YXX24S')
        sequence.Pulses = {-Y,X,-X,Y,-X,-X,Y,-X,X,-Y,X,X,Y,-X,X,-Y,X,X,-Y,X,-X,Y,-X,-X,X,X,-Y,X,-X,Y,-X,-X,Y,-X,X,-Y,-X,-X,Y,-X,X,-Y,X,X,-Y,X,-X,Y};
        sequence.Taus = [0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
    
    elseif strcmp(sequenceName,'AS10M')
        sequence.Pulses = {-X,Y,2*Y,-Y,-X,X,Y,2*Y,-Y,X};
        sequence.Taus = ones(11,1);
        sequence.Taus(6) = 2;
        
    elseif strcmp(sequenceName,'AS10E')
        sequence.Pulses = {X, -Y, Y, -X, 2*Y, X, Y, -Y, -X, 2*Y}; %(added pi-pulse with 0tau after to fix basis)
        sequence.Taus = ones(11,1);
        sequence.Taus(11) = 0;
        
    elseif strcmp(sequenceName,'AS10E-2cycle')
        sequence.Pulses = {X, -Y, Y, -X, 2*Y, X, Y, -Y, -X, X, -Y, Y, -X, 2*Y, X, Y, -Y, -X};
        sequence.Taus = ones(19,1);
        
    elseif strcmp(sequenceName,'WHHWHH')
        sequence.Pulses = {X, -Y, Y, -X, X, -Y, Y, -X};
        sequence.Taus = [1 1 2 1 1 1 1 2 1 1];
            
    elseif strcmp(sequenceName,'ML10')
        sequence.Pulses = {X,X,Y,X,X,-X,-X,-Y,-X,-X};
        sequence.Taus = [1 1 1 1 1 2 1 1 1 1 1];
        
    end
end

% fidelity metric
function fidelity = metric(Utarget, Uexp, N)

fidelity = abs(trace(Utarget' * Uexp)/2^N);

end

% Spectral Norm (sqrt maximum eigenvalue of A'A)
function norm = specnorm(A) 
    mat = A'*A;
    evals = eig(mat);
    norm = sqrt(max(evals));
end

 % calculate U_experiment for a given pulse sequence   
function testUnitaries = getTestUnitaries(sequence,Hsys,pulse,tau,f1)
    global dim getNextUsM
    
    UtauCell = {speye(dim,dim),expm(-1i*Hsys*2*pi*(tau-pulse)), expm(-1i*Hsys*2*pi*(2*tau-pulse))};
        
    testUnitary = speye(dim,dim);
    deltaUnitary = speye(dim,dim);
    
    for p=1:length(sequence.Pulses)
        Utau = UtauCell{sequence.Taus(p)+1};
        nextUs = getNextUsM(sequence,Hsys,pulse,f1,p);
        nextU = nextUs{1};
        nextUD = nextUs{2};
        testUnitary = nextU * Utau * testUnitary;
        deltaUnitary = nextUD *Utau * deltaUnitary;
    end
    testUnitary = UtauCell{sequence.Taus(length(sequence.Taus))+1} * testUnitary;
    deltaUnitary = UtauCell{sequence.Taus(length(sequence.Taus))+1} * deltaUnitary;
    
    testUnitaries = {testUnitary,deltaUnitary};
end

% used in getTestUnitaries
function nextUs = getNextUs(sequence,Hsys,pulse,f1,p)
    nextU = expm(-1i*2*pi*(Hsys+f1*sequence.Pulses{p})*pulse);
    nextUD = expm(-1i*Hsys*2*pi*pulse/2)*expm(-1i*pi*sequence.Pulses{p}/2)*expm(-1i*Hsys*2*pi*pulse/2);  
    nextUs = {nextU,nextUD};
end

% which phase transient to associate with which pulse 
% (note that alpha canbe negative)
function trP = pulseError(pulseIn,roterror,trans)  
    global X Y

    if pulseIn == X
        trP = roterror*pulseIn+trans*Y;
    elseif pulseIn == Y
        trP = roterror*pulseIn+trans*(-X);
    elseif pulseIn == -X
        trP = roterror*pulseIn+trans*(-Y);
    elseif pulseIn == -Y
        trP = roterror*pulseIn+trans*X;
    else
        trP = roterror*pulseIn; % transients not yet implemented for pi pulses
    end
end

% MAGNUS RECURSION FUNCTIONS
% Omega calculates the actual term by calling Q and P. P does the actual
% computation work while the others are just there for the recursive
% structure.

function omega = Omega(n,toggledHsys,Taus)
    global knownOmegas 
    % Check whether this term has already been computed
    if length(knownOmegas)>=n && ~isempty(knownOmegas{n})
        omega=knownOmegas{n};
    
    % if it hasn't been, then compute it
    else
        omega = P(n,toggledHsys,Taus);
        if n>1
            for k=2:n
                omega = omega - (1/factorial(k))*Q(n,k,toggledHsys,Taus);
            end
        end
        knownOmegas{n}=omega;
    end
end

function pn = P(n,toggledHsys,Taus)
    global dim knownPs %#ok<*GVMIS> 
    % Check whether this term has already been computed
    if length(knownPs)>=n && ~isempty(knownPs{n})
        pn = knownPs{n};
    
    % if it hasn't been, then compute it        
    else
        pn = zeros(dim,dim);
        times = ones(n,1);
        maxTime = length(toggledHsys);
        k = n;

        % This while loop enacts n summations and is equivalent to:
        % "for i1=1:N, for i2=1:i1 ... for i(n)=1:i(n-1), do stuff"
        while times(1)<=maxTime
            times'
            if k==1||times(k)<times(k-1)
                hprod = 1;
                for hnum = 1:n
                    % multiply -iH(t1)*...*-iH(tn)
                    hprod = hprod*(-1i)*toggledHsys{times(hnum)}*Taus(times(hnum));
                end
                
                pn = pn+hprod/indexfactor(times);

                times(k) = times(k)+1;

                if k<length(times)
                    times(k+1:n)=1;
                end
                k=n;
            else
                k = k-1;
            end
        end
        knownPs{n}=pn;
    end
end

function qnk = Q(n,k,toggledHsys,Taus)
    global knownQs %#ok<*GVMIS> 

    % Check whether this term has already been computed
    if length(knownQs)>=n && length(knownQs{n})>=k && ~isempty(knownQs{n}{k})
        qnk = knownQs{n}{k};
        
    % if it hasn't been, then compute it        
    else      
        if k==1
            qnk = Omega(n,toggledHsys,Taus);
            
        elseif k==n
            qnk = (Omega(1,toggledHsys,Taus))^n;
            
        else
            qnk = 0;
            for m=1:(n-k+1) 
                qnk = qnk + Omega(m,toggledHsys,Taus)*Q(n-m,k-1,toggledHsys,Taus);
            end            
        end
        
        knownQs{n}{k}=qnk;
    end
end

% The functions below are not written by me (Ben Alford) ----------------

%function PO=mykron(varargin)
% This little function allows you to pass all of
% the product operators in at once and it goes through
% and makes the right answer.  It is equivalent to
% kron(kron(kron(x,x),x),x)....
% written April 29th 2000, -Evan Fortunato
function PO=mykron(varargin)

    if length(varargin) < 1
        error('Please enter atleast 1 Product Operator!');
    elseif length(varargin) == 1
        PO = varargin{1};
    else
        PO = varargin{1};
        for j=2:length(varargin)
            PO=kron(PO, varargin{j});
        end
    end
end

function Hdip = getHdip(N, dim, x, y, z, a)

    % N: number of spins
    % x/y/z: Pauli matrices
    % a: dipolar interaction strengths (NxN symmetric matrix)

    Hdip = sparse(dim, dim);

    for k=1:N
        for h=k+1:N
            % Hdip = a_{kh} (3Z - X - Y)
            Hdip=Hdip+a(k,h)*(2*mykron(speye(2^(k-1)),z,speye(2^(h-k-1)),z,speye(2^(N-h)))-... 
                mykron(speye(2^(k-1)),x,speye(2^(h-k-1)),x,speye(2^(N-h)))-...
                mykron(speye(2^(k-1)),y,speye(2^(h-k-1)),y,speye(2^(N-h)))) ;             
        end
    end

end

% end functions not written by me ------------------------------