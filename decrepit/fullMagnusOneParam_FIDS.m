%% Setup
global dim pulse f1 Pulses X Y Z getNextUsM

% Valid sequences are: WHH, MREV8, CORY48, YXX48, YXX24, YXX24S, AZ48,
%    AS10M, AS10E, WHHWHH, ML10, WPW9-2Cycle

%sequenceName = 'WHH' ;  % select sequence to test over.
%testVarName = 'Tau'; % select parameter to test over (Delta, Tau, Coupling, coupling_lo)


% Control Parameters
%testValueCount = 40;
couplingsCount = 6; % how many different coupling matrices to average over

N = 4;
pulse = 1.4e-6;
tau = 4.4e-6;  % delay spacing
coupling = 2000;
Delta = 0;
overRotation = 0; %0 is a perfect pulse, +/- 0.01 is a 1% over/under rotation, etc.


%Derived parameters
dim = 2^N;
f1 = 1/4/pulse; %Can adjust f1 and w1 by changing 'pulse' variable
w1 = 2*pi*f1;
rotationError = 1+overRotation; %1 is a perfect pulse

%Configure Test Variable
if strcmp(testVarName,'tau')||strcmp(testVarName,'Tau')
    testValueMax = 10e-6;
elseif strcmp(testVarName,'Delta')||strcmp(testVarName,'delta')
    testValueMax = 5000;
elseif strcmp(testVarName,'Coupling')||strcmp(testVarName,'coupling')
    testValueMax = 50000;
elseif strcmp(testVarName,'coupling_lo')
    testValueMax = 5000;
elseif strcmp(testVarName,'Overrotations')||strcmp(testVarName,'overrotations')
    testValueMax = 0.03;
elseif strcmp(testVarName,'tau1')
    testValueMax = 2e-6;
    testVarName = 'tau';
elseif strcmp(testVarName,'tau2')
    testValueMax = 8e-6;
    testVarName = 'tau';
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

%  Hdips = cell(couplingsCount,1);
%  % Generate [couplingsCount] dipole Hamiltonians, each with a different coupling matrix
%  for j=1:couplingsCount
%      dip = abs(randn(N));
%      dip = triu(dip,1) + triu(dip,1)';
%      Hdips{j} = getHdip(N, dim, x, y, z, dip);
%  end

%% Iterate over different parameter values to see how term magnitude changes

testVars = zeros(testValueCount,1);

results_h0 = zeros(length(testVars),1);
results_h1 = zeros(length(testVars),1);
results_h2 = zeros(length(testVars),1);
results_h3 = zeros(length(testVars),1);
results_h4 = zeros(length(testVars),1);

results_f0 = zeros(length(testVars),1);
results_f2 = zeros(length(testVars),1);
results_f4 = zeros(length(testVars),1);

results_Df0 = zeros(length(testVars),1);
results_Df2 = zeros(length(testVars),1);
results_Df4 = zeros(length(testVars),1);

raw_results_h0 = zeros(length(testVars),couplingsCount);
raw_results_h1 = zeros(length(testVars),couplingsCount);
raw_results_h2 = zeros(length(testVars),couplingsCount);
raw_results_h3 = zeros(length(testVars),couplingsCount);
raw_results_h4 = zeros(length(testVars),couplingsCount);

raw_f0 = zeros(length(testVars),couplingsCount);
raw_f2 = zeros(length(testVars),couplingsCount);
raw_f4 = zeros(length(testVars),couplingsCount);

raw_Df0 = zeros(length(testVars),couplingsCount);
raw_Df2 = zeros(length(testVars),couplingsCount);
raw_Df4 = zeros(length(testVars),couplingsCount);

raw_fTS = zeros(length(testVars),couplingsCount);
raw_DfTS = zeros(length(testVars),couplingsCount);

results_fTS = zeros(length(testVars),1);
results_DfTS = zeros(length(testVars),1);

raw_f2Only = zeros(length(testVars),couplingsCount);
f2Only= zeros(length(testVars),1);

raw_h02_closeness = zeros(length(testVars),couplingsCount);
h02_closeness= zeros(length(testVars),1);


for d=1:length(testVars)
    
    % Configure Test Variable and Test Sequence
    
    if strcmp(testVarName,'tau')||strcmp(testVarName,'Tau')
        tau = pulse + d*(testValueMax/testValueCount); % ADJUST FOR TEST VAR 
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
  
    sequence = getSequence(sequenceName);
        
    for editPulse = 1:length(sequence.Pulses)
        sequence.Pulses{editPulse}=sequence.Pulses{editPulse}*rotationError;
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
            toggledHsys{p+1} = getURF(p)'*Hsys*getURF(p);
        end

        %calcMagnus
        %Magnus Calculation is done here ------------------------------------
        
        %Calculate 0th order Magnus term
        H0 = zeros(dim,dim);
        for k=0:length(Pulses)
            H0 = H0 + (Taus(k+1))*toggledHsys{k+1};
        end

        H0 = (1/tCyc)*H0;
        h0 = matOrder(H0)/hsys;


        %Calculate 1st order Magnus term
        H1 = zeros(dim,dim);
        for k=1:length(Pulses)
            for j=0:k-1
                Hk = toggledHsys{k+1};
                Hj = toggledHsys{j+1};
                H1 = H1 + comm(Hk,Hj)*Taus(k+1)*Taus(j+1);
            end
        end

        H1 = (1/(2*1i*tCyc))*H1;
        h1 = matOrder(H1)/hsys;

        %Calculate 2nd order Magnus term
        H2 = zeros(dim,dim);
        for l=0:length(Pulses)
            for k=0:l
                for j=0:k
                    Hl = toggledHsys{l+1};
                    Hk = toggledHsys{k+1};
                    Hj = toggledHsys{j+1};

                    Hterm = comm(Hl,comm(Hk,Hj))+comm(comm(Hl,Hk),Hj);
                    H2 = H2 + Hterm*Taus(l+1)*Taus(k+1)*Taus(j+1);
                end
            end
        end

        H2 = (-1/(6*tCyc))*H2;
        h2 = matOrder(H2)/hsys;

        % Calculate 3rd order Magnus term
%         H3 = zeros(dim,dim);
%         for m=0:length(Pulses)
%             for l=0:m
%                 for k=0:l
%                     for j=0:k
%                         Hm = toggledHsys{m+1};
%                         Hl = toggledHsys{l+1};
%                         Hk = toggledHsys{k+1};
%                         Hj = toggledHsys{j+1};
% 
%                         term1 = comm(comm(comm(Hm,Hl),Hk),Hj);
%                         term2 = comm(Hm,comm(comm(Hl,Hk),Hj));
%                         term3 = comm(Hm,comm(Hl,comm(Hk,Hj)));
%                         term4 = comm(Hl,comm(Hk,comm(Hj,Hm)));
% 
%                         Hterm = term1+term2+term3+term4;
%                         tauProd = Taus(m+1)*Taus(l+1)*Taus(k+1)*Taus(j+1);
% 
%                         H3 = H3 + Hterm*tauProd;
%                     end
%                 end
%             end
%         end
% 
%         H3 = (-1/(12*1i*tCyc))*H3;
%         h3 = matOrder(H3)/hsys;

%        Calculate the 4th term of the Magnus Expansion
%         H4 = zeros(dim,dim);
%         for mm=0:length(Pulses)
%             for m=0:mm
%                 for l=0:m
%                     for k=0:l
%                         for j=0:k
%                             % Express Hsys(t) in the interaction frame
%                             Hmm = toggledHsys{mm+1};
%                             Hm = toggledHsys{m+1};
%                             Hl = toggledHsys{l+1};
%                             Hk = toggledHsys{k+1};
%                             Hj = toggledHsys{j+1};
% 
%                             % Calculate commutators
%                             term1 = (-1/30)*comm(Hmm,comm(Hm,comm(Hl,comm(Hk,Hj))));
%                             term2 = (2/15)*comm(Hj,comm(Hmm,comm(Hm,comm(Hk,Hl))));
%                             term3 = (1/15)*comm(comm(Hmm,Hj),comm(Hm,comm(Hk,Hl)));
%                             term4 = (1/15)*comm(comm(Hm,Hj),comm(Hmm,comm(Hk,Hl)));
%                             term5 = (-1/60)*comm(comm(Hk,Hl),comm(Hmm,comm(Hm,Hj)));
%                             term6 = (1/15)*comm(comm(Hl,Hj),comm(Hmm,comm(Hk,Hm)));
%                             term7 = (-1/60)*comm(comm(Hk,Hm),comm(Hmm,comm(Hl,Hj)));
%                             term8 = (-1/60)*comm(comm(Hk,Hmm),comm(Hm,comm(Hl,Hj)));
%                             term9 = (-1/60)*comm(comm(Hl,Hm),comm(Hmm,comm(Hk,Hj)));
%                             term10 = (-1/60)*comm(comm(Hl,Hm),comm(Hj,comm(Hk,Hmm)));
%                             term11 = (-1/60)*comm(comm(Hmm,Hj),comm(Hl,comm(Hk,Hm)));
%                             term12 = (-1/60)*comm(comm(Hm,Hj),comm(Hl,comm(Hk,Hmm)));
%                             term13 = (-1/60)*comm(comm(Hl,Hmm),comm(Hm,comm(Hk,Hj)));
%                             term14 = (-1/60)*comm(comm(Hl,Hmm),comm(Hj,comm(Hk,Hm)));
%                             term15 = (-1/30)*comm(Hj,comm(Hm,comm(Hl,comm(Hk,Hmm))));
%                             term16 = (-1/60)*comm(comm(Hm,Hmm),comm(Hj,comm(Hk,Hl)));
%                             term17 = (-1/60)*comm(comm(Hk,Hl),comm(Hj,comm(Hm,Hmm)));
%                             term18 = (-1/60)*comm(comm(Hk,Hm),comm(Hj,comm(Hl,Hmm)));
%                             term19 = (-1/60)*comm(comm(Hk,Hj),comm(Hm,comm(Hl,Hmm)));
%                             term20 = (-1/60)*comm(comm(Hm,Hmm),comm(Hl,comm(Hk,Hj)));
%                             term21 = (-1/60)*comm(comm(Hl,Hj),comm(Hm,comm(Hk,Hmm)));
%                             term22 = (-1/30)*comm(Hj,comm(Hmm,comm(Hl,comm(Hk,Hm))));
% 
%                             % Add commutators
%                             Hterm = term1+term2+term3+term4+term5+term6+term7+term8+term9+term10+term11+term12+term13+term14+term15+term16+term17+term18+term19+term20+term21+term22;
%                             tauProd = Taus(mm+1)*Taus(m+1)*Taus(l+1)*Taus(k+1)*Taus(j+1);
% 
%                             H4 = H4 + Hterm*tauProd;
%                         end
%                     end
%                 end
%             end
%         end
% 
%         H4 = (1/tCyc)*H4;
%         h4 = matOrder(H4)/hsys;


        % Magnus Calculation ends here ---------------------------------------


        raw_results_h0(d,c) = h0;
        raw_results_h1(d,c) = h1;
        raw_results_h2(d,c) = h2;
%        raw_results_h3(d,c) = h3;
%        raw_results_h4(d,c) = h4;
        
       
        % FIDELITY CALCULATION HERE
        % get experimental unitary
        testUnitaries = getTestUnitaries(sequence,Hsys,pulse,tau,f1);
        testUnitary = testUnitaries{1};
        deltaUnitary = testUnitaries{2};

        U0 = expm(-1i*H0*2*pi*tCyc);
        U2 = expm(-1i*(H0+H1+H2)*2*pi*tCyc);
        U2Only = expm(-1i*H2*2*pi*tCyc);
%        U4Only = expm(-1i*H4*2*pi*tCyc);
%        U4 = expm(-1i*(H0+H1+H2+H3+H4)*2*pi*tCyc);

        raw_f0(d,c) = metric(testUnitary, U0, N);
        raw_f2(d,c) = metric(testUnitary, U2, N);
%        raw_f4(d,c) = metric(testUnitary, U4, N);

        raw_Df0(d,c) = metric(deltaUnitary, U0, N);
        raw_Df2(d,c) = metric(deltaUnitary, U2, N);
%        raw_Df4(d,c) = metric(deltaUnitary, U4, N); % "fully corrected" fidelity

        raw_fTS(d,c) = metric(testUnitary,speye(dim,dim),N);
        raw_DfTS(d,c) = metric(deltaUnitary,speye(dim,dim),N);

        raw_f2Only(d,c) = metric(testUnitary,U2Only,N);
        %f4Only(d) = metric(testUnitary,U4Only,N);
        
        raw_h02_closeness(d,c) = metric(U0,U2Only,N);
    end
    
    results_h0(d)=mean(raw_results_h0(d,:));
    results_h1(d)=mean(raw_results_h1(d,:));
    results_h2(d)=mean(raw_results_h2(d,:));
%    results_h3(d)=mean(raw_results_h3(d,:));
%    results_h4(d)=mean(raw_results_h4(d,:));
    
    results_f0(d)=mean(raw_f0(d,:));
    results_f2(d)=mean(raw_f2(d,:));
%    results_f4(d)=mean(raw_f4(d,:));
    
    results_Df0(d)=mean(raw_Df0(d,:));
    results_Df2(d)=mean(raw_Df0(d,:));
%    results_Df4(d)=mean(raw_Df0(d,:));

    %time suspension fidelities
    results_fTS(d)=mean(raw_fTS(d,:));
    results_DfTS(d)=mean(raw_DfTS(d,:));
    
    f2Only(d)=mean(raw_f2Only(d,:));
    h02_closeness(d)=mean(raw_h02_closeness(d,:));
    testVarName,sequenceName,d % progress tracker
end


%% Save Result Output
fileDescriptor = strcat(sequenceName,'_',testVarName,'_02D-fids_',date,'.mat');
save(fileDescriptor, 'sequenceName','results_h0', 'results_h1', 'results_h2', 'results_h3', 'results_h4', 'testVars','tau','coupling','Delta','N','couplingsCount','results_f0','results_f2','results_f4','results_Df0','results_Df2','results_Df4','pulse','results_fTS','results_DfTS','Hdips')

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

% SEQUENCES (new code)
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

function fidelity = metric(Utarget, Uexp, N)

fidelity = abs(trace(Utarget' * Uexp)/2^N);

end

function norm = specnorm(A) % Spectral Norm (sqrt maximum eigenvalue of A'A)
    mat = A'*A;
    evals = eig(mat);
    norm = sqrt(max(evals));
end

    
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

function nextUs = getNextUs(sequence,Hsys,pulse,f1,p)
    nextU = expm(-1i*2*pi*(Hsys+f1*sequence.Pulses{p})*pulse);
    nextUD = expm(-1i*Hsys*2*pi*pulse/2)*expm(-1i*pi*sequence.Pulses{p}/2)*expm(-1i*Hsys*2*pi*pulse/2);  
    nextUs = {nextU,nextUD};
end

