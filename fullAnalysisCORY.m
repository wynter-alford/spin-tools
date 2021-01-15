%% Indicate test sequence
testSequence = 'CORY';


%% Initialize Test Variables

% Potential Test Vars
Nin = 4;
dim = 2^Nin ;
cutoffTime = 10e-3;
%deltas = [-10000 -7500 -5000 -2500 -1000 -500 -100 0 100 500 1000 2500 5000 7500 10000];
%couplings = [1000 5000 10000 20000 40000 60000];

% Setup 4-parameter grid
meshing=30; %number of values of each param. to test over

maxTau = 12;
maxPulse = 0.5;
maxDelta = 650;
maxCoupling = 6500;

taus=zeros(meshing,1); 
pulses=zeros(1,1);
deltas = zeros(meshing,1);
couplings=zeros(meshing,1);

for i=1:meshing
    taus(i)=(maxTau/meshing)*i*10^(-6);
    %pulses(i)=(maxPulse/meshing)*i*10^(-6);
    %deltas(i)=(maxDelta/meshing)*i;
    deltas(i)=(-1) * maxDelta + 2*(maxDelta/meshing)*i; % FIX TO AVOID NEGATIVE DELTA
    couplings(i)=(maxCoupling/meshing)*i;
end
pulses(1)=0.7;

%% Generate Coupling Matrices

as = generateCouplingF(Nin);


%% Result Storage Matrices
resultsP0=zeros(meshing+1,meshing+1, 15); % update with number of 'a' matrices
resultsP1=zeros(meshing+1,meshing+1);  % Temporary arrays

resultArray=zeros(meshing+1,meshing+1, length(deltas), length(couplings));  % Final results stored here


%% Setup Code Normally Run by spinSimNoPlots

% initialize Pauli matrices and transformation operators
z=0.5*sparse([1 0; 0 -1]);x=0.5*sparse( [ 0 1;1 0]); y=1i*0.5*sparse([0 -1;1 0]);
ep=sparse([1 0; 0 0]);
em=sparse([0 0; 0 1]);
id=speye(2);p=sparse([0 1;0 0]);m=sparse([0 0; 1 0]);


% collective spin observables

Z=sparse(dim,dim);
X=sparse(dim,dim);
Y=sparse(dim,dim);

for k=1:Nin
      Z = Z + mykron(speye(2^(k-1)),z,speye(2^(Nin-k)));
      X = X + mykron(speye(2^(k-1)),x,speye(2^(Nin-k)));
      Y = Y + mykron(speye(2^(k-1)),y,speye(2^(Nin-k)));
end



%% Run Simulations
for couplingCount=1:length(couplings)
    couplingIn = couplings(couplingCount);
    for deltaCount=1:length(deltas)
        resultsP0=zeros(meshing+1,meshing+1);
        deltaIn = deltas(deltaCount);
        for tauCount=1:meshing
            tauIn=taus(tauCount);
            
            for pulseCount=1:meshing
                pulseIn=(pulses(pulseCount));
                for j=1:length(as)
                    dipIn = as{j};

        % RUN THE SPIN SIMULATION IN HERE - - - - - - - - - - - - - - -            
                    
                    
% define system parameters
                    N = Nin;

                    %symb = strcat('b','m','g','k','r','c');
                    % TODO what does this do?

                    pulse = pulseIn; %1 is pi/2, 0.25e-6 is default
                    tau = tauIn;   % minimum delay between pulses (3e-6 default)
                    Tcyc = 6*tau + 4*pulse; % cycle time

                    coupling = couplingIn; %5 kHz coupling  DON'T CHANGE / eh ok to change
                    f1 = 1/4/pulse; %Can adjust f1 and w1 by changing 'pulse' variable
                    w1 = 2*pi*f1;

                    Delta = deltaIn;    %500 Hz Chemical shift (default)
                    
                    dip = dipIn;

                    % define system observables and Hamiltonians

                    Hdip = getHdip(N, dim, x, y, z, dip);

                    % define random polarization axis and collective observable
                    theta = pi*rand(); phi = 2*pi*rand();
                    NHAT = getNHAT(N, dim, theta, phi);

                    %physical Hamiltonian
                    Hint = Hdip*coupling + Z*Delta; 

                    %Lowest order WAHUHA and MREV hamiltonian
                    tc = 6*tau + 4*pulse;
                    a = 3 * (pulse/tc)*(4/pi - 1);

                    tcM = 12*tau + 8*pulse;
                    aM = 3 * (pulse/tcM)*(4/pi - 1);

                    % Define Unitaries

                    %experimental unitary operators
                    Utau=(expm(-1i*Hint*2*pi*tau));
                    UhalfTau = (expm(-1i*Hint*pi*tau));
                    Ux = (expm(-1i*2*pi*(Hint+f1*X)*pulse));
                    Uy = (expm(-1i*2*pi*(Hint+f1*Y)*pulse));
                    Uxbar = (expm(-1i*2*pi*(Hint-f1*X)*pulse));
                    Uybar = (expm(-1i*2*pi*(Hint-f1*Y)*pulse));

                    %delta pulse unitary operators
                    U1 = (expm(-1i*2*pi*(f1*X)*pulse));
                    U2 = (expm(-1i*2*pi*(f1*Y)*pulse));
                    U3 = (expm(-1i*2*pi*(-f1*X)*pulse));
                    U4 = (expm(-1i*2*pi*(-f1*Y)*pulse));
                    
                    % SIMULATE SPIN DYNAMICS - - - - - - - - - - - -

                    % (Choose desired sequence to test)
                    if strcmp(testSequence, 'WHH')
                        testUnitary = Utau*Uxbar*Utau*Uy*Utau*Utau*Uybar*Utau*Ux*Utau;
                        Havg = (Delta / 3) * (X + Y + Z) * (1 + a); 
                        testCyc = 12;
                        cycleTime = 4 * pulse + 6 * tau;
                        U0 = expm(-1i*Havg*2*pi*cycleTime);
                        
                    elseif strcmp(testSequence, 'MREV8')
                        testUnitary = Utau*Uxbar*Utau*Uy*Utau*Utau*Uybar*Utau*Ux*Utau*Utau*Ux*Utau*Uy*Utau*Utau*Uybar*Utau*Uxbar*Utau;
                        Havg = (Delta / 3) * (X + Z) * (1 + 2*aM);
                        testCyc = 6;
                        cycleTime = 8 * pulse + 12 * tau;
                        U0 = expm(-1i*Havg*2*pi*cycleTime);
                        
                    %elseif strcmp(testSequence, 'MREV16')
                    % TODO: ADD THIS    
                        
                    elseif strcmp(testSequence, 'WK1')
                        testUnitary = UhalfTau*Uxbar*Uxbar*Uy*Utau*Uxbar*Utau*Uy*Uxbar*UhalfTau;
                        Havg = (Delta / 3) * (X + Y + Z) * (1 + a); % average Hamiltonian for WHH-4
                        testCyc = 1; % FIX THIS
                        cycleTime = 3 * tau + 6 * pulse;
                        U0 = expm(-1i*Havg*2*pi*cycleTime);
                        
                    elseif strcmp(testSequence, 'CORY') % Not currently functional
                        testUnitary = Utau*Ux*Utau*Uybar*Utau*Utau*Ux*Utau*Uy*Utau*Utau*Ux*Utau*Uybar*Utau*Utau*Uxbar*Utau*Uy*Utau*Utau*Ux*Utau*Uy*Utau*Utau*Uxbar*Utau*Uy*Utau*Utau*Uybar*Utau*Ux*Utau*Utau*Uybar*Utau*Uxbar*Utau*Utau*Uybar*Utau*Ux*Utau*Utau*Uy*Utau*Uxbar*Utau*Utau*Uybar*Utau*Uxbar*Utau*Utau*Uy*Utau*Uxbar*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Ux*Utau*Uybar*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Uxbar*Utau*Uy*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Uy*Utau*Ux*Utau*Utau*Uybar*Utau*Ux*Utau*Utau*Uy*Utau*Ux*Utau*Utau*Uy*Utau*Ux*Utau*Utau*Uy*Utau*Uxbar*Utau*Utau*Uy*Utau*Ux*Utau;
                        Havg = 0;
                        testCyc = 1;
                        U0 = 1;
                        cycleTime = 72 * tau + 48 * pulse;
                    
                    elseif strcmp(testSequence, 'CORY24')
                        testUnitary = Utau*Ux*Utau*Uybar*Utau*Utau*Ux*Utau*Uy*Utau*Utau*Ux*Utau*Uybar*Utau*Utau*Uy*Utau*Uxbar*Utau*Utau*Uybar*Utau*Uxbar*Utau*Utau*Uy*Utau*Uxbar*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Ux*Utau*Uybar*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Uy*Utau*Ux*Utau*Utau*Uy*Utau*Uxbar*Utau*Utau*Uy*Utau*Ux*Utau;
                    end
                    
                    %Ncyc = testCyc + 1;
                    
           % SPIN SIMULATION - - - - - - - - -

                    % define initial states
                    % % X as initial state
                     rho0=X; % finite-width WHH
                    % NHAT as initial state
                    %rho0=NHAT; % finite-width WHH
                    rho1=rho0; % ideal WHH
                    rho2=rho0; % zeroth order WHH approximation

                    rhoinit = rho0;
                    normD = trace(rho0*rho0);

                    Ncyc = ceil(cutoffTime / cycleTime);
                    sig0 = zeros(Ncyc,1);
                    sig2 = zeros(Ncyc,1);
                    fidelity0 = zeros(Ncyc, 1);

                    %Wahuha Sequence
                    Ucum0 = testUnitary;
                    Ucum2 = U0;

                    simTime = 0;
                    cycleCount = 0;
                    
                    while simTime < cutoffTime
                        
                        cycleCount = cycleCount + 1;
     
                        % evolution of collective magnetization
                        sig0(cycleCount) = trace(rho0*rhoinit)/normD; %experimental
                        sig2(cycleCount) = trace(rho2*rhoinit)/normD; %first order


                        % fidelity metric for unitary operators
                        fidelity0(cycleCount) = metric(Ucum0, Ucum2, N); %exp compared to FO
                        Ucum0 = Ucum0 * testUnitary;
                        Ucum2 = Ucum2 * U0;

                        % Unitary evolution
                        rho0 = testUnitary*rho0*testUnitary';    %experimental
                        rho2 = U0*rho2*U0';                      %average hamiltonian

                        % Adjust time
                        simTime = simTime + cycleTime;
                    end
 
        % SPIN SIMULATION ENDS HERE - - - - - - - - - - - - - - - - - - 
                    
                    resultsP0(pulseCount, tauCount, j) = fidelity0(cycleCount);
                end
            end

        end
        
        %compile average data for this pulse/tau pair
        for i = 1:length(as)
            for pulseCount = 1:meshing
                for tauCount = 1:meshing
                    resultArray(pulseCount,tauCount,deltaCount,couplingCount) = resultArray(pulseCount, tauCount, deltaCount, couplingCount) + resultsP0(pulseCount, tauCount, j);
                end
            end
        end
        
        % close loops
    end
end
resultArray = resultArray / length(as);

%% Save Result Matrix

filename = strcat('test_results_', testSequence,'_fixedTime_waves?','.mat');
save(filename, 'resultArray') 


%% FUNCTION DEFINITIONS

function as = generateCouplingF(Nin)

    a1 = abs(randn(Nin));
    a1 = triu(a1,1) + triu(a1,1)';
    a2 = abs(randn(Nin));
    a2 = triu(a2,1) + triu(a2,1)';
    a3 = abs(randn(Nin));
    a3 = triu(a3,1) + triu(a3,1)';
    a4 = abs(randn(Nin));
    a4 = triu(a4,1) + triu(a4,1)';
    a5 = abs(randn(Nin));
    a5 = triu(a5,1) + triu(a5,1)';
    a6 = abs(randn(Nin));
    a6 = triu(a6,1) + triu(a6,1)';
    a7 = abs(randn(Nin));
    a7 = triu(a7,1) + triu(a7,1)';
    a8 = abs(randn(Nin));
    a8 = triu(a8,1) + triu(a8,1)';
    
    as = {a1,a2,a3,a4,a5,a6,a7, a8};
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

function NHAT = getNHAT(N, dim, theta, phi)

    nhat = [cos(theta/2)^2, 0.5*sin(theta)*exp(-1i*phi); ...
            0.5*sin(theta)*exp(1i*phi), sin(theta/2)^2];

    NHAT=sparse(dim,dim);

    for k=1:N
        NHAT = NHAT + mykron(speye(2^(k-1)),nhat,speye(2^(N-k)));
    end

end

function fidelity = metric(Utarget, Uexp, N)

fidelity = abs(trace(Utarget' * Uexp)/2^N);

end

function PO=mykron(varargin)
    %function PO=mykron(varargin)
    % This little function allows you to pass all of
    % the product operators in at once and it goes through
    % and makes the right answer.  It is equivalent to
    % kron(kron(kron(x,x),x),x)....
    % written April 29th 2000, -Evan Fortunato

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
