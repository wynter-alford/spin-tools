%% Indicate test sequence
testSequence = 'WHH';
% CURRENTLY USES TWO TERMS IN MAGNUS FOR WHH AND MREV8

%% Initialize Test Variables

% Potential Test Vars
N = 4;
global dim Pulses
dim = 2^N ;
Pulses = {};
cutoffTime = 1e-3;
%deltas = [-10000 -7500 -5000 -2500 -1000 -500 -100 0 100 500 1000 2500 5000 7500 10000];
%couplings = [1000 5000 10000 20000 40000 60000];

% Setup 4-parameter grid
meshing=30; %number of values of each param. to test over

maxTau = 14;
maxPulse = 1.4; % low is 0.5; normal is 1.4
maxDelta = 2000; % low is 300; normal is 10000
maxCoupling = 4000; % low is 8000; normal is 60000

taus=zeros(meshing,1); 
pulses=zeros(meshing,1);
deltas = zeros(meshing,1);
couplings=zeros(meshing,1);

for i=1:meshing
    taus(i)=(maxTau/meshing)*i*10^(-6);
    pulses(i)=(maxPulse/meshing)*i*10^(-6);
    deltas(i)=(maxDelta/meshing)*i;  % POSITIVE DELTA ONLY
    %deltas(i)=(-1) * maxDelta + 2*(maxDelta/meshing)*i; % NEGATIVE DELTA
    couplings(i)=(maxCoupling/meshing)*i;
end

%% Generate Coupling Matrices

as = generateCouplingF(N);


%% Setup Code Normally Run by spinSimNoPlots

% initialize Pauli matrices and transformation operators
z=0.5*sparse([1 0; 0 -1]);x=0.5*sparse( [ 0 1;1 0]); y=1i*0.5*sparse([0 -1;1 0]);
ep=sparse([1 0; 0 0]);
em=sparse([0 0; 0 1]);
id=speye(2);p=sparse([0 1;0 0]);m=sparse([0 0; 1 0]);


% collective spin observables
global X Y 

Z=sparse(dim,dim);
X=sparse(dim,dim);
Y=sparse(dim,dim);

for k=1:N
      Z = Z + mykron(speye(2^(k-1)),z,speye(2^(N-k)));
      X = X + mykron(speye(2^(k-1)),x,speye(2^(N-k)));
      Y = Y + mykron(speye(2^(k-1)),y,speye(2^(N-k)));
end

%% Initialize Pulse Sequence

sequence = getSequence(testSequence);
Pulses = sequence.Pulses;


%% Setup Result Storage

raw_f_results3 = zeros(length(as), meshing, meshing, meshing, meshing);
raw_f_results1 = zeros(length(as), meshing, meshing, meshing, meshing);
raw_h1s = zeros(length(as),meshing,meshing,meshing);
raw_h2s = zeros(length(as),meshing,meshing,meshing);


%% Run Simulations

for j=1:length(as)
    dip = as{j};
    Hdip = getHdip(N, dim, x, y, z, dip); %dipole Hamiltonian

    for couplingCount=1:meshing
        coupling = couplings(couplingCount);
        
        for deltaCount=1:meshing
            Delta = deltas(deltaCount);
            Hint = Hdip*coupling + Z*Delta; % system Hamiltonian
            
            % Compute all possible toggled Hamiltonians
            toggledHsys = {};
            for p = 0:length(Pulses)
                toggledHsys{p+1} = getURF(p)'*Hint*getURF(p);
            end
            
            for tauCount=1:meshing
                tau=taus(tauCount);  % lowercase 'taus' are the different values to test             
                Taus = tau * sequence.Taus;   % uppercase 'Taus' are the times between pulses in the sequence (usually either tau or 2*tau)
                Utau=(expm(-1i*Hint*2*pi*tau));
                UhalfTau = (expm(-1i*Hint*pi*tau));

% ------------- CALCULATE AVERAGE HAMILTONIAN (UP TO Hbar(2)) -------------

                if strcmp(testSequence, 'WHH')
                    H0 = (Delta / 3) * (X + Y + Z); 
                    H2 = zeros(dim,dim);
                    
                    for lll=0:length(sequence.Pulses)
                        for qqq=0:lll
                            for jjj=0:qqq
                                Hl = toggledHsys{lll+1};
                                Hk = toggledHsys{qqq+1};
                                Hj = toggledHsys{jjj+1};

                                Hterm = comm(Hl,comm(Hk,Hj))+comm(comm(Hl,Hk),Hj);
                                H2 = H2 + Hterm*Taus(lll+1)*Taus(qqq+1)*Taus(jjj+1);
                            end
                        end
                    end

                    H2 = (-1/(6*sum(Taus)))*H2;    
                    Havg = H0 + H2;
                    h2 = matOrder(H2);
                    raw_h2s(j, tauCount, deltaCount, couplingCount) = h2;

                elseif strcmp(testSequence, 'MREV8')
                    H0 = (Delta / 3) * (X + Z); % * (1 + 2*aM);
                    
                    H1 = zeros(dim,dim);
                    for lll=1:length(sequence.Pulses)
                        for jjj=0:lll-1
                            Hk = toggledHsys{lll+1};
                            Hj = toggledHsys{jjj+1};
                            H1 = H1 + comm(Hk,Hj)*Taus(lll+1)*Taus(jjj+1);
                        end
                    end

                    H1 = (1/(2*1i*sum(Taus)))*H1;
                    h1 = matOrder(H1);
                    raw_h1s(j, tauCount, deltaCount, couplingCount) = h1;

                    %Calculate 2nd order Magnus term
                    H2 = zeros(dim,dim);
                    for lll=0:length(Pulses)
                        for qqq=0:lll
                            for jjj=0:qqq
                                Hl = toggledHsys{lll+1};
                                Hk = toggledHsys{qqq+1};
                                Hj = toggledHsys{jjj+1};

                                Hterm = comm(Hl,comm(Hk,Hj))+comm(comm(Hl,Hk),Hj);
                                H2 = H2 + Hterm*Taus(lll+1)*Taus(qqq+1)*Taus(jjj+1);
                            end
                        end
                    end

                    H2 = (-1/(6*sum(Taus)))*H2;  
                    h2 = matOrder(H2);
                    raw_h2s(j, tauCount, deltaCount, couplingCount) = h2;

                    Havg = H0+H1+H2;  
                    %Havg = H0;
                end
            
% ------------- END AVERAGE HAMILTONIAN CALCULATIONS ---------------------

                for pulseCount=1:meshing
                    pulse=(pulses(pulseCount));

                    f1 = 1/4/pulse; %Can adjust f1 and w1 by changing 'pulse' variable
                    w1 = 2*pi*f1;

                    % Define Unitaries

                    %experimental unitary operators
                    Ux = (expm(-1i*2*pi*(Hint+f1*X)*pulse));
                    Uy = (expm(-1i*2*pi*(Hint+f1*Y)*pulse));
                    Uxbar = (expm(-1i*2*pi*(Hint-f1*X)*pulse));
                    Uybar = (expm(-1i*2*pi*(Hint-f1*Y)*pulse));

                    %delta pulse unitary operators
                    %Ux = (expm(-1i*2*pi*(f1*X)*pulse));
                    %Uy = (expm(-1i*2*pi*(f1*Y)*pulse));
                    %Uxbar = (expm(-1i*2*pi*(-f1*X)*pulse));
                    %Uybar = (expm(-1i*2*pi*(-f1*Y)*pulse));
               
                

       % - - - - - - - - - - SIMULATE SPIN DYNAMICS - - - - - - - - - - - -

                    % (Choose desired sequence to test)
                    if strcmp(testSequence, 'WHH')
                        testUnitary = Utau*Ux*Utau*Uybar*Utau*Utau*Uy*Utau*Uxbar*Utau;
                        cycleTime = 4 * pulse + 6 * tau;
                        U3 = expm(-1i*Havg*2*pi*cycleTime);
                        U1 = expm(-1i*H0*2*pi*cycleTime);
                        
                    elseif strcmp(testSequence, 'MREV8')
                        testUnitary = Utau*Uxbar*Utau*Uy*Utau*Utau*Uybar*Utau*Ux*Utau*Utau*Ux*Utau*Uy*Utau*Utau*Uybar*Utau*Uxbar*Utau;
                        cycleTime = 8 * pulse + 12 * tau;            
                        U3 = expm(-1i*Havg*2*pi*cycleTime); 
                        U1 = expm(-1i*H0*2*pi*cycleTime);

                        
                    elseif strcmp(testSequence, 'WK1')
                        testUnitary = UhalfTau*Uxbar*Uxbar*Uy*Utau*Uxbar*Utau*Uy*Uxbar*UhalfTau;
                        Havg = (Delta / 3) * (X + Y + Z) * (1 + a); % average Hamiltonian for WHH-4
                        cycleTime = 3 * tau + 6 * pulse;
                        U1 = expm(-1i*Havg*2*pi*cycleTime);
                        
                    elseif strcmp(testSequence, 'CORY48')
                        testUnitary = Utau*Ux*Utau*Uybar*Utau*Utau*Ux*Utau*Uy*Utau*Utau*Ux*Utau*Uybar*Utau*Utau*Uxbar*Utau*Uy*Utau*Utau*Ux*Utau*Uy*Utau*Utau*Uxbar*Utau*Uy*Utau*Utau*Uybar*Utau*Ux*Utau*Utau*Uybar*Utau*Uxbar*Utau*Utau*Uybar*Utau*Ux*Utau*Utau*Uy*Utau*Uxbar*Utau*Utau*Uybar*Utau*Uxbar*Utau*Utau*Uy*Utau*Uxbar*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Ux*Utau*Uybar*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Uxbar*Utau*Uy*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Uy*Utau*Ux*Utau*Utau*Uybar*Utau*Ux*Utau*Utau*Uy*Utau*Ux*Utau*Utau*Uy*Utau*Ux*Utau*Utau*Uy*Utau*Uxbar*Utau*Utau*Uy*Utau*Ux*Utau;
                        Havg = 0;
                        U1 = 1;
                        cycleTime = 72 * tau + 48 * pulse;
                    
                    elseif strcmp(testSequence, 'CORY24')
                        testUnitary = Utau*Ux*Utau*Uybar*Utau*Utau*Ux*Utau*Uy*Utau*Utau*Ux*Utau*Uybar*Utau*Utau*Uy*Utau*Uxbar*Utau*Utau*Uybar*Utau*Uxbar*Utau*Utau*Uy*Utau*Uxbar*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Ux*Utau*Uybar*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Uy*Utau*Ux*Utau*Utau*Uy*Utau*Uxbar*Utau*Utau*Uy*Utau*Ux*Utau;
                        Havg = 0;
                        U1 = 1;
                        cycleTime = 36 * tau + 24* pulse;
                        
                    %elseif strcmp(testSequence, 'YXX24')
                        % Do this
                    %elseif strcmp(testSequence, 'YXX48')
                        % Do this
                    end
                    
                    %Ncyc = testCyc + 1;
                    
           % SPIN SIMULATION - - - - - - - - -

                    % define initial states
                    % % X as initial state
            %        rho0=X; % finite-width WHH
                    % NHAT as initial state
                    %rho0=NHAT; % finite-width WHH
            %        rho1=rho0; % ideal WHH
            %        rho2=rho0; % zeroth order WHH approximation

            %        rhoinit = rho0;
            %        normD = trace(rho0*rho0);

                    %Ncyc = ceil(cutoffTime / cycleTime);
                    Ncyc = 1;
                    testCyc = 1;
                    fidelity3 = zeros(Ncyc, 1);
                    fidelity1 = zeros(Ncyc, 1);

                    %Wahuha Sequence
                    Ucum0 = testUnitary; % experimental
                    Ucum1 = U1; % 1 term
                    Ucum3 = U3; % 3 terms (currently only WHH, MREV8)

                    %simTime = 0;
                    cycleCount = 0;
                    
                    %while simTime < cutoffTime
                    while cycleCount < testCyc    
                    
                        cycleCount = cycleCount + 1;
                        
                        % fidelity metric for unitary operators
                        fidelity3(cycleCount) = metric(Ucum0, Ucum3, N); %exp compared to 3terms
                        fidelity1(cycleCount) = metric(Ucum0, Ucum1, N); %exp compared to 1terms
                        Ucum0 = Ucum0 * testUnitary;
                        Ucum1 = Ucum1 * U1; % 1 term
                        Ucum3 = Ucum3 * U3; % 3 terms

                        % Unitary evolution
                %        rho0 = testUnitary*rho0*testUnitary';    %experimental
                %        rho1 = U1*rho1*U1';                      %average hamiltonian (1 term)
                %        rho2 = U3*rho2*U3';                      %average hamiltonian (3 terms)

                        % Adjust time
                        %simTime = simTime + cycleTime;
                    end
                    raw_f_results3(j, pulseCount, tauCount, deltaCount, couplingCount) = fidelity3(cycleCount);
                    raw_f_results1(j, pulseCount, tauCount, deltaCount, couplingCount) = fidelity1(cycleCount);
        % SPIN SIMULATION ENDS HERE - - - - - - - - - - - - - - - - - - 
                    % COMPILE RESULTS
                end
            end
        end

    end
    %save(strcat(string(j),'_',testSequence,'.mat'),'raw_f_results3','raw_f_results1','raw_h1s','raw_h2s')
end

%% Save Results

filename = strcat(testSequence,'_H0+H1+H2_meshing=',string(meshing),'_',date,'.mat');
save(filename, 'raw_f_results3', 'raw_f_results1', 'taus','pulses','deltas','couplings','meshing', 'raw_h1s', 'raw_h2s') 

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
%    a7 = abs(randn(Nin));
%    a7 = triu(a7,1) + triu(a7,1)';
%    a8 = abs(randn(Nin));
%    a8 = triu(a8,1) + triu(a8,1)';
    
    as = {a1,a2,a3,a4,a5,a6};%,a7, a8};
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

function ct = comm(A,B) % calculates the commutator of a pair of matrices
    ct = A*B-B*A;
end

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

function orh = matOrder(A) % calculates the magnitude of a matrix
    orh = log10(sqrt(trace(A*A)));
end

function sequence = getSequence(sequenceName)
    global X Y
    %WAHUHA
    if strcmp(sequenceName, 'WHH')
        sequence.Pulses = {-X, Y, -Y, X};
        sequence.Taus = [1 1 2 1 1];

    %MREV-8
    elseif strcmp(sequenceName, 'MREV8')
        sequence.Pulses = {-X, -Y, Y, X, X, -Y, Y, -X}; % Check this
        sequence.Taus = [1 1 2 1 2 1 2 1 1];

    %CORY 48
    elseif strcmp(sequenceName, 'CORY48')
        sequence.Pulses = {X, Y, -X, Y, X, Y, X, Y, X, -Y, X, Y, -Y, -X, Y, -X, -Y, -X, -Y, -X, -Y, X, -Y, -X, -X, Y, -X, -Y, -X, Y, X, -Y, -X, -Y, X, -Y, Y, -X, Y, X, Y, -X, -Y, X, Y, X, -Y, X};
        sequence.Taus = [1 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 1];
    end
end
