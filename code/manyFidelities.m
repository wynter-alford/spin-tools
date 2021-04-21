%sqs = {'CORY48','YXX48','AZ48'};
sqs = {'WHH','WHH','WHH'};
vars = {'Delta','Tau','coupling','coupling_lo'};
global Pulses dim

N = 4;
dim = 2^N;


testVarCount = 40;

% Setup Code Normally Run by spinSimNoPlots

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

%results are stored as: Sequence, Test Parameter, (avg) Value, (raw value)
results = zeros(3,4,testVarCount);
raw_results = zeros(3,4,testVarCount,8);

results_2 = zeros(3,4,testVarCount);
raw_results_2 = zeros(3,4,testVarCount,8);

for sq = 1:3
    testSequence = sqs{sq};
    
    for var = 1:4
        testVarName = vars{var};
        Ncyc = 1;

        tau = 3e-6;
        pulse = 1.5e-6;
        Delta = 0;
        coupling = 5000;

        
        if strcmp(testVarName, 'coupling_lo')
            maxCoupling = 5000;
        else
            maxCoupling = 50000;
        end

        maxTau = 10e-6;
        maxDelta = 1000;
        
        % Set up test parameter
        if strcmp(testVarName,'coupling')||strcmp(testVarName,'coupling_lo')
            couplings = zeros(testVarCount,1);
            for i=1:testVarCount
                couplings(i) = (i/testVarCount)*maxCoupling;
            end
        elseif strcmp(testVarName,'Tau')
            taus = zeros(testVarCount,1);
            for i=1:testVarCount
                taus(i) = (i/testVarCount)*maxTau;
            end
        elseif strcmp(testVarName,'Delta')
            deltas = zeros(testVarCount,1);
            for i=1:testVarCount
                deltas(i) = 2*(i-(testVarCount/2))*(maxDelta/testVarCount);

            end
        end
        
        % Generate Coupling Matrices

        as = generateCouplingF(N);

        % Simulation

        for dipVal=1:length(as)
            dip = as{dipVal};
            for testVal=1:length(deltas)
                Delta = deltas(testVal);
                Tcyc = 6*tau + 4*pulse; % cycle time

                f1 = 1/4/pulse; %Can adjust f1 and w1 by changing 'pulse' variable
                w1 = 2*pi*f1;

                % define system observables and Hamiltonians

                Hdip = getHdip(N, dim, x, y, z, dip);

                %physical Hamiltonian
                Hint = Hdip*coupling + Z*Delta; 

                %Lowest order WAHUHA and MREV hamiltonian
                tc = 6*tau + 4*pulse;

                tcM = 12*tau + 8*pulse;

                % Define Unitaries

                %experimental unitary operators
                Utau=(expm(-1i*Hint*2*pi*tau));
                UhalfTau = (expm(-1i*Hint*pi*tau));
                Ux = (expm(-1i*2*pi*(Hint+f1*X)*pulse));
                Uy = (expm(-1i*2*pi*(Hint+f1*Y)*pulse));
                Uxbar = (expm(-1i*2*pi*(Hint-f1*X)*pulse));
                Uybar = (expm(-1i*2*pi*(Hint-f1*Y)*pulse));

                % SIMULATE SPIN DYNAMICS - - - - - - - - - - - -

                % (Choose desired sequence to test)
                if strcmp(testSequence, 'CORY48') % Not currently functional
                    testUnitary = Utau*Ux*Utau*Uybar*Utau*Utau*Ux*Utau*Uy*Utau*Utau*Ux*Utau*Uybar*Utau*Utau*Uxbar*Utau*Uy*Utau*Utau*Ux*Utau*Uy*Utau*Utau*Uxbar*Utau*Uy*Utau*Utau*Uybar*Utau*Ux*Utau*Utau*Uybar*Utau*Uxbar*Utau*Utau*Uybar*Utau*Ux*Utau*Utau*Uy*Utau*Uxbar*Utau*Utau*Uybar*Utau*Uxbar*Utau*Utau*Uy*Utau*Uxbar*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Ux*Utau*Uybar*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Uxbar*Utau*Uy*Utau*Utau*Uxbar*Utau*Uybar*Utau*Utau*Uy*Utau*Ux*Utau*Utau*Uybar*Utau*Ux*Utau*Utau*Uy*Utau*Ux*Utau*Utau*Uy*Utau*Ux*Utau*Utau*Uy*Utau*Uxbar*Utau*Utau*Uy*Utau*Ux*Utau;
                    Havg = 0;
                    testCyc = 1;
                    U0 = 1;
                    cycleTime = 72 * tau + 48 * pulse;

                elseif strcmp(testSequence, 'YXX48')
                    testUnitary = Utau*Ux*Utau*Ux*Utau*Uybar*Utau*Uxbar*Utau*Uxbar*Utau*Uy*Utau*Ux*Utau*Ux*Utau*Uybar*Utau*Uxbar*Utau*Uxbar*Utau*Uy*Utau*Ux*Utau*Ux*Utau*Uybar*Utau*Ux*Utau*Ux*Utau*Uybar*Utau*Uxbar*Utau*Uxbar*Utau*Uy*Utau*Ux*Utau*Ux*Utau*Uybar*Utau*Uxbar*Utau*Uxbar*Utau*Uy*Utau*Uxbar*Utau*Uxbar*Utau*Uy*Utau*Ux*Utau*Ux*Utau*Uybar*Utau*Ux*Utau*Ux*Utau*Uybar*Utau*Uxbar*Utau*Uxbar*Utau*Uy*Utau*Ux*Utau*Ux*Utau*Uybar*Utau*Uxbar*Utau*Uxbar*Utau*Uy*Utau*Uxbar*Utau*Uxbar*Utau*Uy*Utau;
                    Havg = 0;
                    testCyc = 1;
                    U0 = 1;
                    cycleTime = 49*tau + 48*pulse;

                elseif strcmp(testSequence, 'AZ48')
                    testUnitary = Utau*Uxbar*Utau*Uxbar*Utau*Uybar*Utau*Uxbar*Utau*Uxbar*Utau*Uybar*Utau*Ux*Utau*Ux*Utau*Uybar*Utau*Ux*Utau*Ux*Utau*Uybar*Utau*Ux*Utau*Uybar*Utau*Uybar*Utau*Ux*Utau*Uybar*Utau*Uy*Utau*Uy*Utau*Uxbar*Utau*Uy*Utau*Uy*Utau*Uxbar*Utau*Uxbar*Utau*Uybar*Utau*Ux*Utau*Ux*Utau*Uybar*Utau*Ux*Utau*Uybar*Utau*Ux*Utau*Ux*Utau*Uybar*Utau*Ux*Utau*Ux*Utau*Uy*Utau*Ux*Utau*Ux*Utau*Uybar*Utau*Ux*Utau*Ux*Utau*Uybar*Utau*Uy*Utau*Uy*Utau*Ux*Utau*Uy*Utau*Uy*Utau*Uxbar*Utau;
                    Havg = 0;
                    testCyc = 1;
                    U0 = 1;
                    cycleTime = 49*tau + 48*pulse;
                    
                elseif strcmp(testSequence,'WHH')
                    testUnitary = Utau*Uxbar*Utau*Uy*Utau*Utau*Uybar*Utau*Ux*Utau;
                    Havg = (Delta / 3) * (X + Y + Z); 
                    testCyc = 12;
                    cycleTime = 4 * pulse + 6 * tau;
                    U0 = expm(-1i*Havg*2*pi*cycleTime);
                end

        % CALCULATE AVERAGE HAMILTONIAN TERMS

                sequence = getSequence(testSequence);
                Pulses = sequence.Pulses;
                Taus = tau * sequence.Taus;
                tCyc = sum(Taus);

                toggledHsys = {};
                for p = 0:length(Pulses)
                    toggledHsys{p+1} = getURF(p)'*Hint*getURF(p);
                end

                %Calculate 0th order Magnus term
                H0 = zeros(dim,dim);
                for k=0:length(Pulses)
                    URF = getURF(k);
                    Htilde = (URF')*Hint*URF;
                    Taus(k+1);
                    H0 = H0 + (Taus(k+1))*Htilde;
                end

                H0 = (1/tCyc)*H0;


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

                Havg = H0 + H1 + H2;

                Uavg = expm(-1i*Havg*2*pi*cycleTime);


        % SPIN SIMULATION - - - - - - - - -

                UcumExp = testUnitary;  % experimental
                UcumH0 = U0;            % target (H0 for spectroscopic, id for time-suspension)
                UcumHavg = Uavg;        % first 3 terms

                for cycleCount=1:Ncyc
                    fidelity0 = zeros(Ncyc, 1);
                    fidelity3 = zeros(Ncyc, 1);

                    %Unitary Evolution

                    simTime = 0;

                    % fidelity metric for unitary operators
                    fidelity0(cycleCount) = metric(UcumExp, UcumH0, N); %exp compared to order 0
                    fidelity3(cycleCount) = metric(UcumExp, UcumHavg, N); %exp compared to orders 0-2

                    UcumExp = UcumExp * testUnitary;
                    UcumH0 = UcumH0 * U0;
                    UcumHavg = UcumHavg * Uavg;

                end      

                raw_results(sq,var,dipVal,testVal)=fidelity0(Ncyc);
                raw_results_2(sq,var,dipVal,testVal)=fidelity3(Ncyc);

            end
        end

        for d=1:40
            for aa=1:8
                results(sq,var,d) = results(sq,var,d)+(1/8)*raw_results(sq,var,d,aa);
                results_2(sq,var,d) = results_2(sq,var,d)+(1/8)*raw_results_2(sq,var,d,aa);
            end
        end
    end
end



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
        sequence.Taus = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
        sequence.Taus(1)=0;
    end
    
end

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

as={a1,a2,a3,a4,a5,a6,a7, a8};

end

function fidelity = metric(Utarget, Uexp, N)

fidelity = abs(trace(Utarget' * Uexp)/2^N);

end
