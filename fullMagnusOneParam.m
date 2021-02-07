%% Setup Code (original Spin Sim code)

global dim pulse f1 Pulses X Y Z

sequenceName = 'MREV8';  % select sequence to test over
testVarName = 'Tau'; % only affects file name currently
testValueCount = 30;
testValueMax = 20e-6;

couplingsCount = 6; % how many different coupling matrices to average over

N = 4;
dim = 2^N;
Ncyc = 1;
pulse = 0;
tau = 3e-6;  % delay spacing
coupling = 5000;
f1 = 1/4/pulse; %Can adjust f1 and w1 by changing 'pulse' variable
w1 = 2*pi*f1;
Delta = 500;

% initvars
z=0.5*sparse([1 0; 0 -1]);x=0.5*sparse( [ 0 1;1 0]); y=1i*0.5*sparse([0 -1;1 0]);
ep=sparse([1 0; 0 0]);
em=sparse([0 0; 0 1]);
id=speye(2);p=sparse([0 1;0 0]);m=sparse([0 0; 1 0]);


%initCollectiveObs;
Z=sparse(dim,dim);
X=sparse(dim,dim);
Y=sparse(dim,dim);
for k=1:N
      Z = Z + mykron(speye(2^(k-1)),z,speye(2^(N-k)));
      X = X + mykron(speye(2^(k-1)),x,speye(2^(N-k)));
      Y = Y + mykron(speye(2^(k-1)),y,speye(2^(N-k)));
end



%% Initialize Hamiltonians (original Spin Sim code)

Hdips = cell(couplingsCount,1);
% Generate 8 dipole Hamiltonians, each with a different coupling matrix
for j=1:couplingsCount
    dip = abs(randn(N));
    dip = triu(dip,1) + triu(dip,1)';
    Hdips{j} = getHdip(N, dim, x, y, z, dip);
end

%% Iterate over different Delta values to see how term magnitude changes

testVars = zeros(testValueCount,1);

results_h0 = zeros(length(testVars),1);
results_h1 = zeros(length(testVars),1);
results_h2 = zeros(length(testVars),1);
results_h3 = zeros(length(testVars),1);
results_h4 = zeros(length(testVars),1);

raw_results_h0 = zeros(length(testVars),couplingsCount);
raw_results_h1 = zeros(length(testVars),couplingsCount);
raw_results_h2 = zeros(length(testVars),couplingsCount);
raw_results_h3 = zeros(length(testVars),couplingsCount);
raw_results_h4 = zeros(length(testVars),couplingsCount);


for d=1:length(testVars)
    tau = d*(testValueMax/testValueCount); % ADJUST FOR TEST VAR 
    testVars(d)=tau; % ADJUST FOR TEST VAR
    
    sequence = getSequence(sequenceName);
    Pulses = sequence.Pulses;
    Taus = tau * sequence.Taus;
    tCyc = sum(Taus);
    
    for c=1:couplingsCount
        
        Hdip = Hdips{c};
        Hsys = Hdip*coupling + Z*Delta;
        
        toggledHsys = {};
        for p = 0:length(Pulses)
            toggledHsys{p+1} = getURF(p)'*Hsys*getURF(p);
        end

        %calcMagnus
        %Magnus Calculation is done here ------------------------------------

        %Calculate 0th order Magnus term
        H0 = zeros(dim,dim);
        for k=0:length(Pulses)
            URF = getURF(k);
            Htilde = (URF')*Hsys*URF;
            Taus(k+1);
            H0 = H0 + (Taus(k+1))*Htilde;
        end

        H0 = (1/tCyc)*H0;
        h0 = real(matOrder(H0));


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
        h1 = real(matOrder(H1));

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
        h2 = real(matOrder(H2));

        % Calculate 3rd order Magnus term
        H3 = zeros(dim,dim);
        for m=0:length(Pulses)
            for l=0:m
                for k=0:l
                    for j=0:k
                        Hm = toggledHsys{m+1};
                        Hl = toggledHsys{l+1};
                        Hk = toggledHsys{k+1};
                        Hj = toggledHsys{j+1};

                        term1 = comm(comm(comm(Hm,Hl),Hk),Hj);
                        term2 = comm(Hm,comm(comm(Hl,Hk),Hj));
                        term3 = comm(Hm,comm(Hl,comm(Hk,Hj)));
                        term4 = comm(Hl,comm(Hk,comm(Hj,Hm)));

                        Hterm = term1+term2+term3+term4;
                        tauProd = Taus(m+1)*Taus(l+1)*Taus(k+1)*Taus(j+1);

                        H3 = H3 + Hterm*tauProd;
                    end
                end
            end
        end

        H3 = (-1/(12*1i*tCyc))*H3;
        h3 = real(matOrder(H3));

        % Calculate the 4th term of the Magnus Expansion
        H4 = zeros(dim,dim);
        for mm=0:length(Pulses)
            for m=0:mm
                for l=0:m
                    for k=0:l
                        for j=0:k
                            % Express Hsys(t) in the interaction frame
                            Hmm = toggledHsys{mm+1};
                            Hm = toggledHsys{m+1};
                            Hl = toggledHsys{l+1};
                            Hk = toggledHsys{k+1};
                            Hj = toggledHsys{j+1};

                            % Calculate commutators
                            term1 = (-1/30)*comm(Hmm,comm(Hm,comm(Hl,comm(Hk,Hj))));
                            term2 = (2/15)*comm(Hj,comm(Hmm,comm(Hm,comm(Hk,Hl))));
                            term3 = (1/15)*comm(comm(Hmm,Hj),comm(Hm,comm(Hk,Hl)));
                            term4 = (1/15)*comm(comm(Hm,Hj),comm(Hmm,comm(Hk,Hl)));
                            term5 = (-1/60)*comm(comm(Hk,Hl),comm(Hmm,comm(Hm,Hj)));
                            term6 = (1/15)*comm(comm(Hl,Hj),comm(Hmm,comm(Hk,Hm)));
                            term7 = (-1/60)*comm(comm(Hk,Hm),comm(Hmm,comm(Hl,Hj)));
                            term8 = (-1/60)*comm(comm(Hk,Hmm),comm(Hm,comm(Hl,Hj)));
                            term9 = (-1/60)*comm(comm(Hl,Hm),comm(Hmm,comm(Hk,Hj)));
                            term10 = (-1/60)*comm(comm(Hl,Hm),comm(Hj,comm(Hk,Hmm)));
                            term11 = (-1/60)*comm(comm(Hmm,Hj),comm(Hl,comm(Hk,Hm)));
                            term12 = (-1/60)*comm(comm(Hm,Hj),comm(Hl,comm(Hk,Hmm)));
                            term13 = (-1/60)*comm(comm(Hl,Hmm),comm(Hm,comm(Hk,Hj)));
                            term14 = (-1/60)*comm(comm(Hl,Hmm),comm(Hj,comm(Hk,Hm)));
                            term15 = (-1/30)*comm(Hj,comm(Hm,comm(Hl,comm(Hk,Hmm))));
                            term16 = (-1/60)*comm(comm(Hm,Hmm),comm(Hj,comm(Hk,Hl)));
                            term17 = (-1/60)*comm(comm(Hk,Hl),comm(Hj,comm(Hm,Hmm)));
                            term18 = (-1/60)*comm(comm(Hk,Hm),comm(Hj,comm(Hl,Hmm)));
                            term19 = (-1/60)*comm(comm(Hk,Hj),comm(Hm,comm(Hl,Hmm)));
                            term20 = (-1/60)*comm(comm(Hm,Hmm),comm(Hl,comm(Hk,Hj)));
                            term21 = (-1/60)*comm(comm(Hl,Hj),comm(Hm,comm(Hk,Hmm)));
                            term22 = (-1/30)*comm(Hj,comm(Hmm,comm(Hl,comm(Hk,Hm))));

                            % Add commutators
                            Hterm = term1+term2+term3+term4+term5+term6+term7+term8+term9+term10+term11+term12+term13+term14+term15+term16+term17+term18+term19+term20+term21+term22;
                            tauProd = Taus(mm+1)*Taus(m+1)*Taus(l+1)*Taus(k+1)*Taus(j+1);

                            H4 = H4 + Hterm*tauProd;
                        end
                    end
                end
            end
        end

        H4 = (1/tCyc)*H4;
        h4 = real(matOrder(H4));


        % Magnus Calculation ends here ---------------------------------------


        raw_results_h0(d,c) = h0;
        raw_results_h1(d,c) = h1;
        raw_results_h2(d,c) = h2;
        raw_results_h3(d,c) = h3;
        raw_results_h4(d,c) = h4;
    end
    results_h0(d)=mean(raw_results_h0(d,:));
    results_h1(d)=mean(raw_results_h1(d,:));
    results_h2(d)=mean(raw_results_h2(d,:));
    results_h3(d)=mean(raw_results_h3(d,:));
    results_h4(d)=mean(raw_results_h4(d,:));
end

%% Save Result Output
fileDescriptor = strcat(sequenceName,'_',testVarName,'_magnus_results_',date,'.mat');
save(strcat(testVarName,fileDescriptor), 'results_h0', 'results_h1', 'results_h2', 'results_h3', 'results_h4', 'testVars')
% save(strcat('h0_',fileDescriptor),'results_h0');
% save(strcat('h1_',fileDescriptor),'results_h1');
% save(strcat('h2_',fileDescriptor),'results_h2');
% save(strcat('h3_',fileDescriptor),'results_h3');
% save(strcat('h4_',fileDescriptor),'results_h4');
% save(strcat(testVarName,'s_',fileDescriptor),'testVars');

%% FUNCTION DEFINITIONS
function ct = comm(A,B) % calculates the commutator of a pair of matrices
    ct = A*B-B*A;
end

function orh = matOrder(A) % calculates the magnitude of a matrix
    orh = log10(sqrt(trace(A*A)));
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
    end
end

