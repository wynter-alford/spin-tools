if curse == 1
    sequenceNames = {'WHH','MREV8','BR24','CORY48'};
    loc = 'airResults';
elseif curse == 0.1
    sequenceNames = {"WHH"};
elseif curse == 4
    sequenceNames = {'MREV8','ML10'};
    loc = 'winResults';
else
    sequenceNames = {'ML10','YXX24','YXX48'};
    loc = 'winResults';
end

%% Setup
global dim pulse f1 Pulses X Y Z getNextUsM %#ok<GVMIS,NUSED> 

% Valid sequences are: WHH, MREV8, CORY48, YXX48, YXX24, YXX24S, AZ48,
%    BR24, AS10M, AS10E, WHHWHH, ML10, WPW9-2Cycle, XPulse, YPulse

%sequenceName = 'MREV8' ;  % select sequence to test over.
testVarName = 'coupling'; % select parameter to test over (Delta, Tau, Coupling, coupling_lo, phaseTrans, delta_lo, 'overrot_hi')
outerVarName = 'tau';
maxTerm = 6; % highest Magnus series term to compute (can do 8+ for WHH; no higher than 4 for 48 pulses)

% Control Parameters
testValueCount = 20;
outerValueCount = 20;
couplingsCount = 4; % how many different coupling matrices to average over

N = 4;
dim = 2^N;
pulse = 1.4e-6;
tau = 5.4e-6;  % delay spacing
coupling = 5000;
Delta = 500;
overRotation = 0; %0 is a perfect pulse, +/- 0.01 is a 1% over/under rotation, etc.
phaseTrans = 0; %0 is no phase transient, 1 is a pi/2 phase transient

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

%% Initialize Hamiltonians (modified code)

 Hdips = cell(couplingsCount,1);
 % Generate [couplingsCount] dipole Hamiltonians, each with a different coupling matrix
 for j=1:couplingsCount
     dip = abs(randn(N));
     dip = triu(dip,1) + triu(dip,1)';
     Hdips{j} = getHdip(N, dim, x, y, z, dip);
 end

%% Loop over sequences and Deltas

for uloop = 1:length(sequenceNames)
    sequenceName = sequenceNames{uloop};
    maxTerm = getMaxTerm(sequenceName);
    for udloop = 1:3
        Delta = 500*(udloop-1);
        superloop;
    end
end