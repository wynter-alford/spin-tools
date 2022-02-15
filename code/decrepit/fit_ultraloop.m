% if curse == 1
%     sequenceNames = {'WHH','MREV8','BR24','CORY48'};
%     loc = 'airResults';
% elseif curse == 0.1
%     sequenceNames = {"WHH"};
% elseif curse == 4
%     sequenceNames = {'MREV8','ML10'};
%     loc = 'winResults';
% else
%     sequenceNames = {'ML10','YXX24','YXX48'};
%     loc = 'winResults';
% end

loc = '2022-01-31';
sequenceNames = {'WHH','MREV8','MG8','ML10','BR24','CORY48','YXX24','YXX48'};

%% Setup
global dim pulse f1 Pulses %#ok<GVMIS,NUSED> 

% Valid sequences are: WHH, MREV8, CORY48, YXX48, YXX24, YXX24S, AZ48,
%    BR24, AS10M, AS10E, WHHWHH, ML10, WPW9-2Cycle, XPulse, YPulse

%sequenceName = 'MREV8' ;  % select sequence to test over.
testVarName = 'tau'; % select parameter to test over (Delta, Tau, Coupling, coupling_lo, phaseTrans, delta_lo, 'overrot_hi')
outerVarName = 'tau';
maxTerm = 0; % highest Magnus series term to compute (can do 8+ for WHH; no higher than 4 for 48 pulses)

% Control Parameters
testValueCount = 30;
outerValueCount = 1;
couplingsCount = 4; % how many different coupling matrices to average over

N = 4;
dim = 2^N;
pulse = 1.4e-6;
tau = 5.4e-6;  % delay spacing
coupling = 5000;
Delta = 0;
overRotation = 0; %0 is a perfect pulse, +/- 0.01 is a 1% over/under rotation, etc.
phaseTrans = 0; %0 is no phase transient, 1 is a pi/2 phase transient

initVars
initCollectiveObs

outerValueMax = getVarMax(outerVarName);
testValueMax = getVarMax(testVarName,1);

% Derived parameters
dim = 2^N;
f1 = 1/4/pulse; %Can adjust f1 and w1 by changing 'pulse' variable
w1 = 2*pi*f1;
rotationError = 1+overRotation; %1 is a perfect pulse
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
    %maxTerm = getMaxTerm(sequenceName);
    maxTerm = 2;
    for udloop = 1:3
        Delta = 500*(udloop-1);
        fit_superloop;
    end
end
