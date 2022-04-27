%clear

%% Control Parameters (change these)

global dim Pulses knownOmegas knownQs knownPs 

% Valid sequences are: WHH, MREV8, CORY48, YXX48, YXX24, YXX24S, AZ48,
%    AS10M, AS10E, WHHWHH, ML10, WPW9-2Cycle, MG8

sequenceName = 'WHH' ;  % select sequence to test over.
testVarName = 'coupling_lo'; % select parameter to test over (Delta, Tau, Coupling, coupling_lo, phaseTrans, delta_lo, 'overrot_hi')
varMaxMod = 1; % factor by which to multiply the default testVarMax

mode = 'max'; %can be 'max' to get up to maxTerm terms, or 'time' to compute for a certain amount of time
maxTerm = 6; % highest Magnus series term to compute 
computationTime = 120; % once elapsed time reaches this many seconds, no further terms are computed for this loop.

% Control Parameters
testValueCount = 20;
couplingsCount = 4; % how many different coupling matrices to average over

N = 8;
pulse = 1.4e-6;
tau = 7.4e-6;  % delay spacing
coupling = 5000;
Delta = 00;
overRotation = 0; %0 is a perfect pulse, +/- 0.01 is a 1% over/under rotation, etc.
phaseTrans = 0; %0 is no phase transient, 1 is a pi/2 phase transient

fileNameAdd = '_6Spins';
saveRuns = true;
%% Analyze Sequence

runAnalysis;