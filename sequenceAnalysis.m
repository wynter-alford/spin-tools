%% sequenceAnalysis.m
% Wynter Alford
% Original version from April 2020. Last updated July 2022. 
%
% Analyze pulse sequences in terms of the Magnus Expansion. WhileThe actual 
% analysis is found in the file runAnalysis.m, this is the setup code.
% The numbers in this script can be easily modified to run different
% simulations or computations.  This script is all that any user of the
% code need interact with except when adding new sequences; new sequences
% are added by modifying the function getSequence.m

clear all %#ok<CLALL>
global dim Pulses knownOmegas knownQs knownPs  %#ok<NUSED>

%% Control Parameters (change these)

load('dij_1_Dipole_Hamiltonians(4).mat')
% Valid sequences are: WHH, MREV8, CORY48, YXX48, YXX24, YXX24S, AZ48,
%    AS10M, AS10E, WHHWHH, ML10, WPW9-2Cycle, MG8, SED24, SED28, SED96

sequenceName = 'MREV8' ;  % select sequence to test over.
testVarName = 'WT1'; % select parameter to test over (tau, tau_lo, coupling, coupling_lo, delta, delta_lo, phaseTrans, overrot, overrot_hi) or use 'WT1' to obtain a single data point for fixed parameters
varMaxMod = 1; % factor by which to multiply the default testVarMax

mode = 'max'; %can be 'max' to get up to maxTerm terms, or 'time' to compute for a certain amount of time
maxTerm = 6; % highest Magnus series term to compute 
computationTime = 120; % once elapsed time reaches this many seconds, no further terms are computed for this loop.

% Control Parameters
testValueCount = 1;
couplingsCount = 1; % how many different coupling matrices to average over
pulseDivs = 2; % how many intervals to divide pulses into for finite pulse AHT

N = 4;
pulse = 1.4e-6;
tau = 4e-6;  % delay spacing
coupling = 4000; % Hz (not rad/s). omega_D=coupling*2*pi
Delta = 0;
overRotation = 0; %0 is a perfect pulse, +/- 0.01 is a 1% over/under rotation, etc.
phaseTrans = 0; %0 is no phase transient, 1 is a pi/2 phase transient

fileNameAdd = '';
saveRuns = false; % if true, the workspace will be saved as a .mat file when the script concludes.
%% Analyze Sequence

runAnalysis