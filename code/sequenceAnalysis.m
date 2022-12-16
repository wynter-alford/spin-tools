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
tic
%% Control Parameters (change these)

% Valid sequences are: WHH, MREV8, CORY48, YXX48, YXX24, YXX24S, AZ48,
%    AS10M, AS10E, WHHWHH, ML10, WPW9-2Cycle, MG8, SED24, SED28, SED96

sequenceName = 'WHH' ;  % select sequence to test over.
testVarName = 'WT1'; % select parameter to test over (tau, tau_lo, coupling, coupling_lo, delta, delta_lo, phaseTrans, overrot, overrot_hi) or use 'WT1' to obtain a single data point for fixed parameters
varMaxMod = 1; % factor by which to multiply the default testVarMax


% Procedure Parameters
do_Magnus = true;
do_Floquet = true;

% Magnus Parameters
maxTerm = 0; % highest Magnus series term to compute 

% Floquet Parameters
maxFloqN = 30; % highest Fourier harmonic to include


% Iteration Parameters
testValueCount = 1;
couplingsCount = 1; % how many different coupling matrices to average over
pulseDivs = 1; % how many intervals to divide pulses into for finite pulse AHT


% System Parameters
N = 3;
pulse = 1.4e-6;
tau = 4e-6;  % delay spacing
coupling = 8000; % Hz (not rad/s). omega_D=coupling*2*pi
Delta = 000;
overRotation = 0; %0 is a perfect pulse, +/- 0.01 is a 1% over/under rotation, etc.
phaseTrans = 0; %0 is no phase transient, 1 is a pi/2 phase transient


% File Saving Parameters
fileNameAdd = '';
saveRuns = false; % if true, the workspace will be saved as a .mat file when the script concludes.
%% Analyze Sequence

runAnalysis
toc