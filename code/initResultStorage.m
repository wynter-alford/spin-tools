%% initResultStorage.m
% Wynter Alford
% Updated December 2022
%
% Constructs result storage arrays for sequenceAnalysis

%% Naming Convention for Result Storage Arrays
% Array naming convention is results_ABC with a 3-letter code.
% 
% The first letter is the AVERAGING METHOD used, which is either the Magnus
% series (AHT) or a truncated Floquet matrix. The averaging method options
% are:
%
% M: Magnus
% Q: Floquet
% S: True Time-Suspension Fidelity (Magnus/Floquet terms not computed)
%
% The second letter is the METRIC in question.  The metric options are:
% O: overlap(Uexp, AHTUnitaries{termInd}) is the propagator overlap of the
% experimental and AHT propagators.
% D: uDist(Uexp, AHTUnitaries{termInd}) is the unitary distance (currently
% using the spectral norm) between the experimental and AHT propagators.
% H: specnorm(MagnusTerms{termInd}) is the norm of the average Hamiltonian
% term
%
% The third letter is the PULSE TREATMENT for constructing Uexp, the
% EXPERIMENTAL propagator. The options are:
% I: Instantaneous pulses (deltaUnitary)
% F: Finite pulses (expUnitary)
% H: used in AHT if the Metric is H, since in that case treatment of Uexp 
% is irrelevant. 
% C: used in Floquet if comparing correction terms, since in that case
% treatment of Uexp is irrelevant.
%
% The fourth letter is the PULSE METHOD for AHT or Floquet, indicating 
% whether the Average Hamiltonian or Floquet computations used pulse 
% divisions. The options are:
% T: Traditional analysis assuming instantaneous pulses. Floquet
% computations use the interaction picture.
% P: Pulse length incorporated into the toggledHsys computations.
% S: Uexp compared only to the identity matrix; AHT/Floquet terms not 
% computed (the true Time-Suspension fidelity) (Only occurs with S first)
%
% For example, results_MOFT is the overlap between testUnitary (Uexp for
% finite pulses) and AHTUnitaries{termInd} (U_AHT,n not accounting for
% pulse widths in the AHT computations).

%% Storage Array Creation

testVars = zeros(testValueCount,1);

% Size of Magnus Terms
results_MHHT = zeros(length(testVars),maxTerm+1);
results_MHHP = zeros(length(testVars),maxTerm+1);

raw_MHHT = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_MHHP = zeros(length(testVars),couplingsCount,maxTerm+1);

% Magnus Overlaps
results_MOIT = zeros(length(testVars),maxTerm+1);
results_MOFT = zeros(length(testVars),maxTerm+1);
results_MOFP = zeros(length(testVars),maxTerm+1);

raw_MOIT = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_MOFT = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_MOFP = zeros(length(testVars),couplingsCount,maxTerm+1);

% Magnus Unitary Distances
results_MDIT = zeros(length(testVars),maxTerm+1);
results_MDFT = zeros(length(testVars),maxTerm+1);
results_MDFP = zeros(length(testVars),maxTerm+1);

raw_MDFT = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_MDIT = zeros(length(testVars),couplingsCount,maxTerm+1);
raw_MDFP = zeros(length(testVars),couplingsCount,maxTerm+1);

% Floquet Correction Term Sizes
results_QOCT = zeros(length(testVars),maxFloqN+1);
results_QDCT = zeros(length(testVars),maxFloqN+1);
results_QOCP = zeros(length(testVars),maxFloqN+1);
results_QDCP = zeros(length(testVars),maxFloqN+1);

raw_QOCT = zeros(length(testVars),couplingsCount,maxFloqN+1);
raw_QDCT = zeros(length(testVars),couplingsCount,maxFloqN+1);
raw_QOCP = zeros(length(testVars),couplingsCount,maxFloqN+1);
raw_QDCP = zeros(length(testVars),couplingsCount,maxFloqN+1);

% Floquet Overlaps
results_QOIT = zeros(length(testVars),maxFloqN+1);
results_QOFT = zeros(length(testVars),maxFloqN+1);
results_QOIP = zeros(length(testVars),maxFloqN+1);
results_QOFP = zeros(length(testVars),maxFloqN+1);

raw_QOIT = zeros(length(testVars),couplingsCount,maxFloqN+1);
raw_QOFT = zeros(length(testVars),couplingsCount,maxFloqN+1);
raw_QOIP = zeros(length(testVars),couplingsCount,maxFloqN+1);
raw_QOFP = zeros(length(testVars),couplingsCount,maxFloqN+1);

% Floquet Unitary Distances
results_QDIT = zeros(length(testVars),maxFloqN+1);
results_QDFT = zeros(length(testVars),maxFloqN+1);
results_QDIP = zeros(length(testVars),maxFloqN+1);
results_QDFP = zeros(length(testVars),maxFloqN+1);

raw_QDFT = zeros(length(testVars),couplingsCount,maxFloqN+1);
raw_QDIT = zeros(length(testVars),couplingsCount,maxFloqN+1);
raw_QDFP = zeros(length(testVars),couplingsCount,maxFloqN+1);
raw_QDIP = zeros(length(testVars),couplingsCount,maxFloqN+1);

% Time-Suspension Overlaps
results_SOFS = zeros(length(testVars),1);
results_SOIS = zeros(length(testVars),1);

raw_SOFS = zeros(length(testVars),couplingsCount);
raw_SOIS = zeros(length(testVars),couplingsCount);

% Time-Suspension Unitary Distances
results_SDFS = zeros(length(testVars),1);
results_SDIS = zeros(length(testVars),1);

raw_SDFS = zeros(length(testVars),couplingsCount);
raw_SDIS = zeros(length(testVars),couplingsCount);