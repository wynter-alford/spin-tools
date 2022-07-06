%% overlap.m
% author, date not listed
%
% Computes the overlap between two unitary matrices, scaled by system size
% Formerly called metric.m

function fidelity = overlap(Utarget, Uexp, N)
    fidelity = abs(trace(Utarget' * Uexp)/2^N);
end