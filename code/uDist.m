%% uDist.m
% author, date not listed
%
% Computes the Hilbert-Schmidt distance between two unitary matrices,
% scaled by the size of the system.

function fidelity = uDist(Utarget, Uexp, N)
    fidelity = specnorm(Utarget-Uexp)/2^N;
end