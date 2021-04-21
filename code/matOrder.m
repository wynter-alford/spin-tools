% Ben Alford, November 2020
% Calculates the order of magnitude of a matrix, as measured by the log10
% of the square root of the trace of A*A.

function orh = matOrder(A) % calculates the magnitude of a matrix
    orh = log10(sqrt(trace(A*A)));
end