%% tracenorm.m
% Wynter Alford
% November 2020
%
% Computes the trace norm (1-norm) of a matrix A

function result = tracenorm(A)
    result = sqrt(trace(A*A));
end