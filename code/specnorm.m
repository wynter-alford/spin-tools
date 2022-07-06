%% specnorm.m
% Wynter Alford
% April 2021
%
% Computes the spectral norm (2-norm) of a matrix.

function norm = specnorm(A) % Spectral Norm (sqrt maximum eigenvalue of A'A)
    mat = A'*A;
    evals = eig(mat);
    norm = sqrt(max(evals));
end