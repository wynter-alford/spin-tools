% Spectral norm of a matrix
% Ben Alford
% April 2021

function norm = specnorm(A) % Spectral Norm (sqrt maximum eigenvalue of A'A)
    mat = A'*A;
    evals = eig(mat);
    norm = sqrt(max(evals));
end