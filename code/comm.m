% Ben Alford, November 2020
% Calculates the commutator of two matrices (or other multipliable objects)

function ct = comm(A,B) % calculates the commutator of a pair of matrices
    ct = A*B-B*A;
end