%% anticomm.m
% Wynter Alford
% January 2022
%
% Computes the anticommutator of two matrices A and B

function ac = anticomm(A,B)
    ac = A*B+B*A;
end