function Hdip = getHdip(N, dim, x, y, z, a)

% N: number of spins
% x/y/z: Pauli matrices
% a: dipolar interaction strengths (NxN symmetric matrix)

Hdip = sparse(dim, dim);

for k=1:N
    for h=k+1:N
        % Hdip = a_{kh} (3Z - X - Y)
        Hdip=Hdip+a(k,h)*(2*mykron(speye(2^(k-1)),z,speye(2^(h-k-1)),z,speye(2^(N-h)))-... 
            mykron(speye(2^(k-1)),x,speye(2^(h-k-1)),x,speye(2^(N-h)))-...
            mykron(speye(2^(k-1)),y,speye(2^(h-k-1)),y,speye(2^(N-h)))) ;             
    end
end

end