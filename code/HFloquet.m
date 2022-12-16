%% HFloquet.m
% Wynter Alford
% December 2022
%
% Computes the Floquet Hamiltonian (in the expanded Floquet space) for the
% piecewise time-independent, periodic Hamiltonian described by
% toggledHsys, Taus. Truncates Fourier modes at Nmax (technically, modes up
% to 2*Nmax remain relevant).

function lambda0 = HFloquet(Nmax,Taus,toggledHsys,dim)
    tCyc = sum(Taus);    
    hf = (2*pi/tCyc)*kron(NFloq(Nmax),eye(dim,dim));
    
    for n=-Nmax:Nmax
        hf = hf + kron(FFloq(n,Nmax),HFourier(n,Taus,toggledHsys));
    end

    [~,diagHF]=eig(hf);    
   
    lambda0 = diagHF(Nmax*dim+1:(Nmax+1)*dim,Nmax*dim+1:(Nmax+1)*dim);

end

function nf = NFloq(Nmax) % COMPLETE
    nf = zeros(2*Nmax+1,2*Nmax+1);    
    for j=-Nmax:Nmax
        nf(j+Nmax+1,j+Nmax+1)=j;
    end    
end

function ff = FFloq(n, Nmax)

    ff = zeros(2*Nmax+1,2*Nmax+1);
    if n==0
        ff = speye(2*Nmax+1,2*Nmax+1);
        
    elseif n<0
        for j=-Nmax-n:Nmax
            ff(j+Nmax+1,j+Nmax+1+n)=1;
        end
        
    elseif n>0
        for j=-Nmax+n:Nmax
            ff(j+Nmax+1-n,j+Nmax+1)=1;
        end
        
    end
end