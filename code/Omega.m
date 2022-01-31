%% Omega.m
% Wynter Alford
% January 2022
%
% Part of the recursive Magnus expansion calculation method. The terms are
% such that H^(n-1)=Omega(n)*(1i/tCyc).

function omega = Omega(n,toggledHsys,Taus)
    global knownOmegas %#ok<GVMIS> 

    % Check whether this term has already been computed
    if length(knownOmegas)>=n && ~isempty(knownOmegas{n})
        omega=knownOmegas{n};

    % if it hasn't been, then compute it
    else
        omega = P(n,toggledHsys,Taus);
        if n>1
            for k=2:n
                omega = omega - (1/factorial(k))*Q(n,k,toggledHsys,Taus);
            end
        end
        knownOmegas{n}=omega;
    end
end


