%% P.m
% Wynter Alford
% January 2022
%
% Part of the recursive Magnus Expansion calculation. This is the function
% that does actual computation, evaluating an integral (as a discrete sum)

function pn = P(n,toggledHsys,Taus)
    global dim knownPs %#ok<*GVMIS> 

    % Check whether this term has already been computed
    if length(knownPs)>=n && ~isempty(knownPs{n})
        pn = knownPs{n};
    
    % if it hasn't been, then compute it        
    else
        pn = zeros(dim,dim);
        times = ones(n,1);
        maxTime = length(toggledHsys);
        k = n;

        % This while loop enacts n summations and is equivalent to:
        % "for i1=1:N, for i2=1:i1 ... for i(n)=1:i(n-1), do stuff"
        while times(1)<=maxTime

            if k==1||times(k)<times(k-1)
                hprod = 1;
                for hnum = 1:n
                    % multiply -iH(t1)*...*-iH(tn)
                    hprod = hprod*(-1i)*toggledHsys{times(hnum)}*Taus(times(hnum));
                end
                
                pn = pn+hprod/indexfactor(times);

                times(k) = times(k)+1;

                if k<length(times)
                    times(k+1:n)=1;
                end
                k=n;
            else
                k = k-1;
            end
        end
        knownPs{n}=pn;
    end
end