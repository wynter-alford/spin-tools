%% Q.m
% Wynter Alford
% January 2022
%
% Part of the recursive Magnus Expansion calculation.  This is an auxilary
% function with a role in calculating the Omega(n) terms.

function qnk = Q(n,k,toggledHsys,Taus)
    global knownQs %#ok<*GVMIS> 

    % Check whether this term has already been computed
    if length(knownQs)>=n && length(knownQs{n})>=k && ~isempty(knownQs{n}{k})
        qnk = knownQs{n}{k};
        
    % if it hasn't been, then compute it        
    else      
        if k==1
            qnk = Omega(n,toggledHsys,Taus);
            
        elseif k==n
            qnk = (Omega(1,toggledHsys,Taus))^n;
            
        else
            qnk = 0;
            for m=1:(n-k+1) 
                qnk = qnk + Omega(m,toggledHsys,Taus)*Q(n-m,k-1,toggledHsys,Taus);
            end            
        end
        
        knownQs{n}{k}=qnk;
    end
end