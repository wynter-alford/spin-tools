%% indexfactor.m
% Wynter Alford
% January 2022
% 
% I'm not quite sure what this function does (I'm writing this in July
% 2022 as I forgot to in January when I made this function) but it is an
% important part of the Magnus Series computations (see P.m for its usage).

% Note that indexList must be SORTED in decreasing order for this function
% to work properly.

function indf = indexfactor(indexList)

    matched = 1;
    indf = 1;
    
    for i=2:length(indexList)
        if indexList(i)==indexList(i-1)
            matched = matched + 1;
            indf = indf*matched;
        else
            matched = 1;
        end        
    end 
    
end