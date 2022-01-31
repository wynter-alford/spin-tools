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