%% getEvenTerms.m
% Wynter Alford
% April 2022
%
% Extracts a sequence of every other item of an array.  If `start` is set
% to 0, the first term extracted is array(1), followed by array(3) and so
% on. If `start` is set to 1, the first term extracted is array(2), and so
% on.

function evenTerms = getEvenTerms(array, start)
    if start ~=1 && start ~= 0
        error("'start' parameter must be 1 or 0")
    end

    evenTerms=zeros(ceil(length(array)/2)-start,1);
    
    for i=1:length(evenTerms)
        evenTerms(i)=array(2*i-1+start);
    end
    
end