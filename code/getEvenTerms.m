function evenTerms = getEvenTerms(array, start)
    if start ~=1 && start ~= 0
        error("'start' parameter must be 1 or 0")
    end

    evenTerms=zeros(ceil(length(array)/2)-start,1);
    
    for i=1:length(evenTerms)
        evenTerms(i)=array(2*i-1+start);
    end
    
end