function modRes = timeAdjust(T, taus, result, tauCount)
    if length(taus)~=length(result)
        error('taus and result must have the same length')
    elseif T<0
        error('T must be positive')
    else
        modRes=zeros(length(result),1);
        cycleTime = taus*tauCount;        
        for j=1:length(result)
            modRes(j)=result(j).^(T/cycleTime(j));
        end
    end
end