function df = mdiff(Utarget,Uexp,N)
    df = abs(trace(Utarget' * Uexp)/2^N) - real(trace(Utarget' * Uexp)/2^N);
end