function fidelity = metric(Utarget, Uexp, N)

fidelity = abs(trace(Utarget' * Uexp)/2^N);

end