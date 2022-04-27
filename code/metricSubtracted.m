function fidelity = metricSubtracted(U1,U2,N)
    fidelity = tracenorm(U1-U2)/2^N;
end