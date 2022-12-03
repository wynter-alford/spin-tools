%% HFourier.m
% Wynter Alford
% December 2022
% 
function hf = HFourier(n,Taus,toggledHsys)
    tCyc = sum(Taus);
    intFn = @(v) exp(1i*n*v*(2*pi/tCyc));
    hf = 0;
    intStart = 0;

    for t=1:length(toggledHsys)
        hf = hf + toggledHsys{t}*integral(intFn,intStart,intStart+Taus(t));
        intStart = intStart + Taus(t);
    end
    hf = (1/tCyc)*hf;
end