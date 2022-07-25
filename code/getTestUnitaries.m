%% getTestUnitaries.m
% Wynter Alford
% February 2022
% 
% Computes the experimental unitary operator for a pulse sequence, with
% both finite pulses (testUnitary) and instantaneous pulses (deltaUnitary)

UtauCell = {speye(dim,dim), Utau, U2tau};
UDtauCell = {speye(dim,dim), UDtau, UD2tau};

expUnitary = speye(dim,dim); 
deltaUnitary = speye(dim,dim);

for p=1:length(Pulses)

    Uint = UtauCell{Taus(p)/tau+1};
    UDint = UDtauCell{Taus(p)/tau+1};

    if sequence.Pulses{p}==X
        nextU = Ux;
        nextUD = UDx;
    elseif sequence.Pulses{p}==Y
        nextU = Uy;
        nextUD = UDy;
    elseif sequence.Pulses{p}==-X
        nextU = Uxbar;
        nextUD = UDxbar;
    elseif sequence.Pulses{p}==-Y
        nextU = Uybar;
        nextUD = UDybar;
    else
        nextU = expm(-1i*Pulses{p}*rotationError*pi/2);
        nextUD = expm(-1i*2*pi*(Hsys+f1*Pulses{p}*rotationError)*pulse);
    end

    expUnitary = nextU * Uint * expUnitary;
    deltaUnitary = nextUD * UDint * deltaUnitary;
end

expUnitary = UtauCell{Taus(length(Taus))/tau+1} * expUnitary;
deltaUnitary = UDtauCell{Taus(length(Taus))/tau+1} * deltaUnitary;