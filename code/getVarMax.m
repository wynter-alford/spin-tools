%% getVarMax.m
% Wynter Alford
% February 2022
%
% Returns the standard testing range for parameter testVarName. 'lo' options
% exist for the coupling and delta (offset) parameters to test over a
% smaller range.  coupling_lo has become the default testing range for the
% coupling strength.

function testValueMax = getVarMax(testVarName,rangeMod)
    if ~exist("rangeMod",'var')
        rangeMod = 1;
    end
    if strcmp(testVarName,'tau')||strcmp(testVarName,'Tau')
        testValueMax = 10e-6;
    elseif strcmp(testVarName,'Delta')||strcmp(testVarName,'delta')
        testValueMax = 1000;
    elseif strcmp(testVarName,'Coupling')||strcmp(testVarName,'coupling')
        testValueMax = 50000;
    elseif strcmp(testVarName,'coupling_lo')
        testValueMax = 8000;
    elseif strcmp(testVarName,'Overrotations')||strcmp(testVarName,'overrotations')
        testValueMax = 0.03;
    elseif strcmp(testVarName,'phaseTrans')
        testValueMax = 0.15;
    elseif strcmp(testVarName,'delta_lo')
        testValueMax = 400;
    elseif strcmp(testVarName,'overrot_hi')
        testValueMax = 0.2;
    end
    testValueMax = testValueMax * rangeMod;
end