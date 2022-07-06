sequenceNames = {'WHH','MREV8','MG8','ML10','BR24','YXX24','YXX48','CORY48','I24'};
%maxTerms = [10 8 8 8 6 6 4 4];
maxTerms = ones(length(sequenceNames));

% Only the two vars above and tau need be commented out in sequenceAnalysis
% for full functionality.

for config = 1:3
    if config == 1
        testVarName = 'tau';
    elseif config == 2
        testVarName = 'coupling_lo';
        tau = 7.4e-6;
    elseif config == 3
        testVarName = 'coupling_lo';
        tau = 3.4e-6;
    end
    for superind = 1:length(sequenceNames)
        sequenceName = sequenceNames{superind};
        maxTerm = maxTerms(superind);
        sequenceAnalysis
    end
end
