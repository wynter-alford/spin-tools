clear
sequenceNames = {'BR24'};%,'BR24','CORY48','SED48'};

for seqInd=1:length(sequenceNames)
    sequenceName = sequenceNames{seqInd};
    for varInd = 1:10
        tau = (0.2+varInd/6) * 10^(-6);
        sequenceAnalysis
        f=fit(log10(testVars(3:20)),-log10(1-results_Df(3:20,1)),'poly1');
        slopes(varInd)=f.p1;
        inters(varInd)=f.p2;
        fitTaus(varInd)=tau;
    end
    save(strcat(date,'_',sequenceName,'_TauSlopeTest.mat'))
end
