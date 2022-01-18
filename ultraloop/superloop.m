
if strcmp(outerVarName,'tau')||strcmp(outerVarName,'Tau')
    outerValueMax = 10e-6;
elseif strcmp(outerVarName,'Delta')||strcmp(outerVarName,'delta')
    outerValueMax = 1000;
elseif strcmp(outerVarName,'Coupling')||strcmp(outerVarName,'coupling')
    outerValueMax = 50000;
elseif strcmp(outerVarName,'coupling_lo')
    outerValueMax = 15000;
elseif strcmp(outerVarName,'Overrotations')||strcmp(outerVarName,'overrotations')
    outerValueMax = 0.03;
elseif strcmp(outerVarName,'phaseTrans')
    outerValueMax = 0.15;
elseif strcmp(outerVarName,'delta_lo')
    outerValueMax = 400;
elseif strcmp(outerVarName,'overrot_hi')
    outerValueMax = 0.2;
end

tauTestVars = zeros(outerValueCount,1);
fitparameters = struct;
fidelities = struct;

for superind = 1:outerValueCount
    tau = (superind/outerValueCount)*outerValueMax;
    tauTestVars(superind)=tau;
    sequenceAnalysisRepeatable
    fitf0 = fitAW0(testVars,results_f(:,1));
    fitDf0 = fitAW0(testVars,results_Df(:,1));
    fitDfn = fitAW0(testVars,results_Df(:,maxTerm+1));
    
    fitparameters.W0(superind,:) = [fitf0.W0,fitDf0.W0,fitDfn.W0];
    fitparameters.a(superind,:) = [fitf0.a,fitDf0.a,fitDfn.a];
    fitparameters.b(superind,:) = [fitf0.b,fitDf0.b,fitDfn.b];

    fidelities.results_f{superind} = results_f;
    fidelities.results_Df{superind} = results_Df;
end

Date = date;
filename = strcat(loc,'/',date,"_",sequenceName,"_",string(Delta),"_superloop_couplingANDtauFitParameters.mat");
save(filename,'fitparameters','Delta','pulse','sequenceName','testValueCount','outerValueCount','testVars','tauTestVars','maxTerm','N','Date');