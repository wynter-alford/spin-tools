%% Colors
myColors = {[230 25 75],[245 130 48], [210 245 60], [60 180 75],  [0 130 200], [145 30 180]};
% colors are:   red          orange      lime          green           blue         violet
for i=1:length(myColors)
    myColors{i} = myColors{i}/255;
end

%% Plot Fidelity vs Tau
hold on
plot(testVars,results_f(:,1),'Color',myColors{1})
plot(testVars,results_Df(:,1),'Color',myColors{4})
plot(testVars,results_Df(:,maxTerm+1),'Color',myColors{5})
plot(fitAW0(testVars,results_Df(:,maxTerm+1)))

legend('Uncorrected','Pulse Correction','Pulse and HOTCorrection')
h = xlabel('Tau Spacing (s)');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('Fidelity');
set(h,'interpreter','latex','fontsize',14);
h = title(strcat(sequenceName,'_',string(Delta)));
set(h,'interpreter','latex','fontsize', 16);


%% Plot Infidelity vs Tau
hold on
%plot(testVars,-log10(1-results_f(:,1)),'Color',myColors{1})
plot(testVars,-log10(1-results_Df(:,1)),'Color',myColors{4})
plot(testVars,-log10(1-results_Df(:,maxTerm+1)),'Color',myColors{5})
legend('Pulse Correction','Pulse and HOTCorrection')
h = xlabel('Tau Spacing (s)');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('Infidelity');
set(h,'interpreter','latex','fontsize',14);
h = title(strcat(sequenceName,'_',string(Delta)));
set(h,'interpreter','latex','fontsize', 16);


%% Plot Fidelity vs Coupling
hold on
plot(testVars,results_f(:,1),'Color',myColors{1})
plot(testVars,results_Df(:,1),'Color',myColors{4})
plot(testVars,results_Df(:,maxTerm+1),'Color',myColors{5})
%plot(fitAW0(testVars,results_Df(:,maxTerm+1)))
legend('Uncorrected','Pulse Correction','Pulse and HOTCorrection')
h = xlabel('Coupling Strength ($s^{-1}$)');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('Fidelity');
set(h,'interpreter','latex','fontsize',14);
h = title(strcat(sequenceName,'_',string(Delta)));
set(h,'interpreter','latex','fontsize', 16);

%% Plot Infidelity vs Coupling
hold on
plot(testVars,-log10(1-results_f(:,1)),'Color',myColors{1})
plot(testVars,-log10(1-results_Df(:,1)),'Color',myColors{4})
plot(testVars,-log10(1-results_Df(:,maxTerm+1)),'Color',myColors{5})
legend('Uncorrected','Pulse Correction',strcat("Pulse and order ",string(maxTerm), ' correction'))
h = xlabel('Coupling Strength ($s^{-1}$)');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('Fidelity');
set(h,'interpreter','latex','fontsize',14);
h = title(strcat(sequenceName,'_',string(Delta)));
set(h,'interpreter','latex','fontsize', 16);

%% Fit ^^ to 1-(x/w0)^a 
aW0 = fittype('b-(W/W0)^a','independent','W');
ff0 = fit(testVars,results_f(:,1),aW0,'Lower',[0 1 0],'StartPoint',[2500,3 1]) %#ok<NOPTS>
hold on
plot(testVars,results_f(:,1))
plot(ff0)

%% Plot Time-Adjusted Infidelity vs Test Var
tauCount = length(getSequence(sequenceName,X,Y).Taus);
cutoff = tauCount * 10e-6;
taf0 = timeAdjust(cutoff,testVars,results_f(:,1),tauCount);
taDf0 = timeAdjust(cutoff,testVars,results_Df(:,1),tauCount);
taDfn = timeAdjust(cutoff,testVars,results_Df(:,maxTerm+1),tauCount);

hold on
plot(testVars,taf0,'Color',myColors{1})
plot(testVars,taDf0,'Color',myColors{4})
plot(testVars,taDfn,'Color',myColors{5})
plot(fitAW0(testVars,taf0))

%% Plot Term Sizes
xEvens = getEvenTerms(0:maxTerm,0);
yEvens = getEvenTerms(results_hsizes,0);
xOdds = getEvenTerms(0:maxTerm,1);
yOdds = getEvenTerms(results_hsizes,1);

hold on
plot(xEvens,yEvens,'Color',myColors{1})
plot(xOdds,yOdds,'Color',myColors{5})
h = xlabel('Term');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('$||\bar H^{(n)}||$ (log-trace method)');
set(h,'interpreter','latex','fontsize',14);
h = title(strcat(sequenceName,',  ',string(coupling),'/',string(Delta)));
set(h,'interpreter','latex','fontsize', 16);
legend('Even Terms','Odd Terms')

%% Plot Commutator Sizes
xEvens = getEvenTerms(2:maxTerm,0);
yEvens = getEvenTerms(results_CS(3:maxTerm+1),0);
xOdds = getEvenTerms(0:maxTerm,1);
yOdds = getEvenTerms(results_CS,1);

hold on
plot(xEvens,yEvens,'Color',myColors{1})
plot(xOdds,yOdds,'Color',myColors{5})
h = xlabel('Term (only even terms included)');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('Term Commutator Size (log-trace method)');
set(h,'interpreter','latex','fontsize',14);
h = title(strcat(sequenceName,',  ',string(coupling),'/',string(Delta)));
set(h,'interpreter','latex','fontsize', 16);
legend('Even Terms','Odd Terms')

%% Plot Fidelity by Increasing Terms

plot(0:maxTerm,results_Df(ceil(length(testVars)/2),:))
h = xlabel('Number of Included AHT Terms');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('Fidelity');
set(h,'interpreter','latex','fontsize',14);
h = title(strcat(sequenceName,',  ',string(coupling),'/',string(Delta)));
set(h,'interpreter','latex','fontsize', 16);