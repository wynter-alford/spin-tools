%% plotCorrections.m
% Wynter Alford
% February 2022
%
% Code to plot various types of results.  Use MATLAB's Run Section feature
% and select the plot type you'd like to make.

%% Colors

myColors = {[250 190 212],[230 25 75],[245 130 48],  [210 245 60], [60 180 75], [70 240 240], [0 130 200], [145 30 180], [0 0 0]};
% colors are:     pink          red      orange           lime          green        cyan         blue         violet       black  


for i=1:length(myColors)
    myColors{i} = myColors{i}/255;    
end

%% Plot Fidelity vs Tau
hold on
%plot(testVars,results_f(:,1),'Color',myColors{1})
plot(testVars,results_Df(:,1),'Color',myColors{1},'LineWidth',1.4)
plot(testVars,results_Df(:,3),'Color',myColors{4},'LineWidth',1.4,'LineStyle','--')
plot(testVars,results_Df(:,maxTerm+1),'Color',myColors{6},'LineWidth',1.4,'LineStyle',':')
%plot(fitAW0(testVars,results_Df(:,maxTerm+1)))

legend('H0','H0+...+H2',strcat('H0+...+H',string(maxTerm)))
h = xlabel('Tau Spacing (s)');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('Fidelity');
set(h,'interpreter','latex','fontsize',14);
%h = title(strcat(sequenceName,'_',string(Delta)));
h = title(sequenceName);
set(h,'interpreter','latex','fontsize', 16);
ax = gca;
ax.XAxis.Exponent = -6;

%% Plot Infidelity vs Tau
hold on
fExp = 72/sum(getSequence(sequenceName,1,1).Taus);
defColors
plot(testVars,-log10(1-results_Df(:,1).^fExp),'Color',myColors{9},'LineWidth',3)
plot(testVars,-log10(1-results_Df(:,3).^fExp),'Color',myColors{2},'LineWidth',3, 'LineStyle','--')
plot(testVars,-log10(1-results_Df(:,maxTerm+1).^fExp),'Color',myColors{6},'LineStyle',':','LineWidth',2.5)

plot(testVars,-log10((results_Dd(:,1)).^fExp),'LineWidth',3)
plot(testVars,-log10((results_Dd(:,3)).^fExp),'LineWidth',3)
plot(testVars,-log10((results_Dd(:,maxTerm+1)).^fExp),'LineWidth',3)

legend('f0','f2',strcat('f',string(maxTerm)))
h = xlabel('Tau Spacing (s)');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('Adjusted Log Metric');
set(h,'interpreter','latex','fontsize',14);
h = title(sequenceName);
set(h,'interpreter','latex','fontsize', 16);
ax = gca;
ax.XAxis.Exponent = -6;
ax.FontSize=16;

%% Plot Fidelity vs Coupling
hold on
%plot(testVars,results_f(:,1),'Color',myColors{1})
plot(testVars,results_Df(:,1),'Color',myColors{1},'LineWidth',1.4)
plot(testVars,results_Df(:,maxTerm+1),'Color',myColors{5},'LineStyle','--','LineWidth',1.4)
%plot(fitAW0(testVars,results_Df(:,maxTerm+1)))
legend('H0',strcat('H0+...+H',string(maxTerm)))
h = xlabel('Coupling Strength ($s^{-1}$)');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('Fidelity');
set(h,'interpreter','latex','fontsize',14);
h = title(sequenceName);
%h = title(strcat(sequenceName,'_',string(Delta)));
set(h,'interpreter','latex','fontsize', 16);

%% Plot Infidelity vs Coupling
hold on
fExp = 72/sum(getSequence(sequenceName,1,1).Taus);
defColors
plot(testVars,-log10(1-results_Df(:,1).^fExp),'Color',myColors{9},'LineWidth',3)
plot(testVars,-log10(1-results_Df(:,3).^fExp),'Color',myColors{2},'LineWidth',3, 'LineStyle','--')
plot(testVars,-log10(1-results_Df(:,maxTerm+1).^fExp),'Color',myColors{6},'LineStyle',':','LineWidth',2.5)
ylabel('Adjusted Log-Overlap'); % FIX FOR DIFFERENT NORMS

yyaxis right
plot(testVars,-log10((results_Dd(:,1)).^fExp),'LineWidth',3,'Color',myColors{7})
plot(testVars,-log10((results_Dd(:,3)).^fExp),'LineWidth',3,'Color',myColors{1},'LineStyle','--')
plot(testVars,-log10((results_Dd(:,maxTerm+1)).^fExp),'LineWidth',2.5,'Color',myColors{8},'LineStyle',':')
ylabel('Adjusted Log-Distance');

legend('f0','f2',strcat('f',string(maxTerm)),'d0','d2',strcat('d',string(maxTerm)))
ax = gca;
ax.FontSize = 16;
h = xlabel('$\omega_D$ ($s^{-1}$)');
set(h,'interpreter','latex','FontSize',21);
title(sequenceName);

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
yEvens = getEvenTerms(results_hsizes(1:maxTerm+1),0);
xOdds = getEvenTerms(0:maxTerm,1);
yOdds = getEvenTerms(results_hsizes(1:maxTerm+1),1);

hold on
plot(xEvens,yEvens,'Color',myColors{2},'LineWidth',2)
plot(xOdds,yOdds,'Color',myColors{7},'LineWidth',2,'LineStyle','--')
h = xlabel('Term (n)');
h = ylabel('$\log||\bar H^{(n)}||$)');
set(h,'interpreter','latex','fontsize',14);
h = title(sequenceName);
legend('Even Terms','Odd Terms')
ax = gca;
ax.FontSize = 16; 

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

plot(0:maxTerm,results_Df(ceil(length(testVars)/2),:),'LineWidth',2,'Color',myColors{2})
xlabel('Number of Included AHT Terms (n)');
ylabel('Fidelity');
%set(h,'interpreter','latex','fontsize', 16);
%h = title(strcat(sequenceName,',  ',string(coupling),'/',string(Delta)));
%h = title(strcat(sequenceName, " Fidelity by terms"));
title(sequenceName)

ax = gca;
ax.FontSize = 16; 
%% Plot Fidelity of Multiple Sequences
%toLoad = 'coupling_lo3,4';
%toLoad = 'coupling_lo7,4';
toLoad = 'tau5000';

plotSeqNames = {'WHH','MREV8','BR24','MG8','ML10','YXX24','CORY48','YXX48'};
termNs = [8 6 4 6 6 4 2 2];
hold on

for plotind = 1:length(plotSeqNames)
    load(strcat("31-Jan-2022_", plotSeqNames{plotind}, '_', toLoad, ']0_REC_',string(termNs(plotind)),'.mat'))
    plot(testVars,results_f(:,1),'Color',myColors{plotind})
end
legend(plotSeqNames);

%% Plot Infidelity of Multiple Sequences
%toLoad = 'coupling_lo3,4';
toLoad = 'coupling_lo7,4';
xLabelName = "$\omega_D$ (Hz)";
%toLoad = 'tau5000';
%xLabelName = '$\tau$ Spacing (s)';

load('MarchMeetingPlots.mat')

xmode = 'log'; % can use 'lin' or 'log'
% plotSeqNames = {'WHH','MG8','BR24','CORY48','SED48'};
% dateInds = [1 1 1 1 2];
% dataDates = {'13-Feb-2022_', '18-Feb-2022_', '19-Feb-2022_'};
% styleInds = [1 1 2 2 3];
% styles = ["-","--",":"];
% markers = ['v' '^' 'o' 'o' 'p'];
% colorInds = [1 2 5 7 9];
% termNs = [1 1 1 1 0];
hold on

for plotind = 1:length(plotSeqNames)
    load(strcat(dataDates(dateInds(plotind)), plotSeqNames{plotind}, '_', toLoad, ']0_REC_',string(termNs(plotind)),'.mat'))
    defColors
    fExp = 72/sum(getSequence(plotSeqNames{plotind},1,1).Taus);
    if strcmp(xmode,'lin')
        plot(testVars,-log10(1-(results_Df(:,1).^fExp)),'Color',myColors{colorInds(plotind)},'lineStyle',styles(styleInds(plotind)),'Marker',markers(plotind),'LineWidth',2)
    elseif strcmp(xmode,'log')
        plot(log10(testVars),-log10(1-(results_Df(:,1).^fExp)),'Color',myColors{colorInds(plotind)},'lineStyle',styles(styleInds(plotind)),'Marker',markers(plotind),'LineWidth',2)
    end
end

legend(plotSeqNames);
h = xlabel(xLabelName);
set(h,'interpreter','Latex');
ylabel("Adjusted Log-Infidelity");

ax = gca;
ax.FontSize = 16;
if strcmp(toLoad(1:3),'tau')
ax.XAxis.Exponent = -6;
end

%% Plot Term Size of Multiple Sequences

toLoad = 'coupling_lo3,4';
dateLd = '08-Feb-2022';
%toLoad = 'coupling_lo7,4';
xLabelName = "$\omega_D$ (Hz)";
%toLoad = 'tau5000';
%xLabelName = '$\tau$ Spacing (s)';

plotSeqNames = {'WHH','MG8'};
termNs = [70 25];
hold on

for plotind = 1:length(plotSeqNames)
    load(strcat(dateLd,"_", plotSeqNames{plotind}, '_', toLoad, ']0_REC_',string(termNs(plotind)),'.mat'))
    fExp = 72/sum(getSequence(plotSeqNames{plotind},1,1).Taus);
    plot(0:termNs(plotind),-log10(1-(results_Df(:,1).^fExp)),'Color',myColors{plotind})
end
legend('WHH4',plotSeqNames{2:length(plotSeqNames)});

h = xlabel(xLabelName);
set(h,'interpreter','Latex');
h = ylabel("Log-Infidelity");
set(h,'interpreter','Latex');

if strcmp(toLoad(1:3),'tau')
ax = gca;
ax.XAxis.Exponent = -6;
end

%% Plot mDiff vs Coupling
hold on
plot(testVars,results_fDiff(:,1),'Color',myColors{5})
plot(testVars,results_fDiff(:,3),'Color',myColors{3},'LineWidth',1.4)
plot(testVars,results_fDiff(:,maxTerm+1),'Color',myColors{7},'LineStyle','--','LineWidth',1.4)
%plot(fitAW0(testVars,results_Df(:,maxTerm+1)))
legend('H0','H2',strcat('H0+...+H',string(maxTerm)))
h = xlabel('Coupling Strength ($s^{-1}$)');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('AbsTr-ReTr');
set(h,'interpreter','latex','fontsize',14);
h = title(sequenceName);
%h = title(strcat(sequenceName,'_',string(Delta)));
set(h,'interpreter','latex','fontsize', 16);

%% Plot Infidelity vs LogCoupling
hold on
fExp = 72/sum(getSequence(sequenceName,1,1).Taus);
defColors
plot(log10(testVars),-log10(1-results_Df(:,1)),'Color',myColors{9},'LineWidth',3)
%plot(log10(testVars),-log10(1-results_Df(:,3).^fExp),'Color',myColors{2},'LineWidth',3, 'LineStyle','--')
%plot(log10(testVars),-log10(1-results_Df(:,maxTerm+1).^fExp),'Color',myColors{6},'LineStyle',':','LineWidth',2.5)
legend('f0','f2',strcat('f',string(maxTerm)))
ax = gca;
ax.FontSize = 16;
h = xlabel('$\omega_D$ ($s^{-1}$)');
set(h,'interpreter','latex','FontSize',21);
ylabel('Adjusted Log-Infidelity');
title(sequenceName);

f=fit(log10(testVars(3:40)),-log10(1-results_Df(3:40,1)),'poly1')
plot(f)