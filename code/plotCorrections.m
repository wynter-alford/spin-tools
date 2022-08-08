%% plotCorrections.m
% Wynter Alford
% February 2022
%
% Code to plot various types of results.  Use MATLAB's Run Section feature
% and select the plot type you'd like to make.

%% Colors
defColors

%% Plot Infidelity of Multiple Sequences
loadTag = '-finite';

%toLoad = 'coupling_lo3,4';
toLoad = 'coupling_lo7,4';
xLabelName = "Dipolar Coupling Strength (Hz)";%, Log Scale)";
%toLoad = 'tau5000';
%toLoad = 'tau_lo5000'
%xLabelName = 'Tau Spacing (s)';

%load('MarchMeetingPlots.mat')

xmode = 'lin'; % can use 'lin' or 'log'
fidType = 'Dd'; % can use 'Df', 'f', 'Dd','d' (delta pulse overlap, finite pulse overlap, delta pulse distance, finite pulse distance)
 plotSeqNames = {'WHH','MREV8','BR24','CORY48','SED48'};%,'YXX24'};
 dateInds = [1 1 1 1 1 1 1 1];
 dataDates = {'2022-07-25/25-Jul-2022_'};
 styleInds = [1 1 1 1 1 1 2 2];
 styles = ["-","--",":"];
 markers = ['o' 'o' 'o' 'o' 'o' 'o' 'o' 'o'];
 colorInds = [1 2 5 7 8];
 termNs = [0 0 0 0 0 0 0 0];
hold on

for plotind = 1:length(plotSeqNames)
    load(strcat(dataDates(dateInds(plotind)), plotSeqNames{plotind}, '_', toLoad, ']0_REC_',string(termNs(plotind)),loadTag,'.mat'))
    defColors
    fExp = 72/sum(getSequence(plotSeqNames{plotind},1,1).Taus);
    
    if strcmp(fidType,'Df')
        plotResult = results_OIT(:,1);
        ylabel("Adjusted Log-Overlap");
    elseif strcmp(fidType,'f')
        plotResult = results_OFT(:,1);
        ylabel("Adjusted Log-Overlap");
    elseif strcmp(fidType,'Dd')
        plotResult = 1-results_DIT(:,1);
        ylabel("Adjusted Log-Distance");
    elseif strcmp(fidType,'d')
        plotResult = 1-results_OFT(:,1);
        ylabel("Adjusted Log-Distance");
    else
        error("Invalid fidType. Choose one of 'Df', 'f', 'Dd','d' (delta pulse overlap, finite pulse overlap, delta pulse distance, finite pulse distance)")
    end
    
    if strcmp(xmode,'lin')
        plot(testVars,-log10(1-(plotResult.^fExp)),'Color',myColors{colorInds(plotind)},'lineStyle',styles(styleInds(plotind)),'Marker',markers(plotind),'LineWidth',2)
        xlabel(xLabelName);
    elseif strcmp(xmode,'log')
        plot(log10(testVars),-log10(1-(plotResult.^fExp)),'Color',myColors{colorInds(plotind)},'lineStyle',styles(styleInds(plotind)),'Marker',markers(plotind),'LineWidth',2)
        xlabel(strcat(xLabelName, " (Log Scale)"))
    end
end

legend(plotSeqNames);
%yticks([3 6 9 12 15])
%ylim([0 15])
%xlim([2.5,4])
ax = gca;
ax.FontSize = 16;
if strcmp(toLoad(1:3),'tau')&&strcmp(xmode,'lin')
    ax.XAxis.Exponent = -6;
end

%% Plot Infidelity vs Coupling
hold on
fExp = 72/sum(getSequence(sequenceName,1,1).Taus);
defColors
plot(testVars,-log10(1-results_OIT(:,1).^fExp),'Color',myColors{9},'LineWidth',3)
plot(testVars,-log10(1-results_OIT(:,3).^fExp),'Color',myColors{2},'LineWidth',3, 'LineStyle','--')
plot(testVars,-log10(1-results_OIT(:,maxTerm+1).^fExp),'Color',myColors{6},'LineStyle',':','LineWidth',2.5)
ylabel('Adjusted Log-Overlap'); % FIX FOR DIFFERENT NORMS

yyaxis right
plot(testVars,-log10(1-(1-results_DIT(:,1)).^fExp),'LineWidth',3,'Color',myColors{7})
plot(testVars,-log10(1-(1-results_DIT(:,3)).^fExp),'LineWidth',3,'Color',myColors{1},'LineStyle','--')
plot(testVars,-log10(1-(1-results_DIT(:,maxTerm+1)).^fExp),'LineWidth',2.5,'Color',myColors{8},'LineStyle',':')
ylabel('Adjusted Log-Distance');

legend('f0','f2',strcat('f',string(maxTerm)),'d0','d2',strcat('d',string(maxTerm)))
ax = gca;
ax.FontSize = 16;
h = xlabel('$\omega_D$ ($s^{-1}$)');
set(h,'interpreter','latex','FontSize',21);
title(sequenceName);

%% Plot Infidelity vs Tau
hold on
fExp = 72/sum(getSequence(sequenceName,1,1).Taus);
defColors
plot(testVars,-log10(1-results_OIT(:,1).^fExp),'Color',myColors{9},'LineWidth',3)
plot(testVars,-log10(1-results_OIT(:,3).^fExp),'Color',myColors{2},'LineWidth',3, 'LineStyle','--')
plot(testVars,-log10(1-results_OIT(:,maxTerm+1).^fExp),'Color',myColors{6},'LineStyle',':','LineWidth',2.5)

% yyaxis right
% plot(testVars,-log10((results_DIT(:,1)).^fExp),'LineWidth',3)
% plot(testVars,-log10((results_DIT(:,3)).^fExp),'LineWidth',3)
% plot(testVars,-log10((results_DIT(:,maxTerm+1)).^fExp),'LineWidth',3)

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

%% Plot Fidelity by Increasing Terms

plot(0:maxTerm,results_OIT(ceil(length(testVars)/2),:),'LineWidth',2,'Color',myColors{2})
xlabel('Number of Included AHT Terms (n)');
ylabel('Overlap (f_n)');

yyaxis right
plot(0:maxTerm,1-results_DIT(ceil(length(testVars)/2),:),'LineWidth',2,'Color',myColors{7},'LineStyle','--')
ylabel('Propagator Closeness (1-d_n)')
%set(h,'interpreter','latex','fontsize', 16);
%h = title(strcat(sequenceName,',  ',string(coupling),'/',string(Delta)));
%h = title(strcat(sequenceName, " Fidelity by terms"));
title(strcat(sequenceName,", \omega_D\tau = ",string(coupling*tau)));

ax = gca;
ax.XColor = 'k';
ax.YColor = myColors{7};
legend('Overlap','Propagator Closeness')
ax.FontSize = 13; 

%% Plot Term Sizes
xEvens = getEvenTerms(0:maxTerm,0);
yEvens = getEvenTerms(results_HHT(1:maxTerm+1),0);
xOdds = getEvenTerms(0:maxTerm,1);
yOdds = getEvenTerms(results_HHT(1:maxTerm+1),1);

hold on
plot(xEvens,yEvens+log10(tCyc),'Color',myColors{5},'LineWidth',2)
plot(xOdds,yOdds+log10(tCyc),'Color',myColors{8},'LineWidth',2,'LineStyle','--')
xlabel('Term (n)');
h = ylabel('$\log(t_c||\bar H^{(n)}||_2 )$');
set(h,'interpreter','latex','fontsize',14);
title(strcat(sequenceName,", \omega_D\tau = ",string(coupling*tau)));
legend('Even Terms','Odd Terms')
ax = gca;
ax.FontSize = 16; 

%% Plot Fidelity vs Tau
hold on
%plot(testVars,results_OFT(:,1),'Color',myColors{1})
plot(testVars,results_OIT(:,1),'Color',myColors{1},'LineWidth',1.4)
plot(testVars,results_OIT(:,3),'Color',myColors{4},'LineWidth',1.4,'LineStyle','--')
plot(testVars,results_OIT(:,maxTerm+1),'Color',myColors{6},'LineWidth',1.4,'LineStyle',':')
%plot(fitAW0(testVars,results_OIT(:,maxTerm+1)))

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

%% Plot Fidelity vs Coupling
hold on
%plot(testVars,results_OFT(:,1),'Color',myColors{1})
plot(testVars,results_OIT(:,1),'Color',myColors{1},'LineWidth',1.4)
plot(testVars,results_OIT(:,maxTerm+1),'Color',myColors{5},'LineStyle','--','LineWidth',1.4)
%plot(fitAW0(testVars,results_OIT(:,maxTerm+1)))
legend('H0',strcat('H0+...+H',string(maxTerm)))
h = xlabel('Coupling Strength ($s^{-1}$)');
set(h,'interpreter','latex','fontsize',14);
h = ylabel('Fidelity');
set(h,'interpreter','latex','fontsize',14);
h = title(sequenceName);
%h = title(strcat(sequenceName,'_',string(Delta)));
set(h,'interpreter','latex','fontsize', 16);

%% Fit ^^ to 1-(x/w0)^a 
aW0 = fittype('b-(W/W0)^a','independent','W');
ff0 = fit(testVars,results_OFT(:,1),aW0,'Lower',[0 1 0],'StartPoint',[2500,3 1]) %#ok<NOPTS>
hold on
plot(testVars,results_OFT(:,1))
plot(ff0)

%% Plot Time-Adjusted Infidelity vs Test Var
tauCount = length(getSequence(sequenceName,X,Y).Taus);
cutoff = tauCount * 10e-6;
taf0 = timeAdjust(cutoff,testVars,results_OFT(:,1),tauCount);
taDf0 = timeAdjust(cutoff,testVars,results_OIT(:,1),tauCount);
taDfn = timeAdjust(cutoff,testVars,results_OIT(:,maxTerm+1),tauCount);

hold on
plot(testVars,taf0,'Color',myColors{1})
plot(testVars,taDf0,'Color',myColors{4})
plot(testVars,taDfn,'Color',myColors{5})
plot(fitAW0(testVars,taf0))

%% Plot Fidelity of Multiple Sequences
%toLoad = 'coupling_lo3,4';
%toLoad = 'coupling_lo7,4';
toLoad = 'tau5000';

plotSeqNames = {'WHH','MREV8','BR24','MG8','ML10','YXX24','CORY48','YXX48'};
termNs = [8 6 4 6 6 4 2 2];
hold on

for plotind = 1:length(plotSeqNames)
    load(strcat("31-Jan-2022_", plotSeqNames{plotind}, '_', toLoad, ']0_REC_',string(termNs(plotind)),'.mat'))
    plot(testVars,results_OFT(:,1),'Color',myColors{plotind})
end
legend(plotSeqNames);


%% Plot Term Size vs N

even_sizes = getEvenTerms(results_HHT,0);
odd_sizes = getEvenTerms(results_HHT,1);

hold on
plot(2*(0:floor(maxTerm/2)),even_sizes,'LineWidth',2,'Color',myColors{2})
plot(2*(1:floor(maxTerm/2))-1,odd_sizes,'LineWidth',2,'Color',myColors{7},'LineStyle','--')
h=xlabel('$n$');
set(h,'interpreter','Latex')
h=ylabel('$\log||\bar H^{(n)}||_2$');
set(h,'interpreter','Latex')
title(sequenceName)

ax = gca;
legend('Even Terms','Odd Terms')
ax.FontSize = 16; 

%% Plot Term Size of Multiple Sequences

toLoad = 'coupling_lo7,4';
dateLd = '12-Jul-2022';
%toLoad = 'coupling_lo7,4';
xLabelName = "\omega_D (Hz)";
%toLoad = 'tau5000';
%xLabelName = '$\tau$ Spacing (s)';

plotSeqNames = {'WHH'};
termNs = [95];
hold on

for plotind = 1:length(plotSeqNames)
    load(strcat(dateLd,"_", plotSeqNames{plotind}, '_', toLoad, ']0_REC_',string(termNs(plotind)),'.mat'))
    plot(0:termNs(plotind),results_HHT,'Color',myColors{plotind+1})
end
legend('WHH4',plotSeqNames{2:length(plotSeqNames)});

h = xlabel(xLabelName);
h = ylabel("log||H^{(n)}||_2");

if strcmp(toLoad(1:3),'tau')
ax = gca;
ax.XAxis.Exponent = -6;
end

%% Plot Infidelity vs LogCoupling
hold on
fExp = 72/sum(getSequence(sequenceName,1,1).Taus);
defColors
plot(log10(testVars),-log10(1-results_OIT(:,1)),'Color',myColors{9},'LineWidth',3)
%plot(log10(testVars),-log10(1-results_OIT(:,3).^fExp),'Color',myColors{2},'LineWidth',3, 'LineStyle','--')
%plot(log10(testVars),-log10(1-results_OIT(:,maxTerm+1).^fExp),'Color',myColors{6},'LineStyle',':','LineWidth',2.5)
legend('f0','f2',strcat('f',string(maxTerm)))
ax = gca;
ax.FontSize = 16;
h = xlabel('$\omega_D$ ($s^{-1}$)');
set(h,'interpreter','latex','FontSize',21);
ylabel('Adjusted Log-Overlap');
title(sequenceName);

f=fit(log10(testVars(3:40)),-log10(1-results_OIT(3:40,1)),'poly1')
plot(f)