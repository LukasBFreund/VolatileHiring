%%=========================================================================
% Replication code for den Haan, Freund and Rendahl (2021)
%
% Plotting of results for lower entry probability
% Uses input files created based on the master code (see documentation)
%
% Run on Matlab R2019b, Dynare 4.4.3
% Last updated: November 2020
% For any questions please email me at lukas.beat.freund@gmail.com
%%=========================================================================

%% Housekeeping
%--------------------------------------------------------------------------
clear;
close all;
clc;
load inputSettings

%% IRFs
%--------------------------------------------------------------------------
% Load data
sPlotting.numModels = 2; 

%% P=0.2
load(fullfile('.', 'Inputs/', 'IRFs_SaMOptionValue_Uniform_p02_sigmaa0001_Recalib'));
mIRFProp_1=mIRFProp_zUncertainty_EMAS; 
% For rates (u, h, f, p), switch to ppt 
uPos = strmatch('u',vNames,'exact');
hPos = strmatch('h',vNames,'exact');
pPos = strmatch('p',vNames,'exact');
JU_dPos = strmatch('JU_d',vNames,'exact'); % and here use abs value

mIRFProp_1(:,uPos,:) = mIRFProp_1(:,uPos,:)*vEMAS(uPos);
mIRFProp_1(:,hPos,:) = mIRFProp_1(:,hPos,:)*vEMAS(hPos);
mIRFProp_1(:,pPos,:) = mIRFProp_1(:,pPos,:)*vEMAS(pPos);
mIRFProp_1(:,JU_dPos,:) = mIRFProp_1(:,JU_dPos,:)*vEMAS(JU_dPos);

%% P=0.5
load(fullfile('.', 'Inputs\', 'IRFs_SaMOptionValue_Uniform_p05_sigmaa0001_Recalib'));
mIRFProp_2=mIRFProp_zUncertainty_EMAS; 

% Transform as well
uPos = strmatch('u',vNames,'exact');
hPos = strmatch('h',vNames,'exact');
pPos = strmatch('p',vNames,'exact');
JU_dPos = strmatch('JU_d',vNames,'exact'); % and here use abs value
mIRFProp_2(:,uPos,:) = mIRFProp_2(:,uPos,:)*vEMAS(uPos);
mIRFProp_2(:,hPos,:) = mIRFProp_2(:,hPos,:)*vEMAS(hPos);
mIRFProp_2(:,pPos,:) = mIRFProp_2(:,pPos,:)*vEMAS(pPos);
mIRFProp_2(:,JU_dPos,:) = mIRFProp_2(:,JU_dPos,:)*vEMAS(JU_dPos);

% What to plot and how
sPlotting.IRFPeriods = 10;
sPlotting.vVariables = {'JU_d','u','p','h'};
sPlotting.vVNames = {'Value of unmatched entrepreneur','Unemployment rate', 'Entry probability','Hiring rate'};
sPlotting.numSubplotV = 2;
sPlotting.numSubplotH = 2;
sPlotting.legend.show = 1;
sPlotting.legend.labels{1} = 'p = 0.2';
sPlotting.legend.labels{2} = 'p = 0.5';
sPlotting.legend.position = 'northeast';
   
% Put together
if sPlotting.numModels == 1
sResults.aIRF = reshape([mIRFProp_1],[size(mIRFProp_1,1),size(mIRFProp_1,2),sPlotting.numModels]);
elseif sPlotting.numModels == 2
sResults.aIRF = reshape([mIRFProp_1,mIRFProp_2],[size(mIRFProp_1,1),size(mIRFProp_1,2),sPlotting.numModels]);
elseif sPlotting.numModels == 3
sResults.aIRF = reshape([mIRFProp_1,mIRFProp_2,mIRFProp_3],[size(mIRFProp_1,1),size(mIRFProp_1,2),sPlotting.numModels]);
elseif sPlotting.numModels == 4
sResults.aIRF = reshape([mIRFProp_1,mIRFProp_2,mIRFProp_3,mIRFProp_4],[size(mIRFProp_1,1),size(mIRFProp_1,2),sPlotting.numModels]);
end

% Adjust values that are virtually zero (not exactly due to use of a solver) to exactly zero to avoid confusing graphs
sResults.aIRF(abs(sResults.aIRF)<1e-10)=0;
 
% Plot 
fig1=figure;
for iV = 1:numel(sPlotting.vVariables)
    subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,iV)
for iN = 1:sPlotting.numModels
    box on 
    hold on
    plot(0:sPlotting.IRFPeriods,100*sResults.aIRF(1:sPlotting.IRFPeriods+1,strmatch(sPlotting.vVariables{iV},vNames,'exact'),iN),sSettings.lines.list{iN},'LineWidth',sSettings.lines.width,'Color',sSettings.colors.list{iN});
    hold on
    plot(0:sPlotting.IRFPeriods,zeros(sPlotting.IRFPeriods+1,1),sSettings.lines.zeros,'HandleVisibility','off','Color','k','Linewidth',0.1);
    xlim([0 sPlotting.IRFPeriods]);
    set(gca,'XTick',[0:2:sPlotting.IRFPeriods],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
    ax = gca;     ax.YAxis.Exponent = 0; 
end
    title(sPlotting.vVNames{iV},'FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
    
    if ismember(sPlotting.vVariables(iV),{'u','f','h','p'}) == 1
        ylabel('Deviation (ppts.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    elseif ismember(sPlotting.vVariables(iV),{'aHat','JU_d','JaHat_aHatSS','JU1'}) == 1
        ylabel('Deviation (abs.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    else        
       ylabel('Deviation (pct.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    end    
    if iV >= (sPlotting.numSubplotV*sPlotting.numSubplotH)-1
       xlabel('Time (quarters)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    end
end

if sPlotting.legend.show == 1
subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,1)
legend1 = legend(sPlotting.legend.labels);
set(legend1,'fontname',sSettings.font.name,'Location',sPlotting.legend.position,'FontSize',sSettings.font.size.legend,'interpreter','none')
end

xSize = 2*8.75; ySize = 2*6.25;  xCut = 1; yCut = 0.5;
set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-xCut ySize-yCut],'PaperPositionMode','auto')

print(fig1,'Output\fig_App_P02vsP05_Recalib_sigmaa001','-dpdf','-painters')

%% Figure 2
%--------------------------------------------------------------------------

% Choose variables to plot
sPlotting.variables = {'JU_d','u'}; 
sPlotting.names = {'Value of unmatched entrepreneur','Unemployment rate'};
sPlotting.numSubplotV = 1;
sPlotting.numSubplotH = 2;

% Load data for P=0.2
load('Inputs\resultsOptionValue_p02_Recalib_sigmaa0To002')
vSigma_a = sResults.sigmaass;
vNames = sResults.names;
mResults = sResults.aStore(1,:,1,:); % impact period, pure uncertainty
mResults = reshape(mResults,size(mResults,2),numel(vSigma_a));

% Comparison: effect from sigma_a = 0.001 with P = 0.5
load('Inputs\resultsOptionValue_p05_Recalib_sigmaa0001')
mResults_Comp = sResults.aStore(1,:,1,:); % impact period, pure uncertainty
mResults_Comp = reshape(mResults_Comp,size(mResults_Comp,2),1);

% Plot           
fig2 = figure; 
 for iV = 1:numel(sPlotting.variables)
    subplot(sPlotting.numSubplotV,sPlotting.numSubplotH,iV)
    p1=plot(vSigma_a,100*mResults(strmatch(sPlotting.variables{iV},vNames,'exact'),:),sSettings.lines.list{1},'LineWidth',sSettings.lines.width,'Color',sSettings.colors.list{1});  
    hold on 
    p2=plot(vSigma_a,100*ones(1,numel(vSigma_a))*mResults_Comp(strmatch(sPlotting.variables{iV},vNames,'exact'),:),sSettings.lines.list{2},'LineWidth',1,'Color',sSettings.colors.list{2});  
    xlabel('$\sigma_a$','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name,'interpreter','latex');
    set(gca,'XTick',[vSigma_a(1):((vSigma_a(end)-vSigma_a(1))/4):vSigma_a(end)],'FontSize',sSettings.font.size.axisticks,'fontname',sSettings.font.name);
    ax = gca; ax.YAxis.Exponent = 0; 

    if ismember(sPlotting.variables(iV),{'u','f','h','p'}) == 1
        ylabel('Deviation (ppts.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    elseif ismember(sPlotting.variables(iV),{'aHat','JU_d','JaHat_aHatSS'}) == 1
        ylabel('Deviation (abs.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    else        
       ylabel('Deviation (pct.)','FontSize',sSettings.font.size.axis,'fontname',sSettings.font.name);
    end    
    box on
    title(sPlotting.names{iV},'FontSize',sSettings.font.size.default,'fontname',sSettings.font.name,'FontWeight','normal');
    axis tight
 end    
 
 
xSize = 2*8.75; ySize = 1*6.25;  xCut = 1; yCut = 0;
set(gcf,'Units','centimeters','Position',[0 0 xSize ySize],'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 xSize ySize],'PaperSize',[xSize-xCut ySize-yCut],'PaperPositionMode','auto')

print(fig2,'Output\fig_App_Recalib_Loop0to002_WithComparison','-dpdf','-painters')
